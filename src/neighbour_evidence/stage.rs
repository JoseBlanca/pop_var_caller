//! **The pass as a stage in front of the writer** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §2, §5).
//!
//! The stage takes finished lines one at a time, in the order the writer would write them, and
//! gives them back rewritten, in the same order, once each one's window has been solved. It holds
//! the current contig's lines from the left margin of the window being filled to the last line
//! read; when the lines reaching past a window's right margin have arrived it solves the window,
//! hands back its core and forgets what no later window will read. A new contig, or the end of
//! the input, finishes the contig's last windows.
//!
//! It runs on the thread that writes, so the threads that call loci wait while a window is solved
//! and the stage never holds more than a window: the bounded buffer of §2.1, at its tightest.
//! Whether letting calling run ahead would pay is a measurement (§8.6).
//!
//! **Whether it runs at all is decided once, before the header is written** ([`decide`]): on by
//! default for a diploid cohort of at least ten samples whose median inbreeding is at least 0.5,
//! and the decision, with its reason, goes into the header either way.

use super::fit::SampleFit;
use super::line::{RewriteTally, rewrite_line, site_of};
use super::model::GenotypePairs;
use super::records::SiteRecord;
use super::refine::{MISMATCH, NoReadRule, refine_window};
use super::templates::MAX_TEMPLATES;
use super::windows::{WindowSizes, plan_window};
use crate::types::{ContigId, Ploidy};
use crate::vcf::RecordPlace;

/// The fewest samples the pass runs on by default (`N_min`, §5) — a starting value for §8.
pub const MIN_SAMPLES: usize = 10;

/// The lowest median inbreeding the pass runs on by default (`F_min`, §5) — a starting value.
pub const MIN_MEDIAN_INBREEDING: f64 = 0.5;

/// **`--neighbour-evidence`**: whether the run decides, or the operator did (§5, ruled §9 Q3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NeighbourEvidenceSwitch {
    /// Run where the cohort is large and inbred enough.
    Auto,
    /// Run wherever the pass can: two or more samples, diploid.
    On,
    /// Never run.
    Off,
}

/// **Whether the pass runs in this run, and why** — what the header states (§4.1).
#[derive(Clone, Debug, PartialEq)]
pub struct NeighbourEvidenceDecision {
    /// Whether it runs.
    pub ran: bool,
    /// Why, in words: the switch, the sample count, the median inbreeding, the ploidy.
    pub why: String,
    /// Each sample's inbreeding coefficient from the parameters, in the run's sample order — each
    /// walk's starting `F` (§3.5).
    pub inbreeding_by_sample: Vec<f64>,
    /// The window sizes the run plans with.
    pub sizes: WindowSizes,
}

impl NeighbourEvidenceDecision {
    /// The value of the `##neighbourEvidence=` header line.
    #[must_use]
    pub fn as_header_value(&self) -> String {
        let why = self.why.replace('"', "'");
        if self.ran {
            format!(
                "<ran=yes,why=\"{why}\",templates={MAX_TEMPLATES},mismatch={MISMATCH},core={},margin={}>",
                self.sizes.core_bases, self.sizes.margin_bases
            )
        } else {
            format!("<ran=no,why=\"{why}\">")
        }
    }
}

/// **Decides whether the pass runs** (§5): `switch` from the command line, each sample's
/// inbreeding coefficient in the run's sample order, the run's ploidy.
#[must_use]
pub fn decide(
    switch: NeighbourEvidenceSwitch,
    inbreeding_by_sample: &[f64],
    ploidy: Ploidy,
) -> NeighbourEvidenceDecision {
    let samples = inbreeding_by_sample.len();
    let mut sorted = inbreeding_by_sample.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = match samples {
        0 => 0.0,
        _ if samples % 2 == 1 => sorted[samples / 2],
        _ => f64::midpoint(sorted[samples / 2 - 1], sorted[samples / 2]),
    };
    let cohort = format!("{samples} samples, median inbreeding {median:.2}");
    let (ran, why) = if switch == NeighbourEvidenceSwitch::Off {
        (false, format!("switched off; {cohort}"))
    } else if ploidy.get() != 2 {
        (
            false,
            format!("ploidy {}; the pass models two copies", ploidy.get()),
        )
    } else if samples < 2 {
        (
            false,
            "one sample: no other sample to copy from".to_string(),
        )
    } else if switch == NeighbourEvidenceSwitch::On {
        (true, format!("switched on; {cohort}"))
    } else if samples < MIN_SAMPLES {
        (false, format!("{cohort}; fewer than {MIN_SAMPLES} samples"))
    } else if median < MIN_MEDIAN_INBREEDING {
        (
            false,
            format!("{cohort}; median inbreeding below {MIN_MEDIAN_INBREEDING}"),
        )
    } else {
        (true, cohort)
    };
    NeighbourEvidenceDecision {
        ran,
        why,
        inbreeding_by_sample: inbreeding_by_sample.to_vec(),
        sizes: WindowSizes::default(),
    }
}

/// One line the stage holds.
struct HeldLine {
    place: RecordPlace,
    line: Vec<u8>,
    /// The line's walk record, as an index counted from the run's first site.
    site: Option<usize>,
}

/// **The pass, fed one line at a time.**
pub struct NeighbourEvidenceStage {
    samples: usize,
    /// Each contig's length, by contig id, for §4.2's rule at a contig's end.
    contig_lengths: Vec<u64>,
    sizes: WindowSizes,
    fits: Vec<SampleFit>,
    contig: Option<ContigId>,
    held: Vec<HeldLine>,
    /// The held lines' positions, kept beside them for the window planner.
    positions: Vec<u64>,
    /// The held sites' walk records; `sites[0]` is site number `first_site`.
    sites: Vec<SiteRecord>,
    first_site: usize,
    next_site: usize,
    /// The first held line not yet written.
    core_start: usize,
    pairs: GenotypePairs,
    tally: RewriteTally,
}

impl NeighbourEvidenceStage {
    /// A stage for a run that decided to run the pass, over contigs of these lengths.
    #[must_use]
    pub fn new(decision: &NeighbourEvidenceDecision, contig_lengths: Vec<u64>) -> Self {
        Self {
            samples: decision.inbreeding_by_sample.len(),
            contig_lengths,
            sizes: decision.sizes,
            fits: decision
                .inbreeding_by_sample
                .iter()
                .map(|&inbreeding| SampleFit::new(inbreeding))
                .collect(),
            contig: None,
            held: Vec::new(),
            positions: Vec::new(),
            sites: Vec::new(),
            first_site: 0,
            next_site: 0,
            core_start: 0,
            pairs: GenotypePairs::default(),
            tally: RewriteTally::default(),
        }
    }

    /// What the pass has done so far.
    #[must_use]
    pub fn tally(&self) -> RewriteTally {
        self.tally
    }

    /// **Takes one line**, its place already checked against the line before it by the writer,
    /// and appends to `ready` every line whose window this one completed.
    pub fn push(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) {
        if self.contig.is_some_and(|contig| contig != place.at.contig) {
            self.finish(ready);
        }
        self.contig = Some(place.at.contig);
        let site = site_of(&line, self.samples).map(|record| {
            self.sites.push(record);
            self.next_site += 1;
            self.next_site - 1
        });
        self.positions.push(place.at.position.0);
        self.held.push(HeldLine { place, line, site });
        self.release(false, ready);
    }

    /// **Finishes the current contig**: solves its remaining windows and appends their lines to
    /// `ready`. Called at a contig change and at the end of the input.
    pub fn finish(&mut self, ready: &mut Vec<(RecordPlace, Vec<u8>)>) {
        self.release(true, ready);
        self.held.clear();
        self.positions.clear();
        self.sites.clear();
        self.first_site = self.next_site;
        self.core_start = 0;
        self.contig = None;
    }

    /// Solves every window whose lines have all arrived (or all of them, when the contig is
    /// complete), hands back their cores, and forgets what no later window reads.
    fn release(&mut self, contig_complete: bool, ready: &mut Vec<(RecordPlace, Vec<u8>)>) {
        while self.core_start < self.held.len() {
            let Some(plan) = plan_window(
                &self.positions,
                self.core_start,
                contig_complete,
                self.samples,
                &self.sizes,
            ) else {
                return;
            };
            let window_sites = self.site_span(plan.context_start..plan.context_end);
            let core_sites = self.site_span(plan.core_start..plan.core_end);
            let refined = if core_sites.is_empty() {
                Vec::new()
            } else {
                let contig_length = self
                    .contig
                    .and_then(|contig| self.contig_lengths.get(contig.0 as usize).copied());
                refine_window(
                    &self.sites
                        [window_sites.start - self.first_site..window_sites.end - self.first_site],
                    core_sites.start - window_sites.start..core_sites.end - window_sites.start,
                    contig_length,
                    &mut self.fits,
                    MAX_TEMPLATES,
                    NoReadRule::default(),
                )
                .records
            };
            for held in &mut self.held[plan.core_start..plan.core_end] {
                let line = std::mem::take(&mut held.line);
                match held.site {
                    Some(site) => {
                        if let Some(rewritten) = rewrite_line(
                            &line,
                            &refined[site - core_sites.start],
                            &mut self.pairs,
                            &mut self.tally,
                        ) {
                            ready.push((held.place, rewritten));
                        }
                    }
                    None => ready.push((held.place, line)),
                }
            }
            self.core_start = plan.core_end;
            self.forget_what_no_window_reads();
        }
    }

    /// The sites among held lines `lines`, as a range of site numbers.
    fn site_span(&self, lines: std::ops::Range<usize>) -> std::ops::Range<usize> {
        let first = self.held[lines.clone()].iter().find_map(|held| held.site);
        let last = self.held[lines].iter().rev().find_map(|held| held.site);
        match (first, last) {
            (Some(first), Some(last)) => first..last + 1,
            _ => 0..0,
        }
    }

    /// Drops held lines further back than a margin before the next core's first line.
    fn forget_what_no_window_reads(&mut self) {
        let next = self
            .positions
            .get(self.core_start)
            .or(self.positions.last())
            .copied()
            .unwrap_or(0);
        let keep_from = self.positions[..self.core_start]
            .partition_point(|&position| position < next.saturating_sub(self.sizes.margin_bases));
        if keep_from == 0 {
            return;
        }
        let sites_dropped = self.held[..keep_from]
            .iter()
            .filter(|held| held.site.is_some())
            .count();
        self.held.drain(..keep_from);
        self.positions.drain(..keep_from);
        self.sites.drain(..sites_dropped);
        self.first_site += sites_dropped;
        self.core_start -= keep_from;
    }
}

#[cfg(test)]
mod tests;
