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

use rayon::prelude::*;

use super::fit::SampleFit;
use super::line::{RewriteTally, rewrite_line, site_of};
use super::model::GenotypePairs;
use super::phase::{
    DISAGREEMENT_COST, EXTRA_PATHS_PER_ROUND, JoinRecord, LEAN_MIN, LastBlock, PassReadCounts,
    PathAgreement, SHARED_RECORDS_TO_JOIN, SWITCH_COST, WrittenPhase,
};
use super::records::SiteRecord;
use super::refine::{
    FillCheck, HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW, MAX_HIDDEN_SHARE, MISMATCH, NoReadRule,
    PHASING_ROUNDS, Phasing, RefinedRecord, record_window_bytes, refine_window_with_reads,
    window_key,
};
use super::templates::MAX_TEMPLATES;
use super::windows::{WindowBytes, WindowSizes, plan_window};
use crate::read_phase::fragments::{PackedRecordFragments, RecordFragments};
use crate::types::{ContigId, Ploidy};
use crate::vcf::RecordPlace;

/// The fewest samples the pass runs on by default (`N_min`, §5) — a starting value for §8.
pub const MIN_SAMPLES: usize = 10;

/// The lowest median inbreeding the pass runs on by default (`F_min`, §5) — a starting value.
pub const MIN_MEDIAN_INBREEDING: f64 = 0.5;

/// **`--refine-genotypes`**: whether the run decides, or the operator did (§5, ruled §9 Q3).
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
    /// The share of drawn paths that must agree at a phase link (`--phase-min-path-agreement`,
    /// phase output spec §9); [`decide`] sets [`PathAgreement::DEFAULT`].
    pub path_agreement: PathAgreement,
    /// Whether each sample's fragments orient its haplotypes in the rounds' draws
    /// (`neighbour_evidence_reads_in_draw.md`); [`decide`] sets it. Off only where a test needs
    /// the link test's reads alone.
    pub(crate) reads_in_draw: bool,
}

impl NeighbourEvidenceDecision {
    /// The same decision, a phase link kept where `agreement` of the drawn paths agree.
    #[must_use]
    pub fn with_path_agreement(self, agreement: PathAgreement) -> Self {
        Self {
            path_agreement: agreement,
            ..self
        }
    }

    /// The same decision, its windows planned within `bytes` (`--refine-window-bytes`).
    #[must_use]
    pub fn with_window_bytes(self, bytes: u64) -> Self {
        Self {
            sizes: WindowSizes {
                budget_bytes: bytes,
                ..self.sizes
            },
            ..self
        }
    }

    /// The value of the `##neighbourEvidence=` header line.
    #[must_use]
    pub fn as_header_value(&self) -> String {
        let why = self.why.replace('"', "'");
        if self.ran {
            format!(
                "<ran=yes,why=\"{why}\",templates={MAX_TEMPLATES},mismatch={MISMATCH},core={},margin={},window_bytes={},\
                 rounds={PHASING_ROUNDS},haplotype_templates={MAX_TEMPLATES},\
                 hidden_per_window={HIDDEN_RECORDS_PER_SAMPLE_PER_WINDOW},\
                 hidden_share={MAX_HIDDEN_SHARE},extra_paths={EXTRA_PATHS_PER_ROUND},\
                 lean_min={LEAN_MIN},path_agreement={},join_records={SHARED_RECORDS_TO_JOIN},\
                 alignment_costs={DISAGREEMENT_COST}/{SWITCH_COST},\
                 read_link_min={:.2},read_break_min={:.2},read_error_floor={},read_gq_min={},\
                 read_reach={},read_cap={}>",
                self.sizes.core_bases,
                self.sizes.margin_bases,
                self.sizes.budget_bytes,
                self.path_agreement.get(),
                crate::read_phase::link::LINK_MIN,
                crate::read_phase::link::BREAK_MIN,
                crate::read_phase::link::ERROR_FLOOR,
                crate::read_phase::link::GQ_MIN,
                crate::read_phase::fragments::READ_REACH,
                crate::read_phase::capture::READ_ALLELES_KEPT,
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
        path_agreement: PathAgreement::DEFAULT,
        reads_in_draw: true,
    }
}

/// **A window the pass cannot solve within its budget** (pass spec §6): at the 500-kb floor it
/// still holds more bytes than `--refine-window-bytes` allows. The run stops here rather than
/// switch refinement off or go on past the budget (owner, 2026-11-13).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WindowOverBudget {
    /// The window's contig, counted from 0 in the header's order, and its core's first position.
    pub contig: u32,
    pub position: u64,
    /// The records the window spans, and the samples each holds.
    pub records: usize,
    pub samples: usize,
    /// The bytes it holds, and the budget.
    pub bytes: u64,
    pub budget: u64,
}

/// **How many held lines are compressed together**, on the pool, once that many have arrived
/// since the last batch: compressing each as it arrives would put the work on the one thread every
/// line passes through.
const LINES_COMPRESSED_TOGETHER: usize = 1_024;

/// **The zstd level held lines are compressed at.** On 800 tomato samples level 1 shrank a line
/// about 5.4 times at about 440 MB/s a thread; levels −3 and 3 gave about the same.
const LINE_COMPRESSION_LEVEL: i32 = 1;

/// **A held line's text**: as it arrived, or compressed. A window holds every line it spans until
/// its core is written, and at 45 bytes a sample a line, uncompressed text was a tenth of a window's
/// memory (pass spec §6).
enum HeldText {
    Plain(Vec<u8>),
    Compressed { bytes: Vec<u8>, length: usize },
}

/// **One thread's compressor and the buffer it compresses into**, made once a batch a thread
/// rather than once a line.
struct LineCompressor {
    compressor: zstd::bulk::Compressor<'static>,
    scratch: Vec<u8>,
}

impl LineCompressor {
    fn new() -> Self {
        Self {
            // PANIC-FREE: a context at a valid level fails only where allocation does.
            compressor: zstd::bulk::Compressor::new(LINE_COMPRESSION_LEVEL)
                .expect("a zstd context at level 1"),
            scratch: Vec::new(),
        }
    }
}

/// One thread's decompression context, made once a core a thread rather than once a line.
fn line_decompressor() -> zstd::bulk::Decompressor<'static> {
    // PANIC-FREE: making a context fails only where allocation does.
    zstd::bulk::Decompressor::new().expect("a zstd context")
}

impl HeldText {
    /// Compresses a plain line; a compressed line, or one already taken, stays as it is.
    ///
    /// **The compressed bytes are copied out at their own length.** zstd compresses into a buffer
    /// of its worst case, at least the line's own length, and a vector keeps its capacity: held as
    /// zstd returns it, a compressed line would hold as much memory as the plain one.
    fn compress(&mut self, with: &mut LineCompressor) {
        let Self::Plain(line) = self else {
            return;
        };
        if line.is_empty() {
            return;
        }
        with.scratch.clear();
        with.scratch
            .reserve(zstd::zstd_safe::compress_bound(line.len()));
        // PANIC-FREE: the buffer holds zstd's bound for the line, so compressing fails only where
        // allocation does.
        with.compressor
            .compress_to_buffer(line.as_slice(), &mut with.scratch)
            .expect("a line in memory compresses");
        *self = Self::Compressed {
            bytes: with.scratch.as_slice().to_vec(),
            length: line.len(),
        };
    }

    /// The line as it arrived, leaving the text empty.
    fn take(&mut self, with: &mut zstd::bulk::Decompressor<'static>) -> Vec<u8> {
        match std::mem::replace(self, Self::Plain(Vec::new())) {
            Self::Plain(line) => line,
            // PANIC-FREE: the bytes are this stage's own compression of a line of that length.
            Self::Compressed { bytes, length } => with
                .decompress(&bytes, length)
                .expect("a line this stage compressed inflates"),
        }
    }
}

/// One line the stage holds.
struct HeldLine {
    place: RecordPlace,
    text: HeldText,
    /// The line's walk record, as an index counted from the run's first site.
    site: Option<usize>,
}

/// **A sample's last written block** (phase output spec §2.4): its phase set, and its records
/// written within a margin of its last one, ascending.
#[derive(Clone, Debug)]
struct WrittenBlock {
    phase_set: u64,
    records: Vec<WrittenRecord>,
}

/// One record of a [`WrittenBlock`]: its run-wide site number, its position, and its phase.
#[derive(Clone, Copy, Debug)]
struct WrittenRecord {
    site: usize,
    position: u64,
    phase: WrittenPhase,
}

/// **The pass, fed one line at a time.**
pub struct NeighbourEvidenceStage {
    samples: usize,
    /// Each contig's length, by contig id, for §4.2's rule at a contig's end.
    contig_lengths: Vec<u64>,
    sizes: WindowSizes,
    path_agreement: PathAgreement,
    reads_in_draw: bool,
    fits: Vec<SampleFit>,
    contig: Option<ContigId>,
    held: Vec<HeldLine>,
    /// The held lines' positions, kept beside them for the window planner.
    positions: Vec<u64>,
    /// What the held lines cost a window, as running totals from the first held line, one longer
    /// than the lines ([`WindowBytes`]): over every line a window spans, and again over its core.
    held_bytes: Vec<u64>,
    core_bytes: Vec<u64>,
    /// The held sites, each with its fragments; the first is site number `first_site`.
    sites: HeldSites,
    first_site: usize,
    next_site: usize,
    /// The first held line not yet written.
    core_start: usize,
    /// Held lines before this one are compressed or written.
    compressed_up_to: usize,
    tally: RewriteTally,
    /// Each sample's fill check, summed over the cores solved so far.
    fill_checks: Vec<FillCheck>,
    /// What reads did to each sample's phase, summed over the cores solved so far.
    read_counts: Vec<PassReadCounts>,
    /// Each sample's last written block on the current contig, kept across windows that write no
    /// phase for it (phase output spec §2.4).
    last_blocks: Vec<Option<WrittenBlock>>,
    /// The first position of every core solved, for the tests that score what crosses a core's
    /// edge.
    #[cfg(test)]
    core_starts: Vec<u64>,
}

impl NeighbourEvidenceStage {
    /// A stage for a run that decided to run the pass, over contigs of these lengths.
    #[must_use]
    pub fn new(decision: &NeighbourEvidenceDecision, contig_lengths: Vec<u64>) -> Self {
        Self {
            samples: decision.inbreeding_by_sample.len(),
            contig_lengths,
            sizes: decision.sizes,
            path_agreement: decision.path_agreement,
            reads_in_draw: decision.reads_in_draw,
            fits: decision
                .inbreeding_by_sample
                .iter()
                .map(|&inbreeding| SampleFit::new(inbreeding))
                .collect(),
            contig: None,
            held: Vec::new(),
            positions: Vec::new(),
            held_bytes: vec![0],
            core_bytes: vec![0],
            sites: HeldSites::default(),
            first_site: 0,
            next_site: 0,
            core_start: 0,
            compressed_up_to: 0,
            tally: RewriteTally::default(),
            fill_checks: vec![FillCheck::default(); decision.inbreeding_by_sample.len()],
            read_counts: vec![PassReadCounts::default(); decision.inbreeding_by_sample.len()],
            last_blocks: vec![None; decision.inbreeding_by_sample.len()],
            #[cfg(test)]
            core_starts: Vec::new(),
        }
    }

    /// What the pass has done so far.
    #[must_use]
    pub fn tally(&self) -> RewriteTally {
        self.tally
    }

    /// Each sample's fill check so far, in the VCF's sample order.
    #[must_use]
    pub fn fill_checks(&self) -> &[FillCheck] {
        &self.fill_checks
    }

    /// What reads did to each sample's phase so far, in the VCF's sample order (read phase spec
    /// §2.7, §2.8).
    #[must_use]
    pub fn read_counts(&self) -> &[PassReadCounts] {
        &self.read_counts
    }

    /// **Takes one line**, its place already checked against the line before it by the writer,
    /// and appends to `ready` every line whose window this one completed. `fragments` are the
    /// line's reads that reach another record (read phase spec §2.4).
    ///
    /// # Errors
    ///
    /// If a window this line completed holds more than the budget even at the smallest size.
    pub fn push(
        &mut self,
        place: RecordPlace,
        line: Vec<u8>,
        fragments: RecordFragments,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) -> Result<(), WindowOverBudget> {
        if self.contig.is_some_and(|contig| contig != place.at.contig) {
            self.finish(ready)?;
        }
        self.contig = Some(place.at.contig);
        // Counted as the line arrived, though most lines are held compressed: what a record costs
        // must not depend on when its batch was compressed.
        let mut held = line.len() as u64;
        let mut core = 0;
        let site = site_of(&line, self.samples).map(|record| {
            let fragments = PackedRecordFragments::from(fragments);
            let (record_held, record_core) = record_window_bytes(&record);
            held += record_held + fragments.heap_bytes() as u64;
            core = record_core;
            self.sites.push(record, fragments);
            self.next_site += 1;
            self.next_site - 1
        });
        let last = |totals: &[u64]| totals.last().copied().unwrap_or(0);
        self.held_bytes.push(last(&self.held_bytes) + held);
        self.core_bytes.push(last(&self.core_bytes) + core);
        self.positions.push(place.at.position.0);
        self.held.push(HeldLine {
            place,
            text: HeldText::Plain(line),
            site,
        });
        self.compressed_up_to = self.compressed_up_to.max(self.core_start);
        if self.held.len() - self.compressed_up_to >= LINES_COMPRESSED_TOGETHER {
            self.held[self.compressed_up_to..]
                .par_iter_mut()
                .for_each_init(LineCompressor::new, |with, held| held.text.compress(with));
            self.compressed_up_to = self.held.len();
        }
        self.release(false, ready)
    }

    /// **Finishes the current contig**: solves its remaining windows and appends their lines to
    /// `ready`. Called at a contig change and at the end of the input.
    ///
    /// # Errors
    ///
    /// As [`Self::push`].
    pub fn finish(
        &mut self,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) -> Result<(), WindowOverBudget> {
        self.release(true, ready)?;
        self.held.clear();
        self.held_bytes = vec![0];
        self.core_bytes = vec![0];
        self.compressed_up_to = 0;
        self.positions.clear();
        self.sites.clear();
        self.first_site = self.next_site;
        self.core_start = 0;
        self.contig = None;
        // A block never crosses a contig.
        self.last_blocks.fill(None);
        Ok(())
    }

    /// Solves every window whose lines have all arrived (or all of them, when the contig is
    /// complete), hands back their cores, and forgets what no later window reads.
    fn release(
        &mut self,
        contig_complete: bool,
        ready: &mut Vec<(RecordPlace, Vec<u8>)>,
    ) -> Result<(), WindowOverBudget> {
        while self.core_start < self.held.len() {
            let Some(plan) = plan_window(
                &self.positions,
                WindowBytes {
                    held: &self.held_bytes,
                    core: &self.core_bytes,
                },
                self.core_start,
                contig_complete,
                &self.sizes,
            ) else {
                return Ok(());
            };
            if plan.over_budget > 0 {
                return Err(WindowOverBudget {
                    contig: self.contig.map_or(0, |contig| contig.0),
                    position: self.positions[plan.core_start],
                    records: plan.records(),
                    samples: self.samples,
                    bytes: plan.bytes,
                    budget: self.sizes.budget_bytes,
                });
            }
            #[cfg(test)]
            self.core_starts.push(self.positions[plan.core_start]);
            let window_sites = self.site_span(plan.context_start..plan.context_end);
            let core_sites = self.site_span(plan.core_start..plan.core_end);
            let refined = if core_sites.is_empty() {
                Vec::new()
            } else {
                let contig_length = self
                    .contig
                    .and_then(|contig| self.contig_lengths.get(contig.0 as usize).copied());
                // The window's contig and its first core site's position key its draws: they depend
                // on the window alone, not on what else the run called.
                let first_core = &self.sites.records[core_sites.start - self.first_site];
                let key = window_key(
                    self.contig.map_or(0, |contig| contig.0),
                    first_core.position,
                );
                // Each sample's last written block, at its records in this window's left margin.
                let joins: Vec<Option<LastBlock>> = self
                    .last_blocks
                    .iter()
                    .map(|last| {
                        last.as_ref().map(|last| LastBlock {
                            phase_set: last.phase_set,
                            records: last
                                .records
                                .iter()
                                .filter(|record| {
                                    (window_sites.start..core_sites.start).contains(&record.site)
                                })
                                .map(|record| JoinRecord {
                                    record: record.site - window_sites.start,
                                    copy_0_allele: record.phase.copy_0_allele,
                                    copy_1_allele: record.phase.copy_1_allele,
                                })
                                .collect(),
                        })
                    })
                    .collect();
                let in_window =
                    window_sites.start - self.first_site..window_sites.end - self.first_site;
                let (records, fragments) = self.sites.window(in_window);
                let solved = refine_window_with_reads(
                    records,
                    Some(fragments),
                    core_sites.start - window_sites.start..core_sites.end - window_sites.start,
                    contig_length,
                    &mut self.fits,
                    MAX_TEMPLATES,
                    NoReadRule::default(),
                    Phasing::rounds_for(key)
                        .with_path_agreement(self.path_agreement)
                        .with_reads_in_draw(self.reads_in_draw),
                    &joins,
                );
                for (sum, check) in self.fill_checks.iter_mut().zip(solved.checks) {
                    *sum += check;
                }
                for (sum, counts) in self.read_counts.iter_mut().zip(solved.read_counts) {
                    *sum += counts;
                }
                self.remember_last_blocks(core_sites.start, &solved.records);
                solved.records
            };
            // **The core's lines are rewritten on the pool and handed back in their own order.**
            // Each rewrite reads only its line and its refined record, so the lines are
            // independent; done on this thread it was a quarter of the writing pass on 59 tomato
            // accessions over 32 Mb, with every other thread waiting for the next window. The
            // tallies are counts, so summing them in any order gives the serial total.
            let rewritten: Vec<(RecordPlace, Option<Vec<u8>>, RewriteTally)> = self.held
                [plan.core_start..plan.core_end]
                .par_iter_mut()
                .map_init(
                    || (GenotypePairs::default(), line_decompressor()),
                    |(pairs, decompressor), held| {
                        let line = held.text.take(decompressor);
                        let mut tally = RewriteTally::default();
                        let written = match held.site {
                            Some(site) => rewrite_line(
                                &line,
                                &refined[site - core_sites.start],
                                pairs,
                                &mut tally,
                            ),
                            None => Some(line),
                        };
                        (held.place, written, tally)
                    },
                )
                .collect();
            for (place, written, tally) in rewritten {
                self.tally.add(tally);
                if let Some(line) = written {
                    ready.push((place, line));
                }
            }
            self.core_start = plan.core_end;
            self.forget_what_no_window_reads();
        }
        Ok(())
    }

    /// **Updates each sample's last written block** from a core's records, the first being site
    /// number `first_site`: the block of its last phased record, extended where it continues the
    /// one before, and cut to the records within a margin of that record — further back, no later
    /// window's margin reaches.
    fn remember_last_blocks(&mut self, first_site: usize, core: &[RefinedRecord]) {
        for (sample, last) in self.last_blocks.iter_mut().enumerate() {
            let phased: Vec<WrittenRecord> = core
                .iter()
                .enumerate()
                .filter_map(|(offset, record)| {
                    record.phases[sample].map(|phase| {
                        let site = first_site + offset;
                        WrittenRecord {
                            site,
                            position: self.sites.records[site - self.first_site].position,
                            phase,
                        }
                    })
                })
                .collect();
            let Some(&final_record) = phased.last() else {
                continue;
            };
            let phase_set = final_record.phase.phase_set;
            let mut records = match last.take() {
                Some(previous) if previous.phase_set == phase_set => previous.records,
                _ => Vec::new(),
            };
            records.extend(
                phased
                    .iter()
                    .filter(|record| record.phase.phase_set == phase_set),
            );
            let reach = final_record
                .position
                .saturating_sub(self.sizes.margin_bases);
            records.retain(|record| record.position >= reach);
            *last = Some(WrittenBlock { phase_set, records });
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
        for totals in [&mut self.held_bytes, &mut self.core_bytes] {
            let dropped = totals[keep_from];
            totals.drain(..keep_from);
            for total in totals.iter_mut() {
                *total -= dropped;
            }
        }
        self.compressed_up_to = self.compressed_up_to.saturating_sub(keep_from);
        self.positions.drain(..keep_from);
        self.sites.drop_first(sites_dropped);
        self.first_site += sites_dropped;
        self.core_start -= keep_from;
    }
}

/// **The held sites' walk records and, beside each, its fragments** (read phase spec §2.7): the
/// reads the link test reads. Two vectors, so a window hands the solve two slices without a copy,
/// but changed only together, here, so one cannot drift a site from the other. A line that is no
/// site keeps none, since the pass phases only sites.
#[derive(Debug, Default)]
struct HeldSites {
    records: Vec<SiteRecord>,
    fragments: Vec<PackedRecordFragments>,
}

impl HeldSites {
    fn push(&mut self, record: SiteRecord, fragments: PackedRecordFragments) {
        self.records.push(record);
        self.fragments.push(fragments);
    }

    fn clear(&mut self) {
        self.records.clear();
        self.fragments.clear();
    }

    /// Forget the first `count` sites.
    fn drop_first(&mut self, count: usize) {
        self.records.drain(..count);
        self.fragments.drain(..count);
    }

    /// The sites in `range`, as the window solve takes them.
    fn window(&self, range: std::ops::Range<usize>) -> (&[SiteRecord], &[PackedRecordFragments]) {
        (&self.records[range.clone()], &self.fragments[range])
    }
}

#[cfg(test)]
mod tests;
