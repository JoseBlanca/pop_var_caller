//! **Why a run decided what it did at chosen loci** — `--explain-loci`
//! (`doc/devel/ng/spec/explain_loci.md`).
//!
//! A run asked to explain the loci overlapping a BED records, at each such locus, the decisions
//! it took between the evidence and the VCF: which candidate alleles it kept and why the others
//! went, which samples it set aside, every genotype's likelihood, the called genotypes, how the
//! partial reads were weighed, the site quality and its penalties, and what became of the
//! record. **It selects what is explained, never what is called**: nothing here changes,
//! reorders or recomputes a decision, and a run that explains nothing pays one question per
//! locus — [`ExplainRegions::covers`] is never even asked when the option is off.
//!
//! # How the rows get here
//!
//! Each step's rows are built where that step's values exist, from values the run computes
//! anyway: selection's at selection, the likelihoods off the calling scratch straight after the
//! call, the site quality off the record the run assembled. They are collected per locus and
//! handed back with the locus's outcome, and [`write_explanations`] sorts and writes them once
//! the run ends — the only point at which every locus's final outcome is known.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

use crate::calling::allele_candidates::{
    CandidateSelectionConfig, LocusSelection, SelectionVerdict,
};
use crate::calling::{
    CallingScratch, CandidateAlleles, GenotypeIdx, GenotypeTable, LocusInference,
    ReferenceBesideLocus, SampleGenotypeCall, partial_row_fits_allele,
};
use crate::fasta::ContigList;
use crate::locus_generation::{LocusKind, LocusLen};
use crate::regions::{BedError, ContigBounds, RegionSet};
use crate::run::cohort_merge::build::{CohortObservation, LocusBorder};
use crate::types::{ContigId, GenomeRegion, Phred, Ploidy};
use crate::vcf::VcfRecord;

/// The layout of the TSV [`write_explanations`] writes. **Changes whenever a column is added,
/// removed or moved, or a step's meaning changes**, so a script can refuse a layout it does not
/// know.
pub const COLUMNS_VERSION: u32 = 1;

/// The columns, in order (`explain_loci.md` §4).
pub const COLUMNS: [&str; 8] = [
    "contig", "start", "end", "step", "sample", "subject", "value", "detail",
];

/// **The loci a run was asked to explain**: the regions of a BED, resolved against the run's
/// contigs.
#[derive(Debug, Clone)]
pub struct ExplainRegions {
    regions: RegionSet,
}

impl ExplainRegions {
    /// Read the BED against the run's contig list.
    ///
    /// # Errors
    ///
    /// The BED cannot be read, names a contig the reference lacks, or holds no region; or a
    /// contig is longer than a BED coordinate can address.
    pub fn from_bed_path(path: &Path, contigs: &ContigList) -> Result<Self, BedError> {
        let bounds: Vec<ContigBounds<'_>> = contigs
            .entries
            .iter()
            .map(|entry| ContigBounds {
                name: &entry.name,
                // A contig past `u32` cannot be named by a BED line either; clamping it means
                // its tail cannot be asked about, which the reader would refuse anyway.
                length: u32::try_from(entry.length).unwrap_or(u32::MAX),
            })
            .collect();
        Ok(Self {
            regions: RegionSet::from_bed_path(path, &bounds)?,
        })
    }

    /// **Whether a locus overlaps a region to explain** — the one question a locus is asked.
    #[must_use]
    pub fn covers(&self, locus: GenomeRegion) -> bool {
        let on_its_contig = self.regions.regions_for(locus.contig.0);
        let first_reaching_it =
            on_its_contig.partition_point(|region| u64::from(region.end) < locus.start.get());
        on_its_contig
            .get(first_reaching_it)
            .is_some_and(|region| u64::from(region.start) <= locus.end.get())
    }
}

/// Which decision a row is about, **in the order a locus's rows are written**
/// (`explain_loci.md` §4.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ExplainStep {
    Locus,
    Allele,
    Reads,
    Partial,
    Callable,
    Likelihood,
    Genotype,
    Site,
    Paralog,
    Outcome,
}

impl ExplainStep {
    /// The word the TSV's `step` column carries.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Locus => "locus",
            Self::Allele => "allele",
            Self::Reads => "reads",
            Self::Partial => "partial",
            Self::Callable => "callable",
            Self::Likelihood => "likelihood",
            Self::Genotype => "genotype",
            Self::Site => "site",
            Self::Paralog => "paralog",
            Self::Outcome => "outcome",
        }
    }
}

/// One value in a row's `detail` — text, or a run sample to be named when the file is written.
#[derive(Debug, Clone, PartialEq)]
pub enum DetailValue {
    Text(String),
    Sample(usize),
}

/// **One fact about one locus.**
#[derive(Debug, Clone, PartialEq)]
pub struct ExplainRow {
    /// The locus, as the run built it.
    pub region: GenomeRegion,
    pub step: ExplainStep,
    /// The run sample the fact is about, or `None` for the whole locus.
    pub sample: Option<usize>,
    pub subject: String,
    pub value: String,
    pub detail: Vec<(&'static str, DetailValue)>,
    /// **Where the record was written**, on a written locus's `outcome` row only: the position
    /// the hidden-paralog filter knows a record by, one base before the locus when a padding
    /// base was put in front of it.
    pub written_at: Option<u64>,
}

impl ExplainRow {
    fn of(region: GenomeRegion, step: ExplainStep, value: impl Into<String>) -> Self {
        Self {
            region,
            step,
            sample: None,
            subject: ".".to_owned(),
            value: value.into(),
            detail: Vec::new(),
            written_at: None,
        }
    }

    fn for_sample(mut self, sample: usize) -> Self {
        self.sample = Some(sample);
        self
    }

    fn about(mut self, subject: impl Into<String>) -> Self {
        self.subject = subject.into();
        self
    }

    fn with(mut self, key: &'static str, value: impl ToString) -> Self {
        self.detail
            .push((key, DetailValue::Text(value.to_string())));
        self
    }

    fn with_sample(mut self, key: &'static str, sample: usize) -> Self {
        self.detail.push((key, DetailValue::Sample(sample)));
        self
    }

    fn written_at(mut self, position: u64) -> Self {
        self.written_at = Some(position);
        self
    }
}

/// An allele's bases as a column can carry them: `-` for the empty allele.
fn spelled(bases: &[u8]) -> String {
    if bases.is_empty() {
        "-".to_owned()
    } else {
        String::from_utf8_lossy(bases).into_owned()
    }
}

/// **The locus, its alleles and what each sample's reads showed** — the rows selection can
/// write, before anything is called.
///
/// An allele the merge assembled is `kept` as the candidate selection numbered it, or it went
/// for one of two reasons on an ordinary locus: **no sample's reads reached the support bar for
/// it** (`below_support`), or **some did and the allele cap cut it** (`cut_by_allele_cap`). That
/// split is recomputed from the reads with selection's own bar
/// (`MinAltReads::reached_by`), because
/// selection keeps it in a scratch buffer. A repeat tract's selection nominates lengths by a
/// different route, so a tract's dropped allele is `not_selected`, or
/// `refused_not_periodic` when the whole tract was.
pub(crate) fn selection_rows(
    rows: &mut Vec<ExplainRow>,
    observation: &CohortObservation,
    selection: &LocusSelection,
    config: &CandidateSelectionConfig,
    reference_beside: ReferenceBesideLocus<'_>,
) {
    let region = observation.region;
    let kind = match observation.kind {
        LocusKind::Generic => "generic",
        LocusKind::Ssr(_) => "repeat_tract",
        LocusKind::SsrBundle => "repeat_bundle",
    };
    rows.push(
        ExplainRow::of(region, ExplainStep::Locus, kind)
            .with("alleles", observation.alleles.len())
            .with("covering_samples", observation.per_sample.len())
            .with("selection", verdict_word(selection.verdict())),
    );

    let compared_reads: Vec<u32> = observation
        .per_sample
        .iter()
        .map(|sample| {
            sample
                .supported
                .iter()
                .map(|row| row.support.num_reads)
                .sum()
        })
        .collect();
    let reads_of = |sample_at: usize, allele: usize| -> u32 {
        observation.per_sample[sample_at]
            .supported
            .iter()
            .filter(|row| row.allele == allele)
            .map(|row| row.support.num_reads)
            .sum()
    };

    for (allele, bases) in observation.alleles.iter().enumerate() {
        let mut cohort_reads = 0_u32;
        let mut best: Option<(usize, u32, f64)> = None;
        let mut cleared_the_bar = false;
        for (sample_at, sample) in observation.per_sample.iter().enumerate() {
            let reads = reads_of(sample_at, allele);
            cohort_reads += reads;
            if reads == 0 {
                continue;
            }
            let share = f64::from(reads) / f64::from(compared_reads[sample_at].max(1));
            if best.is_none_or(|(_, _, best_share)| share > best_share) {
                best = Some((sample.sample, reads, share));
            }
            cleared_the_bar |= config
                .min_allele_support
                .reached_by(reads, compared_reads[sample_at]);
        }
        let value = match selection.remap().candidate_for(allele) {
            Some(candidate) => format!("kept_as_candidate_{}", candidate.0),
            None if matches!(selection.verdict(), SelectionVerdict::NotPeriodic) => {
                "refused_not_periodic".to_owned()
            }
            None if !matches!(observation.kind, LocusKind::Generic) => "not_selected".to_owned(),
            None if cleared_the_bar => "cut_by_allele_cap".to_owned(),
            None => "below_support".to_owned(),
        };
        let mut row = ExplainRow::of(region, ExplainStep::Allele, value)
            .about(spelled(bases))
            .with("cohort_reads", cohort_reads);
        if let Some((sample, reads, share)) = best {
            row = row
                .with_sample("best_sample", sample)
                .with("best_sample_reads", reads)
                .with("best_sample_share", format!("{share:.3}"));
        }
        rows.push(row);
    }

    for sample in &observation.per_sample {
        for (allele, bases) in observation.alleles.iter().enumerate() {
            let groups: Vec<_> = sample
                .supported
                .iter()
                .filter(|row| row.allele == allele)
                .collect();
            if groups.is_empty() {
                continue;
            }
            let reads: u32 = groups.iter().map(|row| row.support.num_reads).sum();
            let forward: u32 = groups.iter().map(|row| row.support.num_fwd).sum();
            rows.push(
                ExplainRow::of(region, ExplainStep::Reads, reads.to_string())
                    .for_sample(sample.sample)
                    .about(spelled(bases))
                    .with("forward", forward)
                    .with("read_groups", groups.len()),
            );
        }
    }

    // **How each partial read was weighed** (`read_likelihoods.md` §5.3): it counts for the
    // genotypes carrying a candidate it fits, and is charged as a sequencing error under the
    // rest. The list of candidates it fits is the whole of what decides which.
    if matches!(observation.kind, LocusKind::Generic) {
        let candidates = selection.alleles();
        let locus_len = LocusLen::from_positions(candidates.reference().len() as u64);
        for sample in &observation.per_sample {
            for partial in &sample.partials {
                let fits: Vec<String> = candidates
                    .iter()
                    .enumerate()
                    .filter(|(_, allele)| {
                        partial_row_fits_allele(partial, allele, reference_beside, locus_len)
                    })
                    .map(|(candidate, _)| candidate.to_string())
                    .collect();
                let witnessed: Vec<String> = partial
                    .witnessed_in_locus
                    .runs()
                    .map(|(start, end)| format!("{start}-{end}"))
                    .collect();
                let row =
                    ExplainRow::of(region, ExplainStep::Partial, partial.num_reads.to_string())
                        .for_sample(sample.sample)
                        .about(spelled(&partial.bases))
                        .with("witnessed", witnessed.join(","));
                // Only on the reads re-read after the merge: they covered the whole locus, and
                // this is what says they were weighed as reads flush to one border.
                let row = match partial.sequence_may_run_on_past {
                    Some(LocusBorder::Right) => row.with("may_run_on_past", "right"),
                    Some(LocusBorder::Left) => row.with("may_run_on_past", "left"),
                    None => row,
                };
                rows.push(row.with(
                    "fits_candidates",
                    if fits.is_empty() {
                        "none".to_owned()
                    } else {
                        fits.join(",")
                    },
                ));
            }
        }
    }

    for (sample_at, sample) in observation.per_sample.iter().enumerate() {
        let leftover = &selection.unmatched()[sample_at];
        let row = if leftover.genotype_must_be_missing() {
            ExplainRow::of(region, ExplainStep::Callable, "set_aside")
                .with(
                    "why",
                    "the allele cap cut an allele this sample's own reads earned",
                )
                .with("earned_reads_cut", leftover.earned_reads_cut_by_the_cap)
        } else {
            ExplainRow::of(region, ExplainStep::Callable, "callable")
        };
        rows.push(
            row.for_sample(sample.sample)
                .with("reads_on_no_candidate", leftover.num_reads),
        );
    }
}

fn verdict_word(verdict: SelectionVerdict) -> String {
    match verdict {
        SelectionVerdict::Selected => "selected".to_owned(),
        SelectionVerdict::Truncated { dropped } => format!("truncated_by_cap_dropping_{dropped}"),
        SelectionVerdict::NotPeriodic => "not_periodic".to_owned(),
    }
}

/// A genotype as a VCF spells it, from how many copies of each candidate it carries.
fn genotype_from_counts(counts: &[u32]) -> String {
    let alleles: Vec<String> = counts
        .iter()
        .enumerate()
        .flat_map(|(allele, &copies)| std::iter::repeat_n(allele.to_string(), copies as usize))
        .collect();
    alleles.join("/")
}

/// **Every sample's genotype likelihoods and its call**, read off the calling scratch straight
/// after the call, and the site row's half the call knows.
///
/// The likelihoods are those the final pass scored with — the rows the scratch holds once
/// `call_locus` returns, one per **callable** sample, over the locus's candidates in the
/// genotype table the call built (rebuilt here the same way; the final pass prunes nothing).
pub(crate) fn calling_rows<S>(
    rows: &mut Vec<ExplainRow>,
    inference: &LocusInference,
    scratch: &CallingScratch<S>,
    ploidy: Ploidy,
) {
    let region = inference.region;
    let candidates: &CandidateAlleles = inference.alleles();
    let table = GenotypeTable::build(ploidy, candidates.len());
    let genotypes = table.view();
    for (row, &run_sample) in scratch.run_sample_of_each_row().iter().enumerate() {
        for (at, likelihood) in scratch.sample_genotype_likelihoods(row).iter().enumerate() {
            let named = genotypes
                .allele_counts_of(GenotypeIdx(at as u32))
                .map_or_else(|| format!("#{at}"), genotype_from_counts);
            rows.push(
                ExplainRow::of(
                    region,
                    ExplainStep::Likelihood,
                    format!("{:.4}", likelihood.0),
                )
                .for_sample(run_sample)
                .about(named),
            );
        }
    }
    for (run_sample, call) in inference.per_sample.iter().enumerate() {
        let row = match call {
            SampleGenotypeCall::Called {
                genotype,
                genotype_quality,
                reads_were_uninformative,
            } => {
                let spelled: Vec<String> = genotype
                    .alleles()
                    .iter()
                    .map(|allele| allele.0.to_string())
                    .collect();
                ExplainRow::of(region, ExplainStep::Genotype, spelled.join("/"))
                    .with("GQ", format!("{:.1}", genotype_quality.get()))
                    .with("reads_were_uninformative", reads_were_uninformative)
            }
            SampleGenotypeCall::Missing => ExplainRow::of(region, ExplainStep::Genotype, "./.")
                .with("why", "set aside, or no reads"),
        };
        rows.push(row.for_sample(run_sample));
    }
    rows.push(
        ExplainRow::of(region, ExplainStep::Site, ".")
            .with(
                "uncorrected_qual",
                format!("{:.1}", inference.uncorrected_site_quality().get()),
            )
            .with("converged", inference.converged)
            .with("passes", inference.passes),
    );
}

/// What became of a locus the caller met: its `outcome` row, and the written record's half of
/// its `site` row.
pub(crate) fn outcome_rows(
    rows: &mut Vec<ExplainRow>,
    region: GenomeRegion,
    outcome: LocusEnd<'_>,
) {
    let (value, why) = match outcome {
        LocusEnd::Written(record) => {
            if let Some(site) = rows
                .iter_mut()
                .find(|row| row.step == ExplainStep::Site && row.region == region)
            {
                site.value = format!("{:.1}", record.site_quality.get());
                if let Some(penalties) = record.artifact_penalties {
                    site.detail.push((
                        "ABPEN",
                        DetailValue::Text(format!("{:.1}", penalties.allele_balance.get())),
                    ));
                    site.detail.push((
                        "SPPEN",
                        DetailValue::Text(format!(
                            "{:.1}",
                            penalties.strand_and_read_position.get()
                        )),
                    ));
                }
                site.detail.push((
                    "filter",
                    DetailValue::Text(record.filter().as_str().to_owned()),
                ));
            }
            let at = crate::vcf::encode::written_position(record);
            rows.push(
                ExplainRow::of(region, ExplainStep::Outcome, "written")
                    .with("why", "handed to the VCF")
                    .written_at(at),
            );
            return;
        }
        LocusEnd::NotWritten => (
            "not_written",
            "no sample was called carrying an alternative allele, or those that were had \
             uninformative reads",
        ),
        LocusEnd::BelowMinimumSiteQuality(quality) => {
            rows.push(
                ExplainRow::of(region, ExplainStep::Outcome, "below_min_site_quality").with(
                    "why",
                    format!(
                        "the site quality after the artifact correction, {:.1}, is below \
                         --min-site-quality",
                        quality.get()
                    ),
                ),
            );
            return;
        }
        LocusEnd::NobodyToCall => (
            "nobody_to_call",
            "every sample was set aside because the allele cap cut an allele its reads earned",
        ),
        LocusEnd::BundleSetAside => (
            "bundle_set_aside",
            "a cluster of repeats without clean flanks is not called",
        ),
        LocusEnd::TractWithoutWholeRepeats => (
            "tract_without_whole_repeats",
            "a candidate tract length holds no whole repeat, which the tract model cannot score",
        ),
        LocusEnd::Failed => ("failed", "the run stopped at this locus"),
    };
    rows.push(ExplainRow::of(region, ExplainStep::Outcome, value).with("why", why));
}

/// How a locus the caller met ended, as [`outcome_rows`] reads it.
#[derive(Debug, Clone, Copy)]
pub(crate) enum LocusEnd<'a> {
    Written(&'a VcfRecord),
    NotWritten,
    BelowMinimumSiteQuality(Phred),
    NobodyToCall,
    BundleSetAside,
    TractWithoutWholeRepeats,
    Failed,
}

/// Why the merge dropped a locus before calling it (`explain_loci.md` §4.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeDrop {
    /// No sample showed enough reads disagreeing with the reference.
    TooQuiet,
    /// Some sample's read group was deeper than `--max-read-group-depth`.
    OverDepthCeiling,
}

/// **The explained loci the merge dropped as too quiet or too deep** — kept beside the merge's
/// observation cache, which every merge path holds, and noted by the merge where it drops them.
///
/// Noting asks [`ExplainRegions::covers`] first, so a run explaining a few regions locks only
/// at the loci it explains; the lock is there because the round driver builds regions on
/// several threads.
#[derive(Debug)]
pub struct DroppedLoci {
    regions: ExplainRegions,
    found: Mutex<Vec<(GenomeRegion, MergeDrop)>>,
}

impl DroppedLoci {
    /// A log for the loci `regions` explains.
    #[must_use]
    pub fn over(regions: ExplainRegions) -> Self {
        Self {
            regions,
            found: Mutex::new(Vec::new()),
        }
    }

    /// Note a dropped locus, if it is explained.
    pub fn note(&self, locus: GenomeRegion, why: MergeDrop) {
        if self.regions.covers(locus) {
            self.found
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .push((locus, why));
        }
    }

    /// The noted loci as rows: a `locus` row and an `outcome` row each.
    #[must_use]
    pub fn into_rows(self) -> Vec<ExplainRow> {
        let found = self
            .found
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let mut rows = Vec::with_capacity(2 * found.len());
        for (locus, why) in found {
            let (value, why) = match why {
                MergeDrop::TooQuiet => (
                    "too_quiet",
                    "no sample showed enough reads disagreeing with the reference for a locus \
                     to be built",
                ),
                MergeDrop::OverDepthCeiling => (
                    "over_depth_ceiling",
                    "some sample's read group is deeper here than --max-read-group-depth, so \
                     the locus is dropped for every sample",
                ),
            };
            rows.push(ExplainRow::of(locus, ExplainStep::Locus, "not_built"));
            rows.push(ExplainRow::of(locus, ExplainStep::Outcome, value).with("why", why));
        }
        rows
    }
}

/// The rows of the explained loci the merge refused for being wider than
/// `--max-cohort-locus-span` — read off the run's own list of them.
#[must_use]
pub fn too_wide_rows(regions: &ExplainRegions, too_wide: &[GenomeRegion]) -> Vec<ExplainRow> {
    let mut rows = Vec::new();
    for &locus in too_wide.iter().filter(|&&locus| regions.covers(locus)) {
        rows.push(ExplainRow::of(locus, ExplainStep::Locus, "not_built"));
        rows.push(
            ExplainRow::of(locus, ExplainStep::Outcome, "too_wide").with(
                "why",
                "the locus is wider than --max-cohort-locus-span, so it is not assembled",
            ),
        );
    }
    rows
}

/// **The hidden-paralog filter's verdicts at the explained loci** — asked of the filter's last
/// pass, which knows a record only by where it was written.
#[derive(Debug, Default)]
pub struct ParalogExplanations {
    wanted: HashSet<(ContigId, u64)>,
    found: Vec<ParalogVerdictAt>,
}

/// What the filter decided about one explained record.
#[derive(Debug, Clone, PartialEq)]
pub struct ParalogVerdictAt {
    pub contig: ContigId,
    /// Where the record was written — its `POS`.
    pub written_at: u64,
    /// The log likelihood ratio, collapsed paralog against ordinary variant; not finite where the
    /// record was not scored.
    pub ratio: f64,
    /// The posterior that the record is a collapsed paralog, where the run's fitted rate gives
    /// one.
    pub posterior: Option<f64>,
    /// Whether the filter flagged it.
    pub flagged: bool,
    /// Whether it was dropped, rather than tagged or kept.
    pub dropped: bool,
}

impl ParalogExplanations {
    /// Ask for the verdicts of every record an explained locus wrote.
    #[must_use]
    pub fn wanting(rows: &[ExplainRow]) -> Self {
        Self {
            wanted: rows
                .iter()
                .filter_map(|row| row.written_at.map(|at| (row.region.contig, at)))
                .collect(),
            found: Vec::new(),
        }
    }

    /// Whether the record written at `position` on `contig` is one an explained locus wrote.
    #[must_use]
    pub fn wants(&self, contig: ContigId, position: u64) -> bool {
        self.wanted.contains(&(contig, position))
    }

    /// Keep one verdict.
    pub fn record(&mut self, verdict: ParalogVerdictAt) {
        self.found.push(verdict);
    }

    /// **Add the verdicts to the rows of the loci that wrote them**: a `paralog` row each, and a
    /// dropped record's `outcome` rewritten to say the filter dropped it.
    pub fn add_to(self, rows: &mut Vec<ExplainRow>) {
        for verdict in self.found {
            let Some(outcome) = rows.iter_mut().find(|row| {
                row.written_at == Some(verdict.written_at) && row.region.contig == verdict.contig
            }) else {
                continue;
            };
            let region = outcome.region;
            if verdict.dropped {
                outcome.value = "dropped_by_paralog_filter".to_owned();
                outcome.detail = vec![(
                    "why",
                    DetailValue::Text(
                        "the hidden-paralog filter flagged it, and flagged records are dropped \
                         unless --paralog-filter-tag is given"
                            .to_owned(),
                    ),
                )];
            }
            let value = match (verdict.ratio.is_finite(), verdict.flagged, verdict.dropped) {
                (false, _, _) => "unscored",
                (true, true, true) => "dropped",
                (true, true, false) => "tagged",
                (true, false, _) => "kept",
            };
            let mut row = ExplainRow::of(region, ExplainStep::Paralog, value)
                .with("LR", format!("{:.4}", verdict.ratio));
            if let Some(posterior) = verdict.posterior {
                row = row.with("posterior", format!("{posterior:.6}"));
            }
            rows.push(row);
        }
    }
}

/// **Sort a run's rows into genome order and write them**, with sample numbers turned into
/// names.
///
/// # Errors
///
/// The file cannot be written.
pub fn write_explanations(
    path: &Path,
    rows: &mut [ExplainRow],
    contigs: &ContigList,
    sample_names: &[String],
) -> io::Result<()> {
    // Stable, so the rows of one step keep the order they were built in.
    rows.sort_by_key(|row| {
        (
            row.region.contig,
            row.region.start,
            row.region.end,
            row.step,
        )
    });
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "#explain_loci_columns_version={COLUMNS_VERSION}")?;
    writeln!(
        out,
        "#coordinates=1-based, both ends included, as in a VCF; candidate numbers are the \
         allele numbers a written record's GT uses before alleles no genotype names are dropped"
    )?;
    writeln!(
        out,
        "#likelihood=natural log of the probability of the sample's reads given the genotype"
    )?;
    writeln!(out, "{}", COLUMNS.join("\t"))?;
    let name_of = |sample: usize| {
        sample_names
            .get(sample)
            .cloned()
            .unwrap_or_else(|| format!("#{sample}"))
    };
    for row in rows.iter() {
        let contig = contigs
            .entries
            .get(row.region.contig.0 as usize)
            .map_or("?", |entry| entry.name.as_str());
        let detail = if row.detail.is_empty() {
            ".".to_owned()
        } else {
            row.detail
                .iter()
                .map(|(key, value)| match value {
                    DetailValue::Text(text) => format!("{key}={text}"),
                    DetailValue::Sample(sample) => format!("{key}={}", name_of(*sample)),
                })
                .collect::<Vec<_>>()
                .join(";")
        };
        writeln!(
            out,
            "{contig}\t{}\t{}\t{}\t{}\t{}\t{}\t{detail}",
            row.region.start.get(),
            row.region.end.get(),
            row.step.as_str(),
            row.sample.map_or_else(|| ".".to_owned(), name_of),
            row.subject,
            row.value,
        )?;
    }
    out.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_genotype_is_spelled_from_its_copy_counts() {
        assert_eq!(genotype_from_counts(&[2, 0]), "0/0");
        assert_eq!(genotype_from_counts(&[1, 1]), "0/1");
        assert_eq!(genotype_from_counts(&[0, 1, 1]), "1/2");
    }

    use crate::types::Position;

    fn regions(bed: &str) -> ExplainRegions {
        let bounds = [
            ContigBounds {
                name: "chr1",
                length: 1_000,
            },
            ContigBounds {
                name: "chr2",
                length: 1_000,
            },
        ];
        ExplainRegions {
            regions: RegionSet::from_bed_reader(bed.as_bytes(), &bounds).expect("the BED reads"),
        }
    }

    fn locus(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **A locus is explained when it shares a base with a region**: one reaching into a region
    /// from before it is, one ending a base short of it is not, and a region on another contig
    /// explains nothing here.
    #[test]
    fn a_locus_is_explained_when_it_overlaps_a_region() {
        // BED 100-200 is bases 101..=200.
        let explain = regions("chr1\t100\t200\n");
        assert!(explain.covers(locus(0, 101, 101)));
        assert!(explain.covers(locus(0, 200, 200)));
        assert!(
            explain.covers(locus(0, 95, 101)),
            "a deletion reaching in from before"
        );
        assert!(!explain.covers(locus(0, 95, 100)));
        assert!(!explain.covers(locus(0, 201, 210)));
        assert!(!explain.covers(locus(1, 150, 150)), "another contig");
    }

    /// **A dropped locus is noted only where it is explained**, and comes back as a `locus` row
    /// and an `outcome` row saying why.
    #[test]
    fn the_merge_s_drops_are_noted_only_where_explained() {
        let dropped = DroppedLoci::over(regions("chr1\t100\t200\n"));
        dropped.note(locus(0, 150, 150), MergeDrop::TooQuiet);
        dropped.note(locus(0, 500, 500), MergeDrop::TooQuiet);
        dropped.note(locus(0, 160, 170), MergeDrop::OverDepthCeiling);
        let rows = dropped.into_rows();
        let outcomes: Vec<&str> = rows
            .iter()
            .filter(|row| row.step == ExplainStep::Outcome)
            .map(|row| row.value.as_str())
            .collect();
        assert_eq!(outcomes, ["too_quiet", "over_depth_ceiling"]);
    }

    /// **A record the filter dropped says so in its outcome**, and every verdict lands beside
    /// the locus that wrote the record — matched by where it was written, which is a base before
    /// the locus when a padding base was put in front.
    #[test]
    fn the_paralog_verdicts_join_the_loci_that_wrote_the_records() {
        let deletion = locus(0, 150, 155);
        let snp = locus(0, 180, 180);
        let mut rows = vec![
            ExplainRow::of(deletion, ExplainStep::Outcome, "written").written_at(149),
            ExplainRow::of(snp, ExplainStep::Outcome, "written").written_at(180),
        ];
        let mut verdicts = ParalogExplanations::wanting(&rows);
        assert!(verdicts.wants(ContigId(0), 149));
        assert!(!verdicts.wants(ContigId(0), 150));
        verdicts.record(ParalogVerdictAt {
            contig: ContigId(0),
            written_at: 149,
            ratio: 30.0,
            posterior: Some(0.999),
            flagged: true,
            dropped: true,
        });
        verdicts.record(ParalogVerdictAt {
            contig: ContigId(0),
            written_at: 180,
            ratio: -5.0,
            posterior: Some(0.001),
            flagged: false,
            dropped: false,
        });
        verdicts.add_to(&mut rows);

        let of = |region: GenomeRegion, step: ExplainStep| {
            rows.iter()
                .find(|row| row.region == region && row.step == step)
                .map(|row| row.value.clone())
        };
        assert_eq!(
            of(deletion, ExplainStep::Outcome).as_deref(),
            Some("dropped_by_paralog_filter")
        );
        assert_eq!(
            of(deletion, ExplainStep::Paralog).as_deref(),
            Some("dropped")
        );
        assert_eq!(of(snp, ExplainStep::Outcome).as_deref(), Some("written"));
        assert_eq!(of(snp, ExplainStep::Paralog).as_deref(), Some("kept"));
    }

    #[test]
    fn the_steps_sort_in_the_order_a_locus_is_read() {
        assert!(ExplainStep::Locus < ExplainStep::Allele);
        assert!(ExplainStep::Site < ExplainStep::Paralog);
        assert!(ExplainStep::Paralog < ExplainStep::Outcome);
    }
}
