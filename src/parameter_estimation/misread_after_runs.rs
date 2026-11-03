//! **Counting, per read group, how much more often reads misread a base just after a long run of
//! one base when they were sequenced through it, and fitting the table the STR path charges from
//! those counts** (`doc/devel/ng/spec/str_context_error_fit.md`).
//!
//! The pre-pass counts at every generic position that has an entry of the table in one read
//! direction and none in the other ([`MisreadCounter`]): the reads of the first direction crossed a
//! run to reach the base, those of the second reached it from the other side. A true variant at
//! the base shows in both groups; a misread caused by the run only in the first. So the fit
//! ([`fit_misread_table`]) takes the difference of the two rates, and needs no truth set.
//!
//! The table itself, and which entry a base belongs to, are the scoring side's
//! ([`crate::calling::likelihood::run_context`]): one rule decides an entry for counting and for
//! charging.

use std::collections::BTreeMap;

use crate::calling::likelihood::run_context::{
    ContextAfterRun, DistanceAfterRun, ENTRY_COUNT, MisreadExcess, MisreadTable, RunEntry,
    contexts_along,
};
use crate::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::parameter_estimation::Provenance;
use crate::ref_seq::RefSeq;
use crate::types::{ContigId, ReadGroupId};

/// **The seven totals one entry keeps for one read group** (`str_context_error_fit.md` §3.3).
///
/// "Through" reads were sequenced through the run before reaching the base; "arriving first"
/// reads reached the same base from the other side. At a further-on entry the two "as the
/// run" totals stay zero and "other wrong" is any wrong base.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MisreadCounts {
    pub through_reads: u64,
    pub through_as_the_run: u64,
    pub through_other_wrong: u64,
    pub arriving_first_reads: u64,
    pub arriving_first_as_the_run: u64,
    pub arriving_first_other_wrong: u64,
    /// **What the base qualities of the reads in `through_as_the_run` claimed**: their summed
    /// natural-log error, in the steps of
    /// [`SummedLogError`](crate::types::SummedLogError). For the check of §4.2 only.
    pub through_as_the_run_log_error_steps: i64,
}

impl MisreadCounts {
    fn add(&mut self, other: &Self) {
        self.through_reads += other.through_reads;
        self.through_as_the_run += other.through_as_the_run;
        self.through_other_wrong += other.through_other_wrong;
        self.arriving_first_reads += other.arriving_first_reads;
        self.arriving_first_as_the_run += other.arriving_first_as_the_run;
        self.arriving_first_other_wrong += other.arriving_first_other_wrong;
        self.through_as_the_run_log_error_steps += other.through_as_the_run_log_error_steps;
    }
}

/// **One read group's totals at every entry of the table**, in [`RunEntry::all`]'s order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MisreadCountsByEntry {
    entries: [MisreadCounts; ENTRY_COUNT],
}

impl Default for MisreadCountsByEntry {
    fn default() -> Self {
        Self {
            entries: [MisreadCounts::default(); ENTRY_COUNT],
        }
    }
}

impl MisreadCountsByEntry {
    /// From every entry's totals, in table order.
    #[must_use]
    pub fn from_entries(entries: [MisreadCounts; ENTRY_COUNT]) -> Self {
        Self { entries }
    }

    /// Every entry's totals, in table order.
    #[must_use]
    pub fn entries(&self) -> &[MisreadCounts; ENTRY_COUNT] {
        &self.entries
    }

    /// One entry's totals.
    #[must_use]
    pub fn at(&self, entry: RunEntry) -> &MisreadCounts {
        &self.entries[entry.index()]
    }

    /// Add another set of totals for the same read group into these.
    pub fn add(&mut self, other: &Self) {
        for (mine, theirs) in self.entries.iter_mut().zip(&other.entries) {
            mine.add(theirs);
        }
    }

    fn entry_mut(&mut self, entry: RunEntry) -> &mut MisreadCounts {
        &mut self.entries[entry.index()]
    }
}

/// **How many reference bases around a position the rule needs**: the 20 bases an entry reaches
/// past a run, and enough of the run before them to put it in its longest class, 23 bases.
const CONTEXT_MARGIN: u64 = 20 + 23;

/// How many positions one cached stretch of reference covers.
const CACHED_STRETCH: u64 = 1 << 16;

/// **The reference around a stretch of positions, with each position's entry for a forward and
/// for a reverse read** — computed once a stretch, so the count is linear in the ground walked.
struct ReferenceContexts {
    contig: ContigId,
    /// 1-based position of `forward[0]` and `reverse[0]`.
    first: u64,
    forward: Vec<Option<ContextAfterRun>>,
    reverse: Vec<Option<ContextAfterRun>>,
}

impl ReferenceContexts {
    fn covers(&self, contig: ContigId, position: u64) -> bool {
        self.contig == contig
            && position >= self.first
            && position < self.first + self.forward.len() as u64
    }
}

/// **The pre-pass's count of misreads after runs, per read group**, fed every generic locus a
/// sample's walk hands the census (`str_context_error_fit.md` §3).
pub struct MisreadCounter {
    reference: Box<dyn RefSeq + Send>,
    cached: Option<ReferenceContexts>,
    counts: BTreeMap<ReadGroupId, MisreadCountsByEntry>,
    /// Scratch, cleared and refilled per stretch.
    fetched: Vec<u8>,
    walked: Vec<u8>,
    contexts: Vec<Option<ContextAfterRun>>,
}

impl std::fmt::Debug for MisreadCounter {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("MisreadCounter")
            .field("read_groups", &self.counts.len())
            .finish_non_exhaustive()
    }
}

impl MisreadCounter {
    /// A counter reading the walk's reference through `reference`.
    #[must_use]
    pub fn new(reference: Box<dyn RefSeq + Send>) -> Self {
        Self {
            reference,
            cached: None,
            counts: BTreeMap::new(),
            fetched: Vec::new(),
            walked: Vec::new(),
            contexts: Vec::new(),
        }
    }

    /// **Count one locus**: generic loci only, and at each position only the observations whose
    /// reads covered the whole locus with no insertion or deletion, compared base by base with
    /// the reference.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations) {
        if !matches!(locus.kind, LocusKind::Generic) {
            return;
        }
        let span = locus.reference_bases.len();
        let first = locus.region.start.get();
        for offset in 0..span {
            let position = first + offset as u64;
            let Some((forward, reverse)) = self.contexts_at(locus.region.contig, position) else {
                continue;
            };
            // **One direction crossed a run and the other did not**, or the position says
            // nothing: with an entry each way, neither group is a clean comparison for the other.
            let (context, through_is_forward) = match (forward, reverse) {
                (Some(context), None) => (context, true),
                (None, Some(context)) => (context, false),
                _ => continue,
            };
            let reference_base = locus.reference_bases[offset];
            let next_base = context.entry.distance() == DistanceAfterRun::NextBase;
            for observation in &locus.observations {
                if !matches!(observation.read_witness, ReadWitness::Complete)
                    || observation.bases.len() != span
                    || observation.num_obs == 0
                {
                    continue;
                }
                let forward_reads = u64::from(observation.num_fwd.min(observation.num_obs));
                let reverse_reads = u64::from(observation.num_obs) - forward_reads;
                let (through, arriving_first) = if through_is_forward {
                    (forward_reads, reverse_reads)
                } else {
                    (reverse_reads, forward_reads)
                };
                let base = observation.bases[offset];
                let as_the_run = next_base && base == context.run_base;
                let other_wrong = base != reference_base && !as_the_run;
                let counts = self
                    .counts
                    .entry(observation.read_group)
                    .or_default()
                    .entry_mut(context.entry);
                counts.through_reads += through;
                counts.arriving_first_reads += arriving_first;
                if as_the_run {
                    counts.through_as_the_run += through;
                    counts.arriving_first_as_the_run += arriving_first;
                    // The observation's summed error, its through reads' share of it: one
                    // sequence's reads, so a share by count is as good as the store allows.
                    let steps = observation.q_sum.steps() as i128 * i128::from(through)
                        / i128::from(observation.num_obs);
                    counts.through_as_the_run_log_error_steps += steps as i64;
                } else if other_wrong {
                    counts.through_other_wrong += through;
                    counts.arriving_first_other_wrong += arriving_first;
                }
            }
        }
    }

    /// Every read group's totals.
    #[must_use]
    pub fn finish(self) -> BTreeMap<ReadGroupId, MisreadCountsByEntry> {
        self.counts
    }

    /// A position's entry for a forward and for a reverse read, computed on the reference.
    /// `None` where the reference cannot be read there.
    fn contexts_at(
        &mut self,
        contig: ContigId,
        position: u64,
    ) -> Option<(Option<ContextAfterRun>, Option<ContextAfterRun>)> {
        if !self
            .cached
            .as_ref()
            .is_some_and(|cached| cached.covers(contig, position))
        {
            self.cached = self.contexts_from(contig, position);
        }
        let cached = self.cached.as_ref()?;
        let at = (position - cached.first) as usize;
        Some((cached.forward[at], cached.reverse[at]))
    }

    /// The contexts of the stretch starting at `first`, from the reference around it. A stretch
    /// near a contig's end is shorter, which is what the fetch answers with.
    fn contexts_from(&mut self, contig: ContigId, first: u64) -> Option<ReferenceContexts> {
        let fetch_from = first.saturating_sub(CONTEXT_MARGIN).max(1);
        let before = first - fetch_from;
        let wanted = before + CACHED_STRETCH + CONTEXT_MARGIN;
        // **A fetch past the contig's end is refused, so near the end the stretch is cut at the
        // end exactly**: the longest fetch that succeeds, found by halving the gap between one
        // that does and one that does not. Anything shorter than the end would leave its last
        // positions without the bases after them, and the positions past it uncached.
        let reached_the_end = self
            .reference
            .fetch_into(contig, fetch_from, wanted, &mut self.fetched)
            .is_err();
        if reached_the_end {
            let fits = |reference: &dyn RefSeq, length: u64, into: &mut Vec<u8>| {
                reference
                    .fetch_into(contig, fetch_from, length, into)
                    .is_ok()
            };
            let (mut fitting, mut too_long) = (before + 1, wanted);
            if !fits(self.reference.as_ref(), fitting, &mut self.fetched) {
                return None;
            }
            while too_long - fitting > 1 {
                let middle = fitting + (too_long - fitting) / 2;
                if fits(self.reference.as_ref(), middle, &mut self.fetched) {
                    fitting = middle;
                } else {
                    too_long = middle;
                }
            }
            self.reference
                .fetch_into(contig, fetch_from, fitting, &mut self.fetched)
                .ok()?;
        }
        let fetched_len = self.fetched.len() as u64;
        let last_covered = if reached_the_end {
            fetched_len
        } else {
            fetched_len - CONTEXT_MARGIN
        };
        let stretch = last_covered.saturating_sub(before) as usize;

        contexts_along(&self.fetched, &mut self.contexts);
        let forward = self.contexts[before as usize..][..stretch].to_vec();

        self.walked.clear();
        self.walked.extend(self.fetched.iter().rev());
        contexts_along(&self.walked, &mut self.contexts);
        let len = self.fetched.len();
        let reverse = (0..stretch)
            .map(|at| self.contexts[len - 1 - (before as usize + at)])
            .collect();
        Some(ReferenceContexts {
            contig,
            first,
            forward,
            reverse,
        })
    }
}

/// **The fewest reads in either group for an entry to be fitted** (`str_context_error_fit.md`
/// §4.1): at 1,000 reads an excess of 1 in 100 has a standard error of about 3 in 1,000.
pub const FEWEST_READS_FITTED: u64 = 1_000;

/// **One read group's fitted entry**: the excess, its standard error, the reads it rests on, and
/// whether it was fitted here or fell back to the default.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedMisreadEntry {
    pub excess: MisreadExcess,
    pub standard_error: MisreadExcess,
    pub through_reads: u64,
    pub arriving_first_reads: u64,
    pub provenance: Provenance,
}

/// How many standard errors above zero an excess must be to be charged.
pub const STANDARD_ERRORS_TO_CHARGE: f64 = 2.0;

/// **The difference of two rates, and its binomial standard error** — zero where the difference is
/// less than [`STANDARD_ERRORS_TO_CHARGE`] standard errors above zero, since such an excess cannot
/// be told from none, and charging it would cost every read near such a run a second scoring
/// for nothing the data showed.
fn excess_of(through: u64, through_reads: u64, first: u64, first_reads: u64) -> (f64, f64) {
    let p_through = through as f64 / through_reads as f64;
    let p_first = first as f64 / first_reads as f64;
    let variance = p_through * (1.0 - p_through) / through_reads as f64
        + p_first * (1.0 - p_first) / first_reads as f64;
    let difference = p_through - p_first;
    let standard_error = variance.sqrt();
    let charged = if difference >= STANDARD_ERRORS_TO_CHARGE * standard_error && difference > 0.0 {
        difference
    } else {
        0.0
    };
    (charged, standard_error)
}

/// **One read group's table from its totals** (`str_context_error_fit.md` §4.1): per entry, the
/// rate in reads through the run minus the rate in reads arriving first, never below zero; an
/// entry with fewer than [`FEWEST_READS_FITTED`] reads in either group takes `defaults`' value,
/// marked [`Provenance::Defaulted`].
#[must_use]
pub fn fit_misread_table(
    counts: &MisreadCountsByEntry,
    defaults: &MisreadTable,
) -> Vec<FittedMisreadEntry> {
    RunEntry::all()
        .map(|entry| {
            let totals = counts.at(entry);
            if totals.through_reads < FEWEST_READS_FITTED
                || totals.arriving_first_reads < FEWEST_READS_FITTED
            {
                return FittedMisreadEntry {
                    excess: defaults.at(entry),
                    standard_error: MisreadExcess::default(),
                    through_reads: totals.through_reads,
                    arriving_first_reads: totals.arriving_first_reads,
                    provenance: Provenance::Defaulted,
                };
            }
            let (as_the_run, as_the_run_error) = excess_of(
                totals.through_as_the_run,
                totals.through_reads,
                totals.arriving_first_as_the_run,
                totals.arriving_first_reads,
            );
            let (other_wrong, other_wrong_error) = excess_of(
                totals.through_other_wrong,
                totals.through_reads,
                totals.arriving_first_other_wrong,
                totals.arriving_first_reads,
            );
            FittedMisreadEntry {
                excess: MisreadExcess {
                    as_the_run,
                    other_wrong,
                },
                standard_error: MisreadExcess {
                    as_the_run: as_the_run_error,
                    other_wrong: other_wrong_error,
                },
                through_reads: totals.through_reads,
                arriving_first_reads: totals.arriving_first_reads,
                provenance: Provenance::FittedHere,
            }
        })
        .collect()
}

/// The table a fit's entries make, in table order.
#[must_use]
pub fn table_of(fitted: &[FittedMisreadEntry]) -> MisreadTable {
    let mut entries = [MisreadExcess::default(); ENTRY_COUNT];
    for (slot, entry) in entries.iter_mut().zip(fitted) {
        *slot = entry.excess;
    }
    MisreadTable::from_entries(entries)
}

/// The fewest misreads at an entry for the base-quality check to speak (§4.2).
pub const FEWEST_MISREADS_CHECKED: u64 = 100;

/// How far the base qualities' claim and the fitted rate may differ before the check warns.
pub const QUALITY_CHECK_FACTOR: f64 = 3.0;

/// **The check of the fit against what the base qualities claimed**, at one next-base entry
/// (`str_context_error_fit.md` §4.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualityCheck {
    pub entry: RunEntry,
    /// The misread reads' base qualities' own error rate: `exp(mean ln error)`.
    pub rate_the_qualities_claim: f64,
    /// **The rate the reads through the run show for reading the base as the run** — the
    /// background both directions share plus the excess. Not the excess alone: the misread bases
    /// include the background ones, and at a short run, where the excess is near zero, the
    /// qualities would then always look worse than the fit.
    pub rate_the_through_reads_show: f64,
    pub misreads: u64,
}

impl QualityCheck {
    /// Whether the two rates differ by more than [`QUALITY_CHECK_FACTOR`].
    #[must_use]
    pub fn disagrees(&self) -> bool {
        let ratio = self.rate_the_qualities_claim / self.rate_the_through_reads_show;
        !(1.0 / QUALITY_CHECK_FACTOR..=QUALITY_CHECK_FACTOR).contains(&ratio)
    }
}

/// **Every fitted next-base entry with enough misreads, checked against its base qualities.**
/// Never changes a value: it exists to make a broken count visible.
#[must_use]
pub fn quality_checks(
    counts: &MisreadCountsByEntry,
    fitted: &[FittedMisreadEntry],
) -> Vec<QualityCheck> {
    RunEntry::all()
        .zip(fitted)
        .filter(|(entry, fit)| {
            entry.distance() == DistanceAfterRun::NextBase
                && fit.provenance == Provenance::FittedHere
                && counts.at(*entry).through_as_the_run >= FEWEST_MISREADS_CHECKED
        })
        .map(|(entry, _)| {
            let totals = counts.at(entry);
            let mean_ln = totals.through_as_the_run_log_error_steps as f64
                / crate::types::SummedLogError::STEPS_PER_NAT as f64
                / totals.through_as_the_run as f64;
            QualityCheck {
                entry,
                rate_the_qualities_claim: crate::float::exp(mean_ln),
                rate_the_through_reads_show: totals.through_as_the_run as f64
                    / totals.through_reads as f64,
                misreads: totals.through_as_the_run,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locus_generation::SequenceObservation;
    use crate::ref_seq::RefSeqError;
    use crate::types::{GenomeRegion, Position, SummedLogError};

    /// One contig's bases, served as the walk's reference is.
    struct OneContig(Vec<u8>);

    impl RefSeq for OneContig {
        fn fetch_into(
            &self,
            _contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            let from = (start_1based - 1) as usize;
            let to = from + length as usize;
            if to > self.0.len() {
                return Err(RefSeqError::OutOfBounds {
                    contig: ContigId(0),
                    contig_length: self.0.len() as u64,
                    start: start_1based,
                    end: start_1based + length,
                });
            }
            dst.clear();
            dst.extend_from_slice(&self.0[from..to]);
            Ok(())
        }
    }

    fn observation(base: u8, reads: u32, forward: u32, group: u32) -> SequenceObservation {
        SequenceObservation {
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            held_flank: None,
            bases: [base][..].into(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(group),
            num_obs: reads,
            num_fwd: forward,
            q_sum: SummedLogError::from_nats(f64::from(reads) * crate::float::ln(0.1)),
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            num_low_quality: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus_at(
        reference: &[u8],
        position: u64,
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        let base = reference[(position - 1) as usize];
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(position),
                end: Position(position),
            },
            reference_bases: [base][..].into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// GGCGTC, then `A`×12, then `C` at position 19, then non-repetitive sequence.
    fn genome() -> Vec<u8> {
        [
            &b"GGCGTC"[..],
            &[b'A'; 12],
            b"CTGACGTTGCAGTCAGCTAGCTGACTGCA",
        ]
        .concat()
    }

    /// **The forward reads at the base after `A`×12 count as through reads**, the reverse ones as
    /// arriving first; a forward read showing `A` is a misread as one more base of the run.
    #[test]
    fn the_base_after_a_run_counts_forward_reads_as_through() {
        let reference = genome();
        let mut counter = MisreadCounter::new(Box::new(OneContig(reference.clone())));
        counter.add_locus(&locus_at(
            &reference,
            19,
            vec![
                observation(b'C', 10, 4, 0),
                observation(b'A', 3, 3, 0),
                observation(b'G', 1, 0, 0),
            ],
        ));
        let counts = counter.finish();
        let entry = RunEntry::of(12, 1).unwrap();
        let totals = counts[&ReadGroupId(0)].at(entry);
        assert_eq!(totals.through_reads, 4 + 3);
        assert_eq!(totals.through_as_the_run, 3);
        assert_eq!(totals.arriving_first_reads, 6 + 1);
        assert_eq!(totals.arriving_first_other_wrong, 1);
        assert_eq!(totals.arriving_first_as_the_run, 0);
        // Three reads at quality 10: three times ln(0.1).
        let claimed =
            totals.through_as_the_run_log_error_steps as f64 / SummedLogError::STEPS_PER_NAT as f64;
        assert!(
            (claimed - 3.0 * crate::float::ln(0.1)).abs() < 1e-3,
            "{claimed}"
        );
    }

    /// A position with a run behind it in both directions is not counted, and neither is an
    /// observation of the wrong length or a position far from any run.
    #[test]
    fn positions_without_a_clean_comparison_are_not_counted() {
        // `A`×12 `G` `T`×12: the `G` follows a charged run from either side.
        let reference = [&b"GCGC"[..], &[b'A'; 12], b"G", &[b'T'; 12], b"CGCA"].concat();
        let mut counter = MisreadCounter::new(Box::new(OneContig(reference.clone())));
        counter.add_locus(&locus_at(&reference, 17, vec![observation(b'A', 5, 2, 0)]));
        // Far from any run.
        let plain = b"GACGTCAGCTAGGTCACGTAGCATGCAGTC".to_vec();
        let mut far = MisreadCounter::new(Box::new(OneContig(plain.clone())));
        far.add_locus(&locus_at(&plain, 10, vec![observation(b'G', 5, 2, 0)]));
        // An observation with an insertion: the wrong length for the locus.
        let reference2 = genome();
        far.reference = Box::new(OneContig(reference2.clone()));
        far.cached = None;
        let mut wrong_length = observation(b'C', 5, 2, 0);
        wrong_length.bases = b"CA"[..].into();
        far.add_locus(&locus_at(&reference2, 19, vec![wrong_length]));
        for counts in [counter.finish(), far.finish()] {
            assert!(
                counts
                    .values()
                    .all(|by_entry| by_entry.entries().iter().all(|t| t.through_reads == 0)),
                "{counts:?}"
            );
        }
    }

    /// The reverse direction mirrors: a run to the right of a base makes its reverse reads the
    /// through reads.
    #[test]
    fn a_run_to_the_right_makes_reverse_reads_the_through_reads() {
        let reference = [&b"GACGTCAGCTAGGTC"[..], &[b'T'; 12], b"GCAGCTA"].concat();
        let mut counter = MisreadCounter::new(Box::new(OneContig(reference.clone())));
        // Position 15 is the `C` just left of the twelve `T`s.
        counter.add_locus(&locus_at(&reference, 15, vec![observation(b'T', 4, 1, 0)]));
        let counts = counter.finish();
        let totals = counts[&ReadGroupId(0)].at(RunEntry::of(12, 1).unwrap());
        assert_eq!(totals.through_reads, 3);
        assert_eq!(totals.through_as_the_run, 3);
        assert_eq!(totals.arriving_first_reads, 1);
        assert_eq!(totals.arriving_first_as_the_run, 1);
    }

    fn counts_at(entry: RunEntry, totals: MisreadCounts) -> MisreadCountsByEntry {
        let mut entries = [MisreadCounts::default(); ENTRY_COUNT];
        entries[entry.index()] = totals;
        MisreadCountsByEntry::from_entries(entries)
    }

    #[test]
    fn the_fit_is_the_difference_never_below_zero_and_sparse_entries_default() {
        let entry = RunEntry::of(16, 1).unwrap();
        let counts = counts_at(
            entry,
            MisreadCounts {
                through_reads: 2_000,
                through_as_the_run: 200,
                through_other_wrong: 4,
                arriving_first_reads: 4_000,
                arriving_first_as_the_run: 40,
                arriving_first_other_wrong: 40,
                through_as_the_run_log_error_steps: 0,
            },
        );
        let defaults = MisreadTable::defaults();
        let fitted = fit_misread_table(&counts, &defaults);
        let this = fitted[entry.index()];
        assert_eq!(this.provenance, Provenance::FittedHere);
        assert!((this.excess.as_the_run - (0.1 - 0.01)).abs() < 1e-12);
        // Through 2 in 1,000 against 1 in 100 arriving first: below zero, so zero.
        assert_eq!(this.excess.other_wrong, 0.0);
        assert!(this.standard_error.as_the_run > 0.006 && this.standard_error.as_the_run < 0.007);
        // Every other entry saw nothing and takes the default.
        let other = RunEntry::of(12, 1).unwrap();
        assert_eq!(fitted[other.index()].provenance, Provenance::Defaulted);
        assert_eq!(fitted[other.index()].excess, defaults.at(other));
        assert_eq!(table_of(&fitted).at(entry), this.excess);
    }

    /// **An excess within two standard errors of zero is charged as zero**: 1.1 in 100 through
    /// against 1.0 in 100 arriving first, over a few thousand reads, is a difference the counts
    /// cannot tell from none.
    #[test]
    fn an_excess_the_counts_cannot_tell_from_none_is_zero() {
        let entry = RunEntry::of(20, 1).unwrap();
        let counts = counts_at(
            entry,
            MisreadCounts {
                through_reads: 2_000,
                through_as_the_run: 22,
                through_other_wrong: 0,
                arriving_first_reads: 4_000,
                arriving_first_as_the_run: 40,
                arriving_first_other_wrong: 0,
                through_as_the_run_log_error_steps: 0,
            },
        );
        let fitted = fit_misread_table(&counts, &MisreadTable::defaults())[entry.index()];
        assert_eq!(fitted.provenance, Provenance::FittedHere);
        assert_eq!(fitted.excess.as_the_run, 0.0);
        assert!(fitted.standard_error.as_the_run > 0.002);
    }

    #[test]
    fn misreads_claiming_far_better_quality_than_the_fit_are_flagged() {
        let entry = RunEntry::of(16, 1).unwrap();
        let misreads = 200_u64;
        let at = |quality: f64| {
            let ln_error = crate::float::ln(crate::float::powf(10.0, -quality / 10.0));
            counts_at(
                entry,
                MisreadCounts {
                    through_reads: 2_000,
                    through_as_the_run: misreads,
                    through_other_wrong: 0,
                    arriving_first_reads: 2_000,
                    arriving_first_as_the_run: 0,
                    arriving_first_other_wrong: 0,
                    through_as_the_run_log_error_steps: (ln_error
                        * misreads as f64
                        * SummedLogError::STEPS_PER_NAT as f64)
                        as i64,
                },
            )
        };
        let defaults = MisreadTable::defaults();
        for (quality, flagged) in [(10.0, false), (30.0, true)] {
            let counts = at(quality);
            let fitted = fit_misread_table(&counts, &defaults);
            let checks = quality_checks(&counts, &fitted);
            assert_eq!(checks.len(), 1);
            assert_eq!(
                checks[0].disagrees(),
                flagged,
                "quality {quality}: {:?}",
                checks[0]
            );
        }
    }
}
