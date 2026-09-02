//! **The repeat-tract path: narrowing one tract's allele table to the tract sequences worth
//! calling over** (`doc/devel/ng/spec/candidate_alleles_ssr.md`, and
//! `doc/devel/ng/arch/candidate_alleles_ssr.md` for the shapes).
//!
//! **Why this is a sibling of [`generic`](super::generic) and not a configuration of it.** At an
//! ordinary locus the merge's alleles are an unordered set and selection is a bar and a cap over
//! it. At a repeat tract they are not unordered: a tract of 11 repeats is adjacent to one of 10
//! and one of 12 and far from one of 4, and both the genotype prior
//! (`doc/devel/ng/spec/calling_priors.md` §5) and the stutter model
//! (`doc/devel/ng/spec/read_likelihoods.md` §4.2) are written on that ordering. So this path
//! groups the table's sequences by how many repeats they carry, lets each sample nominate which
//! groups are worth promoting, and asks the shared support rule of the sequences on a promoted
//! group. The bar, the cap, the leftover and the remapping are the parent module's and are
//! called, not rewritten.
//!
//! **Under construction, in the order `doc/devel/ng/impl_plan/candidate_alleles_ssr.md` sets.**
//! This file holds the ladder (its step B1), the path's configuration and the per-sample length
//! histogram (B2). Nomination, admission, periodicity and the entry point `select_ssr` are the
//! steps after those, and nothing outside this module calls anything here yet.
//!
//! **Which is why a non-test build expects everything here to be dead.** The expectation rather
//! than an `allow`, so that the first real caller — `select_ssr`, which the steps after this one
//! build — turns the line below into a compile error and deletes it. Under `cfg(test)` nothing is
//! exempted: the tests in this file are what exercise the ladder today, and they must keep
//! failing when it is wrong.
#![cfg_attr(not(test), expect(dead_code))]

use super::{CandidateSelectionConfig, MaxCandidateAlleles, SelectionScratch};
use crate::ng::run::cohort_merge::MinAltReadShare;
use crate::ng::run::cohort_merge::build::{CohortObservation, SampleSupport};
use crate::ng::types::{GenomeRegion, Motif, Ploidy};

/// **One locus's observed tract sequences, grouped by how many repeats they carry.**
///
/// The key is `bases.len() / motif.period()`, floored — **the same integer the genotype prior
/// indexes its length spectrum by** (`doc/devel/ng/spec/candidate_alleles_ssr.md` §3), which is
/// the reason the grouping exists at all and the reason both must compute it one way. Production
/// keys its own ladder exactly so (`src/ssr/cohort/rung_ladder.rs:291-320`) and ng keeps that.
///
/// **Sequences inside a rung stay distinct.** The grouping is for nomination and for the prior,
/// never a merge of evidence: two tract sequences of one length differing by an interior base are
/// a real allele class — an interrupted repeat — and the read likelihood separates them by about
/// 28 Phred per distinguishing base (spec §1.1, and `read_likelihoods.md` §4.6). A rung therefore
/// holds *indices into the merge's table*, and owns no bases.
///
/// # It is a scratch buffer, so it means nothing until it is filled
///
/// One of these lives in [`SelectionScratch`] and is refilled at every locus, so a locus costs no
/// allocation for the grouping. Between [`SelectionScratch::reset_for`] and the next
/// [`build_ladder`] it is empty, and after a build it describes **that** locus until the next
/// reset. Every accessor below panics on an empty ladder rather than answering with a number
/// from nowhere.
///
/// # The shape, and where it departs from the architecture sketch
///
/// Arch §2.1 sketches `rungs: Vec<Rung>` with each `Rung` owning a `Vec<u32>` of table indices.
/// That is one heap allocation per rung per locus, which is exactly what a scratch buffer exists
/// to avoid, so the indices live in **one** buffer ordered by `(repeat count, table index)` and a
/// rung names its slice of it. Same information, same order, no per-locus allocation.
#[derive(Default, Debug)]
pub(super) struct RepeatLadder {
    /// Every allele of the merge's table exactly once, ordered by `(repeat count, table index)`.
    /// Each rung's [`Rung::table_indices_from`] and [`Rung::table_indices_len`] name a slice of
    /// it, so within a rung the indices are ascending.
    table_indices_by_rung: Vec<u32>,
    /// The occupied rungs, **ascending by repeat count**. Empty exactly when the ladder has not
    /// been built for a locus yet.
    rungs: Vec<Rung>,
    /// **The other direction: for each allele of the merge's table, in that table's index order,
    /// which rung it sits on.**
    ///
    /// [`table_indices_by_rung`](Self::table_indices_by_rung) answers *which sequences are on this
    /// rung*, which is what admission walks; this answers *which rung is this sequence on*, which
    /// is what a walk over one sample's own support rows needs — and those rows name merge-table
    /// indices, in the merge's order, not the ladder's. Without it a per-sample fold would have to
    /// re-derive each sequence's repeat count from its bases and search the rungs for it, which is
    /// a second producer of the one integer this whole type exists to have one producer of.
    ///
    /// One `u32` per allele, refilled per locus like the rest.
    rung_of_table_index: Vec<u32>,
    /// The cohort's modal repeat count — the most-supported rung, ties to the shorter. Filled by
    /// [`build_ladder`], and meaningless while `rungs` is empty, which is why
    /// [`modal_repeat_count`](Self::modal_repeat_count) asserts rather than reads.
    modal_repeat_count: u32,
}

/// **One repeat count and the merge-table indices of the sequences observed at it.**
///
/// `Copy`, and it names a slice of the ladder's one index buffer rather than owning a `Vec` —
/// see [`RepeatLadder`]'s note on the architecture sketch.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Rung {
    /// How many whole motif units the sequences on this rung carry.
    repeat_count: u32,
    /// Where this rung's indices start in [`RepeatLadder::table_indices_by_rung`].
    table_indices_from: u32,
    /// How many indices this rung holds. At least one — a rung with no sequences is never
    /// created.
    table_indices_len: u32,
    /// This rung's reads over the whole cohort: the sum of its sequences' own cohort totals from
    /// the fold ([`AlleleSummary`](super::AlleleSummary)), which is **every covering sample's
    /// reads on them**, including samples that did not clear the support rule.
    ///
    /// Saturating, like every other read total in this module: it needs 18 quintillion reads on
    /// one rung of one locus to bind.
    cohort_reads: u64,
}

impl RepeatLadder {
    /// Empty the ladder without releasing what earlier loci reserved.
    ///
    /// **`clear` on all three buffers, and the mode reset with them.** A ladder left holding the
    /// previous locus's mode while its rungs are empty would answer
    /// [`modal_repeat_count`](Self::modal_repeat_count) with a neighbour's number instead of
    /// asserting, which is the one failure the assertion exists to prevent.
    #[inline]
    pub(super) fn clear(&mut self) {
        let Self {
            table_indices_by_rung,
            rungs,
            rung_of_table_index,
            modal_repeat_count,
        } = self;
        table_indices_by_rung.clear();
        rungs.clear();
        rung_of_table_index.clear();
        *modal_repeat_count = 0;
    }

    /// How many distinct repeat counts the locus's sequences occupy. Zero only before a build.
    #[inline]
    pub(super) fn rung_count(&self) -> usize {
        self.rungs.len()
    }

    /// The repeat count of the `rung`-th occupied rung, counting from the shortest.
    ///
    /// # Panics
    ///
    /// On a rung this ladder does not hold, which is a caller walking a buffer built for another
    /// locus.
    #[inline]
    pub(super) fn repeat_count_at(&self, rung: usize) -> u32 {
        self.rung_at(rung).repeat_count
    }

    /// The merge-table indices of the sequences on the `rung`-th rung, **ascending**.
    ///
    /// # Panics
    ///
    /// On a rung this ladder does not hold.
    #[inline]
    pub(super) fn table_indices_at(&self, rung: usize) -> &[u32] {
        let rung = self.rung_at(rung);
        let from = rung.table_indices_from as usize;
        &self.table_indices_by_rung[from..from + rung.table_indices_len as usize]
    }

    /// The `rung`-th rung's reads over every covering sample — see [`Rung::cohort_reads`].
    ///
    /// # Panics
    ///
    /// On a rung this ladder does not hold.
    #[inline]
    pub(super) fn cohort_reads_at(&self, rung: usize) -> u64 {
        self.rung_at(rung).cohort_reads
    }

    /// **Which rung a repeat count sits on, or `None` where the merge's table holds no sequence
    /// of that length.** A binary search, because the rungs are ascending by repeat count.
    ///
    /// **This is not quite production's `occupied` test and the difference is one rung.** That
    /// test asks whether the cohort's reads reached a length — `cohort_support(length) > 0`
    /// (`src/ssr/cohort/candidate_set.rs:221`) — and the `±1` rescue asks it so that nothing
    /// invents a length (spec §4). Every rung here holds a sequence some sample's reads showed,
    /// with **one exception**: the merge interns the reference tract at index 0 whether or not a
    /// read landed on it, so the reference's rung can exist carrying zero reads. A caller that
    /// wants production's question therefore asks this one **and**
    /// [`cohort_reads_at`](Self::cohort_reads_at) — which is the rescue's business, in the step
    /// that builds it, not this accessor's.
    ///
    /// # Panics
    ///
    /// On an unbuilt ladder — the question asked of no locus would answer `None` for every
    /// length, which reads as "the table holds none of them".
    #[inline]
    pub(super) fn rung_of_repeat_count(&self, repeat_count: u32) -> Option<usize> {
        self.assert_built();
        self.rungs
            .binary_search_by_key(&repeat_count, |rung| rung.repeat_count)
            .ok()
    }

    /// **Which rung one of the merge table's sequences sits on** — the direction a walk over a
    /// sample's own support rows needs, since those rows name merge-table indices.
    ///
    /// # Panics
    ///
    /// On a `table_index` outside the merge table this ladder was built for — a support row
    /// naming an allele the locus does not hold, which is the merge bug spec §8 makes an assertion
    /// rather than a value.
    #[inline]
    pub(super) fn rung_of_table_index(&self, table_index: usize) -> usize {
        assert!(
            table_index < self.rung_of_table_index.len(),
            "a support row named allele {table_index} of a repeat tract whose table holds {}",
            self.rung_of_table_index.len()
        );
        self.rung_of_table_index[table_index] as usize
    }

    /// **The cohort's modal repeat count: the rung with the most reads, ties to the shorter.**
    ///
    /// **Not the reference tract's count** (arch §2.1). At a tract where the panel carries a
    /// length the reference does not, the reference's own rung can hold a handful of reads and
    /// the mode is elsewhere — which is the case this quantity exists for.
    ///
    /// **What reads it, and what no longer does.** The periodicity test measures each read's
    /// offset from this (spec §7). The genotype prior does *not*: its seed was re-indexed on
    /// 2026-08-27 by offset from the **reference** tract length, which every locus already knows
    /// (`doc/devel/ng/spec/population_diversity.md` §4.2), and handing it the mode instead moves
    /// 0.595 of its mass off the reference length onto 0.091 on its own fixture.
    ///
    /// # Panics
    ///
    /// On an unbuilt ladder. The alternative is answering 0, which is a legal repeat count — a
    /// tract a deletion removed entirely — so a caller could not tell the two apart.
    #[inline]
    pub(super) fn modal_repeat_count(&self) -> u32 {
        self.assert_built();
        self.modal_repeat_count
    }

    /// The `rung`-th rung, with the message a caller reading another locus's buffer needs.
    #[inline]
    fn rung_at(&self, rung: usize) -> Rung {
        assert!(
            rung < self.rungs.len(),
            "rung {rung} of a ladder holding {}: the ladder is a reused buffer, so this is a \
             reader walking one locus's rungs over another locus's ladder",
            self.rungs.len()
        );
        self.rungs[rung]
    }

    /// Refuse a question about a locus this ladder was never built for.
    #[inline]
    fn assert_built(&self) {
        assert!(
            !self.rungs.is_empty(),
            "the ladder holds no rung, so it has not been built for a locus: every merge table \
             holds at least its reference allele, so a built ladder has at least one rung"
        );
    }
}

/// **Group one repeat tract's alleles by repeat count, and find the cohort's modal count**
/// (spec §3; arch §2.1) — the ladder every later step of this path reads.
///
/// **Run it after [`summarise_alleles`](super::summarise_alleles) on the same locus**, which is
/// what fills the per-allele read totals a rung's own total is summed from and what sizes the
/// scratch. Doing the grouping here rather than in a walk of its own is the reuse the parent
/// module's shape asks for: the fold already visits every sample's rows once, and a rung's reads
/// are a sum over its sequences' totals.
///
/// **The key is floor division and nothing else.** A sequence whose length is not a whole number
/// of motif units lands on the rung below and is counted there once — it is not dropped and not
/// counted twice. Whether such a sample is periodic at all is a separate question asked later and
/// per sample (spec §7); the ladder's job is to place every observed sequence somewhere.
///
/// **The mode is the most-supported rung and ties break toward the shorter**, which is the tie
/// rule production uses and is kept for determinism rather than for a measured reason: the output
/// must be byte-identical at any worker count (`doc/devel/ng/spec/candidate_alleles.md` §8), so
/// no tie may fall through to the order the samples were walked in.
///
/// # Panics
///
/// - On an empty allele table. Every merge locus interns its reference at index 0
///   ([`CohortObservation::alleles`]), so this is a caller bug rather than an input, and without
///   the assertion it surfaces as an unbuilt ladder several steps later.
/// - On a scratch whose fold is sized for a different allele table. Summing another locus's read
///   totals onto this locus's rungs is a wrong mode and a wrong nomination, and nothing
///   downstream would see it. **The check is the buffer's width and nothing more**, so it catches
///   a fold of a different-sized locus and cannot catch a fold of a same-sized one — what makes
///   that unreachable is that the only caller folds and then builds, in that order, which is why
///   the order is stated above rather than left to the assertion.
/// - On a ladder that already holds rungs — this tract reached the build without being folded.
///   Appending to it would put two loci's sequences on one ladder.
/// - **On a merge table of more than four billion alleles**, where the `u32` conversions refuse
///   rather than wrap. Not reachable through any input this project can build; stated because the
///   conversions are visible in the code.
pub(super) fn build_ladder(
    observation: &CohortObservation,
    motif: &Motif,
    scratch: &mut SelectionScratch,
) {
    let table_len = observation.alleles.len();
    assert!(
        table_len > 0,
        "a cohort locus always holds at least its reference allele, and the repeat tract at {} \
         holds none",
        observation.region
    );
    assert_eq!(
        scratch.per_allele.len(),
        table_len,
        "the ladder sums each rung's reads from the fold's per-allele totals, so the scratch's \
         fold must be this locus's — it is sized for a table of {} where this tract at {} holds \
         {table_len}",
        scratch.per_allele.len(),
        observation.region
    );

    // Destructured rather than field-accessed, for the reason `reset_for` gives: a buffer added
    // to the scratch later has to be answered for at every site that takes one apart, instead of
    // being silently carried from one locus into the next.
    let SelectionScratch {
        per_allele,
        ranked_table_indices: _,
        ladder,
        sample_reads_per_rung: _,
    } = scratch;
    // **Asserted empty rather than emptied here.** `reset_for` already clears the ladder, and the
    // fold calls it, so a `clear()` on this line would be a second owner of one rule — and being
    // the redundant one, it can only hide the case where the fold did not run. Measured: deleting
    // that `clear()` broke no test, which is what a line with no failing state looks like.
    assert!(
        ladder.rung_count() == 0
            && ladder.table_indices_by_rung.is_empty()
            && ladder.rung_of_table_index.is_empty(),
        "the ladder already holds {} rung(s) over {} index/indices, with {} allele(s) placed, at \
         {}: it is emptied by the fold's own `SelectionScratch::reset_for`, so a full one here \
         means this tract was never folded — and appending to it would put two loci's sequences \
         on one ladder. **All three buffers are named**, because a `clear` that emptied one and \
         not the others would leave the rungs naming slices of another locus's indices",
        ladder.rung_count(),
        ladder.table_indices_by_rung.len(),
        ladder.rung_of_table_index.len(),
        observation.region
    );

    // The period is at least one base — `Motif::new` refuses an empty unit — so the division
    // below cannot be by zero.
    let period = motif.period();
    let repeat_count_of =
        |table_index: u32| (observation.alleles[table_index as usize].len() / period) as u32;

    ladder.table_indices_by_rung.extend(
        (0..table_len)
            .map(|index| u32::try_from(index).expect("a merge table narrower than four billion")),
    );
    // `(repeat count, table index)` rather than the repeat count alone: an unstable sort may
    // reorder equal keys, and the merge table's own order is what a rung's indices must keep so
    // that two runs of the same locus admit its sequences in the same order.
    ladder
        .table_indices_by_rung
        .sort_unstable_by_key(|&index| (repeat_count_of(index), index));

    // One entry per allele, overwritten below — every allele lands on exactly one rung, so no
    // entry survives as the placeholder.
    ladder.rung_of_table_index.resize(table_len, u32::MAX);
    for (position, &table_index) in ladder.table_indices_by_rung.iter().enumerate() {
        let repeat_count = repeat_count_of(table_index);
        let reads = per_allele[table_index as usize].cohort_reads;
        match ladder.rungs.last_mut() {
            Some(rung) if rung.repeat_count == repeat_count => {
                rung.table_indices_len += 1;
                rung.cohort_reads = rung.cohort_reads.saturating_add(reads);
            }
            _ => ladder.rungs.push(Rung {
                repeat_count,
                table_indices_from: u32::try_from(position)
                    .expect("a merge table narrower than four billion"),
                table_indices_len: 1,
                cohort_reads: reads,
            }),
        }
        ladder.rung_of_table_index[table_index as usize] = u32::try_from(ladder.rungs.len() - 1)
            .expect("a merge table narrower than four billion");
    }

    // Strictly greater, walking the rungs shortest-first, so a tie keeps the shorter rung — the
    // tie rule this function's own documentation states.
    let mut modal = ladder.rungs[0];
    for &rung in &ladder.rungs[1..] {
        if rung.cohort_reads > modal.cohort_reads {
            modal = rung;
        }
    }
    ladder.modal_repeat_count = modal.repeat_count;
}

/// **The share of one sample's spanning reads that may sit off the motif's grid before that
/// sample is judged non-periodic** (spec §7).
///
/// A validated fraction of one rather than the bare `f64` arch §2.2 sketches, for the reason
/// [`MinAltReadShare`] gives about its own range and for one this module has already been bitten
/// by: a mistyped share here does not crash, it deletes the gate. **Above one, no sample is ever
/// non-periodic and the verdict can never be reached; below zero, every sample is, and every
/// repeat tract in the run comes back as the reference alone.** Both are runs whose output looks
/// entirely ordinary.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct MaxOffGridShare(f64);

impl MaxOffGridShare {
    /// The default share, [`DEFAULT_MAX_OFF_GRID_SHARE`].
    pub const DEFAULT: Self = DEFAULT_MAX_OFF_GRID_SHARE;

    /// The share, or `None` if it is not a fraction of one — negative, above one, or not a
    /// number. Refusing rather than clamping, for the reason on the type.
    pub fn new(share: f64) -> Option<Self> {
        MinAltReadShare::new(share).map(|share| Self(share.get()))
    }

    /// The same share, for a `const` that has to name one — and it **panics** where
    /// [`new`](Self::new) returns `None`, exactly as [`MinAltReadShare::new_or_panic`] does and
    /// with the same caveat: called at runtime it aborts, so a share an operator typed or a file
    /// carried goes through [`new`](Self::new).
    pub const fn new_or_panic(share: f64) -> Self {
        Self(MinAltReadShare::new_or_panic(share).get())
    }

    /// The share, as a fraction of 1.
    #[inline]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl Default for MaxOffGridShare {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// **One read in ten off the motif grid** — production's `max_out_of_frame_frac`
/// ([`candidate_set.rs:85`](../../../../src/ssr/cohort/candidate_set.rs)), **inherited and never
/// measured, by them or by us** (spec §12, Q1).
///
/// **Soft in the strongest sense available**: every number spec §4.1 and §5 report was taken with
/// this gate switched off entirely, so nothing measured constrains it, and nothing in the shape of
/// the code depends on the value.
pub const DEFAULT_MAX_OFF_GRID_SHARE: MaxOffGridShare = MaxOffGridShare(0.10);

/// **32 tract sequences including the reference, where the ordinary path allows six**
/// ([`DEFAULT_MAX_CANDIDATE_ALLELES`](super::DEFAULT_MAX_CANDIDATE_ALLELES)) — the owner's
/// decision of 2026-08-24 on spec §12's Q2.
///
/// **A repeat tract carries more real alleles than a SNP does**, six was inherited from a
/// SNP/indel setting with nothing behind it for tracts, and HipSTR — the nearest comparator — has
/// no allele limit at all: it admits every sequence clearing a per-sample test and abandons the
/// locus only if the haplotype product exceeds 1,000.
///
/// **Soft, and still unmeasured**: 32 is a judgement bounded by cost. A locus over `A` alleles has
/// `A(A+1)/2` diploid genotypes and the loop scores every sample against every one — **528 here
/// against six's 21**, and 2,080 at 64. What would settle it is the tomato panel's tracts through
/// the merge, histogrammed at 1, 4, 16 and 63 accessions.
pub const DEFAULT_MAX_CANDIDATE_ALLELES_SSR: MaxCandidateAlleles =
    MaxCandidateAlleles::new_or_panic(32);

/// **The repeat-tract path's settings** — the shared support rule and cap, the periodicity gate's
/// share, and the copy number a sample is nominated at (arch §2.2).
///
/// **No `Default`, and that is the design.** [`ploidy`](Self::ploidy) is the run's, from
/// `FrozenParameters`, and a constant here would be a diploid assumption written where a polyploid
/// crop is in scope. [`at_ploidy`](Self::at_ploidy) is the only constructor of the defaults, so
/// there is no way to reach one without naming a ploidy.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SsrSelectionConfig {
    /// The support rule and the cap — [`CandidateSelectionConfig`], reused rather than restated,
    /// so that a sweep of the bar is a sweep of both paths.
    ///
    /// **Only the cap differs by default**, at [`DEFAULT_MAX_CANDIDATE_ALLELES_SSR`]. The support
    /// rule is the ordinary path's, unchanged — see [`at_ploidy`](Self::at_ploidy) for the one
    /// place that departs from what the spec's text says.
    pub shared: CandidateSelectionConfig,
    /// How much of one sample's spanning reads may sit off the motif grid before that sample is
    /// judged non-periodic (spec §7). Default [`DEFAULT_MAX_OFF_GRID_SHARE`].
    pub max_off_grid_share: MaxOffGridShare,
    /// Copies per genome — **how many rungs a sample promotes** (arch §3.1). From the run's
    /// parameters, never a constant in this module: production hard-asserts diploid at its own
    /// genotyping step and ng does not, and a polyploid run changes only this number.
    pub ploidy: Ploidy,
}

impl SsrSelectionConfig {
    /// The defaults, at the run's ploidy.
    ///
    /// # The support share is the ordinary path's 10 in 100, where spec §5 writes 5
    ///
    /// **The spec's own reason for its number is what points at this one.** It sets the tract
    /// share to 5 in 100 *"so that one number governs both paths"* (§5), citing the ordinary
    /// path's value at the time it was written. That value moved to **10 in 100** on the same day,
    /// by the owner's decision, taken against what recall alone would say and for a cost nothing
    /// had yet measured — the candidate count carried in every genotype table, at every sample,
    /// for the life of the locus (see
    /// [`DEFAULT_MIN_ALLELE_SUPPORT`](super::DEFAULT_MIN_ALLELE_SUPPORT)). Setting 5 here would
    /// create the second number the spec's sentence exists to avoid.
    ///
    /// **What it costs, from the spec's own sweep at 300× on HG002** (§5): on the class this path
    /// was designed for — a heterozygote whose two copies are the same length spelled differently,
    /// 296 of HG002's 695 heterozygous tracts — 5 in 100 offers both spellings at **86.1%** of
    /// them at 1.26 candidate sequences per tract, and 10 in 100 at **85.8%** at 1.22. **Three
    /// tracts in a thousand, and fewer candidates.** At 30× and below the two rules are the same
    /// rule, because the floor of two reads decides.
    ///
    /// **This shifts the numbers Milestone E of the plan checks against** — 85.8% and 1.22 where
    /// the plan quotes 86.1% and 1.26 — and it is recorded so that the difference is read as this
    /// decision rather than as the defect that plan asks to be traced.
    pub fn at_ploidy(ploidy: Ploidy) -> Self {
        Self {
            shared: CandidateSelectionConfig {
                max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES_SSR,
                ..CandidateSelectionConfig::DEFAULT
            },
            max_off_grid_share: MaxOffGridShare::DEFAULT,
            ploidy,
        }
    }
}

/// **Count one sample's spanning reads onto the ladder's rungs** — the per-sample length histogram
/// nomination reads (arch §3.1; production's `sample_histogram`,
/// [`rung_ladder.rs:262`](../../../../src/ssr/cohort/rung_ladder.rs), over the merge's rows
/// instead of its own sequence counts).
///
/// `sample_reads_per_rung` comes back one entry per rung of `ladder`, in the ladder's own
/// shortest-first order, holding that sample's reads at each. It is emptied and refilled here, so
/// one buffer serves every sample of every locus.
///
/// **The read-group rows are pooled**, through the same [`one_run_per_allele`] the ordinary path's
/// fold uses. A read is a read whichever lane produced it, and asking the rule of each row
/// separately would be a stricter rule applied to exactly the samples carrying more than one
/// library — 157 of 1,707 in a surveyed tomato archive (`doc/devel/ng/spec/read_groups.md` §1).
/// Sharing the grouping with the ordinary path is what stops the two from coming to disagree.
///
/// **Only spanning reads are counted, and that is what the merge's support rows already are**: a
/// read that ran out inside the tract produced a partial, which names no length and is scored on
/// its own axis (spec §1.3, §8). So the histogram's total is exactly
/// [`compared_reads_of`](super::compared_reads_of) — the denominator nomination divides by — and
/// this function deliberately does not return it, because two spellings of one number is how two
/// rules become different rules.
///
/// # Panics
///
/// On a support row naming an allele the ladder was not built for, and on a sample's rows out of
/// ascending allele order — both merge bugs, both held in release (spec §8). In debug it also
/// checks the total against [`compared_reads_of`](super::compared_reads_of), which is the claim
/// the paragraph above makes.
pub(super) fn fill_sample_reads_per_rung(
    sample: &SampleSupport,
    locus: GenomeRegion,
    ladder: &RepeatLadder,
    sample_reads_per_rung: &mut Vec<u32>,
) {
    sample_reads_per_rung.clear();
    sample_reads_per_rung.resize(ladder.rung_count(), 0);
    for rows in super::one_run_per_allele(sample, locus) {
        let rung = ladder.rung_of_table_index(rows[0].allele);
        let pooled_reads = rows.iter().fold(0_u32, |total, row| {
            total.saturating_add(row.support.num_reads)
        });
        sample_reads_per_rung[rung] = sample_reads_per_rung[rung].saturating_add(pooled_reads);
    }
    debug_assert_eq!(
        sample_reads_per_rung.iter().copied().sum::<u32>(),
        super::compared_reads_of(sample),
        "every support row's allele is on some rung, so a sample's reads over the rungs are its \
         compared reads at {locus} — and nomination divides by the second while counting the first"
    );
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::*;
    use super::super::summarise_alleles;
    use super::*;
    use crate::ng::types::ReadGroupId;

    /// A dinucleotide unit, so that a sequence of an odd number of bases has a rung to floor
    /// onto — the case a homopolymer cannot produce.
    fn dinucleotide() -> Motif {
        Motif::new(b"AT").expect("a two-base motif")
    }

    /// Two copies of the genome — the ordinary case, named rather than repeated, since
    /// [`SsrSelectionConfig`] has no default ploidy on purpose.
    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("two copies")
    }

    /// Fold `observation` and build its ladder, returning the scratch that holds both.
    ///
    /// The fold runs first because the ladder sums its per-allele read totals — the coupling
    /// [`build_ladder`] documents, exercised here rather than worked around, so a test cannot
    /// pass against a ladder no fold ever filled.
    fn ladder_over(observation: &CohortObservation) -> SelectionScratch {
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(observation, &dinucleotide(), &mut scratch);
        scratch
    }

    /// Every rung's repeat count and its table indices, shortest rung first.
    fn rungs_of(scratch: &SelectionScratch) -> Vec<(u32, Vec<u32>)> {
        (0..scratch.ladder.rung_count())
            .map(|rung| {
                (
                    scratch.ladder.repeat_count_at(rung),
                    scratch.ladder.table_indices_at(rung).to_vec(),
                )
            })
            .collect()
    }

    /// **The plan's first oracle: two sequences of one length land on one rung and stay
    /// distinct inside it.**
    ///
    /// The two five-repeat sequences differ by their last base — an interrupted repeat, which
    /// spec §1.1 names as a real allele class the read likelihood separates. A ladder that
    /// merged them would hand nomination one sequence where the locus has two, and the second
    /// would never reach the candidate table however many reads showed it.
    #[test]
    fn two_sequences_of_one_length_share_a_rung_and_stay_distinct() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATATAG", b"ATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 3, -3.0), row(2, 3, -3.0)],
            )],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(
            rungs_of(&scratch),
            vec![(4, vec![2]), (5, vec![0, 1])],
            "rungs ascend by repeat count, and the two ten-base sequences share rung 5 as two \
             separate table indices"
        );
        assert_eq!(scratch.ladder.rung_count(), 2);
    }

    /// **The plan's second oracle: a length that is not a whole number of units lands on the
    /// floored rung, and is counted there once.**
    ///
    /// Seven bases at a dinucleotide is three whole units and a base left over, so it joins the
    /// six-base sequence on rung 3. Its reads count toward that rung's total exactly once —
    /// which is what makes the mode below a count of reads rather than of sequences.
    #[test]
    fn a_length_that_is_not_a_whole_number_of_units_lands_on_the_floored_rung() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATA", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 5, -5.0), row(1, 4, -4.0), row(2, 2, -2.0)],
            )],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(
            rungs_of(&scratch),
            vec![(3, vec![0, 1]), (5, vec![2])],
            "the seven-base sequence floors onto rung 3 beside the six-base one, and rung 4 is \
             unoccupied rather than invented"
        );
        assert_eq!(
            scratch.ladder.cohort_reads_at(0),
            9,
            "rung 3 carries both sequences' reads, each counted once"
        );
    }

    /// **The plan's third oracle, first half: the mode is the rung with the most reads — and it
    /// is not the reference's.**
    ///
    /// The reference tract is five repeats and two reads; the panel's own length is four repeats
    /// and thirty. A mode keyed on the reference would centre the periodicity test a whole unit
    /// away from where the cohort's reads actually sit.
    #[test]
    fn the_mode_is_the_rung_with_the_most_reads_and_not_the_reference_s() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 2, -2.0), row(1, 30, -30.0)])],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(scratch.ladder.modal_repeat_count(), 4);
    }

    /// **The plan's third oracle, second half: a tie breaks toward the shorter rung.**
    ///
    /// Both rungs carry twelve reads. The rule is not a measured preference — it is what keeps
    /// the answer independent of the order the samples were walked in, which spec §8 requires
    /// because the run's output must be byte-identical at any worker count.
    #[test]
    fn the_mode_breaks_a_tie_toward_the_shorter_rung() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 12, -12.0), row(1, 12, -12.0)],
            )],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(
            scratch.ladder.modal_repeat_count(),
            4,
            "twelve reads each, so the shorter rung takes it"
        );
    }

    /// A rung's reads are **every covering sample's**, including a sample whose own reads never
    /// reached the support rule.
    ///
    /// The second sample shows one read at the ten-base sequence against a bar of two, so it
    /// clears nothing — and its read is still on that rung, because the mode asks where the
    /// cohort's reads are and not which sequences were worth calling over. Getting this wrong
    /// costs the mode at exactly the loci where one deep sample carries a length the rest of the
    /// panel shows thinly.
    #[test]
    fn a_rungs_reads_include_a_sample_that_cleared_no_bar() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 6, -6.0), row(1, 6, -6.0)]),
                sample_showing(1, vec![row(1, 1, -1.0)]),
            ],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(
            scratch.ladder.cohort_reads_at(1),
            7,
            "six reads from the sample that cleared the bar and one from the sample that did not"
        );
        assert_eq!(
            scratch.ladder.modal_repeat_count(),
            5,
            "seven reads against six, decided by the read the bar refused"
        );
    }

    /// **A repeat count the table holds answers with its rung; one it does not answers `None`** —
    /// the lookup the `±1` rescue is built on, so that nothing invents a length.
    ///
    /// Rung 4 is missing here between rungs 3 and 5, which is the shape that makes the answer
    /// mean something: a search that returned the nearest rung instead of `None` would let the
    /// rescue promote a length nothing in the cohort ever showed.
    #[test]
    fn a_repeat_count_the_table_does_not_hold_has_no_rung() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 4, -4.0), row(1, 4, -4.0)])],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(scratch.ladder.rung_of_repeat_count(3), Some(0));
        assert_eq!(scratch.ladder.rung_of_repeat_count(5), Some(1));
        assert_eq!(
            scratch.ladder.rung_of_repeat_count(4),
            None,
            "no sequence of four repeats is in the table, so the rescue has no length to promote \
             there"
        );
    }

    /// **The reference tract's rung exists with no reads on it, which is where this lookup and
    /// production's `occupied` test come apart.**
    ///
    /// The merge interns the reference at index 0 whether or not a read landed on it, so at a
    /// tract where every read shows the panel's own length the reference still has a rung — one
    /// carrying zero reads. Production's rescue asks `cohort_support(length) > 0` and would
    /// refuse to promote that length; a caller here that asked only for the rung would promote
    /// it. Pinned so that the step which builds the rescue has the difference in front of it
    /// rather than in a doc comment.
    #[test]
    fn the_reference_rung_exists_with_no_reads_on_it() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(1, 8, -8.0)])],
        );
        let scratch = ladder_over(&observation);

        assert_eq!(
            scratch.ladder.rung_of_repeat_count(5),
            Some(1),
            "the reference's own length has a rung, because the merge interns the reference"
        );
        assert_eq!(
            scratch.ladder.cohort_reads_at(1),
            0,
            "and no read reached it — which is what production's occupancy test asks"
        );
    }

    /// A homopolymer puts every sequence on its own rung, and the ladder still ascends.
    ///
    /// The period-one case is not a corner: it is the commonest tract class in both benchmarks,
    /// and it is the one where floor division is the identity, so a keying bug that survives
    /// here would survive genome-wide on homopolymers.
    #[test]
    fn a_homopolymer_keys_every_length_to_its_own_rung() {
        let observation = locus_of(
            &[b"AAAA", b"AAA", b"AAAAA"],
            vec![sample_showing(
                0,
                vec![row(0, 3, -3.0), row(1, 9, -9.0), row(2, 3, -3.0)],
            )],
        );
        let mut scratch = SelectionScratch::new();
        summarise_alleles(&observation, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(
            &observation,
            &Motif::new(b"A").expect("a one-base motif"),
            &mut scratch,
        );

        assert_eq!(
            rungs_of(&scratch),
            vec![(3, vec![1]), (4, vec![0]), (5, vec![2])],
            "three lengths, three rungs, ascending"
        );
        assert_eq!(scratch.ladder.modal_repeat_count(), 3);
    }

    /// Building a second locus's ladder over the first's leaves nothing of the first behind.
    ///
    /// The scratch is reused across every locus a worker sees, and a rung, an index or a mode
    /// carried between two loci is a wrong nomination that nothing downstream can detect. The
    /// second locus here is deliberately **smaller** than the first, which is the direction that
    /// catches a buffer emptied by `truncate` to its own length.
    #[test]
    fn a_second_locus_leaves_no_rung_of_the_first() {
        let first = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 2, -2.0), row(1, 9, -9.0), row(2, 2, -2.0)],
            )],
        );
        let mut scratch = ladder_over(&first);
        assert_eq!(scratch.ladder.rung_count(), 3);
        assert_eq!(scratch.ladder.modal_repeat_count(), 4);

        let second = locus_of(
            &[b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 5, -5.0)])],
        );
        summarise_alleles(&second, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(&second, &dinucleotide(), &mut scratch);

        assert_eq!(rungs_of(&scratch), vec![(5, vec![0])]);
        assert_eq!(scratch.ladder.modal_repeat_count(), 5);
        assert_eq!(
            scratch.ladder.rung_of_repeat_count(4),
            None,
            "the first locus's rung 4 is gone, not merely unreachable"
        );
    }

    /// A ladder nobody built refuses the two questions whose honest answers are indistinguishable
    /// from real ones: a mode of 0 is a legal repeat count, and `None` from the occupancy test
    /// reads as "no read reached that length".
    #[test]
    #[should_panic(expected = "has not been built for a locus")]
    fn an_unbuilt_ladder_refuses_its_mode() {
        RepeatLadder::default().modal_repeat_count();
    }

    /// The companion of the test above, on the other accessor that would answer plausibly.
    #[test]
    #[should_panic(expected = "has not been built for a locus")]
    fn an_unbuilt_ladder_refuses_an_occupancy_question() {
        RepeatLadder::default().rung_of_repeat_count(4);
    }

    /// **A tract that reached the build without being folded is refused, not appended to.**
    ///
    /// The ladder is emptied by [`SelectionScratch::reset_for`], which the fold calls, so the
    /// build does not empty it a second time. That leaves one failing state — a build with no
    /// fold before it — and this is what makes it loud instead of a ladder carrying two loci's
    /// sequences at once.
    #[test]
    #[should_panic(expected = "the ladder already holds 2 rung(s) over 2 index/indices")]
    fn building_a_ladder_twice_without_a_fold_between_is_refused() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 3, -3.0), row(1, 5, -5.0)])],
        );
        let mut scratch = ladder_over(&observation);
        build_ladder(&observation, &dinucleotide(), &mut scratch);
    }

    /// A ladder built for one locus refuses a rung another locus would have had.
    #[test]
    #[should_panic(expected = "rung 2 of a ladder holding 1")]
    fn a_ladder_refuses_a_rung_it_does_not_hold() {
        let observation = locus_of(
            &[b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 5, -5.0)])],
        );
        let scratch = ladder_over(&observation);
        scratch.ladder.repeat_count_at(2);
    }

    /// The ladder sums the fold's per-allele totals, so a scratch folded on another locus is
    /// refused rather than summed onto this locus's rungs.
    #[test]
    #[should_panic(expected = "the scratch's fold must be this locus's")]
    fn a_ladder_refuses_a_fold_from_another_locus() {
        let folded = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 3, -3.0), row(1, 3, -3.0)])],
        );
        let mut scratch = SelectionScratch::new();
        summarise_alleles(&folded, support_rule_of(2, 0.0), &mut scratch);

        let other = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 3, -3.0), row(1, 3, -3.0), row(2, 3, -3.0)],
            )],
        );
        build_ladder(&other, &dinucleotide(), &mut scratch);
    }

    // ---- B2: the configuration, and the per-sample length histogram ----

    /// The two caps, side by side. **The whole of the tract path's departure from the ordinary
    /// path's cap is this number**, and stating both in one assertion is what stops a later edit
    /// from quietly making them the same again.
    #[test]
    fn the_tract_cap_is_thirty_two_where_the_ordinary_paths_is_six() {
        assert_eq!(DEFAULT_MAX_CANDIDATE_ALLELES_SSR.get(), 32);
        assert_eq!(super::super::DEFAULT_MAX_CANDIDATE_ALLELES.get(), 6);
        assert_eq!(
            SsrSelectionConfig::at_ploidy(diploid())
                .shared
                .max_candidate_alleles,
            DEFAULT_MAX_CANDIDATE_ALLELES_SSR,
            "the config must carry the tract cap and not the shared default"
        );
    }

    /// **The support rule is the ordinary path's, unchanged** — including the share, which spec §5
    /// writes as 5 in 100 and the ordinary path ships at 10.
    ///
    /// Read `SsrSelectionConfig::at_ploidy`'s own documentation for why: the spec's stated reason
    /// for its number is that one number should govern both paths, and the ordinary path's moved.
    /// This test exists so that the choice is a line someone has to edit rather than a default
    /// that drifts.
    #[test]
    fn the_support_rule_is_the_ordinary_paths_share_and_floor() {
        let config = SsrSelectionConfig::at_ploidy(diploid());
        assert_eq!(
            config.shared.min_allele_support,
            super::super::DEFAULT_MIN_ALLELE_SUPPORT
        );
        assert_eq!(config.shared.min_allele_support.share.get(), 0.10);
    }

    /// The ploidy is the caller's and reaches the config unchanged — a triploid run promotes three
    /// rungs a sample, not two.
    #[test]
    fn the_ploidy_is_the_callers() {
        let triploid = Ploidy::try_new(3).expect("three copies");
        assert_eq!(SsrSelectionConfig::at_ploidy(triploid).ploidy, triploid);
        assert_eq!(SsrSelectionConfig::at_ploidy(diploid()).ploidy, diploid());
    }

    /// The off-grid share refuses everything that is not a fraction of one, because neither
    /// failure crashes: above one the periodicity verdict can never be reached, below zero every
    /// tract returns the reference alone.
    #[test]
    fn the_off_grid_share_refuses_anything_that_is_not_a_fraction_of_one() {
        assert_eq!(
            MaxOffGridShare::new(0.0).map(MaxOffGridShare::get),
            Some(0.0)
        );
        assert_eq!(
            MaxOffGridShare::new(1.0).map(MaxOffGridShare::get),
            Some(1.0)
        );
        assert!(MaxOffGridShare::new(-0.001).is_none());
        assert!(MaxOffGridShare::new(1.001).is_none());
        assert!(MaxOffGridShare::new(f64::NAN).is_none());
        assert!(MaxOffGridShare::new(f64::INFINITY).is_none());
        assert_eq!(
            DEFAULT_MAX_OFF_GRID_SHARE.get(),
            0.10,
            "production's max_out_of_frame_frac, inherited and never measured"
        );
    }

    /// **A sample's reads land on the rungs its sequences sit on**, and a rung it showed nothing
    /// at is a zero rather than a missing entry.
    ///
    /// Three rungs at 3, 4 and 5 repeats; the sample shows the outer two and not the middle one.
    /// The histogram is parallel to the ladder's rungs, so the middle zero has to be *there* —
    /// nomination walks the rungs by index and a compacted histogram would shift every count onto
    /// the wrong length.
    #[test]
    fn a_samples_reads_land_on_the_rungs_its_sequences_sit_on() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 4, -4.0), row(2, 6, -6.0)]),
                sample_showing(1, vec![row(1, 5, -5.0)]),
            ],
        );
        let scratch = ladder_over(&observation);
        assert_eq!(
            rungs_of(&scratch),
            vec![(3, vec![0]), (4, vec![1]), (5, vec![2])]
        );

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[0],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(
            histogram,
            vec![4, 0, 6],
            "four reads at three repeats, none at four, six at five"
        );
    }

    /// **One sample's two read groups are one sample**, pooled onto one rung.
    ///
    /// The same rule as the ordinary path's fold, and pooled through the same helper — asking it
    /// of each row separately would be a stricter rule applied to exactly the samples carrying
    /// more than one library.
    #[test]
    fn one_samples_two_read_groups_land_on_one_rung() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 2, -2.0),
                    row_from_group(1, ReadGroupId(0), 3, -3.0),
                    row_from_group(1, ReadGroupId(1), 4, -4.0),
                ],
            )],
        );
        let scratch = ladder_over(&observation);

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[0],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(
            histogram,
            vec![2, 7],
            "three reads from one lane and four from the other are seven reads at four repeats"
        );
    }

    /// **Two spellings of one length are two counts on one rung, added.**
    ///
    /// This is the case the rung's count exists for and the one a per-allele overwrite would get
    /// wrong while every other fixture stayed green: an interrupted repeat gives a sample two
    /// distinct tract sequences of the same length, and nomination asks whether *that length* has
    /// enough reads. Counting only the last-seen spelling would refuse a length the sample plainly
    /// carries.
    #[test]
    fn two_spellings_of_one_length_add_on_their_shared_rung() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATATAG", b"ATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 5, -5.0), row(2, 2, -2.0)],
            )],
        );
        let scratch = ladder_over(&observation);
        assert_eq!(rungs_of(&scratch), vec![(4, vec![2]), (5, vec![0, 1])]);

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[0],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(
            histogram,
            vec![2, 9],
            "four reads and five reads at five repeats are nine reads at five repeats"
        );
    }

    /// A covering sample whose reads all stopped inside the tract counts nothing at any rung.
    ///
    /// Partials say the sample carries *at least* this much of the tract, not what length it
    /// carries, so they name no rung and are scored on their own axis. A histogram that counted
    /// them would nominate a length from reads that never crossed it.
    #[test]
    fn a_sample_with_only_partials_counts_at_no_rung() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 3, -3.0), row(1, 3, -3.0)]),
                sample_with_only_partials(1, 40),
            ],
        );
        let scratch = ladder_over(&observation);

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[1],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(
            histogram,
            vec![0, 0],
            "forty partial reads name no length, so every rung is zero"
        );
    }

    /// **Refilling for a second sample leaves none of the first's counts** — the buffer is one per
    /// worker, walked sample by sample, and a count carried across would nominate a length for a
    /// sample that never showed it.
    #[test]
    fn refilling_for_a_second_sample_leaves_none_of_the_firsts_counts() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 9, -9.0), row(1, 9, -9.0), row(2, 9, -9.0)]),
                sample_showing(1, vec![row(1, 2, -2.0)]),
            ],
        );
        let scratch = ladder_over(&observation);

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[0],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(histogram, vec![9, 9, 9]);
        fill_sample_reads_per_rung(
            &observation.per_sample[1],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(
            histogram,
            vec![0, 2, 0],
            "the second sample's two reads, and nothing of the first's twenty-seven"
        );
    }

    /// **The histogram's total is the sample's compared reads** — the denominator nomination
    /// divides by, which is why this function does not return a total of its own.
    ///
    /// Asserted here against `compared_reads_of` on a sample carrying two lanes and three rungs,
    /// so that the two counts are compared on an input where a per-row rule and a per-sample rule
    /// would differ.
    #[test]
    fn the_histogram_totals_the_samples_compared_reads() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 2, -2.0),
                    row_from_group(1, ReadGroupId(0), 3, -3.0),
                    row_from_group(1, ReadGroupId(1), 4, -4.0),
                    row(2, 5, -5.0),
                ],
            )],
        );
        let scratch = ladder_over(&observation);

        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(
            &observation.per_sample[0],
            observation.region,
            &scratch.ladder,
            &mut histogram,
        );
        assert_eq!(histogram.iter().sum::<u32>(), 14);
        assert_eq!(
            histogram.iter().sum::<u32>(),
            super::super::compared_reads_of(&observation.per_sample[0]),
            "one producer for the numerator's rungs and the denominator's total"
        );
    }

    /// A support row naming an allele the tract's table does not hold is refused, rather than
    /// counted onto whichever rung the index happens to reach.
    #[test]
    #[should_panic(expected = "named allele 5 of a repeat tract whose table holds 2")]
    fn a_row_naming_an_allele_the_tract_does_not_hold_is_refused() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 3, -3.0), row(1, 3, -3.0)])],
        );
        let scratch = ladder_over(&observation);

        let stray = sample_showing(0, vec![row(5, 3, -3.0)]);
        let mut histogram = Vec::new();
        fill_sample_reads_per_rung(&stray, observation.region, &scratch.ladder, &mut histogram);
    }
}
