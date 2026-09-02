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
//! histogram (B2), nomination — each sample's own repeat counts, the `±1` rescue and the cohort's
//! union (C1, C2) — the admission of the sequences on a promoted rung (D1), and the periodicity
//! verdict with the entry point [`select_ssr`] it gates (D2). The two numbers the genotype prior
//! takes are the step after those, and nothing outside this module calls anything here yet.
//!
//! **Which is why a non-test build expects everything here to be dead.** The expectation rather
//! than an `allow`, so that the first real caller — `select_ssr`, which the steps after this one
//! build — turns the line below into a compile error and deletes it. Under `cfg(test)` nothing is
//! exempted: the tests in this file are what exercise the ladder today, and they must keep
//! failing when it is wrong.
#![cfg_attr(not(test), expect(dead_code))]

use super::summarise_alleles;
use super::{
    AlleleRemap, CandidateSelectionConfig, LocusSelection, MaxCandidateAlleles, SelectionScratch,
    SelectionVerdict,
};
use crate::ng::calling::CandidateAlleles;
use crate::ng::locus_generation::{LocusKind, SsrDetail};
use crate::ng::run::cohort_merge::MinAltReadShare;
use crate::ng::run::cohort_merge::build::{CohortObservation, SampleSupport};
use crate::ng::types::{AlleleId, GenomeRegion, Motif, Ploidy};

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
        promoted_rungs: _,
        rung_is_promoted: _,
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

/// **Which repeat counts one sample puts forward** — the rungs whose reads cleared the shared
/// support rule for this sample, cut to the best `ploidy` of them (spec §4; arch §3.1).
///
/// `promoted` comes back holding **rung indices, ascending**, at most `ploidy` of them.
///
/// # The rule is the shared one, asked of a rung
///
/// A repeat count is nominated when this sample's reads at it reach
/// `max(2 reads, ceil(share × this sample's spanning reads))` — [`MinAltReads::reached_by`], the
/// same predicate the merge asks of a sample's non-reference reads and the ordinary path asks of
/// one sequence, with the rung's read total as the numerator. **Not a second predicate**: a second
/// spelling of one rule is how two rules become different rules, which is why neither the
/// numerator nor the denominator is computed here.
///
/// **The denominator is the sample's own spanning reads**, so nothing about which other samples
/// are in the run can change what this one nominates. That is what makes a cohort of one and a
/// cohort of a thousand give this sample the same answer.
///
/// # Top `ploidy`, ties to the shorter length
///
/// A diploid sample carries two copies, so at most two lengths are real and the rest are stutter
/// and error; `ploidy` comes from the run's parameters, so a triploid region promotes three. Ties
/// break toward the **shorter** repeat count, which is production's rule
/// ([`candidate_set.rs:231-238`](../../../../src/ssr/cohort/candidate_set.rs)) and is kept for
/// determinism rather than for a measured reason: at a heterozygote whose two lengths carry
/// exactly equal reads the answer must not depend on the order the rungs were walked, since the
/// run's output has to be byte-identical at any worker count.
///
/// # What is not ported: production's clear-peak test
///
/// Production nominates a length only if its reads exceed **both** neighbouring lengths by more
/// than three (`is_clear_peak`, [`rung_ladder.rs:274`](../../../../src/ssr/cohort/rung_ladder.rs)).
/// A heterozygote whose two copies differ by one repeat is then invisible: neither length is a
/// peak, because each has the other beside it. Spec §4.1 measures what that costs — both alleles
/// offered at 33–78% of such tracts against 97–100% for this rule, and with **fewer** candidates,
/// not more. Nothing here reads a neighbour.
///
/// # Panics
///
/// On a `sample_reads_per_rung` that is not this ladder's — one entry per rung is what makes a
/// rung index mean the same thing in both. A histogram of another locus's width would nominate
/// lengths this tract does not hold.
pub(super) fn promote_rungs_for_sample(
    sample_reads_per_rung: &[u32],
    compared_reads: u32,
    config: &SsrSelectionConfig,
    ladder: &RepeatLadder,
    promoted: &mut Vec<u32>,
) {
    assert_eq!(
        sample_reads_per_rung.len(),
        ladder.rung_count(),
        "the histogram runs parallel to the ladder's rungs, one entry each"
    );
    promoted.clear();
    promoted.extend(
        (0..ladder.rung_count())
            .filter(|&rung| {
                config
                    .shared
                    .min_allele_support
                    .reached_by(sample_reads_per_rung[rung], compared_reads)
            })
            .map(|rung| u32::try_from(rung).expect("a tract with fewer than four billion rungs")),
    );

    // Two orders in turn, as the ordinary path's cap does: rank to decide *which* rungs survive,
    // then back to ascending rung order, because that is the order everything downstream reads
    // them in and it is the merge's own order underneath.
    let copies = usize::from(config.ploidy.get());
    if promoted.len() > copies {
        promoted.sort_unstable_by_key(|&rung| {
            (
                std::cmp::Reverse(sample_reads_per_rung[rung as usize]),
                rung,
            )
        });
        promoted.truncate(copies);
        promoted.sort_unstable();
    }
}

/// **When a sample resolved fewer lengths than it has copies, put forward the neighbours of what
/// it did resolve — but only lengths some sample's reads actually reached** (spec §4; production's
/// `occupied` test, [`candidate_set.rs:239-258`](../../../../src/ssr/cohort/candidate_set.rs)).
///
/// A diploid sample that resolved one length has a copy unaccounted for, and what it most likely
/// carries there is a neighbour of what it did resolve — a second allele one repeat away, hidden
/// under the first by stutter. This is the one part of production's nomination ng keeps.
///
/// **"Occupied" is what stops it inventing a length.** A neighbour is put forward only where the
/// cohort's reads reached it; without that test the rescue would offer a repeat count nothing in
/// the run has ever seen, at every under-resolved sample of every tract.
///
/// **And it fires only on an under-resolved sample.** Firing it on a sample that did resolve its
/// ploidy would widen every locus by up to two rungs a sample, which is extra candidates rather
/// than a crash — every one of them a column in every genotype table for the life of the locus.
///
/// **The rescue is not itself capped at `ploidy`**, which is production's behaviour ported
/// unchanged: a sample that resolved one length of three occupied ones can come back with three.
/// The cap that does bind is the shared one over *sequences*, applied at admission.
///
/// # Occupied means reads, not a rung
///
/// The test is `the ladder has a rung at that count **and** its cohort reads are non-zero`, which
/// is production's `cohort_support(length) > 0`. The two come apart at exactly one rung: the merge
/// interns the reference tract at index 0 whether or not a read landed on it, so the reference's
/// length always has a rung and may have no reads. Promoting it would put forward a length the
/// cohort never showed — and would gain nothing even if it did, since the reference sequence is
/// admitted first and exempt from the bar regardless of whether its rung was promoted.
pub(super) fn rescue_occupied_neighbours(
    ladder: &RepeatLadder,
    ploidy: Ploidy,
    promoted: &mut Vec<u32>,
) {
    if promoted.len() >= usize::from(ploidy.get()) {
        return;
    }
    let resolved = promoted.len();
    for index in 0..resolved {
        let repeat_count = ladder.repeat_count_at(promoted[index] as usize);
        for neighbour in [repeat_count.checked_sub(1), repeat_count.checked_add(1)]
            .into_iter()
            .flatten()
        {
            let Some(rung) = ladder.rung_of_repeat_count(neighbour) else {
                continue;
            };
            if ladder.cohort_reads_at(rung) > 0 {
                promoted.push(rung as u32);
            }
        }
    }
    // Ascending and once each: a length can be the neighbour of two resolved lengths, and the
    // rungs a sample resolved are already in the list.
    promoted.sort_unstable();
    promoted.dedup();
}

/// **The repeat counts the whole cohort puts forward: the union of what each sample does**
/// (spec §4; arch §3.1).
///
/// `rung_is_promoted` comes back one flag per rung of the ladder, true where **some** sample
/// nominated that length or the rescue reached it for that sample.
///
/// **A union and not a vote**, for the same reason one sample reaching the bar admits a sequence
/// for the whole cohort on the ordinary path: an allele one accession of sixty-three carries is
/// still an allele, and a rule that needed two would delete exactly the rare variation a cohort is
/// sequenced to find. The cost of the union is candidates, and that is what the cap is for.
///
/// **Every sample is asked its own question against its own reads**, so this loop is the only
/// place the cohort appears at all — nothing in the bar or the cut reads it.
///
/// # Panics
///
/// Through the functions it calls: on a scratch whose ladder is not this locus's, and on a
/// sample's support rows out of ascending allele order.
pub(super) fn promote_rungs_for_cohort(
    observation: &CohortObservation,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) {
    let SelectionScratch {
        ladder,
        sample_reads_per_rung,
        promoted_rungs,
        rung_is_promoted,
        ..
    } = scratch;
    rung_is_promoted.clear();
    rung_is_promoted.resize(ladder.rung_count(), false);

    for sample in &observation.per_sample {
        fill_sample_reads_per_rung(sample, observation.region, ladder, sample_reads_per_rung);
        promote_rungs_for_sample(
            sample_reads_per_rung,
            super::compared_reads_of(sample),
            config,
            ladder,
            promoted_rungs,
        );
        rescue_occupied_neighbours(ladder, config.ploidy, promoted_rungs);
        for &rung in promoted_rungs.iter() {
            rung_is_promoted[rung as usize] = true;
        }
    }
}

/// **Admit the sequences on the promoted rungs that some sample's reads earned** — the last of
/// nomination's three passes, and the one that turns rungs back into sequences (spec §5;
/// arch §3.2).
///
/// # Every spelling stands on its own reads
///
/// A promoted rung says *this length is worth calling over*; it does not say which sequences at
/// that length are. Each one faces the shared support rule asked of the **sequence** — the same
/// `max(2 reads, share × a sample's spanning reads)` the ordinary path asks, already folded into
/// the per-allele summaries by [`summarise_alleles`](super::summarise_alleles).
///
/// **No representative is privileged and no recurrence term applies.** Production promotes the
/// rung's best-supported sequence unconditionally and makes any sibling clear three further
/// gates — 8 reads, **3 distinct samples**, and a tenth of the rung's reads
/// ([`candidate_set.rs:169-191`](../../../../src/ssr/cohort/candidate_set.rs)). The three-sample
/// term has no cohort-size clamp, so **below three samples no second spelling can ever be
/// promoted** and the mechanism is simply absent. Spec §5 measures what that costs at the class
/// it matters for: at a heterozygote whose two copies are the same length spelled differently —
/// 296 of HG002's 695 heterozygous tracts — production offers both at 35.1% of them and this rule
/// at 86.1%, against a ceiling of 93.6% set by what some read actually carried.
///
/// # The reference tract is admitted first and asked nothing
///
/// It is exempt from the bar *and* from the promotion test (spec §5, §7), so a tract whose reads
/// all sit somewhere else still yields a table the caller can use. It is seeded structurally
/// rather than by reading whether it passed, which it may well have.
///
/// # The cap, and the leftover the tract likelihood does not read
///
/// Above the cap the list is cut to the best-ranked by [`compare_best_first`](super::compare_best_first)
/// and the locus is still called — the shared rule, and here the cap is
/// [`DEFAULT_MAX_CANDIDATE_ALLELES_SSR`] rather than the ordinary path's six.
///
/// **The leftover is filled although this path's likelihood never reads its pool.** Spec §8 is
/// explicit that a read no candidate explains is already carried by the junk term, spread over the
/// tract lengths the stutter model can reach — so `q_sum` here is computed and unread. What *is*
/// read is the other count: [`UnmatchedSupport::earned_reads_cut_by_the_cap`](super::UnmatchedSupport::earned_reads_cut_by_the_cap),
/// which says this sample carried a length the locus is no longer called over and must therefore
/// be emitted as missing. That rule is the same on both paths, so the leftover is built by the
/// same [`leftover_of`](super::leftover_of) and not by a second walk.
///
/// # Panics
///
/// On an empty allele table, and on a promotion flag list that is not this ladder's. Both are
/// caller bugs whose symptom is a wrong candidate list rather than a crash (spec §8).
pub(super) fn admit_promoted_sequences(
    observation: &CohortObservation,
    detail: &SsrDetail,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> LocusSelection {
    assert!(
        !observation.alleles.is_empty(),
        "a cohort locus always holds at least its reference allele, and the repeat tract at {} \
         holds none",
        observation.region
    );
    let SelectionScratch {
        per_allele,
        ranked_table_indices,
        ladder,
        rung_is_promoted,
        ..
    } = scratch;
    assert_eq!(
        rung_is_promoted.len(),
        ladder.rung_count(),
        "the promotion flags run parallel to the ladder's rungs, one each, at {}",
        observation.region
    );

    // Every alternative that sits on a promoted rung *and* some sample's reads earned, in the
    // merge table's own order — which is the order it is admitted in, so nothing has to sort it
    // unless the cap reorders it below.
    ranked_table_indices.clear();
    ranked_table_indices.extend(
        (1..observation.alleles.len())
            .filter(|&index| {
                rung_is_promoted[ladder.rung_of_table_index(index)]
                    && per_allele[index].cleared_the_bar()
            })
            .map(|index| {
                u32::try_from(index).expect("a merge table narrower than four billion alleles")
            }),
    );

    let allowed_alternatives = usize::from(config.shared.max_candidate_alleles.alternatives());
    let verdict = if ranked_table_indices.len() <= allowed_alternatives {
        SelectionVerdict::Selected
    } else {
        let alternative_of = |table_index: u32| super::RankedAlternative {
            summary: per_allele[table_index as usize],
            bases: &observation.alleles[table_index as usize],
        };
        ranked_table_indices.sort_unstable_by(|&left, &right| {
            super::compare_best_first(alternative_of(left), alternative_of(right))
        });
        let dropped = ranked_table_indices.len() - allowed_alternatives;
        ranked_table_indices.truncate(allowed_alternatives);
        ranked_table_indices.sort_unstable();
        SelectionVerdict::Truncated {
            dropped: u32::try_from(dropped).expect("a merge table narrower than four billion"),
        }
    };

    let mut alleles = CandidateAlleles::new(
        observation.alleles[0].clone(),
        LocusKind::Ssr(detail.clone()),
    );
    let mut remap = AlleleRemap::with_all_dropped(observation.alleles.len());
    remap.admit(0, AlleleId::REFERENCE);
    for &table_index in ranked_table_indices.iter() {
        let table_index = table_index as usize;
        let candidate = alleles.admit(observation.alleles[table_index].clone());
        remap.admit(table_index, candidate);
    }

    let covering_samples = observation.per_sample.len();
    let leftovers = observation
        .per_sample
        .iter()
        .map(|sample| {
            super::leftover_of(
                sample,
                observation.region,
                &remap,
                config.shared.min_allele_support,
            )
        })
        .collect();
    LocusSelection::new(alleles, verdict, leftovers, remap, covering_samples)
}

/// **Whether this tract's reads actually vary in whole motif units** — the one verdict this path
/// adds, and the gate in front of everything above (spec §7; arch §3.3).
///
/// A stretch the catalog called a repeat tract whose reads sit at lengths the motif cannot
/// explain is not a tract this caller's model describes: the stutter distribution is written on
/// whole-repeat and part-repeat regimes and the prior's ladder is written on repeat counts.
/// Genotyping it against that model would produce a confident answer from a model that does not
/// apply.
///
/// # The grid is anchored on the reference tract's own length
///
/// A read is **off the grid** when the difference between its tract length and the reference
/// tract's is not a whole number of motif units. A sample is non-periodic when more than
/// `max_off_grid_share` of its spanning reads are off the grid, and **the locus is judged
/// non-periodic only when no sample is periodic** — the same "one sample suffices" shape as
/// everything else on this path.
///
/// **The anchor is a decision, taken by the owner on 2026-09-02, and it is not what the design
/// documents say.** Arch §3.3 and spec §7 anchor the grid on the ladder's mode, and spec §3 adds
/// that ng measures it in units where production measures it in bases. Taken literally those two
/// give a grid anchored at **zero** — the mode is a repeat count, so its length in bases is a whole
/// number of units and cancels out of the subtraction — and a zero-anchored grid refuses a real
/// class of tract:
///
/// - the catalog trims every tract back to whole motif copies **at both ends**, but a
///   length-changing interruption inside puts the two ends out of phase with each other, so the
///   tract's own reference length is then **not** a multiple of the period;
/// - such a tract is admitted whenever the interruption is late enough to clear the catalog's
///   purity floor of 0.8. Measured through the catalog's own `minimal_trim` and `recompute_purity`
///   on 2026-09-02: 49 bases of an `AT` repeat with one extra base 40 bases in trims to 49 bases at
///   **purity 0.816**, and 49 is odd;
/// - at that tract every read at the reference length is off a zero-anchored grid, so every sample
///   is non-periodic and **the locus is refused and never called**.
///
/// Production avoids that by anchoring on the commonest observed length in bases
/// ([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs)) — its own comment
/// says so: *"the grid is anchored to the modal length, not to zero, so an interrupted repeat
/// sitting at an odd reference length stays periodic"*. **The reference tract's length is the
/// third anchor and the one taken here**, because it keeps those tracts *and* is a property of the
/// locus rather than of the reads: the commonest observed length moves with depth, so a shallow
/// sample could shift which lengths count as on-grid. It is also the quantity the genotype prior
/// was re-indexed onto on 2026-08-27 — offset from the reference tract length
/// (`doc/devel/ng/spec/population_diversity.md` §4.2) — so the periodicity grid and the prior's
/// index are now the same number.
///
/// # What the counting includes
///
/// **Spanning reads only.** A read that ran out inside the tract names no length, so it is on no
/// grid and in no denominator (spec §8) — it is not in the merge's support rows at all.
///
/// **A sample with no spanning reads is not asked**, rather than counted as periodic. Counting it
/// as periodic would make one silent sample enough to save every locus, and a silent sample is
/// exactly what a tract too long for a read to span produces. A locus where **no** sample has a
/// spanning read is periodic by default: there is nothing to judge it on, and refusing it would be
/// a verdict about coverage rather than about the tract.
///
/// **A homopolymer can never be off the grid**, since every difference is a whole number of
/// one-base units. That falls out of the arithmetic rather than needing production's explicit
/// short-circuit.
///
/// # Panics
///
/// On an empty allele table, and on a sample's rows out of ascending allele order.
pub(super) fn locus_is_periodic(
    observation: &CohortObservation,
    motif: &Motif,
    max_off_grid_share: MaxOffGridShare,
) -> bool {
    assert!(
        !observation.alleles.is_empty(),
        "a cohort locus always holds at least its reference allele, and the repeat tract at {} \
         holds none",
        observation.region
    );
    let period = motif.period();
    let reference_len = observation.alleles[0].len();
    let mut any_sample_judged = false;

    for sample in &observation.per_sample {
        let spanning_reads = super::compared_reads_of(sample);
        if spanning_reads == 0 {
            continue;
        }
        any_sample_judged = true;
        let mut off_grid_reads = 0_u32;
        for rows in super::one_run_per_allele(sample, observation.region) {
            let length = observation.alleles[rows[0].allele].len();
            if !length.abs_diff(reference_len).is_multiple_of(period) {
                let pooled_reads = rows.iter().fold(0_u32, |total, row| {
                    total.saturating_add(row.support.num_reads)
                });
                off_grid_reads = off_grid_reads.saturating_add(pooled_reads);
            }
        }
        if f64::from(off_grid_reads) <= max_off_grid_share.get() * f64::from(spanning_reads) {
            return true;
        }
    }
    // Nobody had a spanning read to judge, so there is no evidence the tract is not a tract.
    !any_sample_judged
}

/// **Narrow one repeat tract's allele table to the tract sequences worth calling over**
/// (spec §4, §5, §7) — the entry point, and the repeat-tract sibling of
/// [`select_generic`](super::generic::select_generic).
///
/// **Four passes over the locus, in this order**, because each needs the one before: fold the
/// rows into per-allele summaries; build the ladder over them; nominate rungs per sample and union
/// them; then admit the sequences on nominated rungs that cleared the bar, apply the cap and fill
/// the leftover.
///
/// **The periodicity verdict runs first of all**, and a tract that fails it returns the reference
/// tract alone with [`SelectionVerdict::NotPeriodic`] — a usable table rather than a refusal, so
/// what the run does with such a locus stays emission's decision (spec §7). Its leftover is one
/// zeroed entry per covering sample: nothing was cut by the cap, so no sample is uncallable for
/// that reason, and the tract likelihood does not read the pool at all (spec §8).
///
/// # Panics
///
/// On a cohort observation that is not a repeat tract. Which path a locus takes is decided by its
/// kind at the driver's dispatch, so reaching here with a `Generic` or bundle locus is a routing
/// bug, and the alternative — scoring a SNP against a stutter model — is a confident wrong answer
/// rather than a crash.
pub(super) fn select_ssr(
    observation: &CohortObservation,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> LocusSelection {
    let LocusKind::Ssr(detail) = &observation.kind else {
        panic!(
            "the repeat-tract path was handed a {:?} locus at {}: which path a locus takes is \
             decided by its kind at the dispatch, and scoring one kind against the other's model \
             is a confident wrong answer rather than a failure",
            observation.kind, observation.region
        );
    };

    if !locus_is_periodic(observation, &detail.motif, config.max_off_grid_share) {
        return reference_tract_alone(observation, detail, SelectionVerdict::NotPeriodic);
    }

    summarise_alleles(observation, config.shared.min_allele_support, scratch);
    build_ladder(observation, &detail.motif, scratch);
    promote_rungs_for_cohort(observation, config, scratch);
    admit_promoted_sequences(observation, detail, config, scratch)
}

/// A table holding the reference tract and nothing else, under `verdict` — what a tract this path
/// declines to narrow still yields.
fn reference_tract_alone(
    observation: &CohortObservation,
    detail: &SsrDetail,
    verdict: SelectionVerdict,
) -> LocusSelection {
    let alleles = CandidateAlleles::new(
        observation.alleles[0].clone(),
        LocusKind::Ssr(detail.clone()),
    );
    let mut remap = AlleleRemap::with_all_dropped(observation.alleles.len());
    remap.admit(0, AlleleId::REFERENCE);
    let covering_samples = observation.per_sample.len();
    LocusSelection::new(
        alleles,
        verdict,
        vec![super::UnmatchedSupport::default(); covering_samples],
        remap,
        covering_samples,
    )
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

    // ---- C1: which repeat counts one sample puts forward ----

    /// Nominate for one sample of `observation`, at `ploidy` copies and a bar of `floor` reads or
    /// `share` of that sample's own spanning reads.
    ///
    /// It runs the whole chain the shipped path will — fold, ladder, histogram, nomination —
    /// rather than hand-building a histogram, so a test cannot pass against numbers no ladder
    /// produced.
    fn promoted_for(
        observation: &CohortObservation,
        sample: usize,
        ploidy: u8,
        floor: u32,
        share: f64,
    ) -> Vec<u32> {
        let rule = support_rule_of(floor, share);
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, rule, &mut scratch);
        build_ladder(observation, &dinucleotide(), &mut scratch);

        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: rule,
                ..SsrSelectionConfig::at_ploidy(diploid()).shared
            },
            ploidy: Ploidy::try_new(ploidy).expect("at least one copy"),
            ..SsrSelectionConfig::at_ploidy(diploid())
        };
        let SelectionScratch {
            ladder,
            sample_reads_per_rung,
            promoted_rungs,
            ..
        } = &mut scratch;
        let support = &observation.per_sample[sample];
        fill_sample_reads_per_rung(support, observation.region, ladder, sample_reads_per_rung);
        promote_rungs_for_sample(
            sample_reads_per_rung,
            super::super::compared_reads_of(support),
            &config,
            ladder,
            promoted_rungs,
        );
        promoted_rungs.clone()
    }

    /// **The test spec §13 names as the one production cannot pass: a sample with 150 reads at
    /// ten repeats and 150 at eleven nominates both.**
    ///
    /// Production requires a length's reads to exceed *both* neighbours by more than three before
    /// it is nominated at all, so at a heterozygote whose two copies differ by one repeat neither
    /// length is a peak — each has the other beside it — and the sample resolves nothing. Spec
    /// §4.1 measures the cost at 33–78% of such tracts against 97–100%. This asserts the
    /// difference from production rather than a value, so it stays meaningful when the constants
    /// move.
    #[test]
    fn a_sample_with_equal_reads_at_adjacent_lengths_nominates_both() {
        let observation = locus_of(
            &[b"ATATATATATATATATATAT", b"ATATATATATATATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 150, -150.0), row(1, 150, -150.0)],
            )],
        );
        assert_eq!(
            promoted_for(&observation, 0, 2, 2, 0.05),
            vec![0, 1],
            "ten repeats and eleven, both promoted — no neighbour is consulted"
        );
    }

    /// **Only the best `ploidy` rungs survive.** A diploid sample carries two copies, so a third
    /// length is stutter or error however many reads it has relative to the bar.
    #[test]
    fn a_diploid_sample_promotes_its_two_best_supported_rungs() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT", b"ATATATATATAT"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row(1, 30, -30.0),
                    row(2, 40, -40.0),
                    row(3, 6, -6.0),
                ],
            )],
        );
        assert_eq!(
            promoted_for(&observation, 0, 2, 2, 0.0),
            vec![1, 2],
            "thirty reads at four repeats and forty at five, not the four- and six-read rungs — \
             and **in ascending rung order**, where the cut ranked them 2 then 1"
        );
    }

    /// The same locus at three copies promotes three rungs — the ploidy is the run's, and it is
    /// the only thing that changes here.
    #[test]
    fn a_triploid_sample_promotes_three() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT", b"ATATATATATAT"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row(1, 30, -30.0),
                    row(2, 40, -40.0),
                    row(3, 6, -6.0),
                ],
            )],
        );
        assert_eq!(promoted_for(&observation, 0, 3, 2, 0.0), vec![1, 2, 3]);
    }

    /// **A tie between two rungs goes to the shorter repeat count.**
    ///
    /// Both middle rungs carry twenty reads and only one fits under a haploid sample's single
    /// copy. The rule is production's, and it is kept for determinism rather than for a measured
    /// reason: the run's output must be byte-identical at any worker count, so no tie may fall
    /// through to the order the rungs were walked.
    #[test]
    fn a_tie_in_support_goes_to_the_shorter_repeat_count() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0)],
            )],
        );
        assert_eq!(
            promoted_for(&observation, 0, 1, 2, 0.0),
            vec![0],
            "rung 0 is three repeats, rung 1 is four"
        );
    }

    /// A rung whose reads do not reach the bar is not nominated, even where the sample has copies
    /// to spare — the bar and the `ploidy` cut are two different questions.
    #[test]
    fn a_rung_below_the_bar_is_not_nominated_even_with_copies_to_spare() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 9, -9.0), row(1, 1, -1.0)])],
        );
        assert_eq!(
            promoted_for(&observation, 0, 2, 2, 0.0),
            vec![0],
            "one read against a floor of two, and the second copy simply goes unused"
        );
    }

    /// **The share is of this sample's own spanning reads.** Two reads out of ten clear a bar of
    /// 20 in 100; the same two out of a hundred do not.
    ///
    /// The floor is set to 1 so that the share is what decides — at the shipped floor of two reads
    /// both cases would pass on the floor alone and the test would assert nothing.
    #[test]
    fn the_bar_is_a_share_of_the_samples_own_spanning_reads() {
        let shallow = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 8, -8.0), row(1, 2, -2.0)])],
        );
        assert_eq!(
            promoted_for(&shallow, 0, 2, 1, 0.20),
            vec![0, 1],
            "two reads in ten is a fifth of them"
        );

        let deep = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 98, -98.0), row(1, 2, -2.0)])],
        );
        assert_eq!(
            promoted_for(&deep, 0, 2, 1, 0.20),
            vec![0],
            "the same two reads in a hundred are not, and the bar asks for twenty"
        );
    }

    /// **What one sample nominates does not depend on who else is in the run.**
    ///
    /// The same sample, alone and beside a second sample that is deep at a length it never showed.
    /// The denominator is the sample's own spanning reads and no term of the bar reads the cohort,
    /// so a run of one and a run of a thousand give this sample the same answer — the property
    /// that has to hold across the cohort-size range.
    #[test]
    fn a_samples_nomination_is_the_same_alone_and_in_a_cohort() {
        let alone = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 5, -5.0), row(1, 5, -5.0)])],
        );
        let in_cohort = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 5, -5.0), row(1, 5, -5.0)]),
                sample_showing(1, vec![row(2, 400, -400.0)]),
            ],
        );
        assert_eq!(promoted_for(&alone, 0, 2, 2, 0.10), vec![0, 1]);
        assert_eq!(
            promoted_for(&in_cohort, 0, 2, 2, 0.10),
            vec![0, 1],
            "the neighbour's four hundred reads at a third length change nothing here"
        );
    }

    /// A sample whose reads all stopped inside the tract nominates nothing, and does not divide by
    /// its zero denominator on the way.
    #[test]
    fn a_sample_with_only_partials_nominates_nothing() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 5, -5.0), row(1, 5, -5.0)]),
                sample_with_only_partials(1, 40),
            ],
        );
        assert!(promoted_for(&observation, 1, 2, 2, 0.10).is_empty());
    }

    /// **Nominating for a second sample leaves none of the first's rungs**, since the buffer is
    /// one per worker and walked sample by sample.
    ///
    /// The first sample here promotes two rungs and the second promotes one, so a buffer that was
    /// appended to rather than emptied would hand the second sample a length it never showed —
    /// and the union C2 builds would then carry it for the whole cohort.
    #[test]
    fn nominating_for_a_second_sample_leaves_none_of_the_firsts_rungs() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 9, -9.0), row(1, 9, -9.0)]),
                sample_showing(1, vec![row(2, 9, -9.0)]),
            ],
        );
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = ladder_over(&observation);
        let SelectionScratch {
            ladder,
            sample_reads_per_rung,
            promoted_rungs,
            ..
        } = &mut scratch;

        for (sample, expected) in [(0_usize, vec![0_u32, 1]), (1, vec![2])] {
            let support = &observation.per_sample[sample];
            fill_sample_reads_per_rung(support, observation.region, ladder, sample_reads_per_rung);
            promote_rungs_for_sample(
                sample_reads_per_rung,
                super::super::compared_reads_of(support),
                &config,
                ladder,
                promoted_rungs,
            );
            assert_eq!(*promoted_rungs, expected);
        }
    }

    /// A histogram of another locus's width is refused rather than nominating lengths this tract
    /// does not hold.
    #[test]
    #[should_panic(expected = "one entry each")]
    fn a_histogram_that_is_not_this_ladders_is_refused() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![sample_showing(0, vec![row(0, 5, -5.0), row(1, 5, -5.0)])],
        );
        let mut scratch = ladder_over(&observation);
        let mut promoted = Vec::new();
        promote_rungs_for_sample(
            &[3, 3, 3],
            9,
            &SsrSelectionConfig::at_ploidy(diploid()),
            &scratch.ladder,
            &mut promoted,
        );
        scratch.ladder.clear();
    }

    // ---- C2: the ±1 rescue, and the cohort's union ----

    /// Nominate for one sample **with the rescue applied**, at `ploidy` copies.
    fn promoted_with_rescue(
        observation: &CohortObservation,
        sample: usize,
        ploidy: u8,
        floor: u32,
        share: f64,
    ) -> Vec<u32> {
        let rule = support_rule_of(floor, share);
        let copies = Ploidy::try_new(ploidy).expect("at least one copy");
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, rule, &mut scratch);
        build_ladder(observation, &dinucleotide(), &mut scratch);

        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: rule,
                ..SsrSelectionConfig::at_ploidy(copies).shared
            },
            ..SsrSelectionConfig::at_ploidy(copies)
        };
        let SelectionScratch {
            ladder,
            sample_reads_per_rung,
            promoted_rungs,
            ..
        } = &mut scratch;
        let support = &observation.per_sample[sample];
        fill_sample_reads_per_rung(support, observation.region, ladder, sample_reads_per_rung);
        promote_rungs_for_sample(
            sample_reads_per_rung,
            super::super::compared_reads_of(support),
            &config,
            ladder,
            promoted_rungs,
        );
        rescue_occupied_neighbours(ladder, config.ploidy, promoted_rungs);
        promoted_rungs.clone()
    }

    /// The cohort's union of promoted rungs, as repeat counts.
    fn cohort_promoted_counts(observation: &CohortObservation, ploidy: u8) -> Vec<u32> {
        let copies = Ploidy::try_new(ploidy).expect("at least one copy");
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(observation, &dinucleotide(), &mut scratch);
        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: support_rule_of(2, 0.0),
                ..SsrSelectionConfig::at_ploidy(copies).shared
            },
            ..SsrSelectionConfig::at_ploidy(copies)
        };
        promote_rungs_for_cohort(observation, &config, &mut scratch);
        (0..scratch.ladder.rung_count())
            .filter(|&rung| scratch.rung_is_promoted[rung])
            .map(|rung| scratch.ladder.repeat_count_at(rung))
            .collect()
    }

    /// **A sample that resolved one length puts forward its occupied neighbours, and not its
    /// unoccupied ones.**
    ///
    /// The sample shows four repeats only. Three repeats is occupied — another sample's reads
    /// reached it — and five is a rung of the table with no reads on it at all, since only the
    /// reference sits there. So the rescue offers three and refuses five: **the occupancy test is
    /// what stops it inventing a length**.
    #[test]
    fn a_sample_resolving_one_length_gains_its_occupied_neighbour_only() {
        let observation = locus_of(
            &[b"ATATATATAT", b"ATATATAT", b"ATATAT"],
            vec![
                sample_showing(0, vec![row(1, 20, -20.0)]),
                sample_showing(1, vec![row(2, 9, -9.0)]),
            ],
        );
        assert_eq!(
            promoted_with_rescue(&observation, 0, 2, 2, 0.0),
            vec![0, 1],
            "rung 0 is three repeats — occupied by the other sample — and rung 2, five repeats, \
             holds only the reference and no reads"
        );
    }

    /// **A sample that resolved its full ploidy gains no neighbour at all.**
    ///
    /// Firing the rescue on a resolved sample would widen every locus by up to two rungs a
    /// sample, which is extra candidates rather than a crash — every one a column in every
    /// genotype table for the life of the locus.
    #[test]
    fn a_sample_resolving_two_lengths_gains_no_neighbour() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 9, -9.0), row(1, 9, -9.0), row(2, 9, -9.0)],
            )],
        );
        // All three lengths clear the bar; the cut keeps two, and rung 2 is a live neighbour of
        // rung 1 — so a rescue that fired regardless of resolution would put it back.
        assert_eq!(
            promoted_with_rescue(&observation, 0, 2, 2, 0.0),
            vec![0, 1],
            "two copies resolved, so five repeats stays cut"
        );
    }

    /// The neighbour is a repeat count away, not a rung away — an unoccupied length between two
    /// occupied ones breaks the chain rather than being stepped over.
    #[test]
    fn the_neighbour_is_one_repeat_away_and_not_one_rung_away() {
        let observation = locus_of(
            &[b"ATATATATATAT", b"ATATATAT", b"ATATAT"],
            vec![
                sample_showing(0, vec![row(1, 20, -20.0)]),
                sample_showing(1, vec![row(0, 9, -9.0), row(2, 9, -9.0)]),
            ],
        );
        // Rungs are 3, 4 and 6 repeats. The resolved length is 4; 3 is occupied and promoted, and
        // 5 does not exist — so rung 2 at six repeats, which is *adjacent as a rung*, is not.
        assert_eq!(promoted_with_rescue(&observation, 0, 2, 2, 0.0), vec![0, 1]);
    }

    /// A length at the bottom of the ladder has no lower neighbour, and asking for one does not
    /// wrap around.
    ///
    /// Zero repeats is a real tract length — a deletion that removed the tract entirely — and it
    /// is the one length whose lower neighbour would underflow. The sample resolves it alone, so
    /// the rescue does fire and reaches upward only.
    #[test]
    fn a_zero_repeat_length_has_no_lower_neighbour() {
        let observation = locus_of(
            &[b"", b"AT"],
            vec![
                sample_showing(0, vec![row(0, 20, -20.0)]),
                sample_showing(1, vec![row(1, 9, -9.0)]),
            ],
        );
        assert_eq!(
            promoted_with_rescue(&observation, 0, 2, 2, 0.0),
            vec![0, 1],
            "zero repeats resolved, one repeat rescued, and nothing below zero asked for"
        );
    }

    /// **The cohort's promoted set is the union across samples** — an allele one accession
    /// carries is still an allele.
    ///
    /// Three samples, each resolving a different single length. A rule that needed two samples to
    /// agree would return one length or none; the union returns all three.
    #[test]
    fn the_cohorts_promoted_set_is_the_union_across_samples() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![
                sample_showing(0, vec![row(0, 20, -20.0)]),
                sample_showing(1, vec![row(1, 20, -20.0)]),
                sample_showing(2, vec![row(2, 20, -20.0)]),
            ],
        );
        assert_eq!(
            cohort_promoted_counts(&observation, 2),
            vec![3, 4, 5],
            "three repeats, four and five — one sample each"
        );
    }

    /// A cohort of one is the same rule as a cohort of many: the union of one sample's rungs is
    /// that sample's rungs.
    #[test]
    fn a_cohort_of_one_promotes_exactly_what_that_sample_does() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0), row(2, 1, -1.0)],
            )],
        );
        assert_eq!(
            cohort_promoted_counts(&observation, 2),
            vec![3, 4],
            "the single read at five repeats reaches no bar, and no sample under-resolved"
        );
    }

    /// A sample that nominates nothing contributes nothing to the union, and does not stop the
    /// samples after it being asked.
    #[test]
    fn a_silent_sample_neither_contributes_nor_blocks() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT"],
            vec![
                sample_with_only_partials(0, 40),
                sample_showing(1, vec![row(0, 9, -9.0), row(1, 9, -9.0)]),
            ],
        );
        assert_eq!(cohort_promoted_counts(&observation, 2), vec![3, 4]);
    }

    /// Building the union for a second locus leaves no flag of the first — the flags are a
    /// per-worker buffer like every other.
    #[test]
    fn a_second_locus_union_leaves_no_flag_of_the_first() {
        let first = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 9, -9.0), row(1, 9, -9.0), row(2, 9, -9.0)],
            )],
        );
        let config = SsrSelectionConfig::at_ploidy(diploid());
        let mut scratch = SelectionScratch::new();
        summarise_alleles(&first, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(&first, &dinucleotide(), &mut scratch);
        promote_rungs_for_cohort(&first, &config, &mut scratch);
        assert_eq!(scratch.rung_is_promoted.len(), 3);

        let second = locus_of(
            &[b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 9, -9.0)])],
        );
        summarise_alleles(&second, support_rule_of(2, 0.0), &mut scratch);
        build_ladder(&second, &dinucleotide(), &mut scratch);
        promote_rungs_for_cohort(&second, &config, &mut scratch);
        assert_eq!(
            scratch.rung_is_promoted,
            vec![true],
            "one rung, and none of the first locus's three left behind"
        );
    }

    // ---- D1: which sequences on a promoted rung are admitted ----

    /// A tract detail with `motif` and short flanks — everything `LocusKind::Ssr` needs.
    fn detail_of(motif: &[u8]) -> SsrDetail {
        SsrDetail {
            motif: Motif::new(motif).expect("a motif within the period range"),
            left_flank: Box::from(b"CCCCC".as_slice()),
            right_flank: Box::from(b"GGGGG".as_slice()),
        }
    }

    /// Run the whole tract narrowing on `observation`: fold, ladder, nomination, admission.
    fn admitted_over(
        observation: &CohortObservation,
        ploidy: u8,
        floor: u32,
        share: f64,
        cap: u16,
    ) -> LocusSelection {
        let rule = support_rule_of(floor, share);
        let copies = Ploidy::try_new(ploidy).expect("at least one copy");
        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: rule,
                max_candidate_alleles: MaxCandidateAlleles::new(cap).expect("a cap of two or more"),
            },
            ..SsrSelectionConfig::at_ploidy(copies)
        };
        let mut scratch = SelectionScratch::new();
        summarise_alleles(observation, rule, &mut scratch);
        build_ladder(observation, &dinucleotide(), &mut scratch);
        promote_rungs_for_cohort(observation, &config, &mut scratch);
        admit_promoted_sequences(observation, &detail_of(b"AT"), &config, &mut scratch)
    }

    /// The surviving sequences, in candidate-id order.
    fn admitted_bases(selection: &LocusSelection) -> Vec<Vec<u8>> {
        selection.alleles().iter().map(<[u8]>::to_vec).collect()
    }

    /// **Both spellings of one promoted length are admitted, each on its own reads** — the class
    /// production cannot reach below three samples, and the reason this rule replaced its
    /// sibling gates.
    ///
    /// One sample, two sequences of five repeats differing by an interior base, plus the
    /// reference at four. Production would promote the rung's best-supported sequence and make
    /// the other clear 8 reads, **3 distinct samples** and a tenth of the rung — the three-sample
    /// term has no cohort-size clamp, so at one sample the second spelling can never be promoted.
    #[test]
    fn both_spellings_of_a_promoted_length_are_admitted() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT", b"ATATATATAG"],
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 9, -9.0), row(2, 9, -9.0)],
            )],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 32);
        assert_eq!(
            admitted_bases(&selection),
            vec![
                b"ATATATAT".to_vec(),
                b"ATATATATAT".to_vec(),
                b"ATATATATAG".to_vec()
            ],
            "the reference, and both five-repeat spellings on their own reads"
        );
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
    }

    /// **A sequence on a rung nobody promoted is not admitted, however many reads it has.**
    ///
    /// The rungs decide which lengths are in play and the bar decides which spellings on them
    /// are; a sequence has to pass both. Here the sample resolves its two copies at three and
    /// four repeats, so the six-repeat sequence stays out even though its reads clear the bar
    /// comfortably.
    #[test]
    fn a_sequence_on_an_unpromoted_rung_is_not_admitted() {
        let observation = locus_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATATAT"],
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0), row(2, 9, -9.0)],
            )],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 32);
        assert_eq!(
            admitted_bases(&selection),
            vec![b"ATATAT".to_vec(), b"ATATATAT".to_vec()],
            "nine reads at six repeats, on a rung two copies had already accounted for"
        );
    }

    /// A sequence on a promoted rung that no sample's reads earned is not admitted — the rung is
    /// a length, not a licence for every spelling at it.
    #[test]
    fn a_sequence_below_the_bar_on_a_promoted_rung_is_not_admitted() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT", b"ATATATATAG"],
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 20, -20.0), row(2, 1, -1.0)],
            )],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 32);
        assert_eq!(
            admitted_bases(&selection),
            vec![b"ATATATAT".to_vec(), b"ATATATATAT".to_vec()],
            "one read against a floor of two, on a rung that was promoted"
        );
    }

    /// **The reference tract is admitted first and is asked nothing** — not the bar, not the
    /// promotion test.
    ///
    /// Here every read sits at five repeats and the reference at four is promoted by nobody and
    /// earned by no one. It is still candidate 0, so the locus yields a table the caller can use.
    #[test]
    fn the_reference_tract_is_admitted_although_no_read_reached_it() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(0, vec![row(1, 20, -20.0)])],
        );
        let selection = admitted_over(&observation, 1, 2, 0.0, 32);
        assert_eq!(
            admitted_bases(&selection),
            vec![b"ATATATAT".to_vec(), b"ATATATATAT".to_vec()]
        );
        assert_eq!(
            selection.remap().candidate_for(0),
            Some(AlleleId::REFERENCE)
        );
    }

    /// Above the cap the list is cut to the best-ranked and the locus is **still called**, with
    /// the count of cut alternatives reported.
    #[test]
    fn above_the_cap_the_worst_ranked_sequences_are_cut_and_the_locus_is_still_called() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT", b"ATATATATAG", b"ATATATATAC"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row(1, 30, -30.0),
                    row(2, 20, -20.0),
                    row(3, 10, -10.0),
                ],
            )],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 3);
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 },
            "three alleles counting the reference, so one of the three spellings is cut"
        );
        assert_eq!(
            admitted_bases(&selection),
            vec![
                b"ATATATAT".to_vec(),
                b"ATATATATAT".to_vec(),
                b"ATATATATAG".to_vec()
            ],
            "the ten-read spelling is the one cut, and the survivors keep the merge's order"
        );
    }

    /// **A sample whose own sequence the cap cut must be emitted as missing**, and the count that
    /// says so is filled here — the one part of the shared leftover this path's likelihood reads.
    ///
    /// Three samples, each homozygous for a different spelling of five repeats, at 60, 40 and 10
    /// reads. **All three are homozygous on purpose**: the cap's first ranking key is the largest
    /// share of one sample's reads, so with a heterozygote in the fixture the ranking turns on
    /// that share rather than on the read totals, and which sample loses its allele stops being
    /// obvious from the numbers. Homozygous everywhere, all three shares are 1.0 and the cohort
    /// read total decides — so the ten-read sample is the one cut, and the one flagged.
    #[test]
    fn a_sample_whose_earned_sequence_the_cap_cut_is_marked_uncallable() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT", b"ATATATATAG", b"ATATATATAC"],
            vec![
                sample_showing(0, vec![row(1, 60, -60.0)]),
                sample_showing(1, vec![row(2, 40, -40.0)]),
                sample_showing(2, vec![row(3, 10, -10.0)]),
            ],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 3);
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        let leftovers = selection.unmatched();
        assert!(!leftovers[0].genotype_must_be_missing());
        assert!(!leftovers[1].genotype_must_be_missing());
        assert!(
            leftovers[2].genotype_must_be_missing(),
            "the third sample's own ten reads were on the sequence the cap cut"
        );
    }

    /// The candidate table is the tract's own kind, carrying the motif and flanks the read model
    /// needs — not `Generic`, which would send the locus down the SNP/indel scoring arm.
    #[test]
    fn the_candidate_table_carries_the_tracts_kind() {
        let observation = locus_of(
            &[b"ATATATAT", b"ATATATATAT"],
            vec![sample_showing(0, vec![row(0, 9, -9.0), row(1, 9, -9.0)])],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 32);
        match selection.alleles().kind() {
            LocusKind::Ssr(detail) => assert_eq!(detail.motif.as_bytes(), b"AT"),
            other => panic!("a repeat tract's table must be Ssr, not {other:?}"),
        }
    }

    /// A tract where every rung was promoted and every sequence earned still admits them in the
    /// **merge table's** order, which is the order that reaches the VCF's ALT column.
    #[test]
    fn the_admitted_order_is_the_merge_tables_and_not_the_ladders() {
        // Merge order is 6, 3, 4 repeats; the ladder's order is 3, 4, 6.
        let observation = locus_of(
            &[b"ATATATATATAT", b"ATATAT", b"ATATATAT"],
            vec![
                sample_showing(0, vec![row(1, 20, -20.0), row(2, 20, -20.0)]),
                sample_showing(1, vec![row(0, 20, -20.0)]),
            ],
        );
        let selection = admitted_over(&observation, 2, 2, 0.0, 32);
        assert_eq!(
            admitted_bases(&selection),
            vec![
                b"ATATATATATAT".to_vec(),
                b"ATATAT".to_vec(),
                b"ATATATAT".to_vec()
            ],
            "the reference first, then the merge's indices 1 and 2 in that order"
        );
    }

    // ---- D2: the periodicity verdict, and the entry point ----

    /// A tract locus of `alleles` covered by `per_sample`, carrying `motif`.
    fn tract_of(
        alleles: &[&[u8]],
        motif: &[u8],
        per_sample: Vec<SampleSupport>,
    ) -> CohortObservation {
        let mut observation = locus_of(alleles, per_sample);
        observation.kind = LocusKind::Ssr(detail_of(motif));
        observation
    }

    /// Narrow `observation` through the entry point, at the shipped defaults but with the bar
    /// spelled out so a handful of reads can decide it.
    fn selected(observation: &CohortObservation, ploidy: u8) -> LocusSelection {
        let rule = support_rule_of(2, 0.0);
        let copies = Ploidy::try_new(ploidy).expect("at least one copy");
        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: rule,
                ..SsrSelectionConfig::at_ploidy(copies).shared
            },
            ..SsrSelectionConfig::at_ploidy(copies)
        };
        select_ssr(observation, &config, &mut SelectionScratch::new())
    }

    /// **A tract whose reads are all a whole number of units from the reference is periodic,
    /// whatever their lengths in bases.**
    ///
    /// The reference here is 49 bases — an `AT` repeat with one extra base in it, which is what a
    /// length-changing interruption leaves behind and what the catalog admits at purity 0.816.
    /// Every read sits at 49, 51 or 47 bases, all an even number of bases from the reference and
    /// none of them an even number of bases from zero. **Anchored on zero this locus is entirely
    /// off the grid and refused**; anchored on the reference it is periodic, which is the whole of
    /// the owner's decision of 2026-09-02.
    #[test]
    fn a_tract_at_an_odd_reference_length_is_periodic_about_its_own_reference() {
        let reference = [b"AT".repeat(20), b"G".to_vec(), b"AT".repeat(4)].concat();
        let shorter = [b"AT".repeat(19), b"G".to_vec(), b"AT".repeat(4)].concat();
        let longer = [b"AT".repeat(21), b"G".to_vec(), b"AT".repeat(4)].concat();
        assert_eq!(reference.len(), 49, "the measured tract, 49 bases");
        let observation = tract_of(
            &[&reference, &shorter, &longer],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0), row(2, 20, -20.0)],
            )],
        );
        assert!(locus_is_periodic(
            &observation,
            &Motif::new(b"AT").expect("a two-base motif"),
            MaxOffGridShare::DEFAULT
        ));
        assert_ne!(
            selected(&observation, 2).verdict(),
            SelectionVerdict::NotPeriodic
        );
    }

    /// A tract whose reads mostly sit at lengths the motif cannot explain is refused, and the
    /// refusal still yields a table: the reference tract alone.
    #[test]
    fn a_tract_whose_reads_are_off_the_motif_grid_is_refused_but_still_yields_the_reference() {
        let observation = tract_of(
            &[b"ATATATATAT", b"ATATATATATA", b"ATATATATATAAA"],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 1, -1.0), row(1, 20, -20.0), row(2, 20, -20.0)],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(selection.verdict(), SelectionVerdict::NotPeriodic);
        assert_eq!(admitted_bases(&selection), vec![b"ATATATATAT".to_vec()]);
        assert_eq!(
            selection.unmatched().len(),
            1,
            "one zeroed leftover per covering sample"
        );
        assert!(!selection.unmatched()[0].genotype_must_be_missing());
    }

    /// **One periodic sample saves a locus every other sample fails.**
    ///
    /// Two samples of forty reads each: the first is entirely off the grid, the second entirely
    /// on it. The verdict is the cohort's, and it takes one sample to carry it — the same shape as
    /// the support bar, and for the same reason.
    #[test]
    fn one_periodic_sample_saves_a_locus_every_other_sample_fails() {
        let alleles: [&[u8]; 3] = [b"ATATATATAT", b"ATATATATATA", b"ATATATATATAT"];
        let all_off = tract_of(
            &alleles,
            b"AT",
            vec![sample_showing(0, vec![row(1, 40, -40.0)])],
        );
        assert_eq!(
            selected(&all_off, 2).verdict(),
            SelectionVerdict::NotPeriodic
        );

        let one_on = tract_of(
            &alleles,
            b"AT",
            vec![
                sample_showing(0, vec![row(1, 40, -40.0)]),
                sample_showing(1, vec![row(2, 40, -40.0)]),
            ],
        );
        assert_ne!(
            selected(&one_on, 2).verdict(),
            SelectionVerdict::NotPeriodic,
            "the second sample is wholly on the grid, and that is enough"
        );
    }

    /// **One read in ten off the grid is allowed; more is not** — the share, at the boundary.
    ///
    /// Forty reads a sample: four off the grid is exactly a tenth and passes, five does not.
    #[test]
    fn the_off_grid_share_decides_at_one_read_in_ten() {
        let alleles: [&[u8]; 2] = [b"ATATATATAT", b"ATATATATATA"];
        let at_the_share = tract_of(
            &alleles,
            b"AT",
            vec![sample_showing(0, vec![row(0, 36, -36.0), row(1, 4, -4.0)])],
        );
        assert!(locus_is_periodic(
            &at_the_share,
            &Motif::new(b"AT").expect("a two-base motif"),
            MaxOffGridShare::DEFAULT
        ));

        let above_it = tract_of(
            &alleles,
            b"AT",
            vec![sample_showing(0, vec![row(0, 35, -35.0), row(1, 5, -5.0)])],
        );
        assert!(!locus_is_periodic(
            &above_it,
            &Motif::new(b"AT").expect("a two-base motif"),
            MaxOffGridShare::DEFAULT
        ));
    }

    /// **A homopolymer can never be off the grid**, because every difference is a whole number of
    /// one-base units. Production short-circuits period 1 explicitly; here it falls out of the
    /// arithmetic, so this test is what says the arithmetic really does it.
    #[test]
    fn a_homopolymer_is_always_periodic() {
        let observation = tract_of(
            &[b"AAAAA", b"AAAAAA", b"AAA"],
            b"A",
            vec![sample_showing(
                0,
                vec![row(0, 1, -1.0), row(1, 20, -20.0), row(2, 20, -20.0)],
            )],
        );
        assert!(locus_is_periodic(
            &observation,
            &Motif::new(b"A").expect("a one-base motif"),
            MaxOffGridShare::DEFAULT
        ));
    }

    /// **A sample with no spanning reads is not asked**, rather than counted as periodic.
    ///
    /// A tract too long for a read to span produces exactly such a sample, and counting it as
    /// periodic would make one of them enough to save every locus in the run — the verdict would
    /// become unreachable wherever coverage is thin, which is where it is most needed.
    #[test]
    fn a_sample_with_no_spanning_reads_does_not_vote() {
        let observation = tract_of(
            &[b"ATATATATAT", b"ATATATATATA"],
            b"AT",
            vec![
                sample_with_only_partials(0, 40),
                sample_showing(1, vec![row(1, 40, -40.0)]),
            ],
        );
        assert!(
            !locus_is_periodic(
                &observation,
                &Motif::new(b"AT").expect("a two-base motif"),
                MaxOffGridShare::DEFAULT
            ),
            "the silent sample does not carry the locus for the one that is off the grid"
        );
    }

    /// A locus **no** sample spanned is periodic by default: there is nothing to judge it on, and
    /// refusing it would be a verdict about coverage rather than about the tract.
    #[test]
    fn a_locus_no_sample_spanned_is_periodic_by_default() {
        let observation = tract_of(
            &[b"ATATATATAT"],
            b"AT",
            vec![sample_with_only_partials(0, 40)],
        );
        assert!(locus_is_periodic(
            &observation,
            &Motif::new(b"AT").expect("a two-base motif"),
            MaxOffGridShare::DEFAULT
        ));
    }

    /// The entry point refuses a locus of the wrong kind rather than scoring a SNP against a
    /// stutter model.
    #[test]
    #[should_panic(expected = "the repeat-tract path was handed a")]
    fn the_entry_point_refuses_a_locus_that_is_not_a_tract() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 9, -9.0), row(1, 9, -9.0)])],
        );
        selected(&observation, 2);
    }

    /// End to end through the entry point: a heterozygote whose two copies are the same length
    /// spelled differently comes back with both spellings and the reference.
    #[test]
    fn the_entry_point_offers_both_spellings_of_one_length() {
        let observation = tract_of(
            &[b"ATATATAT", b"ATATATATAT", b"ATATATATAG"],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 9, -9.0), row(2, 9, -9.0)],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(
            admitted_bases(&selection),
            vec![
                b"ATATATAT".to_vec(),
                b"ATATATATAT".to_vec(),
                b"ATATATATAG".to_vec()
            ]
        );
    }
}
