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
//! **Built in the order `doc/devel/ng/impl_plan/candidate_alleles_ssr.md` set.**
//! This file holds the ladder (its step B1), the path's configuration and the per-sample length
//! histogram (B2), nomination — each sample's own repeat counts, the `±1` rescue and the cohort's
//! union (C1, C2) — the admission of the sequences on a promoted rung (D1), the entry point
//! [`select_ssr`] (D2, whose periodicity verdict was removed on 2026-10-22, spec §7), and the
//! repeat counts that entry point hands the genotype prior (D3). **`select_ssr` is complete** and the run driver dispatches to it
//! (`src/run/callers.rs`). The STR loop plan's Milestone E added one more pass to it: the
//! **discovery pre-pass** ([`nominate_discovered_sequences`]), off by default, which admits the
//! tract sequences one sample showed too often for the model to book as slippage — wired here
//! rather than around the calling loop because the eligible set is a function of the
//! observations and the candidate table alone
//! (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5).
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
use crate::calling::CandidateAlleles;
use crate::calling::likelihood::ssr_emission::respelled_over;
use crate::locus_generation::{FlankDeletion, LocusKind, SsrBundleDetail, SsrDetail};
use crate::region_typing::run_parse::{RunParse, parse_into_runs};
use crate::run::cohort_merge::MinAltReads;
use crate::run::cohort_merge::build::{CohortObservation, SampleSupport};
use crate::types::{AlleleId, GenomeRegion, Motif, Ploidy};

/// **One locus's observed tract sequences, grouped by how many repeats they carry.**
///
/// The key is `bases.len() / motif.period()`, floored — **the same integer the genotype prior
/// indexes its allele length variation by** (`doc/devel/ng/spec/candidate_alleles_ssr.md` §3), which is
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
    /// **What no longer reads it.** The periodicity test measured each read's offset from this
    /// until 2026-09-02, then from the reference tract's length, and was removed on 2026-10-22
    /// (spec §7). The genotype prior does not either: its seed was re-indexed on
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
/// counted twice. The ladder's job is to place every observed sequence somewhere.
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
    build_ladder_keyed(
        observation,
        |table_index| repeat_count_of_bases(&observation.alleles[table_index as usize], motif),
        scratch,
    );
}

/// [`build_ladder`] with each allele's rung chosen by `rung_key_of` rather than by its repeat
/// count — the one change an STR bundle needs, whose rung is a vector of run counts
/// (`str_bundle_loci.md` §7.2). The bundle path numbers its distinct vectors in ascending order
/// and hands those numbers here; everything below is the segment's grouping, unchanged.
pub(super) fn build_ladder_keyed(
    observation: &CohortObservation,
    rung_key_of: impl Fn(u32) -> u32,
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
        rung_earned_by_reads: _,
        cap_cut_table_indices: _,
        discovery_nominated: _,
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

    let repeat_count_of = rung_key_of;

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

/// **The repeat-tract path's settings** — the shared support rule and cap, and the copy number a
/// sample is nominated at (arch §2.2).
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
    /// Copies per genome — **how many rungs a sample promotes** (arch §3.1). From the run's
    /// parameters, never a constant in this module: production hard-asserts diploid at its own
    /// genotyping step and ng does not, and a polyploid run changes only this number.
    pub ploidy: Ploidy,
    /// **The discovery pre-pass's evidence bar, or `None` — the shipped default — for no
    /// pre-pass at all.** `None` is *off*, not *no bar*: with it absent
    /// [`nominate_discovered_sequences`] never runs and the selection is byte-identical to what
    /// it was before discovery existed.
    ///
    /// When set, [`select_ssr`] runs the discovery decision of
    /// `inference::discovery` as a pre-pass over the
    /// merge's rows: a sequence one sample showed with reads clearing this bar is admitted as a
    /// candidate even where nomination's top-`ploidy` cut left its length behind. The bar's
    /// numbers are discovery's own ([`DEFAULT_DISCOVERY_BAR`](crate::calling::inference::DEFAULT_DISCOVERY_BAR),
    /// 2 reads **and** 15% of the sample's tract-spanning reads), handed down from the run's
    /// [`DiscoveryConfig`](crate::calling::inference::DiscoveryConfig) by the driver —
    /// never a constant here, so a sweep of discovery's bar cannot fork into two spellings.
    ///
    /// **A pre-pass and not a round around the loop**, which is the E1 finding
    /// (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5): the tract locus
    /// generator realigns every read before the caller sees it, so the eligible set is a
    /// function of the observations and the candidate table alone, and a second look with the
    /// first's admissions on the table finds nothing.
    pub discovery: Option<MinAltReads>,
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
            ploidy,
            // Off. The driver switches it on from the run's `DiscoveryConfig`
            // (`NG_TRACT_DISCOVERY=1` today), and off must stay byte-identical to the
            // selection that shipped before discovery was wired.
            discovery: None,
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
/// **Only the support rows are counted**: the merge's crossing reads, and at an STR locus the
/// partial reads that single out one sequence, which selection adds to those rows before it
/// starts (spec §8.1). Any other partial read names no length and is scored on its own axis
/// (spec §1.3, §8). So the histogram's total is exactly
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
    rescue_occupied_neighbours_by(
        ladder,
        ploidy,
        promoted,
        &RungNeighbours::OneRepeatEitherSide,
    );
}

/// **Which rung keys are one unit from a rung key** — the `±1` of the rescue
/// ([`rescue_occupied_neighbours`]). On a segment a key is a repeat count and its neighbours are
/// one repeat either side. On an STR bundle a key numbers a vector of run counts, and its
/// neighbours are the vectors one copy more or one copy fewer **in each run separately**
/// (`str_bundle_loci.md` §7.2).
pub(super) enum RungNeighbours<'a> {
    /// A segment's: the repeat count one less and one more.
    OneRepeatEitherSide,
    /// A bundle's: the distinct run-count vectors of the locus's alleles, ascending, whose
    /// positions are the rung keys.
    OneCopyEitherSideInEachRun(&'a [Box<[u32]>]),
}

impl RungNeighbours<'_> {
    /// Every key one unit from `key`, written into `out` (cleared first). A key no allele of the
    /// locus carries is not written: the ladder has no rung for it, so the rescue could not
    /// offer it anyway.
    fn fill(&self, key: u32, out: &mut Vec<u32>) {
        out.clear();
        match self {
            Self::OneRepeatEitherSide => {
                out.extend(
                    [key.checked_sub(1), key.checked_add(1)]
                        .into_iter()
                        .flatten(),
                );
            }
            Self::OneCopyEitherSideInEachRun(vectors) => {
                let Some(counts) = vectors.get(key as usize) else {
                    return;
                };
                let mut neighbour = counts.to_vec();
                for run in 0..counts.len() {
                    for changed in [counts[run].checked_sub(1), counts[run].checked_add(1)]
                        .into_iter()
                        .flatten()
                    {
                        neighbour[run] = changed;
                        if let Ok(at) =
                            vectors.binary_search_by(|vector| (**vector).cmp(&neighbour))
                        {
                            out.push(at as u32);
                        }
                    }
                    neighbour[run] = counts[run];
                }
            }
        }
    }
}

/// [`rescue_occupied_neighbours`], with the neighbours of a rung given by `neighbours`.
pub(super) fn rescue_occupied_neighbours_by(
    ladder: &RepeatLadder,
    ploidy: Ploidy,
    promoted: &mut Vec<u32>,
    neighbours: &RungNeighbours<'_>,
) {
    if promoted.len() >= usize::from(ploidy.get()) {
        return;
    }
    let resolved = promoted.len();
    let mut keys = Vec::new();
    for index in 0..resolved {
        let repeat_count = ladder.repeat_count_at(promoted[index] as usize);
        neighbours.fill(repeat_count, &mut keys);
        for &neighbour in &keys {
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
    promote_rungs_for_cohort_by(
        observation,
        config,
        scratch,
        &RungNeighbours::OneRepeatEitherSide,
    );
}

/// [`promote_rungs_for_cohort`], with the rescue's neighbours given by `neighbours`.
pub(super) fn promote_rungs_for_cohort_by(
    observation: &CohortObservation,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
    neighbours: &RungNeighbours<'_>,
) {
    let SelectionScratch {
        ladder,
        sample_reads_per_rung,
        promoted_rungs,
        rung_is_promoted,
        rung_earned_by_reads,
        ..
    } = scratch;
    rung_is_promoted.clear();
    rung_is_promoted.resize(ladder.rung_count(), false);
    rung_earned_by_reads.clear();
    rung_earned_by_reads.resize(ladder.rung_count(), false);

    for sample in &observation.per_sample {
        fill_sample_reads_per_rung(sample, observation.region, ladder, sample_reads_per_rung);
        promote_rungs_for_sample(
            sample_reads_per_rung,
            super::compared_reads_of(sample),
            config,
            ladder,
            promoted_rungs,
        );
        // Before the rescue widens the list: these are the lengths this sample's own reads put
        // forward, which admission may keep a spelling of even when no single spelling cleared
        // the bar.
        for &rung in promoted_rungs.iter() {
            rung_earned_by_reads[rung as usize] = true;
        }
        rescue_occupied_neighbours_by(ladder, config.ploidy, promoted_rungs, neighbours);
        for &rung in promoted_rungs.iter() {
            rung_is_promoted[rung as usize] = true;
        }
    }
}

/// **Which of the merge's sequences one sample showed too often for slippage to explain** — the
/// discovery decision of [`inference::discovery`](crate::calling::inference::discovery), run
/// as a pre-pass of selection over the merge's own rows.
///
/// `nominated` comes back one flag per allele of the merge's table, true where **some** sample's
/// reads on that sequence cleared `bar` against that sample's own spanning reads. Admission ORs
/// these flags into its filter, so a nominated sequence is admitted even where the top-`ploidy`
/// cut left its rung unpromoted — which is exactly the sequence discovery exists for: a sample's
/// third length, invisible to nomination by construction because
/// [`promote_rungs_for_sample`] truncates to `ploidy`.
///
/// # This is the E1 decision, on selection's data
///
/// [`discover_tract_alleles`](crate::calling::inference::discovery::discover_tract_alleles)
/// walks a sample's realigned observations; at selection those same spellings **are** the
/// merge's alleles, and a sample's complete observations are its support rows (a partial names
/// no length and has no row — spec §8). So the decision's three clauses carry over unchanged,
/// and each is asked through the one producer this crate has for it:
///
/// - **the bar**: [`MinAltReads::reached_by`], reads pooled across read groups by
///   [`one_run_per_allele`] against the denominator [`compared_reads_of`] — the sample's
///   support rows, so a read that ran out inside the STR segment counts in both halves if it
///   singles out one sequence (spec §8.1) and in neither otherwise;
/// - **whole repeats**: a sequence shorter than one copy of the unit
///   ([`repeat_count_of_bases`] of zero) is refused here, where the reason is — admitted, it
///   would reach the driver's conversion to the evidence's non-zero counts
///   (`repeat_counts_the_tract_model_can_take`, `src/run/callers.rs`) and stop the **whole
///   locus** one step later, a strictly worse outcome than leaving the sequence out;
/// - **not already a candidate**: not asked, deliberately. The decision's third clause exists
///   so a round cannot re-admit what the table holds; here admission is a union, so flagging a
///   sequence the ordinary filter admits anyway changes nothing — and it is *why* a second
///   pre-pass over the same evidence admits nothing, the E1 property the tests hold through
///   this wiring.
///
/// **Per sample, never summed across the cohort**: one sample must clear the bar alone, so
/// three samples showing a sequence once each nominate nothing (E1's
/// `a_cohort_sum_does_not_clear_a_bar_no_single_sample_clears`, held here by
/// `discovery_does_not_pool_samples_against_the_bar`).
///
/// # Panics
///
/// On a support row naming an allele the merge's table does not hold, and on a sample's rows
/// out of ascending allele order — the same merge bugs every other walk in this module holds
/// as release assertions (spec §8).
pub(super) fn nominate_discovered_sequences(
    observation: &CohortObservation,
    motif: &Motif,
    bar: MinAltReads,
    nominated: &mut Vec<bool>,
) {
    nominated.clear();
    nominated.resize(observation.alleles.len(), false);
    for sample in &observation.per_sample {
        // The sample's spanning reads alone — a partial names no length, so the merge holds it
        // off the support rows and it reaches neither half of the bar.
        let spanning = super::compared_reads_of(sample);
        for rows in super::one_run_per_allele(sample, observation.region) {
            let allele = rows[0].allele;
            assert!(
                allele < nominated.len(),
                "sample {}'s support row named allele {allele} of a repeat tract whose table \
                 holds {} at {}",
                sample.sample,
                nominated.len(),
                observation.region
            );
            if nominated[allele] {
                continue;
            }
            let pooled_reads = rows.iter().fold(0_u32, |total, row| {
                total.saturating_add(row.support.num_reads)
            });
            if !bar.reached_by(pooled_reads, spanning) {
                continue;
            }
            // No rung below the ladder's first: the stutter model is written in whole repeats,
            // so a sequence carrying none is not offered at all rather than offered and refused
            // downstream — where the refusal would take the whole locus with it.
            if repeat_count_of_bases(&observation.alleles[allele], motif) == 0 {
                continue;
            }
            nominated[allele] = true;
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
/// **A sequence the discovery pre-pass nominated enters on discovery's bar instead**, promoted
/// rung or none — that is the pre-pass's whole point, since the length it finds is exactly the
/// one the top-`ploidy` cut left behind ([`nominate_discovered_sequences`]). Past this filter it
/// is an alternative like any other: ranked by [`compare_best_first`](super::compare_best_first),
/// cut by the cap when it ranks last, and booked in the leftover like the rest when it is. With
/// discovery off — the flags empty, the shipped default — this clause is dead and the filter is
/// byte-identical to what it was before discovery existed.
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
/// # A length the reads earned keeps one spelling (owner, 2026-10-10)
///
/// **Where some sample's own reads earned a repeat length but sequencing errors split them across
/// spellings so that none cleared the bar, the length's best-supported spelling is admitted** —
/// the most cohort reads; among equals, the one closest to a pure run of the motif
/// ([`bases_off_the_motif`]), because when the reads are split one apiece the spellings that
/// differ from the pure run are most likely the ones carrying an error, and the admitted spelling
/// is what the record writes; then the lexically smaller, for determinism. Without it the length is lost
/// although the nomination said it was real: on the stuttering benchmark at one sample and 15×, a
/// 5-base repeat whose sample had lost two copies showed it on three reads, two of them with one
/// wrong base each, so three spellings of one read apiece and no candidate. The reads carrying an
/// error still count for the admitted spelling, through the substitution term that scores every
/// read against every candidate of its length (`read_likelihoods.md` §4.3).
///
/// Two limits keep it to that case: a length the `±1` rescue added — which no sample's reads
/// earned — gets no such spelling; and a length the reference tract or an admitted spelling
/// already stands for gets nothing more. A spelling with no whole repeat is never chosen, since
/// the tract model cannot score it.
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
    merge_rows: &[SampleSupport],
    detail: &SsrDetail,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> LocusSelection {
    admit_promoted_sequences_on(
        observation,
        merge_rows,
        &SelectionGround::Segment(detail),
        config,
        scratch,
    )
}

/// **The repeat ground a selection is made over** — an STR segment, or an STR bundle with the
/// run parse of every allele of the merge's table (`str_bundle_loci.md` §7.2). It answers the
/// three questions admission asks of the ground, each of which a segment answers from its one
/// motif and a bundle from its runs.
pub(super) enum SelectionGround<'a> {
    /// An STR segment.
    Segment(&'a SsrDetail),
    /// An STR bundle, and each merge allele's run parse, in the merge table's order.
    Bundle {
        detail: &'a SsrBundleDetail,
        parses: &'a [RunParse],
    },
}

impl SelectionGround<'_> {
    /// Whether the read model can score this allele at all: a segment's must hold a whole copy
    /// of the motif; any bundle allele can be scored, a run at no copies simply cannot slip.
    fn can_score(&self, observation: &CohortObservation, table_index: usize) -> bool {
        match self {
            Self::Segment(detail) => {
                repeat_count_of_bases(&observation.alleles[table_index], &detail.motif) > 0
            }
            Self::Bundle { .. } => true,
        }
    }

    /// How far this allele is from pure repeat — bases off the motif on a segment, the run
    /// parse's edited bases on a bundle — the tie-break among spellings of one rung.
    fn bases_off_the_repeat(&self, observation: &CohortObservation, table_index: usize) -> usize {
        match self {
            Self::Segment(detail) => {
                bases_off_the_motif(&observation.alleles[table_index], &detail.motif)
            }
            Self::Bundle { parses, .. } => parses[table_index].edited_bases() as usize,
        }
    }

    /// The reference flanks either side of the ground.
    fn flanks(&self) -> (&[u8], &[u8]) {
        match self {
            Self::Segment(detail) => (&detail.left_flank, &detail.right_flank),
            Self::Bundle { detail, .. } => (&detail.left_flank, &detail.right_flank),
        }
    }

    /// The kind the candidate table is stamped with.
    fn locus_kind(&self) -> LocusKind {
        match self {
            Self::Segment(detail) => LocusKind::Ssr((*detail).clone()),
            Self::Bundle { detail, .. } => LocusKind::SsrBundle((*detail).clone()),
        }
    }
}

/// [`admit_promoted_sequences`] over either kind of repeat ground.
///
/// `observation` is what selection reads — the merge's locus with the deciding partial reads
/// counted ([`count_deciding_partial_reads`](super::deciding_partial_reads::count_deciding_partial_reads))
/// — and `merge_rows` is the merge's own per-sample rows, parallel to `observation.per_sample`,
/// which the leftover's read count and error mass are taken from ([`leftover_of`](super::leftover_of)).
pub(super) fn admit_promoted_sequences_on(
    observation: &CohortObservation,
    merge_rows: &[SampleSupport],
    ground: &SelectionGround<'_>,
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
        rung_earned_by_reads,
        cap_cut_table_indices,
        discovery_nominated,
        ..
    } = scratch;
    cap_cut_table_indices.clear();
    assert_eq!(
        rung_is_promoted.len(),
        ladder.rung_count(),
        "the promotion flags run parallel to the ladder's rungs, one each, at {}",
        observation.region
    );
    // Empty is *discovery off* — the shipped default, under which this function must admit
    // byte-identically to what it did before discovery was wired. Anything else must be one
    // flag per merge allele; a buffer of another locus's width would nominate by a
    // neighbour's evidence, silently.
    assert!(
        discovery_nominated.is_empty() || discovery_nominated.len() == observation.alleles.len(),
        "the discovery flags hold {} entries over a merge table of {} at {}: they are either \
         empty (discovery off) or one per allele, and any other width is another locus's buffer",
        discovery_nominated.len(),
        observation.alleles.len(),
        observation.region
    );

    // Every alternative that sits on a promoted rung *and* some sample's reads earned — or that
    // the discovery pre-pass nominated on its own bar, which is how a length the top-`ploidy`
    // cut hid still reaches the table — in the merge table's own order for now; admission sorts
    // it by bases.
    ranked_table_indices.clear();
    ranked_table_indices.extend(
        (1..observation.alleles.len())
            .filter(|&index| {
                (rung_is_promoted[ladder.rung_of_table_index(index)]
                    && per_allele[index].cleared_the_bar())
                    || (!discovery_nominated.is_empty() && discovery_nominated[index])
            })
            // **A sequence with no whole repeat is left out here, not refused downstream**
            // (owner, 2026-10-10): the tract model cannot score it, and admitted it reaches the
            // driver's conversion to non-zero repeat counts and stops the whole locus —
            // `repeat_counts_the_tract_model_can_take`, `src/run/callers.rs`. On the stuttering
            // benchmark that dropped a segment three soft-clipped reads had spelled as empty,
            // and with it the sample's real allele.
            .filter(|&index| ground.can_score(observation, index))
            .map(|index| {
                u32::try_from(index).expect("a merge table narrower than four billion alleles")
            }),
    );

    // **A length some sample's reads earned, which no admitted spelling and not the reference
    // stands for, keeps its best-supported spelling** — the doc's "a length the reads earned".
    // Empty flags are a selection that never ran nomination (a test driving admission alone);
    // anything else is one flag per rung.
    assert!(
        rung_earned_by_reads.is_empty() || rung_earned_by_reads.len() == ladder.rung_count(),
        "the earned-length flags hold {} entries over a ladder of {} rungs at {}",
        rung_earned_by_reads.len(),
        ladder.rung_count(),
        observation.region
    );
    let reference_rung = ladder.rung_of_table_index(0);
    for rung in (0..rung_earned_by_reads.len()).filter(|&rung| rung_earned_by_reads[rung]) {
        let stood_for = rung == reference_rung
            || ranked_table_indices
                .iter()
                .any(|&index| ladder.rung_of_table_index(index as usize) == rung);
        if stood_for {
            continue;
        }
        let best = (1..observation.alleles.len())
            .filter(|&index| ladder.rung_of_table_index(index) == rung)
            .filter(|&index| ground.can_score(observation, index))
            .max_by(|&left, &right| {
                per_allele[left]
                    .cohort_reads
                    .cmp(&per_allele[right].cohort_reads)
                    .then_with(|| {
                        ground
                            .bases_off_the_repeat(observation, right)
                            .cmp(&ground.bases_off_the_repeat(observation, left))
                    })
                    .then_with(|| observation.alleles[right].cmp(&observation.alleles[left]))
            });
        if let Some(index) = best {
            ranked_table_indices.push(
                u32::try_from(index).expect("a merge table narrower than four billion alleles"),
            );
        }
    }
    // Nothing below reads the buffer's order: the misread fold keeps or drops each entry on its
    // own, the cap ranks on a total order, and admission sorts by bases.

    // **An allele over a flank extent one misread base away from an allele already standing is
    // not a candidate of its own** (2026-10-11; `str_edge_deletion_alleles.md` §5.5). Respelled
    // over the wider extent, such an allele differs from the reference or an admitted plain
    // allele by a single substituted base — the shape a systematic read error beside a run
    // takes: on HG002 at 30×, four reads in 36 showing the `G` after 25 `T`s misread, the same
    // way, earned a candidate and a false heterozygote. Its reads still count for the allele
    // they are one base from, at the substitution rate. **Only where, in every sample showing
    // it, it is a small minority beside that allele** — under [`MISREAD_SHARE_CEILING`] of the
    // two together: a misread base takes a share of the reads of the allele it came from, while
    // a carrier of a real allele shows it at about half its reads or all of them. Folding every
    // allele one base from another removed the commonest real edge deletion on the stuttering
    // benchmark — one flank base, which once written out is one base from a plain contraction —
    // and lost 36 to 77 carrier pairs a run at 50 samples.
    let standing: Vec<usize> = std::iter::once(0)
        .chain(ranked_table_indices.iter().map(|&index| index as usize))
        .filter(|&index| observation.flank_deletion_of(index).is_none())
        .collect();
    ranked_table_indices.retain(|&index| {
        let index = index as usize;
        let extent = observation.flank_deletion_of(index);
        extent.is_none()
            || !standing.iter().any(|&other| {
                one_base_apart(
                    &observation.alleles[index],
                    extent,
                    &observation.alleles[other],
                    observation.flank_deletion_of(other),
                    ground.flanks(),
                ) && a_minority_beside_it_everywhere(observation, index, other)
            })
    });

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
        // Kept, not discarded: on this path an allele is lost three ways — the bar, the
        // per-sample top-`ploidy` cut, and the cap — and only the third may no-call a sample.
        cap_cut_table_indices.extend_from_slice(&ranked_table_indices[allowed_alternatives..]);
        cap_cut_table_indices.sort_unstable();
        ranked_table_indices.truncate(allowed_alternatives);
        SelectionVerdict::Truncated {
            dropped: u32::try_from(dropped).expect("a merge table narrower than four billion"),
        }
    };
    // **Admitted in the order of their bases, not the merge table's**, for the generic path's
    // reason (`generic.rs`): the table's order is the order the samples were given, and two
    // spellings of one length swapped `ALT` places when the cohort was listed the other way.
    ranked_table_indices.sort_unstable_by(|&left, &right| {
        observation.alleles[left as usize].cmp(&observation.alleles[right as usize])
    });

    let mut alleles = CandidateAlleles::new(observation.alleles[0].clone(), ground.locus_kind());
    let mut remap = AlleleRemap::with_all_dropped(observation.alleles.len());
    remap.admit(0, AlleleId::REFERENCE);
    for &table_index in ranked_table_indices.iter() {
        let table_index = table_index as usize;
        let candidate = alleles.admit_with_flank_deletion(
            observation.alleles[table_index].clone(),
            observation.flank_deletion_of(table_index),
        );
        remap.admit(table_index, candidate);
    }

    assert_eq!(
        merge_rows.len(),
        observation.per_sample.len(),
        "the merge's rows and the rows selection read cover different samples at {}",
        observation.region
    );
    let covering_samples = observation.per_sample.len();
    let leftovers = merge_rows
        .iter()
        .zip(&observation.per_sample)
        .map(|(sample, selection_rows)| {
            super::leftover_of(
                sample,
                selection_rows,
                observation.region,
                &remap,
                config.shared.min_allele_support,
                cap_cut_table_indices,
            )
        })
        .collect();
    LocusSelection::new(alleles, verdict, leftovers, remap, covering_samples)
}

/// **How many of a tract sequence's bases differ from a pure run of the motif**, at whichever
/// rotation of the motif fits it best — the tie-break that picks, among spellings of one length
/// with equal reads, the one least likely to carry a sequencing error
/// ([`admit_promoted_sequences`]). Zero for a pure run; an interrupted reference tract scores its
/// interruptions, which is fine, because it only ever ranks spellings of one length against each
/// other.
fn bases_off_the_motif(bases: &[u8], motif: &Motif) -> usize {
    let unit = motif.as_bytes();
    let period = unit.len();
    (0..period)
        .map(|rotation| {
            bases
                .iter()
                .enumerate()
                .filter(|&(at, &base)| base != unit[(at + rotation) % period])
                .count()
        })
        .min()
        .unwrap_or(0)
}

/// **How many whole motif copies a tract sequence carries** — `bases.len() / motif.period()`,
/// floored.
///
/// **The one producer of this integer in ng's calling path, and that is the point.** The ladder
/// keys its rungs by it, [`SsrLocusSelection::repeat_counts`] hands it to the genotype prior, and
/// the prior's allele length variation is indexed by it (`doc/devel/ng/spec/candidate_alleles_ssr.md` §3).
/// Two spellings of a floor division is how two modules come to disagree about which rung an
/// allele is on, and a disagreement there puts a candidate's prior mass on the wrong length with
/// nothing failing.
///
/// The period is at least one base — [`Motif::new`] refuses an empty unit — so this cannot divide
/// by zero.
///
/// **`pub(crate)` for the discovery round**, which asks the same question of a sequence no
/// candidate holds ([`inference::discovery`](crate::calling::inference::discovery)) and must
/// get the ladder's answer rather than a second floor division.
#[inline]
pub(crate) fn repeat_count_of_bases(bases: &[u8], motif: &Motif) -> u32 {
    (bases.len() / motif.period()) as u32
}

/// **What the repeat-tract path returns: the shared selection, plus the repeat count of every
/// candidate it kept** (arch §2.3).
///
/// # Why the repeat counts travel with the selection rather than being recomputed
///
/// The genotype prior's `fill_ssr_seed` takes exactly this slice. It could derive the same
/// integers from the candidates' bases and the motif — and that is the failure this field exists
/// to prevent: a second floor division is a second producer, and the two only have to disagree
/// once for a candidate's prior mass to land on the wrong length, silently
/// (spec §3). Everything here comes from [`repeat_count_of_bases`].
///
/// # What is *not* here, and why
///
/// **The ladder's modal repeat count is not returned.** Arch §2.3 lists it, and arch §4 and §5
/// then retract its reason: on 2026-08-27 the prior's seed was re-indexed by offset from the
/// **reference** tract length, which every locus already knows, so the cohort's commonest length
/// stopped being an input to it — *"nothing consumes it now"*. **This step's own D2 removed the
/// last remaining consumer inside this module**, by anchoring the periodicity grid on the
/// reference length too (the grid itself went on 2026-10-22, spec §7). Returning it would be a
/// field with no reader, and the architecture's own
/// open question — whether selection should carry the mode at all — is left where it is rather
/// than answered by shipping one. The ladder still computes it, so restoring it is one line.
#[derive(Clone, PartialEq, Debug)]
pub struct SsrLocusSelection {
    /// The narrowed table, the verdict, the per-sample leftover and the remapping — the same
    /// bundle the ordinary path returns.
    pub selection: LocusSelection,
    /// **Parallel to `selection.alleles()`**, the reference at index 0: each surviving
    /// candidate's repeat count.
    pub repeat_counts: Vec<u32>,
}

/// **The largest share of a sample's reads, over a misread allele and the allele it is one base
/// from, at which the first is still taken for a misreading of the second** (2026-10-11). A
/// heterozygous carrier of a real allele shows it near a half; the systematic misreading this is
/// for showed at 4 reads against 13 on HG002, under a quarter. A judgement between the two, not
/// a fitted number.
const MISREAD_SHARE_CEILING: f64 = 1.0 / 3.0;

/// **Is allele `minor` a small minority beside allele `major` in every sample that shows it?**
/// — under [`MISREAD_SHARE_CEILING`] of that sample's reads of the two together.
fn a_minority_beside_it_everywhere(
    observation: &CohortObservation,
    minor: usize,
    major: usize,
) -> bool {
    observation.per_sample.iter().all(|sample| {
        let reads_of = |allele: usize| -> u32 {
            sample
                .supported
                .iter()
                .filter(|row| row.allele == allele)
                .map(|row| row.support.num_reads)
                .sum()
        };
        let (minor_reads, major_reads) = (reads_of(minor), reads_of(major));
        minor_reads == 0
            || f64::from(minor_reads) < MISREAD_SHARE_CEILING * f64::from(minor_reads + major_reads)
    })
}

/// **Do two alleles, respelled over the wider of their flank extents, differ by exactly one
/// substituted base?** — the test that keeps a systematic read error beside a run from standing
/// as an allele of its own (`admit_promoted_sequences`).
fn one_base_apart(
    one: &[u8],
    one_extent: FlankDeletion,
    other: &[u8],
    other_extent: FlankDeletion,
    flanks: (&[u8], &[u8]),
) -> bool {
    let wider = FlankDeletion {
        left: one_extent.left.max(other_extent.left),
        right: one_extent.right.max(other_extent.right),
    };
    let (mut one_over, mut other_over) = (Vec::new(), Vec::new());
    respelled_over(one, one_extent, wider, flanks, &mut one_over).is_some()
        && respelled_over(other, other_extent, wider, flanks, &mut other_over).is_some()
        && one_over.len() == other_over.len()
        && one_over
            .iter()
            .zip(&other_over)
            .filter(|(a, b)| a != b)
            .count()
            == 1
}

/// **Narrow one repeat tract's allele table to the tract sequences worth calling over**
/// (spec §4, §5, §7) — the entry point, and the repeat-tract sibling of
/// [`select_generic`](super::generic::select_generic).
///
/// **Four passes over the locus, in this order**, because each needs the one before: fold the
/// rows into per-allele summaries; build the ladder over them; nominate rungs per sample and union
/// them; then admit the sequences on nominated rungs that cleared the bar, apply the cap and fill
/// the leftover. **Before them**, each partial read that fits exactly one sequence of the table is
/// counted as a read of it, on a copy the four passes read (spec §8.1, the private
/// `deciding_partial_reads` module); the leftover's read count and error mass are still the
/// merge's own rows'.
///
/// **With discovery on there is a fifth, between nomination and admission**: the discovery
/// pre-pass ([`nominate_discovered_sequences`]) flags the sequences some single sample showed
/// too often for slippage to explain, and admission ORs them in — so their reads land in the
/// record's `AD` through the remap rather than in the "no written allele explains these"
/// leftover, and the cap and the verdict treat them like any other alternative. It is a
/// pre-pass and not a round around the calling loop, which is the E1 finding
/// (`doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5): the eligible set is a
/// function of the observations and the candidate table alone, so a second look after these
/// admissions finds nothing. Off — [`SsrSelectionConfig::discovery`] `None`, the shipped
/// default — nothing of it executes.
///
/// # Panics
///
/// On a cohort observation that is not a repeat tract. Which path a locus takes is decided by its
/// kind at the driver's dispatch, so reaching here with a `Generic` or bundle locus is a routing
/// bug, and the alternative — scoring a SNP against a stutter model — is a confident wrong answer
/// rather than a crash.
/// **Public because the measurement calls it.** `examples/ng_candidate_selection_probe.rs`'s
/// repeat-tract arm runs the shipped selector over real reads, so a number it prints is a number
/// the caller itself would produce rather than a re-derivation of one (spec §10). The run
/// driver's dispatch will be the second caller.
pub fn select_ssr(
    observation: &CohortObservation,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> SsrLocusSelection {
    let LocusKind::Ssr(detail) = &observation.kind else {
        panic!(
            "the repeat-tract path was handed a {:?} locus at {}: which path a locus takes is \
             decided by its kind at the dispatch, and scoring one kind against the other's model \
             is a confident wrong answer rather than a failure",
            observation.kind, observation.region
        );
    };

    // Partial reads that single out one sequence count toward the bar (spec §8.1).
    let merge_rows = &observation.per_sample;
    let counted = super::deciding_partial_reads::count_deciding_partial_reads(
        observation,
        (&detail.left_flank, &detail.right_flank),
    );
    let observation = counted.as_ref().unwrap_or(observation);

    summarise_alleles(observation, config.shared.min_allele_support, scratch);
    build_ladder(observation, &detail.motif, scratch);
    promote_rungs_for_cohort(observation, config, scratch);
    if let Some(bar) = config.discovery {
        nominate_discovered_sequences(
            observation,
            &detail.motif,
            bar,
            &mut scratch.discovery_nominated,
        );
    }
    let selection = admit_promoted_sequences(observation, merge_rows, detail, config, scratch);
    let repeat_counts = selection
        .alleles()
        .iter()
        .map(|bases| repeat_count_of_bases(bases, &detail.motif))
        .collect();
    SsrLocusSelection {
        selection,
        repeat_counts,
    }
}

/// **What the STR bundle path returns: the shared selection, plus the run parse of every candidate
/// it kept** (`str_bundle_loci.md` §7.1, §7.2) — the bundle's counterpart of
/// [`SsrLocusSelection`], whose one repeat count per candidate becomes a run count per run.
#[derive(Clone, PartialEq, Debug)]
pub struct SsrBundleSelection {
    /// The narrowed table, the verdict, the per-sample leftover and the remapping.
    pub selection: LocusSelection,
    /// **Parallel to `selection.alleles()`**, the reference at index 0: each candidate read
    /// against the bundle's reference runs — its run counts, and each run's place in it.
    pub parses: Vec<RunParse>,
}

/// **Narrow one STR bundle's allele table to the sequences worth calling over** — [`select_ssr`]
/// with each single-motif rule generalised to the bundle's runs (`str_bundle_loci.md` §7.2):
///
/// - **the bar is unchanged** — at least two crossing reads of one sample, and a tenth of them;
/// - **a rung is a vector of run counts**: every allele of the merge's table is parsed against the
///   reference runs once ([`parse_into_runs`]), and two alleles with the same run counts are two
///   spellings on one rung;
/// - **no periodicity check**: a length change that is no whole number of run units is an indel
///   in the fixed sequence as often as anything (owner, 2026-10-12), so nothing is refused for it;
/// - **the `±1` rescue** offers one copy more and one fewer in each run separately
///   ([`RungNeighbours::OneCopyEitherSideInEachRun`]);
/// - **the tie-break toward the purest spelling** counts the parse's edited bases.
///
/// **The discovery pre-pass is not run on bundles**, which is off by default for segments too; its
/// bar is written in a segment's repeat counts and nothing has measured it on a bundle.
///
/// # Panics
///
/// On a cohort observation that is not an STR bundle, for the reason [`select_ssr`] gives.
pub fn select_ssr_bundle(
    observation: &CohortObservation,
    config: &SsrSelectionConfig,
    scratch: &mut SelectionScratch,
) -> SsrBundleSelection {
    let LocusKind::SsrBundle(detail) = &observation.kind else {
        panic!(
            "the STR bundle path was handed a {:?} locus at {}: which path a locus takes is \
             decided by its kind at the dispatch",
            observation.kind, observation.region
        );
    };
    // Partial reads that single out one sequence count toward the bar (`candidate_alleles_ssr.md`
    // §8.1), at bundles as at segments.
    let merge_rows = &observation.per_sample;
    let counted = super::deciding_partial_reads::count_deciding_partial_reads(
        observation,
        (&detail.left_flank, &detail.right_flank),
    );
    let observation = counted.as_ref().unwrap_or(observation);
    let hull = &observation.alleles[0];
    // No periodicity check (owner, 2026-10-12): the fixed sequence between runs carries real
    // indels of any length, so a change that is no whole number of run units is not a sign the
    // locus is not a repeat.
    let parses: Vec<RunParse> = observation
        .alleles
        .iter()
        .map(|allele| parse_into_runs(allele, hull, &detail.runs))
        .collect();
    // The rung keys: each distinct vector of run counts, numbered in ascending order.
    let mut run_count_vectors: Vec<Box<[u32]>> = parses
        .iter()
        .map(|parse| parse.copies().collect())
        .collect();
    run_count_vectors.sort_unstable();
    run_count_vectors.dedup();
    // Each allele's rung key, its vector's place in that numbering, looked up once: the ladder's
    // sort asks for an allele's key at every comparison it makes.
    let rung_keys: Vec<u32> = parses
        .iter()
        .map(|parse| {
            run_count_vectors
                .binary_search_by(|vector| vector.iter().copied().cmp(parse.copies()))
                .expect("every allele's run counts are among the vectors numbered")
                as u32
        })
        .collect();

    summarise_alleles(observation, config.shared.min_allele_support, scratch);
    build_ladder_keyed(
        observation,
        |table_index| rung_keys[table_index as usize],
        scratch,
    );
    promote_rungs_for_cohort_by(
        observation,
        config,
        scratch,
        &RungNeighbours::OneCopyEitherSideInEachRun(&run_count_vectors),
    );
    let selection = admit_promoted_sequences_on(
        observation,
        merge_rows,
        &SelectionGround::Bundle {
            detail,
            parses: &parses,
        },
        config,
        scratch,
    );
    let parses = selection
        .alleles()
        .iter()
        .map(|bases| parse_into_runs(bases, hull, &detail.runs))
        .collect();
    SsrBundleSelection { selection, parses }
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::*;
    use super::super::summarise_alleles;
    use super::*;
    use crate::calling::inference::DEFAULT_DISCOVERY_BAR;
    use crate::types::ReadGroupId;

    /// **The record's `REPCN` and the genotype prior's repeat count are the same integer.**
    ///
    /// The VCF derives a called allele's repeat copies from its own bases and the motif
    /// ([`TractAnnotation::repeat_copies_of`](crate::vcf::TractAnnotation::repeat_copies_of));
    /// this module derives the same count for the same allele and hands it to the prior. **Two
    /// derivations of one number is how two modules come to disagree about which rung an allele
    /// sits on**, and a disagreement there puts a candidate's prior mass at one length while the
    /// file reports another, with nothing failing. They agree because both are the same floor
    /// division — and this is what holds them to it, across the two cases where a floor division
    /// and a rounded one part company.
    #[test]
    fn the_records_repeat_count_is_the_one_the_prior_is_given() {
        let motif = dinucleotide();
        let annotation = crate::vcf::TractAnnotation::new(motif);
        for bases in [
            b"AT".as_slice(),
            b"ATAT",
            b"ATATATATAT",
            // A part unit at the end, where flooring is what the two have to share.
            b"ATATA",
            // And shorter than one whole unit, which both count as none.
            b"A",
            b"",
        ] {
            assert_eq!(
                annotation.repeat_copies_of(bases) as u32,
                repeat_count_of_bases(bases, &motif),
                "the record and the prior count {} differently",
                String::from_utf8_lossy(bases),
            );
        }
    }

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
    /// and thirty. A mode keyed on the reference would sit a whole unit away from where the
    /// cohort's reads actually are.
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
        admit_promoted_sequences(
            observation,
            &observation.per_sample,
            &detail_of(b"AT"),
            &config,
            &mut scratch,
        )
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
                b"ATATATATAG".to_vec(),
                b"ATATATATAT".to_vec()
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
                b"ATATATATAG".to_vec(),
                b"ATATATATAT".to_vec()
            ],
            "the ten-read spelling is the one cut, and the survivors are in the order of their bases"
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

    // ---- D2: the entry point ----

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
        select_ssr(observation, &config, &mut SelectionScratch::new()).selection
    }

    /// The same narrowing, keeping the repeat counts beside the selection.
    fn selected_with_counts(observation: &CohortObservation, ploidy: u8) -> SsrLocusSelection {
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

    /// **A tract is genotyped however many of its reads sit off the motif grid** (owner,
    /// 2026-10-22; spec §7). GIAB HG002's `chr22:19052415` in miniature: a `CA` repeat of 13
    /// units, 13 reads at the reference length, 12 one unit shorter, and 3 a part of a unit longer
    /// — the 5-base insertion the aligner placed at the reads' end. The periodicity check refused
    /// this tract and called it over the reference alone; now the heterozygous deletion is a
    /// candidate beside the reference, and the three reads are left to the stutter model.
    #[test]
    fn a_tract_with_reads_off_the_motif_grid_is_still_genotyped() {
        let reference = b"CA".repeat(13);
        let one_unit_shorter = b"CA".repeat(12);
        let five_bases_longer = [b"CA".repeat(13), b"CACAC".to_vec()].concat();
        let observation = tract_of(
            &[&reference, &one_unit_shorter, &five_bases_longer],
            b"CA",
            vec![sample_showing(
                0,
                vec![row(0, 13, -13.0), row(1, 12, -12.0), row(2, 3, -3.0)],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(
            admitted_bases(&selection),
            vec![reference, one_unit_shorter]
        );
    }

    /// **A length the reads earned, split across spellings by sequencing errors, keeps its purest
    /// spelling** (owner, 2026-10-10). The stuttering benchmark's case: four reads at the
    /// reference length and three a unit shorter, one pure and two with one wrong base each —
    /// three spellings of one read apiece, none clearing the bar of two, while the length's three
    /// reads do. The pure spelling is admitted; the two carrying an error are not.
    #[test]
    fn a_length_split_across_spellings_keeps_its_purest_spelling() {
        let observation = tract_of(
            &[b"ATATATATAT", b"ATACATAT", b"ATATATAT", b"ATATAGAT"],
            b"AT",
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row(1, 1, -1.0),
                    row(2, 1, -1.0),
                    row(3, 1, -1.0),
                ],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(
            admitted_bases(&selection),
            vec![b"ATATATATAT".to_vec(), b"ATATATAT".to_vec()]
        );
    }

    /// **The rule keeps one spelling for a length nothing else stands for, and no more.** Where a
    /// spelling of the earned length cleared the bar on its own, the length is already stood for
    /// and its one-read siblings stay out; and a length only the `±1` rescue reached — one read,
    /// earned by nobody — gets no spelling from this rule.
    #[test]
    fn a_length_already_stood_for_or_only_rescued_gets_nothing_more() {
        let stood_for = tract_of(
            &[b"ATATATATAT", b"ATATATAT", b"ATACATAT"],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 3, -3.0), row(2, 1, -1.0)],
            )],
        );
        assert_eq!(
            admitted_bases(&selected(&stood_for, 2)),
            vec![b"ATATATATAT".to_vec(), b"ATATATAT".to_vec()]
        );

        // Twenty reads at the reference length and one a unit longer: the sample resolved one
        // length, so the rescue promotes the occupied neighbour, but no sample's reads earned it.
        let only_rescued = tract_of(
            &[b"ATATATATAT", b"ATATATATATAT"],
            b"AT",
            vec![sample_showing(0, vec![row(0, 20, -20.0), row(1, 1, -1.0)])],
        );
        assert_eq!(
            admitted_bases(&selected(&only_rescued, 2)),
            vec![b"ATATATATAT".to_vec()]
        );
    }

    /// The purity tie-break reads the motif at its best rotation: a pure run scores zero
    /// whichever base of the unit it starts on, and each wrong base costs one.
    #[test]
    fn bases_off_the_motif_counts_wrong_bases_at_the_best_rotation() {
        let motif = Motif::new(b"AAGAG").expect("a five-base motif");
        assert_eq!(bases_off_the_motif(b"AAGAGAAGAG", &motif), 0);
        assert_eq!(bases_off_the_motif(b"GAGAAGAGAA", &motif), 0);
        assert_eq!(bases_off_the_motif(b"AAGAGAAAAG", &motif), 1);
        assert_eq!(bases_off_the_motif(b"AATAGAAAAG", &motif), 2);
    }

    /// **A sequence with no whole repeat is never admitted, however many reads it has** (owner,
    /// 2026-10-10): admitted, it would stop the whole locus downstream. Three reads spelled the
    /// tract as empty at a 6-base repeat on the stuttering benchmark; the other lengths stay.
    #[test]
    fn a_sequence_with_no_whole_repeat_is_left_out_and_the_rest_stays() {
        let observation = tract_of(
            &[b"AACAATAACAATAACAAT", b"", b"AACAATAACAAT"],
            b"AACAAT",
            vec![sample_showing(
                0,
                vec![row(0, 5, -5.0), row(1, 3, -3.0), row(2, 5, -5.0)],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(
            admitted_bases(&selection),
            vec![b"AACAATAACAATAACAAT".to_vec(), b"AACAATAACAAT".to_vec()]
        );
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
                b"ATATATATAG".to_vec(),
                b"ATATATATAT".to_vec()
            ]
        );
    }

    // ---- D3: the repeat counts the genotype prior takes ----

    /// **The repeat counts run parallel to the surviving candidates, reference at index 0.**
    ///
    /// Three lengths survive here — four, five and six repeats — and the reference is the
    /// four-repeat one. A count list that was parallel to the *merge's* table instead, or that
    /// dropped the reference's entry, would hand the prior a length for the wrong allele at every
    /// locus where anything was dropped, and nothing downstream could see it.
    #[test]
    fn the_repeat_counts_run_parallel_to_the_surviving_candidates() {
        let observation = tract_of(
            &[b"ATATATAT", b"ATATATATATAT", b"ATATATATAT", b"ATATATATA"],
            b"AT",
            vec![
                sample_showing(0, vec![row(0, 20, -20.0), row(2, 20, -20.0)]),
                sample_showing(1, vec![row(1, 20, -20.0), row(3, 1, -1.0)]),
            ],
        );
        let narrowed = selected_with_counts(&observation, 2);
        assert_eq!(
            admitted_bases(&narrowed.selection),
            vec![
                b"ATATATAT".to_vec(),
                b"ATATATATAT".to_vec(),
                b"ATATATATATAT".to_vec()
            ],
            "the one-read nine-base sequence is the only one dropped"
        );
        assert_eq!(
            narrowed.repeat_counts,
            vec![4, 5, 6],
            "each surviving candidate's own repeat count, in candidate-id order"
        );
        assert_eq!(
            narrowed.repeat_counts.len(),
            narrowed.selection.alleles().len()
        );
    }

    /// A repeat count is the floor of the length over the period, so a sequence carrying a part
    /// copy takes the count of the rung it sits on — the same integer the ladder keyed it by, from
    /// the same derivation. Here the reference is 7 bases, an odd number for a two-base motif,
    /// and the alternative at 9 bases is one whole copy longer.
    #[test]
    fn a_part_copy_carries_the_floored_count_the_ladder_keyed_it_by() {
        let observation = tract_of(
            &[b"ATATATA", b"ATATATATA"],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0)],
            )],
        );
        let narrowed = selected_with_counts(&observation, 2);
        assert_eq!(
            narrowed.repeat_counts,
            vec![3, 4],
            "seven bases floor to three copies of a two-base motif and nine to four"
        );
    }

    /// A homopolymer's counts are its lengths, which is the case where a wrong period would be
    /// invisible in the ladder and visible here.
    #[test]
    fn a_homopolymers_repeat_counts_are_its_lengths() {
        let observation = tract_of(
            &[b"AAAA", b"AAAAAA"],
            b"A",
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 20, -20.0)],
            )],
        );
        let narrowed = selected_with_counts(&observation, 2);
        assert_eq!(narrowed.repeat_counts, vec![4, 6]);
    }

    /// **The counts handed to the prior agree with the ladder's own rung keys**, allele by allele.
    ///
    /// This is the coupling the field exists for, asserted rather than described: the ladder keys
    /// its rungs by one floor division and the prior is indexed by another, and they have to be
    /// the same integer. The test walks every surviving candidate back to the merge index it came
    /// from and compares the two.
    #[test]
    fn every_candidates_count_is_the_rung_the_ladder_put_it_on() {
        let observation = tract_of(
            &[b"ATATATAT", b"ATATATATATAT", b"ATATATATAT", b"ATATATATA"],
            b"AT",
            vec![
                sample_showing(0, vec![row(0, 20, -20.0), row(2, 20, -20.0)]),
                sample_showing(1, vec![row(1, 20, -20.0), row(3, 20, -20.0)]),
            ],
        );
        let rule = support_rule_of(2, 0.0);
        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: rule,
                ..SsrSelectionConfig::at_ploidy(diploid()).shared
            },
            ..SsrSelectionConfig::at_ploidy(diploid())
        };
        let mut scratch = SelectionScratch::new();
        let narrowed = select_ssr(&observation, &config, &mut scratch);

        for table_index in 0..observation.alleles.len() {
            let Some(candidate) = narrowed.selection.remap().candidate_for(table_index) else {
                continue;
            };
            let from_the_ladder = scratch
                .ladder
                .repeat_count_at(scratch.ladder.rung_of_table_index(table_index));
            assert_eq!(
                narrowed.repeat_counts[usize::from(candidate.get())],
                from_the_ladder,
                "candidate {candidate:?} came from merge allele {table_index}"
            );
        }
    }

    /// **A sample that earned a sequence the *cut* removed must still be callable** — only the
    /// **cap** makes a sample uncallable, and on this path an allele can be dropped for a third
    /// reason the ordinary path does not have.
    ///
    /// The shared leftover infers "the cap cut it" from "this sample earned it and it is gone",
    /// and on the ordinary path that inference holds: an allele a sample earned always cleared
    /// the bar, so it was always a candidate for the cap. **Here it does not.** The top-`ploidy`
    /// cut removes whole rungs, so a sample can have earned a sequence that admission dropped
    /// because its length was not promoted — with the cap never binding at all.
    ///
    /// This locus has one sample with four lengths, all earned, and two copies. The cap is 32 and
    /// cuts nothing. The two unpromoted lengths must not make the sample uncallable.
    #[test]
    fn a_sequence_dropped_by_the_ploidy_cut_does_not_make_a_sample_uncallable() {
        let observation = tract_of(
            &[b"ATATAT", b"ATATATAT", b"ATATATATAT", b"ATATATATATAT"],
            b"AT",
            vec![sample_showing(
                0,
                vec![
                    row(0, 30, -30.0),
                    row(1, 30, -30.0),
                    row(2, 5, -5.0),
                    row(3, 5, -5.0),
                ],
            )],
        );
        let selection = selected(&observation, 2);
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Selected,
            "the cap of 32 cuts nothing here"
        );
        assert_eq!(
            admitted_bases(&selection),
            vec![b"ATATAT".to_vec(), b"ATATATAT".to_vec()],
            "two copies, so the five-read lengths are not promoted"
        );
        assert!(
            !selection.unmatched()[0].genotype_must_be_missing(),
            "no allele was cut by the cap, so no sample is uncallable — the ploidy cut is not \
             the cap, and a genotype thrown away here is thrown away for a reason that did not \
             happen"
        );
    }

    // ---- the discovery pre-pass (STR loop plan Milestone E; research doc §6.5) ----

    /// `n` whole copies of the `AT` unit.
    fn at_repeats(n: usize) -> Vec<u8> {
        b"AT".repeat(n)
    }

    /// The path's defaults at `ploidy` and `cap`, with the bar spelled out and **the discovery
    /// pre-pass on at its shipped bar** — 2 reads and 15% of the sample's spanning reads
    /// ([`DEFAULT_DISCOVERY_BAR`]).
    fn discovering_config(ploidy: u8, cap: u16) -> SsrSelectionConfig {
        let copies = Ploidy::try_new(ploidy).expect("at least one copy");
        SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: support_rule_of(2, 0.0),
                max_candidate_alleles: MaxCandidateAlleles::new(cap).expect("a cap of two or more"),
            },
            discovery: Some(DEFAULT_DISCOVERY_BAR),
            ..SsrSelectionConfig::at_ploidy(copies)
        }
    }

    /// **E1's flagship shape, as one merge locus**: a diploid sample showing 12 reads at the
    /// reference's 8 repeats, 10 at 7 — its contraction slip — and 8 at 5. Nomination's two
    /// peaks are 8 and 7, so the 5 is invisible to selection by construction, however many
    /// reads carry it.
    fn tract_hiding_a_third_length() -> CohortObservation {
        tract_of(
            &[&at_repeats(8), &at_repeats(7), &at_repeats(5)],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 12, -12.0), row(1, 10, -10.0), row(2, 8, -8.0)],
            )],
        )
    }

    /// **Discovery ships off**: the constructor of the defaults runs no pre-pass, so a run that
    /// never asks for it selects byte-identically to the path before discovery was wired.
    #[test]
    fn the_defaults_run_no_discovery_pre_pass() {
        assert!(SsrSelectionConfig::at_ploidy(diploid()).discovery.is_none());
    }

    /// **The case the pre-pass exists for**: the third length one sample carries clears the bar
    /// — 8 of its 30 spanning reads, against `max(2, ceil(0.15 × 30)) = 5` — and is admitted,
    /// where the same locus with discovery off drops it at the top-`ploidy` cut.
    ///
    /// And the bookkeeping moves with it: the merge index gains a candidate id in the remap, so
    /// the sample's 8 reads leave the "no written allele explains these" leftover — which is
    /// what puts them in the record's `AD` — and the table the loop is handed carries all three
    /// candidates with their repeat counts.
    #[test]
    fn a_third_length_clearing_the_bar_is_admitted_and_its_reads_leave_the_leftover() {
        let observation = tract_hiding_a_third_length();

        let off = selected_with_counts(&observation, 2);
        assert_eq!(
            admitted_bases(&off.selection),
            vec![at_repeats(8), at_repeats(7)],
            "with discovery off the two peaks are all a diploid sample puts forward"
        );
        assert_eq!(
            off.selection.unmatched()[0].num_reads,
            8,
            "the third length's reads sit in the leftover: no written allele explains them"
        );

        let on = select_ssr(
            &observation,
            &discovering_config(2, 32),
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            admitted_bases(&on.selection),
            vec![at_repeats(8), at_repeats(5), at_repeats(7)],
            "the discovered length joins the table, in the order of the bases"
        );
        assert_eq!(
            on.repeat_counts,
            vec![8, 5, 7],
            "and the loop sees three candidates, each with the repeat count the prior indexes by"
        );
        assert_eq!(
            on.selection.remap().candidate_for(2),
            Some(AlleleId(1)),
            "the merge's third sequence now has a candidate id"
        );
        assert_eq!(
            on.selection.unmatched()[0].num_reads,
            0,
            "its 8 reads left the leftover — they are explained by a written allele now, which \
             is what books them in AD"
        );
        assert_eq!(on.selection.verdict(), SelectionVerdict::Selected);
    }

    /// **A length below either half of the bar is not admitted.** The floor half: one read is
    /// below 2, however shallow the sample. The share half: 4 reads of 40 spanning is 10 in
    /// 100, below 15 — and 6 of 40 clears it, so the refusal is the bar and not the wiring.
    #[test]
    fn a_length_below_either_half_of_the_bar_is_not_discovered() {
        let config = discovering_config(2, 32);

        let below_the_floor = tract_of(
            &[&at_repeats(8), &at_repeats(7), &at_repeats(5)],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 4, -4.0), row(1, 3, -3.0), row(2, 1, -1.0)],
            )],
        );
        assert_eq!(
            admitted_bases(
                &select_ssr(&below_the_floor, &config, &mut SelectionScratch::new()).selection
            ),
            vec![at_repeats(8), at_repeats(7)],
            "one read never mints an allele, whatever the share asks at 8 spanning reads"
        );

        let below_the_share = tract_of(
            &[&at_repeats(8), &at_repeats(7), &at_repeats(5)],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 16, -16.0), row(2, 4, -4.0)],
            )],
        );
        assert_eq!(
            admitted_bases(
                &select_ssr(&below_the_share, &config, &mut SelectionScratch::new()).selection
            ),
            vec![at_repeats(8), at_repeats(7)],
            "four reads of forty spanning is 10 in 100, below the bar's 15"
        );

        let clearing_the_share = tract_of(
            &[&at_repeats(8), &at_repeats(7), &at_repeats(5)],
            b"AT",
            vec![sample_showing(
                0,
                vec![row(0, 20, -20.0), row(1, 14, -14.0), row(2, 6, -6.0)],
            )],
        );
        assert_eq!(
            admitted_bases(
                &select_ssr(&clearing_the_share, &config, &mut SelectionScratch::new()).selection
            ),
            vec![at_repeats(8), at_repeats(5), at_repeats(7)],
            "six of forty is 15 in 100 and clears it"
        );
    }

    /// **One sample has to clear the bar alone** — three samples showing a length once each is
    /// one read anywhere, and one read is what a stutter product looks like.
    #[test]
    fn discovery_does_not_pool_samples_against_the_bar() {
        let observation = tract_of(
            &[&at_repeats(8), &at_repeats(5)],
            b"AT",
            vec![
                sample_showing(0, vec![row(0, 5, -5.0), row(1, 1, -1.0)]),
                sample_showing(1, vec![row(0, 5, -5.0), row(1, 1, -1.0)]),
                sample_showing(2, vec![row(0, 5, -5.0), row(1, 1, -1.0)]),
            ],
        );
        let on = select_ssr(
            &observation,
            &discovering_config(2, 32),
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            admitted_bases(&on.selection),
            vec![at_repeats(8)],
            "three reads across the cohort is one read in each sample, and no single sample \
             clears the floor of two"
        );
    }

    /// **A sequence carrying no whole motif copy is refused by the pre-pass, where the reason
    /// is** — admitted, it would reach the driver's non-zero conversion and stop the whole
    /// locus one step later. Its reads clear the bar comfortably; the refusal is the missing
    /// rung, not the evidence.
    ///
    /// The sequence is a whole-tract deletion at a homopolymer.
    #[test]
    fn discovery_does_not_admit_a_sequence_below_one_whole_repeat() {
        let observation = tract_of(
            &[b"AAAA", b"AAA", b""],
            b"A",
            vec![sample_showing(
                0,
                vec![row(0, 12, -12.0), row(1, 10, -10.0), row(2, 8, -8.0)],
            )],
        );
        let on = select_ssr(
            &observation,
            &discovering_config(2, 32),
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            admitted_bases(&on.selection),
            vec![b"AAAA".to_vec(), b"AAA".to_vec()],
            "zero copies has no rung on the ladder the stutter model is written over, so the \
             deleted tract is not offered however many reads showed it"
        );
        assert_eq!(
            on.repeat_counts,
            vec![4, 3],
            "so the counts handed to the loop stay non-zero and the locus stays callable"
        );
    }

    /// **A second pass over the same evidence admits nothing — the E1 property, held through
    /// the wiring.** The eligible set is a function of the observations and the candidate table
    /// alone: after one `select_ssr` with discovery on, every sequence the decision nominates
    /// already holds a candidate id, so a second round would have nothing left to add — which
    /// is why discovery is a pre-pass here and not a loop around the caller.
    #[test]
    fn a_second_discovery_pass_over_the_same_evidence_admits_nothing() {
        let observation = tract_hiding_a_third_length();
        let mut scratch = SelectionScratch::new();
        let on = select_ssr(&observation, &discovering_config(2, 32), &mut scratch);

        let motif = Motif::new(b"AT").expect("a two-base motif");
        let mut nominated_again = Vec::new();
        nominate_discovered_sequences(
            &observation,
            &motif,
            DEFAULT_DISCOVERY_BAR,
            &mut nominated_again,
        );
        let eligible_for_a_second_round: Vec<usize> = nominated_again
            .iter()
            .enumerate()
            .filter(|&(index, &nominated)| {
                nominated && on.selection.remap().candidate_for(index).is_none()
            })
            .map(|(index, _)| index)
            .collect();
        assert!(
            nominated_again.contains(&true),
            "the decision still nominates — what has changed is the table, not the evidence"
        );
        assert!(
            eligible_for_a_second_round.is_empty(),
            "every nominated sequence is already a candidate, so a round after the pre-pass \
             finds nothing: merge indices {eligible_for_a_second_round:?} would contradict E1"
        );
    }

    /// **A discovered allele the cap cuts follows the same truncation rules as any other**: it
    /// is ranked by [`compare_best_first`](super::super::compare_best_first) with the rest,
    /// reported in the verdict's count, and the sample whose reads earned it is marked
    /// uncallable — exactly what the cap does to an allele admitted the ordinary way.
    #[test]
    fn a_discovered_allele_the_cap_cuts_follows_the_shared_truncation_rules() {
        let observation = tract_hiding_a_third_length();
        // Two alleles counting the reference, so of the two alternatives — the 7 on its 10
        // reads and the discovered 5 on its 8 — the ranking keeps one.
        let on = select_ssr(
            &observation,
            &discovering_config(2, 2),
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            on.selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 },
            "the discovered allele entered the ranking and the cap's count says one was cut"
        );
        assert_eq!(
            admitted_bases(&on.selection),
            vec![at_repeats(8), at_repeats(7)],
            "the 7 outranks the discovered 5 — 10 reads of 30 against 8 of 30 on the first key"
        );
        assert!(
            on.selection.unmatched()[0].genotype_must_be_missing(),
            "the sample's own 8 reads were on the sequence the cap cut, so it is uncallable \
             here — a discovered allele is not exempt from the rule that makes truncation \
             defensible"
        );
    }

    // ── STR bundles (`str_bundle_loci.md` §7.2) ──

    /// `GG (CA)8 GT (TTG)5 C`: two runs, of periods 2 and 3, fixed bases around and between.
    fn bundle_hull(ca_copies: usize, ttg_copies: usize) -> Vec<u8> {
        [
            &b"GG"[..],
            &b"CA".repeat(ca_copies),
            b"GT",
            &b"TTG".repeat(ttg_copies),
            b"C",
        ]
        .concat()
    }

    fn bundle_detail() -> SsrBundleDetail {
        use crate::region_typing::reference_runs::ReferenceRun;
        SsrBundleDetail {
            left_flank: Box::from(&b"ACGTACGTAC"[..]),
            right_flank: Box::from(&b"GTCAGTCAGT"[..]),
            runs: vec![
                ReferenceRun::new(2, 16, Motif::new(b"CA").unwrap()).unwrap(),
                ReferenceRun::new(20, 15, Motif::new(b"TTG").unwrap()).unwrap(),
            ]
            .into_boxed_slice(),
        }
    }

    fn bundle_of(
        alleles: &[&[u8]],
        detail: SsrBundleDetail,
        per_sample: Vec<SampleSupport>,
    ) -> CohortObservation {
        let mut observation = locus_of(alleles, per_sample);
        observation.kind = LocusKind::SsrBundle(detail);
        observation
    }

    fn bundle_selected(observation: &CohortObservation) -> SsrBundleSelection {
        let config = SsrSelectionConfig {
            shared: CandidateSelectionConfig {
                min_allele_support: support_rule_of(2, 0.1),
                ..SsrSelectionConfig::at_ploidy(diploid()).shared
            },
            ..SsrSelectionConfig::at_ploidy(diploid())
        };
        select_ssr_bundle(observation, &config, &mut SelectionScratch::new())
    }

    /// A heterozygote one `CA` copy longer than the reference is called over both, with the run
    /// counts of each; a stutter spelling below a tenth of the sample's reads is not a candidate.
    #[test]
    fn a_bundle_heterozygote_keeps_both_alleles_and_their_run_counts() {
        let reference = bundle_hull(8, 5);
        let longer = bundle_hull(9, 5);
        let stutter = bundle_hull(8, 4);
        let observation = bundle_of(
            &[&reference, &longer, &stutter],
            bundle_detail(),
            vec![sample_showing(
                0,
                vec![row(0, 12, -12.0), row(1, 11, -11.0), row(2, 2, -2.0)],
            )],
        );
        let narrowed = bundle_selected(&observation);
        assert_eq!(narrowed.selection.alleles().len(), 2);
        assert_eq!(
            narrowed.selection.alleles().iter().nth(1),
            Some(&longer[..])
        );
        let counts: Vec<Vec<u32>> = narrowed
            .parses
            .iter()
            .map(|parse| parse.copies().collect())
            .collect();
        assert_eq!(counts, vec![vec![8, 5], vec![9, 5]]);
        assert!(matches!(
            narrowed.selection.alleles().kind(),
            LocusKind::SsrBundle(_)
        ));
    }

    /// **The rescue's neighbours are one copy either side in each run separately** — never two
    /// runs at once — and only vectors some allele of the locus carries.
    #[test]
    fn a_bundle_rungs_neighbours_differ_by_one_copy_in_one_run() {
        let vectors: Vec<Box<[u32]>> = [[7, 5], [8, 4], [8, 5], [8, 6], [9, 5], [9, 6]]
            .iter()
            .map(|vector| Box::from(&vector[..]))
            .collect();
        let mut out = Vec::new();
        RungNeighbours::OneCopyEitherSideInEachRun(&vectors).fill(2, &mut out);
        out.sort_unstable();
        assert_eq!(
            out,
            vec![0, 1, 3, 4],
            "[7,5], [8,4], [8,6], [9,5] — not [9,6]"
        );
        RungNeighbours::OneRepeatEitherSide.fill(0, &mut out);
        assert_eq!(out, vec![1]);
    }

    /// **A bundle is never refused for length changes that are not whole units of its runs**
    /// (owner, 2026-10-12; `str_bundle_loci.md` §7.2): the fixed sequence between runs carries
    /// real indels of any length. A one-base insertion carried by half the reads is a candidate
    /// at a bundle whose runs are both of period 2, and at one with no run at all.
    #[test]
    fn a_bundle_keeps_a_length_change_that_is_no_whole_number_of_run_units() {
        use crate::region_typing::reference_runs::ReferenceRun;
        let samples = || {
            vec![sample_showing(
                0,
                vec![row(0, 10, -10.0), row(1, 10, -10.0)],
            )]
        };
        let period_two_hull =
            [&b"GG"[..], &b"CA".repeat(8), b"GG", &b"AT".repeat(6), b"C"].concat();
        let period_two = SsrBundleDetail {
            runs: vec![
                ReferenceRun::new(2, 16, Motif::new(b"CA").unwrap()).unwrap(),
                ReferenceRun::new(20, 12, Motif::new(b"AT").unwrap()).unwrap(),
            ]
            .into_boxed_slice(),
            ..bundle_detail()
        };
        let no_run = SsrBundleDetail {
            runs: Box::new([]),
            ..bundle_detail()
        };
        for (reference, detail) in [
            (bundle_hull(8, 5), bundle_detail()),
            (period_two_hull, period_two),
            (bundle_hull(8, 5), no_run),
        ] {
            // A base inserted in the fixed sequence between the runs.
            let one_base_longer = [&reference[..18], b"T", &reference[18..]].concat();
            let runs = detail.runs.len();
            let observation = bundle_of(&[&reference, &one_base_longer], detail, samples());
            let narrowed = bundle_selected(&observation);
            assert_eq!(
                narrowed.selection.verdict(),
                SelectionVerdict::Selected,
                "{runs} runs"
            );
            assert_eq!(narrowed.selection.alleles().len(), 2, "{runs} runs");
            assert_eq!(narrowed.parses.len(), 2, "{runs} runs");
        }
    }

    /// **The purest spelling of a rung is the one with the fewest edited bases in its run parse**:
    /// a run-count change is no edit, a substituted base is one.
    #[test]
    fn a_bundle_spellings_impurity_is_its_parses_edited_bases() {
        let reference = bundle_hull(8, 5);
        let mut substituted = bundle_hull(9, 5);
        substituted[5] = b'T';
        let observation = bundle_of(
            &[&reference, &bundle_hull(9, 5), &substituted],
            bundle_detail(),
            Vec::new(),
        );
        let parses: Vec<RunParse> = observation
            .alleles
            .iter()
            .map(|allele| parse_into_runs(allele, &reference, &bundle_detail().runs))
            .collect();
        let detail = bundle_detail();
        let ground = SelectionGround::Bundle {
            detail: &detail,
            parses: &parses,
        };
        assert_eq!(ground.bases_off_the_repeat(&observation, 1), 0);
        assert_eq!(ground.bases_off_the_repeat(&observation, 2), 1);
        assert!(ground.can_score(&observation, 0));
    }
}
