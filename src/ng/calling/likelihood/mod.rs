//! Step 7 — how probable this sample's reads are, given a genotype.
//!
//! One call answers *given that this sample's copies of this locus are such-and-such, how
//! probable are the reads we actually saw here?* — one number per candidate genotype, in
//! log space. The caller adds those numbers to the genotype prior
//! ([`super::genotype_prior`], step 8) and normalises; the sum is the posterior ng emits
//! (`doc/devel/ng/spec/read_likelihoods.md` §1).
//!
//! ## Two paths, one formula
//!
//! Both paths compute the same thing — a copy-weighted mixture over the genotype's
//! alleles, wrapped in one logarithm per observation — and differ in a single term, the
//! **emission**: the probability that one copy of one allele produced one observed
//! sequence. At an ordinary site a read can show a base the individual does not carry, and
//! one rate covers it. At a repeat tract a read can also show a whole repeat more or fewer
//! than the individual carries, because the copying enzyme slipped, and that happens
//! thousands of times more often than a base is misread (spec §2.1, §2.2).
//!
//! ## What this module never does
//!
//! **It fits nothing.** Every number it reads is fitted before calling starts, by
//! [`crate::ng::parameter_estimation`]: the per-read-group error rate and the scale that
//! calibrates it, the contamination fraction, and the four numbers that say how a repeat
//! tract slips — which are fitted per **stratum**, one *(motif length, repeat count)* cell,
//! because how much a tract slips depends on how many repeats it has more than on anything
//! else. This module does not choose the candidate alleles, does not run the caller's loop,
//! does not compute a posterior or a genotype, and does not decide whether a site is
//! emitted (spec §1.2).
//!
//! ## Vocabulary a reader needs before the types below
//!
//! An **observation** is one distinct sequence the reads showed at a locus, with the
//! number of reads that showed it — not one read, and not one allele. Its identity is
//! `(bases, read witness, read group)`, so the same sequence seen by a read that spanned
//! the whole locus and by one that ran out inside it is two observations. A **complete**
//! observation pins what the sample carries there; a **partial** one proves only that the
//! locus is *at least* what the read got through, which the statistician calls
//! **censored** (spec §1.3).
//!
//! A **read group** is one `@RG` — one lane, one chemistry. It is part of an observation's
//! identity here for a reason that is easy to lose: the formula puts a logarithm outside a
//! sum, so an observation's reads may be pooled into one term **only if every one of them
//! would have got the same number**, and two reads showing the same bases from two lanes
//! have different error rates (spec §2.3).
//!
//! An allele's **projection** is the bases it spells over one locus's span — what the
//! merge unified the samples' reads onto, so that two samples showing the same sequence
//! land on one allele. It is a whole-span thing, which is why a read that saw only part of
//! the span cannot be compared against it (spec §5.3).
//!
//! ## Two allele tables, and they are not the same numbering
//!
//! **This is the trap this module's types are shaped around.** The cohort merge unifies
//! every distinct sequence the whole cohort showed into one table, and its rows index that
//! table. Candidate *selection* — step 6, which does not exist yet — then keeps some of
//! those alleles and drops the rest, and **dropping allele *k* renumbers every allele above
//! it**. [`AlleleId`] names an index into the surviving table, the one a genotype is built
//! over, and nothing in the two numberings makes them agree.
//!
//! So there is no conversion from a merge row to an [`AlleleId`] that this module could
//! perform on its own, and it does not offer one: the mapping arrives as an argument
//! ([`GenericObservation::fill_from_supported_alleles`]), and the reads on alleles the
//! mapping drops are exactly the reads whose quality the formula pools into
//! [`GenericSampleEvidence::unmatched_q_sum`].

use crate::ng::locus_generation::{ReadWitness, SequenceObservation, SsrDetail};
use crate::ng::run::cohort_merge::build::{PartialObservation, SupportedAllele};
use crate::ng::types::{AlleleId, ReadGroupId};

/// The merge's fold of every read that showed one allele from one read group — one term of
/// the SNP/indel formula's first two sums.
///
/// **Four numbers where the merge's own row carries eight** — [`SupportedAllele`]'s allele
/// and read group, and the six of [`AlleleSupport`] — and the four that are dropped are
/// dropped on purpose: the forward-strand count, the two mapping-quality moments and the
/// read-position count belong to the site filters that run *after* genotyping, and none of
/// them enters a likelihood (spec §1.4). What is kept is what the formula reads, and the
/// type is what stops a filter statistic reaching it.
///
/// **`q_sum` is a sum of logarithms, not a probability.** It is Σ `ln P(this read is
/// wrong)` over the reads folded here, which is exactly the shape the formula needs: a read
/// the genotype cannot explain is charged its own log error, so summing the logs and
/// charging the sum are the same arithmetic as charging each read separately. That is what
/// makes the likelihood of a list of reads and the likelihood of the merge's fold of those
/// same reads agree to the last bit.
///
/// `Copy`, because it is four scalars the row builder indexes in an inner loop and nothing
/// it owns needs dropping. `tests::the_observation_stays_four_scalars_wide` pins the size, so
/// that a field which would make copying expensive — anything owning heap — has to be argued
/// for rather than arriving unnoticed.
///
/// [`AlleleSupport`]: crate::ng::run::cohort_merge::build::AlleleSupport
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct GenericObservation {
    /// Which allele these reads showed — `a(o)` in the formula, an index into the locus's
    /// [`CandidateAlleles`](super::CandidateAlleles).
    ///
    /// **A candidate index and not the merge's own**, which is why nothing here can mint
    /// one: see the module's *Two allele tables* note.
    pub allele: AlleleId,
    /// Which read group they came from — `r(o)`, and **part of the identity**, never a
    /// label to fold away. It selects the calibration scale and the contamination fraction
    /// this term is charged under (spec §2.3).
    pub read_group: ReadGroupId,
    /// How many reads showed it — `n_o`.
    pub num_reads: u32,
    /// Σ `ln P(error)` over those reads — `q_sum_o`, straight off the merge's row.
    pub q_sum: f64,
}

impl GenericObservation {
    /// The four numbers of one merge row, under the candidate id selection gave that row's
    /// allele.
    ///
    /// **The id is an argument because this module cannot compute it.** The row's own
    /// `allele` field indexes the merge's unification table and [`AlleleId`] indexes the
    /// candidate table, and only selection knows how one maps onto the other (see the
    /// module's *Two allele tables* note).
    #[must_use]
    pub fn of_supported_allele(row: &SupportedAllele, allele: AlleleId) -> Self {
        Self {
            allele,
            read_group: row.read_group,
            num_reads: row.support.num_reads,
            q_sum: row.support.q_sum,
        }
    }

    /// Narrow one sample's merge rows into `out` — cleared first — under selection's
    /// mapping, and return the pooled quality of the rows it dropped.
    ///
    /// `candidate_of_merge_allele` is indexed by the merge's own allele index and says which
    /// candidate that allele became, or `None` where selection dropped it.
    ///
    /// **The return value is not the whole of [`GenericSampleEvidence::unmatched_q_sum`],
    /// and calling it that would be wrong.** It is the part this call discarded — Σ `q_sum`
    /// over the rows whose allele has no candidate. Selection pools its own leftovers too,
    /// and the specification that says how is not written; a caller adds this to whatever
    /// selection hands it. What matters here is that a function that drops rows says what it
    /// dropped, because those reads are evidence and the data likelihood is compared between
    /// loci (spec §3.3).
    ///
    /// **`out` is caller scratch, cleared and refilled**, so a worker holds one buffer
    /// across every sample of every locus and the row function still allocates nothing
    /// (spec §8). There is no borrow to be had instead: a merge row is 48 bytes and this one
    /// is 24, so no reinterpretation of the merge's `Vec` exists and somebody has to fill a
    /// parallel buffer.
    ///
    /// **The output order is the input order**, so `out` is ascending on
    /// `(candidate allele, read group)` exactly when selection's mapping keeps alleles in
    /// order — which the renumbering does, since it drops alleles and slides the rest down.
    /// [`GenericSampleEvidence::new`] is where that promise is checked.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a row whose merge allele index is past the end of
    /// `candidate_of_merge_allele`. That is a mapping built against a different locus's
    /// allele table, and the alternative to panicking is scoring reads against whatever
    /// allele happens to sit at that index (spec §8; `per_group_merger.rs:1963` is the
    /// precedent for holding a structural assertion in release).
    pub fn fill_from_supported_alleles(
        rows: &[SupportedAllele],
        candidate_of_merge_allele: &[Option<AlleleId>],
        out: &mut Vec<Self>,
    ) -> f64 {
        out.clear();
        let mut dropped_q_sum = 0.0;
        for row in rows {
            let mapped = candidate_of_merge_allele
                .get(row.allele)
                .unwrap_or_else(|| {
                    panic!(
                        "merge allele {} has no entry in a mapping of {} alleles, so the \
                         mapping was built against a different locus's table",
                        row.allele,
                        candidate_of_merge_allele.len()
                    )
                });
            match mapped {
                Some(allele) => out.push(Self::of_supported_allele(row, *allele)),
                None => dropped_q_sum += row.support.q_sum,
            }
        }
        dropped_q_sum
    }
}

/// One sample's evidence at one cohort locus, as the SNP/indel row consumes it.
///
/// **A view, and the two halves borrow different things.** `partials` borrows the merge's
/// own rows, which is why there is one [`PartialObservation`] type rather than two.
/// `supported` cannot: the merge's row is twice the width of the one the formula reads, so
/// it borrows a staging buffer some caller filled with
/// [`GenericObservation::fill_from_supported_alleles`] and holds across loci. Either way the
/// row function allocates nothing per sample (spec §8).
#[derive(Copy, Clone, Debug)]
pub struct GenericSampleEvidence<'a> {
    /// One entry per `(allele, read group)` this sample's **complete** reads showed, in the
    /// ascending pair order the merge builds them in.
    ///
    /// **The order is a determinism requirement, not tidiness.** Floating-point addition is
    /// not associative, so two runs that summed a sample's observations in different orders
    /// would produce log-likelihoods differing in the last bits, and the run's genotypes
    /// with them (spec §8). The merge sorts on the pair once per sample, on the one path
    /// both its branches converge to (`AlleleTable::assemble` in
    /// `src/ng/run/cohort_merge/build.rs`), and its own
    /// `the_rows_are_ordered_by_allele_then_read_group` pins it against a fixture that
    /// arrives out of order — but that is a test in another module and this view is built
    /// from a buffer, so [`new`](Self::new) checks the promise here as well.
    pub supported: &'a [GenericObservation],
    /// Σ `ln P(error)` pooled over reads matching **no** candidate allele — the support of
    /// the merge's alleles that candidate selection dropped, folded here by that selection.
    ///
    /// **The same number for every genotype, so it cancels in genotyping, and it is kept
    /// anyway**: the data likelihood also feeds emission and QUAL, where an absolute value
    /// is compared between loci (spec §3.3's `q_sum_other`).
    ///
    /// **Nothing in this crate produces the whole of it yet.** It is created by candidate
    /// *selection*, which has no specification — "whoever specifies selection owes the
    /// pool" — so the row takes it as an input and tests supply it from fixtures.
    /// [`GenericObservation::fill_from_supported_alleles`] returns the part it can see, the
    /// quality of the rows the mapping dropped.
    pub unmatched_q_sum: f64,
    /// The partial observations, bases and witnessed run intact.
    ///
    /// **Not folded onto alleles**, because there is no allele: a partial's bases cannot be
    /// compared against a whole-span allele. Padding one out to the span and interning it
    /// would put a sequence in the table that no molecule carried, and it would read as a
    /// *short* allele — the one direction the model must not be biased in. What a candidate
    /// is scored against is its projection *restricted to the positions the read witnessed*,
    /// so the run has to survive to the scoring (spec §5.1, §5.3).
    ///
    /// **The witnessed run and the bases are on different axes, and their lengths are not
    /// interchangeable.** The run counts *locus positions*; the bases are what the read
    /// showed over them, so the two differ by the net indel the read carried — a read
    /// carrying a two-base insertion and a two-base deletion inside the stretch comes back
    /// with as many bases as positions and is still not a positional match for any of them
    /// ([`PartialObservation::bases`]). Scoring indexes the *allele's* projection with the
    /// run, never the partial's own bases.
    pub partials: &'a [PartialObservation],
}

impl<'a> GenericSampleEvidence<'a> {
    /// A sample's evidence from the three parts the merge and candidate selection hand
    /// over.
    ///
    /// A constructor rather than a struct literal at each call site so that a fourth part,
    /// when one arrives, cannot be forgotten at one of them.
    ///
    /// # Panics
    ///
    /// In debug builds, on `supported` rows that are not strictly ascending on
    /// `(allele, read group)`. The merge builds them that way and the determinism of the
    /// sum rests on it (spec §8), but this constructor takes any slice: a call site that
    /// concatenated two samples' rows, or filtered and re-spliced them, or applied a
    /// selection mapping that reordered the alleles, would break it with nothing saying so.
    /// Debug-only because the check is linear in the rows and the property is the merge's to
    /// hold, not this type's to enforce at every call.
    #[must_use]
    pub fn new(
        supported: &'a [GenericObservation],
        unmatched_q_sum: f64,
        partials: &'a [PartialObservation],
    ) -> Self {
        debug_assert!(
            supported
                .windows(2)
                .all(|pair| (pair[0].allele, pair[0].read_group)
                    < (pair[1].allele, pair[1].read_group)),
            "the merge builds one row per (allele, read group) in ascending pair order and \
             the determinism of the sum rests on it (spec §8); these are not ascending: \
             {supported:?}"
        );
        Self {
            supported,
            unmatched_q_sum,
            partials,
        }
    }

    /// A sample that showed nothing at this locus.
    ///
    /// **Not a special case in the row function**, and that is the point of naming it: an
    /// empty sum is zero, so every genotype scores zero and the prior decides alone, which
    /// is the right answer rather than a branch (spec §3.3).
    #[must_use]
    pub fn empty() -> Self {
        Self {
            supported: &[],
            unmatched_q_sum: 0.0,
            partials: &[],
        }
    }
}

/// One sample's evidence at one repeat tract, as the STR row consumes it.
///
/// **The locus generator's own type, unchanged.** [`SequenceObservation`] already keys on
/// `(bases, witness, read group)` and carries a read count, so the aggregation contract
/// holds by construction and there is nothing for this module to re-shape (arch §2.2).
#[derive(Copy, Clone, Debug)]
pub struct SsrSampleEvidence<'a> {
    /// Every observation at the tract, complete and partial alike, in the order the
    /// generator sorted them — which is the fixed order spec §8 requires the sum to run in.
    ///
    /// **Reaching this field directly is how a partial gets scored as though it were
    /// complete**, which mis-scores it as a *short* allele, because its bases are a prefix
    /// of the truth (spec §5.1). [`complete_observations`](Self::complete_observations) and
    /// [`partial_observations`](Self::partial_observations) split it by witness and carry
    /// each observation's position in this slice, which is the key an emission cache is
    /// built on — so the guarded route is also the useful one.
    ///
    /// **Calling them a guard would overstate them, and an earlier version of this comment
    /// did.** The field is public, so nothing stops a consumer walking it — and the
    /// generator's own `SampleLocusObservations::observations` is public too, so the
    /// iterator these are named after was never a guard either. What the two methods buy is
    /// that scoring a partial as complete has to be written rather than fallen into, and
    /// that the split has one spelling instead of one per caller.
    pub observations: &'a [SequenceObservation],
    /// The tract's motif and flanks — what the emission needs to know which repeat it is
    /// scoring against.
    pub detail: &'a SsrDetail,
}

impl<'a> SsrSampleEvidence<'a> {
    /// A sample's observations at one tract, against that tract's repeat detail.
    #[must_use]
    pub fn new(observations: &'a [SequenceObservation], detail: &'a SsrDetail) -> Self {
        Self {
            observations,
            detail,
        }
    }

    /// The observations whose reads spanned the whole tract — the ones an emission may
    /// score directly — each with its position in
    /// [`observations`](Self::observations).
    ///
    /// **The position travels because the STR row's emission cache is keyed by it.** That
    /// cache is what makes a row cost `observations × candidates` rather than
    /// `observations × genotypes`, which is the design and not an optimisation (spec §8);
    /// an iterator yielding only the observation would leave the row re-walking the
    /// unguarded field to recover the key.
    ///
    /// **It splits on the witness and never on the bases.** The same bases seen by a read
    /// that spanned the tract and by one that ran out are two observations, and only the
    /// witness separates them (spec §1.3) — so the plausible shortcut, *is this shorter
    /// than the tract?*, puts both on the same side and is wrong.
    ///
    /// The `use<'a>` is load-bearing: in edition 2024 an `impl Trait` return captures every
    /// lifetime in scope, so without it the iterator borrows `self` and cannot outlive the
    /// view — which the row builder needs it to.
    pub fn complete_observations(
        &self,
    ) -> impl Iterator<Item = (usize, &'a SequenceObservation)> + use<'a> {
        self.observations
            .iter()
            .enumerate()
            .filter(|(_, observation)| observation.read_witness == ReadWitness::Complete)
    }

    /// The observations whose reads ran out inside the tract — the ones the censored term
    /// scores, and nothing else may — each with its position in
    /// [`observations`](Self::observations).
    ///
    /// **This is not a small corner of the evidence.** A tract's locus is as wide as the
    /// tract, so over half the overlapping reads are partial at a 60-base tract, and an
    /// allele longer than a read can only ever be witnessed partially (spec §5.4.1).
    ///
    /// **An exhaustive match rather than `!= Complete`**, so that a third [`ReadWitness`]
    /// variant is a compile error here rather than silently joining the censored term. The
    /// censored term is the one that must never receive an observation it was not written
    /// for, and this is the single place that decides what reaches it.
    pub fn partial_observations(
        &self,
    ) -> impl Iterator<Item = (usize, &'a SequenceObservation)> + use<'a> {
        self.observations
            .iter()
            .enumerate()
            .filter(|(_, observation)| match observation.read_witness {
                ReadWitness::Partial { .. } => true,
                ReadWitness::Complete => false,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::WitnessedLocusPositions;
    use crate::ng::run::cohort_merge::build::AlleleSupport;
    use crate::ng::types::Motif;

    /// A merge row with every field distinguishable, so a field crossed with another one
    /// changes an answer. `num_fwd`, the mapping moments and `placed_left` are deliberately
    /// not zero: a fixture whose uninteresting fields are all zero cannot show that a
    /// mapping of the interesting ones picked them up by accident.
    fn supported_row(
        allele: usize,
        read_group: u32,
        num_reads: u32,
        q_sum: f64,
    ) -> SupportedAllele {
        SupportedAllele {
            allele,
            read_group: ReadGroupId(read_group),
            support: AlleleSupport {
                num_reads,
                num_fwd: 7,
                q_sum,
                mapq_sum: 611,
                mapq_sum_sq: 37_442,
                placed_left: 3,
            },
        }
    }

    /// Selection kept every merge allele, in order — the mapping a locus gets before
    /// anything is pruned.
    fn kept_in_order(alleles: usize) -> Vec<Option<AlleleId>> {
        (0..alleles)
            .map(|index| Some(AlleleId(index as u16)))
            .collect()
    }

    /// A partial whose witnessed stretch and bases are **deliberately different lengths** —
    /// five locus positions against seven bases, the shape a read carrying a two-base
    /// insertion inside the stretch comes back with.
    ///
    /// A fixture where the two happen to be equal is the one a consumer indexing the bases
    /// by a locus offset would pass, and the merge's own documentation says outright that
    /// equality does not license indexing one with the other.
    fn partial_row(num_reads: u32, q_sum: f64) -> PartialObservation {
        PartialObservation {
            witnessed_in_locus: WitnessedLocusPositions::one_run_from_offset_and_length(2, 5)
                .expect("a five-position run from offset two is a witness"),
            read_group: ReadGroupId(1),
            bases: b"ACGTTAA"[..].into(),
            num_reads,
            q_sum,
        }
    }

    /// The tract's repeat detail, shared by the repeat-path fixtures.
    fn ssr_detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(b"AC").expect("AC is a motif"),
            left_flank: b"GGTT"[..].into(),
            right_flank: b"TTGG"[..].into(),
        }
    }

    /// A read that ran out after `positions` positions of the tract.
    fn ran_out_after(positions: u16) -> ReadWitness {
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, positions)
                .expect("a run from offset zero over some positions is a witness"),
        }
    }

    fn ssr_observation(bases: &[u8], witness: ReadWitness, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: witness,
            read_group: ReadGroupId(2),
            num_obs,
            num_fwd: 5,
            q_sum: -13.5,
            mapq_sum: 240,
            mapq_sum_sq: 14_400,
            placed_left: 1,
            chain_ids: Vec::new(),
        }
    }

    #[test]
    fn the_view_of_a_merge_row_keeps_the_four_numbers_the_formula_reads() {
        let observation =
            GenericObservation::of_supported_allele(&supported_row(3, 9, 17, -41.25), AlleleId(3));

        assert_eq!(observation.allele, AlleleId(3));
        assert_eq!(observation.read_group, ReadGroupId(9));
        assert_eq!(observation.num_reads, 17);
        assert_eq!(observation.q_sum, -41.25);
    }

    /// The type's doc promises four scalars and cheap copying. A field owning heap would
    /// make both false, and would also break the no-allocation contract the row function
    /// works under (spec §8), so the width is pinned rather than described.
    #[test]
    fn the_observation_stays_four_scalars_wide() {
        // 24 and not 18: the `f64` wants eight-byte alignment, so the two-byte allele id and
        // the two four-byte counts are padded out to sixteen before it.
        assert_eq!(std::mem::size_of::<GenericObservation>(), 24);
        assert_eq!(std::mem::align_of::<GenericObservation>(), 8);
        assert!(!std::mem::needs_drop::<GenericObservation>());
    }

    /// The merge's allele index is **not** the candidate id, and the mapping between them is
    /// the argument. Here selection has dropped merge allele 1, so merge allele 2 becomes
    /// candidate 1 — and a fill that had assumed the identity would have scored those reads
    /// against the wrong sequence with nothing saying so.
    #[test]
    fn a_dropped_allele_renumbers_the_ones_above_it_and_its_quality_comes_back() {
        let rows = [
            supported_row(0, 1, 11, -3.5),
            supported_row(1, 1, 4, -19.0),
            supported_row(2, 1, 6, -8.25),
        ];
        let mapping = [Some(AlleleId(0)), None, Some(AlleleId(1))];
        let mut out = Vec::new();

        let dropped = GenericObservation::fill_from_supported_alleles(&rows, &mapping, &mut out);

        assert_eq!(dropped, -19.0);
        assert_eq!(
            out.iter().map(|o| o.allele).collect::<Vec<_>>(),
            vec![AlleleId(0), AlleleId(1)]
        );
        assert_eq!(
            out.iter().map(|o| o.num_reads).collect::<Vec<_>>(),
            vec![11, 6]
        );
    }

    /// Nothing is dropped when selection keeps everything, and the pooled quality is then
    /// zero rather than the sum of what was kept — a fill that added on the wrong branch
    /// would return −30.75 here.
    #[test]
    fn keeping_every_allele_drops_no_quality() {
        let rows = [
            supported_row(0, 1, 11, -3.5),
            supported_row(1, 2, 4, -27.25),
        ];
        let mut out = Vec::new();

        let dropped =
            GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);

        assert_eq!(dropped, 0.0);
        assert_eq!(out.len(), 2);
    }

    /// `out` is caller scratch: it is cleared, so a buffer carrying the previous sample's
    /// rows comes back holding this sample's and nothing else.
    #[test]
    fn the_scratch_buffer_holds_no_trace_of_the_previous_sample() {
        let previous = [supported_row(0, 4, 99, -1.0), supported_row(1, 4, 98, -2.0)];
        let current = [supported_row(0, 1, 11, -3.5)];
        let mut out = Vec::new();

        let _ =
            GenericObservation::fill_from_supported_alleles(&previous, &kept_in_order(2), &mut out);
        let _ =
            GenericObservation::fill_from_supported_alleles(&current, &kept_in_order(2), &mut out);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].num_reads, 11);
    }

    /// A mapping too short for the merge's table was built at another locus. Scoring those
    /// reads against whatever allele sat at that index is the failure this refuses, and it
    /// refuses it in every build profile — `[profile.release]` turns debug assertions off.
    #[test]
    #[should_panic(expected = "was built against a different locus's table")]
    fn a_mapping_that_does_not_cover_the_merges_table_is_a_caller_bug() {
        let rows = [supported_row(2, 1, 6, -8.25)];
        let mut out = Vec::new();

        let _ = GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);
    }

    /// The last allele the mapping covers is not itself out of range — a bound written with
    /// `<` where `<=` was meant, or the other way about, would reject it.
    #[test]
    fn the_last_allele_the_mapping_covers_is_accepted() {
        let rows = [supported_row(1, 1, 6, -8.25)];
        let mut out = Vec::new();

        let dropped =
            GenericObservation::fill_from_supported_alleles(&rows, &kept_in_order(2), &mut out);

        assert_eq!(dropped, 0.0);
        assert_eq!(out[0].allele, AlleleId(1));
    }

    /// The constructor carries all three parts through, **and the slices are compared whole
    /// rather than counted**: a length hides a row dropped and another duplicated, or two
    /// reordered, and the order is a determinism requirement rather than tidiness (spec §8).
    ///
    /// The pooled leftover is the part that can go missing without any genotype moving — it
    /// is the same number in every row of the sample — so nothing but a direct comparison
    /// would notice, and it is the term the data likelihood needs to be comparable between
    /// loci.
    #[test]
    fn the_constructor_keeps_the_pooled_leftover_and_both_slices_whole() {
        let supported = [
            GenericObservation::of_supported_allele(&supported_row(0, 1, 11, -3.5), AlleleId(0)),
            GenericObservation::of_supported_allele(&supported_row(1, 2, 4, -19.0), AlleleId(1)),
        ];
        let partials = [partial_row(6, -22.5), partial_row(2, -7.25)];

        let evidence = GenericSampleEvidence::new(&supported, -7.5, &partials);

        assert_eq!(evidence.unmatched_q_sum, -7.5);
        assert_eq!(evidence.supported, &supported[..]);
        assert_eq!(evidence.partials, &partials[..]);
    }

    #[test]
    fn a_sample_that_showed_nothing_carries_no_reads_and_no_pooled_leftover() {
        let evidence = GenericSampleEvidence::empty();

        assert_eq!(evidence.unmatched_q_sum, 0.0);
        assert!(evidence.supported.is_empty());
        assert!(evidence.partials.is_empty());
    }

    /// The merge's order is the view's order, and this constructor takes any slice. A call
    /// site that concatenated two samples' rows, or applied a selection mapping that
    /// reordered the alleles, would break the order with nothing saying so — and
    /// floating-point addition is not associative, so the run's genotypes move with it.
    #[test]
    #[should_panic(expected = "not ascending")]
    fn rows_out_of_pair_order_are_a_caller_bug() {
        let out_of_order = [
            GenericObservation::of_supported_allele(&supported_row(1, 0, 3, -1.0), AlleleId(1)),
            GenericObservation::of_supported_allele(&supported_row(0, 0, 3, -1.0), AlleleId(0)),
        ];

        let _ = GenericSampleEvidence::new(&out_of_order, 0.0, &[]);
    }

    /// The read group is the second half of the key, so rows ascending on the allele alone
    /// are not ascending — a check that compared only the allele would let them through.
    #[test]
    #[should_panic(expected = "not ascending")]
    fn rows_out_of_read_group_order_within_one_allele_are_a_caller_bug() {
        let out_of_order = [
            GenericObservation::of_supported_allele(&supported_row(2, 9, 3, -1.0), AlleleId(2)),
            GenericObservation::of_supported_allele(&supported_row(2, 4, 3, -1.0), AlleleId(2)),
        ];

        let _ = GenericSampleEvidence::new(&out_of_order, 0.0, &[]);
    }

    #[test]
    fn the_repeat_view_tells_complete_reads_from_reads_that_ran_out() {
        let observations = [
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 12),
            ssr_observation(b"ACAC", ran_out_after(4), 5),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        let complete: Vec<(usize, u32)> = evidence
            .complete_observations()
            .map(|(at, o)| (at, o.num_obs))
            .collect();
        let partial: Vec<(usize, u32)> = evidence
            .partial_observations()
            .map(|(at, o)| (at, o.num_obs))
            .collect();

        assert_eq!(complete, vec![(0, 12), (2, 3)]);
        assert_eq!(partial, vec![(1, 5)]);
        // The tract's own detail is the whole of what the emission is told about which
        // repeat it is scoring against, and nothing else here reads it back.
        assert_eq!(evidence.detail, &detail);
        assert_eq!(evidence.observations, &observations[..]);
    }

    /// The same bases seen by a read that spanned the tract and by one that ran out are
    /// **two** observations (spec §1.3), and only the witness separates them. A filter that
    /// read the bases instead — the plausible "is this shorter than the tract?" shortcut —
    /// puts both on the same side, and every other fixture here would still pass, because
    /// in each of them the partial's bases are also the shortest.
    #[test]
    fn two_observations_of_the_same_bases_split_on_the_witness_alone() {
        let observations = [
            ssr_observation(b"ACACAC", ReadWitness::Complete, 9),
            ssr_observation(b"ACACAC", ran_out_after(6), 4),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        assert_eq!(
            evidence
                .complete_observations()
                .map(|(_, o)| o.num_obs)
                .collect::<Vec<_>>(),
            vec![9]
        );
        assert_eq!(
            evidence
                .partial_observations()
                .map(|(_, o)| o.num_obs)
                .collect::<Vec<_>>(),
            vec![4]
        );
    }

    /// Every observation is on exactly one side, and the positions the two yield are the
    /// positions in the slice. That is what lets a complete term and a censored term be
    /// summed without double-counting a read or losing one, and what lets one cache serve
    /// both — and it is a property of the *pair*: each filter's own test passes just as well
    /// if the two overlap.
    #[test]
    fn the_two_repeat_filters_partition_every_observation_and_agree_on_its_position() {
        let observations = [
            ssr_observation(b"ACACACAC", ReadWitness::Complete, 12),
            ssr_observation(b"ACAC", ran_out_after(4), 5),
            ssr_observation(b"ACACAC", ReadWitness::Complete, 3),
        ];
        let detail = ssr_detail();
        let evidence = SsrSampleEvidence::new(&observations, &detail);

        let mut seen: Vec<usize> = evidence
            .complete_observations()
            .chain(evidence.partial_observations())
            .map(|(at, _)| at)
            .collect();
        seen.sort_unstable();

        assert_eq!(seen, vec![0, 1, 2]);
        for (at, observation) in evidence
            .complete_observations()
            .chain(evidence.partial_observations())
        {
            assert_eq!(observation.num_obs, observations[at].num_obs);
        }
    }
}
