//! The STR starting point: the population's own fitted **length spectrum** at this tract's
//! stratum, read onto the locus's candidate lengths.
//!
//! Everything the SNP/indel starting point leans on is false at a repeat tract, and the three
//! failures are separate (`doc/devel/ng/spec/calling_priors.md` §5). The reference accession's
//! length is one draw among several common ones rather than the usual winner, so the reference
//! allele carries no presumption. The alleles are **ordered** — a tract of 11 repeats is
//! adjacent to one of 10 and far from one of 4 — so splitting mass evenly across alternatives
//! throws away the only structure that makes a rare long allele believable. And repeat tracts
//! mutate orders of magnitude faster than bases do, so how variable they are is a separate
//! number, measured separately.
//!
//! ## The two questions, and both are answered by one fitted pair
//!
//! **Where the mass sits** is the stratum's fitted length spectrum: how that stratum's
//! chromosomes are spread over tract lengths, indexed in whole repeat units either side of the
//! **reference** tract length.
//!
//! **How much mass there is** is the concentration fitted beside it: how monomorphic the
//! stratum's tracts are, in chromosomes' worth of belief.
//!
//! Together they are a Dirichlet over a tract's length frequencies, estimated from that
//! stratum's own tracts by the joint repeat fit
//! (`parameter_estimation::joint::ssr_fit::StratumFit`), conditioning on each sample's
//! homozygote excess as it goes so that inbreeding is divided out inside the estimator.
//!
//! ## What this replaced, and why — `population_diversity.md` §4.2
//!
//! Until 2026-08-26 this module built the shape rather than reading it: mass falling off
//! geometrically from the cohort's commonest length at the tract, scaled so that the prior's own
//! implied gene diversity reproduced a separately measured cohort-wide repeat gene diversity.
//! Both halves are gone, and the four reasons are facts about the two constructions:
//!
//! - **The constructed shape had one free parameter and the fitted one has none.** The decay was
//!   a single number per group of loci with a coded fallback of 0.5.
//! - **The constructed version had a failure mode the fitted one does not have.** Scaling a
//!   shape to reproduce a measured diversity is only possible below a ceiling the shape itself
//!   sets, and above it the builder refused. **At one outbred sample that was every tract** — a
//!   single diploid shows at most three lengths, whose shape can imply at most 0.625, against
//!   the ~0.72 repeat diversity HG002 actually has. A fitted Dirichlet asserts no such scaling
//!   and cannot fail this way, which is why there is no refusal left in this module.
//! - **It removed a per-locus input that has no source**: the cohort's commonest length *at this
//!   tract*, which would come from repeat-tract candidate selection, unwritten. The fitted
//!   spectrum is indexed from the reference tract length, which every locus already knows.
//! - **It removed the need for a cohort-wide repeat diversity number entirely**, which nothing
//!   emits.
//!
//! ## Two words this project keeps apart, and they are both "spectrum"
//!
//! A **length spectrum** is this path's: how a stratum's chromosomes are spread over tract
//! **lengths**. A **frequency spectrum** is the ordinary-site path's: how allele **frequencies**
//! are spread across the population, which [`seed_generic`](super::seed_generic) reads. Neither
//! is ever called just "the spectrum" (`population_diversity.md` §2).

use crate::genetics::MIN_ALT_CONCENTRATION;
use crate::ng::parameter_estimation::joint::stratum_fits::LengthSpectrum;

use super::Concentration;

/// The floor on a candidate's share of the prior's mass, so a candidate the fit gives nothing
/// keeps a strictly positive share rather than falling into an absorbing zero.
///
/// **Ported from production's `G0_FLOOR`** (`src/ssr/cohort/allele_freq_prior.rs`). Why a zero
/// would matter is the spec's sentence rather than production's: a masked long heterozygous copy
/// that the candidate set nearly missed has to stay recoverable rather than fall into a prior it
/// can never climb out of (`calling_priors.md` §5).
///
/// **Two candidates reach it, and they are different situations.** One is a candidate *outside*
/// the fitted spectrum's reach — further from the reference tract length than the fit was ever
/// allowed to place mass — which is outside everything the fit saw rather than something it
/// measured as rare. The other is a candidate inside the reach at a length class the fit put no
/// mass on. Both take the floor, and the floor is what makes the shape's total strictly
/// positive, so the shares are always a distribution
/// (`tests::a_candidate_the_fit_reaches_and_gives_nothing_still_keeps_a_share`).
///
/// It has production's value, which happens to equal [`MIN_ALT_CONCENTRATION`], and **the two
/// are not the same quantity**: this floors a dimensionless share of the shape, that floors a
/// count of chromosomes.
const SHAPE_FLOOR: f64 = 1e-12;

/// Fill `out` with the STR seed: the stratum's fitted length spectrum read onto this locus's
/// candidate lengths, at the concentration it was fitted with.
///
/// ```text
/// offset:  d_j = repeat count of j − reference_repeat_count
/// shape:   w_j = max(spectrum(d_j), SHAPE_FLOOR)          (1/K on the stated-flat rung)
/// seed:    α_j = concentration · w_j,  floored at MIN_ALT_CONCENTRATION
/// ```
///
/// ## The scaling is the Dirichlet's own conditioning, not a renormalisation
///
/// A Dirichlet over the fit's `2·span + 1` length classes, **conditioned on the tract carrying
/// one of this locus's candidate lengths**, is exactly the Dirichlet over those candidates with
/// each class's own `α` kept — that is the distribution's marginal-conditional property, and it
/// is why nothing here divides by the retained mass. Two consequences worth stating because
/// neither is an accident:
///
/// - **The shape is right by construction.** Normalising the seed gives back the fitted spectrum
///   restricted to the candidate lengths, which is what
///   `population_diversity.md` §8's first check asks for
///   (`tests::the_seed_normalises_back_to_the_fitted_spectrum_over_the_candidates`).
/// - **The total is `concentration × (mass the candidates cover)`, which is smaller than the
///   concentration wherever the candidate set misses fitted mass.** A locus whose candidates
///   cover a tenth of what the stratum spreads over is held with a tenth of the conviction, and
///   that is the honest reading: the fit expected lengths this locus is not being called over
///   (`tests::a_candidate_set_covering_less_fitted_mass_is_held_with_less_conviction`).
///
/// ## `candidate_repeat_counts` and `reference_repeat_count`
///
/// `candidate_repeat_counts` is parallel to the locus's
/// [`CandidateAlleles`](crate::ng::calling::CandidateAlleles), so entry 0 is the reference
/// allele's — but **nothing here privileges entry 0**, which is the whole difference from the
/// SNP/indel path: at a repeat tract the reference length is one common length among several
/// (`calling_priors.md` §5).
///
/// `reference_repeat_count` is the **tract's own** reference repeat count, the same number that
/// picked the stratum this spectrum came from. It is what the spectrum's offsets are measured
/// from, and it is not any candidate's — the reference allele usually spells it, but a locus
/// whose reference allele was cut by selection still has one.
///
/// ## Three things a caller can get wrong that nothing here can catch
///
/// The lengths of the two slices are checkable and are checked. These are not:
///
/// - **a candidate's repeat count passed as the reference tract's.** The shape re-centres on
///   that candidate, which is the flattening
///   [`StratumFits::length_spectrum_at`](crate::ng::parameter_estimation::joint::stratum_fits::StratumFits::length_spectrum_at)'s
///   own documentation warns about; measured as a difference within the prior row in
///   `tests::the_reference_repeat_count_is_the_tracts_and_nothing_here_can_check_it`.
/// - **another stratum's spectrum.** Every fitted spectrum has the same shape and the same
///   length, so one stratum's reads as another's with nothing to say so.
/// - **the buffer reused across loci with the call skipped.** The previous locus's row, entry
///   for entry.
///
/// ## Two alleles of the same length land on one class
///
/// The seed is keyed by repeat count, so two candidates that differ by an interior
/// substitution — an interrupted repeat — sit at the same offset and **each receives that
/// class's full weight**, which is production's behaviour. Whether to divide the class instead
/// needs the interrupted-repeat work to say how it should be weighted; the signature takes
/// counts rather than sequences precisely so that change lands in this one function
/// (`calling_priors.md` §5.2, open as that spec's Q3).
///
/// ## The ends of the range
///
/// **At one candidate allele** the locus has one genotype, whose prior probability is 1 at any
/// positive concentration, so what this returns there cannot be wrong — it is
/// `concentration × w`, with no special case
/// (`tests::a_locus_with_one_candidate_length_is_seeded_from_that_length_alone`).
///
/// **At one sample nothing special happens**, which is the change this rewrite bought: the
/// stratum's spectrum is fitted across *tracts*, and a single genome carries the same tracts a
/// panel does (`population_diversity.md` §4.3). The construction this replaced refused every
/// tract there (`tests::a_single_diploids_candidate_sets_are_all_seeded`).
///
/// ## Shape and cost
///
/// Fills the caller's buffer and **allocates nothing**. Three passes over the locus's candidate
/// lengths — write the weights, sum them, scale — and no `lgamma`, against one `lgamma` per
/// allele a genotype carries a copy of in the prior row this feeds. (The sum is the shared
/// filler's, and `fill_ssr_seed` discards it: the seed scales the *unnormalised* weights, which
/// is §4.2's conditioning argument above.)
///
/// # Panics
///
/// **Both slice lengths are checked in release**: `out` is the caller's and is reused across
/// loci, so a short one would leave the previous locus's entries standing in this locus's seed.
#[must_use]
pub fn fill_ssr_seed<'a>(
    candidate_repeat_counts: &[u32],
    reference_repeat_count: u32,
    spectrum: LengthSpectrum<'_>,
    out: &'a mut [f64],
) -> Concentration<'a> {
    let _ = fill_candidate_weights(
        candidate_repeat_counts,
        reference_repeat_count,
        spectrum,
        out,
    );
    let concentration = spectrum.concentration();
    for slot in out.iter_mut() {
        *slot = (concentration * *slot).max(MIN_ALT_CONCENTRATION);
    }
    Concentration::new(out)
}

/// Fill `out` with the share of the prior's starting mass each **candidate** at this locus
/// carries — one entry per candidate, in the candidate set's own order, summing to 1.
///
/// This is the same shape [`fill_ssr_seed`] scales into a concentration, before it is scaled:
/// where a length sits in the stratum's fitted spread, with no claim yet about how much
/// conviction there is. One implementation stands behind both, so the prior's belief about
/// lengths cannot drift between its two consumers
/// (`doc/devel/ng/arch/calling_priors.md` §5).
///
/// **The consumer is the read likelihood's contamination term.** Scoring a read against a
/// candidate allele mixes three sources: the sample's own copies, a uniform outlier share, and a
/// share attributable to DNA from another sample — and that third term needs a distribution over
/// the lengths contaminating DNA might carry (`doc/devel/ng/arch/read_likelihoods.md` §4.1,
/// `SsrContamination::length_distribution`). Absent a measured one, the prior's own belief about
/// which lengths are common at this tract is the stand-in, and it is a stand-in rather than a
/// measurement: it assumes the contaminant was drawn from the same population as the cohort,
/// which is what the sibling spec's §4.5.1 says it assumes.
///
/// ## `OPEN:` this is per **candidate**, and the term it feeds is per observed **length**
///
/// The mixture's third term is `c · seed(o)`, where `o` is an observation — a read. Candidates
/// and observed lengths coincide only when every candidate is a distinct length and every read
/// lands on one of them, and three cases break that. **None is this function's to settle**; they
/// are recorded here so the likelihood step meets them as a decision rather than as a surprise:
///
/// - **Two candidates of one length each take that class's full share.** Read as a claim about
///   *lengths* it double-counts
///   (`tests::two_spellings_of_one_length_carry_that_lengths_share_twice`).
/// - **A read at a length that is not a candidate has no entry.** The candidate set is
///   post-prune, while the mixture's sibling uniform term is spread over every length the stutter
///   model can reach from a candidate — a strictly larger support (sibling spec §4.5).
/// - **A censored read carries no length at all**, only a lower bound.
///
/// ## Cost, and when a caller should not use this
///
/// Two passes over the candidate lengths — write the weights, normalise — built **once per
/// locus** rather than once per scoring call. But a loop that has just called [`fill_ssr_seed`]
/// at this locus already holds the shape and should not rebuild it: it is that concentration
/// divided by its own total. This export is for a caller that wants the shape and not the seed.
///
/// **It does not survive a candidate being added.** A discovery round appends candidates mid
/// locus (`doc/devel/ng/arch/calling_em_loop.md` §5), and a frozen buffer is then one entry short
/// with nothing to raise, because by construction it is not refilled. Discovery is off by
/// default; a loop that turns it on has to rebuild this.
///
/// ## Two things a caller can get wrong that nothing here can catch
///
/// - **A candidate's repeat count passed as the reference tract's**, which re-centres the shape.
/// - **A buffer reused across loci with the call skipped.** The previous locus's shares.
///
/// **Shares, not chromosome counts.** The entries sum to 1 and carry no conviction; nothing here
/// may be handed where a [`Concentration`] belongs, which is why this returns nothing rather than
/// wrapping the buffer.
///
/// # Panics
///
/// As [`fill_ssr_seed`], on an empty candidate set or a buffer that is not one entry per
/// candidate.
pub fn fill_seed_share_per_candidate(
    candidate_repeat_counts: &[u32],
    reference_repeat_count: u32,
    spectrum: LengthSpectrum<'_>,
    out: &mut [f64],
) {
    let total = fill_candidate_weights(
        candidate_repeat_counts,
        reference_repeat_count,
        spectrum,
        out,
    );
    for slot in out.iter_mut() {
        *slot /= total;
    }
}

/// Fill `out` with each candidate's **unnormalised** weight under the spectrum, and return their
/// total.
///
/// **Unnormalised is what the seed wants and normalised is what the share wants**, and the
/// difference between the two is the whole of `population_diversity.md`'s conditioning argument
/// — see [`fill_ssr_seed`]. Returning the total rather than recomputing it keeps the share's
/// second pass from walking the buffer a third time.
///
/// **The stated-flat rung writes `1/K` rather than `1`**, so that the seed's total on that rung
/// is the stated concentration itself rather than the concentration times the candidate count.
/// A flat shape asserts no belief about which length is likelier; it must not also assert more
/// conviction at a locus with more candidates
/// (`tests::the_stated_flat_rung_holds_the_locus_at_the_stated_concentration_whatever_the_candidate_count`).
fn fill_candidate_weights(
    candidate_repeat_counts: &[u32],
    reference_repeat_count: u32,
    spectrum: LengthSpectrum<'_>,
    out: &mut [f64],
) -> f64 {
    assert_lengths(candidate_repeat_counts, out);
    match spectrum.fitted_weights() {
        None => {
            let share = 1.0 / candidate_repeat_counts.len() as f64;
            out.fill(share);
        }
        Some(weights) => {
            let span = i64::from(
                spectrum
                    .allele_span()
                    .expect("a fitted spectrum has a span"),
            );
            for (slot, &repeat_count) in out.iter_mut().zip(candidate_repeat_counts) {
                let offset = i64::from(repeat_count) - i64::from(reference_repeat_count);
                // **Outside the fit's reach takes the floor, not the end class's weight.** The
                // end class holds the mass the fit put at exactly `±span` repeats; handing it to
                // a candidate twelve repeats away would assert that a twelve-repeat departure is
                // as common as a six-repeat one, which the fit never said and its own span is
                // the statement that it declined to say
                // (`tests::a_candidate_past_the_fits_reach_takes_the_floor_and_not_the_end_class`).
                let weight = if offset.abs() <= span {
                    weights[(offset + span) as usize]
                } else {
                    0.0
                };
                *slot = weight.max(SHAPE_FLOOR);
            }
        }
    }
    out.iter().sum()
}

/// The one shape check both fillers make, in release.
fn assert_lengths(candidate_repeat_counts: &[u32], out: &[f64]) {
    assert!(
        !candidate_repeat_counts.is_empty(),
        "every repeat tract has a reference allele, so its candidate set has at least one \
         length — the caller has lost track of which locus it is on"
    );
    assert_eq!(
        out.len(),
        candidate_repeat_counts.len(),
        "the buffer must cover the locus's candidate lengths exactly: a longer one normalises \
         the shape against entries that are not candidates and a shorter one leaves the \
         previous locus's entries behind, and both look like answers"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::joint::stratum_fits::{
        FittedFrom, STATED_FLAT_CONCENTRATION,
    };

    /// **A five-class spectrum leaning long, and every class a different number.** Two
    /// properties are deliberate and both kill a class of mutation: no two classes share a
    /// weight, so reading one class for another shows up; and the two ends differ (`0.04`
    /// against `0.05`), so an offset whose sign is flipped is a different answer rather than
    /// the same one by symmetry.
    ///
    /// Indexed `-2 ..= +2` in whole repeat units from the reference tract length.
    fn leaning_long() -> Vec<f64> {
        vec![0.04, 0.16, 0.50, 0.25, 0.05]
    }

    /// The tract this module's fixtures sit at: 11 whole repeats in the reference.
    const REFERENCE_REPEATS: u32 = 11;

    /// The candidate lengths, in an order no sort would produce and at offsets that are not
    /// symmetric about the reference: `0, +1, −2, +2`.
    fn candidates() -> Vec<u32> {
        vec![11, 12, 9, 13]
    }

    fn own_fit(weights: &[f64], concentration: f64) -> LengthSpectrum<'_> {
        LengthSpectrum::fitted(weights, concentration, FittedFrom::ThisStratum)
    }

    fn seed(
        candidate_repeat_counts: &[u32],
        reference: u32,
        spectrum: LengthSpectrum<'_>,
        out: &mut [f64],
    ) -> Vec<f64> {
        fill_ssr_seed(candidate_repeat_counts, reference, spectrum, out)
            .get()
            .to_vec()
    }

    fn shares(
        candidate_repeat_counts: &[u32],
        reference: u32,
        spectrum: LengthSpectrum<'_>,
    ) -> Vec<f64> {
        let mut out = vec![0.0; candidate_repeat_counts.len()];
        fill_seed_share_per_candidate(candidate_repeat_counts, reference, spectrum, &mut out);
        out
    }

    /// **Relative, with no absolute floor**, and the floor is exactly what went wrong when
    /// there was one. Written `1e-12 * right.abs().max(1.0)`, the tolerance never fell below
    /// `1e-12` — the size of [`MIN_ALT_CONCENTRATION`] itself — so
    /// `close(x, MIN_ALT_CONCENTRATION)` accepted anything in `[0, 2e-12]`, zero included, and
    /// the test that exists to pin the concentration floor passed with the floor deleted.
    ///
    /// An exact zero is compared exactly, because a relative tolerance around zero is no
    /// tolerance at all.
    fn close(left: f64, right: f64) -> bool {
        if right == 0.0 {
            return left == 0.0;
        }
        (left - right).abs() <= 1e-9 * right.abs()
    }

    /// **The whole of what this builder does, in one assertion**: each candidate's starting
    /// chromosomes are the concentration times that candidate's own share of the fitted length
    /// spectrum, at the offset of its repeat count from the **tract's** reference length.
    ///
    /// **Three mutations die here**, which is why the fixture is shaped the way it is. Reading
    /// the offset with its sign flipped swaps `9`'s `0.04` for `13`'s `0.05`. Reading the class
    /// at `offset` rather than `offset + span` shifts every candidate two classes. And
    /// normalising the shares before scaling gives a total of the concentration rather than of
    /// `concentration × 0.84`.
    #[test]
    fn the_seed_is_the_fitted_spectrum_at_the_fitted_concentration() {
        let weights = leaning_long();
        let mut out = vec![0.0; 4];
        let filled = seed(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );

        // offsets 0, +1, −2, +2 → weights 0.50, 0.25, 0.04, 0.05, each × 12 chromosomes.
        for (slot, (candidate, expected)) in filled
            .iter()
            .zip(candidates().iter().zip([6.0, 3.0, 0.48, 0.60]))
        {
            assert!(
                close(*slot, expected),
                "the candidate at {candidate} repeats starts with {slot} chromosomes where the \
                 spectrum's share at its offset from {REFERENCE_REPEATS}, times 12, is {expected}"
            );
        }
    }

    /// **The middle rung is seeded from its own weights, not treated as flat** — the rung that
    /// is the whole new fitting cost of this step, and the one no test pinned a number for.
    ///
    /// Measured, and it is why this test exists: with the pooled rung routed to the flat arm,
    /// all seventeen of this module's tests passed. The one fixture that ever handed a pooled
    /// spectrum was `the_shared_export_is_the_seed_divided_by_its_total`, whose only assertion —
    /// share equals seed entry over seed total — is an identity that holds for any shape,
    /// flat included.
    #[test]
    fn the_pooled_rung_is_seeded_from_its_own_weights() {
        let weights = leaning_long();
        let mut out = vec![0.0; 4];
        let filled = seed(
            &candidates(),
            REFERENCE_REPEATS,
            LengthSpectrum::fitted(&weights, 12.0, FittedFrom::ItsPeriodsPooledTracts),
            &mut out,
        );

        for (slot, (candidate, expected)) in filled
            .iter()
            .zip(candidates().iter().zip([6.0, 3.0, 0.48, 0.60]))
        {
            assert!(
                close(*slot, expected),
                "on the pooled rung the candidate at {candidate} repeats starts with {slot} \
                 chromosomes where its class's share times 12 is {expected}; a flat shape would \
                 give every candidate {}",
                12.0 / 4.0
            );
        }
    }

    /// **A buffer *longer* than the candidate set is refused too**, and it is the half the
    /// message names that no fixture reached: both mis-sized-buffer tests passed a shorter one.
    ///
    /// It is not a harmless overshoot. `fill_candidate_weights` writes one entry per candidate
    /// and then sums the **whole** buffer, so the previous locus's trailing entries enter both
    /// the shares and the seed's total.
    #[test]
    #[should_panic(expected = "must cover the locus's candidate lengths exactly")]
    fn a_buffer_longer_than_the_candidate_set_is_refused() {
        let mut out = vec![0.0; 5];
        let weights = leaning_long();
        let _ = fill_ssr_seed(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );
    }

    /// **The shape floor is production's `G0_FLOOR`, and this pins its value** — every other
    /// assertion about it is written `12.0 * SHAPE_FLOOR`, which tracks the constant and would
    /// survive any change to it.
    #[test]
    fn the_shape_floor_is_the_value_production_uses() {
        assert_eq!(SHAPE_FLOOR, 1e-12, "production's `G0_FLOOR`");
    }

    /// **`population_diversity.md` §8's first check**: normalised, the seed is the fitted
    /// spectrum restricted to the lengths this locus is being called over. Nothing about the
    /// shape is invented by the mapping.
    #[test]
    fn the_seed_normalises_back_to_the_fitted_spectrum_over_the_candidates() {
        let weights = leaning_long();
        let mut out = vec![0.0; 4];
        let filled = seed(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );

        let total: f64 = filled.iter().sum();
        let covered = 0.50 + 0.25 + 0.04 + 0.05;
        for (slot, fitted) in filled.iter().zip([0.50, 0.25, 0.04, 0.05]) {
            assert!(
                close(slot / total, fitted / covered),
                "normalised, the seed gives {} where the fitted spectrum restricted to these \
                 four lengths gives {}",
                slot / total,
                fitted / covered
            );
        }
    }

    /// **A candidate set that covers less of the fitted spread is held with less conviction**,
    /// and that is the Dirichlet's own conditioning rather than a policy: the total is the
    /// concentration times the mass the candidates cover.
    #[test]
    fn a_candidate_set_covering_less_fitted_mass_is_held_with_less_conviction() {
        let weights = leaning_long();
        let mut wide = vec![0.0; 4];
        let wide = seed(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut wide,
        );
        let mut narrow = vec![0.0; 1];
        let narrow = seed(
            &[REFERENCE_REPEATS],
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut narrow,
        );

        let wide_total: f64 = wide.iter().sum();
        let narrow_total: f64 = narrow.iter().sum();
        assert!(
            close(wide_total, 12.0 * 0.84),
            "four candidates covering 0.84 of the fitted mass are held with {wide_total} \
             chromosomes"
        );
        assert!(
            close(narrow_total, 12.0 * 0.50),
            "one candidate covering 0.50 of it is held with {narrow_total}"
        );
    }

    /// **Outside the fit's reach takes the floor, not the end class's weight.** The end class
    /// holds what the fit put at exactly two repeats out; handing it to a candidate five repeats
    /// out would claim the fit measured something it declined to measure.
    #[test]
    fn a_candidate_past_the_fits_reach_takes_the_floor_and_not_the_end_class() {
        let weights = leaning_long();
        let mut out = vec![0.0; 2];
        let filled = seed(
            &[REFERENCE_REPEATS, REFERENCE_REPEATS + 5],
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );

        assert!(close(filled[0], 6.0));
        assert!(
            close(filled[1], 12.0 * SHAPE_FLOOR),
            "the candidate five repeats out starts with {} chromosomes; the floor puts it at \
             {}, and the end class at +2 repeats would put it at {}",
            filled[1],
            12.0 * SHAPE_FLOOR,
            12.0 * 0.05
        );
    }

    /// **The span is read off the spectrum's own class count**, so the same candidate is outside
    /// the reach of a three-class fit and inside a five-class one.
    #[test]
    fn how_far_the_spectrum_reaches_is_read_off_its_own_class_count() {
        let narrow = vec![0.2, 0.5, 0.3];
        let wide = leaning_long();
        let candidate = [REFERENCE_REPEATS + 2];

        let mut out = vec![0.0; 1];
        let at_span_one = seed(
            &candidate,
            REFERENCE_REPEATS,
            own_fit(&narrow, 12.0),
            &mut out,
        );
        let mut out = vec![0.0; 1];
        let at_span_two = seed(
            &candidate,
            REFERENCE_REPEATS,
            own_fit(&wide, 12.0),
            &mut out,
        );

        assert!(
            close(at_span_one[0], 12.0 * SHAPE_FLOOR),
            "two repeats out is past a three-class spectrum's reach, and it gave {}",
            at_span_one[0]
        );
        assert!(
            close(at_span_two[0], 12.0 * 0.05),
            "two repeats out is the end class of a five-class spectrum, and it gave {}",
            at_span_two[0]
        );
    }

    /// **A length the fit reached and put no mass on keeps a share**, so a masked long copy the
    /// candidate set nearly missed stays recoverable rather than falling into a prior it can
    /// never climb out of.
    #[test]
    fn a_candidate_the_fit_reaches_and_gives_nothing_still_keeps_a_share() {
        let weights = vec![0.04, 0.0, 0.66, 0.25, 0.05];
        let mut out = vec![0.0; 2];
        let filled = seed(
            &[REFERENCE_REPEATS, REFERENCE_REPEATS - 1],
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );

        assert!(
            filled[1] > 0.0,
            "the candidate one repeat short sits on a class the fit gave 0.0, and it must not \
             be an absorbing zero; it is {}",
            filled[1]
        );
        assert!(
            close(filled[1], 12.0 * SHAPE_FLOOR),
            "the empty class takes the shape's floor times the concentration, {}, and gave {}",
            12.0 * SHAPE_FLOOR,
            filled[1]
        );
    }

    /// **The concentration floor is what makes the shape floor safe.** At a weak concentration
    /// the shape's floor scales below `MIN_ALT_CONCENTRATION`, and an entry under it is what
    /// [`Concentration::new`] refuses in debug — so this is the check, and without it the
    /// builder panics rather than returning a small number.
    #[test]
    fn a_far_candidate_at_a_weak_concentration_still_clears_the_concentration_floor() {
        let weights = leaning_long();
        let mut out = vec![0.0; 2];
        let filled = seed(
            &[REFERENCE_REPEATS, REFERENCE_REPEATS + 5],
            REFERENCE_REPEATS,
            own_fit(&weights, 1e-3),
            &mut out,
        );

        // **The fixture only tests the floor if the unfloored value is below it**: `1e-3 ×
        // SHAPE_FLOOR` is `1e-15`, against a `MIN_ALT_CONCENTRATION` of `1e-12`. A const block,
        // because both operands are constants and the check is about the fixture rather than
        // about the run.
        const {
            assert!(1e-3 * SHAPE_FLOOR < MIN_ALT_CONCENTRATION);
        }
        assert_eq!(
            filled[1], MIN_ALT_CONCENTRATION,
            "the floor is exactly representable, so this is an exact comparison — through a \
             tolerance the size of the floor itself it would accept zero"
        );
    }

    /// **The bottom rung holds the locus at the stated concentration whatever the candidate
    /// count**, because a flat shape asserts no belief about which length is likelier and must
    /// not assert more conviction at a locus with more candidates.
    #[test]
    fn the_stated_flat_rung_holds_the_locus_at_the_stated_concentration_whatever_the_candidate_count()
     {
        for count in [1_usize, 3, 7] {
            let counts: Vec<u32> = (0..count as u32).map(|i| REFERENCE_REPEATS + i).collect();
            let mut out = vec![0.0; count];
            let filled = seed(
                &counts,
                REFERENCE_REPEATS,
                LengthSpectrum::stated_flat(2.5),
                &mut out,
            );

            let total: f64 = filled.iter().sum();
            assert!(
                close(total, 2.5),
                "{count} candidates on the stated-flat rung are held with {total} chromosomes, \
                 where the rung states 2.5"
            );
            for slot in &filled {
                assert!(close(*slot, 2.5 / count as f64));
            }
        }
    }

    /// A locus with one candidate length has one genotype, whose prior probability is 1 at any
    /// positive concentration — so there is no special case, and this pins what the general rule
    /// gives there.
    #[test]
    fn a_locus_with_one_candidate_length_is_seeded_from_that_length_alone() {
        let weights = leaning_long();
        let mut out = vec![0.0; 1];
        let filled = seed(
            &[REFERENCE_REPEATS],
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );
        assert!(close(filled[0], 6.0));
    }

    /// **The failure this rewrite removed.** A single outbred diploid shows at most three
    /// lengths at a tract, and the construction this replaced refused **every** such locus: it
    /// scaled a geometric shape to reproduce a measured repeat diversity, which is possible only
    /// below a ceiling the shape sets — at most 0.625 over three lengths at the coded fallback
    /// decay, against the ~0.72 HG002 actually has (`population_diversity.md` §4.2).
    ///
    /// The fitted pair asserts no such scaling, so each of those candidate sets is seeded, and
    /// seeded from the stratum's own shape rather than from a shape chosen here.
    #[test]
    fn a_single_diploids_candidate_sets_are_all_seeded() {
        let weights = leaning_long();
        for lengths in [
            &[REFERENCE_REPEATS][..],
            &[REFERENCE_REPEATS, REFERENCE_REPEATS + 1][..],
            &[
                REFERENCE_REPEATS - 2,
                REFERENCE_REPEATS,
                REFERENCE_REPEATS + 1,
            ][..],
        ] {
            let mut out = vec![0.0; lengths.len()];
            let filled = seed(
                lengths,
                REFERENCE_REPEATS,
                own_fit(&weights, 12.0),
                &mut out,
            );
            let total: f64 = filled.iter().sum();

            for (slot, &length) in filled.iter().zip(lengths) {
                let offset = i64::from(length) - i64::from(REFERENCE_REPEATS);
                let fitted = weights[(offset + 2) as usize];
                assert!(
                    close(
                        slot / total,
                        fitted
                            / lengths
                                .iter()
                                .map(|&l| weights
                                    [(i64::from(l) - i64::from(REFERENCE_REPEATS) + 2) as usize])
                                .sum::<f64>()
                    ),
                    "a candidate set of {} lengths is seeded from the stratum's own shape",
                    lengths.len()
                );
            }
        }
    }

    /// **The shared export is the seed divided by its own total**, on every rung — one
    /// implementation stands behind both, so the prior's belief about lengths cannot drift
    /// between its two consumers.
    #[test]
    fn the_shared_export_is_the_seed_divided_by_its_total() {
        let weights = leaning_long();
        for spectrum in [
            own_fit(&weights, 12.0),
            LengthSpectrum::fitted(&weights, 3.0, FittedFrom::ItsPeriodsPooledTracts),
            LengthSpectrum::stated_flat(2.5),
        ] {
            let mut out = vec![0.0; 4];
            let filled = seed(&candidates(), REFERENCE_REPEATS, spectrum, &mut out);
            let total: f64 = filled.iter().sum();
            let exported = shares(&candidates(), REFERENCE_REPEATS, spectrum);

            let sum: f64 = exported.iter().sum();
            assert!(close(sum, 1.0), "the shares sum to {sum}");
            for (share, slot) in exported.iter().zip(&filled) {
                assert!(
                    close(*share, slot / total),
                    "the share {share} against the seed entry {slot} over its total {total}"
                );
            }
        }
    }

    /// **Two spellings of one length carry that length's share twice**, which is production's
    /// behaviour as a concentration and is what the read likelihood's contamination term has to
    /// meet as a decision rather than as a surprise.
    #[test]
    fn two_spellings_of_one_length_carry_that_lengths_share_twice() {
        let weights = leaning_long();
        // The reference length twice — an interrupted repeat spells the same count — plus one
        // repeat above it.
        let exported = shares(
            &[REFERENCE_REPEATS, REFERENCE_REPEATS, REFERENCE_REPEATS + 1],
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
        );

        let at_the_reference_length = exported[0] + exported[1];
        let geometry = 0.50 / (0.50 + 0.25);
        assert!(
            close(
                at_the_reference_length,
                (0.50 + 0.50) / (0.50 + 0.50 + 0.25)
            ),
            "the reference length arrives with {at_the_reference_length} of the mass where the \
             spectrum's own two-class ratio gives it {geometry}"
        );
    }

    /// **A candidate's repeat count passed where the tract's belongs re-centres the shape on
    /// that candidate**, and nothing here can catch it — the measured size of the mistake is
    /// what this test records.
    ///
    /// Measured on the four candidates above, with `9` passed as the reference instead of `11`:
    /// the tract's own reference length falls from **0.595 of the prior's mass to 0.091**, a
    /// factor of 6.5, while the candidate at 9 rises from **0.048 to 0.909**, a factor of 19.
    #[test]
    fn the_reference_repeat_count_is_the_tracts_and_nothing_here_can_check_it() {
        let weights = leaning_long();
        let right = shares(&candidates(), REFERENCE_REPEATS, own_fit(&weights, 12.0));
        let wrong = shares(&candidates(), 9, own_fit(&weights, 12.0));

        assert!(
            (right[0] - 0.5952).abs() < 5e-4 && (wrong[0] - 0.0909).abs() < 5e-4,
            "the tract's own length holds {} of the mass and {} when a candidate's count is \
             passed instead",
            right[0],
            wrong[0]
        );
        assert!(
            (right[2] - 0.0476).abs() < 5e-4 && (wrong[2] - 0.9091).abs() < 5e-4,
            "the candidate at 9 repeats holds {} of the mass and {} when its own count is \
             passed as the tract's",
            right[2],
            wrong[2]
        );
    }

    #[test]
    #[should_panic(expected = "its candidate set has at least one")]
    fn an_empty_candidate_set_is_refused() {
        let mut out: [f64; 0] = [];
        let weights = leaning_long();
        let _ = fill_ssr_seed(&[], REFERENCE_REPEATS, own_fit(&weights, 12.0), &mut out);
    }

    #[test]
    #[should_panic(expected = "must cover the locus's candidate lengths exactly")]
    fn a_mis_sized_buffer_is_refused() {
        let mut out = vec![0.0; 3];
        let weights = leaning_long();
        let _ = fill_ssr_seed(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "must cover the locus's candidate lengths exactly")]
    fn a_mis_sized_buffer_is_refused_by_the_shared_export_too() {
        let mut out = vec![0.0; 3];
        let weights = leaning_long();
        fill_seed_share_per_candidate(
            &candidates(),
            REFERENCE_REPEATS,
            own_fit(&weights, 12.0),
            &mut out,
        );
    }

    /// The bottom rung's stated strength is one chromosome's worth of belief, which is what
    /// lets the reads move the prior from the first read onward.
    #[test]
    fn the_stated_flat_concentration_is_one_chromosomes_worth() {
        assert_eq!(STATED_FLAT_CONCENTRATION, 1.0);
    }
}
