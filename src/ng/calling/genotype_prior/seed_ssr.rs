//! The STR starting point: mass falling off geometrically from the cohort's modal repeat
//! count, totalling what the cohort's measured repeat diversity implies.
//!
//! Everything the SNP/indel starting point leans on is false at a repeat tract, and the three
//! failures are separate (`doc/devel/ng/spec/calling_priors.md` §5). The reference accession's
//! length is one draw among several common ones rather than the usual winner, so the reference
//! allele carries no presumption. The alleles are **ordered** — a tract of 11 repeats is
//! adjacent to one of 10 and far from one of 4 — so splitting mass evenly across alternatives
//! throws away the only structure that makes a rare long allele believable. And repeat tracts
//! mutate orders of magnitude faster than bases do, so how variable they are is a separate
//! number the pre-pass measures separately.
//!
//! ## Two questions, two parameters, and this is where ng departs from production
//!
//! **Where the mass sits** is production's `G₀` — its name for the pseudocount prior it puts
//! over a tract's candidate lengths — ported unchanged: a candidate's weight falls by a fixed
//! factor for every repeat unit it sits away from the cohort's mode at this locus
//! ([`fill_seed_shape`], from `src/ssr/cohort/allele_freq_prior.rs`).
//!
//! **How much mass there is** is new here, and it is new because production records a defect
//! against itself. `G₀` was fitted as a smoothing strength for a plug-in frequency estimate,
//! where its total means nothing else. Reused unchanged as a Dirichlet concentration that
//! total stops being a smoothing knob and becomes a claim about how polymorphic the locus is —
//! and how hard the prior will resist the reads. So ng sets the total from the one measured
//! quantity that answers that question, the cohort's repeat gene diversity, and leaves the
//! decay to do only the job it was fitted for (spec §5.1).
//!
//! ## The units error this module exists to avoid
//!
//! Gene diversity is a **probability** — the chance two copies drawn at random carry different
//! lengths. A concentration is a **count of chromosomes**. Setting the total to the diversity
//! equates the two, and the spec records that it was in the design document until 2026-08-19.
//!
//! What a Dirichlet with total `A` actually implies is `A(1 − c)/(A + 1)`, where `c` is the
//! **Simpson index** of the shape — the chance that two copies drawn from the shape alone land
//! on the same length. So a total of `D` asserts `D(1 − c)/(D + 1)`, which is always less.
//! **Measured on 1,236 polymorphic tomato repeat tracts** at the coded fallback decay, the
//! prior would assert a median 0.40 of what was measured, tenth percentile 0.22 (spec §5.1).
//! The total that does reproduce it is `D / (1 − c − D)`, which is what [`fill_ssr_seed`] uses,
//! and the identity is what `tests::the_seed_implies_the_diversity_that_was_measured` pins.

use crate::genetics::{ALPHA_REF, MIN_ALT_CONCENTRATION};
use crate::ng::types::{RepeatGeneDiversity, SeedDecayPerRepeat};

use super::Concentration;

/// The floor on a shape weight, so a candidate far from the mode keeps a strictly positive
/// share rather than falling into an absorbing zero.
///
/// **Ported from production's `G0_FLOOR`** (`src/ssr/cohort/allele_freq_prior.rs`), which gives
/// one reason for it: `decay^|Δ|` must not underflow to exactly `0.0` for a far candidate over
/// a long tract. At a decay of `0.5` that underflow is complete at an offset of **1,075 repeat
/// units and beyond** (`tests::a_candidate_too_far_for_the_decay_to_reach_keeps_its_share`).
/// Why a zero would matter is the spec's sentence rather than production's: a masked long
/// heterozygous copy that the candidate set nearly missed has to stay recoverable rather than
/// fall into a prior it can never climb out of (spec §5).
///
/// **Here it does a second job production's does not need: it keeps the normalisation
/// possible.** ng divides by the weights' total, and a tract whose every candidate sits past
/// the underflow distance would divide by zero and fill the buffer with `NaN`
/// (`tests::a_tract_whose_every_candidate_is_past_the_underflow_distance_is_still_seeded`).
///
/// It has production's value, which happens to equal [`MIN_ALT_CONCENTRATION`], and **the two
/// are not the same quantity**: this floors a dimensionless share of the shape, that floors a
/// count of chromosomes. Production picked this one as "any representable positive number";
/// `MIN_ALT_CONCENTRATION` is sized rather than arbitrary — `src/genetics.rs` puts it far below
/// any real diversity so it cannot perturb a genuine estimate.
const SHAPE_FLOOR: f64 = 1e-12;

/// What [`fill_ssr_seed`] did, and the buffer it wrote, on whichever of the two terms applies.
///
/// **A refusal hands back no concentration**, so the reads cannot be met with a starting point
/// that silently missed its target — the failure spec §12 test 11 exists to prevent. What it
/// hands back instead is the buffer itself, still holding the shape, because completing a
/// refusal is a policy decision that belongs to the caller and not here (spec Q2).
///
/// **This does not make the mistake unrepresentable, and it is not meant to.**
/// [`Concentration::new`] is public, so a caller can still wrap the shape by hand — and it has
/// to be able to, because the provisional policy `doc/devel/ng/arch/calling_priors.md` §5 names
/// for these loci is exactly that: scale the shape to some total and carry the marker onto the
/// locus's output. What the type does is make that an explicit decision with a name on it,
/// rather than something a caller falls into by reading a buffer it was never handed.
#[must_use]
#[derive(Debug)]
pub enum SsrSeedOutcome<'a> {
    /// The total that reproduces the measured diversity exists, and the buffer holds the shape
    /// scaled to it.
    Seeded(Concentration<'a>),
    /// **No total reproduces the measurement**, because `A(1 − c)/(A + 1)` rises to `1 − c` and
    /// stops: the shape itself cannot express that much diversity over this locus's candidate
    /// lengths, and rescaling is not a repair.
    ///
    /// **How often this fires is a property of the panel, not of the caller, and the two ends
    /// of the range are far apart.** On the 63-accession tomato panel at the coded fallback
    /// decay it is 119 of 1,236 polymorphic tracts, about one in ten — 242 at a decay of `0.3`
    /// and 49 at `0.7`. At **one outbred sample it is every locus**: a single diploid genome
    /// shows at most three lengths at a tract, whose shape can imply at most 0.444 (two
    /// lengths) or 0.625 (three) at the fallback decay, while the pre-pass reports that
    /// genome's own repeat diversity — about 0.72 on the GIAB HG002 benchmark, where 72 tandem
    /// repeats in 100 are heterozygous (spec §5.3). Whichever policy spec Q2 settles on has to
    /// work at ten refusals in ten, not only at one in ten
    /// (`tests::a_single_outbred_genome_is_refused_at_every_tract`).
    ///
    /// It also carries a guard against a total too large to represent. **No shape reaches it**:
    /// the worst a shape can ask for is its ceiling divided by that ceiling's own last unit in
    /// the last place, which for a binary float is at most `2^53` — measured at `2^53 − 1`,
    /// about 9.0e15 chromosomes (`tests::the_largest_total_any_shape_can_ask_for_is_finite`).
    /// The guard stays because the alternative is an infinity reaching a `debug_assert!` that
    /// release compiles out.
    DiversityUnreachable {
        /// The cohort's repeat gene diversity, as the pre-pass measured it.
        measured: f64,
        /// `1 − c`: the most this locus's candidate lengths can imply, at any total.
        ceiling: f64,
        /// The buffer, holding the shape normalised to sum to one — no total asserted. It is
        /// handed back mutably so a caller that has a policy for these loci can scale it in
        /// place and wrap the result, with nothing allocated and nothing recomputed.
        shape: &'a mut [f64],
    },
}

/// Fill `out` with the STR seed: mass falling off geometrically from `modal_repeat_count`,
/// scaled so that the prior's own implied gene diversity is the `gene_diversity` the pre-pass
/// measured.
///
/// ```text
/// shape:   w_j ∝ max(decay ^ |repeat count of j − modal repeat count|, SHAPE_FLOOR),  Σ w = 1
/// total:   A   = D / (1 − c − D),      c = Σ_j w_j²        (the shape's Simpson index)
/// seed:    α_j = A · w_j,  floored at MIN_ALT_CONCENTRATION
/// ```
///
/// **Why that total and not `D` itself:** the two are different quantities in different units,
/// and the substitution makes the prior assert a median 0.40 of what was measured on tomato —
/// see this module's own documentation, and spec §5.1.
///
/// `candidate_repeat_counts` is parallel to the locus's
/// [`CandidateAlleles`](crate::ng::calling::CandidateAlleles), so entry 0 is the reference
/// allele's — but **nothing here privileges entry 0**, which is the whole difference from the
/// SNP/indel path: at a repeat tract the reference length is one common length among several
/// (spec §5).
///
/// ## Three things a caller can get wrong that nothing here can catch
///
/// The lengths of the two slices are checkable and are checked. These are not, and each was
/// measured on one diploid sample at three candidate lengths with an inbreeding coefficient of
/// 0.3, as a difference within the prior row — which is all a prior means:
///
/// - **the reference allele's repeat count passed as the modal count.** The prior flips which
///   homozygote it favours, a swing of 1.44 nats (6.3 phred).
/// - **the SNP/indel diversity passed where the repeat diversity belongs.** The heterozygote is
///   penalised 7.49 nats instead of 2.97 — 4.52 nats (19.6 phred) more hostile to
///   heterozygotes. [`RepeatGeneDiversity`] exists so this needs a deliberate conversion.
/// - **the buffer reused across loci with the call skipped.** The previous locus's row, entry
///   for entry.
///
/// ## Two alleles of the same length land on one rung
///
/// The seed is keyed by repeat count, so two candidates that differ by an interior
/// substitution — an interrupted repeat — sit at the same offset and **each receives the
/// rung's full weight**, which is production's behaviour. That costs more than the extra prior
/// mass spec §5.2 names: a second spelling flattens the shape, which lowers its Simpson index
/// and **raises the ceiling from 0.444 to 0.625** on a two-length tract at the fallback decay
/// (`tests::a_second_spelling_of_one_length_raises_the_ceiling`), so whether a locus is refused
/// at all depends on how many spellings the cohort happened to show. Whether to divide the rung
/// instead needs the interrupted-repeat work to say how it should be weighted; the signature
/// takes counts rather than sequences precisely so that change lands in this one function
/// (spec §5.2, open as spec Q3).
///
/// ## The ends of the range
///
/// **At one candidate allele** the shape is a single entry, its Simpson index is exactly 1, and
/// the ceiling is 0 — so the rule above would refuse every monomorphic tract, whatever the
/// measurement. It is refused as a rule with nothing to refuse: a locus with one allele has one
/// genotype, whose prior probability is 1 at **any** positive concentration. This returns
/// [`SsrSeedOutcome::Seeded`] there, at one chromosome
/// (`tests::a_locus_with_one_candidate_length_is_seeded_rather_than_refused`).
///
/// **At one sample every tract is refused**, and that is the cohort end of the range rather
/// than an edge case — see [`SsrSeedOutcome::DiversityUnreachable`].
///
/// **As the measurement approaches the ceiling the total rises without bound**, because
/// `A = D/(1 − c − D)` has the ceiling as its pole: closing the gap to the ceiling by a factor
/// of ten multiplies the total by about ten
/// (`tests::the_total_climbs_without_bound_as_the_measurement_nears_the_ceiling`). That is the
/// honest answer and it is also an immovable prior — at a gap of one part in a million the seed
/// is worth about a million chromosomes, a hundred times the largest cohort this caller commits
/// to, so the reads cannot move it. Capping it is part of what spec Q2 has to settle, and this
/// function does not cap it.
///
/// ## Shape and cost
///
/// Fills the caller's buffer and **allocates nothing**. Four passes over the locus's candidate
/// lengths — write the raw weights, sum them, normalise and accumulate the Simpson index,
/// scale — and no `lgamma`, against one `lgamma` per allele a genotype carries a copy of in the
/// prior row this feeds.
///
/// **Both slice lengths are checked in release**: `out` is the caller's and is reused across
/// loci, so a short one would leave the previous locus's entries standing in this locus's seed,
/// which is the silent failure this module refuses everywhere.
pub fn fill_ssr_seed<'a>(
    candidate_repeat_counts: &[u32],
    modal_repeat_count: u32,
    decay: SeedDecayPerRepeat,
    gene_diversity: RepeatGeneDiversity,
    out: &'a mut [f64],
) -> SsrSeedOutcome<'a> {
    assert_lengths(candidate_repeat_counts, out);

    if candidate_repeat_counts.len() == 1 {
        // One length, one genotype, and its prior probability is 1 at any positive
        // concentration — so no total can be wrong here, and the ceiling rule below, which
        // would refuse every monomorphic tract, has nothing to refuse. The value is arbitrary
        // for that reason; it is `ALPHA_REF` so that a reader who finds it in a buffer
        // recognises it as one chromosome's worth rather than as a number this file chose.
        out[0] = ALPHA_REF;
        return SsrSeedOutcome::Seeded(Concentration::new(out));
    }

    let simpson_index = fill_seed_shape(candidate_repeat_counts, modal_repeat_count, decay, out);
    let ceiling = 1.0 - simpson_index;
    let measured = gene_diversity.get();

    let headroom = ceiling - measured;
    let total = measured / headroom;
    if headroom <= 0.0 || !total.is_finite() {
        // `out` keeps the normalised shape: the refusal asserts no total, and the shape is what
        // any policy for these loci will need (spec Q2).
        return SsrSeedOutcome::DiversityUnreachable {
            measured,
            ceiling,
            shape: out,
        };
    }

    for slot in out.iter_mut() {
        *slot = (total * *slot).max(MIN_ALT_CONCENTRATION);
    }
    SsrSeedOutcome::Seeded(Concentration::new(out))
}

/// Fill `out` with the seed's shape — the share of the prior's mass each candidate length
/// carries, summing to 1 — and return its **Simpson index**, `Σ w²`: the chance that two copies
/// drawn from the shape alone land on the same length.
///
/// `1 − c` is therefore the diversity the shape itself carries, and it is the ceiling
/// [`fill_ssr_seed`] measures the pre-pass's measurement against. It comes back from here rather
/// than being recomputed because the caller has just written every entry.
///
/// **The lengths are checked here too**, and not only in [`fill_ssr_seed`]: a longer `out`
/// leaves stale entries in the total it normalises against, and a shorter one drops candidates
/// silently. Two integer comparisons per locus, against a wrong shape at every locus after it.
fn fill_seed_shape(
    candidate_repeat_counts: &[u32],
    modal_repeat_count: u32,
    decay: SeedDecayPerRepeat,
    out: &mut [f64],
) -> f64 {
    assert_lengths(candidate_repeat_counts, out);
    let decay = decay.get();
    for (slot, &repeat_count) in out.iter_mut().zip(candidate_repeat_counts) {
        let offset = (i64::from(repeat_count) - i64::from(modal_repeat_count)).abs();
        // A tract cannot sit `i32::MAX` repeat units from the mode on any genome, and the
        // weight has underflowed to zero long before that anyway — the clamp is here so the
        // conversion cannot wrap rather than because the case is reachable. Wrapping would send
        // the offset negative and put nearly all the prior's mass on the *far* candidate
        // (`tests::an_offset_past_the_integer_range_is_clamped_rather_than_wrapped`).
        let offset = i32::try_from(offset).unwrap_or(i32::MAX);
        *slot = decay.powi(offset).max(SHAPE_FLOOR);
    }

    let raw_total: f64 = out.iter().sum();
    let mut simpson_index = 0.0;
    for slot in out.iter_mut() {
        *slot /= raw_total;
        simpson_index += *slot * *slot;
    }
    simpson_index
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

    /// What a written seed actually implies about diversity, computed from the buffer rather
    /// than from the arithmetic that produced it: `A(1 − c)/(A + 1)`, with both `A` and `c`
    /// read back off the entries.
    ///
    /// **Reading it back is the point.** Recomputing it from the total the builder chose would
    /// re-run the builder's own algebra and could only ever agree with itself; this sees what
    /// the caller of `Concentration::get` sees, floors included.
    fn implied_gene_diversity(alpha: &[f64]) -> f64 {
        let total: f64 = alpha.iter().sum();
        let simpson_index: f64 = alpha.iter().map(|a| (a / total) * (a / total)).sum();
        total * (1.0 - simpson_index) / (total + 1.0)
    }

    /// Seed a tract and hand back the entries, or panic naming the refusal.
    fn seeded(counts: &[u32], modal: u32, decay: f64, measured: f64, out: &mut [f64]) -> Vec<f64> {
        match fill_ssr_seed(
            counts,
            modal,
            SeedDecayPerRepeat::try_new(decay).unwrap(),
            RepeatGeneDiversity::try_new(measured).unwrap(),
            out,
        ) {
            SsrSeedOutcome::Seeded(concentration) => concentration.get().to_vec(),
            SsrSeedOutcome::DiversityUnreachable {
                measured, ceiling, ..
            } => panic!("expected a seed, got a refusal: measured {measured}, ceiling {ceiling}"),
        }
    }

    /// Seed a tract and hand back the refusal's two numbers, or panic naming the seed.
    fn refused(
        counts: &[u32],
        modal: u32,
        decay: f64,
        measured: f64,
        out: &mut [f64],
    ) -> (f64, f64) {
        match fill_ssr_seed(
            counts,
            modal,
            SeedDecayPerRepeat::try_new(decay).unwrap(),
            RepeatGeneDiversity::try_new(measured).unwrap(),
            out,
        ) {
            SsrSeedOutcome::DiversityUnreachable {
                measured, ceiling, ..
            } => (measured, ceiling),
            SsrSeedOutcome::Seeded(c) => {
                panic!("expected a refusal, got a seed: {:?}", c.get())
            }
        }
    }

    /// The shape's ceiling, `1 − c`, for a candidate set at a decay.
    fn ceiling_of(counts: &[u32], modal: u32, decay: f64) -> f64 {
        let mut shape = vec![f64::NAN; counts.len()];
        1.0 - fill_seed_shape(
            counts,
            modal,
            SeedDecayPerRepeat::try_new(decay).unwrap(),
            &mut shape,
        )
    }

    /// A locus's candidate lengths as offsets from the mode, which is what every fixture here
    /// actually cares about: `spread(&[0, 1, -2])` is a tract whose candidates sit at the mode,
    /// one repeat above it, and two below.
    fn spread(offsets: &[i32]) -> (Vec<u32>, u32) {
        const MODE: u32 = 20_000;
        let counts = offsets
            .iter()
            .map(|&d| u32::try_from(i64::from(MODE) + i64::from(d)).expect("fixture offset"))
            .collect();
        (counts, MODE)
    }

    /// **Spec §12 test 10, and the reason this module exists.** The seed's own implied gene
    /// diversity is the diversity that was measured — not its concentration total, which is a
    /// different quantity in different units and would make the prior assert about two-fifths
    /// of the measurement.
    ///
    /// Swept over decays from a steep `0.2` to a flat `1.0`, over candidate spreads from two
    /// lengths to eleven, and over diversities from a hundredth to a half — every combination
    /// whose measurement the shape can hold. The tolerance is on the **relative** error,
    /// because the diversities span more than an order of magnitude.
    ///
    /// **No entry in this sweep meets the concentration floor**, so what it measures is the
    /// identity's own arithmetic; the floored path is a different measurement and
    /// `tests::the_concentration_floor_lifts_the_implied_diversity_and_by_how_much` makes it.
    #[test]
    fn the_seed_implies_the_diversity_that_was_measured() {
        let spreads: [&[i32]; 5] = [
            &[0, 1],
            &[0, -1, 1],
            &[0, -1, 1, -2, 2],
            &[-3, -1, 0, 2, 5],
            &[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5],
        ];
        let mut worst_relative_error = 0.0_f64;
        let mut cases = 0_u32;
        let mut refusals = 0_u32;
        let mut floored_entries = 0_u32;

        for spread_offsets in spreads {
            let (counts, mode) = spread(spread_offsets);
            for decay in [0.2, 0.5, 0.7, 0.9, 1.0] {
                for measured in [0.01, 0.05, 0.087, 0.2, 0.5] {
                    let mut out = vec![f64::NAN; counts.len()];
                    let outcome = fill_ssr_seed(
                        &counts,
                        mode,
                        SeedDecayPerRepeat::try_new(decay).unwrap(),
                        RepeatGeneDiversity::try_new(measured).unwrap(),
                        &mut out,
                    );
                    let SsrSeedOutcome::Seeded(concentration) = outcome else {
                        refusals += 1; // the shape cannot hold it; the next test's case
                        continue;
                    };
                    cases += 1;
                    floored_entries += concentration
                        .get()
                        .iter()
                        .filter(|&&a| a == MIN_ALT_CONCENTRATION)
                        .count() as u32;
                    let implied = implied_gene_diversity(concentration.get());
                    worst_relative_error =
                        worst_relative_error.max(((implied - measured) / measured).abs());
                }
            }
        }

        // Pinned rather than bounded: a tolerance that held over three cases would look exactly
        // like one that held over 118, and a sweep that silently shrank is the way this test
        // stops testing anything. 118 + 7 is the whole 5 x 5 x 5 grid, so nothing is skipped.
        assert_eq!((cases, refusals), (118, 7));
        // All seven refusals sit at the sweep's widest diversity, 0.5. Five of them are the
        // two-length tract at every decay in the sweep — two lengths cannot imply more than 0.5
        // however flat the shape is — and the other two are the steepest decay, 0.2, on the
        // three-length and the ragged five-length tract.
        assert_eq!(
            floored_entries, 0,
            "the sweep must stay on the floor-free path"
        );
        // Measured, not chosen: the worst is 1.12e-15 relative, about seven units in the last
        // place of a diversity of 0.087, which is the arithmetic of the identity and not a bias.
        assert!(
            worst_relative_error < 2e-15,
            "worst relative error {worst_relative_error:e} over {cases} cases"
        );
    }

    /// Where the concentration floor lifts an entry, the identity drifts **upward** — the seed
    /// implies slightly more diversity than was asked for — and this measures how far.
    ///
    /// It is one-directional and bounded by about `n · 1e-12 / D`, so it becomes material only
    /// far below any diversity a cohort can resolve: at 50 candidate lengths and a diversity of
    /// 1 in 10,000 it is 3 parts in 10 million, against the 1.12e-15 of the floor-free path.
    #[test]
    fn the_concentration_floor_lifts_the_implied_diversity_and_by_how_much() {
        let far = 1_075; // past the underflow distance at the fallback decay
        let (counts, mode) = spread(&[0, 1, far]);
        let mut out = vec![f64::NAN; counts.len()];
        let alpha = seeded(&counts, mode, 0.5, 0.087, &mut out);
        assert_eq!(alpha[2], MIN_ALT_CONCENTRATION);
        let drift = (implied_gene_diversity(&alpha) - 0.087) / 0.087;
        assert!((0.0..2e-11).contains(&drift), "drift {drift:e}");

        let wide: Vec<i32> = (0..50).collect();
        let (counts, mode) = spread(&wide);
        let mut out = vec![f64::NAN; counts.len()];
        let alpha = seeded(&counts, mode, 0.5, 1e-4, &mut out);
        let floored = alpha
            .iter()
            .filter(|&&a| a == MIN_ALT_CONCENTRATION)
            .count();
        assert_eq!(floored, 23);
        let drift = (implied_gene_diversity(&alpha) - 1e-4) / 1e-4;
        assert!((0.0..1e-6).contains(&drift), "drift {drift:e}");
    }

    /// **Spec §12 test 11.** The refusal fires exactly at the bound and not a step early or
    /// late: a measurement a hair below the shape's ceiling is seeded, one at the ceiling and
    /// one above it are refused, and the refusal reports the ceiling it measured against.
    ///
    /// The ceiling is `1 − c`, so it is a property of the candidate spread and the decay alone
    /// — this walks four spreads at three decays and finds each one's bound from the shape
    /// itself rather than from a hard-coded number.
    #[test]
    fn a_diversity_the_shape_cannot_hold_is_refused_exactly_at_the_bound() {
        let spreads: [&[i32]; 4] = [&[0, 1], &[0, -1, 1], &[0, -1, 1, -2, 2], &[-3, 0, 2, 5]];
        for spread_offsets in spreads {
            let (counts, mode) = spread(spread_offsets);
            for decay in [0.3, 0.5, 0.7] {
                let ceiling = ceiling_of(&counts, mode, decay);
                let mut out = vec![f64::NAN; counts.len()];
                seeded(&counts, mode, decay, ceiling * (1.0 - 1e-9), &mut out);

                for measured in [ceiling, ceiling * (1.0 + 1e-9), (ceiling + 1.0) / 2.0, 1.0] {
                    assert_eq!(
                        refused(&counts, mode, decay, measured, &mut out),
                        (measured, ceiling),
                        "at the ceiling {ceiling}"
                    );
                }
            }
        }
    }

    /// A refusal hands back the **shape** — summing to one, with no total asserted — so
    /// whichever policy spec Q2 settles on has the geometry to work from, and it hands it back
    /// through the outcome rather than leaving the caller to read a buffer it was not given.
    #[test]
    fn a_refusal_hands_back_the_shape_and_asserts_no_total() {
        let (counts, mode) = spread(&[0, -1, 1, -2, 2]);
        let ceiling = ceiling_of(&counts, mode, 0.5);
        let mut expected = vec![f64::NAN; counts.len()];
        fill_seed_shape(&counts, mode, SeedDecayPerRepeat::FALLBACK, &mut expected);

        let mut out = vec![f64::NAN; counts.len()];
        let outcome = fill_ssr_seed(
            &counts,
            mode,
            SeedDecayPerRepeat::FALLBACK,
            RepeatGeneDiversity::try_new((ceiling + 1.0) / 2.0).unwrap(),
            &mut out,
        );
        let SsrSeedOutcome::DiversityUnreachable { shape, .. } = outcome else {
            panic!("expected a refusal");
        };
        assert_eq!(shape, expected.as_slice());
        assert!((shape.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    }

    /// **One outbred genome is refused at every tract**, which is the cohort end of this
    /// caller's committed range and not an edge case.
    ///
    /// A single diploid sample shows at most three lengths at a tract, and the pre-pass reports
    /// that genome's own repeat diversity — about 0.72 on GIAB HG002, where 72 tandem repeats
    /// in 100 are heterozygous. Both candidate sets a single sample can produce fall well short
    /// of it, and no decay rescues them: the ceiling saturates at 0.8148 at the fallback decay
    /// however many lengths a tract carries, so this is a property of the shape rather than of
    /// the candidate set being thin.
    #[test]
    fn a_single_outbred_genome_is_refused_at_every_tract() {
        const HG002_REPEAT_DIVERSITY: f64 = 0.72;
        let (two_lengths, mode) = spread(&[0, 1]);
        let (three_lengths, _) = spread(&[0, -1, 1]);

        assert!((ceiling_of(&two_lengths, mode, 0.5) - 0.4444444444444444).abs() < 1e-15);
        assert!((ceiling_of(&three_lengths, mode, 0.5) - 0.625).abs() < 1e-15);

        for counts in [&two_lengths, &three_lengths] {
            let mut out = vec![f64::NAN; counts.len()];
            let (measured, _) = refused(counts, mode, 0.5, HG002_REPEAT_DIVERSITY, &mut out);
            assert_eq!(measured, HG002_REPEAT_DIVERSITY);
        }

        // And no candidate set reaches it: the ceiling saturates at 1 − 5/27 = 0.8148 as the
        // shape spreads over every offset, so a panel reporting more than that is refused
        // everywhere whatever its tracts look like.
        let wide: Vec<i32> = (-60..=60).collect();
        let (counts, mode) = spread(&wide);
        let saturated = ceiling_of(&counts, mode, 0.5);
        // 1 − 5/27 = 0.814814…, to within the 5 parts in a trillion the shape floor adds by
        // lifting the candidates past offset 40, where `0.5^|Δ|` falls under it.
        assert!(
            (saturated - (1.0 - 5.0 / 27.0)).abs() < 1e-10,
            "{saturated}"
        );
        assert!(saturated < 0.815);
    }

    /// The shape is production's `G₀`, ported: consecutive repeat units away from the mode are
    /// a constant factor apart, and it is the same factor on both sides.
    ///
    /// This is a **ratio** test rather than a value test, because normalising the shape divides
    /// every entry by the same total and the ratios are what survives it.
    #[test]
    fn mass_falls_off_by_the_decay_for_every_repeat_unit_away_from_the_mode() {
        let (counts, mode) = spread(&[0, 1, 2, 3, -1, -2, -3]);
        for decay in [0.2, 0.5, 0.9] {
            let mut shape = vec![f64::NAN; counts.len()];
            fill_seed_shape(
                &counts,
                mode,
                SeedDecayPerRepeat::try_new(decay).unwrap(),
                &mut shape,
            );

            // Entries 0..4 are the mode and three lengths above it; 0, 4, 5, 6 are the mode and
            // three below.
            for pair in [[0, 1], [1, 2], [2, 3]] {
                assert!((shape[pair[1]] / shape[pair[0]] - decay).abs() < 1e-15);
            }
            for pair in [[0, 4], [4, 5], [5, 6]] {
                assert!((shape[pair[1]] / shape[pair[0]] - decay).abs() < 1e-15);
            }
        }
    }

    /// **The shape floor, at the one place it is the thing doing the work.** A candidate far
    /// enough from the mode that `decay^|Δ|` underflows keeps a strictly positive *share*, and
    /// this measures that share rather than the concentration — because on the seeded path
    /// `MIN_ALT_CONCENTRATION` would lift a zero to `1e-12` and the test would pass with the
    /// shape floor deleted.
    ///
    /// At the fallback decay the underflow is complete at an offset of **1,075 repeat units**,
    /// which this finds rather than assumes.
    #[test]
    fn a_candidate_too_far_for_the_decay_to_reach_keeps_its_share() {
        let decay = SeedDecayPerRepeat::FALLBACK;
        let underflow_offset = (1..)
            .find(|&d| decay.get().powi(d) == 0.0)
            .expect("underflow");
        assert_eq!(underflow_offset, 1075);

        let (counts, mode) = spread(&[0, 1, underflow_offset]);
        let mut shape = vec![f64::NAN; counts.len()];
        fill_seed_shape(&counts, mode, decay, &mut shape);
        // The raw weights are 1, 0.5 and the floor, so the far candidate's share is the floor
        // over their total. Exactly, not merely positive: `> 0.0` passes with the floor gone.
        let expected_far_share = SHAPE_FLOOR / (1.0 + 0.5 + SHAPE_FLOOR);
        assert!(
            (shape[2] / expected_far_share - 1.0).abs() < 1e-12,
            "{:e}",
            shape[2]
        );

        // And at a diversity near the ceiling the total is large enough that the share survives
        // into the concentration, above `MIN_ALT_CONCENTRATION` rather than rescued by it.
        let ceiling = ceiling_of(&counts, mode, decay.get());
        let mut out = vec![f64::NAN; counts.len()];
        let alpha = seeded(&counts, mode, decay.get(), ceiling * (1.0 - 1e-9), &mut out);
        assert!(alpha[2] > MIN_ALT_CONCENTRATION, "{:e}", alpha[2]);
    }

    /// A tract whose **every** candidate sits past the underflow distance is still seeded. This
    /// is the shape floor's second job: without it every raw weight is exactly zero, the total
    /// they are normalised against is zero, and the buffer fills with `NaN`.
    #[test]
    fn a_tract_whose_every_candidate_is_past_the_underflow_distance_is_still_seeded() {
        let (counts, mode) = spread(&[2_000, 3_000]);
        let mut out = vec![f64::NAN; counts.len()];
        let alpha = seeded(&counts, mode, 0.5, 0.087, &mut out);
        // Both candidates underflowed to the same floor, so the shape is flat and they share
        // the mass equally.
        assert!(alpha.iter().all(|a| a.is_finite() && *a > 0.0), "{alpha:?}");
        assert_eq!(alpha[0], alpha[1]);
    }

    /// An offset past the integer range is clamped rather than wrapped. A wrapping conversion
    /// sends the offset negative, so `decay.powi` returns a huge weight and nearly all the
    /// prior's mass lands on the **far** candidate instead of the mode — the shape inverted.
    #[test]
    fn an_offset_past_the_integer_range_is_clamped_rather_than_wrapped() {
        let mut shape = [f64::NAN; 2];
        fill_seed_shape(
            &[20, u32::MAX],
            20,
            SeedDecayPerRepeat::FALLBACK,
            &mut shape,
        );
        assert!(shape[0] > 0.999, "the mode must keep the mass: {shape:?}");
        assert!(shape[1] < 1e-11, "the far candidate must not: {shape:?}");
    }

    /// **One candidate length is seeded, not refused.** Its shape has a Simpson index of
    /// exactly 1 and therefore a ceiling of 0, so the rule that governs every other locus would
    /// refuse it whatever the measurement — including a measurement of zero. There is nothing
    /// to refuse: one length means one genotype, whose prior probability is 1 at any positive
    /// concentration.
    #[test]
    fn a_locus_with_one_candidate_length_is_seeded_rather_than_refused() {
        for measured in [0.0, 0.087, 0.9, 1.0] {
            let mut out = [f64::NAN];
            assert_eq!(seeded(&[20], 20, 0.5, measured, &mut out), [ALPHA_REF]);
        }

        // And the ceiling really is zero there, which is what makes the branch necessary rather
        // than defensive.
        let mut shape = [f64::NAN];
        let simpson_index = fill_seed_shape(&[20], 20, SeedDecayPerRepeat::FALLBACK, &mut shape);
        assert_eq!(simpson_index, 1.0);
    }

    /// A cohort with no repeat variation at all seeds every candidate at the floor — the
    /// concentration total goes to zero with the measurement, and what survives is an
    /// effectively-certain homozygote of some length, with the reads left to say which.
    ///
    /// That is the repeat-tract twin of production's `θ = 0` behaviour on the SNP path, where
    /// `MIN_ALT_CONCENTRATION` yields an effectively-certain hom-ref prior
    /// (`src/genetics.rs`). **The shape's own preference for the mode is what is lost**, and
    /// this measures how far up the range that loss reaches: the mode's entry clears the floor
    /// again at a diversity of 1.8 in a trillion, so every diversity a pre-pass could report
    /// keeps the geometry.
    #[test]
    fn a_cohort_with_no_repeat_variation_seeds_every_candidate_at_the_floor() {
        let (counts, mode) = spread(&[0, -1, 1, -2, 2]);
        let mut out = vec![f64::NAN; counts.len()];
        let alpha = seeded(&counts, mode, 0.5, 0.0, &mut out);
        assert!(
            alpha.iter().all(|&a| a == MIN_ALT_CONCENTRATION),
            "{alpha:?}"
        );

        // Where the geometry comes back: the smallest diversity at which the modal candidate
        // sits above the floor, found by bisection so the figure is measured, not recalled.
        let mode_is_above_the_floor = |measured: f64| {
            let mut out = vec![f64::NAN; counts.len()];
            seeded(&counts, mode, 0.5, measured, &mut out)[0] > MIN_ALT_CONCENTRATION
        };
        let (mut low, mut high) = (0.0_f64, 1e-6_f64);
        assert!(mode_is_above_the_floor(high));
        for _ in 0..200 {
            let mid = 0.5 * (low + high);
            if mode_is_above_the_floor(mid) {
                high = mid;
            } else {
                low = mid;
            }
        }
        assert!(
            (1.84e-12..1.86e-12).contains(&high),
            "the geometry returns at a diversity of {high:e}"
        );
    }

    /// As the measurement climbs toward the shape's ceiling the total rises without bound,
    /// because the ceiling is the pole of `A = D/(1 − c − D)`. **That is the honest answer and
    /// it is also an immovable prior**, which is half of what spec Q2 has to settle.
    ///
    /// The sizes, on a five-length tract at the fallback decay: closing the gap to the ceiling
    /// by a factor of ten multiplies the total by about ten, and at a gap of one part in a
    /// million the seed is worth about a million chromosomes — a hundred times the largest
    /// cohort this caller commits to, so no cohort can move it.
    #[test]
    fn the_total_climbs_without_bound_as_the_measurement_nears_the_ceiling() {
        let (counts, mode) = spread(&[0, -1, 1, -2, 2]);
        let ceiling = ceiling_of(&counts, mode, 0.5);

        let mut totals = Vec::new();
        for gap_fraction in [1e-2, 1e-3, 1e-4, 1e-6] {
            let mut out = vec![f64::NAN; counts.len()];
            let alpha = seeded(&counts, mode, 0.5, ceiling * (1.0 - gap_fraction), &mut out);
            totals.push(alpha.iter().sum::<f64>());
        }

        for pair in totals.windows(2).take(2) {
            let growth = pair[1] / pair[0];
            assert!((9.0..11.0).contains(&growth), "total grew by {growth}");
        }
        // Several thousand samples is the top of this caller's committed cohort range, so
        // roughly ten thousand chromosomes; the last seed outweighs all of them a hundred-fold.
        assert!(
            (999_000.0..1_001_000.0).contains(&totals[3]),
            "{}",
            totals[3]
        );
    }

    /// **The largest total any shape can ask for is finite**, so the overflow arm of the
    /// refusal is a guard rather than a case — and this is the measurement that says so.
    ///
    /// The worst input is the largest measurement strictly below a shape's ceiling — one unit
    /// in the last place under it — and the total it asks for is then the ceiling divided by
    /// its own last unit in the last place. For a binary float that ratio is at most `2^53`,
    /// which is where the measured worst lands: **9.007e15 chromosomes**, exactly `2^53 − 1`,
    /// over four spreads at five decays. So the arm is unreachable by construction rather than
    /// merely unreached by this sweep.
    ///
    /// It matters that the guard stays: `out` is written before anything reads it back, and the
    /// finiteness check on [`Concentration`] is a `debug_assert!` that release compiles out, so
    /// an infinity here would travel into a prior row as a `NaN` with nothing raised.
    #[test]
    fn the_largest_total_any_shape_can_ask_for_is_finite() {
        let spreads: [&[i32]; 4] = [&[0, 1], &[0, -1, 1], &[0, -1, 1, -2, 2], &[-3, 0, 2, 5]];
        let mut worst_total = 0.0_f64;
        for spread_offsets in spreads {
            let (counts, mode) = spread(spread_offsets);
            for decay in [0.2, 0.5, 0.7, 0.9, 1.0] {
                let ceiling = ceiling_of(&counts, mode, decay);
                let measured = f64::from_bits(ceiling.to_bits() - 1);
                assert!(measured < ceiling);
                let mut out = vec![f64::NAN; counts.len()];
                let total: f64 = seeded(&counts, mode, decay, measured, &mut out)
                    .iter()
                    .sum();
                assert!(total.is_finite(), "{total} at decay {decay}");
                worst_total = worst_total.max(total);
            }
        }
        assert_eq!(worst_total, (1_u64 << 53) as f64 - 1.0, "worst total");
    }

    /// Two candidates of the same length land on one rung and **each takes its full weight** —
    /// production's behaviour, kept deliberately until the interrupted-repeat work says how to
    /// divide it (spec §5.2, Q3).
    ///
    /// Recorded as a test rather than as a comment so that the day the division lands, the
    /// change is visible as this test failing rather than as a silent shift in every
    /// interrupted tract's prior.
    #[test]
    fn two_candidates_of_one_length_each_take_the_rungs_full_weight() {
        let (one_spelling, mode) = spread(&[0, 1]);
        let (two_spellings, _) = spread(&[0, 1, 1]);
        let decay = SeedDecayPerRepeat::FALLBACK;

        let mut single = vec![f64::NAN; one_spelling.len()];
        fill_seed_shape(&one_spelling, mode, decay, &mut single);
        let mut doubled = vec![f64::NAN; two_spellings.len()];
        fill_seed_shape(&two_spellings, mode, decay, &mut doubled);

        // The second spelling adds a whole extra rung to the raw total, 1.5 to 2, so the mode's
        // share falls from 1/1.5 to 1/2.
        assert!((single[0] - 1.0 / 1.5).abs() < 1e-15);
        assert!((doubled[0] - 0.5).abs() < 1e-15);
        assert_eq!(doubled[1], doubled[2]);
    }

    /// **A second spelling of one length also raises the ceiling**, which is spec Q3 with a
    /// size attached: it flattens the shape, which lowers the Simpson index, so whether a tract
    /// is refused at all depends on how many spellings of one length the cohort happened to
    /// show — not only how much prior mass it collects.
    #[test]
    fn a_second_spelling_of_one_length_raises_the_ceiling() {
        let (one_spelling, mode) = spread(&[0, 1]);
        let (two_spellings, _) = spread(&[0, 1, 1]);
        assert!((ceiling_of(&one_spelling, mode, 0.5) - 0.4444444444444444).abs() < 1e-15);
        assert!((ceiling_of(&two_spellings, mode, 0.5) - 0.625).abs() < 1e-15);

        // So a measurement between the two is refused at one tract and seeded at the other.
        let mut out = vec![f64::NAN; one_spelling.len()];
        refused(&one_spelling, mode, 0.5, 0.5, &mut out);
        let mut out = vec![f64::NAN; two_spellings.len()];
        seeded(&two_spellings, mode, 0.5, 0.5, &mut out);
    }

    /// A buffer that does not match the candidate set is refused **in release**, because `out`
    /// is the loop's and is reused across loci: a short one would leave the previous tract's
    /// entries standing in this one's seed.
    #[test]
    #[should_panic(expected = "cover the locus's candidate lengths exactly")]
    fn a_mis_sized_buffer_is_refused() {
        let mut out = [f64::NAN; 2];
        let _ = fill_ssr_seed(
            &[20, 21, 22],
            20,
            SeedDecayPerRepeat::FALLBACK,
            RepeatGeneDiversity::try_new(0.087).unwrap(),
            &mut out,
        );
    }

    /// The shape filler makes the same check, because plan step E2 makes it public and its
    /// caller's assertion then stops covering it.
    #[test]
    #[should_panic(expected = "cover the locus's candidate lengths exactly")]
    fn a_mis_sized_buffer_is_refused_by_the_shape_filler_too() {
        let mut out = [f64::NAN; 2];
        fill_seed_shape(&[20, 21, 22], 20, SeedDecayPerRepeat::FALLBACK, &mut out);
    }

    /// An empty candidate set is refused: every repeat tract has a reference allele, so it is a
    /// caller that has lost track of which locus it is on rather than a thin one.
    #[test]
    #[should_panic(expected = "at least one length")]
    fn an_empty_candidate_set_is_refused() {
        let mut out: [f64; 0] = [];
        let _ = fill_ssr_seed(
            &[],
            20,
            SeedDecayPerRepeat::FALLBACK,
            RepeatGeneDiversity::try_new(0.087).unwrap(),
            &mut out,
        );
    }

    /// The two scalars refuse what is not their own quantity, and they refuse it with a
    /// `DomainError` rather than a panic — a degenerate fit is a run to abandon with a message,
    /// not a process to abort.
    #[test]
    fn the_scalars_refuse_values_that_are_not_their_own_quantity() {
        for bad in [0.0, -0.5, 1.5, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    SeedDecayPerRepeat::try_new(bad),
                    Err(crate::ng::types::DomainError::SeedDecayPerRepeat(_))
                ),
                "decay {bad}"
            );
        }
        assert_eq!(SeedDecayPerRepeat::try_new(1.0).unwrap().get(), 1.0);
        assert_eq!(SeedDecayPerRepeat::FALLBACK.get(), 0.5);

        for bad in [-0.5, 1.5, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    RepeatGeneDiversity::try_new(bad),
                    Err(crate::ng::types::DomainError::RepeatGeneDiversity(_))
                ),
                "diversity {bad}"
            );
        }
        assert_eq!(RepeatGeneDiversity::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(RepeatGeneDiversity::try_new(1.0).unwrap().get(), 1.0);
    }
}
