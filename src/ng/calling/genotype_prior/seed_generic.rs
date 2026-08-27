//! The SNP/indel starting point: the run's two concentration numbers — the chromosomes'
//! worth of prior belief attached to the reference allele and to the alternatives — built from
//! the two things the pre-pass measured about the population.
//!
//! At an ordinary site most alternative alleles are rare, so the one chromosome the reference
//! records is almost always the common one and the reference's number is the larger — but only
//! just, and *how much* larger is measured rather than fixed
//! (`doc/devel/ng/spec/calling_priors.md` §4.1).
//!
//! **Two functions.** [`seed_from_population_moments`] takes the population's mean
//! alternative-allele frequency and its heterozygosity and returns the run's seed;
//! [`fill_locus_concentration`] spreads that one pair over each locus's own alleles.
//!
//! **⛔ What used to stand between them, and left on 2026-08-27.** The seed's expected frequency
//! was not measured but *searched for*: the joint fit's population curve was evaluated into the
//! `2N + 1` allele-count classes a panel of `N` diploid individuals has, and a two-parameter pair
//! was fitted to those classes. That objective has `2N + 1` terms, so its answer moved with the
//! cohort — measured against the panel's own genotypes it was the worst of three routes in **34 of
//! 36 cells**, and **0.749×** the truth at its worst
//! (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §9.1). It is replaced by two
//! integrals of the same curve, which reproduce its first two moments exactly at every panel size
//! (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §2).
//!
//! **The panel's size appears nowhere in this module.** That is the point of the change: the two
//! numbers are statements about a population.

use crate::genetics::MIN_ALT_CONCENTRATION;
use crate::ng::types::{ExpectedAlternativeFrequency, ExpectedHeterozygosity};

use super::{Concentration, SeedRegime, SpectrumSeed};

/// The reference allele's concentration on a neutral panel — where the fit lands rather than a
/// number anyone chose, and the `1/p` density written as a Dirichlet
/// (`doc/devel/ng/spec/calling_priors.md` §4). It is what holds the heterozygote to
/// homozygous-alternative prior ratio near 2:1, so raising it is the §2.3 trap.
const NEUTRAL_ALPHA_REF: f64 = 1.0;

/// Which kind of variant a run's seed is for.
///
/// **Both classes are handed the same diversity today, and the argument exists so that stops
/// being true without touching a call site** (`doc/devel/ng/spec/calling_priors.md` §4.2, Q1).
/// Production ran different pseudocounts for the two — `0.01` against `0.00125`, an 8:1 ratio
/// inherited from another tool and never measured here — while ng's pre-pass measures one
/// heterozygosity for both, because it sums a windowed histogram that does not separate
/// substitutions from short insertions and deletions. Splitting the prior before the estimate is
/// split would mean inventing the ratio.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum VariantClass {
    /// One base carried in place of another — a single-nucleotide polymorphism.
    Substitution,
    /// A short insertion or deletion.
    InsertionOrDeletion,
}

// **What stood here until 2026-08-27: the blend toward the neutral shape, and its one constant.**
//
// The seed's expected frequency used to be interpolated in log space between the neutral shape's
// `θ / (1 + θ)` and the panel's own fitted one, at a weight `N / (N + N₀)` that rose with the
// panel — three functions and a fitted `N₀` of a quarter of an individual.
//
// **It existed to damp the small-panel noise of a search that is itself being deleted**, and the
// sweep that was to set `N₀` said the blend points the wrong way: the panel's own fitted shape is
// exact at one individual and degrades as the panel grows. **All three arms of that sweep's
// headline table — the one averaged over panel size, depth and population — put the best
// half-weight panel size at zero**; its depth-crossed arms did not agree with each other, putting
// it at 0 on a strong rare-allele pile-up and at **200** on a moderate one, which is the
// two-hundred-fold disagreement that says no single constant is right
// (`doc/devel/reports/implementations/ng_seed_shrinkage_2026-08-26.md` §5.2).
//
// Handed the population exactly, the blended seed came back at **0.62× to 0.92× of the truth at
// one individual** across four populations
// (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §9,
// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §6.1).
//
// **The frequency is now integrated off the fitted curve, and there is nothing to blend toward.**
// The program that ran the sweep, `examples/ng_seed_shape_weight_sweep.rs`, was deleted with the
// search it measured through; its numbers are in
// `doc/devel/reports/implementations/ng_seed_shrinkage_2026-08-26.md`.

/// The largest diversity **any** concentration pair can imply.
///
/// A pair of expected frequency `f` implies `2 f (1 − f) · A / (A + 1)`, and `2 f (1 − f)` is at
/// most a half, at `f = 1/2`. So a measured heterozygosity above this is not a thin estimate of
/// something a seed could carry — it is a fit that did not converge
/// (`doc/devel/ng/spec/ordinary_site_seed.md` §3.1).
const MAX_IMPLIED_DIVERSITY: f64 = 0.5;

/// **Solve for how much conviction a pair of this expected frequency needs in order to imply the
/// diversity that was measured** — `doc/devel/ng/spec/ordinary_site_seed.md` §3's identity.
///
/// A Dirichlet with total `A` and expected frequency `f` makes a diploid drawn from it
/// heterozygous at `2 f (1 − f) · A / (A + 1)` — the Beta-binomial at one alternative copy in
/// two draws. Fixing that to the measured `θ` fixes `A`:
///
/// ```text
///   t = θ / (2 f (1 − f)),        A = t / (1 − t)
/// ```
///
/// `t` is the share of the shape's own ceiling the measurement asks for.
///
/// ## Why there is no "no total reaches it" answer here, and there used to be
///
/// **`2 f (1 − f)` is the largest heterozygosity a pair of expected frequency `f` can imply**, so
/// a measurement at or above it has no answer at all rather than a large one. Until 2026-08-27
/// this function said so, and the seed carried a `DiversityUnreachable` regime that fell back to
/// the neutral pair and reported which of two ways it had got there.
///
/// **On the route that ships, that state cannot arise, and the reason is Jensen's inequality.**
/// Both numbers are integrals of one fitted population curve:
///
/// ```text
///   θ  =  E[2 f (1 − f)]  =  2 E[f] − 2 E[f²]
///   the ceiling            =  2 E[f] (1 − E[f])  =  2 E[f] − 2 E[f]²
/// ```
///
/// and `E[f²] ≥ E[f]²` — that difference *is* the spread of the population's frequencies, taken
/// over the whole density, both point masses included. So the measurement is below the ceiling by
/// exactly twice that spread, and equals it only where the whole population sits at one frequency.
/// **No density the fit can produce does**, and the margin has a closed form rather than an
/// estimate. Where every position segregates and neither point mass carries anything, the density
/// is a bare `Beta(a, b)` and the algebra collapses:
///
/// ```text
///   1 − θ / (2 E[f] (1 − E[f]))  =  1 / (a + b + 1),      so      A  =  a + b
/// ```
///
/// The fit clamps `a` and `b` to `[0.02, 50]` independently, so the tightest that gets is
/// `a = b = 50` — one part in 101 of the ceiling left unused, and a solved total of 100
/// chromosomes. **Adding either point mass widens the margin**, and that is measured rather than
/// argued: swept over the whole box the fit clamps to, the largest share of the ceiling any
/// density asks for is that same 100 in 101, and it is asked for only where the point masses carry
/// nothing ([`seed_tests::no_density_the_fit_can_produce_comes_within_one_part_in_a_hundred_of_its_ceiling`],
/// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §6).
///
/// So this is now a **release assertion** rather than a fall-back: tripping it means the two
/// numbers did not come from one curve, which is a run-assembly defect and not a thin cohort.
///
/// # Panics
///
/// On a measured heterozygosity at or above `2 f (1 − f)`, per the paragraph above.
fn total_for_diversity(expected_frequency: f64, diversity: f64) -> f64 {
    debug_assert!(
        expected_frequency > 0.0 && expected_frequency < 1.0,
        "the expected alternative-allele frequency is strictly inside (0, 1); got \
         {expected_frequency}"
    );
    debug_assert!(
        diversity > 0.0 && diversity <= MAX_IMPLIED_DIVERSITY,
        "a diversity of zero and one above a half are both taken before this point; got \
         {diversity}"
    );
    let shapes_ceiling = 2.0 * expected_frequency * (1.0 - expected_frequency);
    let share_of_ceiling = diversity / shapes_ceiling;
    // **At or above the ceiling rather than above it**, and the difference is not pedantry: at
    // exactly the ceiling the solved total is infinite, which `SpectrumSeed` refuses — so a
    // measurement that lands there would panic three frames later with a message about a
    // concentration. The division can also *round* to one from a measurement just below the
    // ceiling, and that lands in the same place.
    //
    // Held in release, because what it catches is a silently wrong seed rather than a crash: a
    // negative solved total makes `SpectrumSeed::new` refuse, but a total that overflowed to
    // infinity or came back enormous would be a prior no depth of reads could move.
    assert!(
        share_of_ceiling < 1.0,
        "a heterozygosity of {diversity} cannot come from a population whose mean \
         alternative-allele frequency is {expected_frequency}: the most such a population can \
         show is 2 f (1 - f) = {shapes_ceiling}, and it reaches that only if every position sits \
         at exactly that frequency. Two moments of one fitted curve always satisfy this strictly \
         — E[2f(1-f)] is below 2 E[f] (1 - E[f]) by twice the spread of the population's \
         frequencies — so these two numbers did not come from one curve"
    );
    share_of_ceiling / (1.0 - share_of_ceiling)
}

/// **Build the run's two starting numbers from the two things the pre-pass measured about the
/// population**: how often a chromosome drawn at random carries something other than the
/// reference base, and how often two chromosomes drawn at random differ.
///
/// What comes out is the reference allele's concentration and the total shared out across
/// whatever alternative alleles a locus turns out to carry.
///
/// ## What it does, in one sentence each
///
/// **A concentration pair is exactly an expected frequency and a total conviction in other
/// clothes**, so the two measurements determine it. The frequency is used as it stands. The total
/// is whatever makes a pair of that frequency imply the heterozygosity that was measured
/// ([`total_for_diversity`]). Neither number is chosen, and neither is fitted here
/// (`doc/devel/ng/spec/ordinary_site_seed.md` §3,
/// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §2).
///
/// **The panel's size appears nowhere in this function, and that is the point of the change that
/// put it here.** Until 2026-08-27 the frequency was not measured but *searched for*: the fitted
/// population curve was evaluated into the `2N + 1` allele-count classes a panel of `N` diploid
/// individuals has, and a two-parameter pair was fitted to those classes. That objective has
/// `2N + 1` terms, so its answer moved with the cohort — measured against the panel's own
/// genotypes it was the worst of three routes in **34 of 36 cells**, and **0.749×** the truth at
/// its worst, on a population with two frequency peaks at 63 individuals
/// (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §9.1). The frequency is a
/// property of the population, and now nothing about the panel reaches it.
///
/// ## The regimes, and none of them is a branch on cohort size
///
/// - **Both moments arrived** — the frequency is used and the total is pinned to the
///   heterozygosity ([`SeedRegime::FittedCurve`]).
/// - **A heterozygosity but no frequency** — there is no measured shape, so the pair is the
///   neutral `(1, θ)` ([`SeedRegime::NeutralShape`]). This is the middle rung of
///   `population_diversity.md` §3.4's ladder, and the per-sample histogram route is what is
///   meant to reach it.
/// - **No heterozygosity at all** — the neutral pair at the species-range guess, and the run must
///   say so ([`SeedRegime::FallbackDiversity`]). **A frequency without a heterozygosity lands
///   here too**, and it is discarded: the total comes from the heterozygosity, so a run without
///   one has nothing to pin a shape to. In practice the two arrive together — the joint route
///   fits the population's density and reads both moments off it.
/// - **The measured heterozygosity is exactly zero** — a cohort with no variation at all. Every
///   entry of the solved pair goes to zero with it, so the alternative concentration is floored at
///   [`MIN_ALT_CONCENTRATION`] and the run says the diversity was zero
///   ([`SeedRegime::ZeroDiversity`]).
///
/// **There is no fifth regime for a heterozygosity no pair can reach**, and until 2026-08-27 there
/// was. Two moments of one fitted curve always leave the ceiling `2 f (1 − f)` above the
/// measurement, by Jensen's inequality — [`total_for_diversity`] carries the argument and holds it
/// as a release assertion.
///
/// **Nothing here tests how many samples the cohort holds.**
///
/// ## It takes no variant class, and that is a decision rather than an omission
///
/// The design keeps the door open to giving substitutions and short insertions or deletions
/// different diversities (`calling_priors.md` §4.2, Q1). **The end that will apply that split is
/// [`fill_locus_concentration`], not this one** — settled by the owner, 2026-08-22. This function
/// reads how variable the population is, which the pre-pass fits without separating the two
/// classes; a class-specific *scale* belongs where the run's total is shared out over a locus's
/// alleles, and that is the only end able to describe a site carrying one alternative of each
/// kind. Carrying the argument at both ends would apply the ratio twice.
///
/// # Panics
///
/// **On a measured heterozygosity above a half**, which no concentration pair can imply
/// ([`MAX_IMPLIED_DIVERSITY`]). `ExpectedHeterozygosity` admits the whole of `[0, 1]`, so this is
/// expressible; it is refused at the run's assembly rather than carried, because a heterozygosity
/// above a half is a fit that did not converge (`ordinary_site_seed.md` §3.1).
///
/// ## Cost
///
/// **Three multiplications and a divide.** There is no fit here and no spectrum: the two moments
/// arrive already measured, and the identity that turns them into a pair is closed-form. The
/// search this replaced cost 399 predictions of the panel's own spectrum, which was 3.8 seconds at
/// 400 individuals and **11.8 minutes at 3,200**.
#[must_use]
pub fn seed_from_population_moments(
    expected_frequency: Option<ExpectedAlternativeFrequency>,
    diversity: Option<ExpectedHeterozygosity>,
) -> SpectrumSeed {
    let Some(measured_diversity) = diversity.map(ExpectedHeterozygosity::get) else {
        return SpectrumSeed::new(
            NEUTRAL_ALPHA_REF,
            ExpectedHeterozygosity::SPECIES_FALLBACK.get(),
            SeedRegime::FallbackDiversity,
        );
    };
    assert!(
        measured_diversity <= MAX_IMPLIED_DIVERSITY,
        "the run's fitted heterozygosity is {measured_diversity}, and no pair of concentrations \
         can imply more than {MAX_IMPLIED_DIVERSITY}: a diploid drawn from a pair of expected \
         frequency f is heterozygous at most 2 f (1 - f) of the time. A heterozygosity above a \
         half is not a thin estimate — it is a population fit that did not converge, and the run \
         is refused here rather than seeded from it"
    );
    if measured_diversity == 0.0 {
        // Solving for the total gives zero and every entry with it, so the alternative
        // concentration is floored the way every other seed builder in this tree floors one.
        return SpectrumSeed::new(
            NEUTRAL_ALPHA_REF,
            MIN_ALT_CONCENTRATION,
            SeedRegime::ZeroDiversity,
        );
    }
    let Some(expected_frequency) = expected_frequency.map(ExpectedAlternativeFrequency::get) else {
        return SpectrumSeed::new(
            NEUTRAL_ALPHA_REF,
            measured_diversity,
            SeedRegime::NeutralShape,
        );
    };

    let total = total_for_diversity(expected_frequency, measured_diversity);
    // **No floor is applied here, and an earlier draft applied one.** Flooring the alternative
    // concentration at `MIN_ALT_CONCENTRATION` would break the one thing this function
    // guarantees: below a measured diversity of about `2e-12` the floor binds and the pair stops
    // implying the measurement. A diversity of exactly zero is the case that needs the floor and
    // it is taken above; every diversity above zero gives a strictly positive total, and the
    // per-locus expansion floors what it shares out (`fill_locus_concentration`).
    SpectrumSeed::new(
        total * (1.0 - expected_frequency),
        total * expected_frequency,
        SeedRegime::FittedCurve,
    )
}

/// **Spread the run's two numbers over one locus's alleles** — the reference allele's
/// concentration first, then the alternative total shared out evenly across however many
/// alternative alleles this locus turns out to carry.
///
/// ```text
///   out[0]    = α_ref,       floored at MIN_ALT_CONCENTRATION
///   out[1..]  = α_alt_total / (number of alternative alleles),
///                            floored at MIN_ALT_CONCENTRATION
/// ```
///
/// **`allele_count` and `out` are checked against each other, and that is the whole reason the
/// count is an argument at all.** `out` is a slice of scratch the calling loop owns and reuses, so
/// its length is a slicing decision rather than a property of the locus: a locus of three alleles
/// handed an eight-wide buffer at its full length gives every alternative `θ/7` where it should
/// have `θ/2` — 3.5 times too little prior mass on every non-reference genotype, about 5.4 phred,
/// with every value still looking like a value. A buffer *shorter* than the locus is worse: the
/// entries past its end keep the previous locus's concentrations, and every length check
/// downstream still passes.
///
/// **The check is only load-bearing if the two arguments come from different expressions** —
/// `candidate_alleles.len()` for the count against the buffer the worker sliced. Written
/// `fill_locus_concentration(seed, class, n, &mut buffer[..n])` it proves nothing, and no type in
/// this module can tell the difference. It is what `arch/calling_priors.md` §4 asks for, and it
/// costs one integer compare per locus.
///
/// **Sharing the total out rather than giving each allele the whole of it keeps a site's total
/// polymorphism independent of how many alleles it happens to carry** — a triallelic site is not
/// twice as polymorphic as a biallelic one merely for holding a third allele
/// (`doc/devel/ng/spec/calling_priors.md` §4). That is the shape of production's
/// `alpha_from_diversity` ([`genetics.rs`](crate::genetics::alpha_from_diversity)), ported with
/// the fitted pair as input instead of `α_ref = 1` and `α_alt = θ` hard-coded, and filling a
/// caller's buffer instead of returning a fresh `Vec` per locus.
///
/// **The floor is what keeps a rare allele recoverable.** A cohort with no polymorphism at all
/// has an alternative total of zero, and a concentration of exactly zero puts `lgamma(0)` in the
/// prior's row. [`MIN_ALT_CONCENTRATION`] sits at `1e-12`, nine orders of magnitude below the
/// diversities this caller has actually been fitted at — 6 in 10,000 on the tomato panel, 1 in
/// 1,000 on human — so it never perturbs a real estimate.
///
/// `α_ref` is floored by the same constant, and **no measurement can make that bind**: the seed's
/// reference concentration is `A (1 − f)`, and reaching `1e-12` needs either a mean frequency
/// within `1e-12` of one — a population where the reference base is essentially absent — or a
/// heterozygosity below about `2e-12`, which is one difference per 500 thousand million bases. It
/// is there because [`SpectrumSeed::new`] admits any strictly positive
/// reference concentration, so a seed built by hand can sit below it, and
/// [`Concentration`](super::Concentration)'s invariant is that **every** entry clears the floor,
/// not only the alternatives. Production never had to check it: its `α_ref` was the constant 1.
///
/// **A monomorphic locus is an answer, not an error**: a buffer of length one gets the reference
/// concentration and nothing else, which is what `alpha_from_diversity` does with a
/// single-allele shape.
///
/// ## Why it hands back a [`Concentration`] rather than only filling the buffer
///
/// **So that forgetting to call it stops compiling.**
/// [`fill_sample_concentration`](super::fill_sample_concentration) takes the seed as a
/// `Concentration`, and the only ways to get one are this function and
/// [`Concentration::new`](super::Concentration::new) — so a calling loop that skipped this step
/// and passed its raw buffer no longer type-checks.
///
/// It was a bare slice, and the omission was caught by nothing. A buffer of zeros is the right
/// length and its entries are legal floats, so every check downstream passes; the row then
/// reaches `lgamma(0)` and comes back `[NaN, −inf, NaN]`, which is what the
/// [`GenotypePriorModel`](super::GenotypePriorModel) contract forbids by name. Downstream the
/// locus runs to its pass cap and is emitted as unconverged with nothing saying why. **Every
/// mistake about the buffer's shape was caught in release; the omission was caught nowhere**, and
/// `Concentration`'s own per-entry check cannot close it — that one is a `debug_assert!` and
/// release builds compile it out.
///
/// ## The variant class is an argument and is not read
///
/// Both classes take the same seed today, for the reason [`VariantClass`] gives.
///
/// **This is the end that will apply the split when it happens** — settled by the owner,
/// 2026-08-22, and [`seed_from_population_moments`] no longer takes the argument because of it. That
/// function reads the *shape* of variation off the panel's allele counts, which the pre-pass fits
/// without separating the two classes; a class-specific *scale* belongs here, where the run's
/// total is shared out. When the pre-pass measures two diversities, [`SpectrumSeed`] carries both
/// totals and this function picks between them — the run still holds one seed, which is what the
/// calling loop's frozen parameters already assume. Carrying the argument at both ends would
/// apply the ratio twice.
///
/// **One class for the locus, and a locus that mixes them is still Q1's to settle.** Production
/// its pseudocount on each *alternative allele* — `0.01` for a substitution against `0.00125` for
/// an indel ([`DEFAULT_SNP_ALT_PSEUDOCOUNT`](crate::var_calling::posterior_engine::DEFAULT_SNP_ALT_PSEUDOCOUNT),
/// [`DEFAULT_INDEL_ALT_PSEUDOCOUNT`](crate::var_calling::posterior_engine::DEFAULT_INDEL_ALT_PSEUDOCOUNT))
/// — and a generic locus can carry one of each, since
/// [`LocusKind::Generic`](crate::ng::locus_generation::LocusKind) covers both. Taking one class
/// forecloses nothing: an alternative allele's class is readable from the locus's own
/// [`CandidateAlleles`](crate::ng::calling::CandidateAlleles), which hold the bases, so the day
/// Q1 splits the estimate this function can read them without a new array threaded in from the
/// loop.
#[must_use]
pub fn fill_locus_concentration<'a>(
    seed: SpectrumSeed,
    class: VariantClass,
    allele_count: usize,
    out: &'a mut [f64],
) -> Concentration<'a> {
    // Both classes take the same seed today. **This is the end that will apply the split when it
    // happens** — the projection reads the shape of variation off allele counts the pre-pass does
    // not separate by class, so a class-specific scale belongs here (spec §4.2, Q1, settled
    // 2026-08-22).
    let _ = class;

    assert!(
        allele_count > 0,
        "every locus has a reference allele, so its concentration has at least one entry — the \
         caller has lost track of which locus it is on"
    );
    assert_eq!(
        out.len(),
        allele_count,
        "the buffer must cover the locus's alleles exactly: a longer one shares the alternative \
         concentration out too thinly and a shorter one leaves the previous locus's entries \
         behind, and both look like answers"
    );
    out[0] = seed.alpha_ref().max(MIN_ALT_CONCENTRATION);

    let alternative_allele_count = allele_count - 1;
    if alternative_allele_count > 0 {
        // A monomorphic locus does not enter this branch: `out[1..]` is empty there, so the fill
        // would write nothing whatever the division produced. The condition is here to name the
        // case rather than to guard the fill.
        let per_alternative_allele =
            (seed.alpha_alt_total() / alternative_allele_count as f64).max(MIN_ALT_CONCENTRATION);
        out[1..].fill(per_alternative_allele);
    }
    Concentration::new(out)
}

#[cfg(test)]
mod seed_tests {
    use super::*;
    use crate::genetics::lgamma;
    use crate::ng::parameter_estimation::joint::fit::FrequencyDensity;

    /// **What a concentration pair says about how often a diploid drawn from it is
    /// heterozygous** — the Beta-binomial at one alternative copy in two draws, evaluated through
    /// gamma functions.
    ///
    /// **The oracle for `doc/devel/ng/spec/ordinary_site_seed.md` §3, and it shares no line of
    /// arithmetic with what it checks.** The pin solves `A` from `θ = 2 f (1 − f) · A / (A + 1)`;
    /// this evaluates `2 · B(1 + α_alt, 1 + α_ref) / B(α_alt, α_ref)` through `lgamma`, so a test
    /// comparing the two is not comparing a value against its own definition.
    ///
    /// **⚠ It checks the pair's total and is blind to how the total was split.**
    /// `2 α_ref α_alt / (A (A + 1))` is symmetric in the two concentrations, so it returns the
    /// same number for a pair and for its mirror, so a test that reads only this cannot see the
    /// two swapped. Where that matters — in
    /// [`the_seeds_implied_diversity_is_the_measured_one_at_every_shape`], which is the one test
    /// here that goes through the seed builder on a shape whose two concentrations differ — the
    /// expected frequency is asserted separately beside it.
    ///
    /// **It used to read the module's own spectrum machinery at one individual**, which was
    /// deleted with the search on 2026-08-27; this is the same quantity by the same route, written
    /// out here because it is now only a test's oracle.
    fn implied_heterozygosity(seed: SpectrumSeed) -> f64 {
        let (a, b) = (seed.alpha_alt_total(), seed.alpha_ref());
        let ln_beta = |x: f64, y: f64| lgamma(x) + lgamma(y) - lgamma(x + y);
        2.0 * (ln_beta(1.0 + a, 1.0 + b) - ln_beta(a, b)).exp()
    }

    /// The seed the shipped builder returns for a measured mean frequency and a measured
    /// heterozygosity — the whole of what it takes.
    fn seed_at(expected_frequency: f64, diversity: f64) -> SpectrumSeed {
        seed_from_population_moments(
            Some(ExpectedAlternativeFrequency::try_new(expected_frequency).unwrap()),
            Some(ExpectedHeterozygosity::try_new(diversity).unwrap()),
        )
    }

    /// One row of the shape list: a density, and **both of its moments worked out by hand**.
    ///
    /// The two truths are literals rather than expressions, so no formula is shared with what they
    /// check. That is the whole point of the type: without them the list is decoration, and a
    /// review measured exactly that — with the truths absent, deleting the Beta's shape from the
    /// heterozygosity and swapping `a` for `b` in the mean frequency both left every test that
    /// used this list green.
    struct ShapeWithItsMoments {
        name: &'static str,
        density: FrequencyDensity,
        /// `p_fixed_alt + p_segregating · a/(a+b)`, evaluated by hand.
        mean_frequency: f64,
        /// `p_segregating · 2ab/((a+b)(a+b+1))`, evaluated by hand.
        heterozygosity: f64,
    }

    /// The densities the closed-form checks run over: the four shapes
    /// `doc/devel/ng/spec/ordinary_site_seed.md` §1.2 measures, the unit tests' own lopsided
    /// fixture, **and one where the reference base is the rare one at the positions that vary**.
    ///
    /// **That last row is why this is a list rather than one density.** The other five all have
    /// `a ≤ b`, and a formula that read `b / (a + b)` for `a / (a + b)` does one of two things on
    /// them: on the three where `a < b` it comes back too *high*, and on the two symmetric ones,
    /// `Beta(1, 1)` and `Beta(4, 4)`, it comes back **identical** — the swap is not merely
    /// consistent there, it is invisible. `Beta(3, 0.6)` (report §2, the population where the
    /// reference base is the rare one) is the only fixture on which the swap points the other way.
    ///
    /// **Six rows, five distinct mean frequencies.** `Beta(1, 1)` and `Beta(4, 4)` are both
    /// symmetric, so their means are both a half and the two rows give the same answer — 3.000 in
    /// 1,000. They differ in spread, which is what `ordinary_site_seed.md` §1.2 measures them for,
    /// and not in mean. The list is kept whole because it is the same set those measurements used,
    /// and the duplication is named here rather than left for a reader to discover.
    ///
    /// **⚠ The heterozygosity column cannot see an `a`-for-`b` swap on any row, and that is a fact
    /// about the quantity rather than about the fixtures.** `2ab/((a+b)(a+b+1))` is symmetric in
    /// its two arguments. The swap is visible only in the mean frequency, which is why both
    /// columns are here.
    fn shapes_spanning_the_beta_and_both_point_masses() -> [ShapeWithItsMoments; 6] {
        [
            ShapeWithItsMoments {
                name: "tomato-like, strong rare-allele pile-up",
                density: FrequencyDensity {
                    p_invariant: 0.9950,
                    p_fixed_alt: 0.0010,
                    a: 0.20,
                    b: 1.00,
                },
                // 0.0010 + 0.0040 · (1/6); 0.0040 · 0.4/2.64.
                mean_frequency: 0.001_666_666_666_667,
                heterozygosity: 0.000_606_060_606_061,
            },
            ShapeWithItsMoments {
                name: "human-like, moderate pile-up",
                density: FrequencyDensity {
                    p_invariant: 0.9949,
                    p_fixed_alt: 0.0004,
                    a: 0.35,
                    b: 1.20,
                },
                // 0.0004 + 0.0047 · (0.35/1.55); 0.0047 · 0.84/3.9525.
                mean_frequency: 0.001_461_290_322_581,
                heterozygosity: 0.000_998_861_480_076,
            },
            ShapeWithItsMoments {
                name: "flat over what segregates",
                density: FrequencyDensity {
                    p_invariant: 0.9950,
                    p_fixed_alt: 0.0010,
                    a: 1.00,
                    b: 1.00,
                },
                // 0.0010 + 0.0040 · 0.5; 0.0040 · 2/6.
                mean_frequency: 0.003,
                heterozygosity: 0.001_333_333_333_333,
            },
            ShapeWithItsMoments {
                name: "the unit tests' own lopsided fixture",
                density: FrequencyDensity {
                    p_invariant: 0.90,
                    p_fixed_alt: 0.01,
                    a: 0.50,
                    b: 2.00,
                },
                // 0.01 + 0.09 · 0.2; 0.09 · 2.0/8.75.
                mean_frequency: 0.028,
                heterozygosity: 0.020_571_428_571_429,
            },
            ShapeWithItsMoments {
                name: "middling frequencies — the shape the family cannot hold",
                density: FrequencyDensity {
                    p_invariant: 0.9950,
                    p_fixed_alt: 0.0010,
                    a: 4.00,
                    b: 4.00,
                },
                // 0.0010 + 0.0040 · 0.5; 0.0040 · 32/72.
                mean_frequency: 0.003,
                heterozygosity: 0.001_777_777_777_778,
            },
            ShapeWithItsMoments {
                name: "where it varies, the reference base is the rare one",
                density: FrequencyDensity {
                    p_invariant: 0.9950,
                    p_fixed_alt: 0.0010,
                    a: 3.00,
                    b: 0.60,
                },
                // 0.0010 + 0.0040 · (5/6); 0.0040 · 3.6/16.56.
                mean_frequency: 0.004_333_333_333_333,
                heterozygosity: 0.000_869_565_217_391,
            },
        ]
    }

    /// **Both closed forms are what a hand calculation gives, on all six shapes.**
    ///
    /// This is the test that makes the shape list load-bearing: it compares each density's two
    /// integrals against literals worked out from the four fitted numbers, so it shares no
    /// expression with either function. Without it the list's own justification was untrue — a
    /// review measured that deleting `a · b` from the heterozygosity and swapping `a` for `b` in
    /// the mean frequency both left every test that used the list green.
    ///
    /// The tolerance is relative and `1e-12`, which is about the accumulated rounding in
    /// `1 − p_invariant − p_fixed_alt` and far tighter than any of the defects above.
    #[test]
    fn both_closed_forms_are_what_a_hand_calculation_gives() {
        for shape in shapes_spanning_the_beta_and_both_point_masses() {
            let frequency = shape.density.expected_alternative_frequency();
            assert!(
                (frequency / shape.mean_frequency - 1.0).abs() < 1e-12,
                "on {} the mean frequency is {frequency:e} where the hand calculation gives {:e}",
                shape.name,
                shape.mean_frequency
            );
            let heterozygosity = shape.density.expected_heterozygosity();
            assert!(
                (heterozygosity / shape.heterozygosity - 1.0).abs() < 1e-12,
                "on {} the heterozygosity is {heterozygosity:e} where the hand calculation gives \
                 {:e}",
                shape.name,
                shape.heterozygosity
            );
        }
    }

    /// **No fitted mean frequency: the pair is the neutral `(1, θ)` at the diversity the pre-pass
    /// did fit** — exactly, with no arithmetic in between, and the regime says where it came from.
    ///
    /// A run arrives without a mean frequency when its pre-pass fitted no population curve: the
    /// per-sample histogram route, which supplies a diversity and no density, is the case the
    /// design has in mind. **A branch on absence, never on cohort size.**
    #[test]
    fn no_fitted_frequency_is_the_neutral_pair_at_the_fitted_diversity() {
        let theta = ExpectedHeterozygosity::try_new(6e-4).unwrap();
        let seed = seed_from_population_moments(None, Some(theta));
        assert_eq!(seed.alpha_ref(), 1.0);
        assert_eq!(seed.alpha_alt_total(), 6e-4);
        assert_eq!(seed.regime(), SeedRegime::NeutralShape);
    }

    /// **No fitted diversity either: the species-range fallback, and the run must say so.** Two
    /// runs that used different information are otherwise indistinguishable in what they emit.
    #[test]
    fn no_fitted_diversity_falls_back_to_the_species_value_and_says_so() {
        let seed = seed_from_population_moments(None, None);
        assert_eq!(seed.alpha_ref(), 1.0);
        assert_eq!(
            seed.alpha_alt_total(),
            ExpectedHeterozygosity::SPECIES_FALLBACK.get()
        );
        assert_eq!(seed.regime(), SeedRegime::FallbackDiversity);
    }

    /// **A fitted mean frequency with no diversity beside it is discarded**, and the run lands on
    /// the species-range guess and says so.
    ///
    /// The total comes from the heterozygosity, so a frequency alone has nothing to pin. On the
    /// joint route the two always arrive together — both are integrals of one fitted curve — so
    /// this is a shape the type system allows rather than a state a run reaches.
    #[test]
    fn a_frequency_with_no_diversity_falls_back_and_says_so() {
        let seed = seed_from_population_moments(
            Some(ExpectedAlternativeFrequency::try_new(0.2).unwrap()),
            None,
        );
        assert_eq!(seed.alpha_ref(), 1.0);
        assert_eq!(
            seed.alpha_alt_total(),
            ExpectedHeterozygosity::SPECIES_FALLBACK.get()
        );
        assert_eq!(seed.regime(), SeedRegime::FallbackDiversity);
    }

    /// **Two numbers that did not come from one population curve are refused, not answered.**
    ///
    /// A pair of expected frequency `f` makes a diploid heterozygous at most `2 f (1 − f)` of the
    /// time, however much conviction it carries, so a heterozygosity at or above that ceiling has
    /// no answer at all. **Before 2026-08-27 the seed fell back to the neutral pair and reported
    /// which of two ways it had got there.** It no longer can: on the shipped route both numbers
    /// are integrals of one curve and Jensen's inequality puts the measurement strictly below the
    /// ceiling, so reaching this state means a caller assembled two numbers that do not go
    /// together.
    ///
    /// **The fixture is the state a fully invariant panel used to produce**: an expected frequency
    /// of 1 in a thousand million — the old search's own floor on the ratio between the two
    /// concentrations — against a measured heterozygosity of 6 in 10,000. The ceiling there is
    /// about 2 in a thousand million, five orders of magnitude short.
    #[test]
    #[should_panic(expected = "did not come from one curve")]
    fn two_moments_that_cannot_belong_to_one_curve_are_refused() {
        let _ = seed_at(1e-9, 6e-4);
    }

    /// **The refusal is at the ceiling and not above it**, which is the difference between a
    /// reported refusal and a panic three frames later about a concentration: at exactly the
    /// ceiling the solved total is infinite, and `SpectrumSeed` refuses a non-finite one.
    ///
    /// At an expected frequency of a half the ceiling is a half, and this asks for exactly that.
    #[test]
    #[should_panic(expected = "did not come from one curve")]
    fn a_heterozygosity_exactly_at_the_ceiling_is_refused() {
        let _ = seed_at(0.5, 0.5);
    }

    /// **And a hair below the ceiling is not refused** — the bound is at the ceiling rather than
    /// near it, which neither of the two tests above can say on its own.
    ///
    /// 999 parts in a thousand of the ceiling gives a total near a thousand: a pair that carries a
    /// thousand chromosomes' worth of conviction, which is what asking for almost all of a shape's
    /// own maximum heterozygosity costs.
    #[test]
    fn just_below_the_ceiling_still_has_a_total() {
        let frequency = 1e-3;
        let ceiling = 2.0 * frequency * (1.0 - frequency);
        let seed = seed_at(frequency, 0.999 * ceiling);
        let total = seed.alpha_ref() + seed.alpha_alt_total();
        assert!(
            (900.0..1_100.0).contains(&total),
            "999 parts in a thousand of the ceiling needs a total near a thousand; got {total}"
        );
        assert!(matches!(seed.regime(), SeedRegime::FittedCurve));
    }

    /// **The seed's implied diversity is the one that was measured, at every shape** —
    /// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §9's third test, and goal 1 of the whole
    /// change.
    ///
    /// **Both of the seed's inputs are now integrals of one fitted curve**, so this runs over the
    /// six shapes those integrals are checked on rather than over a grid of panel sizes: there is
    /// no panel size left in the seed to vary. Where the old two-parameter fit lost 9.9% of the
    /// diversity at 63 individuals on a tomato-like shape, 18.6% on a human-like one and 53.9% on
    /// a middling one (`ordinary_site_seed.md` §1.2), the pinned pair loses none of it on any of
    /// these.
    ///
    /// **Asserted rather than sampled, and it is not a restatement of the identity that produced
    /// it**: the seed's own implied heterozygosity is read off [`implied_heterozygosity`], which
    /// evaluates the Beta-binomial sum at one individual and shares no algebra with the pin.
    #[test]
    fn the_seeds_implied_diversity_is_the_measured_one_at_every_shape() {
        let mut worst: f64 = 0.0;
        let mut worst_at = "";
        for shape in shapes_spanning_the_beta_and_both_point_masses() {
            let name = shape.name;
            // **The hand-computed moments, not the density's own methods.** Handing the seed
            // builder `density.expected_alternative_frequency()` and then asserting the seed
            // reports that same number back is an identity: it can see the builder swap the pair,
            // and it cannot see the moment function be wrong. These two literals can.
            let seed = seed_at(shape.mean_frequency, shape.heterozygosity);
            assert!(
                matches!(seed.regime(), SeedRegime::FittedCurve),
                "on {name} the regime came back {:?}",
                seed.regime()
            );
            // **The pair's own expected frequency is the one it was handed.**
            // [`implied_heterozygosity`] cannot see this: `2 α_ref α_alt / (A (A + 1))` is
            // symmetric in the two concentrations, so it returns the same number for a pair and
            // for its mirror.
            let total = seed.alpha_ref() + seed.alpha_alt_total();
            let frequency = seed.alpha_alt_total() / total;
            assert!(
                (frequency / shape.mean_frequency - 1.0).abs() < 1e-12,
                "on {name} the seed's expected frequency is {frequency:e} where it was handed {:e}",
                shape.mean_frequency
            );
            let off = (implied_heterozygosity(seed) / shape.heterozygosity - 1.0).abs();
            if off > worst {
                worst_at = name;
            }
            worst = worst.max(off);
        }
        assert!(
            worst < 1e-11,
            "the seed must imply the diversity it was handed, whatever the population's shape; \
             worst departure {worst:.2e}, on {worst_at}"
        );
    }

    /// **No density the fit can produce comes within one part in a hundred of its own ceiling**,
    /// which is what makes [`total_for_diversity`]'s refusal a release assertion about the caller
    /// rather than a state a run can reach.
    ///
    /// A pair of expected frequency `f` implies at most `2 f (1 − f)`. Jensen's inequality puts a
    /// curve's own heterozygosity below that by twice the spread of its frequencies — this sweeps
    /// the box the fit clamps its four parameters to and reports how much of the ceiling the worst
    /// density in it asks for.
    ///
    /// **The answer is 100 in 101, and where it is reached says why.** With every position
    /// segregating and neither point mass carrying anything, the density is a bare `Beta(a, b)`,
    /// the share of the ceiling is exactly `(a + b) / (a + b + 1)`, and the solved total is exactly
    /// `a + b` — so the tightest case in the box is `a = b = 50` and the total never exceeds 100
    /// chromosomes. Every point mass added widens the margin.
    #[test]
    fn no_density_the_fit_can_produce_comes_within_one_part_in_a_hundred_of_its_ceiling() {
        // The fit's own clamps: `a` and `b` to [0.02, 50] independently, and the two masses to
        // shares that leave something segregating (`joint::fit`'s M-step).
        let shapes = [0.02_f64, 0.1, 0.5, 1.0, 4.0, 20.0, 50.0];
        let masses = [0.0_f64, 1e-9, 1e-3, 0.1, 0.5, 0.9, 0.999];
        let mut worst_share = 0.0_f64;
        let mut worst_at = (0.0, 0.0, 0.0, 0.0);
        for &a in &shapes {
            for &b in &shapes {
                for &p_invariant in &masses {
                    for &p_fixed_alt in &masses {
                        if p_invariant + p_fixed_alt >= 1.0 {
                            continue;
                        }
                        let density = FrequencyDensity {
                            p_invariant,
                            p_fixed_alt,
                            a,
                            b,
                        };
                        let frequency = density.expected_alternative_frequency();
                        let theta = density.expected_heterozygosity();
                        let share = theta / (2.0 * frequency * (1.0 - frequency));
                        if share > worst_share {
                            worst_share = share;
                            worst_at = (a, b, p_invariant, p_fixed_alt);
                        }
                    }
                }
            }
        }
        assert!(
            worst_share <= 100.0 / 101.0 + 1e-12,
            "the worst density in the fit's own box asks for {worst_share} of its ceiling, at \
             Beta({}, {}) with masses {} and {} — the closed form says the bare Beta at a = b = 50 \
             is the tightest, at 100/101",
            worst_at.0,
            worst_at.1,
            worst_at.2,
            worst_at.3
        );
        // **And it really does get that close**, so the bound above is tight rather than generous:
        // a test asserting only `share < 1` would pass with the clamps ten times wider.
        assert!(
            worst_share > 100.0 / 101.0 - 1e-9,
            "the bare Beta at a = b = 50 should reach 100/101 exactly; the sweep's worst is \
             {worst_share}"
        );
    }

    /// **The pin's algebra is the Beta-binomial's own answer** — checked against a bisection on
    /// [`implied_heterozygosity`] over a grid of frequencies and diversities.
    ///
    /// [`total_for_diversity`] solves `A = t / (1 − t)` in closed form; the oracle searches for
    /// the `A` whose spectrum at one individual puts the measured mass on the heterozygous class.
    /// **They share no line of arithmetic**, which is what makes the comparison worth making.
    #[test]
    fn the_solved_total_is_what_the_beta_binomial_needs() {
        let mut worst: f64 = 0.0;
        for expected_frequency in [1e-4_f64, 1e-3, 6e-3, 0.05, 0.3] {
            for share_of_ceiling in [1e-3_f64, 0.1, 0.5, 0.9, 0.99] {
                let ceiling = 2.0 * expected_frequency * (1.0 - expected_frequency);
                let diversity = share_of_ceiling * ceiling;
                let total = total_for_diversity(expected_frequency, diversity);
                let seed = SpectrumSeed::new(
                    total * (1.0 - expected_frequency),
                    total * expected_frequency,
                    SeedRegime::NeutralShape,
                );
                worst = worst.max((implied_heterozygosity(seed) / diversity - 1.0).abs());
            }
        }
        assert!(
            worst < 1e-11,
            "the solved total must be the one the Beta-binomial needs; worst {worst:.2e}"
        );
    }

    /// **A run that measured no variation at all is floored and says so.**
    ///
    /// Solving for the total at a diversity of zero gives zero, and every entry of the pair with
    /// it — a pair `SpectrumSeed` would refuse, since a reference concentration of zero reaches
    /// `lgamma(0)` downstream. So the alternative concentration takes the floor every other seed
    /// builder in this tree applies and the regime says the diversity was zero
    /// (`ordinary_site_seed.md` §3.1).
    ///
    /// **A fitted mean frequency makes no difference here**, and that is the point of the variant
    /// carrying none: with no diversity there is nothing for a shape to scale.
    #[test]
    fn a_cohort_with_no_variation_is_floored_and_says_the_diversity_was_zero() {
        for frequency in [ExpectedAlternativeFrequency::try_new(0.2).ok(), None] {
            let seed = seed_from_population_moments(
                frequency,
                Some(ExpectedHeterozygosity::try_new(0.0).unwrap()),
            );
            assert_eq!(seed.regime(), SeedRegime::ZeroDiversity);
            assert_eq!(seed.alpha_ref(), 1.0);
            assert_eq!(seed.alpha_alt_total(), MIN_ALT_CONCENTRATION);
            // The floored pair is one the per-locus expansion accepts and `lgamma` can take.
            assert!(seed.alpha_alt_total() > 0.0);
        }
    }

    /// **A heterozygosity above a half refuses the run rather than seeding from it.**
    ///
    /// A pair of expected frequency `f` makes a diploid heterozygous `2 f (1 − f) · A / (A + 1)`
    /// of the time, and `2 f (1 − f)` is at most a half. `ExpectedHeterozygosity` admits the whole
    /// of `[0, 1]`, so a fit that did not converge can hand the caller 0.9 — and the difference
    /// between refusing at the run's assembly and falling back is that the second names no
    /// culprit (`ordinary_site_seed.md` §3.1).
    #[test]
    #[should_panic(expected = "is not a thin estimate")]
    fn a_heterozygosity_above_a_half_refuses_the_run() {
        let _ =
            seed_from_population_moments(None, Some(ExpectedHeterozygosity::try_new(0.9).unwrap()));
    }

    /// **The refusal sits at exactly a half and not somewhere above it**, which neither of its
    /// two neighbours can say: one refuses at 0.9 and one accepts at exactly 0.5, so a bound
    /// moved to 0.8 keeps both of them green. Measured — it does.
    ///
    /// A millionth above the bound is the smallest step that is unambiguously outside it and
    /// still exactly representable well clear of the comparison's own precision.
    #[test]
    #[should_panic(expected = "is not a thin estimate")]
    fn a_heterozygosity_a_millionth_above_a_half_refuses_the_run() {
        let _ = seed_from_population_moments(
            None,
            Some(ExpectedHeterozygosity::try_new(0.500_001).unwrap()),
        );
    }

    /// **How close a measurement can come to the ceiling and still have an answer** — one bit,
    /// where the solved total is 9.0 × 10¹⁵ and the pair is a prior no depth of reads could move.
    ///
    /// It is not reachable from a fit: it needs the measurement and the shape's own ceiling to
    /// agree to one part in 10¹⁶, where Jensen's inequality puts them apart by twice the
    /// population's spread of frequencies. It is recorded rather than guarded, because clamping it
    /// would break the pin for a case that cannot arise.
    #[test]
    fn one_bit_below_the_ceiling_still_has_a_total() {
        let ceiling = 2.0_f64 * 0.3 * (1.0 - 0.3);
        let a_bit_under = f64::from_bits(ceiling.to_bits() - 1);
        let enormous = total_for_diversity(0.3, a_bit_under);
        assert!(
            enormous > 8e15 && enormous.is_finite(),
            "one bit below the ceiling the total is {enormous:e}"
        );
    }

    /// **The pin holds at every measured diversity above zero, however small** — there is no
    /// floor on this path, and an earlier draft of it had one.
    ///
    /// A floor on the alternative concentration would bind below a diversity of about `2e-12` and
    /// the pair would stop implying the measurement there. The case that genuinely needs a floor
    /// is a diversity of exactly zero, which is taken before any of this.
    #[test]
    fn a_diversity_far_below_the_floor_is_still_pinned_rather_than_floored() {
        let theta = 1e-15;
        let seed = seed_at(6e-4 / (1.0 + 6e-4), theta);
        assert!(
            matches!(seed.regime(), SeedRegime::FittedCurve),
            "got {:?}",
            seed.regime()
        );
        assert!(
            seed.alpha_alt_total() < MIN_ALT_CONCENTRATION,
            "a floored alternative concentration would sit at {MIN_ALT_CONCENTRATION:e}; got {}",
            seed.alpha_alt_total()
        );
        assert!(
            (implied_heterozygosity(seed) / theta - 1.0).abs() < 1e-6,
            "the pair must still imply {theta:e}; it implies {}",
            implied_heterozygosity(seed)
        );
    }
}

#[cfg(test)]
mod locus_concentration_tests {
    use super::*;
    use crate::ng::calling::genotype_prior::Concentration;

    fn neutral_seed(theta: f64) -> SpectrumSeed {
        SpectrumSeed::new(1.0, theta, SeedRegime::NeutralShape)
    }

    fn locus_concentration(seed: SpectrumSeed, allele_count: usize) -> Vec<f64> {
        let mut out = vec![f64::NAN; allele_count];
        let _ = fill_locus_concentration(seed, VariantClass::Substitution, allele_count, &mut out);
        out
    }

    /// **A locus carries the same total polymorphism however many alleles it has** — the property
    /// spec §4 gives as the reason for sharing the total out rather than repeating it.
    ///
    /// Without it a site is more polymorphic for having been given a third allele to consider,
    /// which is a statement about the candidate generator rather than about the genome.
    ///
    /// Held across four allele counts and three diversities. **The bound is relative rather than
    /// absolute**, because the total is not: `SpectrumSeed` admits anything finite and the fit's
    /// own box reaches 1e3, where an absolute `1e-15` would fail on a legal answer — the sum's
    /// worst-case rounding over 9 alternatives at a total of 1e3 is 1.1e-13. At the diversities
    /// in this loop the relative form is four orders *tighter* than `1e-15` would have been, and
    /// measured, all twelve cells come back exactly 0.0.
    #[test]
    fn a_locus_carries_the_same_total_polymorphism_however_many_alleles_it_has() {
        for theta in [1e-4, 6e-4, 1e-2] {
            let seed = neutral_seed(theta);
            for allele_count in [2, 3, 4, 6] {
                let out = locus_concentration(seed, allele_count);
                let alternative_total: f64 = out[1..].iter().sum();
                assert!(
                    (alternative_total - theta).abs() <= 8.0 * f64::EPSILON * theta,
                    "{allele_count} alleles at θ = {theta} carry {alternative_total} of \
                     alternative concentration between them"
                );
                assert_eq!(
                    out[0], 1.0,
                    "{allele_count} alleles at θ = {theta}: the reference entry should be the \
                     run's α_ref of 1.0, got {}",
                    out[0]
                );
            }
        }
    }

    /// **The reference allele takes the first entry**, because the concentration is read in the
    /// same order as the locus's candidate alleles and entry 0 is the reference's
    /// (`doc/devel/ng/spec/calling_priors.md` §4). Reversing the two would leave every genotype
    /// row believing the reference is the rare allele.
    #[test]
    fn the_reference_allele_takes_the_first_entry() {
        let out = locus_concentration(SpectrumSeed::new(1.0, 6e-4, SeedRegime::NeutralShape), 3);
        assert_eq!(out[0], 1.0);
        assert_eq!(out[1], 3e-4);
        assert_eq!(out[2], 3e-4);
    }

    /// **A locus with no alternative allele gets only the reference** — an answer rather than an
    /// error, and what `alpha_from_diversity` does with a single-allele shape.
    #[test]
    fn a_locus_with_no_alternative_allele_gets_only_the_reference() {
        assert_eq!(locus_concentration(neutral_seed(6e-4), 1), vec![1.0]);
    }

    /// **A cohort with no polymorphism still clears the floor**, so the prior's row never reaches
    /// `lgamma(0)` and a rare allele stays recoverable rather than falling into a zero it can
    /// never leave.
    ///
    /// The floor cannot perturb a real estimate: at a `θ` of 1 in 10,000 shared between two
    /// alternatives each gets 5e-5, which is 50 million times `MIN_ALT_CONCENTRATION`.
    #[test]
    fn a_cohort_with_no_polymorphism_still_clears_the_floor() {
        let out = locus_concentration(SpectrumSeed::new(1.0, 0.0, SeedRegime::NeutralShape), 3);
        assert_eq!(out[1], MIN_ALT_CONCENTRATION);
        assert_eq!(out[2], MIN_ALT_CONCENTRATION);
        let real = locus_concentration(neutral_seed(1e-4), 3);
        assert!(real[1] > 1e7 * MIN_ALT_CONCENTRATION);
    }

    /// **Every entry clears the floor [`Concentration`] requires**, the reference's included —
    /// which production never had to check, because its `α_ref` was the constant 1 and only the
    /// alternatives could reach zero.
    ///
    /// **The check is written out here rather than left to `Concentration::new`**, whose
    /// per-entry test is a `debug_assert!` and so is compiled out of the `--release` runs this
    /// suite is gated on (`Cargo.toml` sets no `debug-assertions` on `[profile.release]`).
    /// Constructing the type as well keeps the two in step if the invariant moves.
    ///
    /// The `1e-30` seed here is hand-built: no fit can produce one, since the search box bottoms
    /// out near `α_ref = 1e-5`. What it guards is that
    /// [`SpectrumSeed::new`] admits any strictly positive reference concentration.
    #[test]
    fn every_entry_clears_the_floor_the_concentration_type_requires() {
        for seed in [
            SpectrumSeed::new(1.0, 6e-4, SeedRegime::NeutralShape),
            SpectrumSeed::new(1e-30, 0.0, SeedRegime::NeutralShape),
            SpectrumSeed::new(500.0, 500.0, SeedRegime::NeutralShape),
        ] {
            let out = locus_concentration(seed, 4);
            assert!(
                out.iter()
                    .all(|a| a.is_finite() && *a >= MIN_ALT_CONCENTRATION),
                "a seed of ({}, {}) gave {out:?}, which Concentration refuses",
                seed.alpha_ref(),
                seed.alpha_alt_total()
            );
            let _ = Concentration::new(&out);
        }
    }

    /// **Both variant classes get the same concentration today**, which is what makes the
    /// argument a seam rather than a behaviour — see [`fill_locus_concentration`] for which end
    /// of the pipeline would absorb a split.
    #[test]
    fn both_variant_classes_get_the_same_concentration_today() {
        let seed = neutral_seed(6e-4);
        let mut substitution = vec![f64::NAN; 3];
        let mut indel = vec![f64::NAN; 3];
        let _ = fill_locus_concentration(seed, VariantClass::Substitution, 3, &mut substitution);
        let _ = fill_locus_concentration(seed, VariantClass::InsertionOrDeletion, 3, &mut indel);
        assert_eq!(substitution, indel);
    }

    /// **A locus with no alleles is refused.** Every locus has a reference allele, so a count of
    /// zero is a caller that has lost track of which locus it is on. Without the assertion the
    /// write to `out[0]` panics anyway, on an index — the assertion is what makes the panic name
    /// the mistake rather than the line.
    #[test]
    #[should_panic(expected = "lost track of which locus")]
    fn a_locus_with_no_alleles_is_refused() {
        let _ =
            fill_locus_concentration(neutral_seed(6e-4), VariantClass::Substitution, 0, &mut []);
    }

    /// **A buffer that does not cover the locus's alleles exactly is refused, both ways.**
    ///
    /// This is the mistake the calling loop is most likely to make, because `out` is a slice of a
    /// scratch buffer it owns and reuses: handing over the whole buffer rather than the locus's
    /// prefix. Measured with the check removed, a 3-allele locus in an 8-wide buffer at tomato's
    /// fitted `θ` of 6 in 10,000 gives each alternative 8.571e-5 instead of 3.0e-4 — 3.5 times
    /// too little, about 5.4 phred off every non-reference genotype, with every value still
    /// looking like a value. Too short is worse: the entries past the buffer's end keep the
    /// previous locus's concentrations, and every length check downstream passes.
    #[test]
    #[should_panic(expected = "cover the locus's alleles exactly")]
    fn a_buffer_longer_than_the_locus_is_refused() {
        let mut worker_buffer = vec![f64::NAN; 8];
        let _ = fill_locus_concentration(
            neutral_seed(6e-4),
            VariantClass::Substitution,
            3,
            &mut worker_buffer,
        );
    }

    /// The other half of the pair above — see it for what a short buffer costs.
    #[test]
    #[should_panic(expected = "cover the locus's alleles exactly")]
    fn a_buffer_shorter_than_the_locus_is_refused() {
        let mut too_short = vec![f64::NAN; 2];
        let _ = fill_locus_concentration(
            neutral_seed(6e-4),
            VariantClass::Substitution,
            3,
            &mut too_short,
        );
    }
}
