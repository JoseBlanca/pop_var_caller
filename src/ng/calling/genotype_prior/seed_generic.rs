//! The SNP/indel starting point: the run's two concentration numbers — the chromosomes'
//! worth of prior belief attached to the reference allele and to the alternatives — read
//! off the pre-pass's fitted frequency spectrum.
//!
//! At an ordinary site most alternative alleles are rare, so the one chromosome the reference
//! records is almost always the common one and the reference's number is the larger — but only
//! just, and *how much* larger is fitted rather than fixed
//! (`doc/devel/ng/spec/calling_priors.md` §4.1).
//!
//! **Three functions and the view between them.** [`fill_expected_spectrum`] says what a panel's
//! allele counts would look like if a candidate pair were the truth; [`project_spectrum_seed`]
//! takes the pre-pass's measured spectrum as a [`FittedSpectrum`], searches for the pair whose
//! prediction matches it, and hands back the run's seed; [`fill_locus_concentration`] spreads
//! that one pair over each locus's own alleles.

use crate::genetics::{MIN_ALT_CONCENTRATION, PROBABILITY_FLOOR, lgamma};
use crate::ng::parameter_estimation::fitting::multistart::SearchPrecision;
use crate::ng::types::{ExpectedHeterozygosity, InbreedingF};

use super::{Concentration, SeedRegime, SpectrumMatch, SpectrumSeed};

/// The largest concentration — chromosomes' worth of prior belief — this projection will predict a
/// spectrum for.
///
/// **Not a modelling limit but an arithmetic one.** The sum's terms are ratios of gamma functions
/// at arguments that grow with the concentration, and past this point they stop agreeing with each
/// other: measured at 63 individuals and `F = 0.9`, the classes total 1.0000016 at `1e9`, 1.0022 at
/// `1e12` and 1,107 at `1e15`, every entry still finite and non-negative. A run that reached those
/// values would get a spectrum that looks like one.
///
/// `1e6` chromosomes is five hundred thousand diploid individuals, two orders of magnitude past
/// the several thousand this caller commits to, and the total is within 2e-11 of one there.
pub const MAX_PROJECTION_CONCENTRATION: f64 = 1e6;

/// How far below the most likely one a branch split's probability may fall before it is dropped
/// from the sum.
///
/// **Not an approximation in the usual sense: the error has a bound rather than an estimate.** For
/// a fixed number of inbred individuals the classes are themselves a distribution summing to one,
/// so dropping that split moves no class by more than its own probability, and the whole spectrum
/// by no more than the dropped tail.
///
/// **Measured against the term-by-term sum, at tomato's fitted diversity and `F = 0.8`, this value
/// costs nothing at all**: the worst class-by-class disagreement is 4e-13 at 800 and 1,600
/// individuals and 1e-12 at 3,200, which is the same disagreement the untrimmed recurrence has —
/// floating-point accumulation, not the trim. Loosening it does start to cost, and in proportion:
/// `1e-8` gives 1e-9 and `1e-6` gives 1e-7. It buys 3 to 4-fold over not trimming at all
/// (`examples/ng_spectrum_projection_cost.rs`).
const BRANCH_TAIL_TOLERANCE: f64 = 1e-18;

/// A branch split whose probability is below this cannot move any class by more than itself, so it
/// is skipped whatever [`BRANCH_TAIL_TOLERANCE`] says — it is below what a `f64` total can carry.
const NEGLIGIBLE_BRANCH_WEIGHT: f64 = 1e-300;

/// Fill `out` with **what fraction of sites carry the alternative allele on exactly `j` of a
/// panel's chromosomes**, for `j = 0` up to every chromosome, if the population's allele
/// frequencies were drawn from the concentration pair given.
///
/// This is the *frequency spectrum* — the shape the pre-pass measures on real data, and the thing
/// step D2 matches a candidate pair against. `out` holds one entry per class and they sum to one:
/// entry 0 is the fraction of sites where no chromosome carries the alternative allele, entry 1
/// singletons, entry 2 doubletons, and so on.
///
/// ## Where the shape comes from
///
/// Two draws, in order. First the locus's alternative-allele frequency `p` is drawn from the
/// concentration pair — a Beta distribution with those two numbers. Then the panel's genotypes are
/// drawn at that `p` under the **same two-branch model the genotype prior itself uses**
/// (`doc/devel/ng/spec/calling_priors.md` §3.2): with probability `F` an individual's two copies
/// are one ancestral copy counted twice, otherwise they are two independent draws. Averaging over
/// `p` gives the numbers below.
///
/// **Nothing is simulated.** The whole calculation is a finite sum of Beta-function ratios and
/// binomial coefficients — no draws, no quadrature — so the same inputs give the same answer to
/// the last bit at any panel size. The one departure from exactness is the branch-tail trim
/// described under cost, whose error is bounded and, at the tolerance shipped, not measurable.
///
/// **That is what lets step D2's tests state a target rather than only a tolerance** (spec §12
/// tests 5–7): a neutral panel's spectrum can be written down and the fit asked to return `(1, θ)`
/// from it. Note what this does *not* buy — the same function supplies D2's objective, so the
/// tests are not an independent check of the fit's mathematics, only of its search. The
/// independent checks live here, in this file's own oracles.
///
/// ## Why the two-branch model and not `2N` independent chromosomes
///
/// **Because a selfer's panel does not produce an independent-chromosome spectrum at any
/// diversity, and on tomato that is the dominant feature rather than a correction.** In the cohort
/// VCF over 26 accessions, 10,786 sites carry the alternative allele on exactly two chromosomes
/// against 5,142 on exactly one — doubletons outnumbering singletons 2.1 to 1. No
/// independent-chromosome spectrum can do that at any `θ`, because its classes fall monotonically
/// like `θ/k`. Inbreeding is the shape of that spectrum, not a second-order effect on it (spec
/// §4.1). Treating the chromosomes as independent biases `α_ref` **down** with a fixed sign — 12 to
/// 14% at tomato's fitted `F` of 0.8 to 0.9, and 8.6% as far down as `F = 0.6` — which is the
/// number the whole of step D is about.
///
/// ## How it is computed
///
/// Split the panel by which branch each individual took. With `M` of the `N` individuals
/// identical-by-descent, the panel holds `2N − M` **distinct** chromosomes — the other `M` are
/// copies of one of them — and each distinct chromosome carries the alternative allele
/// independently at `p`. So the class probability is a sum over `M`, over how many distinct
/// chromosomes are alternative, and over how many of the duplicated ones are among them.
///
/// **Every term is a product of non-negative factors, and that is the reason for this shape.**
/// Written instead as a polynomial in `p`, the coefficients alternate in sign and grow about
/// 8.3-fold per individual — measured at `F = 0`, they reach 9.5e9 by twelve individuals — so
/// they would cancel catastrophically against a total that must come out at one. These do not
/// cancel at all.
///
/// ## What it costs
///
/// **Paid once per *fit*, not once per run**: this function is the objective step D2 searches, so a
/// multistart fit evaluates it on the order of a hundred times. Measured in release at tomato's
/// fitted diversity and `F = 0.8`, one prediction:
///
/// ```text
///   samples      400     800    1600    3200
///   this        5.8ms  29.9ms   179ms   960ms
///   term-by-term 43.8ms  340ms   2.1s   12.1s
/// ```
///
/// About `N^2.45`, against `N^2.95` for the straight term-by-term sum — so a fit at the top of the
/// committed cohort range is minutes rather than half an hour. Three things buy that, and none of
/// them changes the model:
///
/// - the log-factorials are tabulated once per call rather than rebuilt from `lgamma` per term;
/// - each term is written as **a beta-binomial weight times a hypergeometric one**, both genuine
///   probabilities, so the hypergeometric can be stepped by an exact ratio instead of exponentiated
///   — one exponential per `(split, draw count)` pair rather than one per term;
/// - branch splits far out in their own tail are dropped, at [`BRANCH_TAIL_TOLERANCE`], whose
///   error is bounded rather than estimated.
///
/// **The hypergeometric walk starts at its mode and goes out both ways, and that is not a
/// refinement.** Started at the low end, the first weight underflows to zero long before the ones
/// at the mode become small, and since the walk is multiplicative the whole row then vanishes.
/// Measured at 1,600 individuals when it was written that way: one class came back 5.7e-16 against
/// its true 6.1e-7, and the spectrum lost 3 parts in 10,000 of its mass — with every entry still
/// finite and non-negative, so nothing downstream could have told.
///
/// ## Preconditions and edges
///
/// `out` must hold exactly `2 × individuals + 1` entries, held in release. The panel is diploid:
/// the two-branch model has no state for an individual with two of four copies identical by
/// descent, which spec §3.3 defers to a spec of its own.
///
/// **An alternative concentration of exactly zero is a real answer**, not a degenerate one — a
/// fully invariant cohort — and it puts all the mass in class 0. It is taken early, because the
/// general sum would otherwise reach `lgamma(0)`.
///
/// **Both concentrations are bounded above**, at [`MAX_PROJECTION_CONCENTRATION`]. Past it the sum
/// stops being computable and says nothing about it: measured at 63 individuals and `F = 0.9`, the
/// classes total 1.0000016 at a reference concentration of `1e9`, 1.0022 at `1e12`, and **1,107 at
/// `1e15`** — with every entry still finite and non-negative, so nothing downstream could tell.
/// Step D2 searches over exactly these two axes, which is why this is a release check and not a
/// remark.
pub fn fill_expected_spectrum(
    alpha_ref: f64,
    alpha_alt: f64,
    individuals: u32,
    inbreeding: InbreedingF,
    out: &mut [f64],
) {
    let n = individuals as usize;
    assert_eq!(
        out.len(),
        2 * n + 1,
        "one entry per allele-count class: a panel of {individuals} diploid individuals has \
         {} classes, and `out` holds {}",
        2 * n + 1,
        out.len()
    );
    assert!(
        alpha_ref.is_finite() && alpha_ref > 0.0 && alpha_ref <= MAX_PROJECTION_CONCENTRATION,
        "the reference concentration must be finite, strictly positive and at most \
         {MAX_PROJECTION_CONCENTRATION:e} chromosomes, got {alpha_ref}"
    );
    assert!(
        alpha_alt.is_finite() && (0.0..=MAX_PROJECTION_CONCENTRATION).contains(&alpha_alt),
        "the alternative concentration must be finite, non-negative and at most \
         {MAX_PROJECTION_CONCENTRATION:e} chromosomes, got {alpha_alt}"
    );

    out.fill(0.0);
    if alpha_alt == 0.0 {
        // No alternative allele exists to be drawn, so every site is monomorphic. Taken here
        // because the general sum would otherwise reach `lgamma(0)`.
        out[0] = 1.0;
        return;
    }

    let f = inbreeding.get();
    let concentration_total = alpha_ref + alpha_alt;
    let log_pair_constant = lgamma(concentration_total) - lgamma(alpha_alt) - lgamma(alpha_ref);
    // `ln k!` for every count this call can reach, filled once. The sum asks for binomial
    // coefficients with heavily repeated arguments, and reading them from here rather than calling
    // `lgamma` three times apiece is bit-identical and several times faster. This and the branch
    // splits below are allocated per *prediction*, of which a fit runs several hundred — about
    // 150 kB a time at 3,200 individuals. Spec §8's no-allocation rule governs the prior's row,
    // which the calling loop runs once per sample per pass; this runs once per run.
    let log_factorial: Vec<f64> = (0..=2 * n + 1).map(|k| lgamma(k as f64 + 1.0)).collect();
    let log_binomial = |top: usize, chosen: usize| {
        log_factorial[top] - log_factorial[chosen] - log_factorial[top - chosen]
    };

    // How likely each split of the panel into inbred and outbred individuals is, and the likeliest
    // of them — the tail is measured against that rather than against zero.
    let splits: Vec<Option<f64>> = (0..=n)
        .map(|identical_by_descent| log_branch_split(n, identical_by_descent, f))
        .collect();
    let tail_floor = splits
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        + BRANCH_TAIL_TOLERANCE.ln();

    for (identical_by_descent, split) in splits.iter().enumerate() {
        let Some(log_branch_weight) = *split else {
            continue;
        };
        if log_branch_weight < tail_floor || log_branch_weight < NEGLIGIBLE_BRANCH_WEIGHT.ln() {
            continue;
        }
        let branch_weight = log_branch_weight.exp();
        // The panel's distinct chromosomes: one for each identical-by-descent individual, two for
        // each of the others.
        let distinct = 2 * n - identical_by_descent;
        let single_chromosomes = distinct - identical_by_descent; // 2 × the outbred individuals
        let log_draw_constant = log_pair_constant - lgamma(concentration_total + distinct as f64);

        for alternative_draws in 0..=distinct {
            // How likely it is that exactly this many of the distinct chromosomes carry the
            // alternative allele — a beta-binomial weight, and a probability in its own right.
            let draw_weight = (log_binomial(distinct, alternative_draws)
                + lgamma(alpha_alt + alternative_draws as f64)
                + lgamma(alpha_ref + (distinct - alternative_draws) as f64)
                + log_draw_constant)
                .exp();
            if draw_weight == 0.0 {
                continue;
            }

            // How many of the duplicated chromosomes are among those draws. Each one contributes a
            // second copy, so it moves the class up by one. Given the draw count this is a
            // hypergeometric weight, walked out from its mode by an exact ratio — see the cost
            // section for why the mode and not the low end.
            let lowest = alternative_draws.saturating_sub(single_chromosomes);
            let highest = identical_by_descent.min(alternative_draws);
            let mode = (((alternative_draws + 1) * (identical_by_descent + 1)) / (distinct + 2))
                .clamp(lowest, highest);
            let at_mode = (log_binomial(identical_by_descent, mode)
                + log_binomial(single_chromosomes, alternative_draws - mode)
                - log_binomial(distinct, alternative_draws))
            .exp();
            let scale = branch_weight * draw_weight;
            out[alternative_draws + mode] += scale * at_mode;

            let step_up = |doubled: usize| {
                let rise =
                    ((identical_by_descent - doubled) * (alternative_draws - doubled)) as f64;
                let fall = ((doubled + 1)
                    * ((single_chromosomes + doubled + 1) - alternative_draws))
                    as f64;
                rise / fall
            };
            let mut climbing = at_mode;
            for doubled in mode..highest {
                climbing *= step_up(doubled);
                if climbing == 0.0 {
                    break;
                }
                out[alternative_draws + doubled + 1] += scale * climbing;
            }
            let mut falling = at_mode;
            for doubled in (lowest..mode).rev() {
                falling /= step_up(doubled);
                if falling == 0.0 {
                    break;
                }
                out[alternative_draws + doubled] += scale * falling;
            }
        }
    }
}

/// `ln` of the chance that exactly `identical_by_descent` of `individuals` took the
/// identical-by-descent branch — a binomial at the panel's `F`.
///
/// `None` where that chance is exactly zero, which is every count but one at `F = 0` and at
/// `F = 1`. Returning `None` rather than `−∞` is what keeps `0 × ln 0` — a `NaN` — out of the sum
/// at both ends of the range, where the caller spends most of its life: `F = 0` is an outbred
/// panel and the pre-pass's own default.
fn log_branch_split(
    individuals: usize,
    identical_by_descent: usize,
    inbreeding: f64,
) -> Option<f64> {
    if inbreeding == 0.0 {
        return (identical_by_descent == 0).then_some(0.0);
    }
    if inbreeding == 1.0 {
        return (identical_by_descent == individuals).then_some(0.0);
    }
    let outbred = individuals - identical_by_descent;
    Some(
        lgamma(individuals as f64 + 1.0)
            - lgamma(identical_by_descent as f64 + 1.0)
            - lgamma(outbred as f64 + 1.0)
            + identical_by_descent as f64 * inbreeding.ln()
            + outbred as f64 * (1.0 - inbreeding).ln(),
    )
}

/// The largest panel this projection will fit a pair to.
///
/// **A time bound, not a modelling one.** One prediction grows about as `N^2.45` and a fit runs
/// several hundred of them: 12.7 minutes at 3,200 individuals, so 10,000 would be about two
/// hours and 100,000 about a fortnight, with nothing on the way saying the run had stopped being
/// a run. Ten thousand diploid individuals is three times the several thousand this caller
/// commits to (`CLAUDE.md`, *What this caller has to work on*), so a panel that trips this is
/// asking for something the cost table does not cover, and should be told rather than left
/// waiting.
const MAX_PROJECTION_INDIVIDUALS: u32 = 10_000;

/// How far the pre-pass's class weights may sit from summing to one before they stop being a
/// spectrum. Wide enough for an accumulation over a few thousand classes, far too narrow to admit
/// counts that were never normalised.
const SPECTRUM_NORMALISATION_TOLERANCE: f64 = 1e-9;

/// The reference allele's concentration on a neutral panel — where the fit lands rather than a
/// number anyone chose, and the `1/p` density written as a Dirichlet
/// (`doc/devel/ng/spec/calling_priors.md` §4). It is what holds the heterozygote to
/// homozygous-alternative prior ratio near 2:1, so raising it is the §2.3 trap.
const NEUTRAL_ALPHA_REF: f64 = 1.0;

/// The lowest and highest **total** concentration the fit will consider, in chromosomes —
/// `α_ref + α_alt`, which is how much conviction the prior carries about the frequency.
///
/// The neutral panel's answer is just above 1, and the range is three decades either side.
/// Bounds and not merely starting points, because the search line-searches the whole range on
/// each axis.
const CONCENTRATION_TOTAL_SEARCH_RANGE: (f64, f64) = (1e-3, 1e3);

/// The lowest and highest **ratio** `α_alt / α_ref` the fit will consider — the odds the prior
/// gives the alternative allele, which for small values is the expected frequency itself.
///
/// The bottom is where polymorphism stops being visible to any panel: at a ratio of `1e-9` the
/// share of segregating sites is about `α_alt · H(2N)`, roughly 8 sites in a thousand million
/// over a thousand chromosomes, so nothing below it can be told from a fully invariant cohort.
/// The top is four decades above the most diverse panel anyone would call — a `θ` of 1 in 100 is
/// already ten times human diversity.
const CONCENTRATION_RATIO_SEARCH_RANGE: (f64, f64) = (1e-9, 1e2);

/// Where each start begins, as `(total concentration, α_alt / α_ref)`.
///
/// **Every start differs from every other on both axes**, which is not decoration: the sibling
/// path's inbreeding fit returned a confident zero from five starts that disagreed about the
/// headline number while sharing one guess at a nuisance axis (`fitting/multistart.rs`). The four
/// cover 3.3 of the total's 6 decades and 7 of the ratio's 11, rather than clustering near the
/// neutral pair, so a search that only ever returns its own starting neighbourhood is visible as
/// disagreement.
const SEARCH_STARTS: [(f64, f64); 4] = [(0.02, 1e-7), (0.3, 1e-5), (3.0, 1e-2), (40.0, 1.0)];

/// The three directions each sweep line-searches along, in the search's own log coordinates —
/// **one for each of the three quantities the two parametrisations name between them.**
///
/// Writing `t = ln(α_ref + α_alt)` and `r = ln(α_alt / α_ref)`, the identities are
/// `ln α_ref = t − ln(1 + e^r)` and `ln α_alt = t + r − ln(1 + e^r)`. While the ratio is small
/// the last term is negligible, so:
///
/// ```text
///   [1,  0]        the total concentration, α_ref and α_alt moving together
///   [0,  1]        α_alt alone — α_ref is unchanged
///   [√½, −√½]      α_ref alone — α_alt is unchanged
/// ```
///
/// So this sweeps the axes of *both* parametrisations, which is what [`fit_pair`] needs and what
/// neither pair of coordinates gives on its own.
///
/// **The fourth direction a rotation suggests, `[√½, √½]`, was tried and removed.** It is not a
/// coordinate of either parametrisation — along it `ln α_alt` moves twice as fast as
/// `ln α_ref` — and it does harm rather than nothing: on a spectrum with all its mass at
/// intermediate frequency, which is the shape spec §4.1 says two parameters cannot hold, the
/// sweep including it ended at `α_ref = 3.59` where the best point in the box is 498, a
/// log-likelihood of −3.410 against −2.232. It also cost up to 314 predictions a fit.
const SEARCH_DIRECTIONS: [[f64; 2]; 3] = [
    [1.0, 0.0],
    [0.0, 1.0],
    [
        std::f64::consts::FRAC_1_SQRT_2,
        -std::f64::consts::FRAC_1_SQRT_2,
    ],
];

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

/// The pre-pass's fitted spectrum, in a wrapper whose invariants a struct literal cannot skip.
///
/// **The nesting is load-bearing, and it is the same device `mod.rs` uses**: a private field is
/// visible to a module's *descendants*, so a type declared directly in `seed_generic` could be
/// built field by field — skipping both checks — from anywhere in this file, its test modules
/// included. One level of nesting makes those siblings instead, and the literal fails with
/// `error[E0451]`.
mod checked {
    use super::{MAX_PROJECTION_INDIVIDUALS, SPECTRUM_NORMALISATION_TOLERANCE};

    /// **What the pre-pass measured about how allele frequencies are spread across the panel**,
    /// held as a borrow of whatever the cohort gather emits.
    ///
    /// One weight per allele-count class — what share of sites carry the alternative allele on no
    /// chromosome, on exactly one, on exactly two, and so on to every chromosome of the panel.
    /// `2N + 1` of them for `N` diploid individuals, summing to one
    /// (`doc/devel/ng/spec/parameter_prepass_cohort.md` §4.1).
    ///
    /// **A view rather than a type of its own**, because the pre-pass's cohort gather is not
    /// built yet: when it lands it will own the weights, and this borrows them.
    ///
    /// **It is not
    /// [`FrequencyDensity`](crate::ng::parameter_estimation::joint::fit::FrequencyDensity), which
    /// is a different object.** That is four numbers describing how the *population's* allele
    /// frequency is distributed — two point masses and a Beta over what segregates. This is a
    /// spread over a *panel's* allele counts. The projection matches class weights, so class
    /// weights are what it takes.
    ///
    /// The two numbers beside the weights are how the run says which information produced its
    /// seed: how many sites' worth of pseudo-counts held the estimate at the neutral shape, and
    /// how many census sites actually came out variable. Two runs that used different information
    /// are otherwise indistinguishable in what they emit
    /// (`doc/devel/ng/spec/calling_priors.md` §4.1).
    #[derive(Copy, Clone, Debug)]
    pub struct FittedSpectrum<'a> {
        class_weights: &'a [f64],
        regularizer_site_weight: f64,
        variable_census_sites: f64,
    }

    impl<'a> FittedSpectrum<'a> {
        /// Wrap the pre-pass's weights, refusing a shape that is not a spectrum.
        ///
        /// **The class count fixes the panel size**, `2N + 1` for `N` diploid individuals, so an
        /// even count is refused rather than rounded: nothing else in this module knows how many
        /// individuals the spectrum was fitted at, and a panel size taken from a second argument
        /// is a second place for it to be wrong.
        ///
        /// The weights are checked to be a distribution — non-negative, summing to one within
        /// [`SPECTRUM_NORMALISATION_TOLERANCE`] — because the projection's objective is a
        /// log-likelihood against them, and weights that are not a distribution score every
        /// candidate wrongly by an amount that varies with the candidate.
        pub fn new(
            class_weights: &'a [f64],
            regularizer_site_weight: f64,
            variable_census_sites: f64,
        ) -> Self {
            assert!(
                class_weights.len() >= 3
                    && class_weights.len() % 2 == 1
                    && class_weights.len() <= 2 * MAX_PROJECTION_INDIVIDUALS as usize + 1,
                "a panel of N diploid individuals has 2N + 1 allele-count classes, so the count \
                 is odd, at least 3, and at most {} for the {MAX_PROJECTION_INDIVIDUALS} \
                 individuals this projection will fit; got {}",
                2 * MAX_PROJECTION_INDIVIDUALS as usize + 1,
                class_weights.len()
            );
            // The two halves are checked apart because they fail differently and a caller who
            // trips one learns nothing from the other's number: weights of `[-0.1, 0.6, 0.5]`
            // total exactly 1.
            if let Some((class, weight)) = class_weights
                .iter()
                .enumerate()
                .find(|(_, weight)| weight.is_nan() || **weight < 0.0)
            {
                panic!(
                    "a class weight is the share of sites in that allele-count class, so it \
                     cannot be negative or NaN: class {class} holds {weight}"
                );
            }
            let total: f64 = class_weights.iter().sum();
            assert!(
                (total - 1.0).abs() <= SPECTRUM_NORMALISATION_TOLERANCE,
                "the class weights must sum to 1 within {SPECTRUM_NORMALISATION_TOLERANCE:e} — \
                 the projection's objective is a log-likelihood against them, so raw counts score \
                 every candidate wrongly by an amount that varies with the candidate; got a total \
                 of {total}"
            );
            assert!(
                regularizer_site_weight.is_finite() && regularizer_site_weight >= 0.0,
                "how many sites' worth of pseudo-counts held the spectrum at the neutral shape \
                 must be finite and non-negative, got {regularizer_site_weight}"
            );
            assert!(
                variable_census_sites.is_finite() && variable_census_sites >= 0.0,
                "the count of variable census sites is non-negative, got {variable_census_sites}"
            );
            Self {
                class_weights,
                regularizer_site_weight,
                variable_census_sites,
            }
        }

        /// The share of sites in each allele-count class, class 0 first.
        #[inline]
        pub fn class_weights(&self) -> &'a [f64] {
            self.class_weights
        }

        /// How many sites' worth of pseudo-counts held the estimate at the neutral shape.
        #[inline]
        pub fn regularizer_site_weight(&self) -> f64 {
            self.regularizer_site_weight
        }

        /// How many census sites came out variable across the panel.
        #[inline]
        pub fn variable_census_sites(&self) -> f64 {
            self.variable_census_sites
        }

        /// How many diploid individuals the spectrum was fitted at, read off the class count.
        #[inline]
        pub fn individuals(&self) -> u32 {
            // The constructor's ceiling puts this far inside a `u32`.
            ((self.class_weights.len() - 1) / 2) as u32
        }
    }
}

pub use checked::FittedSpectrum;

/// **Read the run's two starting numbers off the pre-pass's fitted spectrum**: the reference
/// allele's concentration and the total shared out across whatever alternative alleles a locus
/// turns out to carry.
///
/// ## What it does
///
/// The neutral `1/p` density and the neutral frequency spectrum are the same statement written
/// twice — once at a locus, once across a panel. So the pre-pass's spectrum and this module's
/// concentration pair describe one object at two sample sizes, and getting the pair from the
/// spectrum is **a change of representation rather than a second estimate**: nothing is fitted
/// here that the pre-pass has not already fitted
/// (`doc/devel/ng/spec/calling_priors.md` §4.1).
///
/// The fit takes the pair whose predicted spectrum is most likely to have produced the fitted
/// one — equivalently, the pair closest to it in Kullback–Leibler divergence — **over every
/// allele-count class including the monomorphic one**. The monomorphic class is not a
/// throwaway: it is what pins the alternative concentration against the diversity rather than
/// leaving only the shape of the fall-off identified.
///
/// **The prediction carries the panel's inbreeding** ([`fill_expected_spectrum`]), because a
/// panel's chromosomes are not independent draws once its individuals are inbred. Treating them
/// as independent biases the reference concentration down with a fixed sign — 8.6% at `F = 0.6`,
/// 12.1% at 0.8 and 14.0% at 0.9, on a panel of 26 individuals at tomato's diversity, measured in
/// [`projection_tests::the_projection_returns_one_pair_at_every_inbreeding_coefficient`].
///
/// ## The three regimes, and none of them is a branch on cohort size
///
/// - **A spectrum arrived** — fit it, and report the regularizer beside the answer. **A spectrum
///   makes the diversity moot**: it carries its own scale, so `diversity` is not read and the
///   fallback cannot be reported here.
/// - **No spectrum** — the pair is the neutral `(1, θ)` at the diversity the pre-pass *did* fit.
///   The pre-pass emits the spectrum as absent below a panel-size floor rather than as a thin
///   estimate, so a cohort of five arrives here without one while a single sample arrives with
///   one. A spectrum too thin to emit and a panel with nothing to fit carry the same information
///   about shape.
/// - **No spectrum and no diversity either** — the same neutral pair at the species-range
///   fallback, and the run must say so.
///
/// **Nothing here tests how many samples the cohort holds.** The one branch is on what the
/// pre-pass had.
///
/// ## It takes no variant class, and that is a decision rather than an omission
///
/// The design keeps the door open to giving substitutions and short insertions or deletions
/// different diversities (spec §4.2, Q1). **The end that will apply that split is
/// [`fill_locus_concentration`], not this one** — settled by the owner, 2026-08-22. This function
/// reads the *shape* of variation off the panel's own allele counts, which the pre-pass fits
/// without separating the two classes; a class-specific *scale* belongs where the run's total is
/// shared out over a locus's alleles, and that is the only end able to describe a site carrying
/// one alternative of each kind. Carrying the argument at both ends would apply the ratio twice.
///
/// ## Cost
///
/// One prediction is the expensive step and a fit runs **399 of them, once per run** — the same
/// count at every panel size and every inbreeding coefficient measured, because what the search
/// costs is set by how finely it resolves each direction and not by the panel. The count is
/// asserted in [`projection_tests::a_fit_costs_at_most_450_predictions`]; the wall clock is 3.8 s
/// at 400 individuals, 22 s at 800, 2.2 minutes at 1,600 and **11.8 minutes at 3,200**, measured
/// in [`projection_tests::the_cost_of_one_fit_by_panel_size`].
///
/// The search resolves each concentration to about 1% of itself, which is
/// [`SearchPrecision::fast`]'s reasoning applied here: a concentration shifted by 1% moves a call
/// far less than one more read would, so resolving it further would be spending time on digits no
/// genotype registers.
pub fn project_spectrum_seed(
    spectrum: Option<FittedSpectrum<'_>>,
    diversity: Option<ExpectedHeterozygosity>,
    panel_inbreeding: InbreedingF,
) -> SpectrumSeed {
    let Some(spectrum) = spectrum else {
        return match diversity {
            Some(theta) => {
                SpectrumSeed::new(NEUTRAL_ALPHA_REF, theta.get(), SeedRegime::NeutralShape)
            }
            None => SpectrumSeed::new(
                NEUTRAL_ALPHA_REF,
                ExpectedHeterozygosity::SPECIES_FALLBACK.get(),
                SeedRegime::FallbackDiversity,
            ),
        };
    };

    let fit = fit_pair(&spectrum, panel_inbreeding, SearchPrecision::fast());
    SpectrumSeed::new(
        fit.alpha_ref,
        fit.alpha_alt,
        SeedRegime::FittedSpectrum {
            regularizer_site_weight: spectrum.regularizer_site_weight(),
            // **This is the panel-wide comparison, and spec §4.1 is explicit that a panel-wide
            // ratio is the wrong number to quote as reassurance** — on the panel it measures, the
            // aggregate was 3,100 to 1 while the thinnest class held two sites and was outweighed
            // only 39 to 1. The per-class ratio is the pre-pass's to emit beside its spectrum
            // (arch §4); this module carries the aggregate through and claims nothing more.
            census_sites_outweigh_regularizer: spectrum.variable_census_sites()
                > spectrum.regularizer_site_weight(),
            spectrum_match: fit.spectrum_match,
        },
    )
}

/// What the fit returned: the pair, what it scored, and what it cost.
///
/// **Whether the search settled or ran out of sweeps is deliberately not here.** This surface has
/// no concavity proof, unlike the climb the sibling driver wraps, so running out is a data
/// condition rather than a defect, and what protects the answer is that the best point evaluated
/// is the one returned ([`sweep_from`]).
///
/// **The last two are read only by this module's tests**, which is what the attribute says: the
/// score is how [`projection_tests::the_winning_score_is_the_spectrums_own_entropy`] checks the
/// fit found the maximum, and the count is what
/// [`projection_tests::a_fit_costs_at_most_600_predictions`] holds. Carried on the fit rather
/// than recomputed, because neither is derivable from the pair afterwards.
#[derive(Copy, Clone, Debug)]
#[cfg_attr(not(test), allow(dead_code))]
struct ProjectionFit {
    alpha_ref: f64,
    alpha_alt: f64,
    /// What the winning pair scored under [`spectrum_log_likelihood`]: the negative
    /// Kullback–Leibler divergence from the fitted spectrum to the predicted one, **plus** the
    /// fitted spectrum's own entropy — a constant that does not move with the candidate pair, so
    /// the maximum of this is the minimum of that divergence. On a spectrum the two-parameter
    /// family can reproduce exactly, the winning score is therefore that entropy.
    log_likelihood: f64,
    /// How many predictions the whole fit cost — 399 everywhere measured — and **half of what a
    /// fit costs rather than the whole of it**. A prediction *inside a fit* averages 1.78 s at
    /// 3,200 individuals, against the 0.96 s
    /// `doc/devel/ng/reports/spectrum_projection_cost_2026-08-22.md` measures **at the neutral
    /// pair**, because the search spends most of its predictions away from that pair where the
    /// branch-tail trim drops fewer splits. So one fit there is 11.8 minutes measured, against
    /// the 6.4 the count and that report alone would predict.
    predictions: usize,
    /// Whether the winning pair reproduces the spectrum or is the closest the family reaches.
    spectrum_match: SpectrumMatch,
}

/// Fit the concentration pair to a spectrum: line searches along four directions from several
/// starts, on the log scale of the total concentration and of the ratio between the two.
///
/// **On log scales because a concentration spans decades** — the same argument
/// [`SearchableNoise`](crate::ng::parameter_estimation::fitting::multistart::SearchableNoise)
/// makes for a slippage level: a step is then a fixed *fraction* of the value, and the search
/// cannot walk out of the positive range whatever step it takes.
///
/// **Three directions and not two, because the surface is a ridge and which way it runs depends
/// on the panel size.** Searching `ln α_ref` and `ln α_alt` separately finds the answer at 26 and
/// 63 individuals and fails at one, returning about 0.24 to 0.84 where the answer is 1 depending
/// on the box searched, with the starts spread thousands-fold. Searching the total and the ratio
/// separately does the opposite: 7 parts in ten million high at one individual, 1.3% high at 63.
/// The two parametrisations are a 45° rotation of each other in log space, and between them they
/// name three quantities — the total, `α_ref` alone and `α_alt` alone — which is what
/// [`SEARCH_DIRECTIONS`] sweeps. Neither panel size is then a special case, and **nothing here
/// tests how many individuals there are**; the third direction costs the same at every size.
///
/// **This is not `fit_by_multistart`, and the reason is cost rather than taste.** That driver
/// scores one *cell* at a time through
/// [`NoiseModel::append_genotype_likelihoods`](crate::ng::parameter_estimation::fitting::NoiseModel::append_genotype_likelihoods),
/// which takes `&self` and so cannot cache; the natural cell here is one allele-count class, and
/// every class would rebuild the whole spectrum. At 3,200 individuals that is 6,401 predictions
/// where one is needed — about 1.7 hours per candidate against 0.96 seconds. It also has no
/// notion of a diagonal direction, which the paragraph above needs. The search's *shape* is that
/// driver's, and so is [`SearchPrecision`]; the driver itself does not fit.
fn fit_pair(
    spectrum: &FittedSpectrum<'_>,
    panel_inbreeding: InbreedingF,
    precision: SearchPrecision,
) -> ProjectionFit {
    let mut scorer = SpectrumScorer::new(spectrum, panel_inbreeding);

    // **Every start gets one sweep, and only the best-scoring one is then swept to convergence.**
    // A sweep line-searches all four directions over the whole box, so one of them already looks
    // everywhere along each direction; what the other starts are for is a *second* optimum
    // elsewhere on the surface, and one sweep each is what would find it. Sweeping all four to
    // convergence costs about 2.7 times as many predictions and, on every spectrum measured in
    // this module's tests, reaches the same pair to within 1%.
    let mut best = ScoredPoint {
        point: [SEARCH_STARTS[0].0.ln(), SEARCH_STARTS[0].1.ln()],
        log_likelihood: f64::NEG_INFINITY,
    };
    for (from_total, from_ratio) in SEARCH_STARTS {
        let start = ScoredPoint {
            point: [from_total.ln(), from_ratio.ln()],
            log_likelihood: f64::NEG_INFINITY,
        };
        let (swept, _) = sweep_once(start, precision, &mut |point| scorer.score(point));
        best = best.max_by_score(swept);
    }

    let best = sweep_from(best, precision, &mut |point| scorer.score(point));

    let (alpha_ref, alpha_alt) = concentrations_at(best.point);
    let spectrum_match = scorer.match_at(best.point, precision.tolerance);
    ProjectionFit {
        alpha_ref,
        alpha_alt,
        log_likelihood: best.log_likelihood,
        predictions: scorer.predictions,
        spectrum_match,
    }
}

/// The objective, with the buffer it predicts into and the count of how often it has been asked.
/// One per fit.
///
/// **One value rather than four arguments threaded through every closure**, because all four are
/// the same thing — state that outlives a single candidate — and the alternative is the same
/// six-argument call written out at each of the search's two entry points.
struct SpectrumScorer<'a> {
    spectrum: &'a FittedSpectrum<'a>,
    panel_inbreeding: InbreedingF,
    /// Refilled once per candidate rather than allocated there: a fit runs several hundred
    /// candidates and this is the only thing it builds.
    predicted: Vec<f64>,
    predictions: usize,
}

impl<'a> SpectrumScorer<'a> {
    fn new(spectrum: &'a FittedSpectrum<'a>, panel_inbreeding: InbreedingF) -> Self {
        Self {
            spectrum,
            panel_inbreeding,
            predicted: vec![0.0; spectrum.class_weights().len()],
            predictions: 0,
        }
    }

    /// **Whether the winning pair reproduces the measured spectrum**, run once when the search is
    /// over rather than at every candidate — one more prediction against the fit's 399.
    ///
    /// Two things make a pair a compromise rather than an answer, checked in that order because
    /// the first is the more informative:
    ///
    /// - **it predicts effectively nothing for a class the measurement gives real weight to**, so
    ///   no pair in the family can have produced this spectrum. The comparison is against
    ///   [`PROBABILITY_FLOOR`], which is what the objective charged the pair for that class;
    /// - **it sits on the edge of the range searched**, so a better pair may lie outside it.
    ///   "On the edge" means within the resolution the search was asked for, because a line
    ///   search stops at its bracket's midpoint and cannot land exactly on a bound.
    fn match_at(&mut self, point: [f64; 2], tolerance: f64) -> SpectrumMatch {
        let (alpha_ref, alpha_alt) = concentrations_at(point);
        fill_expected_spectrum(
            alpha_ref,
            alpha_alt,
            self.spectrum.individuals(),
            self.panel_inbreeding,
            &mut self.predicted,
        );
        self.predictions += 1;
        let cannot_produce = self
            .spectrum
            .class_weights()
            .iter()
            .zip(&self.predicted)
            .any(|(weight, prediction)| *weight > 0.0 && *prediction < PROBABILITY_FLOOR);
        if cannot_produce {
            return SpectrumMatch::Unreproducible;
        }
        let at_limit = [
            (
                point[0],
                CONCENTRATION_TOTAL_SEARCH_RANGE.0.ln(),
                CONCENTRATION_TOTAL_SEARCH_RANGE.1.ln(),
            ),
            (
                point[1],
                CONCENTRATION_RATIO_SEARCH_RANGE.0.ln(),
                CONCENTRATION_RATIO_SEARCH_RANGE.1.ln(),
            ),
        ]
        .iter()
        .any(|(at, low, high)| (at - low).abs() <= tolerance || (high - at).abs() <= tolerance);
        if at_limit {
            SpectrumMatch::AtSearchLimit
        } else {
            SpectrumMatch::Reproduced
        }
    }

    /// One candidate pair's score: predict the spectrum it implies and read the objective off it.
    fn score(&mut self, point: [f64; 2]) -> f64 {
        self.predictions += 1;
        let (alpha_ref, alpha_alt) = concentrations_at(point);
        fill_expected_spectrum(
            alpha_ref,
            alpha_alt,
            self.spectrum.individuals(),
            self.panel_inbreeding,
            &mut self.predicted,
        );
        spectrum_log_likelihood(self.spectrum.class_weights(), &self.predicted)
    }
}

/// A point and what it scored, so the two never come apart on the way back.
#[derive(Copy, Clone, Debug)]
struct ScoredPoint {
    point: [f64; 2],
    log_likelihood: f64,
}

impl ScoredPoint {
    /// Whichever of the two scored higher, this one on a tie.
    fn max_by_score(self, other: ScoredPoint) -> ScoredPoint {
        if self.log_likelihood >= other.log_likelihood {
            self
        } else {
            other
        }
    }
}

/// Sweep the four directions from one point until nothing moves further than the resolution asked
/// for, or the sweeps run out. Returns the best point any sweep reached, not merely the last.
///
/// **Running out is not asserted against.** This surface has no concavity proof, unlike the climb
/// the sibling driver wraps, so a capped search is a data condition; what protects the answer is
/// that the best point evaluated is the one returned.
fn sweep_from(
    from: ScoredPoint,
    precision: SearchPrecision,
    score: &mut impl FnMut([f64; 2]) -> f64,
) -> ScoredPoint {
    let mut best = from;
    for _ in 0..precision.max_sweeps {
        let (swept, furthest_move) = sweep_once(best, precision, score);
        best = swept;
        if furthest_move < precision.tolerance {
            break;
        }
    }
    best
}

/// One sweep: a line search along each of the three directions in turn. Returns where it ended
/// and what that scored, and how far the furthest of the three searches moved.
///
/// **A sweep cannot end below where it began**, because no line search in it can — see
/// [`line_search`].
fn sweep_once(
    from: ScoredPoint,
    precision: SearchPrecision,
    score: &mut impl FnMut([f64; 2]) -> f64,
) -> (ScoredPoint, f64) {
    let mut current = from;
    let mut furthest_move: f64 = 0.0;
    for direction in SEARCH_DIRECTIONS {
        let (reached, moved) = line_search(current, direction, precision, score);
        current = reached;
        furthest_move = furthest_move.max(moved);
    }
    (current, furthest_move)
}

/// Golden-section along one direction, over the whole part of the search box that lies on that
/// line. Returns where it ended and what that scored, and how far it moved — **the same order as
/// [`sweep_once`]**, because two `f64`s that mean different things are what the next author
/// copies the wrong way round.
///
/// **It never returns a point worse than the one it started from.** The bracket's midpoint is
/// where the search stopped, not necessarily the best thing it saw: on a line that is not
/// unimodal, or simply at the resolution's own edge, it can sit below the start. Measured on the
/// exact neutral spectrum at 26 individuals, 5 of 80 line searches ended below their start (worst
/// 6.3e-9 nats) and on a flat spectrum 31 of 80 (worst 1.1e-6). Keeping the better of the two
/// costs no prediction, because both scores are already in hand.
///
/// **Golden section keeps one of its two interior evaluations at every step**, so a bracket costs
/// one prediction per step rather than two.
fn line_search(
    from: ScoredPoint,
    direction: [f64; 2],
    precision: SearchPrecision,
    score: &mut impl FnMut([f64; 2]) -> f64,
) -> (ScoredPoint, f64) {
    let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
    let (mut low, mut high) = bounds_along(from.point, direction);
    let at = |t: f64| {
        [
            from.point[0] + t * direction[0],
            from.point[1] + t * direction[1],
        ]
    };

    let mut left = high - inverse_phi * (high - low);
    let mut right = low + inverse_phi * (high - low);
    let mut at_left = score(at(left));
    let mut at_right = score(at(right));
    for _ in 0..precision.max_axis_steps {
        if at_left > at_right {
            high = right;
            right = left;
            at_right = at_left;
            left = high - inverse_phi * (high - low);
            at_left = score(at(left));
        } else {
            low = left;
            left = right;
            at_left = at_right;
            right = low + inverse_phi * (high - low);
            at_right = score(at(right));
        }
        if (high - low) < precision.tolerance {
            break;
        }
    }
    let stopped_at = 0.5 * (low + high);
    let point = at(stopped_at);
    let reached = ScoredPoint {
        point,
        log_likelihood: score(point),
    };
    if reached.log_likelihood > from.log_likelihood {
        (reached, stopped_at.abs())
    } else {
        (from, 0.0)
    }
}

/// How far a point may travel along a direction before it leaves the search box, either way.
///
/// **The box is the two ranges above**, and a diagonal leaves it through whichever side comes
/// first, so both ends are the tightest of the per-axis limits.
fn bounds_along(from: [f64; 2], direction: [f64; 2]) -> (f64, f64) {
    let box_bounds = [
        (
            CONCENTRATION_TOTAL_SEARCH_RANGE.0.ln(),
            CONCENTRATION_TOTAL_SEARCH_RANGE.1.ln(),
        ),
        (
            CONCENTRATION_RATIO_SEARCH_RANGE.0.ln(),
            CONCENTRATION_RATIO_SEARCH_RANGE.1.ln(),
        ),
    ];
    let mut low = f64::NEG_INFINITY;
    let mut high = f64::INFINITY;
    for axis in 0..2 {
        if direction[axis] == 0.0 {
            continue;
        }
        let to_low = (box_bounds[axis].0 - from[axis]) / direction[axis];
        let to_high = (box_bounds[axis].1 - from[axis]) / direction[axis];
        low = low.max(to_low.min(to_high));
        high = high.min(to_low.max(to_high));
    }
    (low, high)
}

/// The two concentrations a search point stands for, `(α_ref, α_alt)`. A point is
/// `[ln(α_ref + α_alt), ln(α_alt / α_ref)]` — the log of the total and the log of the ratio.
///
/// **The search runs on the total and the ratio rather than on the two concentrations directly,
/// and at one individual that is the difference between finding the answer and not.** In
/// `(ln α_ref, ln α_alt)` both axes move both of the things the spectrum is actually sensitive
/// to — how often the alternative allele is expected, and how tightly that is believed — so the
/// surface is a long curved ridge and a coordinate search walks along it a step at a time.
/// Measured at one individual and tomato's diversity, the search in those coordinates returned
/// `α_ref = 0.844` where the answer is 1, with the four starts spread 3,206-fold; in these it
/// returns 1 to nine decimal places with the starts agreeing to 1.000000. At 26 individuals the
/// old coordinates worked, so the failure is specific to the panel size where the spectrum has
/// three classes and almost all the weight is in one of them.
#[inline]
fn concentrations_at(point: [f64; 2]) -> (f64, f64) {
    let total = point[0].exp();
    // `α_alt / (α_ref + α_alt)` from the log odds. The logistic rather than `r / (1 + r)` because
    // it stays exact where the ratio is enormous: at `ln r = 745` the naive form's `1 + r`
    // overflows and this returns a share of exactly 1, so `α_ref` would come out at 0.0 and
    // `fill_expected_spectrum` would refuse it. Unreachable while [`bounds_along`] holds — the
    // box stops at `ln r = ln 1e2` — and it is the sort of thing a widened bound would find.
    let alternative_share = 1.0 / (1.0 + (-point[1]).exp());
    (total * (1.0 - alternative_share), total * alternative_share)
}

/// How likely the fitted spectrum's class weights are under a predicted spectrum:
/// `Σ_k w_k · ln q_k`.
///
/// **This is the maximum-likelihood objective spec §4.1 names, and it is the negative
/// Kullback–Leibler divergence up to a constant** — the fitted spectrum's own entropy, which does
/// not move with the candidate pair, so maximising one minimises the other.
///
/// **A class the fitted spectrum gives no weight contributes nothing**, and is skipped rather
/// than multiplied: `0 × ln 0` is a `NaN`, and predicted zeros are ordinary here — at full
/// inbreeding every odd class is exactly zero.
///
/// **A class the fitted spectrum does give weight to and the candidate predicts at zero is
/// floored rather than allowed to reach `−∞`** ([`PROBABILITY_FLOOR`], the rule spec §8 sets for
/// every logarithm in this module). Two reasons, and neither is tidiness:
///
/// - **`−∞` is not an ordering.** A candidate that predicts `1e-320` for an occupied class and
///   one that predicts exactly zero are equally impossible, but the search has to prefer the
///   first to walk towards anything. Golden section compares `−∞ > −∞` as false and then walks to
///   one end of the line blind. That region is not hypothetical: on the exact neutral spectrum at
///   `F = 0.8`, none of 441 points across the search box scores `−∞` at 150 individuals, 17 do at
///   400, and 28 of 225 do at 1,600 — one point in eight in the middle of the committed cohort
///   range.
/// - **A spectrum no candidate can produce would otherwise have no answer at all.** At `F = 1`
///   every odd class is exactly zero, so a spectrum carrying any weight in an odd class scores
///   `−∞` at *every* pair — measured, 441 of 441 grid points — and the fit would return whichever
///   point it happened to start from. `InbreedingF` still admits `1.0` today; the `[0, 1)`
///   tightening spec §7 requires is another plan's (`calling_prerequisites.md` Milestone A).
///
/// **What the floor does not do is say so.** The pair that comes back is the closest the
/// two-parameter family can reach, and nothing on [`SpectrumSeed`] distinguishes that from a fit
/// that matched. That is the same complaint spec §12 test 11 makes about the STR seed's
/// unreachable diversity, and it is recorded as open at the end of this step.
fn spectrum_log_likelihood(class_weights: &[f64], predicted: &[f64]) -> f64 {
    assert_eq!(
        class_weights.len(),
        predicted.len(),
        "the objective runs over every allele-count class, and `zip` on a short prediction would \
         silently drop the top classes rather than fail"
    );
    let mut total = 0.0;
    for (weight, prediction) in class_weights.iter().zip(predicted) {
        if *weight > 0.0 {
            total += weight * prediction.max(PROBABILITY_FLOOR).ln();
        }
    }
    total
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
/// `α_ref` is floored by the same constant, and **no fit can make that bind**: the search box
/// bottoms out near `α_ref = 1e-5` — a total of at least `1e-3` against a ratio of at most `1e2`
/// ([`CONCENTRATION_TOTAL_SEARCH_RANGE`], [`CONCENTRATION_RATIO_SEARCH_RANGE`]) — ten million
/// times the floor. It is there because [`SpectrumSeed::new`] admits any strictly positive
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
/// 2026-08-22, and [`project_spectrum_seed`] no longer takes the argument because of it. That
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
pub fn fill_locus_concentration<'a>(
    seed: SpectrumSeed,
    class: VariantClass,
    allele_count: usize,
    out: &'a mut [f64],
) -> Concentration<'a> {
    // Both classes take the same seed today; the split is `project_spectrum_seed`'s to absorb
    // (spec §4.2, Q1).
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
        // A monomorphic locus takes neither branch: `out[1..]` is empty there, so this fill would
        // write nothing whatever the division produced. The condition is here to name the case.
        let per_alternative_allele =
            (seed.alpha_alt_total() / alternative_allele_count as f64).max(MIN_ALT_CONCENTRATION);
        out[1..].fill(per_alternative_allele);
    }
    Concentration::new(out)
}

/// The exact expected spectrum of a concentration pair — **every projection target in this
/// module's tests is built by this, and never by writing `θ/k`, and never by drawing sites**.
///
/// `θ/k` is the small-diversity approximation, and its own error is 0.272% at tomato's diversity
/// over 52 chromosomes and 4.4% at a `θ` of 1 in 100
/// ([`tests::the_neutral_shape_appears_in_the_small_diversity_limit`]) — larger than the effect
/// these tests measure, so a wiring bug and the approximation would be indistinguishable. Drawing
/// sites fails differently: Monte-Carlo noise falls as one over the square root of the site count,
/// so a floating-point tolerance would need more sites than anyone can generate
/// (`doc/devel/ng/spec/calling_priors.md` §12 test 5).
///
/// The buffer starts as `NaN` so a class the sum never writes shows up as one rather than as a
/// plausible zero.
#[cfg(test)]
fn exact_spectrum(alpha_ref: f64, alpha_alt: f64, individuals: u32, inbreeding: f64) -> Vec<f64> {
    let mut out = vec![f64::NAN; 2 * individuals as usize + 1];
    fill_expected_spectrum(
        alpha_ref,
        alpha_alt,
        individuals,
        InbreedingF::try_new(inbreeding).unwrap(),
        &mut out,
    );
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ln C(n, k)` for the tests' own oracles, written from `lgamma` so it shares no table with
    /// the function under test.
    fn log_binomial(n: usize, k: usize) -> f64 {
        lgamma(n as f64 + 1.0) - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0)
    }

    /// The file's one builder of an exact expected spectrum, named for what these tests read it
    /// as. Why it is never `θ/k` and never sampled is on [`exact_spectrum`].
    use super::exact_spectrum as spectrum;

    /// **Every site falls in exactly one class, so the classes sum to one** — across panel sizes,
    /// diversities and inbreeding coefficients.
    ///
    /// The property the whole construction rests on. It is not free: the sum runs over three
    /// nested indices with binomial coefficients that reach `10^36.8` at the largest panel here,
    /// and any bound that is off by one, or any branch weight that does not itself sum to one,
    /// shows up as a total away from 1.
    #[test]
    fn the_classes_of_a_panel_sum_to_one() {
        let mut worst = 0.0_f64;
        for individuals in [1_u32, 2, 3, 7, 26, 63] {
            for alpha_alt in [1e-4, 6e-4, 1e-2, 0.5, 2.0] {
                for alpha_ref in [0.5, 1.0, 3.0] {
                    for inbreeding in [0.0, 0.25, 0.6, 0.9, 1.0] {
                        let out = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                        let total: f64 = out.iter().sum();
                        worst = worst.max((total - 1.0).abs());
                        assert!(
                            (total - 1.0).abs() < 1e-9,
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, F \
                             {inbreeding}: the classes carry {total} rather than one unit of \
                             probability"
                        );
                        assert!(
                            out.iter().all(|c| c.is_finite() && *c >= 0.0),
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, F \
                             {inbreeding}: a class came out negative or not finite: {out:?}"
                        );
                    }
                }
            }
        }
        // Recorded so the `1e-9` above is known to be a bound rather than a hope: the worst
        // departure over this whole grid is 1.7e-13, about 6,000 times inside it. Asserted, so
        // that a change which quietly costs accuracy fails here rather than being absorbed by the
        // looser budget.
        assert!(worst < 1e-12, "worst departure from one was {worst}");
    }

    /// **At no inbreeding the panel is `2N` independent chromosomes, and the spectrum is exactly
    /// the beta-binomial** — the closed form written a second way, sharing no code with the sum
    /// under test.
    ///
    /// This is the oracle for the whole construction at one end of the `F` range: the three nested
    /// sums collapse to a single term per class, and that term is a beta-binomial probability,
    /// which is written here from `lgamma` directly.
    #[test]
    fn at_no_inbreeding_the_spectrum_is_the_beta_binomial() {
        for individuals in [1_u32, 3, 12, 26] {
            for alpha_alt in [1e-3, 1e-2, 0.4] {
                for alpha_ref in [0.5, 1.0, 2.5] {
                    let out = spectrum(alpha_ref, alpha_alt, individuals, 0.0);
                    let chromosomes = 2 * individuals as usize;
                    let total = alpha_ref + alpha_alt;
                    for (class, got) in out.iter().enumerate() {
                        let want = (log_binomial(chromosomes, class)
                            + lgamma(alpha_alt + class as f64)
                            + lgamma(alpha_ref + (chromosomes - class) as f64)
                            + lgamma(total)
                            - lgamma(alpha_alt)
                            - lgamma(alpha_ref)
                            - lgamma(total + chromosomes as f64))
                        .exp();
                        assert!(
                            (got - want).abs() < 1e-12,
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, \
                             class {class}: got {got}, the beta-binomial says {want}"
                        );
                    }
                }
            }
        }
    }

    /// **At full inbreeding every individual is one chromosome counted twice**, so no odd class
    /// can hold anything and the even ones are a beta-binomial over `N` draws rather than `2N`.
    ///
    /// The oracle at the other end of the `F` range, and the check that the doubling is a doubling
    /// rather than a reweighting: a panel of 26 selfers carries the same information as 26
    /// chromosomes, not 52.
    #[test]
    fn at_full_inbreeding_only_the_even_classes_can_hold_anything() {
        for individuals in [1_u32, 4, 26] {
            for alpha_alt in [1e-3, 0.2] {
                let out = spectrum(1.0, alpha_alt, individuals, 1.0);
                for (class, weight) in out.iter().enumerate() {
                    if class % 2 == 1 {
                        assert_eq!(
                            *weight, 0.0,
                            "{individuals} individuals at F = 1: class {class} is odd and must be \
                             empty, got {weight}"
                        );
                    }
                }
                let n = individuals as usize;
                let total = 1.0 + alpha_alt;
                for pairs in 0..=n {
                    let want = (log_binomial(n, pairs)
                        + lgamma(alpha_alt + pairs as f64)
                        + lgamma(1.0 + (n - pairs) as f64)
                        + lgamma(total)
                        - lgamma(alpha_alt)
                        - lgamma(1.0)
                        - lgamma(total + n as f64))
                    .exp();
                    assert!(
                        (out[2 * pairs] - want).abs() < 1e-12,
                        "{individuals} individuals at F = 1, class {}: got {}, a beta-binomial \
                         over {n} draws says {want}",
                        2 * pairs,
                        out[2 * pairs]
                    );
                }
            }
        }
    }

    /// **Doubletons outnumbering singletons is the signature no independent-chromosome spectrum
    /// can produce**, and it is what this function exists to get right (spec §4.1).
    ///
    /// An independent-chromosome spectrum falls monotonically — about `θ/k`, so singletons always
    /// beat doubletons at every diversity. An inbred panel does the opposite, because a selfer
    /// that carries the allele at all usually carries it twice. On the tomato cohort VCF over 26
    /// accessions the real counts are 10,786 doubletons against 5,142 singletons, 2.1 to 1.
    ///
    /// Measured here at that panel size and tomato's fitted diversity: the ratio is 0.50 at
    /// `F = 0`, passes 1 between `F = 0.5` and `F = 0.6`, and reaches 4.7 at `F = 0.9`. So the
    /// real panel's 2.1 sits inside the range this model spans, which the independent-chromosome
    /// one cannot reach at any diversity — checked here across four decades of it.
    #[test]
    fn doubletons_beat_singletons_only_once_the_panel_is_inbred() {
        let individuals = 26;
        let alpha_alt = 6e-4;

        // Independent chromosomes: singletons always win, whatever the diversity.
        for alpha_alt in [1e-5, 1e-4, 1e-3, 1e-2, 0.1] {
            let out = spectrum(1.0, alpha_alt, individuals, 0.0);
            assert!(
                out[2] < out[1],
                "at F = 0 and α_alt {alpha_alt} the doubleton class reached {} against the \
                 singleton's {} — an independent-chromosome spectrum cannot do that",
                out[2],
                out[1]
            );
        }

        // Inbred: the ratio rises with F and crosses one.
        let mut previous = 0.0_f64;
        for inbreeding in [0.0, 0.25, 0.5, 0.6, 0.8, 0.9] {
            let out = spectrum(1.0, alpha_alt, individuals, inbreeding);
            let ratio = out[2] / out[1];
            assert!(
                ratio > previous,
                "the doubleton-to-singleton ratio fell from {previous} to {ratio} when F rose to \
                 {inbreeding}"
            );
            previous = ratio;
        }
        assert!(
            previous > 4.0,
            "at F = 0.9 doubletons should outnumber singletons several times over, got {previous}"
        );
        // And the panel's own measured 2.1 is reached inside the fitted range, not beyond it.
        let at_half = spectrum(1.0, alpha_alt, individuals, 0.5);
        let at_six = spectrum(1.0, alpha_alt, individuals, 0.6);
        assert!(
            at_half[2] / at_half[1] < 1.0 && at_six[2] / at_six[1] > 1.0,
            "the crossing should sit between F = 0.5 and F = 0.6, got {} and {}",
            at_half[2] / at_half[1],
            at_six[2] / at_six[1]
        );
    }

    /// **The neutral shape is `θ/k`, and this reproduces it in the limit the spec names** — small
    /// diversity, independent chromosomes, `α_ref = 1`.
    ///
    /// The spec is careful that the two are *not* an identity: the gap is the panel's own chance
    /// of being polymorphic, about `θ · H(2N)` with `H` the harmonic number. Measured here at
    /// tomato's fitted diversity of 6 in 10,000 over 52 chromosomes, that predicted gap is 3 in a
    /// thousand, and the largest departure of any class from `θ/k` is 0.272% — inside it.
    ///
    /// **This is why the fit's targets are built from this function and never from `θ/k`**
    /// (spec §12 test 5): at a diversity of 1 in 100 the same departure is 4.4%, which is larger
    /// than the effect step D2 is trying to measure.
    #[test]
    fn the_neutral_shape_appears_in_the_small_diversity_limit() {
        let individuals = 26_u32;
        let chromosomes = 2 * individuals as usize;
        let harmonic: f64 = (1..chromosomes).map(|k| 1.0 / k as f64).sum();

        for (alpha_alt, budget) in [(6e-4, 4e-3), (1e-2, 6e-2)] {
            let out = spectrum(1.0, alpha_alt, individuals, 0.0);
            let mut worst = 0.0_f64;
            for (k, weight) in out.iter().enumerate().take(chromosomes).skip(1) {
                let neutral = alpha_alt / k as f64;
                worst = worst.max((weight - neutral).abs() / neutral);
            }
            assert!(
                worst < budget,
                "at θ {alpha_alt} the worst class departed from θ/k by {worst}, above the {budget} \
                 the panel's own polymorphism rate θ·H(2N) = {} allows",
                alpha_alt * harmonic
            );
            // And the departure is of the size the spec predicts, not merely small: it tracks
            // θ·H(2N) rather than being an order of magnitude off it.
            assert!(
                worst > alpha_alt * harmonic / 10.0,
                "the departure {worst} is far below θ·H(2N) = {}, which would mean the limit is \
                 being reached for the wrong reason",
                alpha_alt * harmonic
            );
        }
    }

    /// **A cohort with no alternative allele at all is a real answer**, and it is every site in
    /// class 0 rather than a division by zero.
    #[test]
    fn a_cohort_with_no_alternative_concentration_puts_every_site_in_the_monomorphic_class() {
        for individuals in [1_u32, 5, 26] {
            for inbreeding in [0.0, 0.7, 1.0] {
                let out = spectrum(1.0, 0.0, individuals, inbreeding);
                assert_eq!(out[0], 1.0);
                assert!(out[1..].iter().all(|c| *c == 0.0), "{out:?}");
            }
        }
    }

    /// **One individual is two chromosomes and the spectrum has three classes**, the low end of
    /// the committed cohort range. At `F = 0` they are the beta-binomial's; at `F = 1` the
    /// heterozygous class is empty, because the individual's two copies are one copy.
    #[test]
    fn a_single_individual_still_has_a_spectrum() {
        let outbred = spectrum(1.0, 1e-3, 1, 0.0);
        assert_eq!(outbred.len(), 3);
        assert!((outbred.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(
            outbred[1] > 0.0,
            "an outbred individual can be heterozygous"
        );

        let selfed = spectrum(1.0, 1e-3, 1, 1.0);
        assert_eq!(selfed[1], 0.0, "at F = 1 there is no heterozygous class");
        assert!((selfed.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// **A mis-sized output buffer is refused in release.** The class count is `2N + 1` and a
    /// caller reusing one buffer across panel sizes is the likeliest way to get it wrong; a short
    /// one would leave the top classes carrying the previous panel's weights.
    #[test]
    #[should_panic(expected = "one entry per allele-count class")]
    fn a_mis_sized_class_buffer_is_refused() {
        let mut out = [f64::NAN; 4];
        fill_expected_spectrum(1.0, 1e-3, 3, InbreedingF::try_new(0.0).unwrap(), &mut out);
    }

    /// **An oracle strictly inside the inbreeding range**, which neither exact test above reaches:
    /// both of those sit at `F = 0` or `F = 1`, where the triple sum collapses to one term per
    /// class and a whole family of wrong models still passes.
    ///
    /// Measured, on this code: replacing `class = draws + doubled` with `class = draws` — deleting
    /// the doubling, which is a different model — still passes the sum-to-one test, the
    /// beta-binomial test *and* the neutral-limit test, because the inner sum is a Vandermonde
    /// identity either way.
    ///
    /// What separates them is the first two moments, which the two-branch model fixes exactly:
    ///
    /// ```text
    /// E[j]  = 2N · E[p]                                   — free of F
    /// E[j²] = 4N²·E[p²] + 2N(1 + F)·(E[p] − E[p²])        — carries F linearly
    /// ```
    ///
    /// The mean is `F`-free because inbreeding rearranges copies between individuals without
    /// creating or destroying any; the second moment is where it shows, because a selfer's two
    /// copies agree. `E[p]` and `E[p²]` are the Beta's own moments, written here from the
    /// concentrations directly.
    ///
    /// **It works at one individual**, where the doubleton test has no two classes to compare —
    /// which is the low end of the committed cohort range.
    #[test]
    fn the_first_two_moments_match_the_two_branch_model_at_every_inbreeding_coefficient() {
        for individuals in [1_u32, 2, 5, 17, 40] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 1e-2), (0.5, 0.5), (3.0, 1.0)] {
                for inbreeding in [1e-6, 0.05, 0.25, 0.5, 0.6, 0.9, 0.999] {
                    let out = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    let n = f64::from(individuals);
                    let total = alpha_ref + alpha_alt;
                    let mean_p = alpha_alt / total;
                    let mean_p_squared = alpha_alt * (alpha_alt + 1.0) / (total * (total + 1.0));

                    let got_mean: f64 = out.iter().enumerate().map(|(j, w)| j as f64 * w).sum();
                    let got_second: f64 = out
                        .iter()
                        .enumerate()
                        .map(|(j, w)| (j * j) as f64 * w)
                        .sum();

                    let want_mean = 2.0 * n * mean_p;
                    let want_second = 4.0 * n * n * mean_p_squared
                        + 2.0 * n * (1.0 + inbreeding) * (mean_p - mean_p_squared);

                    assert!(
                        (got_mean - want_mean).abs() <= 1e-9 * want_mean.max(1.0),
                        "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F {inbreeding}: \
                         mean copies {got_mean}, the model says {want_mean}"
                    );
                    assert!(
                        (got_second - want_second).abs() <= 1e-9 * want_second.max(1.0),
                        "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F {inbreeding}: \
                         second moment {got_second}, the model says {want_second}"
                    );
                }
            }
        }
    }

    /// **The whole model, written a second way and compared class by class** — every assignment of
    /// a branch and a genotype to every individual, enumerated, with each one's chance averaged
    /// over the frequency exactly.
    ///
    /// `5^N` panels, so it is only affordable at a handful of individuals; that is the point of
    /// having it beside the moment identity above, which is cheap at any size but pins two numbers
    /// rather than the whole shape. This shares the Beta-moment formula with the code under test
    /// and shares nothing else — no split by how many individuals were inbred, no hypergeometric
    /// counting, no binomial coefficients at all.
    #[test]
    fn every_panel_enumerated_one_by_one_gives_the_same_spectrum() {
        for individuals in [1_usize, 2, 3, 5] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 0.3), (2.0, 1.5)] {
                for inbreeding in [0.0, 0.3, 0.75, 1.0] {
                    let got = spectrum(alpha_ref, alpha_alt, individuals as u32, inbreeding);
                    let want = enumerated_spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    for (class, (got, want)) in got.iter().zip(&want).enumerate() {
                        assert!(
                            (got - want).abs() <= 1e-12 * want.max(1e-12),
                            "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F \
                             {inbreeding}, class {class}: got {got}, enumeration says {want}"
                        );
                    }
                }
            }
        }
    }

    /// Every `(branch, genotype)` assignment for `individuals` people, each contributing a
    /// monomial in `p` and `1 − p`, averaged over the frequency exactly.
    ///
    /// One individual has five ways to be: identical by descent and alternative (`p`, two copies)
    /// or reference (`1 − p`, none), or outbred with two, one or no alternative copies
    /// (`p²`, `2p(1 − p)`, `(1 − p)²`). The five are walked as base-5 digits.
    fn enumerated_spectrum(
        alpha_ref: f64,
        alpha_alt: f64,
        individuals: usize,
        inbreeding: f64,
    ) -> Vec<f64> {
        // (weight, powers of p, powers of 1 − p, alternative copies)
        let ways: [(f64, usize, usize, usize); 5] = [
            (inbreeding, 1, 0, 2),
            (inbreeding, 0, 1, 0),
            (1.0 - inbreeding, 1 + 1, 0, 2),
            (2.0 * (1.0 - inbreeding), 1, 1, 1),
            (1.0 - inbreeding, 0, 2, 0),
        ];
        let total = alpha_ref + alpha_alt;
        let mut out = vec![0.0; 2 * individuals + 1];
        for panel in 0..5_usize.pow(individuals as u32) {
            let mut weight = 1.0;
            let (mut alt_powers, mut ref_powers, mut class) = (0, 0, 0);
            let mut rest = panel;
            for _ in 0..individuals {
                let (w, a, b, copies) = ways[rest % 5];
                rest /= 5;
                weight *= w;
                alt_powers += a;
                ref_powers += b;
                class += copies;
            }
            if weight == 0.0 {
                continue;
            }
            // The exact average of p^a (1 − p)^b over Beta(α_alt, α_ref).
            let moment = (lgamma(alpha_alt + alt_powers as f64) - lgamma(alpha_alt)
                + lgamma(alpha_ref + ref_powers as f64)
                - lgamma(alpha_ref)
                + lgamma(total)
                - lgamma(total + (alt_powers + ref_powers) as f64))
            .exp();
            out[class] += weight * moment;
        }
        out
    }

    /// **A reference concentration at or below zero is refused in release.** Nothing downstream
    /// could catch it: measured with the check removed, `α_ref = −0.5` at four individuals returns
    /// nine finite, non-negative class weights totalling 1.0097 — a spectrum that looks like a
    /// spectrum. `lgamma`'s own guard is a `debug_assert!`, so release has no other.
    #[test]
    #[should_panic(expected = "the reference concentration must be finite")]
    fn a_reference_concentration_at_zero_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(0.0, 1e-3, 2, InbreedingF::try_new(0.5).unwrap(), &mut out);
    }

    /// **A negative alternative concentration is refused in release**, for the same reason: with
    /// the check removed it returns a plausible spectrum summing to 1.0036.
    #[test]
    #[should_panic(expected = "the alternative concentration must be finite")]
    fn a_negative_alternative_concentration_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(1.0, -1e-3, 2, InbreedingF::try_new(0.5).unwrap(), &mut out);
    }

    /// **A concentration past the point where the sum is computable is refused**, rather than
    /// returning a spectrum that no longer sums to one. Step D2 searches over these two axes, so
    /// the ceiling has to be a check and not a comment.
    #[test]
    #[should_panic(expected = "at most 1e6 chromosomes")]
    fn a_concentration_past_the_computable_range_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(1e15, 1e-3, 2, InbreedingF::try_new(0.9).unwrap(), &mut out);
    }

    /// **The ceiling is where the arithmetic is still good, not where it has already gone.** At the
    /// largest concentration the function accepts, the classes must still sum to one.
    #[test]
    fn the_spectrum_still_sums_to_one_at_the_largest_concentration_accepted() {
        for inbreeding in [0.0, 0.5, 0.9] {
            let out = spectrum(
                MAX_PROJECTION_CONCENTRATION,
                MAX_PROJECTION_CONCENTRATION * 1e-3,
                26,
                inbreeding,
            );
            let total: f64 = out.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-9,
                "at the ceiling and F {inbreeding} the classes carry {total}"
            );
        }
    }

    /// **Skipping the rarest branch splits drops no mass**, which is what makes it a saving rather
    /// than an approximation.
    ///
    /// At an inbreeding coefficient of one in a million, almost every split with more than a
    /// handful of inbred individuals has a probability below `1e-300`, so the skip fires on nearly
    /// the whole range of `M`. Two things must survive that, and neither is about speed:
    ///
    /// - **the classes still sum to one** — a skip that dropped real mass would show here first;
    /// - **the spectrum is within order `F` of the outbred one**, because that is how far apart the
    ///   two genuinely are. Measured at three individuals, the largest class moves by 3.0e-10
    ///   between `F = 0` and `F = 1e-6`, so the bound below is loose by four orders of magnitude
    ///   and is a sanity check rather than the oracle.
    ///
    /// **The oracle for the skip is the moment identity above**, which now runs at `F = 1e-6`: it
    /// pins two numbers exactly, and dropping a branch split that carried real weight moves them.
    #[test]
    fn a_vanishing_inbreeding_coefficient_drops_no_mass() {
        for individuals in [3_u32, 26, 63] {
            let outbred = spectrum(1.0, 6e-4, individuals, 0.0);
            let nearly = spectrum(1.0, 6e-4, individuals, 1e-6);
            let total: f64 = nearly.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-12,
                "{individuals} individuals at F = 1e-6: the classes carry {total}"
            );
            for (class, (a, b)) in outbred.iter().zip(&nearly).enumerate() {
                assert!(
                    (a - b).abs() <= 1e-4 * a.max(1e-12),
                    "{individuals} individuals, class {class}: F = 0 gives {a}, F = 1e-6 gives {b}"
                );
            }
        }
    }

    /// **What one prediction costs, at the panel sizes the caller commits to.**
    ///
    /// Ignored by default because it measures wall-clock and would be flaky in a suite; run it
    /// with `cargo test --release --lib -- --ignored --nocapture seed_generic::tests::cost`. The
    /// numbers it prints are what the function's doc quotes, and re-running it is how they stay
    /// true.
    ///
    /// **Read it as cost per objective evaluation, not per run.** Step D2 searches this function,
    /// so a fit pays it on the order of a hundred times over.
    #[test]
    #[ignore = "measures wall-clock; run explicitly in release"]
    fn cost_of_one_prediction_by_panel_size() {
        for individuals in [26_u32, 63, 200, 400, 800] {
            let mut out = vec![0.0; 2 * individuals as usize + 1];
            let started = std::time::Instant::now();
            fill_expected_spectrum(
                1.0,
                6e-4,
                individuals,
                InbreedingF::try_new(0.8).unwrap(),
                &mut out,
            );
            let elapsed = started.elapsed();
            println!(
                "COST individuals={individuals} terms≈{} elapsed={:?}",
                (individuals as u64).pow(3) / 3,
                elapsed
            );
            assert!((out.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }

    /// The same sum, written the slow obvious way: **one exponential per term, nothing stepped and
    /// nothing skipped**.
    ///
    /// This is what the shipped function computed before it was made fast, and it is kept because
    /// it is the version whose agreement with three outside oracles was established — enumeration,
    /// quadrature and a generating function. Every term is built independently from logarithms, so
    /// a term that underflows was genuinely negligible; nothing here can lose a whole row the way
    /// a multiplicative walk can.
    fn term_by_term_spectrum(
        alpha_ref: f64,
        alpha_alt: f64,
        individuals: u32,
        inbreeding: f64,
    ) -> Vec<f64> {
        let n = individuals as usize;
        let mut out = vec![0.0; 2 * n + 1];
        if alpha_alt == 0.0 {
            out[0] = 1.0;
            return out;
        }
        let concentration_total = alpha_ref + alpha_alt;
        let log_pair_constant = lgamma(concentration_total) - lgamma(alpha_alt) - lgamma(alpha_ref);
        for identical_by_descent in 0..=n {
            let Some(log_branch_weight) = log_branch_split(n, identical_by_descent, inbreeding)
            else {
                continue;
            };
            let distinct = 2 * n - identical_by_descent;
            let singles = distinct - identical_by_descent;
            let log_draw_constant =
                log_pair_constant - lgamma(concentration_total + distinct as f64);
            for alternative_draws in 0..=distinct {
                let log_frequency_weight = lgamma(alpha_alt + alternative_draws as f64)
                    + lgamma(alpha_ref + (distinct - alternative_draws) as f64)
                    + log_draw_constant;
                let lowest = alternative_draws.saturating_sub(singles);
                let highest = identical_by_descent.min(alternative_draws);
                for doubled in lowest..=highest {
                    out[alternative_draws + doubled] += (log_branch_weight
                        + log_binomial(identical_by_descent, doubled)
                        + log_binomial(singles, alternative_draws - doubled)
                        + log_frequency_weight)
                        .exp();
                }
            }
        }
        out
    }

    /// **The fast sum and the slow obvious one agree class by class**, over panel sizes, both
    /// concentrations and the whole inbreeding range.
    ///
    /// The shipped function steps its innermost factor by a ratio, tabulates its factorials and
    /// drops branch splits far out in their own tail. None of that changes the model, and this is
    /// where that claim is checked rather than asserted. Measured worst disagreement on this grid
    /// is 2.5e-14 relative.
    ///
    /// **Panel sizes stop at 120 because the slow sum is cubic and this runs in debug.** The sizes
    /// where the difference between the two actually bites are covered by the release-only test
    /// below and by `examples/ng_spectrum_projection_cost.rs`.
    #[test]
    fn the_fast_sum_matches_the_term_by_term_sum() {
        let mut worst = 0.0_f64;
        for individuals in [1_u32, 2, 5, 26, 63, 120] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 1e-2), (0.5, 0.5), (3.0, 2.0)] {
                for inbreeding in [0.0, 1e-6, 0.25, 0.8, 0.999, 1.0] {
                    let fast = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    let slow = term_by_term_spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    for (class, (a, b)) in fast.iter().zip(&slow).enumerate() {
                        let gap = if *b > 0.0 { (a - b).abs() / b } else { a.abs() };
                        worst = worst.max(gap);
                        assert!(
                            gap < 1e-11,
                            "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F \
                             {inbreeding}, class {class}: fast {a}, term-by-term {b}"
                        );
                    }
                }
            }
        }
        assert!(worst < 1e-12, "worst disagreement was {worst}");
    }

    /// **A panel large enough that the hypergeometric walk's low end underflows.**
    ///
    /// The walk is multiplicative, so a first weight of zero makes the whole row zero. Started at
    /// the low end rather than the mode, that is what happens above about a thousand individuals —
    /// measured at 1,600, one class came back 5.7e-16 against its true 6.1e-7 and the spectrum lost
    /// 3 parts in 10,000 of its mass, with every entry still finite and non-negative.
    ///
    /// Ignored by default because a panel that large costs seconds in a debug build; run it with
    /// `cargo test --release --lib -- --ignored --nocapture seed_generic::tests::the_hypergeometric`.
    #[test]
    #[ignore = "needs a panel of 1,600; run explicitly in release"]
    fn the_hypergeometric_walk_survives_a_panel_that_underflows_its_low_end() {
        let out = spectrum(1.0, 6e-4, 1600, 0.8);
        let mass: f64 = out.iter().sum();
        assert!(
            (mass - 1.0).abs() < 1e-9,
            "the classes carry {mass}; a row lost to underflow shows here first"
        );
        // Every class a panel this size can reach carries something: the walk reaching zero early
        // is what empties one.
        let emptiest = out[..=2 * 1600]
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        assert!(
            emptiest > 1e-30,
            "the emptiest class holds {emptiest}, which is a row the walk failed to fill"
        );
    }
}

#[cfg(test)]
mod projection_tests {
    use super::*;

    fn project(weights: &[f64], inbreeding: f64) -> SpectrumSeed {
        project_spectrum_seed(
            Some(FittedSpectrum::new(weights, 10.0, 3_000.0)),
            Some(ExpectedHeterozygosity::try_new(1e-3).unwrap()),
            InbreedingF::try_new(inbreeding).unwrap(),
        )
    }

    /// **A neutral panel projects to `(1, θ)`** — spec §12 test 5, and the reason the whole step
    /// exists: the seed the caller starts from must be the one the pre-pass's spectrum implies,
    /// not a number anybody chose.
    ///
    /// The target is the exact expected spectrum of `(1, θ)` at each panel size, so the answer is
    /// known to the last bit and what is being tested is the search. Measured worst departure
    /// across every combination below: **0.25% on `α_ref` and 0.31% on `α_alt`**, both at 5
    /// individuals, against a search that resolves each concentration to about 1% of itself
    /// ([`SearchPrecision::fast`]). Asked for a thousand times finer resolution the same search
    /// returns both to within **2 parts in 100,000**, more than a hundred times closer, which is
    /// what says the residue above is the resolution and not a bias — that is the second half of
    /// this test.
    ///
    /// **Panel sizes from one individual to 150, diversities from 1 in 10,000 to 1 in 100.** One
    /// individual is the low end of the committed cohort range and the hardest case here: its
    /// spectrum has three classes, of which one holds 999 parts in a thousand.
    #[test]
    fn a_neutral_panel_projects_to_one_and_theta() {
        let mut worst_alpha_ref: f64 = 0.0;
        let mut worst_alpha_alt: f64 = 0.0;
        let mut worst_at = (0u32, 0.0f64, 0.0f64);
        let mut worst_fine_alpha_ref: f64 = 0.0;
        let mut worst_fine_alpha_alt: f64 = 0.0;
        for individuals in [1, 2, 5, 26, 63, 150] {
            for theta in [1e-4, 6e-4, 1e-2] {
                for inbreeding in [0.0, 0.6] {
                    let weights = exact_spectrum(1.0, theta, individuals, inbreeding);
                    let seed = project(&weights, inbreeding);
                    let off_ref = (seed.alpha_ref() - 1.0).abs();
                    let off_alt = (seed.alpha_alt_total() / theta - 1.0).abs();
                    if off_ref.max(off_alt) > worst_alpha_ref.max(worst_alpha_alt) {
                        worst_at = (individuals, theta, inbreeding);
                    }
                    worst_alpha_ref = worst_alpha_ref.max(off_ref);
                    worst_alpha_alt = worst_alpha_alt.max(off_alt);
                }
            }
        }
        // The finer search costs about three times as much per fit, so it runs at both ends of
        // the panel-size range and the middle rather than over the whole grid: what it is here
        // to say is that the residue above shrinks with the resolution asked for, and that shows
        // wherever it is measured.
        for individuals in [1, 26, 150] {
            for inbreeding in [0.0, 0.6] {
                let theta = 6e-4;
                let weights = exact_spectrum(1.0, theta, individuals, inbreeding);
                let sharp = fit_pair(
                    &FittedSpectrum::new(&weights, 10.0, 3_000.0),
                    InbreedingF::try_new(inbreeding).unwrap(),
                    SearchPrecision::fine(),
                );
                worst_fine_alpha_ref = worst_fine_alpha_ref.max((sharp.alpha_ref - 1.0).abs());
                worst_fine_alpha_alt =
                    worst_fine_alpha_alt.max((sharp.alpha_alt / theta - 1.0).abs());
            }
        }
        assert!(
            worst_alpha_ref < 1e-2 && worst_alpha_alt < 1e-2,
            "the pair must come back inside the 1% the shipped search resolves; worst was \
             {worst_alpha_ref:.2e} on α_ref and {worst_alpha_alt:.2e} on α_alt, at {} \
             individuals, θ = {}, F = {}",
            worst_at.0,
            worst_at.1,
            worst_at.2
        );
        assert!(
            worst_fine_alpha_ref < 5e-5 && worst_fine_alpha_alt < 5e-5,
            "asked for a thousand times finer resolution the same search must reach (1, θ) far \
             closer — otherwise the residue above is a bias, not the resolution; \
             worst was {worst_fine_alpha_ref:.2e} on α_ref and {worst_fine_alpha_alt:.2e} on α_alt"
        );
    }

    /// **On a spectrum the family can reproduce exactly, the winning score is the spectrum's own
    /// entropy** — because the objective is that entropy minus the Kullback–Leibler divergence to
    /// the prediction, and at the true pair the divergence is zero.
    ///
    /// This is the only test that reads the score the fit carries back, and it is what makes that
    /// number a check rather than a value nobody looks at.
    #[test]
    fn the_winning_score_is_the_spectrums_own_entropy() {
        for individuals in [1u32, 26] {
            let weights = exact_spectrum(1.0, 6e-4, individuals, 0.6);
            let entropy: f64 = weights
                .iter()
                .filter(|weight| **weight > 0.0)
                .map(|weight| weight * weight.ln())
                .sum();
            let fit = fit_pair(
                &FittedSpectrum::new(&weights, 10.0, 3_000.0),
                InbreedingF::try_new(0.6).unwrap(),
                SearchPrecision::fine(),
            );
            assert!(
                (fit.log_likelihood - entropy).abs() < 1e-9,
                "at {individuals} individuals the fit scored {} where the spectrum's own entropy \
                 is {entropy}",
                fit.log_likelihood
            );
        }
    }

    /// **One density projects to one pair whatever the panel's inbreeding** — spec §12 test 6,
    /// and the test that holds §4.1's two-branch requirement in place rather than leaving it as
    /// prose.
    ///
    /// A panel's `2N` chromosomes are not `2N` independent draws once its individuals are inbred.
    /// Predicting as though they were biases the reference concentration **down**, with a fixed
    /// sign, and this measures how far on a panel of 26 individuals at tomato's diversity:
    ///
    /// ```text
    ///   F                       0        0.6      0.8      0.9
    ///   two-branch  α_ref     1.0000   1.0000   1.0000   1.0000
    ///   independent α_ref     1.0000   0.9144   0.8793   0.8599
    ///                            —     8.6% low  12.1%    14.0%
    ///   independent α_alt     1.000 θ  0.893 θ  0.848 θ  0.824 θ
    /// ```
    ///
    /// **The `F = 0` column is the comparison's zero**: there the two predictions are the same
    /// model, so a difference between the rows would mean the arms are not being run on one
    /// density.
    ///
    /// **The last row is why spec §4.1's tomato remark is not evidence of anything yet.** A fit
    /// against an independently-called VCF of 18 accessions returned `α_alt = 0.81 θ` and was
    /// read as a hint that a domesticated selfer stretches the two-parameter family. A perfectly
    /// neutral panel run through an independent-chromosome projection returns 0.824 θ at
    /// `F = 0.9` by itself, so that number is consistent with being nothing but this bias.
    #[test]
    fn the_projection_returns_one_pair_at_every_inbreeding_coefficient() {
        let theta = 6e-4;
        let mut two_branch = Vec::new();
        let mut independent = Vec::new();
        for inbreeding in [0.0, 0.6, 0.8, 0.9] {
            let weights = exact_spectrum(1.0, theta, 26, inbreeding);
            let spectrum = FittedSpectrum::new(&weights, 10.0, 3_000.0);
            two_branch.push(fit_pair(
                &spectrum,
                InbreedingF::try_new(inbreeding).unwrap(),
                SearchPrecision::fine(),
            ));
            independent.push(fit_pair(
                &spectrum,
                InbreedingF::try_new(0.0).unwrap(),
                SearchPrecision::fine(),
            ));
        }

        for (fit, inbreeding) in two_branch.iter().zip([0.0, 0.6, 0.8, 0.9]) {
            assert!(
                (fit.alpha_ref - 1.0).abs() < 1e-5 && (fit.alpha_alt / theta - 1.0).abs() < 1e-5,
                "carrying the panel's F, one density must give one pair at every F; at \
                 F = {inbreeding} it gave ({}, {})",
                fit.alpha_ref,
                fit.alpha_alt
            );
        }

        for ((fit, expected), inbreeding) in independent
            .iter()
            .zip([1.0, 0.9144, 0.8793, 0.8599])
            .zip([0.0, 0.6, 0.8, 0.9])
        {
            assert!(
                (fit.alpha_ref - expected).abs() < 5e-4,
                "an independent-chromosome projection must still return its own biased answer — \
                 if this moves, the two-branch numbers above are being compared against \
                 something else; at F = {inbreeding} expected {expected}, got {}",
                fit.alpha_ref
            );
        }
        assert!(
            (independent[3].alpha_alt / theta - 0.8235).abs() < 5e-4,
            "the alternative concentration's own bias at F = 0.9 is 0.824 θ, got {} θ",
            independent[3].alpha_alt / theta
        );
    }

    /// **One individual projects to the neutral pair, reached without any test of the cohort
    /// size** — spec §12 test 7.
    ///
    /// With one sample no census site is variable across the panel, so the pre-pass's spectrum is
    /// its own regularizer at that genome's measured diversity, and projecting it returns the two
    /// numbers §4 sets. The only branch this function has is on the spectrum being **absent**;
    /// nothing here reads a sample count.
    ///
    /// It is also the case the search finds hardest, and it is what put the search on the total
    /// and the ratio rather than on the two concentrations: in the latter it returned
    /// `α_ref = 0.844` here (see [`concentration_pair`]).
    #[test]
    fn at_one_individual_the_projection_is_still_the_neutral_pair() {
        let theta = 6e-4;
        for inbreeding in [0.0, 0.9] {
            let weights = exact_spectrum(1.0, theta, 1, inbreeding);
            assert_eq!(
                weights.len(),
                3,
                "one individual has three allele-count classes"
            );
            let seed = project(&weights, inbreeding);
            assert!(
                (seed.alpha_ref() - 1.0).abs() < 5e-3
                    && (seed.alpha_alt_total() / theta - 1.0).abs() < 5e-3,
                "at F = {inbreeding} one individual must still project to (1, θ); got ({}, {})",
                seed.alpha_ref(),
                seed.alpha_alt_total()
            );
        }
    }

    /// **No spectrum: the pair is the neutral `(1, θ)` at the diversity the pre-pass did fit** —
    /// exactly, with no arithmetic in between, and the regime says where it came from.
    ///
    /// The pre-pass emits the spectrum as *absent* below a panel-size floor rather than as a thin
    /// estimate, so a cohort of five arrives here without one while a single sample arrives with
    /// one. **A branch on absence, never on cohort size.**
    #[test]
    fn an_absent_spectrum_is_the_neutral_pair_at_the_fitted_diversity() {
        let theta = ExpectedHeterozygosity::try_new(6e-4).unwrap();
        let seed = project_spectrum_seed(None, Some(theta), InbreedingF::try_new(0.85).unwrap());
        assert_eq!(seed.alpha_ref(), 1.0);
        assert_eq!(seed.alpha_alt_total(), 6e-4);
        assert_eq!(seed.regime(), SeedRegime::NeutralShape);
    }

    /// **No fitted diversity either: the species-range fallback, and the run must say so.** Two
    /// runs that used different information are otherwise indistinguishable in what they emit.
    #[test]
    fn no_fitted_diversity_falls_back_to_the_species_value_and_says_so() {
        let seed = project_spectrum_seed(None, None, InbreedingF::try_new(0.0).unwrap());
        assert_eq!(seed.alpha_ref(), 1.0);
        assert_eq!(
            seed.alpha_alt_total(),
            ExpectedHeterozygosity::SPECIES_FALLBACK.get()
        );
        assert_eq!(seed.regime(), SeedRegime::FallbackDiversity);
    }

    /// **How strongly the estimate was held toward the neutral shape travels with the seed**, and
    /// so does whether the real sites outweighed that hold. A run whose spectrum was mostly its
    /// own prior and one whose spectrum was measured are otherwise indistinguishable in what they
    /// emit (spec §4.1).
    #[test]
    fn the_regularizer_weight_and_whether_the_census_sites_outweighed_it_travel_with_the_seed() {
        let weights = exact_spectrum(1.0, 6e-4, 5, 0.0);
        for (regularizer, variable_sites, expected_data_won) in
            [(10.0, 3_000.0, true), (10_000.0, 3_000.0, false)]
        {
            let seed = project_spectrum_seed(
                Some(FittedSpectrum::new(&weights, regularizer, variable_sites)),
                Some(ExpectedHeterozygosity::try_new(6e-4).unwrap()),
                InbreedingF::try_new(0.0).unwrap(),
            );
            assert_eq!(
                seed.regime(),
                SeedRegime::FittedSpectrum {
                    regularizer_site_weight: regularizer,
                    census_sites_outweigh_regularizer: expected_data_won,
                    spectrum_match: SpectrumMatch::Reproduced,
                }
            );
        }
    }

    /// **Every start reaches the same pair when each is swept to its own convergence.**
    ///
    /// Only the best-scoring start is swept to convergence in the shipped path, so this is what
    /// says that saving is safe rather than a search that stopped early in a place the other
    /// starts would have left. Measured, `α_ref` across the four starts: identical to 15 digits
    /// at one individual and at 63, and 1.0039 apart at 26 — inside the 1% the search resolves.
    ///
    /// **A spread near one is evidence of nothing on its own**, because a line search over the
    /// whole box overwrites a start's own value before reading it (`fitting/multistart.rs`); it
    /// is paired here with [`a_neutral_panel_projects_to_one_and_theta`], which scores the point
    /// against an answer that is known.
    #[test]
    fn every_start_reaches_the_same_pair_within_one_percent() {
        for (individuals, inbreeding) in [(1u32, 0.9), (26, 0.6), (63, 0.85)] {
            let weights = exact_spectrum(1.0, 6e-4, individuals, inbreeding);
            let spectrum = FittedSpectrum::new(&weights, 10.0, 3_000.0);
            let panel_inbreeding = InbreedingF::try_new(inbreeding).unwrap();
            let mut scorer = SpectrumScorer::new(&spectrum, panel_inbreeding);
            let reached: Vec<f64> = SEARCH_STARTS
                .iter()
                .map(|(total, ratio)| {
                    let start = ScoredPoint {
                        point: [total.ln(), ratio.ln()],
                        log_likelihood: f64::NEG_INFINITY,
                    };
                    let best = sweep_from(start, SearchPrecision::fast(), &mut |point| {
                        scorer.score(point)
                    });
                    concentrations_at(best.point).0
                })
                .collect();
            let highest = reached.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let lowest = reached.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(
                highest / lowest < 1.01,
                "at {individuals} individuals and F = {inbreeding} the starts reached {reached:?}"
            );
        }
    }

    /// **A spectrum the family can reproduce comes back marked as reproduced**, which is the
    /// baseline the two markers below are worth anything against.
    #[test]
    fn a_spectrum_the_family_can_reproduce_says_so() {
        for (individuals, inbreeding) in [(1u32, 0.0), (26, 0.6), (63, 0.9)] {
            let weights = exact_spectrum(1.0, 6e-4, individuals, inbreeding);
            let seed = project(&weights, inbreeding);
            assert_eq!(
                seed.regime(),
                SeedRegime::FittedSpectrum {
                    regularizer_site_weight: 10.0,
                    census_sites_outweigh_regularizer: true,
                    spectrum_match: SpectrumMatch::Reproduced,
                },
                "at {individuals} individuals and F = {inbreeding}"
            );
        }
    }

    /// **A spectrum no pair can produce says so instead of returning a pair that looks fitted.**
    ///
    /// At an inbreeding coefficient of 1 the model puts exactly nothing on an odd number of
    /// chromosomes carrying the allele — every individual's two copies are one copy counted
    /// twice — so a measured spectrum holding any heterozygote cannot have come from any pair.
    /// Measured before this marker existed: all 441 points across the search box scored the same,
    /// and the run returned whichever one it happened to reach.
    #[test]
    fn a_spectrum_no_pair_can_produce_is_marked_rather_than_answered() {
        let weights = [0.90, 0.04, 0.04, 0.01, 0.01];
        let seed = project_spectrum_seed(
            Some(FittedSpectrum::new(&weights, 10.0, 3_000.0)),
            Some(ExpectedHeterozygosity::try_new(6e-4).unwrap()),
            InbreedingF::try_new(1.0).unwrap(),
        );
        assert!(
            matches!(
                seed.regime(),
                SeedRegime::FittedSpectrum {
                    spectrum_match: SpectrumMatch::Unreproducible,
                    ..
                }
            ),
            "got {:?}",
            seed.regime()
        );
        assert!(seed.alpha_ref().is_finite() && seed.alpha_alt_total().is_finite());
    }

    /// **A fit that stopped on the edge of the range searched says so**, because a better pair
    /// may lie outside it and a boundary is not a summit.
    ///
    /// A fully invariant cohort reaches this legitimately and is the cheapest way to show it: its
    /// answer is an alternative concentration of zero, which the search cannot express — the
    /// ratio between the two concentrations floors at `1e-9`. What comes back is that floor, and
    /// without the marker it would read as a fitted answer.
    #[test]
    fn a_fit_that_stopped_at_the_edge_of_its_range_says_so() {
        let mut weights = vec![0.0; 53];
        weights[0] = 1.0;
        let seed = project_spectrum_seed(
            Some(FittedSpectrum::new(&weights, 10.0, 3_000.0)),
            Some(ExpectedHeterozygosity::try_new(6e-4).unwrap()),
            InbreedingF::try_new(0.0).unwrap(),
        );
        assert!(
            matches!(
                seed.regime(),
                SeedRegime::FittedSpectrum {
                    spectrum_match: SpectrumMatch::AtSearchLimit,
                    ..
                }
            ),
            "got {:?}",
            seed.regime()
        );
    }

    /// **The objective is maximised by the truth and by nothing else** — Gibbs' inequality, which
    /// is what makes "maximum likelihood of the class weights" and "closest in Kullback–Leibler
    /// divergence" the same instruction (spec §4.1).
    ///
    /// An independent check of the objective rather than of the search: it uses no spectrum at
    /// all, only distributions built here.
    #[test]
    fn the_objective_is_maximised_by_the_truth() {
        let truth = [0.7, 0.2, 0.06, 0.03, 0.01];
        let at_truth = spectrum_log_likelihood(&truth, &truth);
        for candidate in [
            [0.69, 0.21, 0.06, 0.03, 0.01],
            [0.2, 0.2, 0.2, 0.2, 0.2],
            [0.7, 0.2, 0.05, 0.04, 0.01],
            [0.999, 0.0005, 0.0002, 0.0002, 0.0001],
        ] {
            assert!(
                spectrum_log_likelihood(&truth, &candidate) < at_truth,
                "no candidate may score above the truth; {candidate:?} did"
            );
        }
    }

    /// **A class the fitted spectrum gives no weight contributes nothing, and one it does give
    /// weight to that the candidate predicts at zero is floored rather than sent to `−∞`.**
    ///
    /// Both halves are ordinary here rather than edge cases: at full inbreeding every odd class
    /// is exactly zero. The first would be `0 × ln 0`, a `NaN`, which compares false against
    /// every score and so would make a candidate invisible to the search rather than rejected.
    ///
    /// The second is the property that lets the search move: an impossible candidate must score
    /// **below every possible one and still be ordered against other impossible ones**, so that a
    /// line crossing a region where the prediction underflows still tells the search which way to
    /// go. `−∞` gives no such ordering, and there are spectra on which every candidate is
    /// impossible — see [`spectrum_log_likelihood`].
    #[test]
    fn an_impossible_class_is_floored_rather_than_sent_to_negative_infinity() {
        assert_eq!(
            spectrum_log_likelihood(&[0.5, 0.0, 0.5], &[0.5, 0.0, 0.5]),
            2.0 * 0.5 * 0.5f64.ln(),
            "a class with no weight must be skipped, not multiplied into a NaN"
        );

        let impossible = spectrum_log_likelihood(&[0.5, 0.5, 0.0], &[0.5, 0.0, 0.5]);
        assert_eq!(
            impossible,
            0.5 * 0.5f64.ln() + 0.5 * PROBABILITY_FLOOR.ln(),
            "a candidate that cannot have produced an occupied class pays the floor for it"
        );
        assert!(
            impossible.is_finite(),
            "the score has to be finite for the search to compare it against anything"
        );
        assert!(
            impossible < spectrum_log_likelihood(&[0.5, 0.5, 0.0], &[0.5, 1e-200, 0.5]),
            "an impossible candidate must rank below one that merely makes the class very \
             unlikely, or the search cannot walk out of the region"
        );
        assert!(
            spectrum_log_likelihood(&[0.5, 0.5, 0.0], &[0.5, 0.5, 0.0]) > impossible,
            "and below the candidate that matches"
        );
    }

    /// **An even class count is refused.** `2N + 1` classes is odd; an even count would silently
    /// become a panel one individual smaller, and every class would then be read off by one.
    #[test]
    #[should_panic(expected = "allele-count classes")]
    fn an_even_class_count_is_refused() {
        let _ = FittedSpectrum::new(&[0.5, 0.3, 0.15, 0.05], 10.0, 3_000.0);
    }

    /// **A panel of no individuals is refused**, which an odd count of one otherwise passes for.
    /// Measured with the lower bound removed: a single-class spectrum returns `α_ref = 996.87`,
    /// which is the top corner of the search box and not a fit at all.
    #[test]
    #[should_panic(expected = "at least 3")]
    fn a_spectrum_with_one_class_is_refused() {
        let _ = FittedSpectrum::new(&[1.0], 10.0, 3_000.0);
    }

    /// **A panel past what the projection will fit in a sensible time is refused rather than
    /// started.** At 10,000 individuals a fit is already about two hours; nothing on the way would
    /// say the run had stopped being a run.
    #[test]
    #[should_panic(expected = "at most 20001")]
    fn a_panel_past_the_projections_range_is_refused() {
        let weights = vec![0.0; 2 * MAX_PROJECTION_INDIVIDUALS as usize + 3];
        let mut weights = weights;
        weights[0] = 1.0;
        let _ = FittedSpectrum::new(&weights, 10.0, 3_000.0);
    }

    /// **A negative class weight is refused even when the weights still sum to one.** The
    /// objective skips a class on `weight > 0.0`, so the negative one would simply drop out and
    /// the fit would run against more than a unit of probability. Measured with the check
    /// removed, `[1.5, -0.5, 0.0]` returns `(1.005e-3, 1.008e-12)` — the box's bottom corner,
    /// reported as a fit.
    #[test]
    #[should_panic(expected = "cannot be negative or NaN")]
    fn a_negative_class_weight_is_refused_even_when_the_total_is_one() {
        let _ = FittedSpectrum::new(&[1.5, -0.5, 0.0], 10.0, 3_000.0);
    }

    /// **A regularizer weight that is not a count of sites is refused.** It travels onto the run's
    /// output through [`SeedRegime::FittedSpectrum`], and a `NaN` there is worse than an error:
    /// `NaN > x` is false, so the seed would also report that the prior dominated, which is the
    /// reassuring answer.
    #[test]
    #[should_panic(expected = "sites' worth of pseudo-counts")]
    fn a_regularizer_weight_that_is_not_a_count_of_sites_is_refused() {
        let _ = FittedSpectrum::new(&[0.9, 0.05, 0.05], f64::NAN, 3_000.0);
    }

    /// **A variable-site count that is not a count is refused**, for the same reason: it is the
    /// other half of the comparison the seed reports.
    #[test]
    #[should_panic(expected = "count of variable census sites")]
    fn a_variable_site_count_that_is_not_a_count_is_refused() {
        let _ = FittedSpectrum::new(&[0.9, 0.05, 0.05], 10.0, -1.0);
    }

    /// **A fully inbred panel whose spectrum holds heterozygotes still returns a pair.**
    ///
    /// At `F = 1` the prediction puts exactly zero in every odd allele-count class, so a spectrum
    /// carrying any weight in one of them cannot have come from any pair in this family — 441 of
    /// 441 points across the search box, measured. Before the floor of
    /// [`spectrum_log_likelihood`] every candidate scored `−∞`, no start could beat the sentinel
    /// the search began from, and the run died three frames later complaining that a
    /// concentration was `NaN`.
    ///
    /// `InbreedingF` still admits `1.0`; the `[0, 1)` tightening spec §7 requires belongs to
    /// `calling_prerequisites.md` Milestone A. **What this pins is that the fit survives it**, not
    /// that the answer means anything — see [`spectrum_log_likelihood`] on what is still missing.
    #[test]
    fn a_fully_inbred_panel_whose_spectrum_holds_heterozygotes_still_returns_a_pair() {
        let weights = [0.90, 0.04, 0.04, 0.01, 0.01];
        let seed = project_spectrum_seed(
            Some(FittedSpectrum::new(&weights, 10.0, 3_000.0)),
            Some(ExpectedHeterozygosity::try_new(6e-4).unwrap()),
            InbreedingF::try_new(1.0).unwrap(),
        );
        assert!(
            seed.alpha_ref().is_finite() && seed.alpha_alt_total().is_finite(),
            "got ({}, {})",
            seed.alpha_ref(),
            seed.alpha_alt_total()
        );
    }

    /// **Weights that are not a distribution are refused.** The objective is a log-likelihood
    /// against them, so unnormalised counts score every candidate wrongly by an amount that
    /// varies with the candidate — a fit that returns an answer rather than an error.
    #[test]
    #[should_panic(expected = "must sum to 1 within")]
    fn weights_that_are_not_a_distribution_are_refused() {
        let _ = FittedSpectrum::new(&[500.0, 300.0, 200.0], 10.0, 3_000.0);
    }

    /// **What a fit costs, by panel size** — the number that says whether the projection runs at
    /// the top of the committed cohort range.
    ///
    /// Run with `cargo test --release --lib -- --ignored the_cost_of_one_fit`. Measured:
    ///
    /// ```text
    ///   individuals    400     800    1,600    3,200
    ///   one fit       3.8 s    22 s   2.2 min  11.8 min
    ///   predictions    399     399     399      399
    /// ```
    ///
    /// About `N^2.5`, which is the `N^2.45` one prediction costs plus a constant factor rather
    /// than a power. **The count is printed here and asserted in
    /// [`a_fit_costs_at_most_450_predictions`]**, which runs at panel sizes small enough for the
    /// ordinary test suite; this one exists for the wall clock, which does depend on the machine.
    /// A prediction **inside a fit** averages 1.78 s at 3,200 individuals — not the 0.96 s
    /// `doc/devel/ng/reports/spectrum_projection_cost_2026-08-22.md` measures at the neutral
    /// pair, which the search leaves after its first few steps.
    #[test]
    #[ignore]
    fn the_cost_of_one_fit_by_panel_size() {
        for individuals in [400u32, 800, 1_600, 3_200] {
            let weights = exact_spectrum(1.0, 6e-4, individuals, 0.8);
            let started = std::time::Instant::now();
            let fit = fit_pair(
                &FittedSpectrum::new(&weights, 10.0, 3_000.0),
                InbreedingF::try_new(0.8).unwrap(),
                SearchPrecision::fast(),
            );
            println!(
                "{individuals} individuals: {:?} for {} predictions, α_ref = {:.6}, α_alt = {:.6e}",
                started.elapsed(),
                fit.predictions,
                fit.alpha_ref,
                fit.alpha_alt
            );
        }
    }

    /// **The fit's cost in predictions is flat in the panel size** — 399 at every size and every
    /// inbreeding coefficient measured, here and in [`the_cost_of_one_fit_by_panel_size`] —
    /// which is what lets the wall clock at any size be read off one in-fit per-prediction
    /// measurement. Held here so that a change to the search that makes it cost more has to say
    /// so: adding one more start takes it to 532, adding back the fourth sweep direction takes
    /// it past 600.
    #[test]
    fn a_fit_costs_at_most_450_predictions() {
        let mut most = 0;
        for individuals in [1u32, 5, 26, 63, 150] {
            for inbreeding in [0.0, 0.6, 0.9] {
                let weights = exact_spectrum(1.0, 6e-4, individuals, inbreeding);
                let fit = fit_pair(
                    &FittedSpectrum::new(&weights, 10.0, 3_000.0),
                    InbreedingF::try_new(inbreeding).unwrap(),
                    SearchPrecision::fast(),
                );
                most = most.max(fit.predictions);
            }
        }
        assert!(
            most <= 450,
            "a fit took {most} predictions; at 3,200 individuals a fit is 11.8 minutes measured, \
             so this is the whole run's projection budget"
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
        fill_locus_concentration(seed, VariantClass::Substitution, allele_count, &mut out);
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
        fill_locus_concentration(seed, VariantClass::Substitution, 3, &mut substitution);
        fill_locus_concentration(seed, VariantClass::InsertionOrDeletion, 3, &mut indel);
        assert_eq!(substitution, indel);
    }

    /// **A locus with no alleles is refused.** Every locus has a reference allele, so a count of
    /// zero is a caller that has lost track of which locus it is on. Without the assertion the
    /// write to `out[0]` panics anyway, on an index — the assertion is what makes the panic name
    /// the mistake rather than the line.
    #[test]
    #[should_panic(expected = "lost track of which locus")]
    fn a_locus_with_no_alleles_is_refused() {
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
        fill_locus_concentration(
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
        fill_locus_concentration(
            neutral_seed(6e-4),
            VariantClass::Substitution,
            3,
            &mut too_short,
        );
    }
}
