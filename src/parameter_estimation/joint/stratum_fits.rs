//! **The one borrow of slippage numbers that crosses the calling seam.**
//!
//! The caller scoring a repeat-tract candidate needs three numbers — how often a read reports a
//! tract length other than its allele's, which way the slips go, and how fast multi-repeat slips
//! fall off — for the read group that produced the read and the stratum the candidate sits in
//! (`doc/devel/ng/arch/read_likelihoods.md` §4.2). The pre-pass has already fitted them. What it
//! has not had is a way to ask for one, and this module is that way and nothing more.
//!
//! # What it deliberately does not do
//!
//! **It does not blend, and re-blending here would be a defect rather than a duplication.** The
//! level on a [`StratumOutcome`] is already the emitted one, and its [`LevelProvenance`] says how
//! it was made. The penalised fit (`ssr_fit::penalised_fit::fit_repeat_tracts`) emits a point on
//! its motif length's penalised curve. The climb it replaced (deleted 2026-10-16), whose
//! parameters files are still read, emitted one of three:
//!
//! - a stratum fitted on its own tracts had its level weighed against its period's curve — the
//!   curves were drawn after every stratum had its own answer, then `blend_level` (deleted with
//!   the rest of the curve fitting on 2026-10-16) applied in place;
//! - a stratum too thin to be fitted took its period's curve whole, never through `blend_level`;
//! - a run with curves switched off kept every cell's own answer.
//!
//! Feeding any of those back through `blend_level` would weigh the curve against a number the
//! curve is already inside. `spec/str_slippage_level_curve.md` §5.1 does not name this act — what
//! it forbids in so many words is a *curve* fitted from blended values, "otherwise each round of
//! smoothing fits a curve to the previous round's curve, and the cells stop being evidence". This
//! is the same circularity one step downstream, and the measurement is in the report: re-blending
//! the five blended strata of a small real fit moves their levels by 0.6% to 4.1%, while leaving
//! the `curve_weight` in their provenance unchanged — so a consumer inspecting the provenance
//! would see nothing wrong and only the number would have moved.
//!
//! **It does not decide anything about a stratum with no answer.** Four different absences reach
//! a lookup and [`NoSlippage`] keeps them apart, because a caller that cannot tell *this run
//! never named that read group* from *that library put no read in this stratum* will report a
//! quiet tract as an unsequenced one.
//!
//! # The grain, and why the key has a read group in it
//!
//! Slippage is a property of the chemistry, so it is fitted per read group
//! (`spec/parameter_prepass_joint_fit.md` §4). A run may declare that several of its read groups
//! slip alike by naming them in one **slippage group**, and one group per read group is the
//! **specified grain** — which is not the same as what happens by default: the only builder of
//! that map in this tree pools every read group into one set unless told otherwise
//! (`examples/ng_joint_records_walk.rs`). So a lookup takes the read group the caller has and
//! this module translates, rather than making every caller carry the translation.

use std::collections::BTreeMap;

use crate::parameter_estimation::Provenance;
use crate::types::ReadGroupId;

use super::census::Stratum;
use super::share_curve::ShareSource;
use super::slippage_curve::LevelSource;
use super::ssr_fit::penalised_fit::LOCI_FOR_AN_ERROR;
use super::ssr_fit::{LevelProvenance, SharesProvenance, Slippage, StratumOutcome};

/// What one `(read group, stratum)` cell answers with.
///
/// **The three numbers and where they came from, together.** A level fitted from 8,000 slipped
/// reads and one read off a curve through four cells are the same `f64`, and a consumer that
/// weighs them alike is treating an interpolation as a measurement
/// (`str_slippage_level_curve.md` §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedSlippage {
    /// How often a read reports a tract length other than its allele's, which way the slips go,
    /// and how fast multi-repeat slips fall off — **the emitted numbers**. The level is the one
    /// the fit settled on, which the field beside it says was the cell's, the curve's, or a blend
    /// of the two; it is not always a blend, and this type never makes one.
    pub slippage: Slippage,
    /// Where the level came from: the stratum's own fit, its period's curve, or a blend of the
    /// two with the share the curve carried.
    pub level: LevelProvenance,
    /// Where the direction split and the fall-off came from. Separate from the level's, because
    /// the three numbers are smoothed on their own curves and a stratum can take its level from
    /// a curve while keeping its own shares.
    pub shares: Option<SharesProvenance>,
    /// **The standard error of each of the three numbers, where that number is the stratum's own
    /// fit** (`doc/devel/ng/spec/fit_precision.md` §5.2). Calling never reads them; the parameters
    /// file writes them for the person judging the fit.
    pub own_fit_standard_errors: OwnFitStandardErrors,
    /// **How many samples the stratum's own fit read**, where a large cohort's stratum was fitted
    /// on a subset of its samples (`fit_precision.md` §4.4) — `None` where every sample was read
    /// without a subset being drawn, or nothing was fitted here. The evidence behind the numbers —
    /// [`LevelProvenance::slipped_reads`] among them — is that subset's.
    pub samples_fitted_on: Option<u64>,
}

/// **The standard errors of one slippage group's three numbers in one stratum, each only where the
/// number emitted was fitted with the stratum's own reads**: its own fit, or a point on a
/// penalised curve (`doc/devel/ng/spec/str_fit_to_calling.md` §2.3).
///
/// A number taken whole from its period's parametric curve has none: the own fit's error would
/// describe a number the stratum does not emit. A number blended with the curve has none either,
/// even at a curve weight of 0, where the number is the stratum's own fit: no error is computed for
/// a blend (spec `fit_precision.md` §5.2), the writer and the validator have always refused one
/// beside it, and the climb never emitted a blend at weight 0. Also `None` where the fit gives the
/// number none — not told apart from the stratum's other numbers, not placed by its tracts, or
/// withheld for too few loci.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct OwnFitStandardErrors {
    /// The slippage level's — how often a read reports a tract length other than its allele's.
    pub level: Option<f64>,
    /// The share of slipped reads showing a shorter tract.
    pub shorter_share: Option<f64>,
    /// The fall-off of two-repeat slips against one-repeat slips.
    pub fall_off: Option<f64>,
}

impl OwnFitStandardErrors {
    /// **The fit's errors, kept only for the numbers the stratum emits as fitted with its own
    /// reads** — the level where it came from the cell or a penalised curve, each share where it
    /// came from the stratum or a penalised curve. `errors` is the fit's for this group, `None`
    /// where none was computed.
    fn of_the_numbers_emitted(
        errors: Option<&super::ssr_fit::SlippageErrors>,
        level: Option<&LevelProvenance>,
        shares: Option<&SharesProvenance>,
    ) -> Self {
        let Some(errors) = errors else {
            return Self::default();
        };
        Self {
            level: level
                .filter(|level| a_level_fitted_with_its_strata_reads(level.source))
                .and_then(|_| errors.level.value()),
            shorter_share: shares
                .filter(|shares| a_share_fitted_with_its_strata_reads(shares.shorter_share.source))
                .and_then(|_| errors.shorter_share.value()),
            fall_off: shares
                .filter(|shares| a_share_fitted_with_its_strata_reads(shares.fall_off.source))
                .and_then(|_| errors.fall_off.value()),
        }
    }
}

/// Whether a level from this source was fitted with its stratum's own reads, so that the fit's
/// error describes it. **Every variant is named**, so a new source has to say which side it is
/// on. A blend is not, even at a curve weight of 0: no error is computed for a blend, and the
/// climb, the only fit that emitted blends (deleted 2026-10-16), never emitted one at weight 0.
fn a_level_fitted_with_its_strata_reads(source: LevelSource) -> bool {
    match source {
        LevelSource::Cell | LevelSource::PenalisedCurve(_) => true,
        LevelSource::Curve | LevelSource::Blend { .. } => false,
    }
}

/// The same, for a share.
fn a_share_fitted_with_its_strata_reads(source: ShareSource) -> bool {
    match source {
        ShareSource::Stratum | ShareSource::PenalisedCurve(_) => true,
        ShareSource::Curve | ShareSource::Blend { .. } => false,
    }
}

/// Why a lookup has no numbers.
///
/// **Four absences and not one**, and two of them are ordinary while two say the run is not what
/// it claims. The structural twin is [`NotIdentifiedReason`], which does the same for a
/// contamination fraction: several named reasons behind one absence, because a caller told only
/// "no number" would act on it.
///
/// [`NotIdentifiedReason`]: super::contamination::NotIdentifiedReason
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NoSlippage {
    /// **Ordinary.** No stratum of this period and repeat count is in the fit. Either the cohort
    /// holds no kept tract of that shape, or every one of them was refused. A candidate several
    /// repeats from its reference tract's length can land here on perfectly good data, so a
    /// caller has to have an answer for it.
    NoSuchStratum,
    /// **Ordinary.** The stratum is in the fit, and this read group's slippage group put no read
    /// in it. The group is named because it is what was looked under, and because two read
    /// groups pooled into one group share this answer.
    GroupPutNoReadHere { slippage_group: u32 },
    /// **The run is not what it claims.** This run's slippage fit never named this read group,
    /// so there is no slippage group to look its numbers up under — a library present at calling
    /// time that the pre-pass did not know existed.
    UnknownReadGroup,
    /// **The run is not what it claims.** The read group's slippage group is past the end of the
    /// fit's own rows, so the map this type was built with names more groups than the fit was run
    /// over. Not a quiet library: the map and the fit came from different runs.
    GroupNotInTheFit {
        slippage_group: u32,
        groups_fitted: usize,
    },
}

/// The one type in this module whose fields are private, in a module of its own so that
/// **nothing else here can build one without the check**.
///
/// **The nesting is load-bearing rather than tidy**, and it is the same device
/// `calling::genotype_prior` uses for the five types whose invariants matter: a private field is
/// visible to a module's *descendants*, so a type declared directly in `stratum_fits` could be
/// built field by field — skipping the constructor — from anywhere in this file, its test module
/// included. One level of nesting makes those siblings instead, and the literal fails with
/// `error[E0451]`.
///
/// **It is a struct rather than the enum it reads as, and that is why.** An enum's variant
/// fields carry the enum's own visibility — there is no private field in a public variant — so
/// three variants would have left every check optional.
///
/// **What the checks stop was measured rather than imagined.** Built by hand with an empty
/// `weights`, [`AlleleLengthVariation::allele_span`] computed `(0 - 1) / 2`: a subtract-with-overflow
/// panic in debug, and in **release** a span of `-1`, against which `offset.abs() <= -1` is never
/// true — so every candidate at every tract of that stratum took the shape floor and the prior
/// came back `[1e-12, 1e-12, …]`, degenerate and silent. An even class count did the same more
/// quietly: `[0.1, 0.2, 0.3, 0.4]` gives a span of 1, which puts the *second* class at the
/// reference offset and leaves the fourth unreachable.
mod checked_spectrum {
    use super::AlleleLengthVariationRung;

    /// **Which of the three rungs that carry shares an allele length variation came from.**
    ///
    /// One value fewer than [`AlleleLengthVariationRung`], which also has the flat bottom rung,
    /// so that [`AlleleLengthVariation::fitted`] has no unreachable arm for the rung that has no
    /// weights.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum FittedFrom {
        /// The stratum's own tracts: its shares and its concentration.
        ThisStratum,
        /// **The stratum's own shares, held at the stated concentration**, because the fit
        /// reported no concentration for it (`doc/devel/ng/spec/str_fit_to_calling.md` §3.2).
        ThisStratumsSharesStatedSpread,
        /// Every tract of the stratum's motif period, pooled.
        ItsPeriodsPooledTracts,
    }

    /// **What a repeat tract's genotype prior is seeded from: a shape and a strength**, and
    /// which rung of the tract ladder they came from
    /// (`doc/devel/ng/spec/population_diversity.md` §4.4).
    ///
    /// The shape is a **allele length variation** — how a stratum's chromosomes are spread over tract
    /// lengths, one share per whole-repeat offset from the **reference** tract length,
    /// `-span ..= +span`. It is not the ordinary-site path's **frequency spectrum**, which is how
    /// allele frequencies are spread across the population; the two are separate quantities on
    /// separate paths and this project keeps the two words apart in code as in prose
    /// (`population_diversity.md` §2).
    ///
    /// **Every rung answers**: a run cannot get *no shape* at a tract, only a shape from further
    /// down the ladder, and which rung it was travels with the numbers so that whoever carries
    /// the rung into the run's output can tell a call resting on a measurement from one resting
    /// on a stated constant.
    ///
    /// **⛔ There used to be a second spectrum to cross this with, and there is not any more.**
    /// `population_diversity.md` §8's fifth check asked that a run's *allele-frequency* spectrum
    /// and a tract's *length* spectrum be separate types with no conversion between them, and a
    /// `compile_fail` doctest here handed the first where the second belonged. The first is
    /// deleted: the SNP/indel seed is two integrals of the fitted population curve and evaluates
    /// nothing into allele-count classes
    /// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5).
    ///
    /// **The doctest went with it rather than being left standing**, because it would still have
    /// failed to compile — for the wrong reason. A `compile_fail` test that passes because the
    /// type it names no longer exists proves nothing about the type it was written to protect,
    /// and this repository has shipped that shape of test before.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct AlleleLengthVariation<'a> {
        /// **`None` on the bottom rung and `Some` on the other three**, which is the invariant the
        /// two constructors hold and the reason the field is private.
        weights: Option<&'a [f64]>,
        concentration: f64,
        rung: AlleleLengthVariationRung,
    }

    impl<'a> AlleleLengthVariation<'a> {
        /// Wrap a fitted shape and the strength it is held with, at one of the three rungs that
        /// carry shares.
        ///
        /// `weights` sums to one; `concentration` is the Dirichlet total those shares are the
        /// mean of.
        ///
        /// # Panics
        ///
        /// **On a class count that is not `2 · span + 1` for some span of at least one**, and on
        /// a concentration that is not finite and strictly positive. Both are structural and
        /// cost two comparisons at one lookup a locus, against a prior that is wrong at every
        /// tract of the stratum — see this module's own documentation for what each of them
        /// produced when they were absent.
        ///
        /// **That the shares sum to one is *not* checked here**, and deliberately: it is
        /// `O(classes)` where these are `O(1)`, and it is already checked once per run where the
        /// spectra are gathered — [`StratumFits::over`](super::StratumFits::over) and
        /// [`of_gathered_rows`](super::StratumFits::of_gathered_rows), both
        /// through `checked_allele_length_variation`. A caller that builds one of these from weights
        /// that did not come through the gather owes that check itself: the seed's total is
        /// `concentration × (mass the candidates cover)`, which is a claim about conviction only
        /// if the mass is a share.
        #[must_use]
        pub fn fitted(weights: &'a [f64], concentration: f64, from: FittedFrom) -> Self {
            assert!(
                weights.len() >= 3 && weights.len() % 2 == 1,
                "an allele length variation runs from -span to +span in whole repeat units, so its class \
                 count is odd and at least three; got {}",
                weights.len()
            );
            assert!(
                concentration.is_finite() && concentration > 0.0,
                "an allele length variation is held with a finite, strictly positive number of \
                 chromosomes' worth of belief; got {concentration}"
            );
            Self {
                weights: Some(weights),
                concentration,
                rung: match from {
                    FittedFrom::ThisStratum => AlleleLengthVariationRung::StratumsOwnFit,
                    FittedFrom::ThisStratumsSharesStatedSpread => {
                        AlleleLengthVariationRung::StratumsOwnSharesStatedSpread
                    }
                    FittedFrom::ItsPeriodsPooledTracts => {
                        AlleleLengthVariationRung::PeriodsPooledTracts
                    }
                },
            }
        }

        /// The ladder's bottom rung: a flat shape over whatever lengths the locus offers, at a
        /// stated strength.
        ///
        /// # Panics
        ///
        /// On a concentration that is not finite and strictly positive — a Dirichlet with a
        /// total of zero has no mean for a shape to be.
        #[must_use]
        pub fn stated_flat(concentration: f64) -> Self {
            assert!(
                concentration.is_finite() && concentration > 0.0,
                "the stated-flat rung states a finite, strictly positive number of chromosomes' \
                 worth of belief; got {concentration}"
            );
            Self {
                weights: None,
                concentration,
                rung: AlleleLengthVariationRung::StatedFlat,
            }
        }

        /// The fitted shares, or nothing on the bottom rung — where the shape is flat over the
        /// locus's own candidate lengths and there is no vector to hand out.
        ///
        /// **The three rungs that carry shares are deliberately one answer here.** A consumer
        /// builds the same seed from any of them; what differs is the provenance, which
        /// [`Self::rung`] carries.
        #[inline]
        #[must_use]
        pub fn fitted_weights(&self) -> Option<&'a [f64]> {
            self.weights
        }

        /// The Dirichlet total: how many chromosomes' worth of belief the shape is held with.
        #[inline]
        #[must_use]
        pub fn concentration(&self) -> f64 {
            self.concentration
        }

        /// How far either side of the reference tract length the fitted shares reach, in whole
        /// repeat units — `None` on the bottom rung, which reaches nowhere in particular.
        ///
        /// The two outermost classes hold every allele this far or further, so the seed builder
        /// shares an end class's weight among every candidate at or past it
        /// (`str_allele_window.md` §1.2).
        ///
        /// **The subtraction cannot underflow**, because [`Self::fitted`] is the only door onto
        /// a `Some` and it refuses a class count below three.
        #[inline]
        #[must_use]
        pub fn allele_span(&self) -> Option<i32> {
            self.weights.map(|weights| ((weights.len() - 1) / 2) as i32)
        }

        /// Which rung of the tract ladder this came from — the value a run carries into its
        /// output beside the call.
        #[inline]
        #[must_use]
        pub fn rung(&self) -> AlleleLengthVariationRung {
            self.rung
        }
    }
}

pub use checked_spectrum::{AlleleLengthVariation, FittedFrom};

/// Which rung of the tract ladder a tract's prior shape came from — the value that **belongs**
/// in the run's output beside the call (`doc/devel/ng/spec/population_diversity.md` §4.4 for the
/// ladder, §1's third goal for why it has to travel).
///
/// **The calling loop's driver carries it there**: it looks the spectrum up once per tract and
/// puts the rung on the locus's
/// [`LocusInference`](crate::calling::LocusInference). A SNP/indel locus carries none,
/// because its prior comes from a frequency spectrum instead.
///
/// **Separate from [`AlleleLengthVariation`] because it outlives it.** The spectrum borrows the run's
/// frozen parameters and dies with the locus; the rung is what an output carries.
///
/// **The order the variants are declared in is their rank, strongest first.** The enum derives
/// `Ord`, and a repeat-tract bundle reports the largest rung among its runs as the weakest thing
/// its prior rests on ([`BundlePrior::rung`](crate::calling::genotype_prior::seed_ssr::BundlePrior::rung)).
/// A new rung goes where it falls in the ladder, not at the end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum AlleleLengthVariationRung {
    /// Fitted here, from this stratum's own tracts.
    StratumsOwnFit,
    /// **This stratum's own shares, held at the stated concentration**, because the fit placed
    /// no concentration for it — at one sample, nothing in the data can
    /// (`doc/devel/ng/spec/str_fit_to_calling.md` §3.2, decided by the owner 2026-10-14).
    ///
    /// **Ranked above the period's pool** because the class's own shares were measured to predict
    /// held-out tracts better than a flat shape, and the pool is no longer produced, so its order
    /// against this rung was never measured (same section).
    StratumsOwnSharesStatedSpread,
    /// Borrowed from the stratum's motif period, pooled over every tract of it.
    PeriodsPooledTracts,
    /// Defaulted: no fit for this period at all.
    StatedFlat,
}

/// One fitted allele length variation and the strength it is held with, owned.
///
/// Both a stratum's own fit and a period's pool produce exactly this, which is why they are one
/// type here. The lookup tells them apart, and tells a stratum with a concentration from one
/// without.
///
/// **A period's pool always carries a concentration**: the only door onto the pooled table,
/// [`StratumFits::of_gathered_rows`], refuses one without, because the ladder has no rung for a
/// pooled shape at a stated strength.
///
/// **Public, and its fields with it, so that [`StratumFits::of_gathered_rows`] has a named pair
/// to take** — the alternative was a bare `(Vec<f64>, f64)` in that signature, where the two
/// halves are told apart only by their types. Nothing is weakened by it: the checks that make a
/// spectrum a distribution live in the constructors that take one, not in this type, exactly as
/// they did while it was private.
#[derive(Debug, Clone, PartialEq)]
pub struct FittedAlleleLengthVariation {
    /// The share of this stratum's chromosomes at each whole-repeat offset from the reference
    /// length, running from `-span` to `+span` — an odd count, summing to one.
    pub weights: Vec<f64>,
    /// How many chromosomes' worth of belief the shares are held with — **`None` where the fit
    /// reported none**, as it does at one sample (`doc/devel/ng/spec/str_fit_to_calling.md`
    /// §3.2). Such a stratum is seeded from its own shares at the stated concentration.
    pub concentration: Option<f64>,
    /// The concentration's standard error, on the concentration's own scale; `None` where the
    /// fit gave it none, and always where the concentration itself is `None`
    /// (`str_fit_to_calling.md` §2.1).
    pub concentration_standard_error: Option<f64>,
}

/// **The concentration the tract ladder's bottom rung states where the run fitted no stratum
/// at all.**
///
/// Ten chromosomes' worth of belief, spread over whatever lengths the locus offers — a
/// between-locus share of `F = 1/(1 + 10) = 1/11`, the chance that two chromosomes drawn at a tract
/// carry the same length.
///
/// **Why ten and not one (owner, 2026-10-17).** The total sets how strongly the prior prefers a
/// homozygote: at a flat shape over `K` candidates it makes each homozygote `(c/K + 1)/(2c/K)`
/// times as likely as each heterozygote. At one chromosome — the value until 2026-10-17 — that is
/// 1.5 at two candidates and 2 at three, so every run that reaches this rung, including every
/// single sample, was told its genome is homozygous at three tracts in four. That is a homozygote
/// excess of F = 0.5 layered under the sample's own inbreeding coefficient, which is a separate
/// term of the prior (`calling_priors.md` §3.2) and the one that should carry it. On GIAB HG002,
/// HG003 and HG004 called one at a time with `--defaults`, the biallelic repeat tracts pvc calls
/// at 30x are 74 heterozygous, 49 homozygous for the other length and 5 for the reference. Moving
/// from one to ten made 17 tract genotypes right and 1 wrong across 5x–30x, nearly all of them
/// near ties the prior had tipped to a homozygote, and moved nothing at 300x
/// (`doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`). At two candidates ten gives
/// 0.27 / 0.45 / 0.27; results were flat between ten and thirty.
///
/// **Measured on one human trio, so still soft — and only for a flat shape.** Ten is what makes
/// a *flat* shape behave: the shape is itself wrong at a real tract (it gives the reference
/// length no more weight than any other), and a large total keeps it from pushing near ties
/// to a homozygote. Over a stratum's own fitted shares it is the wrong number, which is why the
/// second rung has its own, [`STATED_OWN_SHARES_CONCENTRATION`].
///
/// **The shape it is spread over is the locus's candidate lengths, where
/// `population_diversity.md` §4.4's table says the *reachable* lengths** — every length the
/// stutter model can produce from a candidate, which is a strictly larger support. The candidate
/// set is what this builder is handed; the reachable lengths are built by the read likelihood
/// (`likelihood::ssr`'s `fill_reachable_lengths`) and are not in the prior's hands. Recorded as a
/// departure rather than taken silently: on the bottom rung the shape is flat either way, so the
/// two differ only in how the *total* is spread — over `K` candidates here against the larger
/// reachable set — and nothing has measured what that costs.
///
/// **It is soft, and naming it is what makes it movable**
/// (`doc/devel/ng/spec/population_diversity.md` §4.4). It is reached only where the run fitted
/// no stratum with a concentration: a run that has any takes its own median instead
/// ([`StratumFits::stated_concentration`]), which is §9's question 4 and its leaning.
pub const STATED_FLAT_CONCENTRATION: f64 = 10.0;

/// **The concentration the tract ladder's second rung holds a stratum's own fitted shares at,
/// where the run fitted no stratum's concentration** — one chromosome, a between-locus share of
/// `F = 1/(1 + 1) = 1/2`. This is every stratum of a one-sample fit: at one sample the fit measures
/// how a stratum's chromosomes spread over lengths but cannot separate how much of that spread lies
/// between loci from the sample's own homozygote excess (`str_fit_to_calling.md` §3.2).
///
/// **Measured, on HG002's assembly truth** (GIAB tandem-repeat benchmark v1.0.1, 50,000 regions):
/// a stratum's concentration is what makes the shares reproduce how often the genome's two copies
/// carry the same length, `A = (1 − P_same)/(P_same − Σ m²)` for shares `m`, and it comes out at
/// 0.9–1.6 over the homopolymer strata with 30 or more tracts and 1.0–2.5 over the dinucleotide
/// ones, while the one-sample fit's shares match the truth's to about 0.02 in most homopolymer
/// strata (`doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md` §8). Holding those shares
/// at the flat rung's ten instead pulled true homozygous-alternative tracts to heterozygous: at 5x,
/// correct genotypes fell from 2,758 to 2,487.
///
/// **One human genome, and an outbred one.** For an inbred sample the same `P_same` includes its
/// own homozygote excess, which is the sample's inbreeding coefficient's to carry, not this
/// number's; nothing here has been measured on one.
pub const STATED_OWN_SHARES_CONCENTRATION: f64 = 1.0;

/// How far a stratum's allele length variation may sum from one before it is refused as not a
/// distribution.
///
/// The fit renormalises after every coordinate move, so a real one is within a few units in the
/// last place; this leaves room for a caller that assembled a `StratumFit` by hand from rounded
/// numbers, and refuses raw counts.
const LENGTH_SPECTRUM_NORMALISATION_TOLERANCE: f64 = 1e-9;

/// One stratum's numbers, per slippage group — the three vectors a [`StratumOutcome`] exposes,
/// held together so a lookup indexes once.
#[derive(Debug, Clone, PartialEq)]
struct StratumRow {
    slippage: Vec<Option<Slippage>>,
    level: Vec<Option<LevelProvenance>>,
    shares: Vec<Option<SharesProvenance>>,
    own_fit_standard_errors: Vec<OwnFitStandardErrors>,
    /// Per group, as [`FittedSlippage::samples_fitted_on`] — one number for the stratum when the
    /// fit produced it, kept a group at a time so a file read back keeps what each row said.
    samples_fitted_on: Vec<Option<u64>>,
}

/// **One stratum's numbers as a parameters file holds them, with no fit behind them**, per
/// slippage group — what a test hands [`StratumFits::with_rows_of_a_file`]. The climb wrote such
/// rows for strata too thin to fit, taking every number from its period's curves.
#[cfg(test)]
#[derive(Debug, Clone)]
pub(crate) struct StratumRowsOfAFile {
    pub stratum: Stratum,
    /// Per slippage group, its three numbers — `None` where that group put no read here.
    pub slippage: Vec<Option<Slippage>>,
    /// Where each group's level came from; `Some` wherever the group has numbers.
    pub level_provenance: Vec<Option<LevelProvenance>>,
    /// Where each group's two shares came from.
    pub shares_provenance: Vec<Option<SharesProvenance>>,
}

/// **Every stratum's slippage numbers, indexed by the key the caller has.**
///
/// Built once per run from what the repeat-tract fit
/// ([`fit_repeat_tracts`](super::ssr_fit::penalised_fit::fit_repeat_tracts)) returned, or from a
/// parameters file's rows, then read
/// unchanged for the whole of calling: slippage is a frozen parameter, and nothing about it is
/// re-estimated per locus (`arch/calling_em_loop.md` §2).
#[derive(Debug, Clone, PartialEq)]
pub struct StratumFits {
    /// Which set of slippage numbers each read group's reads are drawn under — the run's own
    /// declaration, the same map [`gather_strata`](super::ssr_fit::gather_strata) was given.
    slippage_group_of: BTreeMap<ReadGroupId, u32>,
    by_stratum: BTreeMap<Stratum, StratumRow>,
    /// **Per stratum, the allele length variation and concentration its own fit produced** — the tract
    /// ladder's top rung, and its second rung where the fit placed no concentration; harvested
    /// from the same outcomes the slippage numbers come from.
    ///
    /// **Only a stratum fitted on its own tracts has an entry.** A stratum an older parameters
    /// file holds numbers for that the climb furnished whole from its period's slippage curves
    /// carries no allele length variation at all — nothing about it was estimated — and that
    /// absence is exactly what the pooled rung exists to answer.
    allele_length_variation_by_stratum: BTreeMap<Stratum, FittedAlleleLengthVariation>,
    /// **Per motif period, one allele length variation and concentration fitted over every tract of
    /// that period pooled** — the tract ladder's pooled rung. **No fit produces one now**: only
    /// an older parameters file, read by [`Self::of_gathered_rows`], fills it.
    allele_length_variation_by_period: BTreeMap<u8, FittedAlleleLengthVariation>,
    /// The strength the bottom rung states, over a flat shape: the run's own median fitted
    /// concentration, taken on the between-locus share over the strata with enough loci (see
    /// `median_concentration`), where any stratum has a concentration, and
    /// [`STATED_FLAT_CONCENTRATION`] where none has.
    stated_concentration: f64,
    /// The strength the second rung states, over a stratum's own fitted shares: the same median
    /// where there is one, and [`STATED_OWN_SHARES_CONCENTRATION`] where there is none — a
    /// different constant from the flat rung's, because the shape it holds is a different claim.
    stated_concentration_with_own_shares: f64,
    /// Where [`Self::stated_concentration_with_own_shares`] came from, carried as the flat rung's
    /// warrant is: a file may supply one of the two numbers and leave the other defaulted.
    stated_concentration_with_own_shares_warrant: Provenance,
    /// **Where that strength came from, carried rather than re-derived.**
    ///
    /// [`Provenance::FittedHere`] where the run took the median of its own strata's
    /// concentrations and [`Provenance::Defaulted`] where no stratum has one — and
    /// [`Provenance::Supplied`] where a parameters file handed the number over, which is the state
    /// that cannot be worked out from the numbers beside it.
    ///
    /// **It is a field because it is not a function of anything else this type holds.** Until
    /// 2026-08-30 the one caller that needed it — the parameters file's writer — re-derived it
    /// from whether `allele_length_variation_by_stratum` was empty, which is true of a run that fitted
    /// the median itself and false of a run handed one, and those two are exactly what a
    /// warrant is for. A file demoted to `supplied` (spec §2.1) and written back out re-emerged
    /// as `fitted_here`.
    stated_concentration_warrant: Provenance,
}

impl StratumFits {
    /// Gather what the fit produced.
    ///
    /// **A refused stratum contributes nothing rather than an empty row.** Its three per-group
    /// accessors return empty slices, so a row built from one would answer every read group with
    /// [`NoSlippage::GroupPutNoReadHere`] — which would say a library was silent where the truth
    /// is that the stratum has no answer for anybody. Leaving it out makes the lookup say
    /// [`NoSuchStratum`](NoSlippage::NoSuchStratum), which is the honest one of the two.
    ///
    /// # Panics
    ///
    /// **When two outcomes name one stratum**, which would otherwise lose one of them without a
    /// word — a map insert keeps the last and says nothing, and the two levels can differ by a
    /// factor of five.
    ///
    /// **The guarantee that they do not belongs to
    /// [`gather_strata`](super::ssr_fit::gather_strata), not to the fit.** The fit returns one
    /// outcome per *evidence* it was handed; what makes the strata distinct is that
    /// `gather_strata` keys its evidence off a map. The fit is public and the research bench
    /// builds `StratumEvidence` by hand, so a caller that assembled its own list could reach this.
    /// **Release-level, deliberately**: the cost is one comparison per stratum, of which a run has
    /// tens, and the alternative is a caller scoring every tract of one shape against another
    /// shape's polymerase.
    #[must_use]
    pub fn over(
        outcomes: &[StratumOutcome],
        slippage_group_of: BTreeMap<ReadGroupId, u32>,
    ) -> Self {
        let mut by_stratum = BTreeMap::new();
        let mut allele_length_variation_by_stratum = BTreeMap::new();
        // Each fitted class's concentration beside its count of loci with reads: the count is
        // what `median_concentration` chooses the classes it counts by.
        let mut concentration_and_loci_of_each_class = Vec::new();
        for outcome in outcomes {
            if matches!(outcome, StratumOutcome::Refused { .. }) {
                continue;
            }
            let stratum = outcome.stratum();
            // **The own fit's errors and its samples.**
            let (own_fit_errors, samples_fitted_on) = match outcome {
                StratumOutcome::Fitted(fit) => (
                    fit.standard_errors.as_ref(),
                    fit.samples_fitted_on
                        .map(|samples| u64::try_from(samples).expect("a cohort's size fits a u64")),
                ),
                StratumOutcome::Refused { .. } => (None, None),
            };
            let groups = outcome.slippage().len();
            let row = StratumRow {
                slippage: outcome.slippage().to_vec(),
                level: outcome.level_provenance().to_vec(),
                shares: outcome.shares_provenance().to_vec(),
                own_fit_standard_errors: (0..groups)
                    .map(|group| {
                        OwnFitStandardErrors::of_the_numbers_emitted(
                            own_fit_errors
                                .and_then(|errors| errors.slippage.get(group))
                                .and_then(Option::as_ref),
                            outcome
                                .level_provenance()
                                .get(group)
                                .and_then(Option::as_ref),
                            outcome
                                .shares_provenance()
                                .get(group)
                                .and_then(Option::as_ref),
                        )
                    })
                    .collect(),
                samples_fitted_on: vec![samples_fitted_on; groups],
            };
            // **Checked once here rather than at every lookup.** The fit builds the three
            // vectors from one mask in one pass, so they are the same length by
            // construction — but every field of `StratumFit` is public, so
            // a caller that assembled its own outcome could hand over a short one, and the
            // failure would then be an index panic inside a lookup rather than a sentence
            // naming the stratum. Build time is where a caller can act on it.
            assert!(
                row.slippage.len() == row.level.len() && row.slippage.len() == row.shares.len(),
                "the outcome for period {}, {} repeats holds {} slippage groups, {} level \
                 provenances and {} shares provenances, where the fit builds all three from one \
                 mask and they are always the same length",
                stratum.period,
                stratum.reference_repeats,
                row.slippage.len(),
                row.level.len(),
                row.shares.len(),
            );
            let displaced = by_stratum.insert(stratum, row);
            assert!(
                displaced.is_none(),
                "two of the fit's outcomes are for period {}, {} repeats, and one of them would \
                 be lost without a word — look at how the evidence handed to the fit was \
                 assembled, since `gather_strata` cannot produce a repeat",
                stratum.period,
                stratum.reference_repeats,
            );
            // **The stratum's own allele length variation, harvested from the same outcome.**
            // **A fitted stratum with no concentration keeps its shares**
            // (`str_fit_to_calling.md` §3.2); the lookup holds them at the stated concentration.
            if let StratumOutcome::Fitted(fit) = outcome {
                let checked = checked_allele_length_variation(
                    &fit.allele_length_variation,
                    fit.concentration,
                    fit.standard_errors
                        .as_ref()
                        .and_then(|errors| errors.concentration.value()),
                    &format_args!(
                        "the fit of period {}, {} repeats",
                        stratum.period, stratum.reference_repeats
                    ),
                );
                if let Some(concentration) = checked.concentration {
                    concentration_and_loci_of_each_class
                        .push((concentration, fit.tracts_of_its_own));
                }
                allele_length_variation_by_stratum.insert(stratum, checked);
            }
        }
        let median = median_concentration(&concentration_and_loci_of_each_class);
        let stated_warrant = if median.is_some() {
            Provenance::FittedHere
        } else {
            Provenance::Defaulted
        };
        Self {
            slippage_group_of,
            // **The fit's own warrant for it**: a median over this run's strata is a number the
            // run measured, and the flat constant where it measured none is a stated guess.
            // **It asks whether any stratum has a concentration, not whether any has shares**
            // (`str_fit_to_calling.md` §3.5): at one sample every stratum with reads has shares
            // and none has a concentration, and the constant would otherwise be written as
            // fitted, where the file's validator checks a `defaulted` value against the constant
            // and cannot catch a `fitted_here` one.
            //
            // **Both stated strengths are that median where there is one**; only the constants
            // they fall back to differ (`STATED_FLAT_CONCENTRATION`,
            // `STATED_OWN_SHARES_CONCENTRATION`).
            stated_concentration_warrant: stated_warrant,
            stated_concentration_with_own_shares_warrant: stated_warrant,
            by_stratum,
            allele_length_variation_by_stratum,
            allele_length_variation_by_period: BTreeMap::new(),
            stated_concentration: median.unwrap_or(STATED_FLAT_CONCENTRATION),
            stated_concentration_with_own_shares: median.unwrap_or(STATED_OWN_SHARES_CONCENTRATION),
        }
    }

    /// **The same numbers read back out of a run's parameters file**, row by row, rather than
    /// harvested from the fit's own outcomes.
    ///
    /// **Why not [`Self::over`].** That one takes [`StratumOutcome`]s — the fit's raw output,
    /// which carries a refusal and an allele length variation beside the numbers a lookup
    /// answers with. A parameters file carries the answers and not
    /// the fitting: a refused stratum is simply a stratum with no rows, and the median that seeds
    /// the ladder's bottom rung was computed when the fit ran and is written down
    /// ([`Self::stated_concentration`]) rather than recomputed here. Recomputing it would make a
    /// file that was hand-edited disagree with itself in a way no reader could see.
    ///
    /// **The warrant on that median is taken from the file too**, for the same reason and one
    /// more: spec §2.1 demotes every number of a file fitted under another census to `supplied`,
    /// and *supplied* is not a state anything about this run's strata could imply. A caller that
    /// worked the warrant out from whether any stratum was fitted would turn every demoted file
    /// back into a fitted one.
    ///
    /// `slippage_by_stratum` is one entry per stratum the file names, holding that stratum's
    /// slippage groups **densely from zero** — `None` where the pair has no row, which is spec
    /// §5's fifth state and never a zero slip rate.
    ///
    /// # Panics
    ///
    /// On a stratum whose row is empty, which is a stratum no lookup can answer for — the file
    /// writes no row for such a pair rather than an empty one. On an allele length variation that is not
    /// a distribution, and on a concentration that is not finite and positive: the same checks
    /// [`Self::over`] runs, for the same reason, and here they also stand behind
    /// `ParametersFile::validate`, which refuses the same shapes earlier and naming the key. On a
    /// period's pool with no concentration, which no rung of the ladder can hold.
    #[allow(
        clippy::too_many_arguments,
        reason = "a parameters file's repeat-tract section, row table by row table and each stated \
                  value with its warrant; the two stated values are passed apart because a file \
                  may supply one and leave the other defaulted"
    )]
    #[must_use]
    pub fn of_gathered_rows(
        slippage_group_of: BTreeMap<ReadGroupId, u32>,
        slippage_by_stratum: BTreeMap<Stratum, Vec<Option<FittedSlippage>>>,
        allele_length_variation_by_stratum: BTreeMap<Stratum, FittedAlleleLengthVariation>,
        allele_length_variation_by_period: BTreeMap<u8, FittedAlleleLengthVariation>,
        stated_concentration: f64,
        stated_concentration_warrant: Provenance,
        stated_concentration_with_own_shares: f64,
        stated_concentration_with_own_shares_warrant: Provenance,
    ) -> Self {
        let by_stratum = slippage_by_stratum
            .into_iter()
            .map(|(stratum, groups)| {
                assert!(
                    !groups.is_empty(),
                    "the row for period {}, {} repeats holds no slippage group at all; a stratum \
                     nothing was fitted for has no row rather than an empty one, and an empty \
                     one answers every read group with `GroupPutNoReadHere` — a library reported \
                     silent where the truth is that the stratum has no answer for anybody",
                    stratum.period,
                    stratum.reference_repeats
                );
                let mut row = StratumRow {
                    slippage: Vec::with_capacity(groups.len()),
                    level: Vec::with_capacity(groups.len()),
                    shares: Vec::with_capacity(groups.len()),
                    own_fit_standard_errors: Vec::with_capacity(groups.len()),
                    samples_fitted_on: Vec::with_capacity(groups.len()),
                };
                for fitted in groups {
                    row.slippage.push(fitted.map(|fitted| fitted.slippage));
                    row.level.push(fitted.map(|fitted| fitted.level));
                    row.shares.push(fitted.and_then(|fitted| fitted.shares));
                    row.own_fit_standard_errors.push(
                        fitted.map_or_else(OwnFitStandardErrors::default, |fitted| {
                            fitted.own_fit_standard_errors
                        }),
                    );
                    row.samples_fitted_on
                        .push(fitted.and_then(|fitted| fitted.samples_fitted_on));
                }
                (stratum, row)
            })
            .collect();
        fn checked<K: std::fmt::Debug + Ord>(
            spectra: BTreeMap<K, FittedAlleleLengthVariation>,
            what: &str,
        ) -> BTreeMap<K, FittedAlleleLengthVariation> {
            spectra
                .into_iter()
                .map(|(key, spectrum)| {
                    let checked = checked_allele_length_variation(
                        &spectrum.weights,
                        spectrum.concentration,
                        spectrum.concentration_standard_error,
                        &format_args!("the {what} allele length variation filed under {key:?}"),
                    );
                    (key, checked)
                })
                .collect()
        }
        // **A pool without a concentration is refused**: the ladder has a rung for a stratum's
        // own shares at the stated strength and none for a pool's (`str_fit_to_calling.md` §3.2
        // names only the first).
        if let Some(period) = allele_length_variation_by_period
            .iter()
            .find_map(|(period, pool)| pool.concentration.is_none().then_some(period))
        {
            panic!(
                "the period's pooled allele length variation filed under period {period} holds no \
                 concentration; only a stratum's own shares may be held at the stated one"
            );
        }
        Self {
            slippage_group_of,
            by_stratum,
            allele_length_variation_by_stratum: checked(
                allele_length_variation_by_stratum,
                "stratum's own",
            ),
            allele_length_variation_by_period: checked(
                allele_length_variation_by_period,
                "period's pooled",
            ),
            stated_concentration: checked_stated_concentration(stated_concentration),
            stated_concentration_warrant,
            stated_concentration_with_own_shares: checked_stated_concentration(
                stated_concentration_with_own_shares,
            ),
            stated_concentration_with_own_shares_warrant,
        }
    }

    /// **These fits plus rows only a parameters file holds, read the way the file's are**: strata
    /// whose numbers came whole from elsewhere, with no fit and no allele length variation of
    /// their own, and pooled allele length variations a motif period. The climb wrote both and
    /// no fit produces either now, so tests that need them build them here; they go through
    /// [`Self::of_gathered_rows`], the door an older file comes in by. The stated concentration
    /// and its warrant stay `self`'s, and `pools` replaces any pooled table `self` held.
    ///
    /// # Panics
    ///
    /// On a row for a stratum `self` already holds, as [`Self::over`] refuses two outcomes for
    /// one stratum; and on whatever [`Self::of_gathered_rows`] refuses.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_rows_of_a_file(
        mut self,
        rows: impl IntoIterator<Item = StratumRowsOfAFile>,
        pools: BTreeMap<u8, FittedAlleleLengthVariation>,
    ) -> Self {
        let slippage_by_stratum = rows
            .into_iter()
            .map(|row| {
                let groups = row
                    .slippage
                    .iter()
                    .zip(&row.level_provenance)
                    .zip(&row.shares_provenance)
                    .map(|((slippage, level), shares)| {
                        slippage.map(|slippage| FittedSlippage {
                            slippage,
                            level: level.expect(
                                "a test row gives every group with numbers a level provenance",
                            ),
                            shares: *shares,
                            own_fit_standard_errors: OwnFitStandardErrors::default(),
                            samples_fitted_on: None,
                        })
                    })
                    .collect();
                (row.stratum, groups)
            })
            .collect();
        let read = Self::of_gathered_rows(
            self.slippage_group_of.clone(),
            slippage_by_stratum,
            BTreeMap::new(),
            pools,
            self.stated_concentration,
            self.stated_concentration_warrant,
            self.stated_concentration_with_own_shares,
            self.stated_concentration_with_own_shares_warrant,
        );
        for (stratum, row) in read.by_stratum {
            assert!(
                self.by_stratum.insert(stratum, row).is_none(),
                "two rows for period {}, {} repeats",
                stratum.period,
                stratum.reference_repeats
            );
        }
        self.allele_length_variation_by_period = read.allele_length_variation_by_period;
        self
    }

    /// **The shape and the strength a tract's genotype prior is seeded from**, and which rung of
    /// the tract ladder they came from (`doc/devel/ng/spec/population_diversity.md` §4.4).
    ///
    /// **Fill the stratum from the *tract*, not from the candidate — the opposite of
    /// [`Self::at`], and the two calls sit a screen apart in the same assembly.** Slippage is a
    /// property of the tract a read was copied from, so it is looked up per candidate allele.
    /// This is the population's belief about *which lengths this tract can be*, which is one
    /// question per locus: the spectrum is indexed by whole-repeat offset from the **reference**
    /// tract length, so the locus's own reference repeat count is what picks the stratum and
    /// what the offsets are measured from. Looking it up from a candidate would re-centre every
    /// candidate's own shape on itself and flatten the prior.
    ///
    /// **It always answers.** The four rungs are the stratum's own fit; the stratum's own shares
    /// at [`Self::stated_concentration`], where the fit reported no concentration for it; its
    /// motif period's pooled tracts; and a flat shape at the stated concentration. A tract can
    /// land on the last one, never on nothing.
    #[must_use]
    pub fn allele_length_variation_at(
        &self,
        period: u8,
        reference_repeats: u64,
    ) -> AlleleLengthVariation<'_> {
        let stratum = Stratum {
            period,
            reference_repeats,
        };
        if let Some(fitted) = self.allele_length_variation_by_stratum.get(&stratum) {
            let (concentration, from) = match fitted.concentration {
                Some(own) => (own, FittedFrom::ThisStratum),
                None => (
                    self.stated_concentration_with_own_shares,
                    FittedFrom::ThisStratumsSharesStatedSpread,
                ),
            };
            return AlleleLengthVariation::fitted(&fitted.weights, concentration, from);
        }
        if let Some(pooled) = self.allele_length_variation_by_period.get(&period) {
            return AlleleLengthVariation::fitted(
                &pooled.weights,
                pooled
                    .concentration
                    .expect("both doors onto the pooled table refuse a pool with no concentration"),
                FittedFrom::ItsPeriodsPooledTracts,
            );
        }
        AlleleLengthVariation::stated_flat(self.stated_concentration)
    }

    /// **The strength the tract ladder's bottom rung states, over a flat shape**; the second rung
    /// (a stratum's own shares) reads [`Self::stated_concentration_with_own_shares`]. It is
    /// the median of the run's strata's between-locus shares `F = 1 / (1 + c)`, converted back to
    /// a concentration `c`. The median counts the strata with a concentration and at least
    /// [`LOCI_FOR_AN_ERROR`] loci with reads, or every stratum with a concentration where none has
    /// that many (see `median_concentration`). It is [`STATED_FLAT_CONCENTRATION`] where no stratum
    /// has a concentration — which includes a one-sample run, where every stratum with reads is
    /// fitted but none carries a concentration (`str_fit_to_calling.md` §3.5).
    ///
    /// **The run's own median rather than a constant wherever there is one**, because how
    /// monomorphic a species' tracts are is a fact about the species and the run has measured it
    /// at every stratum that could be fitted — where a stated constant is a fact about nothing
    /// (`population_diversity.md` §9, question 4 and its leaning).
    ///
    /// **Taken on `F`, not on the concentration** — a choice the spec makes, and one that matters
    /// only at an even count. There the median is the mean of the two middle `F`s converted back,
    /// so concentrations 4 and 6 give 29/6 ≈ 4.83, not 5. At an odd count it is the middle
    /// stratum's own concentration, exactly as fitted.
    #[must_use]
    pub fn stated_concentration(&self) -> f64 {
        self.stated_concentration
    }

    /// **Where the bottom rung's stated strength came from** — see
    /// [`Self::stated_concentration`] for what the number is.
    ///
    /// **Carried, never inferred.** `Supplied` is the state that says a parameters file handed
    /// this number over, and nothing about a run's own strata distinguishes it from a median the
    /// run took itself.
    #[must_use]
    pub fn stated_concentration_warrant(&self) -> Provenance {
        self.stated_concentration_warrant
    }

    /// **The strength the tract ladder's second rung holds a stratum's own shares at**, where the
    /// fit gave the stratum shares and no concentration — every stratum of a one-sample fit. The
    /// run's own median where it has one, as for [`Self::stated_concentration`]; otherwise
    /// [`STATED_OWN_SHARES_CONCENTRATION`].
    #[must_use]
    pub fn stated_concentration_with_own_shares(&self) -> f64 {
        self.stated_concentration_with_own_shares
    }

    /// **Where the second rung's stated strength came from** — carried, never inferred, as
    /// [`Self::stated_concentration_warrant`] is.
    #[must_use]
    pub fn stated_concentration_with_own_shares_warrant(&self) -> Provenance {
        self.stated_concentration_with_own_shares_warrant
    }

    #[cfg(test)]
    /// How many strata carry an allele length variation of their own — the tract ladder's top rung,
    /// which is **not** [`Self::strata`]: a stratum an older parameters file holds numbers for
    /// that the climb furnished from its period's slippage curves carries slippage numbers and no
    /// spectrum, so the second count can be the larger of the two.
    #[must_use]
    pub fn strata_with_a_allele_length_variation(&self) -> usize {
        self.allele_length_variation_by_stratum.len()
    }

    #[cfg(test)]
    /// How many motif periods carry a pooled allele length variation — zero unless an older
    /// parameters file supplied pools.
    #[must_use]
    pub fn periods_with_a_pooled_allele_length_variation(&self) -> usize {
        self.allele_length_variation_by_period.len()
    }

    /// The numbers for one read group at one stratum.
    ///
    /// **Fill the stratum from the *candidate*, not from the tract.** A read's chance of slipping
    /// is a property of the tract it was copied from, and that is the candidate allele
    /// (`spec/read_likelihoods.md` §4.4): a candidate of 6 repeats and one of 12 at the same
    /// locus are drawn from different strata and slip at measurably different rates — slippage
    /// rises about 1.3-fold per repeat count over the measured range. **So the stutter parameters
    /// cannot be hoisted out of the candidate loop**, and a caller that looked one up per locus
    /// from the reference tract's own length would score every candidate there against one
    /// polymerase model.
    ///
    /// **The two numbers are taken by name rather than as a [`Stratum`], and that is the whole
    /// reason.** `Stratum`'s own field is `reference_repeats` — the right word on the fit's side
    /// of the seam, where a stratum is the bin a *reference* tract was sorted into so that
    /// tracts of one shape could be pooled. A caller handed that type has a `Stratum` for the
    /// locus already in its hand and would pass it, which is the wrong number and nothing would
    /// say so. Naming the argument `candidate_repeats` makes the mistake one somebody has to
    /// type on purpose. The bins are the same bins; only which length picks one differs.
    ///
    /// # Errors
    ///
    /// [`NoSlippage`], which names which of the four absences this is.
    pub fn at(
        &self,
        read_group: ReadGroupId,
        period: u8,
        candidate_repeats: u64,
    ) -> Result<FittedSlippage, NoSlippage> {
        let stratum = Stratum {
            period,
            reference_repeats: candidate_repeats,
        };
        let group = *self
            .slippage_group_of
            .get(&read_group)
            .ok_or(NoSlippage::UnknownReadGroup)?;
        let row = self
            .by_stratum
            .get(&stratum)
            .ok_or(NoSlippage::NoSuchStratum)?;
        let index = group as usize;
        // **A group past the end of the row is not a quiet library**, and answering as though it
        // were would hide the thing worth knowing: the map this type was built with names more
        // groups than the fit was run over, so the two were assembled from different runs. It is
        // the same class of fact as [`NoSlippage::UnknownReadGroup`] and gets its own answer.
        if index >= row.slippage.len() {
            return Err(NoSlippage::GroupNotInTheFit {
                slippage_group: group,
                groups_fitted: row.slippage.len(),
            });
        }
        // Indexed rather than `get`-ed from here on: `over` has checked that the three vectors
        // are the same length, and the bound above is that length.
        fitted_at(row, stratum, index).ok_or(NoSlippage::GroupPutNoReadHere {
            slippage_group: group,
        })
    }

    /// Which slippage group a read group's reads are drawn under, for a caller reporting what it
    /// looked up rather than looking one up.
    #[must_use]
    pub fn slippage_group_of(&self, read_group: ReadGroupId) -> Option<u32> {
        self.slippage_group_of.get(&read_group).copied()
    }

    /// **Every `(stratum × slippage group)` that has numbers**, in stratum order and then in
    /// slippage-group order.
    ///
    /// **The name says *with numbers* because this is not a dense walk.** A pair with no entry is
    /// a slippage group that put no read in that stratum, and it is skipped rather than yielded
    /// as an absence — the same claim [`Self::at`] makes through
    /// [`NoSlippage::GroupPutNoReadHere`], said by omission instead of by a variant. A caller
    /// that zipped this against a per-stratum table of its own would misalign it.
    ///
    /// **For a consumer that must walk every cell rather than look one up**, which is what
    /// writing the run's parameters down is (`doc/devel/ng/spec/parameters_file.md` §3.7).
    /// [`Self::at`] answers one `(read group, candidate length)` question and cannot enumerate.
    ///
    /// # Panics
    ///
    /// On a slippage group that has numbers at a stratum and no level provenance beside them,
    /// naming the stratum and the group. [`Self::over`] checks that the three per-group vectors
    /// are the same *length* and not that they agree cell by cell, and every field of the fit's
    /// outcome types is public.
    pub fn each_stratum_and_group_with_numbers(
        &self,
    ) -> impl Iterator<Item = (Stratum, u32, FittedSlippage)> {
        self.by_stratum.iter().flat_map(|(stratum, row)| {
            (0..row.slippage.len()).filter_map(move |group| {
                let fitted = fitted_at(row, *stratum, group)?;
                let group = u32::try_from(group).expect("a fit's slippage groups fit in a u32");
                Some((*stratum, group, fitted))
            })
        })
    }

    /// **Every stratum that was fitted on its own tracts: its shares by whole-repeat offset, the
    /// concentration they are held with, and that concentration's error**, as the fit reported
    /// them. The concentration and its error are `None` where the fit placed none; such a stratum
    /// is still here, because its shares are its own.
    ///
    /// **A stratum absent here is data and not a hole**: one an older file's climb furnished from
    /// its period's slippage curves was never estimated, which is what the pooled rung exists to
    /// answer. The
    /// name says *fitted* for that reason — this is not one item a stratum.
    ///
    /// **The record is handed out as stored rather than wrapped in an [`AlleleLengthVariation`]**,
    /// which holds an absent concentration already replaced by the stated one: a writer handed
    /// the wrapper would write the stated number as though the stratum had fitted it.
    pub fn fitted_allele_length_variation_of_each_stratum(
        &self,
    ) -> impl Iterator<Item = (Stratum, &FittedAlleleLengthVariation)> {
        self.allele_length_variation_by_stratum
            .iter()
            .map(|(stratum, fitted)| (*stratum, fitted))
    }

    /// **Every motif period the run pooled an allele length variation over** — the tract ladder's
    /// pooled rung, and empty unless the run read an older parameters file that carries pools.
    /// Each yields its shares and the concentration they are
    /// held with, which a pool always has.
    ///
    /// The shares are handed out borrowed rather than wrapped in an [`AlleleLengthVariation`],
    /// whose `fitted_weights` is an `Option` this iterator can never leave empty, so a consumer
    /// would carry an unwrap that cannot fire.
    pub fn pooled_allele_length_variation_of_each_period(
        &self,
    ) -> impl Iterator<Item = (u8, &[f64], f64)> {
        self.allele_length_variation_by_period
            .iter()
            .map(|(period, pool)| {
                (
                    *period,
                    pool.weights.as_slice(),
                    pool.concentration.expect(
                        "both doors onto the pooled table refuse a pool with no concentration",
                    ),
                )
            })
    }

    #[cfg(test)]
    /// How many strata carry an answer — what a run summary reports.
    ///
    /// **It cannot tell a cohort with no repeat tracts from one where every stratum was
    /// refused**, and it is not meant to: both are zero, and both mean the caller gets no
    /// slippage anywhere. A run that needs to tell them apart reads the refusals off the
    /// outcomes, which carry their own reason.
    #[must_use]
    pub fn strata(&self) -> usize {
        self.by_stratum.len()
    }
}

/// The bottom rung's concentration, refused where it is not a Dirichlet total.
///
/// **`over` cannot reach this** — its median is either the middle stratum's own concentration,
/// which [`checked_allele_length_variation`] has already passed, or, at an even count, a ratio of
/// sums and products of two such concentrations that is positive whenever both are (see
/// [`median_concentration`]) — but [`StratumFits::of_gathered_rows`] takes
/// the number straight from a file, where it is one line a person can edit and where a run that
/// fitted no stratum at all still has one. A zero has no mean to be the shape of, and every
/// tract of every stratum that falls to the bottom rung would take it.
fn checked_stated_concentration(concentration: f64) -> f64 {
    assert!(
        concentration.is_finite() && concentration > 0.0,
        "the tract ladder's bottom rung states {concentration} chromosomes' worth of belief; a \
         Dirichlet total is finite and strictly positive, and every tract of every stratum with \
         no fit of its own and no pool for its period is seeded from it"
    );
    concentration
}

/// Check one fitted allele length variation and its concentration, and own them.
///
/// **Checked where they are gathered rather than where they are read**, exactly as the three
/// slippage vectors' lengths are, and for the same reason: every field of
/// [`StratumFit`](super::ssr_fit::StratumFit) and [`FittedAlleleLengthVariation`] is public, and a
/// parameters file is text a person can edit, so a caller that assembled one by hand can hand over raw counts or a negative share. What that costs if it is
/// not caught here is not a crash: a seed built from weights that do not sum to one is a prior
/// that is wrong at **every tract of that stratum**, by an amount nothing downstream can see.
///
/// `what` names the fit in the message, because a run holds tens of these and the number alone
/// says nothing about which.
///
/// # Panics
///
/// On a class count that is not `2 · span + 1` for some span of at least one; on a share that is
/// negative or not finite; on shares that do not sum to one within
/// [`LENGTH_SPECTRUM_NORMALISATION_TOLERANCE`]; on a present concentration that is not finite or
/// not strictly positive — a Dirichlet with a total of zero has no mean; and on a concentration
/// error with no concentration beside it.
///
/// **An absent concentration passes**: the fit reports none where nothing in the data can place
/// it, and the lookup then holds the stratum's shares at the stated concentration
/// (`doc/devel/ng/spec/str_fit_to_calling.md` §3.2). Only a present one is checked.
fn checked_allele_length_variation(
    weights: &[f64],
    concentration: Option<f64>,
    concentration_standard_error: Option<f64>,
    what: &std::fmt::Arguments<'_>,
) -> FittedAlleleLengthVariation {
    assert!(
        weights.len() >= 3 && weights.len() % 2 == 1,
        "{what} holds {} length classes; the spectrum runs from -span to +span in whole repeat \
         units, so the count is odd and at least three",
        weights.len()
    );
    if let Some((class, weight)) = weights
        .iter()
        .enumerate()
        .find(|(_, weight)| weight.is_nan() || **weight < 0.0)
    {
        panic!(
            "{what} gives the offset of {} repeats a share of {weight}; a share of the \
             stratum's chromosomes cannot be negative or NaN",
            class as i64 - (weights.len() as i64 - 1) / 2
        );
    }
    let total: f64 = weights.iter().sum();
    assert!(
        (total - 1.0).abs() <= LENGTH_SPECTRUM_NORMALISATION_TOLERANCE,
        "{what} has length shares totalling {total}, where a spectrum sums to one within \
         {LENGTH_SPECTRUM_NORMALISATION_TOLERANCE:e} — raw counts here would scale the tract \
         prior by the count itself"
    );
    match concentration {
        Some(concentration) => assert!(
            concentration.is_finite() && concentration > 0.0,
            "{what} was held with {concentration} chromosomes' worth of belief; a Dirichlet total \
             is finite and strictly positive, and a zero has no mean to be the shape of"
        ),
        None => assert!(
            concentration_standard_error.is_none(),
            "{what} carries a concentration error of {} and no concentration; an error is the \
             error of a number, and there is none here",
            concentration_standard_error.unwrap_or_default()
        ),
    }
    FittedAlleleLengthVariation {
        weights: weights.to_vec(),
        concentration,
        concentration_standard_error,
    }
}

/// The numbers at one already-bounds-checked `(row, slippage group)`, or **nothing where that
/// group put no read in the stratum**.
///
/// **One place rather than two**, because [`StratumFits::at`] and
/// [`StratumFits::each_stratum_and_group_with_numbers`] build the same value from the same three
/// vectors, and two copies would drift the first time [`FittedSlippage`] gains a field: the
/// lookup would carry it and the enumeration would not, with nothing to fail.
///
/// # Panics
///
/// On a group with slippage numbers and no level provenance beside them, **naming the stratum
/// and the group** — this is reached inside a walk over every stratum of a run, where a message
/// that named neither would leave the reader bisecting.
fn fitted_at(row: &StratumRow, stratum: Stratum, index: usize) -> Option<FittedSlippage> {
    let slippage = row.slippage[index]?;
    let level = row.level[index].unwrap_or_else(|| {
        panic!(
            "period {}, {} repeats, slippage group {index} has slippage numbers and no level \
             provenance beside them; `over` checks the three vectors are the same length and not \
             that they agree cell by cell, and every field of the fit's outcome types is public — \
             so look at how this outcome was assembled",
            stratum.period, stratum.reference_repeats
        )
    });
    Some(FittedSlippage {
        slippage,
        level,
        shares: row.shares[index],
        own_fit_standard_errors: row.own_fit_standard_errors[index],
        samples_fitted_on: row.samples_fitted_on[index],
    })
}

/// The median of the run's strata's between-locus shares, converted back to a concentration, or
/// [`STATED_FLAT_CONCENTRATION`] where no stratum has a concentration.
///
/// **What is averaged.** A stratum's concentration `c` says how many chromosomes' worth of belief
/// its allele length variation carries; its between-locus share `F = 1 / (1 + c)` is the same
/// number on a scale from 0 to 1, and is the one the median is taken on
/// (`str_fit_to_calling.md` §3.5). The median is converted back to a concentration once, by
/// `c = (1 − F) / F` — at an even count in a closed form that gives the same number without the
/// round trip's cancellation (see the comment in the body).
///
/// **Where the scale matters.** `F` falls as `c` rises, so the strata sort into the same middle
/// either way:
///
/// - **At an odd count** the median is the middle stratum's own `F`, and converting it back gives
///   that stratum's own concentration. That concentration is returned as it was fitted rather
///   than round-tripped through `F`, which would change its last digits and nothing else.
/// - **At an even count** the median is the mean of the middle two `F`s, which is not `F` of the
///   mean of the middle two concentrations. Concentrations 4 and 6 have shares 1/5 and 1/7; their
///   mean, 6/35, converts back to 29/6 ≈ 4.83, where the mean of the concentrations is 5.
///
/// **Which strata count** (`str_fit_to_calling.md` §3.5, owner 2026-10-15). `fitted` holds one
/// pair per stratum that has a concentration: the concentration, and the stratum's count of loci
/// with at least one spanning read. A stratum with no concentration is not in it: it has nothing
/// to contribute, and its shares are held at this very median.
///
/// - **The median is over the strata with at least [`LOCI_FOR_AN_ERROR`] loci with reads**, the
///   same count below which the fit withholds a stratum's errors on its F, its slip rate and its
///   allele-length shares. A thin stratum's concentration is mostly its motif length's curve,
///   sometimes extended past the last stratum with data, and would otherwise get the same vote as
///   a stratum of a thousand loci.
/// - **Where no stratum has that many, the median is over every stratum in `fitted`**, so a run
///   of thin strata still states a number it measured.
///
/// The values are finite and positive by [`checked_allele_length_variation`], so the sort is total
/// and every share lies strictly between 0 and 1.
fn median_concentration(fitted: &[(f64, usize)]) -> Option<f64> {
    let with_enough_loci: Vec<f64> = fitted
        .iter()
        .filter(|&&(_, loci_with_reads)| loci_with_reads >= LOCI_FOR_AN_ERROR)
        .map(|&(concentration, _)| concentration)
        .collect();
    let mut concentrations = if with_enough_loci.is_empty() {
        fitted
            .iter()
            .map(|&(concentration, _)| concentration)
            .collect()
    } else {
        with_enough_loci
    };
    if concentrations.is_empty() {
        return None;
    }
    concentrations.sort_by(f64::total_cmp);
    let middle = concentrations.len() / 2;
    Some(if concentrations.len() % 2 == 1 {
        concentrations[middle]
    } else {
        // The mean of the two middle shares, `F̄ = (1/(1 + c1) + 1/(1 + c2)) / 2`, converted back
        // by `(1 − F̄) / F̄`, is algebraically `(c1 + c2 + 2·c1·c2) / (2 + c1 + c2)`: the round
        // trip through the file's conversions, `between_locus_share_of` and `concentration_of`
        // (`calling/parameters_file/between_locus_share.rs`, whose test
        // `the_mean_of_two_shares_converted_back_matches_the_closed_form` checks that round trip
        // against its own copy of this closed form; it would not notice a change here). This
        // file's `the_stated_concentration_is_the_median_of_the_runs_own_fitted_ones` pins this
        // function to the same 29/6 that test expects.
        // They are not called here, because `parameter_estimation` imports nothing from
        // `calling`. The closed form is used because the round trip cancels: below about 1e-16
        // both shares round to 1.0 and `1 − F̄` is exactly zero, where the closed form stays
        // positive for any two positive concentrations.
        let (c1, c2) = (concentrations[middle - 1], concentrations[middle]);
        (c1 + c2 + 2.0 * c1 * c2) / (2.0 + c1 + c2)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::share_curve::ShareSource;
    use crate::parameter_estimation::joint::slippage_curve::{
        LevelSource, RiseShape, SlippageCurve,
    };
    use crate::parameter_estimation::joint::ssr_fit::{
        ShareProvenance, StratumFit, StratumRefusal,
    };

    fn stratum(period: u8, reference_repeats: u64) -> Stratum {
        Stratum {
            period,
            reference_repeats,
        }
    }

    fn slippage(level: f64) -> Slippage {
        Slippage {
            level,
            shorter_share: 0.83,
            fall_off: 0.25,
        }
    }

    /// A curve whose value at 10 repeats is 0.11 — the shape a level provenance carries, so a
    /// test can ask the curve itself what it would have said.
    fn a_curve() -> SlippageCurve {
        SlippageCurve {
            rise_shape: RiseShape::new(1.0).expect("a rise shape of one is on the grid"),
            intercept: 0.01,
            slope: 0.01,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.077,
            cells: 5,
        }
    }

    /// A level provenance carrying a distinguishable `slipped_reads`, so that a test can tell
    /// one slippage group's provenance from another's.
    ///
    /// **The two are told apart by a number rather than by the variant**, because the variant is
    /// what a mutation would most easily preserve: reading group 0's provenance for every group
    /// leaves every `source` right and every count wrong.
    fn from_the_cell(slipped_reads: f64) -> LevelProvenance {
        LevelProvenance {
            source: LevelSource::Cell,
            curve: None,
            reach: None,
            slipped_reads: Some(slipped_reads),
        }
    }

    /// A shares provenance with real content — **the shape the fit actually emits**, where the
    /// helper below used to hand out `None`. The climb's thin strata and its pooled fit (both
    /// deleted 2026-10-16) set a shares provenance wherever they set a slippage number, so a
    /// fixture without one is a shape no file holds.
    fn shares_from_a_curve(slipped_reads: f64) -> SharesProvenance {
        let from_a_curve = ShareProvenance {
            source: ShareSource::Curve,
            curve: None,
            reach: None,
        };
        SharesProvenance {
            slipped_reads: Some(slipped_reads),
            shorter_share: from_a_curve,
            fall_off: from_a_curve,
        }
    }

    /// A stratum nothing was fitted at, whose numbers came from elsewhere — the rows a parameters
    /// file holds for one, built by hand so the tests do not pay for a fit.
    ///
    /// Every group with a slippage number gets a level provenance and a shares provenance beside
    /// it, which is the invariant every fit holds.
    fn derived(
        at: Stratum,
        slippage: Vec<Option<Slippage>>,
        level: Vec<Option<LevelProvenance>>,
    ) -> StratumRowsOfAFile {
        let shares = slippage
            .iter()
            .zip(&level)
            .map(|(numbers, level)| {
                numbers
                    .and(level.map(|level| shares_from_a_curve(level.slipped_reads.unwrap_or(1.0))))
            })
            .collect();
        StratumRowsOfAFile {
            stratum: at,
            slippage,
            level_provenance: level,
            shares_provenance: shares,
        }
    }

    /// The fits a parameters file holding only these rows reads into.
    fn of_rows(
        rows: impl IntoIterator<Item = StratumRowsOfAFile>,
        slippage_group_of: BTreeMap<ReadGroupId, u32>,
    ) -> StratumFits {
        StratumFits::over(&[], slippage_group_of).with_rows_of_a_file(rows, BTreeMap::new())
    }

    /// **A stratum's own-fit errors are kept only for the numbers it emits as its own, and the
    /// samples it was fitted on travel with every group's numbers.** Group 0's level is the
    /// cell's and its fall-off the stratum's own, its shorter share a curve's: it keeps the
    /// level's and the fall-off's errors and not the shorter share's, though the own fit computed
    /// all three. Group 1's level is a blend and its shares a curve's: no error at all. A stratum
    /// derived from its period's curves, fitted on nothing, has neither errors nor a sample count.
    #[test]
    fn own_fit_errors_are_kept_only_for_the_numbers_the_stratum_emits_as_its_own() {
        use crate::parameter_estimation::joint::slippage_curve::CurveReach;
        use crate::parameter_estimation::joint::ssr_fit::{
            SlippageErrors, StratumError, StratumErrors,
        };
        let errors = |level: f64, shorter_share: f64, fall_off: f64| {
            Some(SlippageErrors {
                level: StratumError::Estimated(level),
                shorter_share: StratumError::Estimated(shorter_share),
                fall_off: StratumError::Estimated(fall_off),
            })
        };
        let mut shares_partly_own = shares_from_a_curve(400.0);
        shares_partly_own.fall_off = ShareProvenance {
            source: ShareSource::Stratum,
            curve: None,
            reach: None,
        };
        let blended = LevelProvenance {
            source: LevelSource::Blend { curve_weight: 0.4 },
            curve: Some(a_curve()),
            reach: Some(CurveReach::Inside),
            slipped_reads: Some(90.0),
        };
        let fitted_on_a_subset = StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: Some(StratumErrors {
                slippage: vec![errors(0.001, 0.02, 0.03), errors(0.004, 0.05, 0.06)],
                allele_length_variation: vec![StratumError::NoShare; 3],
                concentration: StratumError::Estimated(0.5),
            }),
            stratum: stratum(2, 10),
            slippage: vec![Some(slippage(0.1)), Some(slippage(0.2))],
            allele_length_variation: leaning_short(),
            concentration: Some(2.0),
            samples_fitted_on: Some(300),
            tracts_of_its_own: 40,
            reads_crossing: 400,
            level_provenance: vec![Some(from_the_cell(400.0)), Some(blended)],
            shares_provenance: vec![Some(shares_partly_own), Some(shares_from_a_curve(90.0))],
        }));
        let fits = StratumFits::over(
            &[fitted_on_a_subset],
            BTreeMap::from([(ReadGroupId(0), 0), (ReadGroupId(1), 1)]),
        )
        .with_rows_of_a_file(
            [derived(
                stratum(2, 12),
                vec![Some(slippage(0.3))],
                vec![Some(from_the_cell(5.0))],
            )],
            BTreeMap::new(),
        );
        let cells: Vec<_> = fits.each_stratum_and_group_with_numbers().collect();
        assert_eq!(cells.len(), 3);
        assert_eq!(
            cells[0].2.own_fit_standard_errors,
            OwnFitStandardErrors {
                level: Some(0.001),
                shorter_share: None,
                fall_off: Some(0.03),
            }
        );
        assert_eq!(
            cells[1].2.own_fit_standard_errors,
            OwnFitStandardErrors::default()
        );
        assert_eq!(
            (cells[0].2.samples_fitted_on, cells[1].2.samples_fitted_on),
            (Some(300), Some(300))
        );
        assert_eq!(cells[2].0, stratum(2, 12));
        assert_eq!(
            (
                cells[2].2.own_fit_standard_errors,
                cells[2].2.samples_fitted_on
            ),
            (OwnFitStandardErrors::default(), None)
        );
    }

    /// **Errors beside a point on a penalised curve are kept**, because the number was fitted
    /// with the stratum's own reads (`str_fit_to_calling.md` §2.3): the level and both shares are
    /// penalised points, and all three keep the errors the fit gave them.
    #[test]
    fn own_fit_errors_are_kept_beside_a_penalised_point() {
        use crate::parameter_estimation::joint::slippage_curve::{CurveReach, PenalisedPoint};
        use crate::parameter_estimation::joint::ssr_fit::{
            SlippageErrors, StratumError, StratumErrors,
        };
        let point = PenalisedPoint {
            strength: 12.5,
            eased_here: false,
        };
        let on_the_curve = ShareProvenance {
            source: ShareSource::PenalisedCurve(point),
            curve: None,
            reach: Some(CurveReach::Inside),
        };
        let fitted = StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: Some(StratumErrors {
                slippage: vec![Some(SlippageErrors {
                    level: StratumError::Estimated(0.002),
                    shorter_share: StratumError::Estimated(0.04),
                    fall_off: StratumError::Estimated(0.07),
                })],
                allele_length_variation: vec![StratumError::NoShare; 3],
                concentration: StratumError::NotPlaced,
            }),
            stratum: stratum(2, 10),
            slippage: vec![Some(slippage(0.1))],
            allele_length_variation: leaning_short(),
            concentration: None,
            samples_fitted_on: None,
            tracts_of_its_own: 40,
            reads_crossing: 400,
            level_provenance: vec![Some(LevelProvenance {
                source: LevelSource::PenalisedCurve(point),
                curve: None,
                reach: Some(CurveReach::Inside),
                slipped_reads: Some(400.0),
            })],
            shares_provenance: vec![Some(SharesProvenance {
                slipped_reads: Some(400.0),
                shorter_share: on_the_curve,
                fall_off: on_the_curve,
            })],
        }));
        let fits = StratumFits::over(&[fitted], one_group());
        let cells: Vec<_> = fits.each_stratum_and_group_with_numbers().collect();
        assert_eq!(cells.len(), 1);
        assert_eq!(
            cells[0].2.own_fit_standard_errors,
            OwnFitStandardErrors {
                level: Some(0.002),
                shorter_share: Some(0.04),
                fall_off: Some(0.07),
            }
        );
    }

    fn one_group() -> BTreeMap<ReadGroupId, u32> {
        BTreeMap::from([(ReadGroupId(0), 0)])
    }

    /// **The numbers come back keyed by the pair the caller has** — the read group that produced
    /// the read, and the candidate's own motif period and repeat count.
    #[test]
    fn a_lookup_answers_with_the_stratum_and_groups_own_numbers() {
        let fits = of_rows(
            [
                derived(
                    stratum(2, 9),
                    vec![Some(slippage(0.08)), Some(slippage(0.11))],
                    vec![Some(from_the_cell(400.0)), Some(from_the_cell(9_999.0))],
                ),
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.10)), Some(slippage(0.13))],
                    vec![Some(from_the_cell(401.0)), Some(from_the_cell(9_998.0))],
                ),
                derived(
                    stratum(3, 9),
                    vec![Some(slippage(0.02)), Some(slippage(0.03))],
                    vec![Some(from_the_cell(402.0)), Some(from_the_cell(9_997.0))],
                ),
            ],
            BTreeMap::from([(ReadGroupId(7), 0), (ReadGroupId(9), 1)]),
        );
        let level_of = |group, period, repeats| {
            fits.at(group, period, repeats)
                .expect("the stratum has numbers")
                .slippage
                .level
        };

        // Both halves of the key are load-bearing: the same read group at two strata, and the
        // same stratum for two read groups, all four differ.
        assert_eq!(level_of(ReadGroupId(7), 2, 9), 0.08);
        assert_eq!(
            level_of(ReadGroupId(7), 2, 10),
            0.10,
            "the same library at one repeat count more",
        );
        assert_eq!(
            level_of(ReadGroupId(7), 3, 9),
            0.02,
            "the same library at the same repeat count of a longer motif",
        );
        assert_eq!(
            level_of(ReadGroupId(9), 2, 9),
            0.11,
            "the other library, which the run put in its own slippage group",
        );

        // **Each group's provenance is its own**, told apart by a count rather than by a
        // variant: reading group 0's provenance for every group would leave every `source`
        // right and every number wrong.
        let answer = fits
            .at(ReadGroupId(9), 2, 9)
            .expect("the second group has numbers");
        assert_eq!(answer.level.slipped_reads, Some(9_999.0));
        assert_eq!(
            answer
                .shares
                .expect("a group with numbers has a shares provenance")
                .slipped_reads,
            Some(9_999.0),
            "and so is its shares provenance",
        );
    }

    /// **Two candidates at one locus with different repeat counts get different numbers** — the
    /// case `spec/read_likelihoods.md` §4.4 names, and the reason the stutter parameters cannot
    /// be hoisted out of the candidate loop.
    ///
    /// A caller that keyed on the tract's *reference* length would score both against one
    /// polymerase model, and at tomato dinucleotides that is the difference between about 6 %
    /// of reads slipping and about 15 %.
    #[test]
    fn two_candidates_at_one_tract_are_two_strata() {
        let fits = of_rows(
            [
                derived(
                    stratum(2, 6),
                    vec![Some(slippage(0.06))],
                    vec![Some(from_the_cell(400.0))],
                ),
                derived(
                    stratum(2, 12),
                    vec![Some(slippage(0.15))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        );

        let six = fits
            .at(ReadGroupId(0), 2, 6)
            .expect("six repeats is fitted");
        let twelve = fits
            .at(ReadGroupId(0), 2, 12)
            .expect("twelve repeats is fitted");
        assert_eq!(six.slippage.level, 0.06);
        assert_eq!(twelve.slippage.level, 0.15);
        assert_ne!(
            six.slippage.level, twelve.slippage.level,
            "one tract, two candidates, two slippage rates",
        );
    }

    /// **Two read groups a run declares alike share one answer**, which is what a slippage group
    /// is for: a run that knows two libraries ran on one machine may pool them, and one that
    /// pools everything is saying it cannot tell them apart.
    #[test]
    fn read_groups_pooled_into_one_slippage_group_get_one_answer() {
        let fits = of_rows(
            [derived(
                stratum(2, 9),
                vec![Some(slippage(0.08))],
                vec![Some(from_the_cell(400.0))],
            )],
            BTreeMap::from([(ReadGroupId(3), 0), (ReadGroupId(4), 0)]),
        );

        // **Unwrapped on both sides**, because two absences also compare equal: written as a
        // comparison of two `Result`s this assertion passed even under a mutation that made
        // every lookup fail.
        let one = fits.at(ReadGroupId(3), 2, 9).expect("the first library");
        let other = fits.at(ReadGroupId(4), 2, 9).expect("the second library");
        assert_eq!(one, other);
        assert_eq!(one.slippage.level, 0.08);
        assert_eq!(fits.slippage_group_of(ReadGroupId(4)), Some(0));
        assert_eq!(
            fits.slippage_group_of(ReadGroupId(5)),
            None,
            "and a read group the run never named has no group at all",
        );
    }

    /// **The four ways a lookup can come back empty are four different answers**, and a caller
    /// that could not tell them apart would read a library the fit never saw as a library that
    /// was quiet here.
    #[test]
    fn the_four_absences_are_told_apart() {
        let fits = StratumFits::over(
            &[StratumOutcome::Refused {
                stratum: stratum(2, 20),
                tracts: 0,
                reason: StratumRefusal::NoSpanningReads,
            }],
            BTreeMap::from([
                (ReadGroupId(7), 0),
                (ReadGroupId(9), 1),
                (ReadGroupId(11), 4),
            ]),
        )
        .with_rows_of_a_file(
            [derived(
                stratum(2, 9),
                vec![Some(slippage(0.08)), None],
                vec![Some(from_the_cell(400.0)), None],
            )],
            BTreeMap::new(),
        );

        assert_eq!(
            fits.at(ReadGroupId(5), 2, 9),
            Err(NoSlippage::UnknownReadGroup),
            "a read group this run's fit never named",
        );
        assert_eq!(
            fits.at(ReadGroupId(7), 2, 11),
            Err(NoSlippage::NoSuchStratum),
            "a candidate repeat count no kept tract of that period occupies",
        );
        assert_eq!(
            fits.at(ReadGroupId(9), 2, 9),
            Err(NoSlippage::GroupPutNoReadHere { slippage_group: 1 }),
            "a stratum with an answer, and a library that put no read in it",
        );
        assert_eq!(
            fits.at(ReadGroupId(11), 2, 9),
            Err(NoSlippage::GroupNotInTheFit {
                slippage_group: 4,
                groups_fitted: 2,
            }),
            "a group past the end of the fit's own rows — the map and the fit disagree, which \
             is not a quiet library",
        );
        assert_eq!(
            fits.at(ReadGroupId(7), 2, 20),
            Err(NoSlippage::NoSuchStratum),
            "a refused stratum has no answer for anybody, which is not the same claim as one \
             library being silent — so it is left out rather than carried as an empty row",
        );
        assert_eq!(
            fits.strata(),
            1,
            "the refusal is not a stratum with numbers"
        );
    }

    /// **The lookup returns the level the fit emitted, not one it recomputed from the curve.**
    ///
    /// This is the property the module exists to hold, and the fixture is built so that it can
    /// fail. The stratum's emitted level is a blend of its own fit with its period's curve, as a
    /// file the climb wrote carries it, so it sits **between** the two and equals neither — a
    /// fixture whose level happened to equal the curve's value would make every recomputation
    /// look correct and this test would prove nothing.
    ///
    /// Measured on a small real fit when the climb still ran: putting the emitted level back
    /// through the blend moved the five blended strata by 0.6 % to 4.1 %, and left the
    /// `curve_weight` in their provenance unchanged — so the number moved and the provenance did
    /// not say so.
    #[test]
    fn the_level_is_the_one_the_fit_emitted_and_not_the_curves_own() {
        let curve = a_curve();
        let repeats = 10;
        // `a_curve` is a straight line in the level (rise shape 1): 0.01 + 0.01 × 10 = 0.11 at
        // 10 repeats, inside its fitted 8 to 12. A cell at 0.08 weighed against that 0.11 sits
        // between the two and equals neither.
        let the_curves_level = 0.11;
        let emitted = 0.084_711_129_860_078_91;
        assert!(
            0.08 < emitted && emitted < the_curves_level,
            "the fixture must not be a level the curve would also produce, or nothing it \
             asserts can fail",
        );
        let fits = of_rows(
            [derived(
                stratum(2, repeats),
                vec![Some(slippage(emitted))],
                vec![Some(LevelProvenance {
                    source: LevelSource::Blend {
                        curve_weight: 0.179_7,
                    },
                    curve: Some(curve),
                    reach: None,
                    slipped_reads: Some(400.0),
                })],
            )],
            one_group(),
        );

        let answer = fits
            .at(ReadGroupId(0), 2, repeats)
            .expect("the stratum has numbers");
        assert_eq!(
            answer.slippage.level, emitted,
            "the emitted level, not the curve's own value and not a second blend of the two",
        );
        assert_eq!(
            answer.level.curve,
            Some(curve),
            "and the curve travels beside it, so a consumer can see what it would have said",
        );
        assert_eq!(
            answer.slippage.shorter_share, 0.83,
            "the two shape numbers are the ones the row carried",
        );
        assert_eq!(answer.slippage.fall_off, 0.25);
    }

    /// **Two outcomes naming one stratum is refused rather than silently halved.** A map insert
    /// keeps the last and says nothing, and the two levels can differ by a factor of five.
    #[test]
    #[should_panic(expected = "two of the fit's outcomes are for period 2, 10 repeats")]
    fn two_outcomes_for_one_stratum_are_refused() {
        let _ = StratumFits::over(
            &[
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0),
                fitted(stratum(2, 10), 0.99, leaning_short(), 4.0),
            ],
            one_group(),
        );
    }

    // ---------------------------------------------------------------
    // The tract ladder: an allele length variation and a concentration per locus
    // ---------------------------------------------------------------

    /// A stratum fitted on its own tracts, carrying the allele length variation and concentration the
    /// fit produced — built by hand so the tests do not pay for a fit.
    ///
    /// **The slippage numbers and the spectrum are given separately and never derived from one
    /// another**, so a fixture can hold a stratum whose slippage says one thing and whose
    /// spectrum says another — which is what tells the two lookups apart.
    fn fitted(
        at: Stratum,
        level: f64,
        allele_length_variation: Vec<f64>,
        concentration: f64,
    ) -> StratumOutcome {
        StratumOutcome::Fitted(Box::new(StratumFit {
            standard_errors: None,
            stratum: at,
            slippage: vec![Some(slippage(level))],
            allele_length_variation,
            concentration: Some(concentration),
            samples_fitted_on: None,
            tracts_of_its_own: 40,
            reads_crossing: 400,
            level_provenance: vec![Some(from_the_cell(400.0))],
            shares_provenance: vec![Some(shares_from_a_curve(400.0))],
        }))
    }

    /// A three-class spectrum leaning short, so that a sign-flipped offset is a different
    /// number rather than the same one.
    fn leaning_short() -> Vec<f64> {
        vec![0.6, 0.3, 0.1]
    }

    /// A pooled period as an older parameters file holds it, built by hand.
    fn pool(weights: Vec<f64>, concentration: f64) -> FittedAlleleLengthVariation {
        FittedAlleleLengthVariation {
            weights,
            concentration: Some(concentration),
            concentration_standard_error: None,
        }
    }

    /// **The top rung.** A stratum fitted on its own tracts answers with its own spectrum and
    /// its own concentration, and says the fit was its own.
    #[test]
    fn a_stratum_fitted_here_answers_from_its_own_allele_length_variation() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        );

        let spectrum = fits.allele_length_variation_at(2, 10);
        assert_eq!(spectrum.rung(), AlleleLengthVariationRung::StratumsOwnFit);
        assert_eq!(spectrum.fitted_weights(), Some(&leaning_short()[..]));
        assert_eq!(spectrum.concentration(), 4.0);
        assert_eq!(spectrum.allele_span(), Some(1));
        assert_eq!(fits.strata_with_a_allele_length_variation(), 1);
    }

    /// **The pooled rung.** A stratum furnished from its period's slippage curves carries
    /// slippage numbers and no spectrum of its own, so it falls to the period's pooled tracts —
    /// and the two are told apart by their numbers, not only by the variant.
    #[test]
    fn a_stratum_with_no_fit_of_its_own_takes_its_periods_pooled_tracts() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        )
        .with_rows_of_a_file(
            [derived(
                stratum(2, 14),
                vec![Some(slippage(0.09))],
                vec![Some(from_the_cell(400.0))],
            )],
            BTreeMap::from([(2, pool(vec![0.2, 0.5, 0.3], 9.0))]),
        );

        // The derived stratum has slippage — it is not absent from the fit — and no spectrum.
        assert!(fits.at(ReadGroupId(0), 2, 14).is_ok());
        assert_eq!(fits.strata_with_a_allele_length_variation(), 1);
        assert_eq!(fits.periods_with_a_pooled_allele_length_variation(), 1);

        let spectrum = fits.allele_length_variation_at(2, 14);
        assert_eq!(
            spectrum.rung(),
            AlleleLengthVariationRung::PeriodsPooledTracts
        );
        assert_eq!(spectrum.fitted_weights(), Some(&[0.2, 0.5, 0.3][..]));
        assert_eq!(spectrum.concentration(), 9.0);

        // …and the stratum that has its own does **not** read the pool.
        let own = fits.allele_length_variation_at(2, 10);
        assert_eq!(own.rung(), AlleleLengthVariationRung::StratumsOwnFit);
        assert_eq!(own.concentration(), 4.0);
    }

    /// **A stratum that is not in the fit at all reads the pool too**, which is the case a
    /// caller meets most often: a locus whose reference tract length no stratum was fitted at.
    #[test]
    fn a_stratum_the_fit_never_saw_reads_its_periods_pool() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        )
        .with_rows_of_a_file([], BTreeMap::from([(2, pool(vec![0.2, 0.5, 0.3], 9.0))]));

        assert_eq!(
            fits.at(ReadGroupId(0), 2, 31),
            Err(NoSlippage::NoSuchStratum)
        );
        assert_eq!(
            fits.allele_length_variation_at(2, 31).rung(),
            AlleleLengthVariationRung::PeriodsPooledTracts
        );
    }

    /// **The bottom rung, and it is a different period's pool that must not be read.** A period
    /// with no pool of its own falls to the stated flat shape rather than borrowing the
    /// dinucleotides'.
    #[test]
    fn a_period_with_no_pool_falls_to_the_stated_flat_shape() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        )
        .with_rows_of_a_file([], BTreeMap::from([(2, pool(vec![0.2, 0.5, 0.3], 9.0))]));

        let spectrum = fits.allele_length_variation_at(3, 10);
        assert_eq!(spectrum.rung(), AlleleLengthVariationRung::StatedFlat);
        assert_eq!(spectrum.fitted_weights(), None);
        assert_eq!(spectrum.allele_span(), None);
        assert_eq!(
            spectrum.concentration(),
            4.0,
            "the run fitted one stratum, at a concentration of 4, so its median is 4"
        );
    }

    /// **A run that fitted no stratum states the constant**, which is the only place
    /// [`STATED_FLAT_CONCENTRATION`] is reached.
    #[test]
    fn a_run_that_fitted_no_stratum_states_the_constant() {
        let fits = StratumFits::over(&[], BTreeMap::new());
        let spectrum = fits.allele_length_variation_at(2, 10);
        assert_eq!(spectrum.rung(), AlleleLengthVariationRung::StatedFlat);
        assert_eq!(spectrum.concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
    }

    /// **The stated concentration is the run's own median, taken on the between-locus share
    /// `F = 1 / (1 + c)`**, at an odd count and at an even one — and the two arms are checked
    /// apart, because a mid-point taken with the wrong rounding is right at one of them and wrong
    /// at the other, and because only the even arm differs from a median taken on `c`.
    #[test]
    fn the_stated_concentration_is_the_median_of_the_runs_own_fitted_ones() {
        let three = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 30.0),
                fitted(stratum(2, 12), 0.05, leaning_short(), 6.0),
            ],
            one_group(),
        );
        // 6 is chosen because it does not survive the round trip: `(1 − F) / F` with
        // `F = 1 / (1 + 6)` is 6.000000000000001 in f64, so exact equality here fails if the
        // odd arm converts the middle share back instead of returning the stored concentration.
        assert_eq!(
            three.stated_concentration(),
            6.0,
            "1, 6 and 30 have shares 1/2, 1/7 and 1/31, whose median is 1/7, which is 6 — not the \
             mean, 12.33, and not the first inserted; the middle stratum's own 6 comes back \
             exactly, not round-tripped through its share"
        );
        assert_eq!(three.stated_concentration_warrant(), Provenance::FittedHere);

        let four = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 30.0),
                fitted(stratum(2, 12), 0.05, leaning_short(), 4.0),
                fitted(stratum(2, 14), 0.05, leaning_short(), 6.0),
            ],
            one_group(),
        );
        // Shares 1/2, 1/5, 1/7 and 1/31; the middle two are 1/5 and 1/7, whose mean is
        // (7 + 5) / 70 = 6/35, and (1 − 6/35) / (6/35) = 29/6. The median taken on the
        // concentrations themselves, (4 + 6) / 2 = 5, is the number this replaced.
        let expected = 29.0 / 6.0;
        assert!(
            (four.stated_concentration() - expected).abs() < 1e-12,
            "1, 4, 6 and 30 have a median share of 6/35, which is a concentration of 29/6 ≈ 4.833, \
             not (4 + 6) / 2 = 5; got {}",
            four.stated_concentration()
        );
        assert_eq!(four.stated_concentration_warrant(), Provenance::FittedHere);
    }

    /// A stratum fitted on its own tracts as [`fitted`] builds it, on `loci_with_reads` loci with
    /// at least one spanning read rather than `fitted`'s 40.
    fn fitted_on_loci(at: Stratum, concentration: f64, loci_with_reads: usize) -> StratumOutcome {
        let StratumOutcome::Fitted(mut fit) = fitted(at, 0.05, leaning_short(), concentration)
        else {
            unreachable!("`fitted` builds a fitted outcome")
        };
        fit.tracts_of_its_own = loci_with_reads;
        StratumOutcome::Fitted(fit)
    }

    /// **The median counts only the strata with at least [`LOCI_FOR_AN_ERROR`] loci with
    /// reads** (`str_fit_to_calling.md` §3.5, owner 2026-10-15). Three thick strata, at 8, 40
    /// and 40 loci, have concentrations 1, 6 and 30, a median of 6. Four thin ones, at 7, 3, 1 and
    /// 1 loci, have concentrations 0.1 to 0.4, between-locus shares near 1; counted, they would
    /// move the median of the seven to 0.4. The stratum at exactly 8 loci counts and the one at 7
    /// does not.
    #[test]
    fn the_stated_concentration_is_the_median_of_the_strata_with_enough_loci() {
        let fits = StratumFits::over(
            &[
                fitted_on_loci(stratum(1, 8), 1.0, LOCI_FOR_AN_ERROR),
                fitted_on_loci(stratum(1, 10), 6.0, 40),
                fitted_on_loci(stratum(2, 8), 30.0, 40),
                fitted_on_loci(stratum(2, 20), 0.1, LOCI_FOR_AN_ERROR - 1),
                fitted_on_loci(stratum(3, 20), 0.2, 3),
                fitted_on_loci(stratum(3, 30), 0.3, 1),
                fitted_on_loci(stratum(4, 20), 0.4, 1),
            ],
            one_group(),
        );
        assert_eq!(
            fits.stated_concentration(),
            6.0,
            "the three strata with 8 or more loci have concentrations 1, 6 and 30, a median of \
             6; over all seven strata it would be 0.4"
        );
        assert_eq!(fits.stated_concentration_warrant(), Provenance::FittedHere);
    }

    /// **Where no stratum has 8 loci with reads, the median is over every stratum with a
    /// concentration**, as it was before the threshold: a run of thin strata still states a
    /// number it measured.
    #[test]
    fn with_no_stratum_of_enough_loci_the_median_is_over_every_stratum_with_a_concentration() {
        let fits = StratumFits::over(
            &[
                fitted_on_loci(stratum(2, 8), 1.0, LOCI_FOR_AN_ERROR - 1),
                fitted_on_loci(stratum(2, 10), 30.0, 3),
                fitted_on_loci(stratum(2, 12), 6.0, 1),
            ],
            one_group(),
        );
        assert_eq!(fits.stated_concentration(), 6.0);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::FittedHere);
    }

    /// **A thick stratum counts only if it has a concentration.** Here the strata with 40 loci
    /// have shares and no concentration, so no stratum has both, and the median falls back to the
    /// thin stratum that has one.
    #[test]
    fn a_thick_stratum_without_a_concentration_does_not_count_towards_the_threshold() {
        let fits = StratumFits::over(
            &[
                fitted_without_a_concentration(stratum(2, 10), leaning_short()),
                fitted_without_a_concentration(stratum(2, 12), leaning_short()),
                fitted_on_loci(stratum(3, 8), 6.0, 2),
            ],
            one_group(),
        );
        assert_eq!(fits.stated_concentration(), 6.0);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::FittedHere);
    }

    /// **At one sample, thick and thin strata alike have no concentration**, so the stated one is
    /// the constant and its warrant `Defaulted`, whatever the loci counts.
    #[test]
    fn at_one_sample_the_threshold_leaves_the_stated_concentration_defaulted() {
        let StratumOutcome::Fitted(mut thin) =
            fitted_without_a_concentration(stratum(3, 8), leaning_short())
        else {
            unreachable!("`fitted_without_a_concentration` builds a fitted outcome")
        };
        thin.tracts_of_its_own = 2;
        let fits = StratumFits::over(
            &[
                fitted_without_a_concentration(stratum(2, 10), leaning_short()),
                StratumOutcome::Fitted(thin),
            ],
            one_group(),
        );
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::Defaulted);
    }

    /// **The even arm moves the median in the direction the share's curvature says**: at
    /// concentrations 1 and 9, shares 1/2 and 1/10, the mean share is 3/10 and the concentration
    /// (1 − 3/10) / (3/10) = 7/3 ≈ 2.33 — well below the mean of the concentrations, 5. Every
    /// high concentration has a share close to 0, so on the share's scale a stratum with little
    /// between-locus variation pulls the median less than it does on the concentration's.
    #[test]
    fn at_an_even_count_the_median_is_the_mean_of_the_middle_two_shares() {
        let two = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 9.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 1.0),
            ],
            one_group(),
        );
        let expected = 7.0 / 3.0;
        assert!(
            (two.stated_concentration() - expected).abs() < 1e-12,
            "shares 1/2 and 1/10 have a mean of 3/10, a concentration of 7/3; got {}",
            two.stated_concentration()
        );
    }

    /// **Two tiny middle concentrations give a tiny median, not zero.** At 1e-17 both shares
    /// round to exactly 1.0 in f64, so converting their mean back by `(1 − F) / F` returns 0.0,
    /// which the ladder's lower rungs refuse at the first tract that reads it. The closed form
    /// the median uses gives 1e-17.
    #[test]
    fn the_stated_concentration_is_positive_when_both_middle_concentrations_are_tiny() {
        let two = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1e-17),
                fitted(stratum(2, 10), 0.05, leaning_short(), 1e-17),
            ],
            one_group(),
        );
        let stated = two.stated_concentration();
        assert!(
            stated > 0.0 && ((stated - 1e-17) / 1e-17).abs() < 1e-12,
            "two concentrations of 1e-17 have a median of 1e-17; got {stated}"
        );
    }

    /// **At one sample the stated concentration is the constant, and says so.** The penalised fit
    /// reports shares and no concentration for every class with reads at one sample
    /// (`str_fit_to_calling.md` §3.5; arch `str_penalised_fit.md` §1.2), built here by hand as
    /// fitted outcomes with `concentration: None`: the classes keep their shares, the median has
    /// nothing to be taken over, and the warrant is `Defaulted`, not `FittedHere` — which it was
    /// while the warrant asked whether any class had shares.
    #[test]
    fn at_one_sample_no_class_has_a_concentration_and_the_stated_one_is_defaulted() {
        let fits = StratumFits::over(
            &[
                fitted_without_a_concentration(stratum(2, 10), leaning_short()),
                fitted_without_a_concentration(stratum(3, 8), vec![0.2, 0.5, 0.3]),
            ],
            one_group(),
        );
        assert_eq!(fits.strata_with_a_allele_length_variation(), 2);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::Defaulted);
        assert_eq!(
            fits.allele_length_variation_at(2, 10).rung(),
            AlleleLengthVariationRung::StratumsOwnSharesStatedSpread
        );
        assert_eq!(
            fits.allele_length_variation_at(3, 8).rung(),
            AlleleLengthVariationRung::StratumsOwnSharesStatedSpread
        );
    }

    /// **One class with a concentration is enough to make the stated one fitted**, however many
    /// beside it have none: the median is over that one, and it was measured.
    #[test]
    fn one_class_with_a_concentration_makes_the_stated_one_fitted() {
        let fits = StratumFits::over(
            &[
                fitted_without_a_concentration(stratum(2, 10), leaning_short()),
                fitted_without_a_concentration(stratum(2, 12), leaning_short()),
                fitted(stratum(3, 8), 0.05, leaning_short(), 6.0),
            ],
            one_group(),
        );
        assert_eq!(fits.stated_concentration(), 6.0);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::FittedHere);
    }

    /// A stratum fitted on its own tracts whose fit reported **its shares and no concentration**,
    /// as the fit does at one sample.
    fn fitted_without_a_concentration(
        at: Stratum,
        allele_length_variation: Vec<f64>,
    ) -> StratumOutcome {
        let StratumOutcome::Fitted(mut fit) = fitted(at, 0.05, allele_length_variation, 1.0) else {
            unreachable!("`fitted` builds a fitted outcome")
        };
        fit.concentration = None;
        StratumOutcome::Fitted(fit)
    }

    /// **The second rung, from the fit's outcomes.** A stratum with its own shares and no
    /// concentration keeps its shares and is held at the stated concentration, the median of the
    /// strata that have one; it does not fall to its period's pool, which ranks below it
    /// (`str_fit_to_calling.md` §3.2).
    #[test]
    fn a_stratum_with_shares_and_no_concentration_keeps_its_shares_at_the_stated_one() {
        let own_shares = vec![0.1, 0.2, 0.7];
        let fits = StratumFits::over(
            &[
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0),
                fitted_without_a_concentration(stratum(2, 12), own_shares.clone()),
            ],
            one_group(),
        )
        .with_rows_of_a_file([], BTreeMap::from([(2, pool(vec![0.2, 0.5, 0.3], 9.0))]));

        let spectrum = fits.allele_length_variation_at(2, 12);
        assert_eq!(
            spectrum.rung(),
            AlleleLengthVariationRung::StratumsOwnSharesStatedSpread
        );
        assert_eq!(spectrum.fitted_weights(), Some(&own_shares[..]));
        assert_eq!(
            spectrum.concentration(),
            4.0,
            "the stated concentration is the median of the one stratum that has one, 4; not the \
             pool's 9"
        );
        assert_eq!(fits.stated_concentration(), 4.0);
        assert_eq!(fits.strata_with_a_allele_length_variation(), 2);
        assert_eq!(
            fits.allele_length_variation_at(2, 10).rung(),
            AlleleLengthVariationRung::StratumsOwnFit
        );
        // The writer is handed the absence, not the stated number in its place.
        assert_eq!(
            fits.fitted_allele_length_variation_of_each_stratum()
                .map(|(at, fitted)| (at, fitted.concentration))
                .collect::<Vec<_>>(),
            vec![(stratum(2, 10), Some(4.0)), (stratum(2, 12), None)]
        );
    }

    /// **The second rung, from a hand-built table**: the same answer when the rows come from a
    /// parameters file instead of the fit, held at the stated concentration the file gives for that
    /// rung — not the flat rung's, which is a different number.
    #[test]
    fn a_hand_built_row_with_no_concentration_is_held_at_the_files_stated_one() {
        let own_shares = vec![0.25, 0.5, 0.25];
        let fits = StratumFits::of_gathered_rows(
            one_group(),
            BTreeMap::new(),
            BTreeMap::from([(
                stratum(3, 7),
                FittedAlleleLengthVariation {
                    weights: own_shares.clone(),
                    concentration: None,
                    concentration_standard_error: None,
                },
            )]),
            BTreeMap::from([(
                3,
                FittedAlleleLengthVariation {
                    weights: vec![0.2, 0.5, 0.3],
                    concentration: Some(9.0),
                    concentration_standard_error: None,
                },
            )]),
            7.0,
            Provenance::Defaulted,
            2.5,
            Provenance::Supplied,
        );

        let spectrum = fits.allele_length_variation_at(3, 7);
        assert_eq!(
            spectrum.rung(),
            AlleleLengthVariationRung::StratumsOwnSharesStatedSpread
        );
        assert_eq!(spectrum.fitted_weights(), Some(&own_shares[..]));
        assert_eq!(spectrum.concentration(), 2.5);
        assert_eq!(spectrum.allele_span(), Some(1));
        // A stratum of the same period with no row of its own still reads the pool.
        assert_eq!(
            fits.allele_length_variation_at(3, 8).rung(),
            AlleleLengthVariationRung::PeriodsPooledTracts
        );
    }

    /// **The rungs rank in ladder order**, strongest first, because a bundle reports the largest
    /// rung among its runs as the weakest thing its prior rests on.
    #[test]
    fn the_rungs_rank_in_ladder_order() {
        use AlleleLengthVariationRung::{
            PeriodsPooledTracts, StatedFlat, StratumsOwnFit, StratumsOwnSharesStatedSpread,
        };
        let mut rungs = [
            StatedFlat,
            PeriodsPooledTracts,
            StratumsOwnSharesStatedSpread,
            StratumsOwnFit,
        ];
        rungs.sort();
        assert_eq!(
            rungs,
            [
                StratumsOwnFit,
                StratumsOwnSharesStatedSpread,
                PeriodsPooledTracts,
                StatedFlat
            ]
        );
    }

    /// **The concentration's error travels from the fit to the table**, and only the estimated
    /// kind: a reason for having none is not a number.
    #[test]
    fn the_concentrations_error_is_taken_from_the_fit() {
        use crate::parameter_estimation::joint::ssr_fit::{StratumError, StratumErrors};
        let with_error = |concentration_error| {
            let StratumOutcome::Fitted(mut fit) =
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)
            else {
                unreachable!("`fitted` builds a fitted outcome")
            };
            fit.standard_errors = Some(StratumErrors {
                slippage: vec![None],
                allele_length_variation: vec![StratumError::NotPlaced; 3],
                concentration: concentration_error,
            });
            let fits = StratumFits::over(&[StratumOutcome::Fitted(fit)], one_group());
            let (_, fitted) = fits
                .fitted_allele_length_variation_of_each_stratum()
                .next()
                .expect("one stratum was fitted");
            fitted.concentration_standard_error
        };
        assert_eq!(with_error(StratumError::Estimated(0.4)), Some(0.4));
        assert_eq!(with_error(StratumError::NotPlaced), None);
    }

    /// An error with no concentration beside it is refused: there is no number for it to be the
    /// error of.
    #[test]
    #[should_panic(expected = "carries a concentration error of 0.3 and no concentration")]
    fn a_concentration_error_without_a_concentration_is_refused() {
        let _ = StratumFits::of_gathered_rows(
            one_group(),
            BTreeMap::new(),
            BTreeMap::from([(
                stratum(2, 10),
                FittedAlleleLengthVariation {
                    weights: leaning_short(),
                    concentration: None,
                    concentration_standard_error: Some(0.3),
                },
            )]),
            BTreeMap::new(),
            1.0,
            Provenance::Supplied,
            1.0,
            Provenance::Supplied,
        );
    }

    /// **A period's pool with no concentration is refused**: the ladder holds a stratum's own
    /// shares at the stated strength, and has no such rung for a pool.
    #[test]
    #[should_panic(expected = "filed under period 2 holds no concentration")]
    fn a_pool_with_no_concentration_is_refused() {
        let _ = StratumFits::of_gathered_rows(
            one_group(),
            BTreeMap::new(),
            BTreeMap::new(),
            BTreeMap::from([(
                2,
                FittedAlleleLengthVariation {
                    weights: leaning_short(),
                    concentration: None,
                    concentration_standard_error: None,
                },
            )]),
            1.0,
            Provenance::Supplied,
            1.0,
            Provenance::Supplied,
        );
    }

    /// **The two lookups are keyed by two different repeat counts, and this is the fixture that
    /// can tell them apart.** The tract sits at 10 repeats and carries a candidate at 14; the
    /// stratum at 10 and the stratum at 14 hold different slippage levels *and* different
    /// spectra. A slippage lookup keyed by the tract, or a spectrum lookup keyed by the
    /// candidate, changes both answers.
    #[test]
    fn slippage_is_keyed_by_the_candidate_and_the_allele_length_variation_by_the_tract() {
        let fits = StratumFits::over(
            &[
                fitted(stratum(2, 10), 0.05, vec![0.6, 0.3, 0.1], 4.0),
                fitted(stratum(2, 14), 0.20, vec![0.1, 0.3, 0.6], 25.0),
            ],
            one_group(),
        );

        // The candidate's own stratum answers the slippage question…
        let at_the_candidate = fits
            .at(ReadGroupId(0), 2, 14)
            .expect("the stratum at 14 repeats is in the fit");
        assert_eq!(at_the_candidate.slippage.level, 0.20);

        // …and the tract's own answers the prior question, at the same locus, differently.
        let at_the_tract = fits.allele_length_variation_at(2, 10);
        assert_eq!(at_the_tract.fitted_weights(), Some(&[0.6, 0.3, 0.1][..]));
        assert_eq!(at_the_tract.concentration(), 4.0);
    }

    /// **The counters count**, and no fixture that reads one had more than a single stratum or
    /// a single pool until this one — so `.len().min(1)` survived both.
    #[test]
    fn the_two_counters_report_how_many_of_each_the_run_holds() {
        let fits = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0),
            ],
            one_group(),
        )
        .with_rows_of_a_file(
            [derived(
                stratum(2, 14),
                vec![Some(slippage(0.09))],
                vec![Some(from_the_cell(400.0))],
            )],
            BTreeMap::from([
                (2, pool(vec![0.2, 0.5, 0.3], 9.0)),
                (3, pool(vec![0.3, 0.4, 0.3], 6.0)),
            ]),
        );

        assert_eq!(
            fits.strata_with_a_allele_length_variation(),
            2,
            "three strata, of which two were fitted here and one was furnished from curves"
        );
        assert_eq!(fits.periods_with_a_pooled_allele_length_variation(), 2);
        assert_eq!(
            fits.strata(),
            3,
            "every one of the three carries slippage, which is the count `strata` is about — it \
             is the larger of the two and they are not interchangeable"
        );
    }

    /// **A shape that misses being a distribution by a little is refused too**, which is what
    /// the tolerance is for — the only other fixture is off by a factor of ten, so a tolerance
    /// loosened from `1e-9` to `1e-2` would have changed nothing.
    ///
    /// A spectrum summing to 1.005 scales the tract prior by half a percent at every tract of
    /// that stratum, silently.
    #[test]
    #[should_panic(expected = "length shares totalling")]
    fn an_allele_length_variation_off_by_half_a_percent_is_refused() {
        let _ = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, vec![0.605, 0.3, 0.1], 4.0)],
            one_group(),
        );
    }

    /// …and one inside the tolerance is accepted, so the check is a tolerance rather than an
    /// equality that only a renormalised vector could pass.
    #[test]
    fn an_allele_length_variation_within_the_tolerance_is_accepted() {
        let fits = StratumFits::over(
            &[fitted(
                stratum(2, 10),
                0.05,
                vec![0.6 + 5e-10, 0.3, 0.1],
                4.0,
            )],
            one_group(),
        );
        assert_eq!(fits.strata_with_a_allele_length_variation(), 1);
    }

    /// **The lookup's own type refuses a shape it cannot index**, and it is a second door: the
    /// gather's checks run once per run, this one runs at every lookup, and only this one stands
    /// between a hand-built `AlleleLengthVariation` and the seed builder.
    ///
    /// An empty `weights` is the case that was silent before the constructor existed: the span
    /// is `(0 − 1) / 2`, which overflows in debug and is `-1` in release, and against `-1` no
    /// candidate is ever in reach.
    #[test]
    #[should_panic(expected = "class count is odd and at least three")]
    fn an_allele_length_variation_with_no_classes_cannot_be_built() {
        let _ = AlleleLengthVariation::fitted(&[], 4.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "class count is odd and at least three")]
    fn an_allele_length_variation_with_an_even_class_count_cannot_be_built() {
        let _ = AlleleLengthVariation::fitted(&[0.25; 4], 4.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "chromosomes' worth of belief; got 0")]
    fn an_allele_length_variation_held_with_no_conviction_cannot_be_built() {
        let _ = AlleleLengthVariation::fitted(&[0.25, 0.5, 0.25], 0.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "the stated-flat rung states a finite")]
    fn a_stated_flat_rung_with_no_conviction_cannot_be_built() {
        let _ = AlleleLengthVariation::stated_flat(f64::NAN);
    }

    /// A stratum's spectrum that does not sum to one is raw counts or a truncated copy, and
    /// seeding from it scales the prior by the total at every tract of that stratum.
    #[test]
    #[should_panic(expected = "the fit of period 3, 17 repeats has length shares totalling")]
    fn an_allele_length_variation_that_is_not_a_distribution_is_refused() {
        // Period 3 at 17 repeats, so that the message's two interpolations cannot be swapped
        // without the expectation failing — at period 2, 10 repeats either order reads
        // plausibly and neither number is the other's.
        let _ = StratumFits::over(
            &[fitted(stratum(3, 17), 0.05, vec![6.0, 3.0, 1.0], 4.0)],
            one_group(),
        );
    }

    /// A negative share is refused apart from the total, because `[-0.1, 0.6, 0.5]` sums to
    /// exactly one.
    #[test]
    #[should_panic(expected = "cannot be negative or NaN")]
    fn a_negative_length_share_is_refused() {
        let _ = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, vec![-0.1, 0.6, 0.5], 4.0)],
            one_group(),
        );
    }

    /// An even class count cannot be `2 · span + 1`, so the offsets it would be read at are not
    /// the offsets it was fitted at.
    ///
    /// **Four classes rather than two, and that is the whole point of the fixture.** `[0.5, 0.5]`
    /// fails the *count* check and the *even* check at once, so dropping `len() % 2 == 1` left it
    /// panicking on the other half and all nineteen tests green. At four, only the even check
    /// stands between the run and a span of `(4 − 1) / 2 = 1`, which puts the second class at the
    /// reference offset and leaves the fourth unreachable.
    #[test]
    #[should_panic(expected = "length classes; the spectrum runs from")]
    fn an_even_length_class_count_is_refused() {
        let _ = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, vec![0.25; 4], 4.0)],
            one_group(),
        );
    }

    /// A Dirichlet with a total of zero has no mean for the shape to be.
    #[test]
    #[should_panic(expected = "chromosomes' worth of belief")]
    fn a_concentration_of_zero_is_refused() {
        let _ = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 0.0)],
            one_group(),
        );
    }

    /// **The pool is checked the way a stratum's own fit is**, and by the same function — a
    /// second door onto the same seam is a second place for raw counts to get in.
    #[test]
    #[should_panic(
        expected = "the period's pooled allele length variation filed under 2 has length shares \
                    totalling"
    )]
    fn a_pooled_allele_length_variation_that_is_not_a_distribution_is_refused_too() {
        let _ = StratumFits::over(&[], BTreeMap::new())
            .with_rows_of_a_file([], BTreeMap::from([(2, pool(vec![6.0, 3.0, 1.0], 4.0))]));
    }

    /// **Each allele length variation comes back naming the rung it came off**, and the two rungs are
    /// different claims: a stratum's own tracts against every tract of its motif period pooled.
    ///
    /// Nothing else in the crate reads which rung these two iterators stamp — the parameters
    /// file writes the spectrum's numbers and places it by which table it is in — so without
    /// this, swapping the two rungs changes no test.
    #[test]
    fn each_allele_length_variation_names_the_rung_it_came_off() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        )
        .with_rows_of_a_file(
            [derived(
                stratum(2, 14),
                vec![Some(slippage(0.09))],
                vec![Some(from_the_cell(40.0))],
            )],
            BTreeMap::from([(2, pool(vec![0.2, 0.5, 0.3], 2.5))]),
        );

        assert_eq!(
            fits.fitted_allele_length_variation_of_each_stratum()
                .map(|(at, fitted)| (at, fitted.clone()))
                .collect::<Vec<_>>(),
            vec![(
                stratum(2, 10),
                FittedAlleleLengthVariation {
                    weights: leaning_short(),
                    concentration: Some(4.0),
                    concentration_standard_error: None,
                }
            )],
            "only the stratum fitted on its own tracts has a spectrum; the derived one carries \
             none by construction"
        );
        assert_eq!(
            fits.pooled_allele_length_variation_of_each_period()
                .map(|(period, shares, concentration)| (period, shares.to_vec(), concentration))
                .collect::<Vec<_>>(),
            vec![(2, vec![0.2, 0.5, 0.3], 2.5)]
        );
        assert_eq!(
            fits.allele_length_variation_at(2, 10).rung(),
            AlleleLengthVariationRung::StratumsOwnFit
        );
        assert_eq!(
            fits.allele_length_variation_at(2, 14).rung(),
            AlleleLengthVariationRung::PeriodsPooledTracts,
            "the derived stratum falls to its period's pool, which is the middle rung"
        );
    }

    /// **A slippage group that put no read in a stratum is skipped rather than yielded**, which
    /// is the same claim [`StratumFits::at`] makes through `GroupPutNoReadHere` — and the reason
    /// the iterator's name says *with numbers*.
    #[test]
    fn each_stratum_and_group_with_numbers_skips_a_pair_with_none() {
        let two_groups = BTreeMap::from([(ReadGroupId(0), 0), (ReadGroupId(1), 1)]);
        let fits = of_rows(
            [
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.05)), None],
                    vec![Some(from_the_cell(400.0)), None],
                ),
                derived(
                    stratum(2, 14),
                    vec![None, Some(slippage(0.09))],
                    vec![None, Some(from_the_cell(40.0))],
                ),
            ],
            two_groups,
        );

        assert_eq!(
            fits.each_stratum_and_group_with_numbers()
                .map(|(at, group, numbers)| (at, group, numbers.slippage.level))
                .collect::<Vec<_>>(),
            vec![(stratum(2, 10), 0, 0.05), (stratum(2, 14), 1, 0.09)],
            "four cells, two of them with numbers, and each on the other slippage group"
        );
        assert_eq!(
            fits.each_stratum_and_group_with_numbers()
                .map(|(_, _, numbers)| numbers.level.slipped_reads)
                .collect::<Vec<_>>(),
            vec![Some(400.0), Some(40.0)],
            "each pair carries its own group's provenance, not the first group's"
        );
    }

    /// **A slippage group with numbers and no level provenance beside them is refused naming the
    /// stratum**, because this fires inside a walk over every stratum of a run.
    #[test]
    #[should_panic(expected = "period 2, 10 repeats, slippage group 0 has slippage numbers")]
    fn a_group_with_numbers_and_no_level_provenance_is_refused_by_name() {
        let StratumOutcome::Fitted(mut fit) = fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)
        else {
            unreachable!("`fitted` builds a fitted outcome")
        };
        fit.level_provenance = vec![None];
        fit.shares_provenance = vec![None];
        let fits = StratumFits::over(&[StratumOutcome::Fitted(fit)], one_group());
        let _ = fits.each_stratum_and_group_with_numbers().count();
    }
}
