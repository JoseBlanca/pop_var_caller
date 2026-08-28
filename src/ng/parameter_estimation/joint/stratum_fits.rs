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
//! level on a [`StratumOutcome`] is already the emitted one, by one of three routes, and its
//! [`LevelProvenance`] says which:
//!
//! - a stratum fitted on its own tracts has had its level weighed against its period's curve —
//!   [`fit_strata`](super::ssr_fit::fit_strata) draws the curves after every stratum has its own
//!   answer, then applies [`blend_level`](super::slippage_curve::blend_level) in place;
//! - a stratum too thin to be fitted takes its period's curve whole, never through `blend_level`;
//! - a run with curves switched off keeps every cell's own answer.
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

use crate::ng::types::ReadGroupId;

use super::census::Stratum;
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
/// `weights`, [`LengthSpectrum::allele_span`] computed `(0 - 1) / 2`: a subtract-with-overflow
/// panic in debug, and in **release** a span of `-1`, against which `offset.abs() <= -1` is never
/// true — so every candidate at every tract of that stratum took the shape floor and the prior
/// came back `[1e-12, 1e-12, …]`, degenerate and silent. An even class count did the same more
/// quietly: `[0.1, 0.2, 0.3, 0.4]` gives a span of 1, which puts the *second* class at the
/// reference offset and leaves the fourth unreachable.
mod checked_spectrum {
    use super::LengthSpectrumRung;

    /// **Which of the two fitted rungs a length spectrum came from.**
    ///
    /// Two-valued where [`LengthSpectrumRung`] is three-valued, so that
    /// [`LengthSpectrum::fitted`] has no unreachable arm for the rung that has no weights.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum FittedFrom {
        /// The stratum's own tracts.
        ThisStratum,
        /// Every tract of the stratum's motif period, pooled.
        ItsPeriodsPooledTracts,
    }

    /// **What a repeat tract's genotype prior is seeded from: a shape and a strength**, and
    /// which rung of the tract ladder they came from
    /// (`doc/devel/ng/spec/population_diversity.md` §4.4).
    ///
    /// The shape is a **length spectrum** — how a stratum's chromosomes are spread over tract
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
    pub struct LengthSpectrum<'a> {
        /// **`None` on the bottom rung and `Some` on the other two**, which is the invariant the
        /// two constructors hold and the reason the field is private.
        weights: Option<&'a [f64]>,
        concentration: f64,
        rung: LengthSpectrumRung,
    }

    impl<'a> LengthSpectrum<'a> {
        /// Wrap a fitted shape and the strength it is held with, at one of the two fitted rungs.
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
        /// [`with_period_length_spectra`](super::StratumFits::with_period_length_spectra), both
        /// through `checked_length_spectrum`. A caller that builds one of these from weights
        /// that did not come through the gather owes that check itself: the seed's total is
        /// `concentration × (mass the candidates cover)`, which is a claim about conviction only
        /// if the mass is a share.
        #[must_use]
        pub fn fitted(weights: &'a [f64], concentration: f64, from: FittedFrom) -> Self {
            assert!(
                weights.len() >= 3 && weights.len() % 2 == 1,
                "a length spectrum runs from -span to +span in whole repeat units, so its class \
                 count is odd and at least three; got {}",
                weights.len()
            );
            assert!(
                concentration.is_finite() && concentration > 0.0,
                "a length spectrum is held with a finite, strictly positive number of \
                 chromosomes' worth of belief; got {concentration}"
            );
            Self {
                weights: Some(weights),
                concentration,
                rung: match from {
                    FittedFrom::ThisStratum => LengthSpectrumRung::StratumsOwnFit,
                    FittedFrom::ItsPeriodsPooledTracts => LengthSpectrumRung::PeriodsPooledTracts,
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
                rung: LengthSpectrumRung::StatedFlat,
            }
        }

        /// The fitted shares, or nothing on the bottom rung — where the shape is flat over the
        /// locus's own candidate lengths and there is no vector to hand out.
        ///
        /// **The two fitted rungs are deliberately one answer here.** A consumer builds the same
        /// seed from either; what differs is the provenance, which [`Self::rung`] carries.
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
        /// A candidate further from the reference than this is outside everything the fit ever
        /// saw, and the seed builder gives it a floor rather than the end class's weight.
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
        pub fn rung(&self) -> LengthSpectrumRung {
            self.rung
        }
    }
}

pub use checked_spectrum::{FittedFrom, LengthSpectrum};

/// Which rung of the tract ladder a tract's prior shape came from — the value that **belongs**
/// in the run's output beside the call (`doc/devel/ng/spec/population_diversity.md` §4.4 for the
/// ladder, §1's third goal for why it has to travel).
///
/// **The calling loop's driver carries it there**: it looks the spectrum up once per tract and
/// puts the rung on the locus's
/// [`LocusInference`](crate::ng::calling::LocusInference). A SNP/indel locus carries none,
/// because its prior comes from a frequency spectrum instead.
///
/// **Separate from [`LengthSpectrum`] because it outlives it.** The spectrum borrows the run's
/// frozen parameters and dies with the locus; the rung is what an output carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LengthSpectrumRung {
    /// Fitted here, from this stratum's own tracts.
    StratumsOwnFit,
    /// Borrowed from the stratum's motif period, pooled over every tract of it.
    PeriodsPooledTracts,
    /// Defaulted: no fit for this period at all.
    StatedFlat,
}

/// One fitted length spectrum and the strength it is held with, owned.
///
/// Both a stratum's own fit and a period's pool produce exactly this, which is why the two
/// rungs are one type here and two variants at the lookup.
#[derive(Debug, Clone, PartialEq)]
struct FittedLengthSpectrum {
    weights: Vec<f64>,
    concentration: f64,
}

/// **The concentration the tract ladder's bottom rung states where the run fitted no stratum
/// at all.**
///
/// One chromosome's worth of belief, spread flat over whatever lengths the locus offers — the
/// same quantity and the same reading [`ALPHA_REF`](crate::genetics::ALPHA_REF) carries on the
/// ordinary-site path, where it is the *count of chromosomes* a neutral seed gives the reference
/// allele. At one chromosome the reads move the prior from the first read onward, which is the
/// honest posture for a number with no measurement behind it.
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
/// no stratum whatever: a run that fitted any takes its own median instead
/// ([`StratumFits::stated_concentration`]), which is §9's question 4 and its leaning.
pub const STATED_FLAT_CONCENTRATION: f64 = 1.0;

/// How far a stratum's length spectrum may sum from one before it is refused as not a
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
}

/// **Every stratum's slippage numbers, indexed by the key the caller has.**
///
/// Built once per run from what [`fit_strata`](super::ssr_fit::fit_strata) returned, then read
/// unchanged for the whole of calling: slippage is a frozen parameter, and nothing about it is
/// re-estimated per locus (`arch/calling_em_loop.md` §2).
#[derive(Debug, Clone, PartialEq)]
pub struct StratumFits {
    /// Which set of slippage numbers each read group's reads are drawn under — the run's own
    /// declaration, the same map [`gather_strata`](super::ssr_fit::gather_strata) was given.
    slippage_group_of: BTreeMap<ReadGroupId, u32>,
    by_stratum: BTreeMap<Stratum, StratumRow>,
    /// **Per stratum, the length spectrum and concentration its own fit produced** — the tract
    /// ladder's top rung, harvested from the same outcomes the slippage numbers come from.
    ///
    /// **Only a stratum fitted on its own tracts has an entry.** One furnished from its
    /// period's slippage curves ([`DerivedStratum`](super::ssr_fit::DerivedStratum)) carries no
    /// length spectrum at all, by construction — nothing about it was estimated — and that
    /// absence is exactly what the middle rung exists to answer.
    length_spectrum_by_stratum: BTreeMap<Stratum, FittedLengthSpectrum>,
    /// **Per motif period, one length spectrum and concentration fitted over every tract of
    /// that period pooled** — the tract ladder's middle rung. Empty unless the run asked for
    /// it, which is what [`Self::with_period_length_spectra`] is.
    length_spectrum_by_period: BTreeMap<u8, FittedLengthSpectrum>,
    /// The strength the bottom rung states: the run's own median fitted concentration where any
    /// stratum was fitted, and [`STATED_FLAT_CONCENTRATION`] where none was.
    stated_concentration: f64,
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
    /// **A fitted stratum and a derived one are kept alike**, exactly as
    /// [`StratumOutcome::slippage`] keeps them: the read likelihood reads three numbers and the
    /// provenance beside them, and where the fit stopped and the curve started is in the
    /// provenance rather than in which door the numbers came through.
    ///
    /// # Panics
    ///
    /// **When two outcomes name one stratum**, which would otherwise lose one of them without a
    /// word — a map insert keeps the last and says nothing, and the two levels can differ by a
    /// factor of five.
    ///
    /// **The guarantee that they do not belongs to
    /// [`gather_strata`](super::ssr_fit::gather_strata), not to `fit_strata`.** `fit_strata`
    /// returns one outcome per *evidence* it was handed, and `derive_thin_strata` rewrites those
    /// in place rather than adding any; what makes the strata distinct is that `gather_strata`
    /// keys its evidence off a map. `fit_strata` is public and three examples and a benchmark
    /// build `StratumEvidence` by hand, so a caller that assembled its own list could reach this.
    /// **Release-level, deliberately**: the cost is one comparison per stratum, of which a run has
    /// tens, and the alternative is a caller scoring every tract of one shape against another
    /// shape's polymerase.
    #[must_use]
    pub fn over(
        outcomes: &[StratumOutcome],
        slippage_group_of: BTreeMap<ReadGroupId, u32>,
    ) -> Self {
        let mut by_stratum = BTreeMap::new();
        let mut length_spectrum_by_stratum = BTreeMap::new();
        for outcome in outcomes {
            if matches!(outcome, StratumOutcome::Refused { .. }) {
                continue;
            }
            let stratum = outcome.stratum();
            let row = StratumRow {
                slippage: outcome.slippage().to_vec(),
                level: outcome.level_provenance().to_vec(),
                shares: outcome.shares_provenance().to_vec(),
            };
            // **Checked once here rather than at every lookup.** Both of the fit's paths build
            // the three vectors from one mask in one pass, so they are the same length by
            // construction — but every field of `StratumFit` and `DerivedStratum` is public, so
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
                 be lost without a word — look at how the evidence handed to `fit_strata` was \
                 assembled, since `gather_strata` cannot produce a repeat",
                stratum.period,
                stratum.reference_repeats,
            );
            // **The stratum's own length spectrum, harvested from the same outcome**, and only
            // where the fit produced one: a stratum furnished from its period's slippage curves
            // carries none, which is the absence the middle rung answers.
            if let StratumOutcome::Fitted(fit) = outcome {
                length_spectrum_by_stratum.insert(
                    stratum,
                    checked_length_spectrum(
                        &fit.length_spectrum,
                        fit.concentration,
                        &format_args!(
                            "the fit of period {}, {} repeats",
                            stratum.period, stratum.reference_repeats
                        ),
                    ),
                );
            }
        }
        let stated_concentration = median_concentration(&length_spectrum_by_stratum);
        Self {
            slippage_group_of,
            by_stratum,
            length_spectrum_by_stratum,
            length_spectrum_by_period: BTreeMap::new(),
            stated_concentration,
        }
    }

    /// **Set the tract ladder's middle rung**: one pooled length spectrum and concentration a
    /// motif period, for the strata whose own fit does not exist.
    ///
    /// **It replaces rather than merges**, as a `with_` builder does: a second call keeps only
    /// the second call's periods.
    ///
    /// **A second call rather than an argument to [`Self::over`], because the pool costs a run a
    /// second pass over its tracts** ([`fit_period_length_spectra`](super::ssr_fit::fit_period_length_spectra)),
    /// where everything else this type carries is already computed. A run that skips it is not
    /// broken: [`Self::length_spectrum_at`] then answers from the ladder's bottom rung, a flat
    /// shape at a stated concentration, and says so — so the omission shows up in the run's own
    /// record as a rung rather than as a wrong number.
    ///
    /// **It does not move [`Self::stated_concentration`].** That median is over the strata's own
    /// fits, and a period's pool is fitted from the very same tracts — counting both would weigh
    /// one period's tracts twice.
    ///
    /// # Panics
    ///
    /// On a pool whose length spectrum is not a distribution, or whose concentration is not
    /// finite and positive — the same checks [`Self::over`] runs on a stratum's own, and for the
    /// same reason: every field of the fit's output types is public.
    ///
    /// On a pool filed under a motif period that is not its own. The period is carried twice —
    /// as the map's key and on [`PeriodLengthSpectrum`](super::ssr_fit::PeriodLengthSpectrum)
    /// itself — and the lookup reads the key, so a disagreement would file a period's tracts
    /// under another period's tracts with nothing to say so.
    #[must_use]
    pub fn with_period_length_spectra(
        mut self,
        pools: BTreeMap<u8, super::ssr_fit::PeriodLengthSpectrum>,
    ) -> Self {
        self.length_spectrum_by_period = pools
            .into_iter()
            .map(|(period, pool)| {
                assert_eq!(
                    period, pool.period,
                    "a pooled fit of motif period {} is filed under period {period}; the lookup \
                     reads the key, so every tract of period {period} would be seeded from \
                     period {}'s spread",
                    pool.period, pool.period
                );
                let checked = checked_length_spectrum(
                    &pool.length_spectrum,
                    pool.concentration,
                    &format_args!("the pooled fit of motif period {period}"),
                );
                (period, checked)
            })
            .collect();
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
    /// **It always answers.** The three rungs are the stratum's own fit, its motif period's
    /// pooled tracts, and a flat shape at [`Self::stated_concentration`]; a tract can land on
    /// the last one, never on nothing.
    #[must_use]
    pub fn length_spectrum_at(&self, period: u8, reference_repeats: u64) -> LengthSpectrum<'_> {
        let stratum = Stratum {
            period,
            reference_repeats,
        };
        if let Some(fitted) = self.length_spectrum_by_stratum.get(&stratum) {
            return LengthSpectrum::fitted(
                &fitted.weights,
                fitted.concentration,
                FittedFrom::ThisStratum,
            );
        }
        if let Some(pooled) = self.length_spectrum_by_period.get(&period) {
            return LengthSpectrum::fitted(
                &pooled.weights,
                pooled.concentration,
                FittedFrom::ItsPeriodsPooledTracts,
            );
        }
        LengthSpectrum::stated_flat(self.stated_concentration)
    }

    /// **The strength the tract ladder's bottom rung states**: the median of the concentrations
    /// this run's strata fitted, or [`STATED_FLAT_CONCENTRATION`] where the run fitted none.
    ///
    /// **The run's own median rather than a constant wherever there is one**, because how
    /// monomorphic a species' tracts are is a fact about the species and the run has measured it
    /// at every stratum that could be fitted — where a stated constant is a fact about nothing
    /// (`population_diversity.md` §9, question 4 and its leaning).
    ///
    /// **The median of an even count is the mean of the middle two**, which is the ordinary
    /// definition and not a choice this file is making.
    #[must_use]
    pub fn stated_concentration(&self) -> f64 {
        self.stated_concentration
    }

    /// How many strata carry a length spectrum of their own — the tract ladder's top rung,
    /// which is **not** [`Self::strata`]: a stratum furnished from its period's slippage curves
    /// carries slippage numbers and no spectrum, so the second count is the larger of the two.
    #[must_use]
    pub fn strata_with_a_length_spectrum(&self) -> usize {
        self.length_spectrum_by_stratum.len()
    }

    /// How many motif periods carry a pooled length spectrum — zero unless the run called
    /// [`Self::with_period_length_spectra`].
    #[must_use]
    pub fn periods_with_a_pooled_length_spectrum(&self) -> usize {
        self.length_spectrum_by_period.len()
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

    /// **Every stratum that was fitted on its own tracts: its shares by whole-repeat offset, and
    /// the strength they are held with** — the tract ladder's top rung.
    ///
    /// **A stratum absent here is data and not a hole**: one furnished from its period's
    /// slippage curves was never estimated, which is what the middle rung exists to answer. The
    /// name says *fitted* for that reason — this is not one item a stratum.
    ///
    /// **The shares are handed out borrowed rather than wrapped in a [`LengthSpectrum`]**, whose
    /// `fitted_weights` is an `Option` this iterator can never leave empty — so a consumer would
    /// carry an unwrap that cannot fire. A consumer that wants the wrapper builds one.
    pub fn fitted_length_spectrum_of_each_stratum(
        &self,
    ) -> impl Iterator<Item = (Stratum, &[f64], f64)> {
        self.length_spectrum_by_stratum
            .iter()
            .map(|(stratum, fitted)| (*stratum, fitted.weights.as_slice(), fitted.concentration))
    }

    /// **Every motif period the run pooled a length spectrum over** — the tract ladder's middle
    /// rung, and empty unless the run called [`Self::with_period_length_spectra`].
    ///
    /// Same shape and same reason as [`Self::fitted_length_spectrum_of_each_stratum`].
    pub fn pooled_length_spectrum_of_each_period(&self) -> impl Iterator<Item = (u8, &[f64], f64)> {
        self.length_spectrum_by_period
            .iter()
            .map(|(period, pool)| (*period, pool.weights.as_slice(), pool.concentration))
    }

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

/// Check one fitted length spectrum and its concentration, and own them.
///
/// **Checked where they are gathered rather than where they are read**, exactly as the three
/// slippage vectors' lengths are, and for the same reason: every field of
/// [`StratumFit`](super::ssr_fit::StratumFit) and
/// [`PeriodLengthSpectrum`](super::ssr_fit::PeriodLengthSpectrum) is public, so a caller that
/// assembled one by hand can hand over raw counts or a negative share. What that costs if it is
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
/// [`LENGTH_SPECTRUM_NORMALISATION_TOLERANCE`]; and on a concentration that is not finite and
/// strictly positive — a Dirichlet with a total of zero has no mean.
fn checked_length_spectrum(
    weights: &[f64],
    concentration: f64,
    what: &std::fmt::Arguments<'_>,
) -> FittedLengthSpectrum {
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
    assert!(
        concentration.is_finite() && concentration > 0.0,
        "{what} was held with {concentration} chromosomes' worth of belief; a Dirichlet total is \
         finite and strictly positive, and a zero has no mean to be the shape of"
    );
    FittedLengthSpectrum {
        weights: weights.to_vec(),
        concentration,
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
    })
}

/// The median of the concentrations the run's strata fitted, or [`STATED_FLAT_CONCENTRATION`]
/// where it fitted none.
///
/// **The mean of the middle two at an even count**, which is the ordinary definition. The values
/// are finite and positive by [`checked_length_spectrum`], so the sort is total.
fn median_concentration(fitted: &BTreeMap<Stratum, FittedLengthSpectrum>) -> f64 {
    if fitted.is_empty() {
        return STATED_FLAT_CONCENTRATION;
    }
    let mut totals: Vec<f64> = fitted.values().map(|fit| fit.concentration).collect();
    totals.sort_by(f64::total_cmp);
    let middle = totals.len() / 2;
    if totals.len() % 2 == 1 {
        totals[middle]
    } else {
        (totals[middle - 1] + totals[middle]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::joint::share_curve::ShareSource;
    use crate::ng::parameter_estimation::joint::slippage_curve::{
        LevelSource, RiseShape, SlippageCurve,
    };
    use crate::ng::parameter_estimation::joint::ssr_fit::{
        DerivedStratum, PeriodLengthSpectrum, ShareProvenance, StratumFit, StratumRefusal,
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
    /// helper below used to hand out `None`. `derive_thin_strata` and the pooled fit both set a
    /// shares provenance wherever they set a slippage number, so a fixture without one is a
    /// shape no run can produce.
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

    /// A stratum nothing was fitted at, whose numbers came from its period's curves — the shape
    /// [`StratumOutcome::Derived`] carries, built by hand so the tests do not pay for a fit.
    ///
    /// Every group with a slippage number gets a level provenance and a shares provenance beside
    /// it, which is the invariant both of the fit's paths hold.
    fn derived(
        at: Stratum,
        slippage: Vec<Option<Slippage>>,
        level: Vec<Option<LevelProvenance>>,
    ) -> StratumOutcome {
        let shares = slippage
            .iter()
            .zip(&level)
            .map(|(numbers, level)| {
                numbers
                    .and(level.map(|level| shares_from_a_curve(level.slipped_reads.unwrap_or(1.0))))
            })
            .collect();
        StratumOutcome::Derived(Box::new(DerivedStratum {
            stratum: at,
            slippage,
            level_provenance: level,
            shares_provenance: shares,
            tracts_of_its_own: 4,
            reads_crossing: 40,
        }))
    }

    fn one_group() -> BTreeMap<ReadGroupId, u32> {
        BTreeMap::from([(ReadGroupId(0), 0)])
    }

    /// **The numbers come back keyed by the pair the caller has** — the read group that produced
    /// the read, and the candidate's own motif period and repeat count.
    #[test]
    fn a_lookup_answers_with_the_stratum_and_groups_own_numbers() {
        let fits = StratumFits::over(
            &[
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
        let fits = StratumFits::over(
            &[
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
        let fits = StratumFits::over(
            &[derived(
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
            &[
                derived(
                    stratum(2, 9),
                    vec![Some(slippage(0.08)), None],
                    vec![Some(from_the_cell(400.0)), None],
                ),
                StratumOutcome::Refused {
                    stratum: stratum(2, 20),
                    tracts: 3,
                    reason: StratumRefusal::BelowTheFloor {
                        tracts: 3,
                        floor: 50,
                    },
                },
            ],
            BTreeMap::from([
                (ReadGroupId(7), 0),
                (ReadGroupId(9), 1),
                (ReadGroupId(11), 4),
            ]),
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
    /// fail. The stratum's emitted level is a blend of its own fit with its period's curve, so
    /// it sits **between** the two and equals neither — `assert_ne!` pins that, because a
    /// fixture whose level happened to equal the curve's value would make every recomputation
    /// look correct and this test would prove nothing.
    ///
    /// Measured on a small real fit: putting the emitted level back through `blend_level` moves
    /// the five blended strata by 0.6 % to 4.1 %, and leaves the `curve_weight` in their
    /// provenance unchanged — so the number moves and the provenance does not say so.
    #[test]
    fn the_level_is_the_one_the_fit_emitted_and_not_the_curves_own() {
        let curve = a_curve();
        let repeats = 10;
        // A cell at 0.08 weighed against this curve's 0.11: between the two, equal to neither.
        let emitted = 0.084_711_129_860_078_91;
        assert_ne!(
            emitted,
            curve.level_at(repeats),
            "the fixture must not be a level the curve would also produce, or nothing it \
             asserts can fail",
        );
        let fits = StratumFits::over(
            &[derived(
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
            answer
                .level
                .curve
                .expect("the curve is carried")
                .level_at(repeats),
            0.11,
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
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.05))],
                    vec![Some(from_the_cell(400.0))],
                ),
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.99))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        );
    }

    // ---------------------------------------------------------------
    // The tract ladder: a length spectrum and a concentration per locus
    // ---------------------------------------------------------------

    /// A stratum fitted on its own tracts, carrying the length spectrum and concentration the
    /// fit produced — built by hand so the tests do not pay for a climb.
    ///
    /// **The slippage numbers and the spectrum are given separately and never derived from one
    /// another**, so a fixture can hold a stratum whose slippage says one thing and whose
    /// spectrum says another — which is what tells the two lookups apart.
    fn fitted(
        at: Stratum,
        level: f64,
        length_spectrum: Vec<f64>,
        concentration: f64,
    ) -> StratumOutcome {
        StratumOutcome::Fitted(Box::new(StratumFit {
            stratum: at,
            slippage: vec![Some(slippage(level))],
            length_spectrum,
            concentration,
            log_likelihood_a_tract: -1.5,
            tracts_fitted: 40,
            borrowed: Vec::new(),
            converged: true,
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

    /// A pooled period, built by hand.
    fn pool(period: u8, length_spectrum: Vec<f64>, concentration: f64) -> PeriodLengthSpectrum {
        PeriodLengthSpectrum {
            period,
            length_spectrum,
            concentration,
            tracts_fitted: 900,
            strata_pooled: 5,
            converged: true,
        }
    }

    /// **The top rung.** A stratum fitted on its own tracts answers with its own spectrum and
    /// its own concentration, and says the fit was its own.
    #[test]
    fn a_stratum_fitted_here_answers_from_its_own_length_spectrum() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        );

        let spectrum = fits.length_spectrum_at(2, 10);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::StratumsOwnFit);
        assert_eq!(spectrum.fitted_weights(), Some(&leaning_short()[..]));
        assert_eq!(spectrum.concentration(), 4.0);
        assert_eq!(spectrum.allele_span(), Some(1));
        assert_eq!(fits.strata_with_a_length_spectrum(), 1);
    }

    /// **The middle rung.** A stratum furnished from its period's slippage curves carries
    /// slippage numbers and no spectrum of its own, so it falls to the period's pooled tracts —
    /// and the two are told apart by their numbers, not only by the variant.
    #[test]
    fn a_stratum_with_no_fit_of_its_own_takes_its_periods_pooled_tracts() {
        let fits = StratumFits::over(
            &[
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0),
                derived(
                    stratum(2, 14),
                    vec![Some(slippage(0.09))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        )
        .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![0.2, 0.5, 0.3], 9.0))]));

        // The derived stratum has slippage — it is not absent from the fit — and no spectrum.
        assert!(fits.at(ReadGroupId(0), 2, 14).is_ok());
        assert_eq!(fits.strata_with_a_length_spectrum(), 1);
        assert_eq!(fits.periods_with_a_pooled_length_spectrum(), 1);

        let spectrum = fits.length_spectrum_at(2, 14);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::PeriodsPooledTracts);
        assert_eq!(spectrum.fitted_weights(), Some(&[0.2, 0.5, 0.3][..]));
        assert_eq!(spectrum.concentration(), 9.0);

        // …and the stratum that has its own does **not** read the pool.
        let own = fits.length_spectrum_at(2, 10);
        assert_eq!(own.rung(), LengthSpectrumRung::StratumsOwnFit);
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
        .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![0.2, 0.5, 0.3], 9.0))]));

        assert_eq!(
            fits.at(ReadGroupId(0), 2, 31),
            Err(NoSlippage::NoSuchStratum)
        );
        assert_eq!(
            fits.length_spectrum_at(2, 31).rung(),
            LengthSpectrumRung::PeriodsPooledTracts
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
        .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![0.2, 0.5, 0.3], 9.0))]));

        let spectrum = fits.length_spectrum_at(3, 10);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::StatedFlat);
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
        let spectrum = fits.length_spectrum_at(2, 10);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::StatedFlat);
        assert_eq!(spectrum.concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
    }

    /// **The stated concentration is the run's own median**, at an odd count and at an even
    /// one — and the two arms are checked apart, because a mid-point taken with the wrong
    /// rounding is right at one of them and wrong at the other.
    #[test]
    fn the_stated_concentration_is_the_median_of_the_runs_own_fitted_ones() {
        let three = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 30.0),
                fitted(stratum(2, 12), 0.05, leaning_short(), 4.0),
            ],
            one_group(),
        );
        assert_eq!(
            three.stated_concentration(),
            4.0,
            "1, 4 and 30 have a median of 4 — not the mean, 11.67, and not the first inserted"
        );

        let four = StratumFits::over(
            &[
                fitted(stratum(2, 8), 0.05, leaning_short(), 1.0),
                fitted(stratum(2, 10), 0.05, leaning_short(), 30.0),
                fitted(stratum(2, 12), 0.05, leaning_short(), 4.0),
                fitted(stratum(2, 14), 0.05, leaning_short(), 6.0),
            ],
            one_group(),
        );
        assert_eq!(
            four.stated_concentration(),
            5.0,
            "1, 4, 6 and 30 have a median of (4 + 6) / 2"
        );
    }

    /// **A period's pool does not move the stated concentration**, because it is fitted from
    /// the very same tracts the strata's own fits read — counting both weighs one period twice.
    #[test]
    fn a_periods_pool_does_not_move_the_stated_concentration() {
        let fits = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        );
        let with_pool = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, leaning_short(), 4.0)],
            one_group(),
        )
        .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![0.2, 0.5, 0.3], 500.0))]));

        assert_eq!(fits.stated_concentration(), 4.0);
        assert_eq!(
            with_pool.stated_concentration(),
            4.0,
            "a pool at 500 chromosomes does not drag the run's median off its strata's own"
        );
    }

    /// **The two lookups are keyed by two different repeat counts, and this is the fixture that
    /// can tell them apart.** The tract sits at 10 repeats and carries a candidate at 14; the
    /// stratum at 10 and the stratum at 14 hold different slippage levels *and* different
    /// spectra. A slippage lookup keyed by the tract, or a spectrum lookup keyed by the
    /// candidate, changes both answers.
    #[test]
    fn slippage_is_keyed_by_the_candidate_and_the_length_spectrum_by_the_tract() {
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
        let at_the_tract = fits.length_spectrum_at(2, 10);
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
                derived(
                    stratum(2, 14),
                    vec![Some(slippage(0.09))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        )
        .with_period_length_spectra(BTreeMap::from([
            (2, pool(2, vec![0.2, 0.5, 0.3], 9.0)),
            (3, pool(3, vec![0.3, 0.4, 0.3], 6.0)),
        ]));

        assert_eq!(
            fits.strata_with_a_length_spectrum(),
            2,
            "three outcomes, of which two were fitted here and one was furnished from curves"
        );
        assert_eq!(fits.periods_with_a_pooled_length_spectrum(), 2);
        assert_eq!(
            fits.strata(),
            3,
            "every one of the three carries slippage, which is the count `strata` is about — it \
             is the larger of the two and they are not interchangeable"
        );
    }

    /// **A pool filed under a period that is not its own is refused.** The period is carried
    /// twice — as the map key and on the pool — and the lookup reads the key, so a disagreement
    /// would seed every tract of one period from another period's spread.
    #[test]
    #[should_panic(expected = "is filed under period 3")]
    fn a_pool_filed_under_another_period_is_refused() {
        let _ = StratumFits::over(&[], BTreeMap::new())
            .with_period_length_spectra(BTreeMap::from([(3, pool(2, vec![0.2, 0.5, 0.3], 9.0))]));
    }

    /// **A shape that misses being a distribution by a little is refused too**, which is what
    /// the tolerance is for — the only other fixture is off by a factor of ten, so a tolerance
    /// loosened from `1e-9` to `1e-2` would have changed nothing.
    ///
    /// A spectrum summing to 1.005 scales the tract prior by half a percent at every tract of
    /// that stratum, silently.
    #[test]
    #[should_panic(expected = "length shares totalling")]
    fn a_length_spectrum_off_by_half_a_percent_is_refused() {
        let _ = StratumFits::over(
            &[fitted(stratum(2, 10), 0.05, vec![0.605, 0.3, 0.1], 4.0)],
            one_group(),
        );
    }

    /// …and one inside the tolerance is accepted, so the check is a tolerance rather than an
    /// equality that only a renormalised vector could pass.
    #[test]
    fn a_length_spectrum_within_the_tolerance_is_accepted() {
        let fits = StratumFits::over(
            &[fitted(
                stratum(2, 10),
                0.05,
                vec![0.6 + 5e-10, 0.3, 0.1],
                4.0,
            )],
            one_group(),
        );
        assert_eq!(fits.strata_with_a_length_spectrum(), 1);
    }

    /// **The lookup's own type refuses a shape it cannot index**, and it is a second door: the
    /// gather's checks run once per run, this one runs at every lookup, and only this one stands
    /// between a hand-built `LengthSpectrum` and the seed builder.
    ///
    /// An empty `weights` is the case that was silent before the constructor existed: the span
    /// is `(0 − 1) / 2`, which overflows in debug and is `-1` in release, and against `-1` no
    /// candidate is ever in reach.
    #[test]
    #[should_panic(expected = "class count is odd and at least three")]
    fn a_length_spectrum_with_no_classes_cannot_be_built() {
        let _ = LengthSpectrum::fitted(&[], 4.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "class count is odd and at least three")]
    fn a_length_spectrum_with_an_even_class_count_cannot_be_built() {
        let _ = LengthSpectrum::fitted(&[0.25; 4], 4.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "chromosomes' worth of belief; got 0")]
    fn a_length_spectrum_held_with_no_conviction_cannot_be_built() {
        let _ = LengthSpectrum::fitted(&[0.25, 0.5, 0.25], 0.0, FittedFrom::ThisStratum);
    }

    #[test]
    #[should_panic(expected = "the stated-flat rung states a finite")]
    fn a_stated_flat_rung_with_no_conviction_cannot_be_built() {
        let _ = LengthSpectrum::stated_flat(f64::NAN);
    }

    /// A stratum's spectrum that does not sum to one is raw counts or a truncated copy, and
    /// seeding from it scales the prior by the total at every tract of that stratum.
    #[test]
    #[should_panic(expected = "the fit of period 3, 17 repeats has length shares totalling")]
    fn a_length_spectrum_that_is_not_a_distribution_is_refused() {
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
    #[should_panic(expected = "the pooled fit of motif period 2 has length shares totalling")]
    fn a_pooled_length_spectrum_that_is_not_a_distribution_is_refused_too() {
        let _ = StratumFits::over(&[], BTreeMap::new())
            .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![6.0, 3.0, 1.0], 4.0))]));
    }

    /// **Each length spectrum comes back naming the rung it came off**, and the two rungs are
    /// different claims: a stratum's own tracts against every tract of its motif period pooled.
    ///
    /// Nothing else in the crate reads which rung these two iterators stamp — the parameters
    /// file writes the spectrum's numbers and places it by which table it is in — so without
    /// this, swapping the two rungs changes no test.
    #[test]
    fn each_length_spectrum_names_the_rung_it_came_off() {
        let fits = StratumFits::over(
            &[
                fitted(stratum(2, 10), 0.05, leaning_short(), 4.0),
                derived(
                    stratum(2, 14),
                    vec![Some(slippage(0.09))],
                    vec![Some(from_the_cell(40.0))],
                ),
            ],
            one_group(),
        )
        .with_period_length_spectra(BTreeMap::from([(2, pool(2, vec![0.2, 0.5, 0.3], 2.5))]));

        assert_eq!(
            fits.fitted_length_spectrum_of_each_stratum()
                .map(|(at, shares, concentration)| (at, shares.to_vec(), concentration))
                .collect::<Vec<_>>(),
            vec![(stratum(2, 10), leaning_short(), 4.0)],
            "only the stratum fitted on its own tracts has a spectrum; the derived one carries \
             none by construction"
        );
        assert_eq!(
            fits.pooled_length_spectrum_of_each_period()
                .map(|(period, shares, concentration)| (period, shares.to_vec(), concentration))
                .collect::<Vec<_>>(),
            vec![(2, vec![0.2, 0.5, 0.3], 2.5)]
        );
        assert_eq!(
            fits.length_spectrum_at(2, 10).rung(),
            LengthSpectrumRung::StratumsOwnFit
        );
        assert_eq!(
            fits.length_spectrum_at(2, 14).rung(),
            LengthSpectrumRung::PeriodsPooledTracts,
            "the derived stratum falls to its period's pool, which is the middle rung"
        );
    }

    /// **A slippage group that put no read in a stratum is skipped rather than yielded**, which
    /// is the same claim [`StratumFits::at`] makes through `GroupPutNoReadHere` — and the reason
    /// the iterator's name says *with numbers*.
    #[test]
    fn each_stratum_and_group_with_numbers_skips_a_pair_with_none() {
        let two_groups = BTreeMap::from([(ReadGroupId(0), 0), (ReadGroupId(1), 1)]);
        let fits = StratumFits::over(
            &[
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
        let fits = StratumFits::over(
            &[StratumOutcome::Derived(Box::new(DerivedStratum {
                stratum: stratum(2, 10),
                slippage: vec![Some(slippage(0.05))],
                level_provenance: vec![None],
                shares_provenance: vec![None],
                tracts_of_its_own: 4,
                reads_crossing: 40,
            }))],
            one_group(),
        );
        let _ = fits.each_stratum_and_group_with_numbers().count();
    }
}
