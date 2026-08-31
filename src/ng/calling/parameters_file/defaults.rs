//! **What a run puts in each of its parameters when nothing measured one** — the compiled-in
//! defaults, in one place, together with the one parameter that has none and why it has none.
//!
//! **The design is `doc/devel/ng/spec/parameters_file.md` §8**, and every bare section number
//! below is that document's unless another one is named.
//!
//! **Why the defaults are in the binary rather than in a shipped file.** A user choosing to run
//! without a fit should not have to find a file on disk first, so *run with defaults* is a flag and
//! never a path (§8, owner's decision of 2026-08-28). What a shipped file would have bought —
//! defaults a person can see and edit — is bought instead by §7: **every run writes out the
//! parameters it used**, so a defaults run still produces the file, and each number below arrives
//! in it with `warrant = "defaulted"` beside it.
//!
//! # The whole list, and it is not one list
//!
//! Seven things a run needs and no fit gave it. §8 sorts the first four by what their default *is*;
//! the last three are the ones §8's three cases have no slot for, and they are why this table has
//! seven rows rather than four.
//!
//! | what a run needs | what it takes with no fit | what that is |
//! |---|---|---|
//! | the base-quality multiplier, per read group | [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER), one | the value at which the model does nothing |
//! | the repeat-tract outlier weight, one per run | [`DEFAULT_OUTLIER_WEIGHT`](crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT), 0.01 | another caller's number, never measured here |
//! | the tract ladder's fallback concentration | [`STATED_FLAT_CONCENTRATION`](crate::ng::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION), one | a stated uninformative prior |
//! | contamination, per read group | **absence** — no `[contamination]` section | a model state, not a guess |
//! | the repeat-tract substitution rate, per (read group × stratum) | [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE), 0.001 | a default taken at the tract, not written in the file |
//! | the slippage numbers, per (stratum × slippage group) | **no row**, and the tract falls back to another caller's shipped model | a default that is owed a measurement |
//! | the inbreeding coefficient, per sample | [`DEFAULT_INBREEDING_COEFFICIENT`], zero | Hardy–Weinberg, which the *fit* may not assume and the *run* may |
//!
//! The prior's seed is deliberately absent from the table: it has a fallback
//! ([`ExpectedHeterozygosity::SPECIES_FALLBACK`](crate::ng::types::ExpectedHeterozygosity::SPECIES_FALLBACK)),
//! but what records it is [`SeedRegime::FallbackDiversity`](crate::ng::calling::genotype_prior::SeedRegime),
//! a rung in the file's `[ordinary_site_prior]` section, and not a `warrant`. So it is marked, and
//! it is not marked the way the four above are — which is a difference §8 does not discuss.
//!
//! # What each row means for a run
//!
//! - **A multiplier of one declines to recalibrate; it does not abstain from a claim.** It leaves
//!   every read's error probability at what the instrument minted, which asserts the instrument was
//!   right — the assumption `read_likelihoods.md` §3.2 says the calibration exists to remove, and
//!   one that fitted multipliers refute routinely (a multiplier above one is common and says the
//!   instrument was optimistic). What is true is that it is the value at which the model does
//!   nothing, so it is the one default that cannot push a call in a particular direction. **It is
//!   also the only per-read-group default**, so a cohort can be half calibrated and half not.
//!
//!   **⚑ And it is the one row of this table whose `defaulted` warrant does not mean the run took
//!   the number in the middle column.** The other three warranted numbers are *taken* from a
//!   constant, so `defaulted` fixes the value; a multiplier is a fitted error **rate** divided by
//!   the geometric mean of that read group's minted error, and
//!   [`from_fitted_rate`](crate::ng::calling::likelihood::ReadGroupCalibration::from_fitted_rate)
//!   copies the **rate's** warrant onto the ratio. The pre-pass's error-rate ladder has a
//!   `Defaulted` bottom rung of its own —
//!   [`DEFAULT_ERROR_RATE`](crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE) at
//!   0.001, taken by a read group with too few sites to fit, no sibling to borrow from and nothing
//!   supplied — so a run can legitimately write a `defaulted` multiplier of `0.001 / that
//!   library's mean minted error`, which is one only by coincidence. **So the file's reader cannot
//!   hold this key to its constant the way it holds the other two**, and `validate` says at length
//!   why not.
//!
//!   **And that is deliberate — owner's ruling, 2026-08-31.** A library's real error rate is never
//!   its reported sequencing quality: the quality scores describe base calling, and the reads also
//!   carry mismapping, chimeras and damage. So a read group the fit could not measure is charged a
//!   stated rate rather than taken at its word, and on any real library that pushes the reads the
//!   *conservative* way — on HG002's mean minted error of 2.9055 × 10⁻⁴ the multiplier is 3.44,
//!   5.4 Phred less confident than the instrument claimed. Spec §5's third row says such a read
//!   group gets "scale 1.0" and is the sentence to correct; `DEFAULT_ERROR_RATE` is itself a
//!   placeholder until it is fitted from GIAB. Recorded in `PROJECT_STATUS.md`.
//! - **The outlier weight is another caller's number**, inherited from the existing caller at 0.01
//!   and never measured here (§3.8). It is the share of a repeat tract's reads the model expects to
//!   have come from somewhere it cannot explain — a chimera, a paralogous tract, a mismapped read.
//!   **Too low and a stray read has nowhere to go but into a genotype**, so a tract with one
//!   aberrant read is called over-confidently; too high and every repeat-tract genotype loses
//!   evidence to a term that explains nothing. Nothing here measures which way 0.01 errs; what a
//!   run can look at is how its repeat-tract calls move when the number is edited, which is why
//!   §3.8 puts it in the file: *marking a number soft is the point of writing it down.*
//! - **The fallback concentration is a stated uninformative prior, not an inherited one.** It is
//!   one chromosome's worth of belief spread flat over a tract's candidate lengths — the same
//!   quantity and the same reading `ALPHA_REF` carries on the ordinary-site path — and at one
//!   chromosome the reads move the prior from the first read onward. It is reached only where a run
//!   fitted no stratum at all; a run that fitted any takes its own median (§3.7).
//! - **Contamination's default is absence, and absence is a real answer.** A run told nothing about
//!   contamination is *scored as* uncontaminated: the read likelihood has no fraction to mix in, so
//!   it computes its plain two-term formula rather than the three-term mixture. That is a modelling
//!   default and not a finding about the samples — the file says the first ("nobody identified any
//!   contamination") and never the second. There is no constant here at all: the file writes no
//!   `[contamination]` section, which is §5's first row and the absence a reader most easily
//!   collapses into a table of zeros.
//! - **The substitution rate's default is taken at the tract and never written down.** The pre-pass
//!   emits the rate as `FittedHere` or not at all, so a `(read group × stratum)` it never
//!   accumulated has no row in the file; the cell that reaches it takes
//!   [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE)
//!   and `Provenance::Defaulted`, and
//!   [`TractScoringFits`](crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits)
//!   counts how many cells did. **So `defaulted` is not a warrant this file's substitution-rate
//!   rows can legitimately carry**, and nothing checks that they do not.
//!
//! # The one parameter with no default at all
//!
//! **The slippage numbers are owed a measurement, and until it exists they are the only row of
//! the table above whose middle column is not a number this caller holds.** The inbreeding
//! coefficient used to stand beside them and no longer does — the owner's ruling of 2026-08-31
//! gave it zero, and the bullet in the table says what taking that costs.
//!
//! - **Slippage, per (stratum × slippage group).** §8's third bullet decides these are to be fitted
//!   from the GIAB HG002 alignments and compiled in like the rest, and §12 question 1 records that
//!   the measurement does not exist. So a run with no slippage fit writes no slippage rows, and
//!   **the gap is filled one level down rather than left open**:
//!   [`repeat_tract_parameters`](crate::ng::calling::inference::repeat_tract_parameters) gives such
//!   a cell [`StutterModel::hipstr_shipped`](crate::ng::alignment::StutterModel::hipstr_shipped)
//!   with `Provenance::Defaulted` and counts it. **Those are HipSTR's shipped constants and not a
//!   fit**: one read in twenty comes back a whole repeat short and one in twenty a whole repeat
//!   long — symmetric, where `StutterModel::hipstr_shipped`'s own documentation records that
//!   HipSTR's *fitted* values are contraction-biased. So a run without slippage numbers scores its
//!   tracts under a symmetric guess taken from another caller, on no organism in particular, and
//!   [`cells_with_no_fitted_slippage`](crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits::cells_with_no_fitted_slippage)
//!   is what says how much of the run that was.
//! - **The inbreeding coefficient, per sample — and this row moved on 2026-08-31.** A run that
//!   fitted nothing has to put *some* coefficient into the genotype prior, and zero is the one that
//!   leaves Hardy–Weinberg in place. **That is a real assumption and not an abstention**, and it is
//!   wrong in a known direction on any selfing or structured cohort: the prior multiplies its
//!   heterozygote branch by `1 − F` (`calling_priors.md` §7), so a landrace whose true coefficient
//!   is 0.9 scored at zero has that branch **ten times** what it should be, and heterozygotes are
//!   over-called. A single sample always lands here, because there is no cohort to fit a
//!   coefficient from.
//!
//!   **The fit is still forbidden from taking it**, and the reason is not that zero is harmless.
//!   `parameter_estimation::generic::fallback`'s header states it: *"The inbreeding coefficient has
//!   one rung and it is not a default … it is the parameter that differs most between an outcrosser
//!   and a selfing landrace, and a cohort's diversity divides by `1 − F`, so a wrong constant would
//!   be amplified rather than absorbed."* **What differs is how far the error travels.** A fitted
//!   diversity divides by `1 − F` and carries the mistake into every downstream number the fit
//!   emits; a defaulted coefficient at calling time carries it into the calls and no further. So
//!   the run may take it and the fit may not — by the owner's ruling of 2026-08-31, which also
//!   gives the user one value for the whole run or a different value for any sample
//!   ([`DeclaredInbreeding`]). **The file marks every such coefficient `defaulted`**, which is how
//!   an operator who knows the crop sees that nobody has said. §3.5 requires at least one row and a
//!   defaults run now has one for every sample.
//!
//! # Where each constant lives, and why not here
//!
//! **Beside the code that reads it, and named there once.** A constant re-declared here would be a
//! second spelling of a number the caller already reads, and the two could then disagree — which is
//! the failure this file exists to make visible rather than one to introduce. So this module
//! documents the set and pins its behaviour; the numbers stay with their readers:
//! [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER)
//! beside [`ReadGroupCalibration`](crate::ng::calling::likelihood::ReadGroupCalibration), whose
//! `defaulted` constructor is the only thing that *reads* it — the projection in from a file builds
//! the struct literally, which is the path the reader's check below guards;
//! [`DEFAULT_OUTLIER_WEIGHT`](crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT) beside
//! [`RepeatTractOutlierWeight`](crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight);
//! [`STATED_FLAT_CONCENTRATION`](crate::ng::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION)
//! beside [`StratumFits`](crate::ng::parameter_estimation::joint::stratum_fits::StratumFits);
//! [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE)
//! beside the assembly that takes it.
//!
//! **The sentence a reader is shown beside each defaulted number is `to_toml`'s `origins`**, which
//! is the list this module reconciles against — one place, so a number's origin and the comment
//! above it in the file cannot drift apart.
//!
//! # Assembling a run that fitted nothing
//!
//! [`RunParameters::of_defaults`] is the door. **Not
//! [`assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)**, which takes the
//! *fit's* raw per-read-group maps and derives the run's read-group axis from them — a run with no
//! fit hands it two empty maps and it refuses a run with no read groups, which is the wrong
//! complaint about the right situation. **Not
//! [`of_gathered_values`](crate::ng::calling::run_parameters::RunParameters::of_gathered_values)
//! directly either**, though that is what this builds on: its nine arguments are the whole of a
//! run's parameters, and a caller assembling them by hand is a caller that can leave one of the
//! defaults out. So the third door takes what a run with no fit actually has — its read groups, its
//! ploidy, and what it was told about inbreeding — and fills in everything else from the list above.
//!
//! # What "the warrant says `defaulted`" is worth
//!
//! **A defaulted value and a fitted one are the same number.** A multiplier of one that nobody
//! fitted and a multiplier of one a fit arrived at multiply every read's error probability
//! identically, so nothing downstream of the arithmetic can tell them apart — the warrant beside
//! the value is the whole of the difference, and §5's third row is the requirement that it survive
//! a write and a read. The tests below hold each of the three constants to its warrant at the point
//! of use, and hold the file's reader to refusing a `defaulted` value that is not the number this
//! caller holds.

use std::collections::BTreeMap;

use crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight;
use crate::ng::calling::likelihood::{ContaminationView, ReadGroupCalibration};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::types::{InbreedingF, Ploidy, ReadGroupId};

/// **How inbred a sample is where nobody has said — zero, an outcrosser.**
///
/// **Zero is the value at which the genotype prior leaves the heterozygote branch alone.** The
/// prior multiplies that branch by `1 − F` (`doc/devel/ng/spec/calling_priors.md` §7), so at zero
/// the factor is exactly one and every genotype is scored under Hardy–Weinberg.
///
/// **That is a claim and not an abstention**, exactly as a base-quality multiplier of one is:
/// Hardy–Weinberg asserts random mating, no selfing and no structure, and on a selfing crop it is
/// large and wrong. What the two defaults share is that they are the value at which the model does
/// nothing, not that they assume nothing.
///
/// **⚑ The *fit* may not take this and the *run* may, and what separates them is how far the
/// error travels.** `parameter_estimation::generic::fallback` refuses to infer a coefficient from
/// data that cannot carry one, because a fitted diversity divides by `1 − F` and would carry the
/// mistake into every number the fit goes on to emit. A defaulted coefficient at calling time
/// carries it into the calls and stops there — which is smaller, not harmless. Owner's ruling,
/// 2026-08-31, which also gives the user both a run-wide value and a per-sample one — see
/// [`DeclaredInbreeding`].
///
/// **A run that knows better must say so**, and the file it writes is where a reader sees which
/// happened: `defaulted` at zero is *nobody said*, and `supplied` is *somebody did*. On a selfing
/// crop the difference is large — a landrace at `F = 0.9` scored at zero is told every homozygous
/// stretch of its genome is a surprise.
pub const DEFAULT_INBREEDING_COEFFICIENT: f64 = 0.0;

/// **What a run was told about how inbred its samples are** — nothing, one number for all of them,
/// or a number for some of them by name.
///
/// **Names, never positions.** A per-sample statement is joined to the run's sample order through
/// the run's own [`ReadGroups`], which is the rule spec §3.5 already settles for this quantity —
/// *"the file writes the name beside the value, because the order is the run's and a file that
/// carried only an order would be silently wrong against a re-ordered sample list"*. Joining on
/// row order makes two swapped rows a silent mis-attribution rather than an error.
///
/// **A sample nothing at all was said about takes [`DEFAULT_INBREEDING_COEFFICIENT`] and is marked
/// `Defaulted`.** *Nothing at all* is the operative phrase: a sample this does not name still takes
/// a run-wide value where one was given, marked `Supplied`, and falls to the default only where
/// there is no run-wide value either. So a cohort where the operator knew about three of its plants
/// and said nothing about the rest writes three `supplied` rows among `defaulted` ones; a cohort
/// where they also declared a run-wide coefficient writes `supplied` on every row.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DeclaredInbreeding {
    /// What every sample not named in `by_sample` takes. `None` is *nobody said anything at all*.
    everyone: Option<InbreedingF>,
    by_sample: BTreeMap<Box<str>, InbreedingF>,
}

impl DeclaredInbreeding {
    /// **Nobody said anything** — every sample takes [`DEFAULT_INBREEDING_COEFFICIENT`],
    /// `Defaulted`.
    #[must_use]
    pub fn nothing_said() -> Self {
        Self::default()
    }

    /// **One coefficient for every sample of the run**, `Supplied` — an operator who knows the
    /// cohort is a selfing crop and says so once.
    #[must_use]
    pub fn one_value_for_every_sample(coefficient: InbreedingF) -> Self {
        Self {
            everyone: Some(coefficient),
            by_sample: BTreeMap::new(),
        }
    }

    /// **This sample's own coefficient**, `Supplied`, overriding whatever the samples around it
    /// take. Naming one sample twice keeps the last statement.
    #[must_use]
    pub fn and_this_sample(mut self, sample: &str, coefficient: InbreedingF) -> Self {
        self.by_sample.insert(sample.into(), coefficient);
        self
    }

    /// One estimate per sample, **in the run's own sample order** — what
    /// [`RunParameters::of_defaults`] hands to the assembly, and what the parameters file writes as
    /// its `[inbreeding]` table.
    ///
    /// **`observations` is zero on every row, whichever warrant it carries**, because neither a
    /// stated constant nor a number an operator typed has any data behind it. The writer then
    /// leaves the key out of the file entirely — for a `defaulted` value because a stated constant
    /// has nothing behind it, and for a `supplied` one because a zero count would claim a
    /// measurement over no genome (`from_run_parameters`'s `warranted_value`, which is the only
    /// place either rule lives).
    ///
    /// **A name this was given that the run does not have is silently unused here**, deliberately:
    /// this type does not know the run's samples until it is handed them, and the caller that took
    /// the name from a command line is the one that can say *no sample of this run is called
    /// that*. [`Self::names_not_in`] is what that caller asks, and a caller that does not ask
    /// scores a plant under a coefficient meant for another.
    #[must_use]
    pub fn of_each_sample(&self, read_groups: &ReadGroups) -> Vec<Estimate<InbreedingF>> {
        read_groups
            .read_groups_per_sample()
            .iter()
            .map(|sample| {
                let stated = self
                    .by_sample
                    .get(sample.sample.as_ref())
                    .or(self.everyone.as_ref());
                match stated {
                    Some(coefficient) => Estimate {
                        value: *coefficient,
                        provenance: Provenance::Supplied,
                        observations: 0,
                    },
                    None => Estimate {
                        value: InbreedingF::try_new(DEFAULT_INBREEDING_COEFFICIENT)
                            .expect("zero is a coefficient in [0, 1)"),
                        provenance: Provenance::Defaulted,
                        observations: 0,
                    },
                }
            })
            .collect()
    }

    /// **Every sample this names that the run does not have, sorted by name**, for a caller that
    /// has to refuse a mistyped name rather than ignore it.
    ///
    /// **Sorted, not in the order they were given**, because the statements are held in a
    /// `BTreeMap` keyed by name and that order is gone by the time this is asked. It is the right
    /// order for a message a person reads, and it is the one thing about this list a caller can
    /// rely on.
    #[must_use]
    pub fn names_not_in(&self, read_groups: &ReadGroups) -> Vec<&str> {
        let of_the_run: std::collections::BTreeSet<&str> = read_groups
            .read_groups_per_sample()
            .iter()
            .map(|sample| sample.sample.as_ref())
            .collect();
        self.by_sample
            .keys()
            .map(std::convert::AsRef::as_ref)
            .filter(|name| !of_the_run.contains(name))
            .collect()
    }
}

impl RunParameters {
    /// **The parameters a run scores with when nothing was fitted and no file was supplied** —
    /// spec §8's defaults, assembled over this run's own read groups and samples.
    ///
    /// What the run brings is `read_groups` (from the alignment headers), `ploidy` (its own
    /// declaration, which §3.2 calls a property of the run rather than of the fit), and
    /// `inbreeding` (what it was told). **Everything else comes from the module list above**, and
    /// the point of the door is that a caller cannot leave one out.
    ///
    /// **The slippage group is declared and empty, which are two different things.** Every read
    /// group is declared into slippage group 0 — the run's own declaration, the same one
    /// `gather_strata` takes and the same default the joint walk uses — and no stratum has any
    /// numbers. So a repeat tract's lookup answers
    /// [`NoSlippage::NoSuchStratum`](crate::ng::parameter_estimation::joint::stratum_fits::NoSlippage::NoSuchStratum),
    /// *ordinary*, rather than
    /// [`UnknownReadGroup`](crate::ng::parameter_estimation::joint::stratum_fits::NoSlippage::UnknownReadGroup),
    /// which means *the run is not what it claims* and is counted apart for that reason. Declaring
    /// nothing would have made every cell of every tract report the second, and a defaults run is
    /// exactly what it claims.
    ///
    /// **The warrants on the inbreeding coefficients do not come back**, because this type does not
    /// hold them — it keeps what *calling* reads, which is the bare coefficient. A caller that also
    /// has to write the run's parameters file asks [`DeclaredInbreeding::of_each_sample`] for them,
    /// with the same two arguments; it is a pure function of those, so the two calls cannot
    /// disagree.
    ///
    /// # Panics
    ///
    /// On a run with no read groups or no samples. **The refusal is
    /// [`SequencingBatches::all_together`](crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches::all_together)'s**,
    /// which is evaluated before [`Self::of_gathered_values`] is entered and so is the message a
    /// caller sees; `of_gathered_values` refuses the same two shapes and never gets the chance.
    #[must_use]
    pub fn of_defaults(
        read_groups: &ReadGroups,
        ploidy: Ploidy,
        inbreeding: &DeclaredInbreeding,
    ) -> Self {
        // **A statement naming a plant this run does not have is a mistyped name, and it scores
        // that plant at the default in silence** — measured on a two-plant cohort where one name
        // was misspelled: both coefficients came back 0.0, and the operator's 0.42 reached
        // nothing. Held in debug only because the caller that took the name from a command line
        // is the one that can produce a message worth reading
        // ([`DeclaredInbreeding::names_not_in`]); this is what stops a test or an example
        // shipping the mistake unnoticed.
        debug_assert!(
            inbreeding.names_not_in(read_groups).is_empty(),
            "these samples were given inbreeding coefficients and are not in this run: {:?}. \
             A name that matches no sample is not applied to anything, and the plant it was \
             meant for is scored at the default instead",
            inbreeding.names_not_in(read_groups)
        );
        let slippage_group_of_each_read_group: BTreeMap<ReadGroupId, u32> = (0..read_groups.len())
            .map(|group| {
                (
                    ReadGroupId(u32::try_from(group).expect(
                        "a run has fewer read groups than a \
                                                             u32 can name",
                    )),
                    0,
                )
            })
            .collect();
        Self::of_gathered_values(
            // **⚑ Every read group takes a multiplier of one, which is the side of the module
            // header's ⚑ paragraph the owner ruled *against* for a read group the fit could not
            // measure** — there, a library nothing was fitted for is charged the pre-pass's
            // stated rate rather than taken at its reported quality. **A defaults run cannot be
            // charged that way**, and the reason is arithmetic rather than a second decision: the
            // multiplier is a rate divided by that library's own mean minted error, and the mean
            // comes from a `MintedReadErrors` accumulator only a pre-pass fills. A defaults run
            // has read no read, so it has no denominator. Spec §8's first bullet says 1.0, and
            // that is what this can do; a run that wanted the stated rate applied would need a
            // mean minted error from somewhere, which is a change to what a defaults run reads.
            vec![ReadGroupCalibration::defaulted(); read_groups.len()],
            // **Absence, not a row of zeros** (spec §5's first row): a run told nothing about
            // contamination is scored on the read likelihood's plain formula.
            Vec::<ContaminationView>::new(),
            SequencingBatches::all_together(read_groups),
            inbreeding
                .of_each_sample(read_groups)
                .into_iter()
                .map(|estimate| estimate.value)
                .collect(),
            // **The seed's own bottom rung**, which carries its origin in its regime rather than
            // in a warrant: no fitted frequency and no fitted heterozygosity gives
            // `SeedRegime::FallbackDiversity` at `ExpectedHeterozygosity::SPECIES_FALLBACK`.
            Self::seed_from_moments(None, None),
            StratumFits::over(&[], slippage_group_of_each_read_group),
            BTreeMap::new(),
            ploidy,
            RepeatTractOutlierWeight::defaulted(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{DEFAULT_INBREEDING_COEFFICIENT, DeclaredInbreeding};
    use crate::ng::calling::genotype_prior::SeedRegime;
    use crate::ng::calling::likelihood::ssr::{DEFAULT_OUTLIER_WEIGHT, RepeatTractOutlierWeight};
    use crate::ng::calling::likelihood::{
        DEFAULT_ERROR_PROBABILITY_MULTIPLIER, MIN_BASE_ERROR, ReadGroupCalibration,
    };
    use crate::ng::calling::parameters_file::ParametersFile;
    use crate::ng::calling::parameters_file::Warrant;
    use crate::ng::calling::parameters_file::tests::{
        THE_REFERENCE_A_RUN_FITTED_AGAINST, a_census_a_run_could_have_fitted_under,
        a_file_using_every_shape,
    };
    use crate::ng::calling::run_parameters::RunParameters;
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::stratum_fits::{
        LengthSpectrumRung, NoSlippage, STATED_FLAT_CONCENTRATION, StratumFits,
    };
    use crate::ng::read::input::read_groups::ReadGroups;
    use crate::ng::types::{ExpectedHeterozygosity, InbreedingF, Ploidy, ReadGroupId};
    use std::collections::BTreeMap;

    /// The cohort every `of_defaults` test below assembles over: **two lanes of one plant and one
    /// lane of another**, so the read-group axis (three) and the sample axis (two) have different
    /// lengths and a projection that worked only because they were equal would show.
    fn a_runs_read_groups() -> ReadGroups {
        ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", "Ailsa Craig", "lib5"),
        ])
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("two copies is a ploidy")
    }

    fn a_coefficient(value: f64) -> InbreedingF {
        InbreedingF::try_new(value).expect("a coefficient in [0, 1)")
    }

    /// **The multiplier a read group with no fitted rate is charged under leaves its reads
    /// exactly as the instrument minted them**, and says it was defaulted.
    ///
    /// The second half is what the assertion on the arithmetic is for: `charged_error` is
    /// `scale · exp(q_sum / n)`, so at a scale of one it is the geometric mean of the reads' own
    /// error probabilities and nothing else. Compared on `to_bits` as the module's other float
    /// assertions are — it is the same verdict as `==` on an ordinary finite double, and the
    /// habit is what step C3 established so that a lost sign on a zero cannot pass.
    #[test]
    fn a_read_group_with_no_fitted_rate_is_charged_what_its_reads_were_minted_with() {
        let defaulted = ReadGroupCalibration::defaulted();
        assert_eq!(defaulted.scale, DEFAULT_ERROR_PROBABILITY_MULTIPLIER);
        assert_eq!(defaulted.provenance, Provenance::Defaulted);

        // Three reads spanning Phred 40 to Phred 13, as the calibration's own fixture does.
        let minted = [1e-4_f64, 10f64.powf(-2.0), 10f64.powf(-1.3)];
        let q_sum: f64 = minted.iter().map(|error| error.ln()).sum();
        let reads = u32::try_from(minted.len()).expect("three reads fit in a u32");
        let geometric_mean = (q_sum / f64::from(reads)).exp();
        assert!(
            geometric_mean > MIN_BASE_ERROR,
            "the fixture must sit above the floor, or this asserts the floor instead"
        );
        assert_eq!(
            defaulted.charged_error(q_sum, reads).to_bits(),
            geometric_mean.to_bits(),
            "a defaulted calibration changes no read's error probability"
        );
        // `ln 1 = 0` exactly, so the log-space form adds nothing either.
        assert_eq!(defaulted.log_scale(), 0.0);
    }

    /// **The outlier weight a run inherits is the existing caller's 0.01 and says it was
    /// defaulted** — the one number here that is a guess at a quantity nothing in this project has
    /// measured.
    #[test]
    fn the_outlier_weight_a_run_inherits_is_the_existing_callers_number() {
        let inherited = RepeatTractOutlierWeight::defaulted();
        assert_eq!(inherited.value(), DEFAULT_OUTLIER_WEIGHT);
        assert_eq!(inherited.provenance(), Provenance::Defaulted);
        assert_eq!(DEFAULT_OUTLIER_WEIGHT, 0.01);
    }

    /// **A run that fitted no stratum states the flat concentration, says it was defaulted, and
    /// seeds every tract from the ladder's bottom rung.**
    ///
    /// All three matter together: the number alone cannot say where it came from — a run whose
    /// own strata's median happened to be 1.0 would carry the same number — and the rung is what
    /// a tract's own record shows.
    #[test]
    fn a_run_that_fitted_no_stratum_seeds_every_tract_from_the_flat_rung() {
        let nothing_fitted = StratumFits::over(&[], BTreeMap::new());
        assert_eq!(
            nothing_fitted.stated_concentration(),
            STATED_FLAT_CONCENTRATION
        );
        assert_eq!(
            nothing_fitted.stated_concentration_warrant(),
            Provenance::Defaulted
        );
        assert_eq!(STATED_FLAT_CONCENTRATION, 1.0);

        let spectrum = nothing_fitted.length_spectrum_at(2, 11);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::StatedFlat);
        assert_eq!(spectrum.concentration(), STATED_FLAT_CONCENTRATION);
        assert!(
            spectrum.fitted_weights().is_none(),
            "the bottom rung has no shape of its own to hand out"
        );
    }

    /// **The two parameters with no default are absent from the file a run with nothing fitted
    /// writes, and the file's own reader accepts that.**
    ///
    /// **Contamination's absence is already held elsewhere** — `to_run_parameters`'s
    /// `an_absent_contamination_table_is_not_a_table_of_zeros` is spec §5's first-row fixture and
    /// covers what the run then scores and reports, and `to_toml` pins that no section is
    /// written. What no test held is the pair above: that a file naming **no slippage row and no
    /// substitution-rate row at all** is a legal file rather than a broken one, which is what a
    /// defaults run will produce and what step E2 depends on.
    #[test]
    fn a_file_with_nothing_fitted_for_repeat_tracts_is_still_a_legal_file() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;
        file.repeat_tracts.slippage_by_stratum_and_group.clear();
        file.repeat_tracts.slippage_group_by_read_group.clear();
        file.repeat_tracts.substitution_rate_by_stratum.clear();
        file.repeat_tracts.length_spectrum_by_stratum.clear();
        file.repeat_tracts.length_spectrum_by_period.clear();
        // The bottom rung is then the only thing left saying anything about a tract, and its
        // warrant has to be the one a run that fitted nothing carries.
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .warrant = Warrant::Defaulted;
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .value = STATED_FLAT_CONCENTRATION;
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .observations = None;

        let projected = file
            .to_run_parameters()
            .expect("a run that fitted no repeat tracts writes a file its own reader accepts");
        let fits = projected.parameters.ssr_slippage_fits();
        assert_eq!(fits.strata(), 0);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::Defaulted);
        assert_eq!(projected.parameters.ssr_substitution_rate().count(), 0);
        assert!(projected.parameters.view().contamination_is_absent());
    }

    /// **The file's reader holds the two keys it can to their own constant**: a value marked
    /// `defaulted` that is not the number this caller holds is refused, naming the key, quoting
    /// both numbers as the file spells them, and saying what to type instead.
    ///
    /// This is the edit spec §7's third bullet invites — copy the file your run wrote and change
    /// one line — made to a number whose warrant the reader forgot to move. **Both refusals share
    /// one closing clause**, which is what the last assertion holds: a reader who has met one of
    /// them can act on the other without re-reading.
    ///
    /// **The base-quality multiplier is deliberately not here**, and the module header says why:
    /// its `defaulted` warrant is copied from the error rate the multiplier was built from, so a
    /// legitimate run writes one at a value that is not
    /// [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER).
    /// The test below is the one that pins it accepted.
    ///
    /// **The concentration's edited value is a whole number on purpose.** `Display` and `Debug`
    /// agree on 3.5 and differ on 3.0, so a fixture at 3.5 leaves the `{:?}` this message is
    /// written with unpinned — E1's review measured the revert to `{}` surviving its whole
    /// mutation suite. The outlier weight has no whole number to use, since a legal weight is
    /// strictly inside zero and one, so its `{:?}` is convention rather than something a test
    /// can hold.
    #[test]
    fn a_defaulted_value_that_is_not_the_binarys_own_number_is_refused() {
        let weight_edited = {
            let mut file = a_file_using_every_shape();
            file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Defaulted;
            file.stated_constants.repeat_tract_outlier_weight.value = DEFAULT_OUTLIER_WEIGHT * 2.0;
            file.stated_constants
                .repeat_tract_outlier_weight
                .observations = None;
            file
        };
        let concentration_edited = {
            let mut file = a_file_using_every_shape();
            let rung = &mut file.repeat_tracts.fallback_length_spectrum_concentration;
            rung.warrant = Warrant::Defaulted;
            rung.value = STATED_FLAT_CONCENTRATION + 2.0;
            rung.observations = None;
            file
        };

        for (edited, key, constant, edited_to) in [
            (
                weight_edited,
                "stated_constants.repeat_tract_outlier_weight",
                DEFAULT_OUTLIER_WEIGHT,
                DEFAULT_OUTLIER_WEIGHT * 2.0,
            ),
            (
                concentration_edited,
                "repeat_tracts.fallback_length_spectrum_concentration",
                STATED_FLAT_CONCENTRATION,
                STATED_FLAT_CONCENTRATION + 2.0,
            ),
        ] {
            let refusal = edited
                .validate()
                .expect_err("a `defaulted` value that is not the constant is refused")
                .to_string();
            assert!(
                refusal.contains(key),
                "the refusal must name the key to edit; got {refusal}"
            );
            // **Both numbers spelled as the file spells them**, which is `Debug` for a float —
            // the writer formats every value with it, so `3.0` is the string a reader can search
            // their own file for and `3` is not.
            assert!(
                refusal.contains(&format!("{constant:?}")),
                "the refusal must quote the number this caller holds, as the file spells it; \
                 got {refusal}"
            );
            assert!(
                refusal.contains(&format!("{edited_to:?}")),
                "the refusal must quote the number in the file, as the file spells it; \
                 got {refusal}"
            );
            // **One closing clause, word for word, on both.** Two sentences that each *mention*
            // `supplied` are two messages a reader has to parse separately; a reader who has
            // acted on one of these has acted on the other. The clause is the outlier weight's,
            // which was the first written and is the one that names the claim and the fix.
            assert!(
                refusal.ends_with(
                    "a number you changed is one the run was handed, so change the warrant \
                     beside it to `supplied`"
                ),
                "the two refusals must close with one clause, so a reader meets one shape of \
                 message rather than two; got {refusal}"
            );
        }
    }

    /// **Each of the three constants is accepted beside a `defaulted` warrant**, so the test
    /// above refuses the edit rather than the shape it was made in — and the multiplier, which
    /// nothing checks, is accepted at its constant too.
    #[test]
    fn each_of_the_three_constants_is_accepted_beside_a_defaulted_warrant() {
        let mut file = a_file_using_every_shape();

        let row = &mut file.base_quality_calibration.by_read_group[0];
        row.error_probability_multiplier.warrant = Warrant::Defaulted;
        row.error_probability_multiplier.value = DEFAULT_ERROR_PROBABILITY_MULTIPLIER;
        row.error_probability_multiplier.observations = None;

        file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Defaulted;
        file.stated_constants.repeat_tract_outlier_weight.value = DEFAULT_OUTLIER_WEIGHT;
        file.stated_constants
            .repeat_tract_outlier_weight
            .observations = None;

        let rung = &mut file.repeat_tracts.fallback_length_spectrum_concentration;
        rung.warrant = Warrant::Defaulted;
        rung.value = STATED_FLAT_CONCENTRATION;
        rung.observations = None;

        file.validate()
            .expect("the three constants are what a `defaulted` warrant claims");
    }

    // -----------------------------------------------------------------
    // E2 — a run with no fit and no supplied file
    // -----------------------------------------------------------------

    /// **A run that fitted nothing assembles, and every number it holds is the defaulted one.**
    ///
    /// This is step E2's whole claim, asserted field by field over the nine `RunParameters`
    /// holds — a run that got one of them from somewhere else would still assemble, and the
    /// genotypes would be wrong in a way nothing says.
    #[test]
    fn a_run_with_no_fit_takes_every_default() {
        let read_groups = a_runs_read_groups();
        let run = RunParameters::of_defaults(
            &read_groups,
            diploid(),
            &DeclaredInbreeding::nothing_said(),
        );

        assert_eq!(run.read_group_count(), 3);
        assert_eq!(run.ploidy(), diploid());

        for calibration in run.calibration_by_read_group() {
            assert_eq!(calibration.scale, DEFAULT_ERROR_PROBABILITY_MULTIPLIER);
            assert_eq!(calibration.provenance, Provenance::Defaulted);
        }

        // **Absence, not a row of zeros.** The empty axis is what `view` reads as the
        // uncontaminated run, which computes the read likelihood's plain formula.
        assert!(run.contamination_by_read_group().is_empty());
        assert!(run.view().contamination_is_absent());

        assert_eq!(run.inbreeding_coefficient_by_sample().len(), 2);
        // **Against the literal zero and not against the constant**, which is what makes this an
        // assertion rather than a restatement: comparing a run's coefficient to
        // `DEFAULT_INBREEDING_COEFFICIENT` moves both sides together. **Measured**: the first
        // draft of this module compared against the constant on every one of these lines, and a
        // build shipping it at 0.5 passed all 184 tests. Three assertions now compare against a
        // literal — this one and the two in
        // `a_stated_coefficient_is_supplied_and_an_unstated_one_is_defaulted` — and that build
        // fails all three.
        assert_eq!(DEFAULT_INBREEDING_COEFFICIENT, 0.0);
        for coefficient in run.inbreeding_coefficient_by_sample() {
            assert_eq!(coefficient.get(), 0.0);
            // What that zero buys: the prior weights its heterozygote branch by `1 − F`
            // (`calling_priors.md` §7, in the default path and in the `hardy_weinberg`
            // comparator alike), so at zero the branch is left alone and every genotype is
            // scored under Hardy–Weinberg. Asserting `1.0 - 0.0 == 1.0` here would restate the
            // line above rather than check anything.
        }

        // The prior's own bottom rung, which records itself in the regime rather than a warrant.
        let seed = run.prior_seed();
        assert_eq!(seed.regime(), SeedRegime::FallbackDiversity);
        assert_eq!(
            seed.alpha_alt_total(),
            ExpectedHeterozygosity::SPECIES_FALLBACK.get()
        );

        let fits = run.ssr_slippage_fits();
        assert_eq!(fits.strata(), 0);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::Defaulted);
        assert_eq!(run.ssr_substitution_rate().count(), 0);

        let weight = run.repeat_tract_outlier_weight();
        assert_eq!(weight.value(), DEFAULT_OUTLIER_WEIGHT);
        assert_eq!(weight.provenance(), Provenance::Defaulted);

        // The batching nobody declared: one batch holding the run.
        assert!(run.sequencing_batches().is_default());
        assert_eq!(run.sequencing_batches().batch_count(), 1);
    }

    /// **⚑ A defaults run's repeat tracts land on the *ordinary* absence, not on the one that
    /// means the parameters and the reads came from different runs.**
    ///
    /// `NoSlippage` has four arms and two of them
    /// [say the run is not what it claims](crate::ng::parameter_estimation::joint::stratum_fits::NoSlippage::UnknownReadGroup);
    /// the calling side counts those apart for exactly that reason
    /// (`TractScoringFits::cells_whose_read_group_the_fit_does_not_describe`). **A defaults run is
    /// precisely what it claims**, so it must not land there — and it would, on every cell of
    /// every tract, if `of_defaults` declared no slippage group. Declaring every read group into
    /// group 0 is what makes the answer `NoSuchStratum` instead, which is *ordinary*.
    #[test]
    fn a_defaults_runs_tracts_find_no_stratum_rather_than_an_unknown_read_group() {
        let read_groups = a_runs_read_groups();
        let run = RunParameters::of_defaults(
            &read_groups,
            diploid(),
            &DeclaredInbreeding::nothing_said(),
        );
        let fits = run.ssr_slippage_fits();

        for group in 0..3 {
            let read_group = ReadGroupId(group);
            assert_eq!(
                fits.slippage_group_of(read_group),
                Some(0),
                "every read group of the run is declared, and into one group"
            );
            assert_eq!(
                fits.at(read_group, 2, 11),
                Err(NoSlippage::NoSuchStratum),
                "the fit holds no stratum, which is ordinary; an undeclared read group would \
                 answer `UnknownReadGroup`, which says the parameters and the reads came from \
                 different runs"
            );
        }
    }

    /// **Nothing said, one value for the run, and a value for one plant — three states the file
    /// keeps apart.**
    ///
    /// The warrant is the whole of the difference and it is per row: a cohort where the operator
    /// knew about one of its plants writes one `supplied` row among `defaulted` ones.
    #[test]
    fn a_stated_coefficient_is_supplied_and_an_unstated_one_is_defaulted() {
        let read_groups = a_runs_read_groups();

        let nobody_said = DeclaredInbreeding::nothing_said().of_each_sample(&read_groups);
        assert_eq!(nobody_said.len(), 2);
        for estimate in &nobody_said {
            assert_eq!(estimate.value.get(), 0.0);
            assert_eq!(estimate.provenance, Provenance::Defaulted);
            assert_eq!(estimate.observations, 0);
        }

        let a_selfing_cohort = DeclaredInbreeding::one_value_for_every_sample(a_coefficient(0.9))
            .of_each_sample(&read_groups);
        for estimate in &a_selfing_cohort {
            assert_eq!(estimate.value.get(), 0.9);
            assert_eq!(estimate.provenance, Provenance::Supplied);
            assert_eq!(
                estimate.observations, 0,
                "a number an operator typed has no data behind it, whatever its warrant"
            );
        }

        // One plant named, the other left alone. The run's sample order is first-seen, so `TS-1`
        // is sample 0 and `Ailsa Craig` sample 1.
        let one_plant = DeclaredInbreeding::nothing_said()
            .and_this_sample("Ailsa Craig", a_coefficient(0.42))
            .of_each_sample(&read_groups);
        assert_eq!(one_plant[0].value.get(), 0.0);
        assert_eq!(one_plant[0].provenance, Provenance::Defaulted);
        assert_eq!(one_plant[1].value.get(), 0.42);
        assert_eq!(one_plant[1].provenance, Provenance::Supplied);
    }

    /// **A per-sample statement overrides the run-wide one, and is joined by name.**
    ///
    /// The join is what the assertion is really about. `TS-1` is the run's *first* sample and
    /// `Ailsa Craig` its second, so a builder that stored statements in the order they were made
    /// and zipped them against the run's samples would put the 0.42 on the wrong plant and pass
    /// every count-based check — which is D2's Blocker one level down.
    #[test]
    fn a_per_sample_statement_lands_on_that_sample_and_overrides_the_run_wide_one() {
        let read_groups = a_runs_read_groups();
        let mixed = DeclaredInbreeding::one_value_for_every_sample(a_coefficient(0.9))
            .and_this_sample("Ailsa Craig", a_coefficient(0.42))
            .of_each_sample(&read_groups);

        assert_eq!(mixed[0].value.get(), 0.9, "TS-1 keeps the run-wide value");
        assert_eq!(mixed[1].value.get(), 0.42, "Ailsa Craig takes its own");
        assert!(
            mixed
                .iter()
                .all(|estimate| estimate.provenance == Provenance::Supplied),
            "both were stated, by different statements"
        );
    }

    /// **A name no sample of the run carries is reported rather than silently ignored**, so a
    /// caller that took it off a command line can refuse a typo instead of scoring a plant under
    /// a coefficient meant for another.
    #[test]
    fn a_statement_naming_a_plant_the_run_does_not_have_is_reported() {
        let read_groups = a_runs_read_groups();
        let with_a_typo = DeclaredInbreeding::nothing_said()
            .and_this_sample("Ailsa Craigg", a_coefficient(0.42))
            .and_this_sample("TS-1", a_coefficient(0.3));

        assert_eq!(with_a_typo.names_not_in(&read_groups), vec!["Ailsa Craigg"]);

        // **Two typos come back sorted by name and not in the order they were typed**, which is
        // what the method promises: the statements are held keyed by name, so the order they
        // arrived in is gone by the time this is asked. `Zed` was named first.
        let two_typos = DeclaredInbreeding::nothing_said()
            .and_this_sample("Zed", a_coefficient(0.1))
            .and_this_sample("Alpha", a_coefficient(0.2));
        assert_eq!(two_typos.names_not_in(&read_groups), vec!["Alpha", "Zed"]);
        // And the run still assembles, with the typo having done nothing.
        let coefficients = with_a_typo.of_each_sample(&read_groups);
        assert_eq!(coefficients[0].value.get(), 0.3);
        assert_eq!(coefficients[1].value.get(), 0.0);
        assert_eq!(coefficients[1].provenance, Provenance::Defaulted);
    }

    /// **The smallest run there is takes the defaults too** — one sample, one library, which
    /// `CLAUDE.md` makes a first-class case. The two axes are the same length here, as they are
    /// in any run of one library a sample, so it is the read-group count of one rather than the
    /// equal lengths that this shape adds.
    #[test]
    fn one_sample_one_library_assembles_from_defaults() {
        let read_groups = ReadGroups::of_lanes(&[("HWI.3", "TS-1", "lib3")]);
        let run = RunParameters::of_defaults(
            &read_groups,
            diploid(),
            &DeclaredInbreeding::nothing_said(),
        );

        assert_eq!(run.read_group_count(), 1);
        assert_eq!(run.inbreeding_coefficient_by_sample().len(), 1);
        assert!(
            run.view().contamination_is_absent(),
            "at one sample there is no panel to be surprised by, so contamination is not \
             estimable at all — and a defaults run has not measured it either"
        );
    }

    /// **A defaults run's parameters write a file, and the file reads back into the same run.**
    ///
    /// This is what step F1 will do unconditionally (§7: every run writes the parameters it used,
    /// whatever the numbers came from), and it is the test that says the defaults are a *state the
    /// artefact can hold* rather than an in-memory shape. It goes the whole way — parameters →
    /// file → TOML → file → parameters — so nothing is proved about a file nobody could write.
    ///
    /// **The writer is handed no rates at all**, which is what a run with no fit has. That is
    /// legal exactly because every calibration is `Defaulted` and a defaulted number writes no
    /// count; `from_run_parameters` refuses a missing rate under any other warrant.
    ///
    /// **⚑ The census is the module's fixture and a direct-mode run has none.** `of_run` takes a
    /// `CensusIdentity` because §3.1 binds a fitted file to the census it was fitted from, and a
    /// run that fitted nothing has no census to name. What a defaults run writes there is F1's
    /// question and the owner's — it is the first of the three Milestone D left open.
    #[test]
    fn a_defaults_run_writes_a_file_that_reads_back_as_the_same_run() {
        let read_groups = a_runs_read_groups();
        let declared =
            DeclaredInbreeding::nothing_said().and_this_sample("TS-1", a_coefficient(0.9));
        let run = RunParameters::of_defaults(&read_groups, diploid(), &declared);

        let file = ParametersFile::of_run(
            &run,
            &read_groups,
            // No rates: nothing was fitted for anybody.
            &BTreeMap::new(),
            &declared.of_each_sample(&read_groups),
            &THE_REFERENCE_A_RUN_FITTED_AGAINST,
            a_census_a_run_could_have_fitted_under(),
        );

        file.validate()
            .expect("a defaults run writes a file its own reader accepts");

        // Through the text, as a run reading one back does.
        let text = file.to_toml();
        let parsed = ParametersFile::from_toml(&text).expect("the file this writer produced");
        let read_back = parsed
            .to_run_parameters()
            .expect("a defaults file means something");

        assert_eq!(
            read_back.parameters.read_group_count(),
            run.read_group_count()
        );
        assert_eq!(read_back.parameters.ploidy(), run.ploidy());
        assert_eq!(
            read_back.parameters.calibration_by_read_group(),
            run.calibration_by_read_group()
        );
        assert!(
            read_back
                .parameters
                .contamination_by_read_group()
                .is_empty()
        );
        assert_eq!(
            read_back.parameters.inbreeding_coefficient_by_sample(),
            run.inbreeding_coefficient_by_sample()
        );
        assert_eq!(read_back.parameters.prior_seed(), run.prior_seed());
        assert_eq!(
            read_back.parameters.repeat_tract_outlier_weight(),
            run.repeat_tract_outlier_weight()
        );
        assert_eq!(
            read_back.parameters.ssr_slippage_fits(),
            run.ssr_slippage_fits()
        );
        assert_eq!(read_back.parameters.ssr_substitution_rate().count(), 0);
        // The ninth field, so that all nine are checked *after* the trip and not eight after and
        // one before: a batching nobody declared has to come back saying nobody declared it.
        assert!(read_back.parameters.sequencing_batches().is_default());
        assert_eq!(
            read_back.parameters.sequencing_batches().batch_count(),
            run.sequencing_batches().batch_count()
        );

        // **The warrants survive the trip, which is the whole point of writing them.** One plant
        // was stated and the other was not, and the file is where a reader sees which.
        let warrants: Vec<_> = read_back
            .inbreeding_by_sample
            .iter()
            .map(|estimate| estimate.provenance)
            .collect();
        assert_eq!(warrants, vec![Provenance::Supplied, Provenance::Defaulted]);
        assert!(
            read_back
                .reads_behind_each_calibration
                .iter()
                .all(Option::is_none),
            "a defaulted multiplier has no count, and the file must not invent one"
        );
    }
}
