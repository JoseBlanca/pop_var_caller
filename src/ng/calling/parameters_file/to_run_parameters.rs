//! **A parameters file, back into the shape a run scores with** — the reverse of
//! `from_run_parameters`, and the half of step C2 that makes the round trip a round trip
//! (`doc/devel/ng/spec/parameters_file.md` §1.2 goal 1).
//!
//! # What it is not
//!
//! **It does not re-run the fit's assembly rules.**
//! [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble) takes
//! the pre-pass's raw outputs and decides things: it refuses a zero error rate and substitutes a
//! defaulted scale of one, it densifies the contamination axis, it decides whether the run has a
//! mixture at all. **A parameters file is the output of those decisions, not their input** — it
//! carries the scale and the warrant the fit settled on. Re-deriving them here would let a file's
//! stated calibration be quietly replaced by one computed from something else, which is the
//! failure spec §5 exists to prevent, one level up.
//! [`RunParameters::of_gathered_values`](crate::ng::calling::run_parameters::RunParameters::of_gathered_values)
//! is the door this uses instead, and the two sibling constructors on `StratumFits` and
//! `SequencingBatches` are there for the same reason.
//!
//! **It does not check the file against the run.** The four bindings of spec §6 — the reference's
//! digest, the sample list, the read-group table, the census — are step D2's, and this module
//! reads none of them except as the file's own axes: the read-group table for its length and the
//! sample list for its order. It hands none of them back, because its caller is already holding
//! the [`ParametersFile`] they are written in. A file that projects here is a file that *means
//! something*; whether it means it about **this** run is the next question.
//!
//! # Two things the file does not carry, and one it carries that memory cannot
//!
//! - **The fitted error rate.** The file writes the *multiplier* a rate produced and the count
//!   behind it, never the rate itself — `RunParameters` does not hold the rate either, so nothing
//!   is lost, but it means a projection cannot hand back an `Estimate<ErrorRate>`. What comes
//!   back instead is [`RunParametersFromFile::reads_behind_each_calibration`].
//! - **Which reads a slippage number was fitted from.** The file records the origin — this
//!   stratum's own fit, its period's curve, or a blend — and the curve itself, which is what a
//!   consumer weighs an interpolation against. That is all the fit kept too.
//! - **Every sample's inbreeding warrant and count**, which the run's own parameters drop: the
//!   seam into calling takes the bare coefficients. They come back in
//!   [`RunParametersFromFile::inbreeding_by_sample`], so a file written from a projection is the
//!   file that was read.
//!
//! # Where a refusal comes from
//!
//! [`ParametersFile::validate`] runs first and covers what a value *means* — a share outside its
//! range, a spectrum that is not a distribution, a contamination table nobody was measured in.
//! What is left here is the places where a value passes every one of those and a newtype's own
//! constructor still turns it down. **Three of those are reachable**: a motif period past
//! `MAX_MOTIF_LEN`, a repeat count past what 32 bits hold, and a slippage curve's rise shape
//! outside the range [`RiseShape`] admits, where `validate` checks only that it is a number. Two
//! more are written and unreachable — an inbreeding coefficient and a substitution rate refuse
//! exactly what `validate` already refuses — and a ploidy refuses only zero, which `validate`
//! catches at the top level and not on a substitution-rate row. They refuse in `validate`'s own
//! words and with the same key path, through the same [`refuse`] helper, so a reader meets one
//! shape of message rather than two.

use std::collections::BTreeMap;

use super::validate::refuse;
use super::{
    ContaminationFittedFrom, ContaminationMeasurement, CurveReach, EvidenceCount, LevelSmoothing,
    ParametersFile, ParametersFileError, SeedRung, ShareCurve, ShareCurveRung, ShareShape,
    ShareSmoothing, SlippageCurve, Warrant,
};
use crate::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
use crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight;
use crate::ng::calling::likelihood::{ContaminationView, ReadGroupCalibration};
use crate::ng::calling::run_parameters::{RunParameters, UNMEASURED_READ_GROUP};
use crate::ng::parameter_estimation::joint::census::Stratum;
use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use crate::ng::parameter_estimation::joint::share_curve::{
    ShareCurve as FittedShareCurve, ShareCurveSource, ShareShape as FittedShape, ShareSource,
};
use crate::ng::parameter_estimation::joint::slippage_curve::{
    CurveReach as FittedCurveReach, LevelSource, RiseShape, SlippageCurve as FittedSlippageCurve,
};
use crate::ng::parameter_estimation::joint::ssr_fit::{
    LevelProvenance, ShareProvenance, SharesProvenance, Slippage,
};
use crate::ng::parameter_estimation::joint::stratum_fits::{
    FittedLengthSpectrum, FittedSlippage, StratumFits,
};
use crate::ng::parameter_estimation::ssr::{RepeatCount, Stratum as SsrStratum, StratumKey};
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::types::{BatchId, ErrorRate, InbreedingF, Ploidy, ReadGroupId, SsrPeriod};

/// **What a parameters file gives a run back** — the parameters themselves, and the two things
/// about them the run's own shape does not keep.
///
/// **Three fields rather than one, because `RunParameters` is lossy on purpose.** It holds what
/// *calling* reads, and calling reads neither the count behind a multiplier nor the warrant on an
/// inbreeding coefficient. The file holds both, and a run that read a file and then wrote one
/// (spec §7 — every run writes the file it used) has to write back what it read.
#[derive(Debug)]
pub struct RunParametersFromFile {
    /// What calling scores with.
    pub parameters: RunParameters,
    /// How much data stood behind each read group's base-quality multiplier, in the run's dense
    /// read-group order — **`None` where the file states none**, which is what a defaulted scale
    /// carries.
    ///
    /// **The count and not the rate.** The fitted error rate is nowhere in the file: what was
    /// written is the multiplier it produced, which is what `RunParameters` keeps as well.
    pub reads_behind_each_calibration: Vec<Option<EvidenceCount>>,
    /// Each sample's inbreeding coefficient with its warrant and its count, in the run's sample
    /// order — the file's `[inbreeding]` table, whole, where `parameters` keeps only the values.
    pub inbreeding_by_sample: Vec<Estimate<InbreedingF>>,
}

impl ParametersFile {
    /// **Read this file back into the parameters a run scores with.**
    ///
    /// [`Self::validate`] runs first, so a file that parses and means nothing is refused here
    /// rather than becoming a plausible run. This is the caller spec §9's refusals exist for and
    /// the reason `validate` is not folded into [`Self::from_toml`]: that entry point answers a
    /// question about the *text*, and a caller reading a file to inspect it should not have to
    /// satisfy the run's constraints.
    ///
    /// # Errors
    ///
    /// [`ParametersFileError::Meaningless`], naming the key's full path in the file's own
    /// spelling — either from `validate` or from the few newtype constructors whose rules no
    /// shape can state (this module's header says which).
    pub fn to_run_parameters(&self) -> Result<RunParametersFromFile, ParametersFileError> {
        self.validate()?;

        let read_group_count = self.fitted_from.read_groups.len();
        let sample_count = self.fitted_from.samples.len();

        // **Every table is read by keyed lookup into the dense axis, never positionally.** The
        // rows may be written in any order — a hand-edited file often is — and `validate` has
        // already refused a table that does not cover `0..n` once each, so every lookup below
        // finds its row.
        let mut calibration_by_read_group = vec![None; read_group_count];
        let mut reads_behind_each_calibration = vec![None; read_group_count];
        for row in &self.base_quality_calibration.by_read_group {
            let at = row.read_group as usize;
            calibration_by_read_group[at] = Some(ReadGroupCalibration {
                scale: row.error_probability_multiplier.value,
                provenance: row.error_probability_multiplier.warrant.into(),
            });
            reads_behind_each_calibration[at] = row.error_probability_multiplier.observations;
        }
        let calibration_by_read_group: Vec<ReadGroupCalibration> = calibration_by_read_group
            .into_iter()
            .map(|found| found.expect("validate refuses a calibration table with a gap in it"))
            .collect();

        // **One index over the sample names, built once and shared.** Both per-sample tables are
        // read by name into the run's sample order, and resolving each row with a linear scan
        // would be a quadratic walk over the file's one cohort-sized axis — 9 million string
        // comparisons at the 3,000 samples `CLAUDE.md` commits to, twice over.
        let sample_at: BTreeMap<&str, usize> = self
            .fitted_from
            .samples
            .iter()
            .enumerate()
            .map(|(at, sample)| (sample.as_str(), at))
            .collect();

        let contamination_by_read_group = self.contamination_views(read_group_count);
        let sequencing_batches = self.batching(read_group_count, sample_count, &sample_at);
        let inbreeding_by_sample = self.inbreeding_estimates(&sample_at)?;
        let prior_seed = SpectrumSeed::new(
            self.ordinary_site_prior.reference_concentration,
            self.ordinary_site_prior.alternative_concentration_total,
            self.ordinary_site_prior.rung.into(),
        );

        Ok(RunParametersFromFile {
            parameters: RunParameters::of_gathered_values(
                calibration_by_read_group,
                contamination_by_read_group,
                sequencing_batches,
                inbreeding_by_sample
                    .iter()
                    .map(|estimate| estimate.value)
                    .collect(),
                prior_seed,
                self.slippage_fits()?,
                self.substitution_rates()?,
                a_ploidy("ploidy", self.ploidy)?,
                self.outlier_weight(),
            ),
            reads_behind_each_calibration,
            inbreeding_by_sample,
        })
    }

    /// **The contamination axis, dense — or empty, which is the uncontaminated run.**
    ///
    /// Spec §5's first row, and the one absence a reader most easily collapses: no
    /// `[contamination]` section at all means nobody identified a fraction, and the read
    /// likelihood then computes its plain formula rather than a mixture with every share zero.
    /// A row with no `measurement` is the *other* absence — this library was not measured inside
    /// a run where others were — and it becomes a view of zero fraction and zero evidence, which
    /// `ContaminationView::was_measured` still tells apart from a library measured and found
    /// clean.
    fn contamination_views(&self, read_group_count: usize) -> Vec<ContaminationView> {
        let Some(table) = &self.contamination else {
            return Vec::new();
        };
        let mut views = vec![None; read_group_count];
        for row in &table.by_read_group {
            views[row.read_group as usize] = Some(match &row.measurement {
                Some(measurement) => a_measured_view(measurement),
                None => UNMEASURED_READ_GROUP,
            });
        }
        views
            .into_iter()
            .map(|found| found.expect("validate refuses a contamination table with a gap in it"))
            .collect()
    }

    /// The two batch columns and the flag that says whether anybody declared them.
    fn batching(
        &self,
        read_group_count: usize,
        sample_count: usize,
        sample_at: &BTreeMap<&str, usize>,
    ) -> SequencingBatches {
        let mut of_each_read_group = vec![BatchId::ALL_TOGETHER; read_group_count];
        for row in &self.sequencing_batches.by_read_group {
            of_each_read_group[row.read_group as usize] = BatchId(row.batch);
        }
        // **The sample column is read by the file's own sample order**, which is
        // `fitted_from.samples` — the order every per-sample axis of the run is indexed in, and
        // the order `validate` has already checked this table names exactly.
        let mut of_each_sample = vec![BatchId::ALL_TOGETHER; sample_count];
        for row in &self.sequencing_batches.by_sample {
            of_each_sample[sample_at[row.sample.as_str()]] = BatchId(row.batch);
        }
        let batch_count = of_each_read_group
            .iter()
            .chain(&of_each_sample)
            .map(|batch| batch.get() as usize + 1)
            .max()
            .unwrap_or(1);
        SequencingBatches::of_gathered_columns(
            of_each_read_group,
            of_each_sample,
            batch_count,
            !self.sequencing_batches.batching_was_declared,
        )
    }

    /// Each sample's coefficient, its warrant and its count, in the file's own sample order.
    fn inbreeding_estimates(
        &self,
        sample_at: &BTreeMap<&str, usize>,
    ) -> Result<Vec<Estimate<InbreedingF>>, ParametersFileError> {
        let mut by_sample = vec![None; self.fitted_from.samples.len()];
        for row in &self.inbreeding.by_sample {
            let at = sample_at[row.sample.as_str()];
            let field = format!(
                "inbreeding.by_sample[{:?}].inbreeding_coefficient",
                row.sample
            );
            let value = InbreedingF::try_new(row.inbreeding_coefficient.value).map_err(|_| {
                refuse(
                    &field,
                    format!(
                        "is {}, and an inbreeding coefficient is a fraction in [0, 1)",
                        row.inbreeding_coefficient.value
                    ),
                )
            })?;
            by_sample[at] = Some(Estimate {
                value,
                provenance: row.inbreeding_coefficient.warrant.into(),
                observations: an_evidence_count(row.inbreeding_coefficient.observations),
            });
        }
        Ok(by_sample
            .into_iter()
            .map(|found| found.expect("validate refuses a per-sample table with a sample missing"))
            .collect())
    }

    /// The whole `[repeat_tracts]` section except the substitution rates, which are keyed
    /// differently.
    fn slippage_fits(&self) -> Result<StratumFits, ParametersFileError> {
        let tracts = &self.repeat_tracts;
        let slippage_group_of: BTreeMap<ReadGroupId, u32> = tracts
            .slippage_group_by_read_group
            .iter()
            .map(|row| (ReadGroupId(row.read_group), row.slippage_group))
            .collect();

        // **One width for every stratum, taken over the whole section.** The fit builds the three
        // per-group vectors from one mask, so every stratum of a run has as many slots as the run
        // has slippage groups — and a lookup reads that width to tell *this group put no read
        // here* from *this group is not in the fit at all*. Taking the maximum over both the
        // declaration and the rows keeps those two answers the ones the fit would have given.
        let groups = tracts
            .slippage_group_by_read_group
            .iter()
            .map(|row| row.slippage_group)
            .chain(
                tracts
                    .slippage_by_stratum_and_group
                    .iter()
                    .map(|row| row.slippage_group),
            )
            .max()
            .map_or(0, |highest| highest as usize + 1);

        let mut slippage_by_stratum: BTreeMap<Stratum, Vec<Option<FittedSlippage>>> =
            BTreeMap::new();
        for row in &tracts.slippage_by_stratum_and_group {
            let stratum = Stratum {
                period: row.period,
                reference_repeats: row.reference_repeats,
            };
            let at = format!(
                "repeat_tracts.slippage_by_stratum_and_group[period = {}, reference_repeats = {}, slippage_group = {}]",
                row.period, row.reference_repeats, row.slippage_group
            );
            let level = a_level_provenance(&at, &row.share_of_reads_that_slip_origin)?;
            let shares = row
                .shorter_share_and_fall_off_origin
                .as_ref()
                .map(|origin| SharesProvenance {
                    slipped_reads: origin.expected_slipped_reads,
                    shorter_share: a_share_provenance(&origin.shorter_share_smoothing),
                    fall_off: a_share_provenance(&origin.fall_off_smoothing),
                });
            slippage_by_stratum
                .entry(stratum)
                .or_insert_with(|| vec![None; groups])[row.slippage_group as usize] =
                Some(FittedSlippage {
                    slippage: Slippage {
                        level: row.share_of_reads_that_slip,
                        shorter_share: row.shorter_share,
                        fall_off: row.fall_off,
                    },
                    level,
                    shares,
                });
        }

        let length_spectrum_by_stratum = tracts
            .length_spectrum_by_stratum
            .iter()
            .map(|row| {
                (
                    Stratum {
                        period: row.period,
                        reference_repeats: row.reference_repeats,
                    },
                    FittedLengthSpectrum {
                        weights: row.shares_by_repeat_offset.clone(),
                        concentration: row.concentration,
                    },
                )
            })
            .collect();
        let length_spectrum_by_period = tracts
            .length_spectrum_by_period
            .iter()
            .map(|row| {
                (
                    row.period,
                    FittedLengthSpectrum {
                        weights: row.shares_by_repeat_offset.clone(),
                        concentration: row.concentration,
                    },
                )
            })
            .collect();

        Ok(StratumFits::of_gathered_rows(
            slippage_group_of,
            slippage_by_stratum,
            length_spectrum_by_stratum,
            length_spectrum_by_period,
            tracts.fallback_length_spectrum_concentration.value,
        ))
    }

    /// The per-`(read group × stratum × ploidy)` substitution rates, keyed as the fit keys them.
    fn substitution_rates(
        &self,
    ) -> Result<BTreeMap<StratumKey, Estimate<ErrorRate>>, ParametersFileError> {
        let mut rates = BTreeMap::new();
        for row in &self.repeat_tracts.substitution_rate_by_stratum {
            let at = format!(
                "repeat_tracts.substitution_rate_by_stratum[read_group = {}, period = {}, reference_repeats = {}, ploidy = {}]",
                row.read_group, row.period, row.reference_repeats, row.ploidy
            );
            let key = StratumKey {
                read_group: ReadGroupId(row.read_group),
                stratum: SsrStratum::new(
                    a_period(&at, row.period)?,
                    a_repeat_count(&at, row.reference_repeats)?,
                ),
                ploidy: a_ploidy(&at, row.ploidy)?,
            };
            let value = ErrorRate::try_new(row.rate.value).map_err(|_| {
                refuse(
                    format!("{at}.rate"),
                    format!(
                        "is {}, and a substitution rate is a probability",
                        row.rate.value
                    ),
                )
            })?;
            rates.insert(
                key,
                Estimate {
                    value,
                    provenance: row.rate.warrant.into(),
                    observations: an_evidence_count(row.rate.observations),
                },
            );
        }
        Ok(rates)
    }

    /// The one stated constant, and which of its two warrants the file gave it.
    ///
    /// **`validate` has already held it to those two**, so this match has no third arm to
    /// invent: `fitted_here` and `borrowed` are claims about a number nothing fits, and a
    /// `defaulted` value that is not the compiled-in constant is a number somebody changed
    /// without saying so.
    fn outlier_weight(&self) -> RepeatTractOutlierWeight {
        match self.stated_constants.repeat_tract_outlier_weight.warrant {
            Warrant::Defaulted => RepeatTractOutlierWeight::defaulted(),
            Warrant::Supplied => RepeatTractOutlierWeight::supplied(
                self.stated_constants.repeat_tract_outlier_weight.value,
            ),
            // **Spelled out rather than caught by a wildcard**, so that a fifth warrant is a
            // compile error here as it is in every other conversion in this file, instead of
            // arriving silently as `supplied`.
            fitted @ (Warrant::FittedHere | Warrant::Borrowed) => unreachable!(
                "`validate` holds this key to `supplied` or `defaulted` — nothing fits it — and \
                 a warrant of {fitted:?} reached the projection, so the two have come apart"
            ),
        }
    }
}

fn a_measured_view(measurement: &ContaminationMeasurement) -> ContaminationView {
    ContaminationView {
        fraction: measurement.share_of_reads_from_another_sample,
        markers_with_reads: measurement.markers_with_reads,
        reads_on_markers: measurement.reads_on_markers,
        source: measurement.fitted_from_reads_of.into(),
    }
}

/// **An absent evidence count is zero**, which is the shape `Estimate<T>` has: its
/// `observations` is a bare `u64` and the file's is an `Option`.
///
/// **Not a loss, because the two absences agree.** The projection out writes no count for a
/// `defaulted` value — a stated constant has nothing behind it — and a count of zero says the
/// same thing in the shape memory has. What the file adds is the ability to say it by leaving
/// the key out, which is what makes a hand-written file legible.
fn an_evidence_count(count: Option<EvidenceCount>) -> u64 {
    match count {
        // **The unit is dropped here and nowhere else.** `Estimate<T>`'s count is a bare `u64`
        // whose unit follows the quantity; the file names the unit because a reader cannot be
        // sent to the source to find it. The projection out puts the right one back, per
        // quantity, so the trip is lossless — the file's own
        // `every_evidence_count_names_the_unit_its_quantity_is_fitted_over` is what holds that.
        Some(
            EvidenceCount::Reads(count)
            | EvidenceCount::CoveredPositions(count)
            | EvidenceCount::BasesCompared(count),
        ) => count,
        None => 0,
    }
}

fn a_period(at: &str, period: u8) -> Result<SsrPeriod, ParametersFileError> {
    SsrPeriod::try_new(usize::from(period)).map_err(|_| {
        refuse(
            format!("{at}.period"),
            format!("is {period}, and this build indexes motif periods this caller can score"),
        )
    })
}

fn a_repeat_count(at: &str, repeats: u64) -> Result<RepeatCount, ParametersFileError> {
    u32::try_from(repeats).map(RepeatCount).map_err(|_| {
        refuse(
            format!("{at}.reference_repeats"),
            format!("is {repeats}, and a tract's repeat count is held in 32 bits"),
        )
    })
}

fn a_ploidy(at: &str, ploidy: u8) -> Result<Ploidy, ParametersFileError> {
    Ploidy::try_new(ploidy).map_err(|_| {
        refuse(
            if at == "ploidy" {
                "ploidy".to_owned()
            } else {
                format!("{at}.ploidy")
            },
            format!("is {ploidy}, and a run calls a number of copies this caller can tabulate"),
        )
    })
}

fn a_level_provenance(
    at: &str,
    origin: &super::LevelOrigin,
) -> Result<LevelProvenance, ParametersFileError> {
    let at = format!("{at}.share_of_reads_that_slip_origin.smoothing");
    let (source, curve, reach) = match &origin.smoothing {
        LevelSmoothing::ThisStratum => (LevelSource::Cell, None, None),
        LevelSmoothing::ThisPeriodsCurve { curve, reach } => (
            LevelSource::Curve,
            Some(a_slippage_curve(&at, curve)?),
            Some((*reach).into()),
        ),
        LevelSmoothing::Blend {
            curve_weight,
            curve,
            reach,
        } => (
            LevelSource::Blend {
                curve_weight: *curve_weight,
            },
            Some(a_slippage_curve(&at, curve)?),
            Some((*reach).into()),
        ),
    };
    Ok(LevelProvenance {
        source,
        curve,
        reach,
        slipped_reads: origin.expected_slipped_reads,
    })
}

/// **Infallible where its level counterpart is not**, because a share curve carries no newtype:
/// every one of its fields is a plain float or count, and `validate` has already refused a
/// weight outside `[0, 1]` and any field that is not a number.
fn a_share_provenance(smoothing: &ShareSmoothing) -> ShareProvenance {
    let (source, curve, reach) = match smoothing {
        ShareSmoothing::ThisStratum => (ShareSource::Stratum, None, None),
        ShareSmoothing::ThisPeriodsCurve { curve, reach } => (
            ShareSource::Curve,
            Some(a_share_curve(curve)),
            Some((*reach).into()),
        ),
        ShareSmoothing::Blend {
            curve_weight,
            curve,
            reach,
        } => (
            ShareSource::Blend {
                curve_weight: *curve_weight,
            },
            Some(a_share_curve(curve)),
            Some((*reach).into()),
        ),
    };
    ShareProvenance {
        source,
        curve,
        reach,
    }
}

/// **The only curve field with a rule no range check states**: the rise shape is a newtype whose
/// constructor refuses what it refuses, and the file writes it as a bare float.
fn a_slippage_curve(
    at: &str,
    curve: &SlippageCurve,
) -> Result<FittedSlippageCurve, ParametersFileError> {
    let rise_shape = RiseShape::new(curve.rise_shape).ok_or_else(|| {
        refuse(
            format!("{at}.curve.rise_shape"),
            format!(
                "is {}, and a slippage curve's rise shape is a number this build can hold",
                curve.rise_shape
            ),
        )
    })?;
    Ok(FittedSlippageCurve {
        rise_shape,
        intercept: curve.intercept,
        slope: curve.slope,
        fitted_from: curve.fitted_from_repeats,
        fitted_to: curve.fitted_to_repeats,
        held_out_error: curve.held_out_error,
        cells: curve.cells as usize,
    })
}

fn a_share_curve(curve: &ShareCurve) -> FittedShareCurve {
    FittedShareCurve {
        shape: curve.shape.into(),
        intercept: curve.intercept,
        slope: curve.slope,
        bend: curve.bend,
        centre: curve.centre_repeats,
        fitted_from: curve.fitted_from_repeats,
        fitted_to: curve.fitted_to_repeats,
        held_out_error: curve.held_out_error,
        strata: curve.strata as usize,
        source: curve.curve_fitted_on.into(),
    }
}

// ---------------------------------------------------------------------
// The six words the file spells, back into the pre-pass's own
// ---------------------------------------------------------------------
//
// **Six, where `from_run_parameters` has eight `From` impls**: two of those eight convert a whole
// curve rather than a word, and this file's counterparts are plain functions because one of them
// can refuse (`a_slippage_curve`). The six here are the vocabularies.
//
// **Every one is an exhaustive match**, as its counterpart going the other way is, so a variant
// added to either vocabulary stops this file compiling rather than being silently dropped.
// **Exhaustiveness cannot catch a variant crossed with another** — `Flat => Sloping` compiles and
// spells correctly — which is what `every_word_the_file_spells_reads_back_to_its_own` names pair
// by pair.

impl From<Warrant> for Provenance {
    fn from(warrant: Warrant) -> Self {
        match warrant {
            Warrant::FittedHere => Self::FittedHere,
            Warrant::Borrowed => Self::Borrowed,
            Warrant::Supplied => Self::Supplied,
            Warrant::Defaulted => Self::Defaulted,
        }
    }
}

impl From<ContaminationFittedFrom> for ContaminationSource {
    fn from(from: ContaminationFittedFrom) -> Self {
        match from {
            ContaminationFittedFrom::ThisReadGroupsOwnReads => Self::ThisReadGroupsReads,
            ContaminationFittedFrom::EveryReadOfThisSample => Self::TheWholeSamplesReads,
        }
    }
}

impl From<SeedRung> for SeedRegime {
    fn from(rung: SeedRung) -> Self {
        match rung {
            SeedRung::FittedCurve => Self::FittedCurve,
            SeedRung::NeutralShape => Self::NeutralShape,
            SeedRung::ZeroDiversity => Self::ZeroDiversity,
            SeedRung::StatedHeterozygosity => Self::FallbackDiversity,
        }
    }
}

impl From<CurveReach> for FittedCurveReach {
    fn from(reach: CurveReach) -> Self {
        match reach {
            CurveReach::InsideTheFittedRange => Self::Inside,
            CurveReach::BelowTheFittedRange => Self::BelowFitted,
            CurveReach::AboveTheFittedRange => Self::AboveFitted,
        }
    }
}

impl From<ShareCurveRung> for ShareCurveSource {
    fn from(rung: ShareCurveRung) -> Self {
        match rung {
            ShareCurveRung::ThisPeriod => Self::ThisPeriod,
            ShareCurveRung::ThisPeriodUnscored => Self::ThisPeriodUnscored,
            ShareCurveRung::OtherPeriods => Self::OtherPeriods,
            ShareCurveRung::BuiltInDefault => Self::BuiltInDefault,
        }
    }
}

impl From<ShareShape> for FittedShape {
    fn from(shape: ShareShape) -> Self {
        match shape {
            ShareShape::Flat => Self::Flat,
            ShareShape::Sloping => Self::Sloping,
            ShareShape::Turning => Self::Turning,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::a_file_using_every_shape;
    use super::super::{ParametersFile, ParametersFileError, Warrant};
    use crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT;
    use crate::ng::parameter_estimation::{Estimate, Provenance};
    use crate::ng::read::input::read_groups::ReadGroups;
    use crate::ng::types::{ErrorRate, ReadGroupId};
    use std::collections::BTreeMap;

    /// The run's read-group table, rebuilt from the identity block the file carries.
    ///
    /// **Step D2's job is to compare this against the run's own**, and this is not that: it is
    /// the table the file *says* the numbers were fitted from, which is the only one a round
    /// trip through the file can use.
    pub(super) fn the_files_read_groups(file: &ParametersFile) -> ReadGroups {
        let lanes: Vec<(&str, &str, &str)> = file
            .fitted_from
            .read_groups
            .iter()
            .map(|row| {
                (
                    row.declared_id.as_str(),
                    row.sample.as_str(),
                    row.library.as_str(),
                )
            })
            .collect();
        ReadGroups::of_lanes(&lanes)
    }

    /// **The rates `of_run` wants, rebuilt from what the file actually carries.**
    ///
    /// **The file does not carry the fitted error rate and never has** — it carries the
    /// *multiplier* that rate produced, which is what `RunParameters` keeps too. `of_run` takes
    /// the rates because it reads two things off them: the warrant, which it checks against the
    /// calibration's, and the observation count, which it writes into the file. Both come back
    /// from the projection; the rate's own value does not, and nothing reads it, so this fills a
    /// placeholder rather than pretending to recover one.
    pub(super) fn the_rates_the_projection_out_reads(
        projected: &super::RunParametersFromFile,
    ) -> BTreeMap<ReadGroupId, Estimate<ErrorRate>> {
        projected
            .parameters
            .calibration_by_read_group()
            .iter()
            .zip(&projected.reads_behind_each_calibration)
            .enumerate()
            .map(|(read_group, (calibration, count))| {
                (
                    ReadGroupId(read_group as u32),
                    Estimate {
                        value: ErrorRate::try_new(1e-3).expect("a placeholder rate"),
                        provenance: calibration.provenance,
                        observations: super::an_evidence_count(*count),
                    },
                )
            })
            .collect()
    }

    /// **A file, read into a run's parameters and written back out, is the file that was read.**
    ///
    /// This is the round trip the whole design rests on (spec §1.2 goal 1), on the fixture that
    /// uses every shape: three read groups over two samples, one library measured for
    /// contamination and one not, a defaulted calibration beside two that were not, four strata
    /// with three different slippage origins between them, a stratum's own length spectrum and a
    /// period's pooled one, and a substitution rate keyed by all three of its axes.
    ///
    /// **It compares the two `ParametersFile`s rather than the two `RunParameters`.** The file
    /// is the artefact; and `RunParameters` has no `PartialEq`, because it holds a borrowed-from
    /// bundle nothing else compares.
    #[test]
    fn a_file_read_into_a_run_and_written_back_is_the_file_that_was_read() {
        let file = a_file_using_every_shape();
        let read_groups = the_files_read_groups(&file);
        let projected = file
            .to_run_parameters()
            .expect("the fixture is a usable file");

        let written = ParametersFile::of_run(
            &projected.parameters,
            &read_groups,
            &the_rates_the_projection_out_reads(&projected),
            &projected.inbreeding_by_sample,
            &file.fitted_from.reference_digest,
            file.fitted_from.census.clone(),
        );
        assert_eq!(written, file);
    }

    /// **And the text survives the whole trip**, which is the loop a run actually makes: a file
    /// on disk, into memory, scored with, and written out beside the VCF (spec §7).
    #[test]
    fn the_text_survives_the_trip_through_memory() {
        let file = a_file_using_every_shape();
        let read_groups = the_files_read_groups(&file);
        let text = file.to_toml();

        let read = ParametersFile::from_toml(&text).expect("the writer's own text parses");
        let projected = read.to_run_parameters().expect("and means something");
        let written = ParametersFile::of_run(
            &projected.parameters,
            &read_groups,
            &the_rates_the_projection_out_reads(&projected),
            &projected.inbreeding_by_sample,
            &read.fitted_from.reference_digest,
            read.fitted_from.census.clone(),
        );
        assert_eq!(written.to_toml(), text);
    }

    /// **Each of the numbers a run scores with arrives where calling looks for it**, rather than
    /// merely surviving the trip back out to a file — the round trip above would pass on a
    /// projection that put every value in the wrong field consistently.
    #[test]
    fn the_projected_parameters_are_what_the_file_says() {
        let file = a_file_using_every_shape();
        let run = file
            .to_run_parameters()
            .expect("the fixture is a usable file")
            .parameters;
        let view = run.view();

        assert_eq!(run.ploidy().get(), file.ploidy);
        assert_eq!(run.read_group_count(), 3);
        // The multiplier and its warrant, per read group — including the defaulted one, whose
        // scale of exactly 1.0 is a legitimate fitted answer as well as the default (spec §3.3).
        assert_eq!(view.calibration_by_read_group()[0].scale, 1.0324);
        assert_eq!(
            view.calibration_by_read_group()[0].provenance,
            Provenance::FittedHere
        );
        assert_eq!(view.calibration_by_read_group()[1].scale, 1.0);
        assert_eq!(
            view.calibration_by_read_group()[1].provenance,
            Provenance::Defaulted
        );

        // The contamination axis is dense and its second row was never measured — a fraction of
        // zero beside two zero counts, which `was_measured` tells from measured-and-clean.
        assert!(!view.contamination_is_absent());
        assert_eq!(view.contamination_by_read_group().len(), 3);
        assert_eq!(view.contamination_by_read_group()[0].fraction, 0.031);
        assert!(view.contamination_by_read_group()[0].was_measured());
        assert!(!view.contamination_by_read_group()[1].was_measured());
        assert!(view.contamination_by_read_group()[2].was_measured());

        assert_eq!(
            run.inbreeding_coefficient_by_sample()
                .iter()
                .map(|f| f.get())
                .collect::<Vec<_>>(),
            vec![0.42, 0.17]
        );
        assert_eq!(
            run.prior_seed().alpha_ref(),
            file.ordinary_site_prior.reference_concentration
        );
        assert_eq!(
            run.repeat_tract_outlier_weight().value(),
            DEFAULT_OUTLIER_WEIGHT
        );
        assert_eq!(
            run.repeat_tract_outlier_weight().provenance(),
            Provenance::Defaulted
        );

        // The slippage lookup answers where the file has a row and says *no such stratum* where
        // it does not — never a zero slip rate (spec §5's fifth row).
        let fits = run.ssr_slippage_fits();
        assert_eq!(fits.slippage_group_of(ReadGroupId(2)), Some(1));
        assert_eq!(fits.stated_concentration(), 3.5);
        assert_eq!(fits.strata_with_a_length_spectrum(), 1);
        assert_eq!(fits.periods_with_a_pooled_length_spectrum(), 1);
        assert_eq!(fits.strata(), 3);
        assert_eq!(run.ssr_substitution_rate().count(), 1);
    }

    /// **A supplied outlier weight comes back supplied**, which is the whole of spec §3.8: the
    /// number is in the file so a person can change it, and a run that read a changed one has to
    /// score under it and say so.
    #[test]
    fn an_edited_outlier_weight_arrives_supplied() {
        let mut file = a_file_using_every_shape();
        file.stated_constants.repeat_tract_outlier_weight.value = 0.04;
        file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Supplied;

        let run = file.to_run_parameters().expect("a usable file").parameters;
        assert_eq!(run.repeat_tract_outlier_weight().value(), 0.04);
        assert_eq!(
            run.repeat_tract_outlier_weight().provenance(),
            Provenance::Supplied
        );
        assert_eq!(run.view().repeat_tract_outlier_weight().value(), 0.04);
    }

    /// **An absent `[contamination]` section is the uncontaminated run**, and the read
    /// likelihood then takes its plain formula rather than a mixture whose every share is zero
    /// (spec §5's first row). Collapsing the two changes what is computed at every locus.
    #[test]
    fn an_absent_contamination_table_gives_an_uncontaminated_run() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;

        let run = file.to_run_parameters().expect("a usable file").parameters;
        assert!(run.contamination_by_read_group().is_empty());
        assert!(
            run.view().contamination_is_absent(),
            "the view calling borrows takes the plain formula, not a mixture of zeros"
        );
    }

    /// **The bottom of the committed input range projects** — one sample, one library, no
    /// contamination, no repeat tracts at all.
    ///
    /// `CLAUDE.md` puts a single low-coverage sample first among the cases a design has to have
    /// an answer for, and it is the shape every other test here is furthest from: the fixture has
    /// three read groups over two samples with four strata between them. What it exercises that
    /// nothing else does is the **empty** repeat-tract section — no slippage groups declared, so
    /// the per-stratum row width is zero and no row is ever allocated; no stratum fitted, so the
    /// ladder's bottom rung is the compiled-in constant; and an absent contamination table, so
    /// the run takes the read likelihood's plain formula.
    #[test]
    fn the_smallest_run_the_caller_commits_to_projects() {
        let mut small = a_file_using_every_shape();
        small.contamination = None;
        small.fitted_from.samples.truncate(1);
        small.fitted_from.read_groups.truncate(1);
        small.base_quality_calibration.by_read_group.truncate(1);
        small.sequencing_batches.by_read_group.truncate(1);
        small.sequencing_batches.by_sample.truncate(1);
        small.inbreeding.by_sample.truncate(1);
        small.repeat_tracts.slippage_group_by_read_group.clear();
        small.repeat_tracts.slippage_by_stratum_and_group.clear();
        small.repeat_tracts.length_spectrum_by_stratum.clear();
        small.repeat_tracts.length_spectrum_by_period.clear();
        small.repeat_tracts.substitution_rate_by_stratum.clear();
        small
            .repeat_tracts
            .fallback_length_spectrum_concentration
            .warrant = Warrant::Defaulted;

        let run = small
            .to_run_parameters()
            .expect("one sample with no repeat tracts is a run this caller commits to")
            .parameters;
        assert_eq!(run.read_group_count(), 1);
        assert_eq!(run.inbreeding_coefficient_by_sample().len(), 1);
        assert!(run.view().contamination_is_absent());
        assert_eq!(run.ssr_slippage_fits().strata(), 0);
        assert_eq!(
            run.ssr_slippage_fits().slippage_group_of(ReadGroupId(0)),
            None,
            "a run with no repeat tracts declares no slippage group, and the lookup says so \
             rather than answering under a group that does not exist"
        );
        assert_eq!(run.ssr_substitution_rate().count(), 0);
    }

    /// **The projection refuses before it projects.** `validate` is not folded into `from_toml`,
    /// so this is the entry point that runs it — and until it existed nothing in a run did.
    #[test]
    fn a_file_that_means_nothing_is_refused_rather_than_projected() {
        let mut file = a_file_using_every_shape();
        file.inbreeding.by_sample[0].inbreeding_coefficient.value = 1.7;

        let error = file
            .to_run_parameters()
            .expect_err("an inbreeding coefficient of 1.7 is not a fraction");
        let ParametersFileError::Meaningless { field, .. } = &error else {
            panic!("a file that parses and means nothing is `Meaningless`, not {error:?}")
        };
        assert!(field.starts_with("inbreeding.by_sample"), "{field}");
    }

    /// **A batching in which one sample's two lanes ran in different batches is refused**, and
    /// this fixture is where that came from: it wrote exactly that state until 2026-08-30.
    ///
    /// `SequencingBatches::declared` refuses it in memory — a sample's libraries all ran in one
    /// batch, because the batch is the population a contaminating read is drawn from — so a file
    /// saying otherwise is one no run could have produced. Its symptom is not a crash: the
    /// sample's contaminant genotype would be drawn against one batch's frequencies while one of
    /// its libraries' reads were scored against another's.
    #[test]
    fn a_sample_sequenced_in_two_batches_is_refused() {
        let mut file = a_file_using_every_shape();
        assert_eq!(
            file.fitted_from.read_groups[0].sample, file.fitted_from.read_groups[1].sample,
            "the fixture's first two read groups are two lanes of one plant, which is the case \
             this refusal exists for"
        );
        file.sequencing_batches.by_read_group[1].batch = 1;

        let error = file
            .to_run_parameters()
            .expect_err("one sample cannot be in two batches");
        let ParametersFileError::Meaningless { field, problem } = &error else {
            panic!("{error:?}")
        };
        assert_eq!(field, "sequencing_batches.by_read_group[read_group = 1]");
        assert!(problem.contains("TS-1"), "it names the sample: {problem}");
    }

    /// **The three newtypes whose rules no range check states**, each refused naming its key
    /// rather than panicking several frames later.
    ///
    /// These are the ones `validate` cannot express: a motif period, a repeat count and a ploidy
    /// are each held in a type with its own bound, and the file writes all three as bare
    /// integers a person can type anything into.
    #[test]
    fn a_value_a_newtype_refuses_is_refused_naming_its_key() {
        let refused = |edit: fn(&mut ParametersFile)| -> (String, String) {
            let mut file = a_file_using_every_shape();
            edit(&mut file);
            match file.to_run_parameters() {
                Err(ParametersFileError::Meaningless { field, problem }) => (field, problem),
                other => panic!("expected a refusal, got {:?}", other.map(|_| "a run")),
            }
        };

        // **A ploidy of zero is `validate`'s** and is refused before this walk runs; every
        // other `u8` is a ploidy this caller can tabulate, so the top-level key has no refusal
        // of its own to reach. The rows' ploidy does: nothing checks it earlier.
        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0].ploidy = 0;
        });
        assert!(field.ends_with(".ploidy"), "{field}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0].period = u8::MAX;
        });
        assert!(field.ends_with(".period"), "{field}");

        let (field, _) = refused(|file| {
            file.repeat_tracts.substitution_rate_by_stratum[0].reference_repeats = u64::MAX / 2;
        });
        assert!(field.ends_with(".reference_repeats"), "{field}");
    }

    /// **Every word the file spells reads back to its own word upstream, pair by pair.**
    ///
    /// The eight conversions each way are exhaustive matches, so a variant *added* to either
    /// vocabulary is a compile error. **Exhaustiveness cannot catch a variant crossed with
    /// another** — mapping `flat` to `Sloping` compiles and spells correctly — and that is what
    /// this names. It is the mirror of `from_run_parameters`'s
    /// `every_pre_pass_word_maps_to_its_own_word_in_the_file`, and it is a separate test rather
    /// than a round trip through both because a round trip cannot see two crossings that undo
    /// each other.
    ///
    /// **It asserts one direction, and the sibling asserts the other.** Every loop below reads
    /// `Upstream::from(the file's word)`; the file-ward direction is `from_run_parameters`'s
    /// `every_pre_pass_word_maps_to_its_own_word_in_the_file`. Two tests rather than one, because
    /// a test that checked both at once could pass on a pair crossed the same way twice.
    #[test]
    fn every_word_the_file_spells_reads_back_to_its_own() {
        use super::super::{
            ContaminationFittedFrom, CurveReach, SeedRung, ShareCurveRung, ShareShape,
        };
        use crate::ng::calling::genotype_prior::SeedRegime;
        use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
        use crate::ng::parameter_estimation::joint::share_curve::{
            ShareCurveSource, ShareShape as FittedShape,
        };
        use crate::ng::parameter_estimation::joint::slippage_curve::CurveReach as FittedCurveReach;

        for (upstream, in_the_file) in [
            (Provenance::FittedHere, Warrant::FittedHere),
            (Provenance::Borrowed, Warrant::Borrowed),
            (Provenance::Supplied, Warrant::Supplied),
            (Provenance::Defaulted, Warrant::Defaulted),
        ] {
            assert_eq!(Provenance::from(in_the_file), upstream, "{in_the_file:?}");
        }
        for (upstream, in_the_file) in [
            (
                ContaminationSource::ThisReadGroupsReads,
                ContaminationFittedFrom::ThisReadGroupsOwnReads,
            ),
            (
                ContaminationSource::TheWholeSamplesReads,
                ContaminationFittedFrom::EveryReadOfThisSample,
            ),
        ] {
            assert_eq!(
                ContaminationSource::from(in_the_file),
                upstream,
                "{in_the_file:?}"
            );
        }
        for (upstream, in_the_file) in [
            (SeedRegime::FittedCurve, SeedRung::FittedCurve),
            (SeedRegime::NeutralShape, SeedRung::NeutralShape),
            (SeedRegime::ZeroDiversity, SeedRung::ZeroDiversity),
            (
                SeedRegime::FallbackDiversity,
                SeedRung::StatedHeterozygosity,
            ),
        ] {
            assert_eq!(SeedRegime::from(in_the_file), upstream, "{in_the_file:?}");
        }
        for (upstream, in_the_file) in [
            (FittedCurveReach::Inside, CurveReach::InsideTheFittedRange),
            (
                FittedCurveReach::BelowFitted,
                CurveReach::BelowTheFittedRange,
            ),
            (
                FittedCurveReach::AboveFitted,
                CurveReach::AboveTheFittedRange,
            ),
        ] {
            assert_eq!(
                FittedCurveReach::from(in_the_file),
                upstream,
                "{in_the_file:?}"
            );
        }
        for (upstream, in_the_file) in [
            (ShareCurveSource::ThisPeriod, ShareCurveRung::ThisPeriod),
            (
                ShareCurveSource::ThisPeriodUnscored,
                ShareCurveRung::ThisPeriodUnscored,
            ),
            (ShareCurveSource::OtherPeriods, ShareCurveRung::OtherPeriods),
            (
                ShareCurveSource::BuiltInDefault,
                ShareCurveRung::BuiltInDefault,
            ),
        ] {
            assert_eq!(
                ShareCurveSource::from(in_the_file),
                upstream,
                "{in_the_file:?}"
            );
        }
        for (upstream, in_the_file) in [
            (FittedShape::Flat, ShareShape::Flat),
            (FittedShape::Sloping, ShareShape::Sloping),
            (FittedShape::Turning, ShareShape::Turning),
        ] {
            assert_eq!(FittedShape::from(in_the_file), upstream, "{in_the_file:?}");
        }
    }
}

#[cfg(test)]
mod the_north_star_round_trip {
    //! **Step C4: a run's parameters, written and read back, equal field for field.**
    //!
    //! Spec §1.2 goal 1 and §13's first test, and the plan calls it "the test the whole design
    //! rests on". What it holds is that direct mode and psp mode score the same reads under the
    //! same numbers: the two-mode oracle compares VCFs, so a parameter that comes back changed
    //! shows up there as a changed genotype at some locus, with nothing to say why.
    //!
    //! # What this is built on, and what is owed
    //!
    //! **The plan asks for a `RunParameters` assembled from the joint fit on real tomato records.
    //! This is not that**, by the owner's ruling of 2026-08-30, and the reason is that no program
    //! in this tree produces one: `examples/ng_joint_records_walk.rs` stops at a `JointFit` and
    //! never calls `RunParameters::assemble`; all thirty-three calls to `assemble` are inside
    //! `#[cfg(test)]` modules; and `assemble` needs the per-read-group minted-error totals, which the joint route
    //! does not fit — they come from the generic calibration pre-pass, whose only real-data
    //! program is `examples/ng_minted_error_means.rs`. Joining the two is a new program.
    //!
    //! **So what is built here is a run shaped like a fit's output**: the maps keyed as the
    //! pre-pass keys them, every one of them through the fit's own doors — `StratumFits::over`
    //! takes `StratumOutcome`s, `RunParameters::assemble` takes the raw per-read-group maps — with
    //! many strata carrying mixed origins and the substitution rates at the
    //! `(read group × stratum × ploidy)` grain §9 prices. **⚑ The real-fit round trip is owed**,
    //! and it is owed a program rather than a test.
    //!
    //! # Why the comparison is on the file and not on `RunParameters`
    //!
    //! `RunParameters` has no `PartialEq` — it holds a borrowed-from bundle nothing else compares
    //! — so "equal field for field" is asserted two ways instead. **The file it writes the second
    //! time is byte-identical to the first**, which covers every number the file carries; and
    //! every accessor the run exposes is compared before against after.
    //!
    //! **What the second of those adds is the lookup**, which a file comparison cannot reach at
    //! all: `StratumFits::at`'s four answers for a cell with no numbers, the three rungs of
    //! `length_spectrum_at`, and the bottom rung's median. None of that is written down anywhere
    //! in the file — the file carries the rows the lookup is built from, and the lookup is what
    //! a locus asks.
    //!
    //! **Two things the trip does not cover, and neither is this step's.** The identity block —
    //! the digest, the census, the sample list and the read-group table — is fed back into the
    //! second write rather than round-tripped, because `to_run_parameters` deliberately carries
    //! none of spec §6's four bindings; they are step D2's. And **the tract ladder's fallback
    //! concentration carries a warrant the writer re-derives**, so a projection that lost it
    //! would be invisible here; `validate` refuses the contradiction instead.

    use std::collections::BTreeMap;

    use super::super::tests::a_file_using_every_shape;
    use super::super::{CensusIdentity, CensusTerm, ParametersFile};
    // **The one helper both modules need**, rather than a second copy of it: the rates `of_run`
    // reads back are rebuilt the same way whichever fixture they came from.
    use super::tests::the_rates_the_projection_out_reads;
    use super::*;
    use crate::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
    use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
    use crate::ng::parameter_estimation::joint::contamination::{
        ContaminationEstimate, ContaminationSource, NotIdentifiedReason,
    };
    use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches as DeclaredBatches;
    use crate::ng::parameter_estimation::joint::share_curve::{
        ShareCurve as FittedShareCurve, ShareCurveSource, ShareShape as FittedShape, ShareSource,
    };
    use crate::ng::parameter_estimation::joint::slippage_curve::{
        LevelSource, RiseShape, SlippageCurve as FittedSlippageCurve,
    };
    use crate::ng::parameter_estimation::joint::ssr_fit::{
        DerivedStratum, LevelProvenance, ShareProvenance, SharesProvenance, Slippage, StratumFit,
        StratumOutcome, StratumRefusal,
    };
    use crate::ng::parameter_estimation::joint::stratum_fits::{LengthSpectrumRung, NoSlippage};
    use crate::ng::read::input::read_groups::ReadGroups;
    use crate::ng::types::{Ploidy, SsrPeriod};

    /// **How many samples, libraries, periods and repeat counts the fixture run has.**
    ///
    /// Small enough to run in a unit test and large enough that every table of the file has more
    /// than one row: 3 samples × 2 libraries is 6 read groups, 3 periods × 12 repeat counts is 36
    /// strata, and the substitution rate is keyed by all three of read group, stratum and ploidy.
    /// **The repeat-tract section comes to 132 rows** — 6 slippage-group, 36 slippage, 12
    /// stratum spectra, 2 period pools and 76 substitution rates — which is measured by
    /// `the_files_two_cohort_sized_rows_are_the_size_the_spec_prices_them_at` rather than
    /// asserted here.
    const SAMPLES: usize = 3;
    const LIBRARIES_A_SAMPLE: usize = 2;
    const PERIODS: [u8; 3] = [1, 2, 3];
    const REPEAT_COUNTS: std::ops::RangeInclusive<u64> = 5..=16;
    const SLIPPAGE_GROUPS: usize = 2;

    /// A number that differs at every index, so no two cells of the file hold the same value.
    ///
    /// **Two rows carrying one number is how a projection that transposed an axis passes**: the
    /// file would be written wrong and read back wrong the same way. So this mixes its indices
    /// rather than adding them — a first draft added them with small weights and the 72
    /// substitution rows then carried only 53 distinct values, because `8Δg + 15Δp + 22Δr = 0`
    /// has a solution inside the fixture's own index ranges.
    ///
    /// **Its range is `[scale, 2·scale)`**, so a caller picks the magnitude and the mixing picks
    /// the digits; every value the fixture uses as a share is called with a scale that keeps it
    /// inside the range its key allows, and no clamp is needed.
    fn a_number(scale: f64, of: &[usize]) -> f64 {
        let mut mixed: u64 = 0x9e37_79b9_7f4a_7c15;
        for index in of {
            mixed ^= *index as u64 + 1;
            mixed = mixed.wrapping_mul(0x1000_0000_01b3);
            mixed ^= mixed >> 29;
        }
        // The low 40 bits as a fraction, so the value is spread over the whole of `[1, 2)` and
        // two different index tuples collide only if their mixes agree in 40 bits.
        scale * (1.0 + (mixed & 0xff_ffff_ffff) as f64 / (1u64 << 40) as f64)
    }

    fn the_runs_read_groups() -> ReadGroups {
        let lanes: Vec<(String, String, String)> = (0..SAMPLES)
            .flat_map(|sample| {
                (0..LIBRARIES_A_SAMPLE).map(move |library| {
                    (
                        format!("HWI.{sample}.{library}"),
                        format!("TS-{sample}"),
                        format!("lib{sample}-{library}"),
                    )
                })
            })
            .collect();
        let borrowed: Vec<(&str, &str, &str)> = lanes
            .iter()
            .map(|(id, sample, library)| (id.as_str(), sample.as_str(), library.as_str()))
            .collect();
        ReadGroups::of_lanes(&borrowed)
    }

    fn read_group_count() -> usize {
        SAMPLES * LIBRARIES_A_SAMPLE
    }

    /// The two per-read-group maps assembly takes, keyed as the pre-pass keys them.
    ///
    /// **One read group's fitted rate is zero**, which is how assembly reaches
    /// `ReadGroupCalibration::defaulted` — a multiplier of exactly 1.0 whose warrant is
    /// `Defaulted`, spec §5's third row and the one a reader must not read as a fitted answer.
    /// A zero rate is a probability and the fit can emit one; what `from_fitted_rate` refuses is
    /// the *scale* it would give, which charges every read of the library the error floor.
    /// **The others' multipliers differ from one another**, because the minted mean is built from
    /// a different rate than the one fitted: a fixture whose minted total is `rate.ln() · reads`
    /// gives every library a scale of exactly one, which is the column where a transposed axis
    /// would hide best.
    #[allow(clippy::type_complexity)]
    fn the_runs_calibration_inputs() -> (
        BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        BTreeMap<ReadGroupId, MintedReadErrors>,
    ) {
        let mut rates = BTreeMap::new();
        let mut minted = BTreeMap::new();
        for group in 0..read_group_count() {
            let id = ReadGroupId(group as u32);
            // **Both maps hold every read group**, because assembly refuses a gap in either:
            // the axis is `0..n` with nothing missing. What makes one library defaulted is the
            // *value*, not an absence.
            let rate = if group == THE_LIBRARY_NOTHING_WAS_FITTED_FOR {
                0.0
            } else {
                a_number(0.002, &[group])
            };
            rates.insert(
                id,
                Estimate {
                    value: ErrorRate::try_new(rate).expect("a rate"),
                    // **Three warrants across the libraries**, so the file carries more than one
                    // and a projection that hard-coded any of them fails.
                    provenance: match group % 3 {
                        0 => Provenance::FittedHere,
                        1 => Provenance::Borrowed,
                        _ => Provenance::Supplied,
                    },
                    observations: 100_000 + group as u64 * 7_919,
                },
            );
            let reads = 1_000 + group as u32 * 13;
            // The instrument's own mean, which is not the fitted rate — so the scale between them
            // is a number of its own rather than exactly one at every library.
            let reported = a_number(0.0025, &[group, 1]);
            minted.insert(
                id,
                MintedReadErrors::of_observation(reported.ln() * f64::from(reads), reads),
            );
        }
        (rates, minted)
    }

    /// **The one library the calibration fit produced nothing for**, so that the run carries a
    /// defaulted scale beside fitted ones.
    const THE_LIBRARY_NOTHING_WAS_FITTED_FOR: usize = 4;

    /// **The one sample whose contamination could not be identified**, and it is a *sample*
    /// rather than a library: the only per-unit refusal the contamination fit has is
    /// `OwnFrequencyIsItsOwnEcho`, whose own documentation says "the refusal is the sample's, and
    /// it refuses every library the sample has". Its other refusal stamps the whole run.
    const THE_SAMPLE_THAT_IS_ITS_OWN_ECHO: usize = 2;

    /// **Contamination measured on some samples and not on others**, and inside the measured ones
    /// one library measured and found clean — a zero fraction with evidence behind it, which is
    /// spec §5's second row and not the same claim as unmeasured.
    fn the_runs_contamination() -> BTreeMap<ReadGroupId, ContaminationEstimate> {
        (0..read_group_count())
            .map(|group| {
                let id = ReadGroupId(group as u32);
                let sample = group / LIBRARIES_A_SAMPLE;
                let estimate = if sample == THE_SAMPLE_THAT_IS_ITS_OWN_ECHO {
                    ContaminationEstimate::NotIdentified {
                        reason: NotIdentifiedReason::OwnFrequencyIsItsOwnEcho,
                    }
                } else {
                    ContaminationEstimate::Estimated {
                        alpha: if group % 3 == 1 {
                            0.0
                        } else {
                            a_number(0.01, &[group])
                        },
                        source: ContaminationSource::ThisReadGroupsReads,
                        panel_markers: 10_000,
                        markers_with_reads: 4_000 + group as u64,
                        reads_on_markers: 90_000 + group as u64 * 11,
                        // Below `MAX_LEVERAGE`, above which the fit refuses the estimate rather
                        // than emitting it.
                        leverage: 0.4,
                    }
                };
                (id, estimate)
            })
            .collect()
    }

    fn the_runs_inbreeding() -> Vec<Estimate<InbreedingF>> {
        (0..SAMPLES)
            .map(|sample| Estimate {
                value: InbreedingF::try_new(a_number(0.1, &[sample])).expect("a coefficient"),
                provenance: if sample % 2 == 0 {
                    Provenance::FittedHere
                } else {
                    Provenance::Borrowed
                },
                observations: 180_000_000 + sample as u64 * 1_009,
            })
            .collect()
    }

    /// **Which periods the run fitted curves for.** A period with no curve is a real state — it
    /// is what `LevelSource::Cell` and `ShareSource::Stratum` mean, in their own words — and it
    /// is also what decides which outcomes a period's strata can have: `derive_thin_strata`
    /// produces a derived stratum **only** where its period has a curve, so a period without one
    /// holds fitted strata and refusals and nothing else.
    fn the_period_has_curves(period: u8) -> bool {
        period != PERIODS[0]
    }

    fn a_slippage_curve_over(from: u64, to: u64, seed: usize) -> FittedSlippageCurve {
        FittedSlippageCurve {
            rise_shape: RiseShape::new(a_number(0.4, &[seed])).expect("a rise shape"),
            intercept: a_number(0.01, &[seed]),
            slope: a_number(0.004, &[seed, 1]),
            fitted_from: from,
            fitted_to: to,
            held_out_error: a_number(0.3, &[seed, 2]),
            cells: (to - from + 1) as usize,
        }
    }

    fn a_share_curve_over(from: u64, to: u64, seed: usize) -> FittedShareCurve {
        FittedShareCurve {
            shape: FittedShape::Turning,
            intercept: a_number(1.4, &[seed]),
            slope: -a_number(0.09, &[seed, 1]),
            bend: a_number(0.006, &[seed, 2]),
            centre: (from + to) as f64 / 2.0,
            fitted_from: from,
            fitted_to: to,
            held_out_error: a_number(0.167, &[seed, 3]),
            strata: (to - from + 1) as usize,
            source: ShareCurveSource::ThisPeriod,
        }
    }

    /// **The fit's own outcomes, mixed, and every one of them a state the fit can reach.**
    ///
    /// The mixture is by period and by repeat count, and the two are kept apart deliberately:
    ///
    /// - **the first period has no curves**, so its strata are fitted on their own tracts or
    ///   refused for want of a spanning read, its levels come from the cell and its shares from
    ///   the stratum — which is what `LevelSource::Cell` and `ShareSource::Stratum` mean;
    /// - **the other two have both curves**, so a stratum with its own evidence blends
    ///   (`LevelSource::Blend`, the fit's ordinary case) and one without is derived whole from
    ///   them (`LevelSource::Curve`, both shares `ShareSource::Curve` and no slipped-read count,
    ///   which is exactly what `derive_thin_strata` writes).
    ///
    /// **A refused stratum is `NoSpanningReads` everywhere**, because the other refusal —
    /// `BelowTheFloor` — is what `derive_thin_strata` converts into a derived stratum wherever a
    /// curve exists, so it cannot survive in a period that has one.
    fn the_runs_outcomes() -> Vec<StratumOutcome> {
        let mut outcomes = Vec::new();
        for (which_period, period) in PERIODS.iter().enumerate() {
            // **One curve a period**, which is what the fit produces and what the file's own
            // section comment says.
            let from = *REPEAT_COUNTS.start() + which_period as u64;
            let to = *REPEAT_COUNTS.end() - which_period as u64;
            let has_curves = the_period_has_curves(*period);
            let level_curve = has_curves.then(|| a_slippage_curve_over(from, to, which_period));
            let share_curve = has_curves.then(|| a_share_curve_over(from, to, which_period));
            for (which_repeats, repeats) in REPEAT_COUNTS.enumerate() {
                let stratum = Stratum {
                    period: *period,
                    reference_repeats: repeats,
                };
                let keys = [which_period, which_repeats];
                // **What this stratum is**: fitted on its own tracts, derived whole from its
                // period's curves — only possible where there are curves — or refused.
                let fitted_here = which_repeats % 3 == 0 || (!has_curves && which_repeats % 3 == 1);
                let refused = which_repeats % 3 == 2;
                let slippage: Vec<Option<Slippage>> = (0..SLIPPAGE_GROUPS)
                    .map(|group| {
                        // **A pair with no row**: the second group is silent at every third
                        // stratum, which is spec §5's fifth state.
                        (group == 0 || which_repeats % 3 != 0).then(|| Slippage {
                            level: a_number(0.04, &[which_period, which_repeats, group]),
                            shorter_share: a_number(0.4, &[which_period, which_repeats, group, 1]),
                            fall_off: a_number(0.2, &[which_period, which_repeats, group, 2]),
                        })
                    })
                    .collect();
                let level: Vec<Option<LevelProvenance>> = slippage
                    .iter()
                    .enumerate()
                    .map(|(group, fitted)| {
                        fitted.map(|_| match (&level_curve, fitted_here) {
                            // Its period had no curve at all, so the cell keeps its own answer.
                            (None, _) => LevelProvenance {
                                source: LevelSource::Cell,
                                curve: None,
                                reach: None,
                                slipped_reads: Some(a_number(8_000.0, &keys)),
                            },
                            // Its own evidence weighed against its period's curve.
                            (Some(curve), true) => LevelProvenance {
                                source: LevelSource::Blend {
                                    curve_weight: a_number(0.3, &[group, which_repeats]),
                                },
                                curve: Some(*curve),
                                reach: Some(curve.reach(repeats)),
                                slipped_reads: Some(a_number(8_000.0, &keys)),
                            },
                            // Nothing of its own: the curve, whole.
                            (Some(curve), false) => LevelProvenance {
                                source: LevelSource::Curve,
                                curve: Some(*curve),
                                reach: Some(curve.reach(repeats)),
                                slipped_reads: None,
                            },
                        })
                    })
                    .collect();
                // **A share's origin, and the curve it names where it names one.** A source of
                // `Stratum` means the period had no curve at all, so it carries neither the curve
                // nor a reach — which is what `blend_share` writes in that arm.
                let a_share = |source: ShareSource| match source {
                    ShareSource::Stratum => ShareProvenance {
                        source,
                        curve: None,
                        reach: None,
                    },
                    _ => ShareProvenance {
                        source,
                        curve: share_curve,
                        reach: share_curve.map(|curve| curve.reach(repeats)),
                    },
                };
                let shares: Vec<Option<SharesProvenance>> = slippage
                    .iter()
                    .enumerate()
                    .map(|(group, fitted)| {
                        // **A row with no shares provenance at all**, at one stratum in four:
                        // the `shorter_share_and_fall_off_origin` key the file leaves out.
                        fitted.filter(|_| which_repeats % 4 != 3).map(|_| {
                            match (&share_curve, fitted_here) {
                                (None, _) => SharesProvenance {
                                    slipped_reads: Some(a_number(
                                        8_000.0,
                                        &[which_period, which_repeats, group],
                                    )),
                                    shorter_share: a_share(ShareSource::Stratum),
                                    fall_off: a_share(ShareSource::Stratum),
                                },
                                (Some(_), true) => SharesProvenance {
                                    slipped_reads: Some(a_number(
                                        8_000.0,
                                        &[which_period, which_repeats, group],
                                    )),
                                    shorter_share: a_share(ShareSource::Blend {
                                        curve_weight: a_number(0.45, &[which_repeats, group]),
                                    }),
                                    fall_off: a_share(ShareSource::Curve),
                                },
                                // What `derive_thin_strata` writes, and it writes nothing else.
                                (Some(_), false) => SharesProvenance {
                                    slipped_reads: None,
                                    shorter_share: a_share(ShareSource::Curve),
                                    fall_off: a_share(ShareSource::Curve),
                                },
                            }
                        })
                    })
                    .collect();

                if refused {
                    // **`NoSpanningReads`, not `BelowTheFloor`**: the second is what a period's
                    // curves turn into a derived stratum, so it cannot survive beside them.
                    outcomes.push(StratumOutcome::Refused {
                        stratum,
                        tracts: 0,
                        reason: StratumRefusal::NoSpanningReads,
                    });
                } else if fitted_here {
                    let span = 1 + which_period;
                    let classes = 2 * span + 1;
                    let mut weights: Vec<f64> = (0..classes)
                        .map(|class| a_number(1.0, &[which_period, which_repeats, class]))
                        .collect();
                    // Normalised exactly, since `StratumFits::over` refuses a spectrum that is
                    // not a distribution.
                    let total: f64 = weights.iter().sum();
                    for weight in &mut weights {
                        *weight /= total;
                    }
                    outcomes.push(StratumOutcome::Fitted(Box::new(StratumFit {
                        stratum,
                        slippage,
                        length_spectrum: weights,
                        concentration: a_number(3.0, &keys),
                        log_likelihood_a_tract: -a_number(2.0, &keys),
                        tracts_fitted: 100 + which_repeats,
                        // **Empty, which is what a stratum that stood on its own tracts has** —
                        // this field names the neighbouring repeat counts it borrowed from.
                        borrowed: Vec::new(),
                        converged: true,
                        tracts_of_its_own: 100 + which_repeats,
                        reads_crossing: 5_000 + which_repeats as u64,
                        level_provenance: level,
                        shares_provenance: shares,
                    })));
                } else {
                    outcomes.push(StratumOutcome::Derived(Box::new(DerivedStratum {
                        stratum,
                        slippage,
                        level_provenance: level,
                        shares_provenance: shares,
                        tracts_of_its_own: 3,
                        reads_crossing: 40,
                    })));
                }
            }
        }
        outcomes
    }

    /// **The tract ladder's middle rung**, one pooled spectrum a period, for the periods whose
    /// curves the fit drew — the rung a run that skips the second pass over its tracts does not
    /// have, and the one the file writes in `length_spectrum_by_period`.
    fn the_runs_period_pools()
    -> BTreeMap<u8, crate::ng::parameter_estimation::joint::ssr_fit::PeriodLengthSpectrum> {
        PERIODS
            .iter()
            .enumerate()
            .filter(|(_, period)| the_period_has_curves(**period))
            .map(|(which_period, period)| {
                let classes = 2 * (1 + which_period) + 1;
                let mut weights: Vec<f64> = (0..classes)
                    .map(|class| a_number(1.0, &[which_period, class, 7]))
                    .collect();
                let total: f64 = weights.iter().sum();
                for weight in &mut weights {
                    *weight /= total;
                }
                (
                    *period,
                    crate::ng::parameter_estimation::joint::ssr_fit::PeriodLengthSpectrum {
                        period: *period,
                        length_spectrum: weights,
                        concentration: a_number(2.5, &[which_period]),
                        tracts_fitted: 400 + which_period,
                        strata_pooled: REPEAT_COUNTS.count(),
                        converged: true,
                    },
                )
            })
            .collect()
    }

    /// The substitution rate at the grain §9 prices it at: one row a
    /// `(read group × stratum × ploidy)` that was fitted, and none where it was not.
    fn the_runs_substitution_rates() -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
        let mut rates = BTreeMap::new();
        for group in 0..read_group_count() {
            for (which_period, period) in PERIODS.iter().enumerate() {
                for (which_repeats, repeats) in REPEAT_COUNTS.enumerate() {
                    // **A row only where one was fitted**, which is what makes this axis's size a
                    // question about the cohort rather than about the shape (spec §9).
                    if (group + which_period + which_repeats) % 3 != 0 {
                        continue;
                    }
                    let key = StratumKey {
                        read_group: ReadGroupId(group as u32),
                        stratum: SsrStratum::new(
                            SsrPeriod::try_new(usize::from(*period)).expect("a period"),
                            RepeatCount(u32::try_from(repeats).expect("a repeat count")),
                        ),
                        ploidy: Ploidy::try_new(2).expect("diploid"),
                    };
                    rates.insert(
                        key,
                        Estimate {
                            value: ErrorRate::try_new(a_number(
                                0.001,
                                &[group, which_period, which_repeats],
                            ))
                            .expect("a rate"),
                            provenance: Provenance::FittedHere,
                            // **A count at the size a real one has**, because the width of this
                            // number is part of what the file's largest axis costs a cohort.
                            observations: 172_000_000 + group as u64 * 1_009,
                        },
                    );
                }
            }
        }
        rates
    }

    /// Everything `of_run` needs, from one place, so the two directions cannot be handed
    /// different inputs.
    struct TheRun {
        parameters: RunParameters,
        read_groups: ReadGroups,
        rates: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        inbreeding: Vec<Estimate<InbreedingF>>,
    }

    fn a_run_shaped_like_a_fits() -> TheRun {
        let read_groups = the_runs_read_groups();
        let (rates, minted) = the_runs_calibration_inputs();
        let inbreeding = the_runs_inbreeding();
        let slippage_group_of: BTreeMap<ReadGroupId, u32> = (0..read_group_count())
            .map(|group| (ReadGroupId(group as u32), (group % SLIPPAGE_GROUPS) as u32))
            .collect();
        let parameters = RunParameters::assemble(
            &rates,
            &minted,
            &the_runs_contamination(),
            DeclaredBatches::all_together(&read_groups),
            inbreeding.iter().map(|estimate| estimate.value).collect(),
            SpectrumSeed::new(1.0, 6e-4, SeedRegime::FittedCurve),
            StratumFits::over(&the_runs_outcomes(), slippage_group_of)
                .with_period_length_spectra(the_runs_period_pools()),
            the_runs_substitution_rates(),
            Ploidy::try_new(2).expect("diploid"),
        );
        TheRun {
            parameters,
            read_groups,
            rates,
            inbreeding,
        }
    }

    fn a_census() -> CensusIdentity {
        CensusIdentity {
            terms: vec![CensusTerm {
                term: "the loci actually kept".into(),
                digest: "fedcba9876543210".into(),
            }],
        }
    }

    const A_REFERENCE_DIGEST: &str = "0123456789abcdef";

    fn written(run: &TheRun) -> ParametersFile {
        ParametersFile::of_run(
            &run.parameters,
            &run.read_groups,
            &run.rates,
            &run.inbreeding,
            A_REFERENCE_DIGEST,
            a_census(),
        )
    }

    /// **The whole trip: parameters, file, text, file, parameters, file — and the two files are
    /// the same file.**
    ///
    /// Six stages rather than three, because three would not see a projection that lost something
    /// on the way *in*: a field dropped by `to_run_parameters` and re-derived by `of_run` gives
    /// the same file after one trip. Writing the file a second time, from the parameters the file
    /// produced, is what makes the middle of the trip visible.
    #[test]
    fn a_run_shaped_like_a_fits_survives_the_whole_trip() {
        let run = a_run_shaped_like_a_fits();
        let first = written(&run);
        first
            .validate()
            .expect("a run this caller can assemble writes a file it can use");

        let text = first.to_toml();
        let read = ParametersFile::from_toml(&text).expect("the writer's own text parses");
        assert_eq!(read, first, "the text is the file");

        let back = read
            .to_run_parameters()
            .expect("and the file is a run's parameters");
        let again = ParametersFile::of_run(
            &back.parameters,
            &run.read_groups,
            &the_rates_the_projection_out_reads(&back),
            &back.inbreeding_by_sample,
            &read.fitted_from.reference_digest,
            read.fitted_from.census.clone(),
        );
        assert_eq!(again, first, "and the parameters are the file again");
        // **Implied by the line above except in one place**: the shape's `PartialEq` compares
        // floats with `==`, and `-0.0 == 0.0` while the two are written differently. That one
        // case is the whole of what this buys.
        assert_eq!(again.to_toml(), text, "byte for byte");
    }

    /// **Every number a locus reads answers the same after the trip**, which the file comparison
    /// cannot see: two strata whose rows were exchanged, numbers and all, write the same file.
    #[test]
    fn every_lookup_a_locus_makes_answers_the_same_after_the_trip() {
        let run = a_run_shaped_like_a_fits();
        let before = &run.parameters;
        let after = written(&run)
            .to_run_parameters()
            .expect("a run's own file is a run's parameters")
            .parameters;

        assert_eq!(after.ploidy(), before.ploidy());
        assert_eq!(after.read_group_count(), before.read_group_count());
        assert_eq!(
            after.calibration_by_read_group(),
            before.calibration_by_read_group()
        );
        assert_eq!(
            after.contamination_by_read_group(),
            before.contamination_by_read_group()
        );
        assert_eq!(
            after.inbreeding_coefficient_by_sample(),
            before.inbreeding_coefficient_by_sample()
        );
        assert_eq!(after.prior_seed(), before.prior_seed());
        assert_eq!(
            after.repeat_tract_outlier_weight(),
            before.repeat_tract_outlier_weight()
        );
        assert_eq!(
            after.ssr_substitution_rate().collect::<Vec<_>>(),
            before.ssr_substitution_rate().collect::<Vec<_>>()
        );

        // **The slippage lookup at every cell a locus can ask about**, which is every read group
        // against every `(period, repeat count)` the run has — and two beyond its range, so that
        // *no such stratum* is compared as well as the answers.
        let mut answered = 0_u32;
        let mut absent = 0_u32;
        for group in 0..read_group_count() {
            let id = ReadGroupId(group as u32);
            for period in PERIODS {
                for repeats in *REPEAT_COUNTS.start() - 1..=*REPEAT_COUNTS.end() + 1 {
                    let was = before.ssr_slippage_fits().at(id, period, repeats);
                    let is = after.ssr_slippage_fits().at(id, period, repeats);
                    assert_eq!(
                        is, was,
                        "read group {group}, period {period}, {repeats} repeats"
                    );
                    match was {
                        Ok(_) => answered += 1,
                        Err(_) => absent += 1,
                    }
                }
            }
        }
        // **A length spectrum at every stratum**, which rides on a different table and does not
        // depend on the read group — so it is walked once rather than once a library.
        let mut rungs: BTreeMap<&str, u32> = BTreeMap::new();
        for period in PERIODS {
            for repeats in REPEAT_COUNTS {
                let was = before
                    .ssr_slippage_fits()
                    .length_spectrum_at(period, repeats);
                let is = after
                    .ssr_slippage_fits()
                    .length_spectrum_at(period, repeats);
                assert_eq!(is.rung(), was.rung(), "period {period}, {repeats} repeats");
                assert_eq!(is.concentration(), was.concentration());
                assert_eq!(is.fitted_weights(), was.fitted_weights());
                *rungs
                    .entry(match was.rung() {
                        LengthSpectrumRung::StratumsOwnFit => "its own",
                        LengthSpectrumRung::PeriodsPooledTracts => "its period's",
                        LengthSpectrumRung::StatedFlat => "the flat one",
                    })
                    .or_default() += 1;
            }
        }
        // **All three rungs of the ladder are walked**, which is what stops this comparing one
        // rung with itself thirty-six times. Measured rather than reasoned about: the 16 strata
        // fitted on their own tracts take the top rung; the 16 of the two curved periods that
        // were not — 8 derived and 8 refused — take their period's pool; and the 4 refused strata
        // of the period that has no curves, and so no pool, fall to the flat one.
        assert_eq!(
            rungs,
            BTreeMap::from([("its own", 16), ("its period's", 16), ("the flat one", 4)])
        );
        // **The counts are asserted so that a lookup answering *absent* everywhere cannot pass.**
        // Measured on this fixture rather than reasoned about: 6 read groups × 3 periods × 14
        // repeat counts is 252 cells, of which **108 carry numbers and 144 do not** — a third of
        // the strata were refused, two of the fourteen repeat counts are outside the run's range
        // altogether, and one slippage group is silent at every third stratum.
        assert_eq!((answered, absent), (108, 144));
        assert_eq!(
            (
                before.ssr_slippage_fits().strata(),
                before.ssr_slippage_fits().strata_with_a_length_spectrum(),
            ),
            (24, 16),
            "a third of the 36 strata were refused and contribute no row; of the 24 that remain, \
             the 16 with evidence of their own were fitted and the 8 without were derived from \
             their period's curves"
        );
    }

    /// **What the file costs at this shape**, which spec §9 prices from rows measured on the
    /// built shape and says nothing has counted on a real cohort.
    ///
    /// **This is a synthetic run and not a cohort**, so what it can say is bytes a row rather
    /// than how many rows a cohort has. It is here because §9's two per-row figures — 146 bytes
    /// an inbreeding row, 146 a substitution-rate row — are the numbers its 0.44 MB and 62 MB
    /// rest on, and nothing had re-measured them since the key names changed.
    #[test]
    fn the_files_two_cohort_sized_rows_are_the_size_the_spec_prices_them_at() {
        let text = written(&a_run_shaped_like_a_fits()).to_toml();
        // **The section as well as the table**, because `by_sample` names two different tables —
        // the batching's and the inbreeding coefficients' — and the batching's comes first.
        let a_row_of = |section: &str, table: &str| -> usize {
            let rows: Vec<&str> = text
                .lines()
                .skip_while(|line| !line.starts_with(section))
                .skip_while(|line| !line.starts_with(table))
                .skip(1)
                .take_while(|line| !line.starts_with(']'))
                // **Rows only.** A `defaulted` row carries a provenance note above it, and a
                // "bytes a row" that averaged prose in would be a different number quietly.
                .filter(|line| !line.trim_start().starts_with('#'))
                .collect();
            assert!(!rows.is_empty(), "{section} has a {table}");
            rows.iter().map(|row| row.len() + 1).sum::<usize>() / rows.len()
        };
        let inbreeding = a_row_of("[inbreeding]", "by_sample = [");
        let substitution = a_row_of("[repeat_tracts]", "substitution_rate_by_stratum = [");
        // **Not an assertion on the spec's own two numbers**, which were measured on a different
        // fixture with different names: what is pinned is the order of magnitude, so that a
        // change which doubled a row's width would be noticed here rather than at 3,000 samples.
        assert!(
            (100..200).contains(&inbreeding),
            "an inbreeding row is {inbreeding} bytes, where spec §9 prices it at 146"
        );
        assert!(
            (100..200).contains(&substitution),
            "a substitution-rate row is {substitution} bytes, where spec §9 prices it at 146"
        );
        // **The two numbers, so a reader has them rather than the band.** ⚑ Both are larger
        // than the 146 bytes spec §9 prices each at, and the substitution rate — the axis §9's
        // 62 MB at 3,000 samples rests on — is larger by a quarter, which would put that figure
        // nearer 79 MB.
        //
        // **The difference is what is in the rows, not the format.** No key either row carries
        // was touched by the 2026-08-30 renames, and neither row's shape has changed since §9 was
        // written; what differs is how wide the numbers and names inside them are. So **185 is a
        // floor and 79 MB with it**, and two things push a real cohort above it: a `read_group`
        // is one digit here and four at 3,000 libraries, and a `bases_compared` count is nine
        // digits here where the fit's own per-read-group total on HG002 is 172,616,054. The
        // inbreeding row carries a sample name as well, and this fixture's are four characters.
        //
        // **§9's own 146 reproduces nothing in the tree**, including the fixture that existed
        // when it was written, so what is compared here is an order of magnitude rather than a
        // measurement against a measurement.
        assert_eq!(
            (inbreeding, substitution),
            (157, 185),
            "the two rows the file's cohort-sized axes are made of"
        );
    }

    /// **A run with no repeat tracts at all makes the file-to-parameters-to-file trip too**,
    /// which is the bottom of the committed input range: `CLAUDE.md` names one sample as the case
    /// a design has to have an answer for, and it is the case every other fixture here is
    /// furthest from.
    ///
    /// **It is a shorter trip than the one above** — file → parameters → file, with no text in it
    /// — because what it exercises is the projection's behaviour on empty tables rather than the
    /// writer's.
    #[test]
    fn the_smallest_run_the_caller_commits_to_makes_the_trip_too() {
        let file = {
            let mut small = a_file_using_every_shape();
            small.contamination = None;
            small.fitted_from.samples.truncate(1);
            small.fitted_from.read_groups.truncate(1);
            small.base_quality_calibration.by_read_group.truncate(1);
            small.sequencing_batches.by_read_group.truncate(1);
            small.sequencing_batches.by_sample.truncate(1);
            small.inbreeding.by_sample.truncate(1);
            small.repeat_tracts.slippage_group_by_read_group.clear();
            small.repeat_tracts.slippage_by_stratum_and_group.clear();
            small.repeat_tracts.length_spectrum_by_stratum.clear();
            small.repeat_tracts.length_spectrum_by_period.clear();
            small.repeat_tracts.substitution_rate_by_stratum.clear();
            small
                .repeat_tracts
                .fallback_length_spectrum_concentration
                .warrant = super::super::Warrant::Defaulted;
            small
        };
        let read_groups = ReadGroups::of_lanes(&[(
            file.fitted_from.read_groups[0].declared_id.as_str(),
            file.fitted_from.read_groups[0].sample.as_str(),
            file.fitted_from.read_groups[0].library.as_str(),
        )]);

        let back = file.to_run_parameters().expect("one sample, no tracts");
        let again = ParametersFile::of_run(
            &back.parameters,
            &read_groups,
            &the_rates_the_projection_out_reads(&back),
            &back.inbreeding_by_sample,
            &file.fitted_from.reference_digest,
            file.fitted_from.census.clone(),
        );
        assert_eq!(again, file);
        assert_eq!(
            back.parameters.ssr_slippage_fits().at(ReadGroupId(0), 2, 6),
            Err(NoSlippage::UnknownReadGroup),
            "a run that declared no slippage group says so, rather than answering under one"
        );
    }
}

#[cfg(test)]
mod the_five_states_survive_the_round_trip {
    //! **Step C5: absent, zero and default are three different claims, and collapsing any two of
    //! them changes an answer.**
    //!
    //! Spec §5 is a table of five rows and one sentence — *`Option<T>` is absence, never a
    //! sentinel, and a warrant is carried rather than inferred from the value*. The plan asks for
    //! a fixture per row **built so that collapsing the two states it separates changes an
    //! answer, not merely so they differ**, because a test that only shows two values are unequal
    //! passes for a reader who has collapsed them into a third thing.
    //!
    //! # The answer is not the same kind of thing for all five rows, and saying which is the
    //! point
    //!
    //! **Two of the five change a number a locus is scored against**, and those are the expensive
    //! ones: a stratum's length spectrum is the prior every tract of it is seeded from, and a
    //! `(stratum × slippage group)` with no row sends the caller to a shipped stutter model where
    //! a zero slip rate would tell it no read of that stratum can ever report a neighbouring
    //! length.
    //!
    //! **Three change what the run says about itself** — the report an output prints, and the
    //! warrant every call resting on that number carries — and change no number it computes. That
    //! is not a weaker finding; it is the finding. A defaulted multiplier of 1.0 and a fitted one
    //! of 1.0 multiply every read's error probability by exactly the same number, and a
    //! contamination fraction of zero gives the three-term read likelihood bit-for-bit what the
    //! two-term one gives (`likelihood::ssr`'s own
    //! `a_contamination_fraction_of_zero_is_the_two_term_form`). So **the warrant and the report
    //! are the only things that tell those states apart, and a reader that inferred either from
    //! the value would have lost the distinction with nothing to notice**. Spec §5's own sentence
    //! says so; these tests are what make it true of the code.
    //!
    //! **What is new here, and what restates a neighbour.** The five tests overlap the module's
    //! other two: the round trips already carry the *shapes* through. What only these have is the
    //! **answer** each collapse changes — the run's report for rows 1 and 2, the warrant a call
    //! folds for row 3, the prior's own numbers for row 4, and the stutter model a read is scored
    //! against for row 5 — and, in row 1, the one file nothing else in the tree builds: a
    //! contamination table of *measured* zeros, which is legal, projects to a mixture, and is the
    //! collapse §5's first row is actually about.

    use super::super::tests::a_file_using_every_shape;
    use super::super::{
        ContaminationFittedFrom, ContaminationMeasurement, ContaminationRow, ParametersFile,
        StratumLengthSpectrumRow, Warrant,
    };
    use super::*;
    // **The one helper both test modules need**, rather than a second copy of it.
    use super::tests::the_files_read_groups;
    use crate::ng::alignment::StutterModel;
    use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::ng::calling::run_report::ContaminationUsed;
    use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;
    use crate::ng::parameter_estimation::joint::stratum_fits::{LengthSpectrumRung, NoSlippage};

    fn a_run_from(file: &ParametersFile) -> RunParameters {
        file.to_run_parameters()
            .unwrap_or_else(|error| panic!("{error}"))
            .parameters
    }

    /// **Row 1. No contamination table at all is an uncontaminated run**, and a table of zeros is
    /// not the same claim.
    ///
    /// **The answer that changes is which formula the read likelihood runs.** An uncontaminated
    /// run computes its plain form; a run carrying a mixture computes the three-term one at every
    /// locus, against a contaminant population drawn from the batching — so a reader that wrote
    /// zeros for the absent table would have every locus of the run scored through a term that
    /// says a share of every read came from somebody else and that share is nothing.
    ///
    /// **And the longhand form is refused outright**, which is the second half: a table in which
    /// no row was measured is the uncontaminated run written out, and `validate` says so rather
    /// than letting it become a run.
    #[test]
    fn an_absent_contamination_table_is_not_a_table_of_zeros() {
        let read_groups = the_files_read_groups(&a_file_using_every_shape());

        let mut absent = a_file_using_every_shape();
        absent.contamination = None;
        let uncontaminated = a_run_from(&absent);
        assert!(uncontaminated.view().contamination_is_absent());
        assert_eq!(
            uncontaminated.report(&read_groups).contamination(),
            &ContaminationUsed::NoneFitted,
            "and the run says so, which is what an output prints beside the calls"
        );

        // The same run with the table written out in full, every fraction zero and every count
        // real — the collapse spec §5's first row forbids.
        let mut longhand = a_file_using_every_shape();
        longhand.contamination = Some(super::super::Contamination {
            by_read_group: absent
                .fitted_from
                .read_groups
                .iter()
                .map(|row| ContaminationRow {
                    read_group: row.read_group,
                    library: row.library.clone(),
                    measurement: Some(ContaminationMeasurement {
                        share_of_reads_from_another_sample: 0.0,
                        markers_with_reads: 4_211,
                        reads_on_markers: 90_233,
                        fitted_from_reads_of: ContaminationFittedFrom::ThisReadGroupsOwnReads,
                    }),
                })
                .collect(),
        });
        let collapsed = a_run_from(&longhand);
        assert!(
            !collapsed.view().contamination_is_absent(),
            "a table of measured zeros is a mixture, and every locus is scored through it"
        );
        assert!(
            matches!(
                collapsed.report(&read_groups).contamination(),
                ContaminationUsed::PerReadGroup(_)
            ),
            "and the run reports it as one"
        );

        // **And the same table with no row measured is refused**, rather than becoming a run:
        // that is the uncontaminated state written longhand, and it is the file a reader who
        // collapsed the two would produce.
        let mut every_row_unmeasured = longhand;
        for row in &mut every_row_unmeasured
            .contamination
            .as_mut()
            .expect("a table")
            .by_read_group
        {
            row.measurement = None;
        }
        let error = every_row_unmeasured
            .to_run_parameters()
            .expect_err("an uncontaminated run written longhand is not a run");
        assert!(
            format!("{error}").contains("leave the section out instead"),
            "{error}"
        );
    }

    /// **Row 2. Measured and found clean is not the same as never measured**, and only the
    /// evidence counts tell them apart.
    ///
    /// **The answer that changes is what the run reports for that library.** Both are a fraction
    /// of zero, so both correct nothing; what differs is whether the run may say the library was
    /// looked at. A cohort where one lane was measured and found clean and another could not be
    /// measured at all is a cohort where one number is a result and the other is a gap, and an
    /// output that printed both as "0" would be reporting a gap as a finding.
    #[test]
    fn a_zero_fraction_with_evidence_is_not_an_unmeasured_one() {
        let read_groups = the_files_read_groups(&a_file_using_every_shape());
        // **The read group, found by its own row rather than by a position.** The fixture's
        // rows happen to be in read-group order, so an index would agree today and would stop
        // agreeing the first time a row moved — which has happened once already in this file.
        let clean_at = a_file_using_every_shape()
            .contamination
            .expect("the fixture has a table")
            .by_read_group
            .iter()
            .find(|row| {
                row.measurement
                    .as_ref()
                    .is_some_and(|found| found.share_of_reads_from_another_sample == 0.0)
            })
            .expect("one lane was measured and found clean")
            .read_group as usize;

        let file = a_file_using_every_shape();
        let measured = a_run_from(&file);
        let view = measured.view().contamination_by_read_group()[clean_at];
        assert_eq!(view.fraction, 0.0);
        assert!(
            view.was_measured(),
            "a zero fraction with evidence behind it is a measurement"
        );

        let mut collapsed = a_file_using_every_shape();
        collapsed
            .contamination
            .as_mut()
            .expect("a table")
            .by_read_group[clean_at]
            .measurement = None;
        let unmeasured = a_run_from(&collapsed);
        let view = unmeasured.view().contamination_by_read_group()[clean_at];
        assert_eq!(
            view.fraction, 0.0,
            "the fraction is the same number in both, which is the whole difficulty"
        );
        assert!(!view.was_measured());

        // **What the run says about that library differs**, and it is the only thing that does.
        let of_the_clean_lane = |run: &RunParameters| {
            let report = run.report(&read_groups);
            let ContaminationUsed::PerReadGroup(rows) = report.contamination() else {
                panic!("this run fitted contamination somewhere")
            };
            rows.iter()
                .find(|row| row.read_group.get() as usize == clean_at)
                .expect("the lane has a row")
                .was_measured()
        };
        assert_eq!(
            (of_the_clean_lane(&measured), of_the_clean_lane(&unmeasured)),
            (true, false)
        );
    }

    /// **Row 3. A defaulted multiplier of 1.0 is not a fitted one**, and nothing a run computes
    /// can tell them apart.
    ///
    /// **That is why this row is in spec §5 at all.** The two multiply every read's reported error
    /// probability by exactly the same number, so no score, no genotype and no quality differs —
    /// the *only* thing that differs is the warrant, which is why §5's rule is that a warrant is
    /// carried and never inferred from the value. **The answer that changes is what the run
    /// writes down**: a run that read a defaulted scale and wrote `fitted_here` would tell its
    /// next reader that a library nothing could be fitted for was calibrated against a
    /// measurement.
    #[test]
    fn a_defaulted_multiplier_of_one_is_not_a_fitted_one() {
        // Likewise found by its own row: the one multiplier the fixture marks `defaulted`.
        let defaulted_at = a_file_using_every_shape()
            .base_quality_calibration
            .by_read_group
            .iter()
            .find(|row| row.error_probability_multiplier.warrant == Warrant::Defaulted)
            .expect("one library's rate could not be fitted")
            .read_group as usize;

        let file = a_file_using_every_shape();
        let run = a_run_from(&file);
        let calibration = run.view().calibration_by_read_group()[defaulted_at];
        assert_eq!(calibration.scale, 1.0);
        assert_eq!(calibration.provenance, Provenance::Defaulted);

        let mut collapsed = a_file_using_every_shape();
        collapsed.base_quality_calibration.by_read_group[defaulted_at]
            .error_probability_multiplier
            .warrant = Warrant::FittedHere;
        let collapsed_run = a_run_from(&collapsed);
        let collapsed_calibration = collapsed_run.view().calibration_by_read_group()[defaulted_at];

        assert_eq!(
            collapsed_calibration.scale, calibration.scale,
            "the number is the same in both — a scale of one is a legitimate fitted answer as \
             well as the default, which is the whole of why the warrant travels with it"
        );
        assert_ne!(collapsed_calibration.provenance, calibration.provenance);

        // **And a score that rests on this library says so, or does not.** Spec §2's rule is
        // that consumers *combine* warrants rather than branching on them: a call resting on one
        // fitted parameter and one defaulted one is a defaulted call. So the answer that changes
        // is the warrant every SNP or indel call touching this library reports — which
        // `summarise_condition`'s fold reads off exactly this field, and which is the only place
        // in calling that reads it at all.
        let as_a_call_would_report_it = |calibration: ReadGroupCalibration| {
            Provenance::FittedHere.weaker_of(calibration.provenance)
        };
        assert_eq!(
            (
                as_a_call_would_report_it(calibration),
                as_a_call_would_report_it(collapsed_calibration)
            ),
            (Provenance::Defaulted, Provenance::FittedHere),
            "a call whose every other parameter was fitted is a defaulted call where this \
             library's calibration was defaulted, and a fitted one where it was not"
        );
    }

    /// **Row 4. A stratum with no length spectrum of its own was furnished from its period's
    /// curves**, and giving it one changes the prior every tract of it is seeded from.
    ///
    /// **This row changes a number.** The length spectrum is the shape a repeat tract's genotype
    /// prior starts from, and the concentration is how many chromosomes' worth of belief it is
    /// held with — so a reader that invented a spectrum for a stratum that has none would seed
    /// every tract of that class from a distribution the fit never produced, at a strength it
    /// never chose.
    #[test]
    fn a_stratum_with_no_spectrum_of_its_own_falls_to_its_periods_and_the_numbers_differ() {
        let file = a_file_using_every_shape();
        // The fixture fits period 2 at 6 repeats and pools period 2; period 2 at 11 repeats has
        // no spectrum of its own and falls to that pool.
        let (period, furnished, fitted) = (2_u8, 11_u64, 6_u64);

        let run = a_run_from(&file);
        let fell_back = run
            .ssr_slippage_fits()
            .length_spectrum_at(period, furnished);
        assert_eq!(fell_back.rung(), LengthSpectrumRung::PeriodsPooledTracts);

        let mut collapsed = a_file_using_every_shape();
        let its_own = StratumLengthSpectrumRow {
            period,
            reference_repeats: furnished,
            concentration: 9.5,
            shares_by_repeat_offset: vec![0.25, 0.5, 0.25],
        };
        collapsed
            .repeat_tracts
            .length_spectrum_by_stratum
            .push(its_own);
        collapsed
            .repeat_tracts
            .length_spectrum_by_stratum
            .sort_by_key(|row| (row.period, row.reference_repeats));
        // **The fallback concentration is left where it was**, and with two fitted strata beside
        // it the median it claims is no longer theirs. Nothing here reads the bottom rung, and
        // `validate` checks that key's warrant rather than its value — but a file with an
        // invented spectrum is a file no run wrote, and this is the one number in it that says
        // so.
        let collapsed_run = a_run_from(&collapsed);
        let invented = collapsed_run
            .ssr_slippage_fits()
            .length_spectrum_at(period, furnished);

        assert_eq!(invented.rung(), LengthSpectrumRung::StratumsOwnFit);
        assert_ne!(
            invented.concentration(),
            fell_back.concentration(),
            "the strength the prior is held with is a different number: 9.5 chromosomes against \
             the {} its period pooled",
            fell_back.concentration()
        );
        assert_ne!(invented.fitted_weights(), fell_back.fitted_weights());
        // **And the stratum that was fitted is untouched**, so this is a difference at the
        // stratum the row is about rather than a change to the whole period.
        assert_eq!(
            collapsed_run
                .ssr_slippage_fits()
                .length_spectrum_at(period, fitted)
                .concentration(),
            run.ssr_slippage_fits()
                .length_spectrum_at(period, fitted)
                .concentration()
        );
    }

    /// **Row 5. A `(stratum × slippage group)` with no row put no read there**, and a slip rate of
    /// zero says something a fit never found.
    ///
    /// **This row changes a number too, and the number is the one a read is scored against.**
    /// `TractScoringFits::gather_for_locus` takes the two answers to different places: a lookup
    /// that fails takes `StutterModel::hipstr_shipped`, and one that succeeds takes
    /// `stutter_model_for` on the numbers it found. So the share of reads the model expects to
    /// come back one repeat short is **0.05 under the shipped model and exactly 0 under a slip
    /// rate of zero**. At zero the emission for such a read collapses to exactly nothing, and the
    /// only term left to explain it is the row's outlier weight — one read in a hundred by
    /// default. **The genotype is not ruled out**: it is charged the junk rate for every read
    /// that slipped, which at a stratum where reads do slip is every read of one allele.
    #[test]
    fn a_pair_with_no_row_is_not_a_slip_rate_of_zero() {
        let file = a_file_using_every_shape();
        // The fixture's read group 2 is the only one in slippage group 1, and period 1 at 30
        // repeats has a row for group 0 alone.
        let (silent, period, repeats) = (ReadGroupId(2), 1_u8, 30_u64);

        let run = a_run_from(&file);
        assert_eq!(
            run.ssr_slippage_fits().at(silent, period, repeats),
            Err(NoSlippage::GroupPutNoReadHere { slippage_group: 1 }),
            "the lookup names the group it looked under, and says it put no read here"
        );

        // **What the caller does with each answer**, which is where the number a read is scored
        // against comes from: a failed lookup takes the shipped model, a successful one takes the
        // numbers it found (`TractScoringFits::gather_for_locus`).
        let shipped = StutterModel::hipstr_shipped();
        assert_eq!(
            shipped.whole_repeat_shorter_share(),
            0.05,
            "the model a cell with no numbers falls to expects one read in twenty to come back a \
             repeat short"
        );
        let never_slips = stutter_model_for(&Slippage {
            level: 0.0,
            shorter_share: 0.8,
            fall_off: 0.3,
        });
        assert_eq!(
            never_slips.whole_repeat_shorter_share(),
            0.0,
            "a slip rate of zero expects none: every direction share is that rate times a share \
             of it, so a zero rate is a model that says this stratum's reads never slip"
        );

        // And what the stratum's own fitted numbers give, which is the answer this file carries
        // for the group that *did* put reads there.
        let fitted = run
            .ssr_slippage_fits()
            .at(ReadGroupId(0), period, repeats)
            .expect("slippage group 0 has numbers at this stratum");
        assert!(
            stutter_model_for(&fitted.slippage).whole_repeat_shorter_share() > 0.0,
            "a fitted stratum expects a real share of its reads to come back short"
        );

        // **The collapse is a row of zeros**, and the lookup then answers where it answered
        // nothing — which is the whole difference.
        let mut collapsed = a_file_using_every_shape();
        let mut zeroed = collapsed.repeat_tracts.slippage_by_stratum_and_group
            [super::super::tests::THE_ROW_WHOSE_SHARES_BLEND]
            .clone();
        zeroed.slippage_group = 1;
        zeroed.share_of_reads_that_slip = 0.0;
        // **And no slipped-read count beside it**, so the invented row models only the collapse
        // this test is about rather than also claiming twelve thousand slipped reads at a
        // stratum where none slips.
        zeroed
            .share_of_reads_that_slip_origin
            .expected_slipped_reads = None;
        if let Some(shares) = &mut zeroed.shorter_share_and_fall_off_origin {
            shares.expected_slipped_reads = None;
        }
        collapsed
            .repeat_tracts
            .slippage_by_stratum_and_group
            .push(zeroed);
        let answered = a_run_from(&collapsed)
            .ssr_slippage_fits()
            .at(silent, period, repeats)
            .expect("the invented row answers");
        assert_eq!(answered.slippage.level, 0.0);
        assert_eq!(
            stutter_model_for(&answered.slippage).whole_repeat_shorter_share(),
            0.0,
            "and the model that library's reads are scored against at this stratum then expects \
             none of them to slip, where the absent row would have taken the shipped 0.05"
        );
    }
}
