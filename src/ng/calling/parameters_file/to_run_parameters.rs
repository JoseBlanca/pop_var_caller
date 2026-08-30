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
    fn the_files_read_groups(file: &ParametersFile) -> ReadGroups {
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
    fn the_rates_the_projection_out_reads(
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
