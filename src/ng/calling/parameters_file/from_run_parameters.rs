//! **What a run used, projected onto the file's shape** — every number
//! [`RunParameters`](crate::ng::calling::run_parameters::RunParameters) holds, mapped to the row
//! that carries it (`doc/devel/ng/spec/parameters_file.md` §3).
//!
//! No TOML here: this step produces a [`ParametersFile`] value and stops. The text is step B2's,
//! and keeping the two apart is what lets the projection be tested on the shape rather than on a
//! rendering of it.
//!
//! What the projection needs beyond the assembled parameters, and the rules it keeps, are
//! documented on [`ParametersFile::of_run`] itself — this module is private, so prose here
//! reaches no rendered page.

use std::collections::BTreeMap;

use super::bindings::hex_digest;
use super::{
    BaseQualityCalibration, BaseQualityCalibrationRow, CensusIdentity, Contamination,
    ContaminationFittedFrom, ContaminationMeasurement, ContaminationRow, CurveReach, EvidenceCount,
    FORMAT_VERSION, Inbreeding, InbreedingRow, InputsFittedFrom, LevelOrigin, LevelSmoothing,
    OrdinarySitePrior, ParametersFile, PeriodLengthSpectrumRow, ReadGroupBatchRow, ReadGroupRow,
    RepeatTracts, SampleBatchRow, SeedRung, SequencingBatches, ShareCurve, ShareCurveRung,
    ShareShape, ShareSmoothing, SharesOrigin, SlippageCurve, SlippageGroupRow, SlippageRow,
    StatedConstants, StratumLengthSpectrumRow, SubstitutionRateRow, Warrant, WarrantedValue,
};
use crate::ng::calling::genotype_prior::SeedRegime;
use crate::ng::calling::likelihood::{ContaminationView, ReadGroupCalibration};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::parameter_estimation::joint::census::Stratum;
use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches as DeclaredBatches;
use crate::ng::parameter_estimation::joint::share_curve::{
    ShareCurve as FittedShareCurve, ShareCurveSource, ShareShape as FittedShape, ShareSource,
};
use crate::ng::parameter_estimation::joint::slippage_curve::{
    CurveReach as FittedCurveReach, LevelSource, SlippageCurve as FittedSlippageCurve,
};
use crate::ng::parameter_estimation::joint::ssr_fit::{LevelProvenance, ShareProvenance};
use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use crate::ng::parameter_estimation::{Estimate, Provenance};
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::types::{ErrorRate, InbreedingF, ReadGroupId};

impl ParametersFile {
    /// **Write down what this run scored its reads under.**
    ///
    /// `run` is the parameters calling read; `read_groups` is the run's own table, which supplies
    /// every name in the file; `base_quality_rate_by_read_group` and `inbreeding_by_sample` are
    /// the two pre-pass estimate sets assembly consumed and did not keep (below);
    /// `reference` and `census` are what the file is bound to (spec §3.1, §6), and neither is
    /// derivable from the numbers.
    ///
    /// # It reads more than `RunParameters`, and that is the point
    ///
    /// **Assembly drops two things the file needs**, and a projection written from the assembled
    /// parameters alone would produce a plausible file that is wrong about both:
    ///
    /// - **the base-quality calibration's evidence count.**
    ///   [`ReadGroupCalibration`](crate::ng::calling::likelihood::ReadGroupCalibration) is a
    ///   multiplier and a warrant with no count on it, and spec §3.3 asks for the count by name.
    ///   It is on the `Estimate<ErrorRate>` that assembly reads and does not store.
    /// - **every sample's inbreeding warrant and count.** The pre-pass fits an
    ///   `Estimate<InbreedingF>` per sample and the seam into calling takes the bare
    ///   coefficients, so a file written from those would mark every sample's number as handed
    ///   over when the run fitted it.
    ///
    /// # Three rules the whole projection keeps
    ///
    /// **A `defaulted` value carries no evidence count.** A stated constant has nothing behind
    /// it, and writing the count of reads the fit *did* see beside a multiplier of one would say
    /// that one rests on them. The rule is applied to every warranted number rather than to the
    /// calibration alone — a defaulted substitution rate carries no count either.
    ///
    /// **A `supplied` value keeps its count, and the asymmetry is deliberate.** Spec §2.1 demotes
    /// every number of a mismatched file to `supplied` wholesale, and those numbers were fitted
    /// on *some* cohort — the evidence behind them is real and stays worth reporting. A defaulted
    /// number has none by construction.
    ///
    /// **Absence is written as absence.** A read group that identified no contamination gets a
    /// row with no `measurement`, never a fraction of zero beside two zero counts; a `(stratum ×
    /// slippage group)` with no reads gets no row; a stratum with no length spectrum of its own
    /// gets no row. Those are three of spec §5's five states, and this is the only code that can
    /// get them wrong on the way out.
    ///
    /// # The two derived bindings, and why only one of them is this run's to choose
    ///
    /// **The reference arrives as the run's own [`ReferenceDigest`] and is spelled here**, so
    /// that the string this writes and the string a later run compares it against come out of
    /// one function — step D2's is what does the comparing. A file naming another reference is
    /// refused (spec §6), so a file this run is allowed to write already names this run's
    /// reference, and the run's own digest is the only value this could write.
    ///
    /// **The census arrives already minted**, by [`CensusIdentity::of`], and there are two
    /// reasons rather than one. A mismatched census *demotes* rather than refusing (spec §6), so
    /// a run that read a file fitted under other terms and writes its parameters out again (spec
    /// §7 — every run writes the file it used) has to write back the terms it read, not its own.
    /// And a **direct-mode run has no census at all** — no pre-pass and no psp
    /// (`run_streaming.md` §2) — so there is no `RecordingTerms` for it to hand over. Which
    /// census a file names is therefore its caller's to say, where which reference it names is
    /// not.
    ///
    /// **⚑ One run cannot satisfy this signature**: `ReferenceDigest::of` refuses a reference
    /// read from a `.fai` alone, which holds no bases to digest, while spec §7 says writing is
    /// unconditional. What such a run writes is step F1's, at the call site.
    ///
    /// **Six arguments rather than a bundle type**, on the same grounds
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)
    /// gives for its nine: the list is the point, a struct naming the same six things would be a
    /// second place for them to go out of step with the run, and no two of the six share a type,
    /// so no pair can be exchanged at a call site.
    ///
    /// # Panics
    ///
    /// Held in release. Each is two tables minted from different inputs and joined positionally,
    /// whose symptom is one library's numbers written under another's name — which looks like an
    /// answer rather than a failure:
    ///
    /// - **the read-group table covers this run's libraries**, and **covers this run's samples**.
    ///   Two checks, whose messages share no opening phrase: a test pinning one of them by a
    ///   shared substring passes when the other fires, which is how the first draft of this
    ///   step's own test could not fail.
    /// - **the fitted rates name exactly the run's read groups** — the same count, none missing,
    ///   and each rate's warrant is the one the run's calibration carries (or that calibration is
    ///   `Defaulted`, which is what assembly substitutes when it refuses a rate). The count a row
    ///   writes comes from that rate, so a rate set from another fit would write another run's
    ///   reads beside this run's multipliers.
    /// - **the inbreeding estimates are the run's own**, one per sample and each carrying the
    ///   coefficient the parameters hold.
    /// - **every repeat-tract substitution rate is keyed to one of the run's read groups.** Every
    ///   section of the file joins on the run's dense index, so a rate naming a read group the
    ///   identity block does not list is a rate no reader can attach to a library.
    /// - **a read group that identified no contamination carries a fraction of zero.** Such a row
    ///   is written with no measurement, so a non-zero fraction there is one the run scored with
    ///   and the file cannot carry.
    /// - **a slippage number that came off a curve carries the curve, and says whether the
    ///   stratum sat inside the curve's fitted range.** The file cannot say *the curve supplied
    ///   this* without saying which curve and how far it reached, by design (spec §5's rule that
    ///   a warrant is carried and never inferred).
    /// - **a slippage group with numbers at a stratum has a level provenance beside them**, which
    ///   `StratumFits::each_stratum_and_group_with_numbers` refuses naming the stratum.
    ///
    /// **Panics rather than a `Result`, and this is the open half of that choice.** These are
    /// wiring bugs in whoever assembles the six arguments, not states any input data can reach,
    /// and the alternative to refusing is writing a false provenance. What a panic costs here is
    /// larger than at assembly, though: this runs *after* the last locus, so it discards a
    /// cohort's calling work where
    /// [`RunParameters::assemble`](crate::ng::calling::run_parameters::RunParameters::assemble)'s
    /// equivalent checks discard a startup. **Nothing calls this yet** — whether the run driver
    /// should instead log and keep its VCF is step F1's to settle, at the call site, where the
    /// order of the two writes is decided.
    #[must_use]
    pub fn of_run(
        run: &RunParameters,
        read_groups: &ReadGroups,
        base_quality_rate_by_read_group: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
        inbreeding_by_sample: &[Estimate<InbreedingF>],
        reference: &ReferenceDigest,
        census: CensusIdentity,
    ) -> Self {
        let read_group_count = run.read_group_count();
        let samples = read_groups.read_groups_per_sample();
        assert_eq!(
            read_groups.len(),
            read_group_count,
            "this run's read-group table covers {} libraries and its parameters cover {}; the two \
             were minted from different inputs, and joining them positionally would write one \
             library's calibration under another's name",
            read_groups.len(),
            read_group_count
        );
        assert_eq!(
            samples.len(),
            run.inbreeding_coefficient_by_sample().len(),
            "the read-group table names {} samples and the run's parameters cover {}; every \
             per-sample row of the file is written in the run's sample order, which is this \
             table's own first-seen order",
            samples.len(),
            run.inbreeding_coefficient_by_sample().len()
        );
        // **A rate for every read group, or none at all** — the same two-state contract the
        // contamination axis carries one level down, and for the same reason: a *short* list is
        // the failure worth refusing, because it writes some other read group's count beside a
        // multiplier. **Empty is the run with no fit** (spec §8), whose calibrations are all
        // `Defaulted` and so write no count; `calibration_rows` refuses a missing rate under any
        // other warrant, read group by read group.
        assert!(
            base_quality_rate_by_read_group.is_empty()
                || base_quality_rate_by_read_group.len() == read_group_count,
            "the fit supplied base-quality rates for {} read groups and the run has {}; the count \
             written beside each multiplier comes from these rates, so a set covering some of \
             them is a different fit — a run that fitted none supplies none",
            base_quality_rate_by_read_group.len(),
            read_group_count
        );

        Self {
            format_version: FORMAT_VERSION,
            ploidy: run.ploidy().get(),
            fitted_from: InputsFittedFrom {
                reference_digest: hex_digest(&reference.0),
                samples: samples
                    .iter()
                    .map(|of_sample| of_sample.sample.to_string())
                    .collect(),
                read_groups: read_groups
                    .iter()
                    .map(|(id, declared)| ReadGroupRow {
                        read_group: id.get(),
                        declared_id: declared.id.to_string(),
                        library: declared.library.value.to_string(),
                        sample: declared.sample.to_string(),
                    })
                    .collect(),
                census,
            },
            base_quality_calibration: BaseQualityCalibration {
                by_read_group: calibration_rows(
                    run.calibration_by_read_group(),
                    base_quality_rate_by_read_group,
                ),
            },
            contamination: contamination_of(run.contamination_by_read_group(), read_groups),
            sequencing_batches: sequencing_batches_of(run.sequencing_batches(), read_groups),
            inbreeding: inbreeding_of(
                run.inbreeding_coefficient_by_sample(),
                read_groups,
                inbreeding_by_sample,
            ),
            ordinary_site_prior: OrdinarySitePrior {
                reference_concentration: run.prior_seed().alpha_ref(),
                alternative_concentration_total: run.prior_seed().alpha_alt_total(),
                rung: run.prior_seed().regime().into(),
            },
            repeat_tracts: repeat_tracts_of(run),
            stated_constants: StatedConstants {
                // **Whatever this run scored under, and its warrant with it.** Until 2026-08-30
                // `RunParameters` had no field for this number and the projection wrote the
                // compiled-in constant marked `defaulted`, so a file whose weight a person had
                // edited round-tripped back to the default. The weight now rides on the run
                // (`RunParameters::repeat_tract_outlier_weight`), and its two reachable
                // warrants are `supplied` — read out of a parameters file — and `defaulted`.
                //
                // **No evidence count either way**, which is the projection's rule for a
                // defaulted value and true here of a supplied one too: nothing counted anything
                // to arrive at it.
                repeat_tract_outlier_weight: WarrantedValue {
                    value: run.repeat_tract_outlier_weight().value(),
                    warrant: run.repeat_tract_outlier_weight().provenance().into(),
                    observations: None,
                },
            },
        }
    }
}

/// One row a read group: the multiplier and its warrant from the run, the count from the rate the
/// multiplier was built from.
fn calibration_rows(
    calibration_by_read_group: &[ReadGroupCalibration],
    rate_by_read_group: &BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
) -> Vec<BaseQualityCalibrationRow> {
    calibration_by_read_group
        .iter()
        .enumerate()
        .map(|(read_group, calibration)| {
            let id = ReadGroupId(
                u32::try_from(read_group).expect("a run's read-group axis fits in a u32"),
            );
            // **A read group with no rate at all is legal exactly where its calibration is
            // `Defaulted`, and nowhere else.** That is the run with no fit (spec §8): nothing was
            // fitted for anybody, so there are no rates to offer, and a `Defaulted` calibration
            // writes no `observations` anyway — the projection's rule that a stated constant has
            // nothing behind it. So the row this builds is the same row whether a rate was
            // offered or not, and requiring one would mean a caller inventing an `Estimate` to
            // satisfy a lookup whose value is then dropped.
            //
            // **The panic stays for every other warrant**, and it is the one worth keeping: a
            // *fitted* calibration whose rate is missing means the rate set and the calibration
            // axis came from different fits, and the count that row is about to write is then
            // the other fit's.
            let rate = rate_by_read_group.get(&id);
            assert!(
                rate.is_some() || calibration.provenance == Provenance::Defaulted,
                "read group {read_group}'s calibration is {:?} and no rate was offered for it; \
                 only a `Defaulted` calibration can have none, because it is the one that writes \
                 no count — so a rate set missing a fitted read group's entry was fitted over a \
                 different set of read groups than the run has",
                calibration.provenance
            );
            // **The scale carries the rate's own warrant** — `from_fitted_rate` copies it — with
            // one legitimate disagreement: a rate assembly *refused* leaves the calibration
            // `Defaulted` beside a rate that is not. Anything else means the two came from
            // different fits, and the count this row is about to write is then the other fit's.
            assert!(
                calibration.provenance == Provenance::Defaulted
                    || rate.is_none_or(|rate| calibration.provenance == rate.provenance),
                "read group {read_group}'s calibration is {:?} and the rate offered for it is \
                 {:?}; the scale carries the rate's own warrant and the count beside it comes \
                 from that rate, so two that disagree are two different fits",
                calibration.provenance,
                rate.map(|rate| rate.provenance)
            );
            BaseQualityCalibrationRow {
                read_group: id.get(),
                error_probability_multiplier: warranted_value(
                    calibration.scale,
                    calibration.provenance,
                    // **The count on the no-rate arm is never read, and the two lines that make
                    // that true are both above.** The assertion admits a missing rate only under
                    // `Defaulted`, and `warranted_value` writes no `observations` at all under
                    // `Defaulted` — so any number here produces the same row. Measured: putting
                    // 7 in its place passes all 184 tests of `ng::calling::parameters_file`,
                    // this module's 33 among them, which is an equivalent mutant rather than a
                    // hole. Zero is written because it is the true count.
                    EvidenceCount::Reads(rate.map_or(0, |rate| rate.observations)),
                ),
            }
        })
        .collect()
}

/// The contamination section, **or its absence**.
///
/// **Absent is the uncontaminated run** — the state assembly records by leaving its list empty,
/// and the one a table of zeros would misreport as every library measured and found clean
/// (spec §5, first row).
fn contamination_of(
    contamination_by_read_group: &[ContaminationView],
    read_groups: &ReadGroups,
) -> Option<Contamination> {
    if contamination_by_read_group.is_empty() {
        return None;
    }
    Some(Contamination {
        by_read_group: read_groups
            .iter()
            .zip(contamination_by_read_group)
            .map(|((id, declared), view)| {
                // **What an unmeasured row drops is a fraction that must be zero.** Such a row
                // carries no measurement, so a non-zero fraction here would be one calling scored
                // with and the file cannot express — psp mode and direct mode would then score
                // the same reads differently, which is what the two-mode oracle exists to catch
                // and would find only by diffing two VCFs. The estimator holds this today
                // (`fit_alpha` returns exactly zero where no marker carries a read); nothing but
                // this states it.
                assert!(
                    view.was_measured() || view.fraction == 0.0,
                    "read group {} identified no contamination and carries a fraction of {}; the \
                     file writes no measurement for such a row, so a non-zero fraction there is \
                     one the run scored with and the file cannot carry",
                    id.get(),
                    view.fraction
                );
                ContaminationRow {
                    read_group: id.get(),
                    library: declared.library.value.to_string(),
                    // **`was_measured` and not the fraction.** A read group that identified
                    // nothing carries a fraction of zero beside two zero counts in memory, and
                    // reading the fraction would write that out as a measurement of zero — the
                    // one thing spec §5's second row says a reader must never do.
                    measurement: view.was_measured().then(|| ContaminationMeasurement {
                        share_of_reads_from_another_sample: view.fraction,
                        markers_with_reads: view.markers_with_reads,
                        reads_on_markers: view.reads_on_markers,
                        fitted_from_reads_of: view.source.into(),
                    }),
                }
            })
            .collect(),
    })
}

/// Who was sequenced beside whom, with each row naming its own axis value.
fn sequencing_batches_of(batches: &DeclaredBatches, read_groups: &ReadGroups) -> SequencingBatches {
    SequencingBatches {
        // **The only thing that tells a declared batching from an assumed one.** The dense rows a
        // run that declared nothing writes and the rows a run that declared one batch holding
        // every library writes are the same rows.
        batching_was_declared: !batches.is_default(),
        by_read_group: read_groups
            .iter()
            .zip(batches.of_each_read_group().0)
            .map(|((id, _), batch)| ReadGroupBatchRow {
                read_group: id.get(),
                batch: batch.get(),
            })
            .collect(),
        by_sample: read_groups
            .read_groups_per_sample()
            .iter()
            .zip(batches.of_each_sample().0)
            .map(|(of_sample, batch)| SampleBatchRow {
                sample: of_sample.sample.to_string(),
                batch: batch.get(),
            })
            .collect(),
    }
}

/// One row a sample: the name from the run's table, the value, warrant and count from the fit.
fn inbreeding_of(
    assembled_coefficients: &[InbreedingF],
    read_groups: &ReadGroups,
    estimates: &[Estimate<InbreedingF>],
) -> Inbreeding {
    assert_eq!(
        estimates.len(),
        assembled_coefficients.len(),
        "the fit supplied {} inbreeding estimates and the run's parameters hold {} coefficients; \
         the warrants written beside these values come from the estimates, so two lists of \
         different lengths are two different fits",
        estimates.len(),
        assembled_coefficients.len()
    );
    Inbreeding {
        by_sample: read_groups
            .read_groups_per_sample()
            .iter()
            .zip(estimates)
            .zip(assembled_coefficients)
            .map(|((of_sample, estimate), coefficient)| {
                // **Exact equality is the intended comparison.** The seam copies the estimate's
                // value into assembly rather than recomputing it, and `InbreedingF` refuses a
                // `NaN`, so two numbers that are not bit-identical are two different fits.
                assert_eq!(
                    estimate.value.get(),
                    coefficient.get(),
                    "sample {} was assembled with an inbreeding coefficient of {} and the \
                     estimate offered here holds {}; the warrant and the count would then \
                     describe a number the run never used",
                    of_sample.sample,
                    coefficient.get(),
                    estimate.value.get()
                );
                InbreedingRow {
                    sample: of_sample.sample.to_string(),
                    inbreeding_coefficient: warranted_value(
                        estimate.value.get(),
                        estimate.provenance,
                        EvidenceCount::CoveredPositions(estimate.observations),
                    ),
                }
            })
            .collect(),
    }
}

/// The whole repeat-tract section: the slippage numbers, both fitted rungs of the length-spectrum
/// ladder, and the substitution rate.
fn repeat_tracts_of(run: &RunParameters) -> RepeatTracts {
    let fits: &StratumFits = run.ssr_slippage_fits();
    let read_group_count = run.read_group_count();
    RepeatTracts {
        // **The run's own warrant, copied rather than re-derived** — fitted where it had strata
        // to take a median over, defaulted where it had none, and *supplied* where a parameters
        // file handed the number over, which is the state no arithmetic over this run's strata
        // could reach. Until 2026-08-30 this worked the warrant out from
        // `strata_with_a_length_spectrum()`, and a file demoted under spec §2.1 came back out
        // saying `fitted_here`.
        //
        // **The value and the warrant can disagree about how interesting they are**: a median
        // over fitted strata can land on exactly the stated 1.0, which is why the warrant
        // travels at all.
        fallback_length_spectrum_concentration: WarrantedValue {
            value: fits.stated_concentration(),
            warrant: fits.stated_concentration_warrant().into(),
            observations: None,
        },
        // **One row a read group the run declared a slippage group for.** A read group the
        // declaration does not name has no row, which is what the fit itself says about it —
        // `StratumFits::at` answers `UnknownReadGroup` and no slippage number is ever looked up
        // under it.
        slippage_group_by_read_group: (0..read_group_count)
            .filter_map(|read_group| {
                let id = ReadGroupId(
                    u32::try_from(read_group).expect("a run's read-group axis fits in a u32"),
                );
                fits.slippage_group_of(id)
                    .map(|slippage_group| SlippageGroupRow {
                        read_group: id.get(),
                        slippage_group,
                    })
            })
            .collect(),
        slippage_by_stratum_and_group: fits
            .each_stratum_and_group_with_numbers()
            .map(|(stratum, slippage_group, fitted)| SlippageRow {
                period: stratum.period,
                reference_repeats: stratum.reference_repeats,
                slippage_group,
                share_of_reads_that_slip: fitted.slippage.level,
                shorter_share: fitted.slippage.shorter_share,
                fall_off: fitted.slippage.fall_off,
                share_of_reads_that_slip_origin: level_origin_of(
                    fitted.level,
                    stratum,
                    slippage_group,
                ),
                shorter_share_and_fall_off_origin: fitted.shares.map(|shares| SharesOrigin {
                    expected_slipped_reads: shares.slipped_reads,
                    shorter_share_smoothing: share_smoothing_of(
                        shares.shorter_share,
                        stratum,
                        slippage_group,
                        "the share of slipped reads showing a shorter tract",
                    ),
                    fall_off_smoothing: share_smoothing_of(
                        shares.fall_off,
                        stratum,
                        slippage_group,
                        "the fall-off",
                    ),
                }),
            })
            .collect(),
        length_spectrum_by_stratum: fits
            .fitted_length_spectrum_of_each_stratum()
            .map(
                |(stratum, shares_by_repeat_offset, concentration)| StratumLengthSpectrumRow {
                    period: stratum.period,
                    reference_repeats: stratum.reference_repeats,
                    concentration,
                    shares_by_repeat_offset: shares_by_repeat_offset.to_vec(),
                },
            )
            .collect(),
        length_spectrum_by_period: fits
            .pooled_length_spectrum_of_each_period()
            .map(
                |(period, shares_by_repeat_offset, concentration)| PeriodLengthSpectrumRow {
                    period,
                    concentration,
                    shares_by_repeat_offset: shares_by_repeat_offset.to_vec(),
                },
            )
            .collect(),
        substitution_rate_by_stratum: run
            .ssr_substitution_rate()
            .map(|(key, rate)| {
                assert!(
                    (key.read_group.get() as usize) < read_group_count,
                    "a repeat-tract substitution rate is keyed to read group {} and the run has \
                     {read_group_count} read groups; every section of the file joins on the run's \
                     dense index, so a row naming one the identity block does not list is a rate \
                     no reader can attach to a library",
                    key.read_group.get()
                );
                SubstitutionRateRow {
                    read_group: key.read_group.get(),
                    period: key.stratum.period.get(),
                    // **One width for a repeat count**, as this module's own conventions say: the
                    // fit's STR key counts repeats in a `u32` and the census's stratum in a
                    // `u64`, and they are the same quantity.
                    reference_repeats: u64::from(key.stratum.repeats.0),
                    // **The key's ploidy and not the run's.** It is the set of genotypes the fit
                    // scored this table's entries against, which is part of what the rate means.
                    ploidy: key.ploidy.get(),
                    rate: warranted_value(
                        rate.value.get(),
                        rate.provenance,
                        EvidenceCount::BasesCompared(rate.observations),
                    ),
                }
            })
            .collect(),
    }
}

/// A number, its warrant, and its count **in the unit that number counts in** — with the count
/// dropped where the warrant says the value is a stated constant.
///
/// **`Estimate<T>`'s count is a bare integer and the file's names a unit**, because the three
/// units in this file differ by orders of magnitude on one cohort and a reader comparing two
/// counts without knowing which is which is not comparing anything. Which unit a quantity counts
/// in is fixed at the call site, where the quantity is known.
///
/// **A count of zero is written rather than dropped.** *Fitted from nothing* and *a stated
/// constant* are spec §5's third row read in both directions, and only the warrant separates them.
fn warranted_value(
    value: f64,
    provenance: Provenance,
    observations: EvidenceCount,
) -> WarrantedValue {
    let warrant = Warrant::from(provenance);
    WarrantedValue {
        value,
        warrant,
        // **A defaulted value has no evidence behind it, whatever the fit saw.** A read group
        // whose fitted rate was refused still has a count of reads on that rate; writing it here
        // would say the multiplier of one rests on them, and it rests on nothing. `Supplied` is
        // deliberately *not* treated this way — see `of_run`'s three rules.
        //
        // **A `supplied` number with nothing behind it writes no count either**, which is a
        // narrower rule than the one above and was added at step E2. `Supplied` keeps its count
        // in general, and must: a file demoted under spec §2.1 marks every number `supplied`,
        // and those numbers were fitted on some cohort. But a coefficient an operator typed
        // (`parameters_file::DeclaredInbreeding`) is `Supplied` with a count of **zero**, and
        // writing `observations = { covered_positions = 0 }` then says the number was measured
        // over no genome at all — while the file's own editing rule, three lines from the top,
        // tells the reader to *delete* `observations` on a value they supplied. A geneticist
        // reading a produced file could not tell whether that row was correct or their own
        // mistake.
        //
        // **A *fitted* count of zero is still written**, and deliberately: *this fit produced a
        // number from no reads at all* is an alarming state and the count is what says so.
        // **The trip stays lossless** — the projection back reads an absent count as zero
        // (`to_run_parameters`'s `an_evidence_count`), which is the number it came from.
        observations: (warrant != Warrant::Defaulted
            && !(warrant == Warrant::Supplied && observations.count() == 0))
            .then_some(observations),
    }
}

/// Where one `(stratum × slippage group)`'s level came from.
fn level_origin_of(
    provenance: LevelProvenance,
    stratum: Stratum,
    slippage_group: u32,
) -> LevelOrigin {
    LevelOrigin {
        smoothing: match provenance.source {
            // **The curve is not read on this arm, and cannot be lost by not reading it.** This
            // source means the stratum's period had no curve at all — it is what `blend_level`
            // returns for a fitted cell with no curve beside it — so there is nothing to record.
            LevelSource::Cell => LevelSmoothing::ThisStratum,
            LevelSource::Curve => LevelSmoothing::ThisPeriodsCurve {
                curve: the_recorded_curve(provenance.curve, stratum, slippage_group, "the level")
                    .into(),
                reach: the_recorded_reach(provenance.reach, stratum, slippage_group, "the level"),
            },
            LevelSource::Blend { curve_weight } => LevelSmoothing::Blend {
                curve_weight,
                curve: the_recorded_curve(provenance.curve, stratum, slippage_group, "the level")
                    .into(),
                reach: the_recorded_reach(provenance.reach, stratum, slippage_group, "the level"),
            },
        },
        expected_slipped_reads: provenance.slipped_reads,
    }
}

/// Where one of a stratum's two slippage shares came from. `which_number` names the quantity, for
/// the message a contradictory provenance would raise.
fn share_smoothing_of(
    provenance: ShareProvenance,
    stratum: Stratum,
    slippage_group: u32,
    which_number: &str,
) -> ShareSmoothing {
    match provenance.source {
        ShareSource::Stratum => ShareSmoothing::ThisStratum,
        ShareSource::Curve => ShareSmoothing::ThisPeriodsCurve {
            curve: the_recorded_curve(provenance.curve, stratum, slippage_group, which_number)
                .into(),
            reach: the_recorded_reach(provenance.reach, stratum, slippage_group, which_number),
        },
        ShareSource::Blend { curve_weight } => ShareSmoothing::Blend {
            curve_weight,
            curve: the_recorded_curve(provenance.curve, stratum, slippage_group, which_number)
                .into(),
            reach: the_recorded_reach(provenance.reach, stratum, slippage_group, which_number),
        },
    }
}

/// The curve a number came off, refusing a provenance that says a curve supplied it and does not
/// say which.
fn the_recorded_curve<T>(
    curve: Option<T>,
    stratum: Stratum,
    slippage_group: u32,
    which_number: &str,
) -> T {
    curve.unwrap_or_else(|| {
        panic!(
            "at period {}, {} repeats, slippage group {slippage_group}, {which_number} came off \
             its period's curve and no curve was recorded beside it — the file says which curve \
             supplied a number or does not say a curve did, because a reader cannot weigh an \
             interpolation it cannot see",
            stratum.period, stratum.reference_repeats
        )
    })
}

/// Whether the stratum sat inside the curve's fitted range, refusing a curve that does not say.
fn the_recorded_reach(
    reach: Option<FittedCurveReach>,
    stratum: Stratum,
    slippage_group: u32,
    which_number: &str,
) -> CurveReach {
    reach
        .unwrap_or_else(|| {
            panic!(
                "at period {}, {} repeats, slippage group {slippage_group}, {which_number} came \
                 off its period's curve and nothing says whether this stratum sat inside the \
                 curve's fitted range — a number held at a fitted end is wrong in a known \
                 direction, and the file has no way to say so without this",
                stratum.period, stratum.reference_repeats
            )
        })
        .into()
}

// ---------------------------------------------------------------------
// The file's spelling of each of the pre-pass's own words
// ---------------------------------------------------------------------
//
// **Every one of these is an exhaustive match, and that is half the drift guard.** The file's
// enums mirror the pre-pass's so that renaming a variant upstream cannot silently re-interpret a
// file on disk — and a mirror nothing checks goes stale. A variant *added* upstream now fails to
// compile here. **It has already caught one**: `SeedRegime` has four states and the file's
// `SeedRung` was built with three, so a run whose cohort turned out to have no variation at all
// had no rung to be written under.
//
// **Exhaustiveness cannot catch a variant crossed with another** — `Flat => Sloping` compiles and
// spells correctly — which is the other half, and
// `every_pre_pass_word_maps_to_its_own_word_in_the_file` names every pair.

impl From<Provenance> for Warrant {
    fn from(provenance: Provenance) -> Self {
        match provenance {
            Provenance::FittedHere => Self::FittedHere,
            Provenance::Borrowed => Self::Borrowed,
            Provenance::Supplied => Self::Supplied,
            Provenance::Defaulted => Self::Defaulted,
        }
    }
}

impl From<ContaminationSource> for ContaminationFittedFrom {
    fn from(source: ContaminationSource) -> Self {
        match source {
            ContaminationSource::ThisReadGroupsReads => Self::ThisReadGroupsOwnReads,
            ContaminationSource::TheWholeSamplesReads => Self::EveryReadOfThisSample,
        }
    }
}

impl From<SeedRegime> for SeedRung {
    fn from(regime: SeedRegime) -> Self {
        match regime {
            SeedRegime::FittedCurve => Self::FittedCurve,
            SeedRegime::NeutralShape => Self::NeutralShape,
            SeedRegime::ZeroDiversity => Self::ZeroDiversity,
            SeedRegime::FallbackDiversity => Self::StatedHeterozygosity,
        }
    }
}

impl From<FittedCurveReach> for CurveReach {
    fn from(reach: FittedCurveReach) -> Self {
        match reach {
            FittedCurveReach::Inside => Self::InsideTheFittedRange,
            FittedCurveReach::BelowFitted => Self::BelowTheFittedRange,
            FittedCurveReach::AboveFitted => Self::AboveTheFittedRange,
        }
    }
}

impl From<FittedSlippageCurve> for SlippageCurve {
    fn from(curve: FittedSlippageCurve) -> Self {
        Self {
            rise_shape: curve.rise_shape.get(),
            intercept: curve.intercept,
            slope: curve.slope,
            fitted_from_repeats: curve.fitted_from,
            fitted_to_repeats: curve.fitted_to,
            held_out_error: curve.held_out_error,
            cells: curve.cells as u64,
        }
    }
}

impl From<FittedShareCurve> for ShareCurve {
    fn from(curve: FittedShareCurve) -> Self {
        Self {
            shape: curve.shape.into(),
            intercept: curve.intercept,
            slope: curve.slope,
            bend: curve.bend,
            centre_repeats: curve.centre,
            fitted_from_repeats: curve.fitted_from,
            fitted_to_repeats: curve.fitted_to,
            held_out_error: curve.held_out_error,
            strata: curve.strata as u64,
            curve_fitted_on: curve.source.into(),
        }
    }
}

impl From<ShareCurveSource> for ShareCurveRung {
    fn from(source: ShareCurveSource) -> Self {
        match source {
            ShareCurveSource::ThisPeriod => Self::ThisPeriod,
            ShareCurveSource::ThisPeriodUnscored => Self::ThisPeriodUnscored,
            ShareCurveSource::OtherPeriods => Self::OtherPeriods,
            ShareCurveSource::BuiltInDefault => Self::BuiltInDefault,
        }
    }
}

impl From<FittedShape> for ShareShape {
    fn from(shape: FittedShape) -> Self {
        match shape {
            FittedShape::Flat => Self::Flat,
            FittedShape::Sloping => Self::Sloping,
            FittedShape::Turning => Self::Turning,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::genotype_prior::SpectrumSeed;
    use crate::ng::calling::likelihood::ssr::{DEFAULT_OUTLIER_WEIGHT, RepeatTractOutlierWeight};
    use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
    use crate::ng::parameter_estimation::joint::contamination::{
        ContaminationEstimate, NotIdentifiedReason,
    };
    use crate::ng::parameter_estimation::joint::slippage_curve::RiseShape;
    use crate::ng::parameter_estimation::joint::ssr_fit::{
        DerivedStratum, PeriodLengthSpectrum, SharesProvenance, Slippage, StratumFit,
        StratumOutcome, StratumRefusal,
    };
    use crate::ng::parameter_estimation::ssr::{RepeatCount, Stratum as SsrStratum, StratumKey};
    use crate::ng::types::{Ploidy, ReadGroupId, SsrPeriod};
    use std::collections::BTreeSet;

    /// A plant whose name needs escaping and is not ASCII — sample names come from `@RG SM` and
    /// are whatever the sequencing centre typed.
    const AWKWARD_SAMPLE: &str = "Ailsa ‘Craig’ \"×2\"";

    // The one reference every fixture here is fitted against, so that a test varying something
    // else does not look as though it varies this too. The module's own, shared with the round
    // trip in `to_run_parameters`, which writes a file this one's fixtures produced.
    use super::super::tests::THE_REFERENCE_A_RUN_FITTED_AGAINST as A_REFERENCE;

    /// The run's read-group table: **two lanes of one plant and one lane of another**, so the
    /// two per-axis lengths differ (three read groups, two samples) and a row that named its
    /// sample where it should name its library is visible.
    fn a_runs_read_groups() -> ReadGroups {
        ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", AWKWARD_SAMPLE, "lib5"),
        ])
    }

    fn a_fitted_rate(value: f64, provenance: Provenance, observations: u64) -> Estimate<ErrorRate> {
        Estimate {
            value: ErrorRate::try_new(value).expect("a legal error rate"),
            provenance,
            observations,
        }
    }

    /// The accumulator's total for a read group whose reads averaged `mean_reported_error` apiece.
    fn a_read_groups_minted_totals(mean_reported_error: f64, reads: u32) -> MintedReadErrors {
        MintedReadErrors::of_observation(mean_reported_error.ln() * f64::from(reads), reads)
    }

    /// **Three fitted rates whose scales all differ**, so a row written under the wrong read
    /// group cannot pass a value check:
    ///
    /// - read group 0: 0.004 over reads averaging 0.008, so a **scale of 0.5**, fitted;
    /// - read group 1: a rate of **zero**, which `ReadGroupCalibration::from_fitted_rate` refuses
    ///   — so its calibration is the honest defaulted one, scale 1.0, and the 4,242 reads behind
    ///   the refused rate are what the file must **not** write beside that 1.0;
    /// - read group 2: 0.001 over reads averaging 0.004, so a **scale of 0.25**, borrowed.
    fn the_runs_fitted_rates() -> BTreeMap<ReadGroupId, Estimate<ErrorRate>> {
        BTreeMap::from([
            (
                ReadGroupId(0),
                a_fitted_rate(0.004, Provenance::FittedHere, 812_344),
            ),
            (
                ReadGroupId(1),
                a_fitted_rate(0.0, Provenance::FittedHere, 4_242),
            ),
            (
                ReadGroupId(2),
                a_fitted_rate(0.001, Provenance::Borrowed, 640_918),
            ),
        ])
    }

    fn the_runs_minted_totals() -> BTreeMap<ReadGroupId, MintedReadErrors> {
        BTreeMap::from([
            (ReadGroupId(0), a_read_groups_minted_totals(0.008, 1_000)),
            (ReadGroupId(1), a_read_groups_minted_totals(0.008, 1_000)),
            (ReadGroupId(2), a_read_groups_minted_totals(0.004, 1_000)),
        ])
    }

    /// A contamination fraction this read group's own reads produced.
    fn identified(alpha: f64, markers: u64, reads: u64) -> ContaminationEstimate {
        ContaminationEstimate::Estimated {
            alpha,
            source: ContaminationSource::ThisReadGroupsReads,
            panel_markers: 10_000,
            markers_with_reads: markers,
            reads_on_markers: reads,
            leverage: 1.0,
        }
    }

    /// **One read group of each of spec §3.4's three states**: contaminated, not measured, and
    /// measured and found clean.
    fn the_runs_contamination() -> BTreeMap<ReadGroupId, ContaminationEstimate> {
        BTreeMap::from([
            (ReadGroupId(0), identified(0.031, 4_211, 90_233)),
            (
                ReadGroupId(1),
                ContaminationEstimate::NotIdentified {
                    reason: NotIdentifiedReason::TooFewMarkers,
                },
            ),
            (
                ReadGroupId(2),
                ContaminationEstimate::Estimated {
                    alpha: 0.0,
                    source: ContaminationSource::TheWholeSamplesReads,
                    panel_markers: 10_000,
                    markers_with_reads: 2_903,
                    reads_on_markers: 64_118,
                    leverage: 1.0,
                },
            ),
        ])
    }

    /// Two batches — and **not the default**, which is the distinction the file's flag carries.
    /// The two lanes of `TS-1` go together and the other plant alone, because a sample's read
    /// groups cannot be split across batches.
    fn a_declared_batching(read_groups: &ReadGroups) -> DeclaredBatches {
        DeclaredBatches::declared(
            read_groups,
            &[
                BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
                BTreeSet::from([ReadGroupId(2)]),
            ],
        )
        .expect("a partition of this run")
    }

    fn an_inbreeding_estimate(
        value: f64,
        provenance: Provenance,
        observations: u64,
    ) -> Estimate<InbreedingF> {
        Estimate {
            value: InbreedingF::try_new(value).expect("a legal coefficient"),
            provenance,
            observations,
        }
    }

    /// One fitted coefficient and one borrowed one, at counts a factor of **19** apart —
    /// 180,600,412 covered positions against 9,411,027.
    fn the_runs_inbreeding() -> Vec<Estimate<InbreedingF>> {
        vec![
            an_inbreeding_estimate(0.42, Provenance::FittedHere, 180_600_412),
            an_inbreeding_estimate(0.17, Provenance::Borrowed, 9_411_027),
        ]
    }

    fn a_slippage_curve() -> FittedSlippageCurve {
        FittedSlippageCurve {
            rise_shape: RiseShape::new(0.55).expect("a rise shape in [0, 1]"),
            intercept: 0.011,
            slope: 0.004,
            fitted_from: 5,
            fitted_to: 19,
            held_out_error: 0.077,
            cells: 23,
        }
    }

    /// **Its fitted range is 4 to 17 where the slippage curve's is 5 to 19**, so a field taken
    /// from the wrong curve is visible and not merely equal.
    fn a_share_curve() -> FittedShareCurve {
        FittedShareCurve {
            shape: FittedShape::Turning,
            intercept: 1.4,
            slope: -0.09,
            bend: 0.006,
            centre: 11.5,
            fitted_from: 4,
            fitted_to: 17,
            held_out_error: 0.167,
            strata: 12,
            source: ShareCurveSource::ThisPeriod,
        }
    }

    fn ssr_stratum(period: u8, repeats: u32) -> SsrStratum {
        SsrStratum::new(
            SsrPeriod::try_new(usize::from(period)).expect("a motif period"),
            RepeatCount(repeats),
        )
    }

    /// **Four strata, and only the first of them was fitted here.**
    ///
    /// - period 2 at 6 repeats — **fitted**, so it carries a length spectrum of its own. Only
    ///   slippage group 0 put reads in it; group 1 has no numbers there, which is the
    ///   `(stratum × slippage group)` the file must leave out rather than write a zero for. Its
    ///   level is a **blend**, its shorter share is its own, and its fall-off came off the
    ///   period's curve — three different origins in one row.
    /// - period 2 at 11 repeats — **derived** from the period's curves, so it has no length
    ///   spectrum at all, and only group 1 has numbers. Its shares provenance is absent, which
    ///   is the `shorter_share_and_fall_off_origin` key the file must leave out.
    /// - period 3 at 9 repeats — **derived**, group 1 only, and its shares came off a **blend**.
    ///   Its group-0 shares are absent where its group-1 shares are not, so a reader that took
    ///   group 0's shares for every group would write this row's `shorter_share_and_fall_off_origin` as missing.
    /// - period 1 at 30 repeats — **refused**, which contributes nothing at all.
    fn the_runs_slippage() -> StratumFits {
        let fitted = StratumOutcome::Fitted(Box::new(StratumFit {
            stratum: Stratum {
                period: 2,
                reference_repeats: 6,
            },
            slippage: vec![
                Some(Slippage {
                    level: 0.0421,
                    shorter_share: 0.83,
                    fall_off: 0.31,
                }),
                None,
            ],
            length_spectrum: vec![0.1, 0.8, 0.1],
            concentration: 3.5,
            log_likelihood_a_tract: -1.25,
            tracts_fitted: 900,
            borrowed: Vec::new(),
            converged: true,
            tracts_of_its_own: 900,
            reads_crossing: 19_000,
            level_provenance: vec![
                Some(LevelProvenance {
                    source: LevelSource::Blend { curve_weight: 0.37 },
                    curve: Some(a_slippage_curve()),
                    reach: Some(FittedCurveReach::Inside),
                    slipped_reads: Some(8_000.5),
                }),
                None,
            ],
            shares_provenance: vec![
                Some(SharesProvenance {
                    slipped_reads: Some(8_000.5),
                    shorter_share: ShareProvenance {
                        source: ShareSource::Stratum,
                        curve: None,
                        reach: None,
                    },
                    fall_off: ShareProvenance {
                        source: ShareSource::Curve,
                        curve: Some(a_share_curve()),
                        reach: Some(FittedCurveReach::AboveFitted),
                    },
                }),
                None,
            ],
        }));
        let derived = StratumOutcome::Derived(Box::new(DerivedStratum {
            stratum: Stratum {
                period: 2,
                reference_repeats: 11,
            },
            slippage: vec![
                None,
                Some(Slippage {
                    level: 0.0913,
                    shorter_share: 0.79,
                    fall_off: 0.28,
                }),
            ],
            level_provenance: vec![
                None,
                Some(LevelProvenance {
                    source: LevelSource::Curve,
                    curve: Some(a_slippage_curve()),
                    reach: Some(FittedCurveReach::BelowFitted),
                    slipped_reads: None,
                }),
            ],
            shares_provenance: vec![None, None],
            tracts_of_its_own: 12,
            reads_crossing: 140,
        }));
        let blended_shares = StratumOutcome::Derived(Box::new(DerivedStratum {
            stratum: Stratum {
                period: 3,
                reference_repeats: 9,
            },
            slippage: vec![
                None,
                Some(Slippage {
                    level: 0.0617,
                    shorter_share: 0.71,
                    fall_off: 0.22,
                }),
            ],
            level_provenance: vec![
                None,
                Some(LevelProvenance {
                    source: LevelSource::Curve,
                    curve: Some(a_slippage_curve()),
                    reach: Some(FittedCurveReach::Inside),
                    slipped_reads: None,
                }),
            ],
            shares_provenance: vec![
                None,
                Some(SharesProvenance {
                    slipped_reads: Some(31.0),
                    shorter_share: ShareProvenance {
                        source: ShareSource::Blend { curve_weight: 0.6 },
                        curve: Some(a_share_curve()),
                        reach: Some(FittedCurveReach::Inside),
                    },
                    fall_off: ShareProvenance {
                        source: ShareSource::Stratum,
                        curve: None,
                        reach: None,
                    },
                }),
            ],
            tracts_of_its_own: 4,
            reads_crossing: 61,
        }));
        let refused = StratumOutcome::Refused {
            stratum: Stratum {
                period: 1,
                reference_repeats: 30,
            },
            tracts: 3,
            reason: StratumRefusal::BelowTheFloor {
                tracts: 3,
                floor: 50,
            },
        };
        StratumFits::over(
            &[fitted, derived, blended_shares, refused],
            BTreeMap::from([
                (ReadGroupId(0), 0),
                (ReadGroupId(1), 0),
                (ReadGroupId(2), 1),
            ]),
        )
        .with_period_length_spectra(BTreeMap::from([(
            2,
            PeriodLengthSpectrum {
                period: 2,
                length_spectrum: vec![0.15, 0.7, 0.15],
                concentration: 2.75,
                tracts_fitted: 1_400,
                strata_pooled: 2,
                converged: true,
            },
        )]))
    }

    /// **Two rates, and every part of the two keys differs** — different read groups, periods,
    /// repeat counts and **ploidies**, so a row that wrote the run's own ploidy, or its period,
    /// in place of the key's would show. The second is defaulted, so the no-count rule is
    /// exercised on a quantity other than the calibration.
    fn the_runs_substitution_rates() -> BTreeMap<StratumKey, Estimate<ErrorRate>> {
        BTreeMap::from([
            (
                StratumKey {
                    read_group: ReadGroupId(0),
                    stratum: ssr_stratum(2, 6),
                    ploidy: diploid(),
                },
                a_fitted_rate(0.0012, Provenance::Borrowed, 40_122),
            ),
            (
                StratumKey {
                    read_group: ReadGroupId(2),
                    stratum: ssr_stratum(3, 9),
                    ploidy: Ploidy::try_new(4).expect("a tetraploid"),
                },
                a_fitted_rate(0.0007, Provenance::Defaulted, 5),
            ),
        ])
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
    }

    fn a_seed() -> SpectrumSeed {
        SpectrumSeed::new(1.0, 0.0006, SeedRegime::FittedCurve)
    }

    /// The whole run, assembled the way a fitted run assembles it.
    fn a_fitted_run(
        read_groups: &ReadGroups,
        contamination: &BTreeMap<ReadGroupId, ContaminationEstimate>,
    ) -> RunParameters {
        a_run_with(
            read_groups,
            contamination,
            a_declared_batching(read_groups),
            the_runs_slippage(),
            the_runs_substitution_rates(),
        )
    }

    /// The same run with the three things a test may want to vary handed in.
    fn a_run_with(
        read_groups: &ReadGroups,
        contamination: &BTreeMap<ReadGroupId, ContaminationEstimate>,
        batching: DeclaredBatches,
        slippage: StratumFits,
        substitution_rates: BTreeMap<StratumKey, Estimate<ErrorRate>>,
    ) -> RunParameters {
        let _ = read_groups;
        RunParameters::assemble(
            &the_runs_fitted_rates(),
            &the_runs_minted_totals(),
            contamination,
            batching,
            the_runs_inbreeding()
                .iter()
                .map(|estimate| estimate.value)
                .collect(),
            a_seed(),
            slippage,
            substitution_rates,
            diploid(),
        )
    }

    /// The census this module's fixtures name — the module's own, so the shape of a real
    /// identity is what the projection is exercised against.
    fn a_census() -> CensusIdentity {
        super::super::tests::a_census_a_run_could_have_fitted_under()
    }

    /// The file this run and this table write — for the tests that vary one of the two.
    fn projected(run: &RunParameters, read_groups: &ReadGroups) -> ParametersFile {
        ParametersFile::of_run(
            run,
            read_groups,
            &the_runs_fitted_rates(),
            &the_runs_inbreeding(),
            &A_REFERENCE,
            a_census(),
        )
    }

    /// The file the fixture run writes.
    fn the_projected_file() -> ParametersFile {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        projected(&run, &read_groups)
    }

    /// **Every name in the file comes from the run's own read-group table**, and each axis is
    /// written in the order that axis is indexed in.
    #[test]
    fn every_name_in_the_file_comes_from_the_runs_own_table() {
        let file = the_projected_file();

        assert_eq!(file.format_version, FORMAT_VERSION);
        assert_eq!(file.ploidy, 2);
        assert_eq!(
            file.fitted_from.reference_digest, "0123456789abcdef0123456789abcdef",
            "the file spells the run's reference as 32 lower-case hex characters"
        );
        assert_eq!(
            file.fitted_from.samples,
            vec!["TS-1".to_owned(), AWKWARD_SAMPLE.to_owned()],
            "the sample axis is the read-group table's first-seen order, and a plant sequenced \
             from two libraries appears once"
        );
        assert_eq!(file.fitted_from.census, a_census());

        let rows = &file.fitted_from.read_groups;
        assert_eq!(rows.len(), 3);
        assert_eq!(rows[0].read_group, 0);
        assert_eq!(rows[0].declared_id, "HWI.3");
        assert_eq!(rows[0].library, "lib3");
        assert_eq!(rows[0].sample, "TS-1");
        assert_eq!(
            (rows[1].library.as_str(), rows[1].sample.as_str()),
            ("lib4", "TS-1"),
            "read groups 0 and 1 are two libraries of one plant, which is the grain the \
             contamination fraction exists at"
        );
        assert_eq!(rows[2].sample, AWKWARD_SAMPLE);

        assert_eq!(
            file.inbreeding
                .by_sample
                .iter()
                .map(|row| row.sample.as_str())
                .collect::<Vec<_>>(),
            vec!["TS-1", AWKWARD_SAMPLE]
        );
        assert_eq!(
            file.sequencing_batches
                .by_sample
                .iter()
                .map(|row| (row.sample.as_str(), row.batch))
                .collect::<Vec<_>>(),
            vec![("TS-1", 0), (AWKWARD_SAMPLE, 1)]
        );
        assert_eq!(
            file.sequencing_batches
                .by_read_group
                .iter()
                .map(|row| (row.read_group, row.batch))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 0), (2, 1)],
            "the two batching axes are different lengths — three read groups over two samples — \
             so an exchange of the two cannot go unnoticed"
        );
        assert!(
            file.sequencing_batches.batching_was_declared,
            "this run declared its batching"
        );
    }

    /// **A run that declared no batching says so**, and it is the only thing that says so: the
    /// rows a defaulted batching writes are the rows a declaration of one batch holding every
    /// library would write.
    #[test]
    fn a_run_that_declared_no_batching_writes_the_flag_false() {
        let read_groups = a_runs_read_groups();
        let run = a_run_with(
            &read_groups,
            &the_runs_contamination(),
            DeclaredBatches::all_together(&read_groups),
            the_runs_slippage(),
            the_runs_substitution_rates(),
        );
        let file = projected(&run, &read_groups);

        assert!(!file.sequencing_batches.batching_was_declared);
        assert_eq!(
            file.sequencing_batches
                .by_read_group
                .iter()
                .map(|row| row.batch)
                .collect::<Vec<_>>(),
            vec![0, 0, 0],
            "these are the same rows a declaration of one batch holding every library would \
             write, which is why the flag is the only thing that separates them"
        );
    }

    /// **The calibration's evidence count is the one thing `RunParameters` cannot supply**, and
    /// the file carries it in reads.
    ///
    /// `ReadGroupCalibration` is a multiplier and a warrant: the count of reads behind the fitted
    /// rate lives on the `Estimate<ErrorRate>` that assembly reads and drops. A projection
    /// written from the assembled parameters alone would have to leave this absent.
    #[test]
    fn the_calibration_carries_the_count_the_assembled_parameters_dropped() {
        let file = the_projected_file();
        let fitted = &file.base_quality_calibration.by_read_group[0];

        assert_eq!(fitted.read_group, 0);
        // **A half, to the accumulator's quantum rather than exactly.** The multiplier is the
        // fitted rate over the reads' own mean reported error, and that denominator is summed in
        // fixed point so that shards merged in different orders give the same number — which
        // leaves the ratio a few parts in ten billion off the real one. Measured on this
        // fixture: 0.4999999998066437, which is 3.9 in 10^10 below a half, against the bound of
        // 4.8 in 10^7 `ReadGroupCalibration` states for the mean it comes from.
        assert!(
            (fitted.error_probability_multiplier.value / 0.5 - 1.0).abs() < 5e-7,
            "a fitted rate of 0.004 over reads averaging 0.008 is a multiplier of a half, and \
             this one is {}",
            fitted.error_probability_multiplier.value
        );
        assert_eq!(
            fitted.error_probability_multiplier.warrant,
            Warrant::FittedHere
        );
        assert_eq!(
            fitted.error_probability_multiplier.observations,
            Some(EvidenceCount::Reads(812_344)),
            "the count is the fitted rate's own, and it names its unit in the file"
        );

        let borrowed = &file.base_quality_calibration.by_read_group[2];
        assert_eq!(
            borrowed.error_probability_multiplier.warrant,
            Warrant::Borrowed,
            "the scale carries the rate's warrant: a rate borrowed from a sibling read group \
             makes a borrowed calibration, and stamping it fitted would launder it"
        );
        assert_eq!(
            borrowed.error_probability_multiplier.observations,
            Some(EvidenceCount::Reads(640_918))
        );
    }

    /// **A defaulted value carries no evidence count, and the fit's count is what it must not
    /// carry.**
    ///
    /// Read group 1's rate fitted to zero and was refused, so its multiplier is the honest
    /// defaulted one. The refused rate still has 4,242 reads on it; writing that count here
    /// would say a multiplier of one rests on 4,242 reads, and it rests on nothing.
    #[test]
    fn a_defaulted_number_carries_no_evidence_count() {
        let file = the_projected_file();
        let defaulted = &file.base_quality_calibration.by_read_group[1];

        assert_eq!(defaulted.error_probability_multiplier.value, 1.0);
        assert_eq!(
            defaulted.error_probability_multiplier.warrant,
            Warrant::Defaulted,
            "a multiplier of exactly one is a legitimate fitted answer as well as the default, \
             which is why the warrant travels beside it rather than being read off the value"
        );
        assert_eq!(
            defaulted.error_probability_multiplier.observations, None,
            "the rate this read group's calibration was refused from carries 4,242 reads, and \
             none of them stand behind the stated one that replaced it"
        );

        let defaulted_rate = file
            .repeat_tracts
            .substitution_rate_by_stratum
            .iter()
            .find(|row| row.rate.warrant == Warrant::Defaulted)
            .expect("the fixture defaults one substitution rate");
        assert_eq!(
            defaulted_rate.rate.observations, None,
            "the rule is the file's and not the calibration's: no defaulted number carries a count"
        );
    }

    /// **A *fitted* count of zero is a count, and the warrant is what makes one absent** — with
    /// the one exception step E2 added.
    ///
    /// *Fitted from nothing* and *a stated constant* are two different claims, and a rule that
    /// dropped every zero count would collapse them: a fit that produced a number from no reads
    /// at all is alarming, and the count is what says so.
    ///
    /// **The exception is a `supplied` number with nothing behind it**, which is what an operator
    /// who declares an inbreeding coefficient produces. Writing
    /// `observations = { covered_positions = 0 }` beside it says the number was measured over no
    /// genome, while the file's own editing rule tells that same reader to delete `observations`
    /// on a value they supplied — so a produced file contradicted its own instructions. A
    /// `supplied` number with a real count still keeps it, which is the demoted-file case spec
    /// §2.1 creates.
    #[test]
    fn an_evidence_count_of_zero_is_written_and_is_not_absence() {
        assert_eq!(
            warranted_value(0.5, Provenance::FittedHere, EvidenceCount::Reads(0)).observations,
            Some(EvidenceCount::Reads(0))
        );
        assert_eq!(
            warranted_value(
                0.5,
                Provenance::Supplied,
                EvidenceCount::CoveredPositions(u64::MAX)
            )
            .observations,
            Some(EvidenceCount::CoveredPositions(u64::MAX)),
            "a supplied number keeps its count: spec §2.1 demotes a whole mismatched file to \
             supplied, and those numbers were fitted on some cohort"
        );
        assert_eq!(
            warranted_value(
                0.9,
                Provenance::Supplied,
                EvidenceCount::CoveredPositions(0)
            )
            .observations,
            None,
            "a coefficient an operator typed has nothing behind it, and a zero count beside it \
             claims a measurement over no genome"
        );
        assert_eq!(
            warranted_value(1.0, Provenance::Defaulted, EvidenceCount::Reads(4_242)).observations,
            None
        );
    }

    /// **Each sample's inbreeding warrant and count come from the fit**, which is the second
    /// thing the seam into calling drops.
    ///
    /// `RunParameters` holds a bare `Vec<InbreedingF>`. A file written from that alone could say
    /// nothing about where a coefficient came from, so every sample would read as handed over
    /// when the run fitted it.
    #[test]
    fn each_samples_inbreeding_carries_the_warrant_the_seam_dropped() {
        let file = the_projected_file();
        let rows = &file.inbreeding.by_sample;

        assert_eq!(rows[0].inbreeding_coefficient.value, 0.42);
        assert_eq!(
            rows[0].inbreeding_coefficient.warrant,
            Warrant::FittedHere,
            "this sample's coefficient was fitted from its own runs of homozygosity"
        );
        assert_eq!(
            rows[0].inbreeding_coefficient.observations,
            Some(EvidenceCount::CoveredPositions(180_600_412)),
            "an inbreeding coefficient is fitted over covered reference positions, not over \
             windows and not over reads"
        );

        assert_eq!(rows[1].inbreeding_coefficient.value, 0.17);
        assert_eq!(rows[1].inbreeding_coefficient.warrant, Warrant::Borrowed);
        assert_eq!(
            rows[1].inbreeding_coefficient.observations,
            Some(EvidenceCount::CoveredPositions(9_411_027))
        );
    }

    /// **A read group that identified nothing gets a row with no measurement**, where in memory
    /// it carries a fraction of zero beside two zero counts.
    ///
    /// This is spec §5's second row and the one a reader is most likely to collapse: *measured
    /// and found clean* and *not measured* are both a fraction near zero, and only the evidence
    /// tells them apart. The fixture holds one of each.
    #[test]
    fn an_unmeasured_read_group_writes_no_measurement_and_a_clean_one_writes_a_zero() {
        let file = the_projected_file();
        let rows = &file
            .contamination
            .as_ref()
            .expect("some read group identified a fraction, so the section is written")
            .by_read_group;

        assert_eq!(
            rows.len(),
            3,
            "where some read group identified a fraction, every read group needs a row"
        );

        let contaminated = rows[0]
            .measurement
            .as_ref()
            .expect("read group 0 identified a fraction");
        assert_eq!(contaminated.share_of_reads_from_another_sample, 0.031);
        assert_eq!(contaminated.markers_with_reads, 4_211);
        assert_eq!(contaminated.reads_on_markers, 90_233);
        assert_eq!(
            contaminated.fitted_from_reads_of,
            ContaminationFittedFrom::ThisReadGroupsOwnReads
        );
        assert_eq!(rows[0].library, "lib3");

        assert_eq!(
            rows[1].measurement, None,
            "read group 1 identified nothing, and in memory that is a fraction of zero beside \
             two zero counts — which written out would claim a measurement of no contamination"
        );
        assert_eq!(
            rows[1].library, "lib4",
            "an unmeasured row still names its library, because it is the run's second lane of \
             the same plant and a reader has to be able to find it"
        );

        let clean = rows[2]
            .measurement
            .as_ref()
            .expect("read group 2 was measured and found clean");
        assert_eq!(clean.share_of_reads_from_another_sample, 0.0);
        assert_eq!(clean.markers_with_reads, 2_903);
        assert_eq!(clean.reads_on_markers, 64_118);
        assert_eq!(
            clean.fitted_from_reads_of,
            ContaminationFittedFrom::EveryReadOfThisSample
        );
    }

    /// **A run where nobody identified anything writes no contamination section at all** — spec
    /// §5's first row, and the only one of the five expressed by a whole section going missing.
    ///
    /// The same run, with the same three read groups, differing only in that no fraction was
    /// identified: the section disappears rather than becoming three unmeasured rows.
    #[test]
    fn an_uncontaminated_run_writes_no_contamination_section() {
        let read_groups = a_runs_read_groups();
        let nothing_identified: BTreeMap<ReadGroupId, ContaminationEstimate> = (0..3)
            .map(|group| {
                (
                    ReadGroupId(group),
                    ContaminationEstimate::NotIdentified {
                        reason: NotIdentifiedReason::NoPanel,
                    },
                )
            })
            .collect();
        let run = a_fitted_run(&read_groups, &nothing_identified);
        let file = projected(&run, &read_groups);

        assert_eq!(
            file.contamination, None,
            "the read likelihood's plain formula is what an absent table asks for; three rows of \
             no measurement would ask for the mixture path with every fraction zero"
        );
        assert_eq!(
            file.sequencing_batches.by_read_group.len(),
            3,
            "the batching is written even where no contamination was fitted, because it is a \
             fact about the run rather than about the fit"
        );
    }

    /// **A `(stratum × slippage group)` with no reads gets no row, and a stratum with no fit of
    /// its own gets no length spectrum** — two more of spec §5's five states.
    ///
    /// The fixture's four strata would fill eight `(stratum × group)` cells if the file wrote a
    /// row for every pair; three of them have numbers.
    #[test]
    fn a_pair_with_no_reads_and_a_stratum_with_no_fit_get_no_row() {
        let file = the_projected_file();
        let tracts = &file.repeat_tracts;

        assert_eq!(
            tracts
                .slippage_by_stratum_and_group
                .iter()
                .map(|row| (row.period, row.reference_repeats, row.slippage_group))
                .collect::<Vec<_>>(),
            vec![(2, 6, 0), (2, 11, 1), (3, 9, 1)],
            "slippage group 1 put no read in the stratum at 6 repeats and group 0 none in the \
             two derived ones; the refused stratum at period 1 has no answer for anybody"
        );
        assert_eq!(
            tracts
                .length_spectrum_by_stratum
                .iter()
                .map(|row| (row.period, row.reference_repeats))
                .collect::<Vec<_>>(),
            vec![(2, 6)],
            "only the stratum fitted on its own tracts has a spectrum; the derived ones carry \
             none by construction, and that absence is what the middle rung answers"
        );
        assert_eq!(tracts.length_spectrum_by_stratum[0].concentration, 3.5);
        assert_eq!(
            tracts.length_spectrum_by_stratum[0].shares_by_repeat_offset,
            vec![0.1, 0.8, 0.1]
        );
        assert_eq!(
            tracts
                .length_spectrum_by_period
                .iter()
                .map(|row| (row.period, row.concentration))
                .collect::<Vec<_>>(),
            vec![(2, 2.75)],
            "the middle rung is present only where the run asked for it"
        );
        assert_eq!(
            tracts.length_spectrum_by_period[0].shares_by_repeat_offset,
            vec![0.15, 0.7, 0.15]
        );
        assert_eq!(
            tracts
                .slippage_group_by_read_group
                .iter()
                .map(|row| (row.read_group, row.slippage_group))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 0), (2, 1)],
            "the map is the run's own declaration and is not the identity"
        );
    }

    /// **A slippage number that came off a curve carries the curve and the reach; one that is
    /// the stratum's own carries neither** — and the three numbers of one row can have three
    /// different origins.
    #[test]
    fn a_number_off_a_curve_carries_the_curve_and_a_stratums_own_carries_neither() {
        let file = the_projected_file();
        let blended = &file.repeat_tracts.slippage_by_stratum_and_group[0];

        assert_eq!(blended.share_of_reads_that_slip, 0.0421);
        assert_eq!(blended.shorter_share, 0.83);
        assert_eq!(blended.fall_off, 0.31);
        assert_eq!(
            blended
                .share_of_reads_that_slip_origin
                .expected_slipped_reads,
            Some(8_000.5)
        );
        match &blended.share_of_reads_that_slip_origin.smoothing {
            LevelSmoothing::Blend {
                curve_weight,
                reach,
                ..
            } => {
                assert_eq!(*curve_weight, 0.37);
                assert_eq!(*reach, CurveReach::InsideTheFittedRange);
            }
            other => panic!("this stratum's level was blended, not {other:?}"),
        }

        let shares = blended
            .shorter_share_and_fall_off_origin
            .as_ref()
            .expect("this stratum has a shares provenance");
        assert_eq!(shares.expected_slipped_reads, Some(8_000.5));
        assert_eq!(
            shares.shorter_share_smoothing,
            ShareSmoothing::ThisStratum,
            "a stratum can keep its own direction split while its level takes a curve"
        );
        match &shares.fall_off_smoothing {
            ShareSmoothing::ThisPeriodsCurve { reach, .. } => {
                assert_eq!(*reach, CurveReach::AboveTheFittedRange);
            }
            other => panic!("this stratum's fall-off came off its period's curve, not {other:?}"),
        }

        let derived = &file.repeat_tracts.slippage_by_stratum_and_group[1];
        match &derived.share_of_reads_that_slip_origin.smoothing {
            LevelSmoothing::ThisPeriodsCurve { reach, .. } => {
                assert_eq!(*reach, CurveReach::BelowTheFittedRange)
            }
            other => panic!("a derived stratum's level is its period's curve, not {other:?}"),
        }
        assert_eq!(
            derived
                .share_of_reads_that_slip_origin
                .expected_slipped_reads,
            None,
            "a stratum that borrowed has reads of its own but no level of its own to say how \
             many of them slipped — and absent is not zero"
        );
        assert_eq!(
            derived.shorter_share_and_fall_off_origin, None,
            "no shares provenance was recorded for this pair, so the file writes no key"
        );

        // The third row's shares are its group's own, where the row above it has none at all —
        // so a reader that took slippage group 0's shares for every group would write this one
        // as missing.
        let blended_shares = &file.repeat_tracts.slippage_by_stratum_and_group[2];
        let shares = blended_shares
            .shorter_share_and_fall_off_origin
            .as_ref()
            .expect("slippage group 1 has a shares provenance at this stratum");
        assert_eq!(shares.expected_slipped_reads, Some(31.0));
        match &shares.shorter_share_smoothing {
            ShareSmoothing::Blend {
                curve_weight,
                reach,
                ..
            } => {
                assert_eq!(
                    *curve_weight, 0.6,
                    "the weight is the share the curve carried, not the stratum's"
                );
                assert_eq!(*reach, CurveReach::InsideTheFittedRange);
            }
            other => panic!("this share was blended, not {other:?}"),
        }
        assert_eq!(shares.fall_off_smoothing, ShareSmoothing::ThisStratum);
    }

    /// **Every field of both curves reaches the file under its own name.**
    ///
    /// The two `From` impls rename field by field — `fitted_from` becomes
    /// `fitted_from_repeats`, `centre` becomes `centre_repeats` — which is where a transposition
    /// happens and where nothing but an assertion can see it. An intercept written under `slope`
    /// evaluates to a different number at every repeat count, silently.
    #[test]
    fn every_field_of_both_curves_reaches_the_file_under_its_own_name() {
        let file = the_projected_file();
        let blended = &file.repeat_tracts.slippage_by_stratum_and_group[0];

        let LevelSmoothing::Blend { curve, .. } =
            &blended.share_of_reads_that_slip_origin.smoothing
        else {
            panic!("this stratum's level was blended")
        };
        assert_eq!(
            curve.rise_shape, 0.55,
            "the rise shape is a checked newtype upstream and a bare number in the file"
        );
        assert_eq!(curve.intercept, 0.011);
        assert_eq!(curve.slope, 0.004);
        assert_eq!(curve.fitted_from_repeats, 5);
        assert_eq!(curve.fitted_to_repeats, 19);
        assert_eq!(curve.held_out_error, 0.077);
        assert_eq!(curve.cells, 23);

        let shares = blended
            .shorter_share_and_fall_off_origin
            .as_ref()
            .expect("a shares provenance");
        let ShareSmoothing::ThisPeriodsCurve { curve, .. } = &shares.fall_off_smoothing else {
            panic!("this stratum's fall-off came off its period's curve")
        };
        assert_eq!(curve.shape, ShareShape::Turning);
        assert_eq!(curve.intercept, 1.4);
        assert_eq!(curve.slope, -0.09);
        assert_eq!(curve.bend, 0.006);
        assert_eq!(curve.centre_repeats, 11.5);
        assert_eq!(
            (curve.fitted_from_repeats, curve.fitted_to_repeats),
            (4, 17),
            "the share curve's fitted range is not the slippage curve's, so a field taken from \
             the wrong curve shows here"
        );
        assert_eq!(curve.held_out_error, 0.167);
        assert_eq!(curve.strata, 12);
        assert_eq!(curve.curve_fitted_on, ShareCurveRung::ThisPeriod);
    }

    /// **A level that is the stratum's own carries no curve** — the arm reached at a motif period
    /// the run drew no curve for, which is the low-coverage corner where a period has too few
    /// strata to fit one.
    #[test]
    fn a_level_that_is_the_stratums_own_carries_no_curve() {
        let origin = level_origin_of(
            LevelProvenance {
                source: LevelSource::Cell,
                curve: None,
                reach: None,
                slipped_reads: Some(412.0),
            },
            Stratum {
                period: 2,
                reference_repeats: 6,
            },
            0,
        );
        assert_eq!(origin.smoothing, LevelSmoothing::ThisStratum);
        assert_eq!(origin.expected_slipped_reads, Some(412.0));
    }

    /// **The stated concentration says whether the run fitted it**, which is the only thing that
    /// can tell the run's own median from the stated constant when the two are the same number.
    #[test]
    fn the_stated_concentration_says_whether_the_run_fitted_it() {
        let file = the_projected_file();
        assert_eq!(
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .value,
            3.5,
            "one stratum was fitted, so the run's median is that stratum's own concentration"
        );
        assert_eq!(
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .warrant,
            Warrant::FittedHere
        );
        assert_eq!(
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .observations,
            None,
            "a median over strata is not an estimate with a sample size, so no count stands \
             behind it whether the run had strata to take a median over or not"
        );

        let read_groups = a_runs_read_groups();
        let with_no_strata = a_run_with(
            &read_groups,
            &the_runs_contamination(),
            a_declared_batching(&read_groups),
            StratumFits::over(&[], BTreeMap::new()),
            BTreeMap::new(),
        );
        let nothing_fitted = projected(&with_no_strata, &read_groups);
        let stated = &nothing_fitted
            .repeat_tracts
            .fallback_length_spectrum_concentration;
        assert_eq!(
            stated.warrant,
            Warrant::Defaulted,
            "a run that fitted no stratum states the compiled-in concentration"
        );
        assert_eq!(
            stated.value, 1.0,
            "and the median of a run that fitted strata can land on exactly this number, which \
             is why the warrant is what separates them"
        );
        assert_eq!(stated.observations, None);
        assert!(
            nothing_fitted
                .repeat_tracts
                .slippage_group_by_read_group
                .is_empty(),
            "no read group was declared into a slippage group, so no row names one"
        );
    }

    /// **The substitution rate names all three parts of its key and counts bases compared.**
    ///
    /// Not reads: a read spanning a stratum contributes one read and as many bases as it crosses.
    #[test]
    fn the_substitution_rate_names_its_key_and_counts_bases_compared() {
        let file = the_projected_file();
        let rows = &file.repeat_tracts.substitution_rate_by_stratum;

        assert_eq!(rows.len(), 2);
        assert_eq!(
            (
                rows[0].read_group,
                rows[0].period,
                rows[0].reference_repeats,
                rows[0].ploidy
            ),
            (0, 2, 6, 2)
        );
        assert_eq!(rows[0].rate.value, 0.0012);
        assert_eq!(rows[0].rate.warrant, Warrant::Borrowed);
        assert_eq!(
            rows[0].rate.observations,
            Some(EvidenceCount::BasesCompared(40_122))
        );
        assert_eq!(
            (
                rows[1].read_group,
                rows[1].period,
                rows[1].reference_repeats,
                rows[1].ploidy
            ),
            (2, 3, 9, 4),
            "the rate is per read group as well as per stratum, because how often a base \
             misreads inside a tract is a property of the chemistry — and the ploidy is the \
             key's, not the run's, which is 2"
        );
    }

    /// **Every regime the seed can be built under has a rung in the file**, including the one an
    /// earlier draft of `SeedRung` had no word for.
    ///
    /// A cohort with no variation at all gets a legal pair of concentrations that says nothing
    /// about how it was arrived at, so the rung is the only place that state is recorded.
    #[test]
    fn every_seed_regime_has_a_rung_in_the_file() {
        assert_eq!(
            SeedRung::from(SeedRegime::FittedCurve),
            SeedRung::FittedCurve
        );
        assert_eq!(
            SeedRung::from(SeedRegime::NeutralShape),
            SeedRung::NeutralShape
        );
        assert_eq!(
            SeedRung::from(SeedRegime::ZeroDiversity),
            SeedRung::ZeroDiversity
        );
        assert_eq!(
            SeedRung::from(SeedRegime::FallbackDiversity),
            SeedRung::StatedHeterozygosity,
            "the file names this rung after what it rests on — a stated species-range \
             heterozygosity taken from human data"
        );

        let file = the_projected_file();
        assert_eq!(file.ordinary_site_prior.rung, SeedRung::FittedCurve);
        assert_eq!(file.ordinary_site_prior.reference_concentration, 1.0);
        assert_eq!(
            file.ordinary_site_prior.alternative_concentration_total,
            0.0006
        );
    }

    /// **Every word the pre-pass uses maps to its own word in the file, and not to a neighbour's.**
    ///
    /// The exhaustive matches catch a variant *added* upstream at compile time and cannot catch
    /// one *crossed* with another: `Flat => Sloping` compiles, spells as `"sloping"`, and passes
    /// every other test in this crate. Half of these pairs are reached by no fixture — among them
    /// `ShareShape::Flat`, which the project's own measurement makes the commonest real answer.
    #[test]
    fn every_pre_pass_word_maps_to_its_own_word_in_the_file() {
        for (upstream, in_the_file) in [
            (Provenance::FittedHere, Warrant::FittedHere),
            (Provenance::Borrowed, Warrant::Borrowed),
            (Provenance::Supplied, Warrant::Supplied),
            (Provenance::Defaulted, Warrant::Defaulted),
        ] {
            assert_eq!(Warrant::from(upstream), in_the_file, "{upstream:?}");
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
                ContaminationFittedFrom::from(upstream),
                in_the_file,
                "{upstream:?}"
            );
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
            assert_eq!(CurveReach::from(upstream), in_the_file, "{upstream:?}");
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
            assert_eq!(ShareCurveRung::from(upstream), in_the_file, "{upstream:?}");
        }
        for (upstream, in_the_file) in [
            (FittedShape::Flat, ShareShape::Flat),
            (FittedShape::Sloping, ShareShape::Sloping),
            (FittedShape::Turning, ShareShape::Turning),
        ] {
            assert_eq!(ShareShape::from(upstream), in_the_file, "{upstream:?}");
        }
    }

    /// **The repeat-tract outlier weight is written out, marked as the guess it is — and a
    /// supplied one is written as supplied.**
    ///
    /// The second half is what stops the round trip losing an edit. Spec §3.8 puts this number
    /// in the file so that a person can change it, and until 2026-08-30 the projection wrote the
    /// compiled-in constant whatever the run held, so an edited file came back `defaulted` at
    /// 0.01.
    #[test]
    fn the_outlier_weight_is_written_as_the_inherited_guess_it_is() {
        let file = the_projected_file();
        assert_eq!(
            file.stated_constants.repeat_tract_outlier_weight.value,
            DEFAULT_OUTLIER_WEIGHT
        );
        assert_eq!(
            file.stated_constants.repeat_tract_outlier_weight.warrant,
            Warrant::Defaulted,
            "it is inherited from the existing caller and was never measured here; a number in \
             a file the user can edit is a number the project has admitted is a guess"
        );
        assert_eq!(
            file.stated_constants
                .repeat_tract_outlier_weight
                .observations,
            None
        );

        let read_groups = a_runs_read_groups();
        let supplied = a_fitted_run(&read_groups, &the_runs_contamination())
            .with_repeat_tract_outlier_weight(RepeatTractOutlierWeight::supplied(0.04));
        let file = projected(&supplied, &read_groups);
        assert_eq!(
            file.stated_constants.repeat_tract_outlier_weight.value,
            0.04
        );
        assert_eq!(
            file.stated_constants.repeat_tract_outlier_weight.warrant,
            Warrant::Supplied,
            "a weight the run was handed is not one it inherited, and the file is where that \
             difference has to survive"
        );
        assert_eq!(
            file.stated_constants
                .repeat_tract_outlier_weight
                .observations,
            None,
            "nothing counted anything to arrive at it, whichever way it got here"
        );
    }

    /// **The smallest legal run projects** — one sample, one library, no strata, no
    /// contamination. `CLAUDE.md` makes this a first-class case, and it is the only shape in
    /// which the read-group and sample axes are the same length, so it is where a projection that
    /// worked only because they differ would show.
    #[test]
    fn of_run_writes_a_single_sample_single_read_group_run() {
        let read_groups = ReadGroups::of_lanes(&[("HWI.3", AWKWARD_SAMPLE, "lib3")]);
        let rates = BTreeMap::from([(
            ReadGroupId(0),
            a_fitted_rate(0.004, Provenance::FittedHere, 812_344),
        )]);
        let totals = BTreeMap::from([(ReadGroupId(0), a_read_groups_minted_totals(0.008, 1_000))]);
        let one_coefficient = vec![an_inbreeding_estimate(
            0.42,
            Provenance::FittedHere,
            180_600_412,
        )];
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            DeclaredBatches::all_together(&read_groups),
            one_coefficient
                .iter()
                .map(|estimate| estimate.value)
                .collect(),
            a_seed(),
            StratumFits::over(&[], BTreeMap::new()),
            BTreeMap::new(),
            diploid(),
        );

        let file = ParametersFile::of_run(
            &run,
            &read_groups,
            &rates,
            &one_coefficient,
            &A_REFERENCE,
            a_census(),
        );

        assert_eq!(file.fitted_from.samples, vec![AWKWARD_SAMPLE.to_owned()]);
        assert_eq!(file.base_quality_calibration.by_read_group.len(), 1);
        assert_eq!(file.inbreeding.by_sample.len(), 1);
        assert_eq!(
            file.contamination, None,
            "one sample has no panel to be surprised by, so contamination is not estimable at all"
        );
        assert!(!file.sequencing_batches.batching_was_declared);
        assert_eq!(
            file.repeat_tracts
                .fallback_length_spectrum_concentration
                .warrant,
            Warrant::Defaulted,
            "a run with no repeat tract fitted states the compiled-in concentration"
        );
        assert!(file.repeat_tracts.slippage_by_stratum_and_group.is_empty());
    }

    /// **⚑ A run whose error rate itself was defaulted writes a `defaulted` multiplier that is
    /// not 1.0, and its own reader accepts the file.**
    ///
    /// **The state, and how a run reaches it.**
    /// [`resolve_error_rates`](crate::ng::parameter_estimation::generic::fallback::resolve_error_rates)
    /// hands a read group with too few sites to fit, no sibling above the floor to borrow from
    /// and nothing supplied the pre-pass's own constant —
    /// [`DEFAULT_ERROR_RATE`](crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE) at
    /// 0.001, marked `Defaulted`. `ReadGroupCalibration::from_fitted_rate` then copies **the
    /// rate's** warrant onto `rate / mean minted error`, so what is written is `defaulted` beside
    /// a multiplier of 0.001 over that library's own average. That is the low-data corner
    /// `CLAUDE.md` commits this caller to: one gene, a panel, a shallow single sample.
    ///
    /// **Ruled intended by the owner, 2026-08-31**, and the fixture's two libraries are chosen to
    /// show the direction that makes it so. A library's real error rate is never its reported
    /// sequencing quality — the quality scores describe base calling, and the reads also carry
    /// mismapping, chimeras and damage — so a read group the fit could not measure is charged a
    /// stated rate rather than taken at its word. **On any real library that is the conservative
    /// direction**: the two here report 2.5 × 10⁻⁴ and 5 × 10⁻⁴, either side of HG002's measured
    /// 2.9055 × 10⁻⁴ (`read_likelihoods.md` §3.2), and get multipliers of 4.0 and 2.0 — every read
    /// charged four times and twice worse than it claimed. **An earlier fixture reported 0.008,
    /// which is Q21 and unlike any library in this project**, and so showed the multiplier below
    /// one and the reads made *more* confident. Spec §5's third row says such a read group gets
    /// "scale 1.0" and is the sentence to correct.
    ///
    /// **What this guards.** Step E1 added a rung to `validate` reading *a `defaulted` multiplier
    /// is [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER)
    /// and nothing else*, on the model of the outlier weight's and the fallback concentration's,
    /// and it **refused the file this test's run writes**. The rung is gone and this is what stops
    /// it coming back.
    #[test]
    fn a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts() {
        let read_groups =
            ReadGroups::of_lanes(&[("HWI.3", "TS-1", "lib3"), ("HWI.4", AWKWARD_SAMPLE, "lib4")]);
        // The pre-pass's bottom rung, for both libraries: nothing could be fitted and nothing was
        // supplied, so each takes the stated constant.
        let defaulted_rate = crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE;
        let rates = BTreeMap::from([
            (
                ReadGroupId(0),
                a_fitted_rate(defaulted_rate, Provenance::Defaulted, 0),
            ),
            (
                ReadGroupId(1),
                a_fitted_rate(defaulted_rate, Provenance::Defaulted, 0),
            ),
        ]);
        // Two libraries at real reported qualities — about Q36 and Q33, either side of HG002's
        // measured 2.9055e-4 — so the two multipliers differ, both are above one, and the
        // fixture cannot be read as saying a defaulted rate makes reads cleaner.
        let totals = BTreeMap::from([
            (ReadGroupId(0), a_read_groups_minted_totals(2.5e-4, 1_000)),
            (ReadGroupId(1), a_read_groups_minted_totals(5e-4, 1_000)),
        ]);
        let coefficients = vec![
            an_inbreeding_estimate(0.42, Provenance::FittedHere, 180_600_412),
            an_inbreeding_estimate(0.11, Provenance::FittedHere, 180_600_412),
        ];
        let run = RunParameters::assemble(
            &rates,
            &totals,
            &BTreeMap::new(),
            DeclaredBatches::all_together(&read_groups),
            coefficients.iter().map(|estimate| estimate.value).collect(),
            a_seed(),
            StratumFits::over(&[], BTreeMap::new()),
            BTreeMap::new(),
            diploid(),
        );

        let file = ParametersFile::of_run(
            &run,
            &read_groups,
            &rates,
            &coefficients,
            &A_REFERENCE,
            a_census(),
        );

        for row in &file.base_quality_calibration.by_read_group {
            assert_eq!(
                row.error_probability_multiplier.warrant,
                Warrant::Defaulted,
                "the rate's warrant travels onto the multiplier built from it"
            );
        }
        // 0.001 over reads reporting 2.5e-4 and 5e-4 — and **above one on both**, which is the
        // half of this that says the ruling is the safe direction rather than merely a decision.
        let multipliers: Vec<f64> = file
            .base_quality_calibration
            .by_read_group
            .iter()
            .map(|row| row.error_probability_multiplier.value)
            .collect();
        // **Compared relatively, against the accumulator's own quantum.** `MintedReadErrors` sums
        // the per-read log error in fixed point in units of 2⁻²⁰ nats, and its documentation
        // bounds the resulting relative miss on the mean at 2⁻²¹ ≈ 4.8 × 10⁻⁷; measured here the
        // multipliers come back 3.9999999984 and 1.9999999992, four parts in ten thousand million.
        // An absolute 1e-9 tolerance rejects them, which is a test asserting the fixed point
        // rather than the calibration.
        for (multiplier, expected) in multipliers.iter().zip([4.0, 2.0]) {
            assert!(
                (multiplier - expected).abs() / expected < 1e-6,
                "a defaulted rate over a library's own minted mean, which is not one: \
                 {multipliers:?}"
            );
        }
        assert!(
            multipliers.iter().all(|multiplier| *multiplier > 1.0),
            "a library reporting better than the stated rate is charged worse than it claimed, \
             not better: {multipliers:?}"
        );

        file.validate()
            .expect("this caller must not refuse a file it has just written");
    }

    /// **What the projection writes is a legal file** — it parses back to the same value through
    /// the shape's own serde derives.
    ///
    /// Not the round trip step C4 owes: that one starts from a `RunParameters` off a real fit and
    /// comes back to a `RunParameters`, and would catch a wrong value. This one only says that a
    /// projected file is expressible as TOML at all, which is what step B2 will render.
    #[test]
    fn the_projected_file_is_expressible_as_toml() {
        let written = the_projected_file();
        let text = toml::to_string(&written).expect("the projection is expressible as TOML");
        let read: ParametersFile = toml::from_str(&text).expect("and parses back");
        assert_eq!(read, written);
    }

    /// **A fitted run writes a file, and it reads back as the same file** — Checkpoint B's claim,
    /// end to end from `RunParameters` through the hand-written writer.
    ///
    /// Not step C4's round trip, which comes back to a `RunParameters` and needs the reader; this
    /// one stops at the shape. What it adds to the writer's own tests is that the text a *run*
    /// produces parses — the every-shape fixture is hand-built, and a projection can emit a row
    /// that fixture has no equivalent of.
    #[test]
    fn a_fitted_run_writes_a_file_that_reads_back() {
        let written = the_projected_file();
        let text = written.to_toml();
        let read: ParametersFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{error}\n\n{text}"));
        assert_eq!(read, written);

        assert!(
            text.contains("[contamination]"),
            "this run had a read group identify a fraction:\n{text}"
        );
        assert_eq!(
            text.lines()
                .filter(|line| line.contains("inbreeding_coefficient"))
                .count(),
            2,
            "one line a sample, which is the form spec §9 prices the per-sample axis in"
        );
    }

    #[test]
    #[should_panic(expected = "libraries and its parameters cover")]
    fn a_read_group_table_with_another_runs_library_count_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        // Four lanes over the run's own two plants: the library counts disagree and the sample
        // counts agree, so only the first guard can fire.
        let four_lanes = ReadGroups::of_lanes(&[
            ("HWI.3", "TS-1", "lib3"),
            ("HWI.4", "TS-1", "lib4"),
            ("HWI.5", AWKWARD_SAMPLE, "lib5"),
            ("HWI.6", AWKWARD_SAMPLE, "lib6"),
        ]);
        let _ = projected(&run, &four_lanes);
    }

    #[test]
    #[should_panic(expected = "samples and the run's parameters cover")]
    fn a_read_group_table_with_another_runs_sample_count_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        // Three lanes as the run has, over three plants rather than two, so only the second
        // guard can fire.
        let three_plants = ReadGroups::of_lanes(&[
            ("HWI.3", "S1", "lib3"),
            ("HWI.4", "S2", "lib4"),
            ("HWI.5", "S3", "lib5"),
        ]);
        let _ = projected(&run, &three_plants);
    }

    /// **The read group whose rate went missing is `Borrowed`, and that is now what makes this
    /// refuse.** Step E2 made a missing rate legal where the calibration is `Defaulted` — the run
    /// with no fit, which has no rates for anybody and writes no counts — so the message this
    /// pins moved with it, from *has a calibration and no fitted rate* to one that names the
    /// warrant. A rate set missing a **fitted or borrowed** read group's entry is still the two
    /// fits coming apart, which is what the refusal is for.
    #[test]
    #[should_panic(expected = "and no rate was offered for it")]
    fn a_rate_set_missing_one_of_the_runs_read_groups_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        let mut one_short = the_runs_fitted_rates();
        one_short.remove(&ReadGroupId(2));
        one_short.insert(
            ReadGroupId(9),
            a_fitted_rate(0.001, Provenance::Borrowed, 640_918),
        );
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &one_short,
            &the_runs_inbreeding(),
            &A_REFERENCE,
            a_census(),
        );
    }

    #[test]
    #[should_panic(expected = "base-quality rates for 5 read groups and the run has 3")]
    fn a_rate_set_over_more_read_groups_than_the_run_has_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        let mut a_wider_fit = the_runs_fitted_rates();
        for group in [3, 4] {
            a_wider_fit.insert(
                ReadGroupId(group),
                a_fitted_rate(0.002, Provenance::FittedHere, 500),
            );
        }
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &a_wider_fit,
            &the_runs_inbreeding(),
            &A_REFERENCE,
            a_census(),
        );
    }

    /// **A rate set covering *some* of the run's read groups is refused at the door, not at the
    /// row.** Step E2 turned that guard from `==` into *empty or complete*, and only the wide case
    /// above pinned it: with the guard relaxed to `<=`, a short set slips past the door and is
    /// caught one frame later by the per-row rule instead, under a message about one read group
    /// rather than about the fit. Both refuse, and which one speaks is the difference between
    /// *this rate set is not this run's* and *read group 2 is missing*.
    #[test]
    #[should_panic(expected = "rates for 2 read groups and the run has 3")]
    fn a_rate_set_covering_some_of_the_runs_read_groups_is_refused_at_the_door() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        let mut one_short = the_runs_fitted_rates();
        one_short.remove(&ReadGroupId(2));
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &one_short,
            &the_runs_inbreeding(),
            &A_REFERENCE,
            a_census(),
        );
    }

    #[test]
    #[should_panic(expected = "two different fits")]
    fn a_rate_set_from_another_fit_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        // The run's own three read groups, another fit's numbers: read group 0's calibration was
        // fitted and this rate says it was supplied, so the count beside the multiplier would be
        // the other fit's.
        let another_fit = BTreeMap::from([
            (ReadGroupId(0), a_fitted_rate(0.09, Provenance::Supplied, 7)),
            (ReadGroupId(1), a_fitted_rate(0.09, Provenance::Borrowed, 8)),
            (
                ReadGroupId(2),
                a_fitted_rate(0.09, Provenance::Defaulted, 9),
            ),
        ]);
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &another_fit,
            &the_runs_inbreeding(),
            &A_REFERENCE,
            a_census(),
        );
    }

    #[test]
    #[should_panic(expected = "the estimate offered here holds")]
    fn an_inbreeding_estimate_that_is_not_the_runs_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        let another_fits_estimates = vec![
            an_inbreeding_estimate(0.31, Provenance::FittedHere, 180_600_412),
            an_inbreeding_estimate(0.17, Provenance::Borrowed, 9_411_027),
        ];
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &the_runs_fitted_rates(),
            &another_fits_estimates,
            &A_REFERENCE,
            a_census(),
        );
    }

    #[test]
    #[should_panic(expected = "inbreeding estimates and the run's parameters hold")]
    fn an_inbreeding_estimate_list_of_another_length_is_refused() {
        let read_groups = a_runs_read_groups();
        let run = a_fitted_run(&read_groups, &the_runs_contamination());
        let one_short = vec![an_inbreeding_estimate(
            0.42,
            Provenance::FittedHere,
            180_600_412,
        )];
        let _ = ParametersFile::of_run(
            &run,
            &read_groups,
            &the_runs_fitted_rates(),
            &one_short,
            &A_REFERENCE,
            a_census(),
        );
    }

    #[test]
    #[should_panic(expected = "is a rate no reader can attach to a library")]
    fn a_substitution_rate_keyed_past_the_runs_read_groups_is_refused() {
        let read_groups = a_runs_read_groups();
        let mut rates = the_runs_substitution_rates();
        rates.insert(
            StratumKey {
                read_group: ReadGroupId(7),
                stratum: ssr_stratum(2, 6),
                ploidy: diploid(),
            },
            a_fitted_rate(0.0009, Provenance::FittedHere, 100),
        );
        let run = a_run_with(
            &read_groups,
            &the_runs_contamination(),
            a_declared_batching(&read_groups),
            the_runs_slippage(),
            rates,
        );
        let _ = projected(&run, &read_groups);
    }

    #[test]
    #[should_panic(expected = "identified no contamination and carries a fraction of")]
    fn a_read_group_measured_at_nothing_may_not_carry_a_fraction() {
        let read_groups = a_runs_read_groups();
        // The estimator returns exactly zero where no marker carries a read; this fixture is what
        // the file could not express if that ever stopped holding.
        let mut contamination = the_runs_contamination();
        contamination.insert(ReadGroupId(1), identified(0.25, 0, 0));
        let run = a_fitted_run(&read_groups, &contamination);
        let _ = projected(&run, &read_groups);
    }

    #[test]
    #[should_panic(expected = "no curve was recorded beside it")]
    fn a_level_that_came_off_a_curve_with_no_curve_recorded_is_refused() {
        level_origin_of(
            LevelProvenance {
                source: LevelSource::Curve,
                curve: None,
                reach: Some(FittedCurveReach::Inside),
                slipped_reads: None,
            },
            Stratum {
                period: 2,
                reference_repeats: 6,
            },
            0,
        );
    }

    #[test]
    #[should_panic(expected = "nothing says whether this stratum sat inside")]
    fn a_level_off_a_curve_that_does_not_say_how_far_it_reached_is_refused() {
        level_origin_of(
            LevelProvenance {
                source: LevelSource::Curve,
                curve: Some(a_slippage_curve()),
                reach: None,
                slipped_reads: None,
            },
            Stratum {
                period: 2,
                reference_repeats: 6,
            },
            0,
        );
    }

    #[test]
    #[should_panic(expected = "the fall-off came off its period's curve")]
    fn a_share_that_came_off_a_curve_with_no_curve_recorded_is_refused() {
        share_smoothing_of(
            ShareProvenance {
                source: ShareSource::Curve,
                curve: None,
                reach: None,
            },
            Stratum {
                period: 2,
                reference_repeats: 6,
            },
            1,
            "the fall-off",
        );
    }
}
