//! **What pass two needs to score a record, built once when the calling pass ends.**
//!
//! Pass one calls every locus and parks each finished record on the spill; the coverage
//! histograms are only complete when that pass has ended, so nothing can be scored until it
//! has. This is what is built in between: one fitted coverage model per sample, the inbreeding
//! coefficient the parameters file carries for each, the σ₀ slice the scorer reads beside them,
//! and the per-pass tables the scorer precomputes from the cohort's size.
//!
//! **A sample can fail to get a model, and which way it failed is kept.** Four ways, and they
//! are not the same news: the pass reached no covered position for that sample at all; every
//! window it did reach was too sparse to trust; the windows had no positive median depth to cut
//! bins at; or the histogram existed and the fit refused it, for any of the six reasons
//! [`CoverageModelError`] names — no usable tiles at all, a single-copy peak in the bottom bin,
//! a mode on the wrong copy-number peak, too much of the sample in the overflow bin, no GC bin
//! dense enough to anchor the curve, or a configuration that is not a configuration. The first
//! three come from [`SampleHistogram`] and the
//! fourth carries the fit's own reason, whichever of the six it was. The run report says how many samples fell to each, because "the
//! coverage model rested on 61 of 63 samples" is a different statement from "on 61 of 63, and
//! the two that dropped out were nearly uncovered".
//!
//! **A sample with no model is absent from every record's score**, not scored as neutral. That
//! is production's behaviour and the copied scorer's: `None` in the observation slice, and the
//! sample simply does not enter the sums.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.1, §3.2.

use thiserror::Error;

use crate::paralog::{
    CoverageFitConfig, CoverageModelError, LocusObservations, ParalogModelParams,
    ParalogScorePrecompute, SampleObservation, SingleCopyCoverageModel, score_locus_for_paralogy,
};
use crate::types::InbreedingF;
use crate::window_coverage::SampleHistogram;

use super::{GenericLocusSample, RepeatTractSample, SpillEntry, SpilledSamples};

/// **A record's sample count disagrees with the run's.**
///
/// Not a property of the data: the sink that filled the spill and this context both claim the
/// run's sample count, so a disagreement is a wiring error between them. It names the record
/// because pass two walks millions of them and a count alone would say nothing about which.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[error(
    "the record at contig {contig} position {position} carries {record} sample(s) and the run \
     has {cohort}; scoring it would silently use a prefix of the cohort"
)]
pub struct CohortSizeMismatch {
    /// The record's contig.
    pub contig: u32,
    /// The record's written position.
    pub position: u64,
    /// How many samples the record carried.
    pub record: usize,
    /// How many the run has.
    pub cohort: usize,
}

/// **The fit configuration is not a usable one**, so no sample could get a model.
///
/// Separate from [`WhyNoCoverageModel`] on purpose: that enum says what a *sample's* data came to,
/// and this says the run was asked to fit with knobs that do not describe a fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
#[error("the paralog coverage fit was configured with {reason}, so no sample could be fitted")]
pub struct CoverageFitConfigRefused {
    /// What the fit said was wrong with it.
    pub reason: &'static str,
}

/// Why a sample has no fitted coverage model, and therefore enters no record's score.
///
/// **Kept per sample rather than counted**, because the run report groups them and because a
/// sample that was never covered and one whose fit was refused are different things to tell an
/// operator about.
/// **Not `Clone`**, because the fit's own error is not: it is a byte-for-byte copy of
/// production's type and gaining a derive would be an edit to it.
/// **`Display` is what spec §3.5's run-report line prints** — "samples whose coverage model was
/// rejected, each with its reason" wants words, and `{:?}` would give an operator a variant name.
#[derive(Debug, PartialEq, Error)]
#[non_exhaustive]
pub enum WhyNoCoverageModel {
    /// The calling pass reached no covered position for this sample, so no window was ever
    /// finalised. Ordinarily this is a sample with no reads in the run's intervals.
    #[error("no window was finalised — the run reached no covered position for this sample")]
    NoWindowFinalised,
    /// Every window this sample finalised held too few covered positions to be trusted, so all
    /// of them were refused and none could train a yardstick.
    #[error("every window this sample finalised held too few covered positions to be trusted")]
    EveryWindowUnderTheFloor,
    /// The windows the depth axis would have been fitted from had no positive median depth, so
    /// there was no bin width to cut the histogram at.
    #[error("the sample's windows had no positive median depth to cut depth bins at")]
    MedianDepthNotPositive,
    /// The histogram existed and the fit refused it. **Carries the fit's own reason**, which
    /// distinguishes an essentially uncovered sample from one whose depth ran off the top of
    /// the histogram — opposite problems that a bare count would merge.
    #[error("the coverage fit refused this sample's histogram: {0}")]
    TheFitRefusedTheHistogram(#[source] CoverageModelError),
}

/// **Everything the scorer reads that does not change from record to record.**
///
/// Built once, between the calling pass and the scoring pass. Holds one entry per sample of the
/// run, in the run's sample order, so an index means the same sample here, in the spill's rows,
/// and in the parameters file.
#[derive(Debug)]
pub struct ParalogScoringContext {
    /// One fitted model per sample, `None` where there is none.
    coverage_models: Vec<Option<SingleCopyCoverageModel>>,
    /// Why, for each sample that has no model. `None` where it has one.
    why_no_model: Vec<Option<WhyNoCoverageModel>>,
    /// The parameters file's inbreeding coefficient per sample.
    inbreeding: Vec<f64>,
    /// σ₀ per sample — the one-copy relative-depth SD — with `NaN` where the sample has no
    /// model. **The scorer takes this as a slice parallel to the observations** and returns a
    /// neutral score if its length disagrees with the cohort's, so it is built here from the
    /// same iteration that builds the models and can never be a different length.
    single_copy_depth_sd: Vec<f64>,
    /// The per-pass tables, built once from the model params and the cohort's coefficients.
    precompute: ParalogScorePrecompute,
}

impl ParalogScoringContext {
    /// Fit one coverage model per sample and build the tables the scorer reuses.
    ///
    /// `histograms` and `inbreeding` are both in the run's sample order and must be the same
    /// length; that length is the cohort size every scored record must also have.
    ///
    /// **The histograms are taken by value on purpose.** Each sample's bins are freed as its fit
    /// finishes rather than at the end of the run, which is the difference between holding one
    /// histogram and holding three thousand. `fit` itself only borrows, so this looks like
    /// unnecessary ownership until you count what a cohort's worth of them costs.
    ///
    /// # Panics
    ///
    /// If the two slices disagree in length. They come from the same run, indexed by the same
    /// sample order, so a mismatch is a wiring error rather than an input the caller can meet —
    /// and the alternative is a scorer that silently drops the tail of the cohort.
    ///
    /// # Errors
    ///
    /// If the fit configuration is not a usable one. **That is an operator's mistake and not a
    /// verdict about the cohort**, and it has to be told apart from one: the copied fit
    /// re-validates its configuration on every call, so folding the refusal into the per-sample
    /// outcome would make one wrong knob report that the coverage model rested on 0 of 63
    /// samples, with 63 identical reasons — a configuration error dressed as a statement about
    /// every sample's coverage.
    pub fn new(
        histograms: Vec<SampleHistogram>,
        inbreeding: &[InbreedingF],
        params: &ParalogModelParams,
        fit: &CoverageFitConfig,
    ) -> Result<Self, CoverageFitConfigRefused> {
        assert_eq!(
            histograms.len(),
            inbreeding.len(),
            "the run has one histogram and one inbreeding coefficient per sample, and these \
             disagree; a scorer built from them would drop the tail of the cohort"
        );

        // **Unwrapped once, here.** The parameters file carries these validated to `[0, 1)`; the
        // copied scorer takes plain `f64`s, because it is production's code unchanged. This is
        // the one boundary between the two, so it is the one place the newtype comes off.
        let inbreeding: Vec<f64> = inbreeding.iter().map(|f| f.get()).collect();

        let mut coverage_models = Vec::with_capacity(histograms.len());
        let mut why_no_model = Vec::with_capacity(histograms.len());
        let mut single_copy_depth_sd = Vec::with_capacity(histograms.len());

        for histogram in histograms {
            let outcome = match histogram {
                SampleHistogram::Fitted(histogram) => SingleCopyCoverageModel::fit(&histogram, fit)
                    .map_err(WhyNoCoverageModel::TheFitRefusedTheHistogram),
                SampleHistogram::NoWindowFinalised => Err(WhyNoCoverageModel::NoWindowFinalised),
                SampleHistogram::EveryWindowUnderTheFloor => {
                    Err(WhyNoCoverageModel::EveryWindowUnderTheFloor)
                }
                SampleHistogram::MedianDepthNotPositive => {
                    Err(WhyNoCoverageModel::MedianDepthNotPositive)
                }
            };

            // A refused configuration is the same refusal for every sample, so it is reported
            // once and construction fails, rather than being counted 63 times as data.
            if let Err(WhyNoCoverageModel::TheFitRefusedTheHistogram(
                CoverageModelError::InvalidConfig { reason },
            )) = outcome
            {
                return Err(CoverageFitConfigRefused { reason });
            }

            match outcome {
                Ok(model) => {
                    // **The three slices are filled in one step each**, so they cannot come out
                    // of step with each other or with the cohort.
                    single_copy_depth_sd.push(model.single_copy_depth_sd());
                    coverage_models.push(Some(model));
                    why_no_model.push(None);
                }
                Err(why) => {
                    single_copy_depth_sd.push(f64::NAN);
                    coverage_models.push(None);
                    why_no_model.push(Some(why));
                }
            }
        }

        Ok(Self {
            precompute: ParalogScorePrecompute::new(params, &inbreeding),
            coverage_models,
            why_no_model,
            inbreeding,
            single_copy_depth_sd,
        })
    }

    /// **What one sample showed at one record, as the scorer reads it** — or `None` where the
    /// sample says nothing about this record.
    ///
    /// Three ways to say nothing, and they are all ordinary rather than errors: the sample has
    /// no coverage model; it has no usable window at this locus, which arrives as the `NaN`
    /// pair; or — **at a generic locus only** — no read reached it.
    ///
    /// **The zero-read skip is production's, and it must not reach a repeat tract** (spec §6
    /// trap 1). Production drops a sample at zero total reads because a sample with no reads
    /// says nothing about allele balance. A tract's rows carry no read counts at all — the
    /// filter declines to read a split that slippage has smeared, which is an abstention and
    /// not an observation of zero — so there is no total for the skip to test, and the type is
    /// what says so rather than a condition someone has to remember.
    #[must_use]
    fn observation_of(&self, entry: &SpillEntry, sample: usize) -> Option<SampleObservation> {
        let model = self.coverage_models.get(sample)?.as_ref()?;

        // **One match yields the window and the counts together**, and both row types are
        // destructured exhaustively — the sibling codec does the same, for the same reason: a
        // field added to either row is a field this, the only reader of them, would otherwise
        // drop in silence. Spec §8's deferred tract-aware allele term is exactly such a field.
        let (window, alt_reads, total_reads) = match &entry.samples {
            SpilledSamples::GenericLocus(samples) => {
                let GenericLocusSample {
                    window,
                    ref_reads,
                    alt_reads,
                } = samples.get(sample)?;
                // A pair that cannot be summed is a corrupt spill row — two counts whose
                // total overflows `u32`, which no run produces and no codec refuses. It is
                // treated as an absence, which is what an unreadable row is: the sample says
                // nothing about this record. **Nothing counts it**, and the only caller is
                // `score`, which cannot tell it from a sample that was simply not covered.
                let total = ref_reads.checked_add(*alt_reads)?;
                if total == 0 {
                    return None;
                }
                (*window, *alt_reads, total)
            }
            // No read counts to read, and no skip. The allele term is
            // `alt·ln(vaf) + (total − alt)·ln(1 − vaf)`, which is zero at zero reads under every
            // genotype and every carrier count, so it cancels and the ratio rests on coverage
            // alone (spec §3.2).
            SpilledSamples::RepeatTract(samples) => {
                let RepeatTractSample { window } = samples.get(sample)?;
                (*window, 0, 0)
            }
        };

        // **Stricter than [`WindowCoverage::is_absent`], deliberately.** That predicate asks
        // whether either field is `NaN`; this rejects any non-finite value, and the difference is
        // `±∞`. An infinite depth would divide to `+∞`, be winsorised by the scorer to
        // `max_relative_copy_number` and read as a *confident four-copy paralog* — the filter's
        // strongest possible coverage signal, from a measurement that never happened. The `NaN`
        // half is load-bearing too, and against a panic rather than a wrong answer: a `NaN` GC
        // fraction reaches `gc_multiplier`, where both its range comparisons are false, the
        // floor saturates to zero, and the interpolation indexes one past the end of the curve.
        if !window.gc_fraction.is_finite() || !window.mean_depth.is_finite() {
            return None;
        }
        let relative_copy_number =
            model.relative_copy_number(f64::from(window.gc_fraction), f64::from(window.mean_depth));

        Some(SampleObservation {
            relative_copy_number,
            alt_reads,
            total_reads,
            inbreeding_coefficient: *self.inbreeding.get(sample)?,
        })
    }

    /// **Every sample's observation for one record**, into a buffer the caller keeps.
    ///
    /// This is what pass two calls; [`Self::observation_of`] is one sample of it. The buffer is
    /// cleared and refilled, so one allocation serves the whole run — the project's
    /// load / use / clear / reload shape.
    ///
    /// **The record's width is checked once here, and that is the point of the method.** Per
    /// sample there is nothing to check against: an index past the record's own rows is
    /// indistinguishable from a sample that was simply not covered, so a record narrower than
    /// the cohort would score on a prefix and a wider one would be truncated —
    /// [`score_locus_for_paralogy`](crate::paralog::score_locus_for_paralogy) answers a
    /// length mismatch with a *neutral score* rather than an error, so neither would fail
    /// loudly. Both are wiring errors between the sink that filled the spill and this context,
    /// and both produce a plausible, quietly weaker score on some records — the defect shape
    /// that survives a whole run because nothing about the output looks wrong.
    ///
    /// # Errors
    ///
    /// If the record carries a different number of samples than the run has.
    fn observations_of(
        &self,
        entry: &SpillEntry,
        out: &mut Vec<Option<SampleObservation>>,
    ) -> Result<(), CohortSizeMismatch> {
        if entry.samples.len() != self.sample_count() {
            return Err(CohortSizeMismatch {
                contig: entry.contig.get(),
                position: entry.position.get(),
                record: entry.samples.len(),
                cohort: self.sample_count(),
            });
        }
        out.clear();
        out.extend((0..self.sample_count()).map(|sample| self.observation_of(entry, sample)));
        Ok(())
    }

    /// **One parked record's likelihood ratio — or `NaN`, which means it was not scored.**
    ///
    /// This is the only way in: the σ₀ slice and the precomputed tables are the scorer's other
    /// two arguments, and handing them out separately would let a caller pair them with
    /// observations they were not built for. The scorer answers that with a *neutral score*
    /// rather than an error (spec §6 trap 3), so the mistake would show up as a run that
    /// quietly flags nothing. Here the three arguments come from one place and cannot be
    /// mismatched.
    ///
    /// **`NaN` where no sample was usable, and never `0.0`.** The scorer's neutral verdict is a
    /// ratio of `0.0` — the value it returns when the two stories are exactly balanced, and also
    /// the value it returns when it had nothing to weigh at all. Those are opposite states and
    /// they must not share a number: a `0.0` folded into the run's histogram is one more record
    /// saying "not a duplication", and a whole run of unscorable records would fit a paralog rate
    /// from evidence that does not exist (spec §6 trap 4). `samples_used` is the scorer's own
    /// count of what it weighed, so this reads it rather than re-deriving the condition.
    ///
    /// `observations` is a scratch buffer the caller keeps across records — cleared and refilled
    /// here, so one allocation serves the run.
    ///
    /// # Errors
    ///
    /// If the record carries a different number of samples than the run has.
    pub fn score(
        &self,
        entry: &SpillEntry,
        observations: &mut Vec<Option<SampleObservation>>,
    ) -> Result<f64, CohortSizeMismatch> {
        self.observations_of(entry, observations)?;
        let score = score_locus_for_paralogy(
            &LocusObservations {
                samples: observations,
            },
            &self.single_copy_depth_sd,
            &self.precompute,
        );
        Ok(if score.samples_used == 0 {
            f64::NAN
        } else {
            score.paralog_log_likelihood_ratio
        })
    }

    /// How many samples the run has — the cohort size every scored record must match.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.coverage_models.len()
    }

    /// How many samples have a fitted coverage model.
    #[must_use]
    pub fn how_many_samples_have_a_coverage_model(&self) -> usize {
        self.coverage_models.iter().flatten().count()
    }

    /// Why each sample without a model has none, in the run's sample order — `None` where it
    /// has one. What the run report groups (spec §3.5).
    #[must_use]
    pub fn why_no_model(&self) -> &[Option<WhyNoCoverageModel>] {
        &self.why_no_model
    }

    /// **What each sample's coverage fit came to**, in the run's sample order — `None` where the
    /// fit was refused.
    ///
    /// **Why this is here rather than being inferred later.** The two numbers below are what the
    /// score divides by and measures against, and a run that prints neither can only be explained
    /// by inverting its own output — which is how D1's first report explained it, and the
    /// explanation was wrong twice over: the depth a record was compared against was taken to be
    /// the median depth of the records the run *wrote*, which is a variant-site subset and not
    /// one copy of anything. Plan step D1 asks for "the fit's outcome"; this is it.
    ///
    /// It is deliberately a small owned summary rather than the models themselves. A caller
    /// holding a [`SingleCopyCoverageModel`] could ask it for a copy number at a GC of its
    /// choosing and index the answer against the wrong sample, which is the shape of spec §6's
    /// trap 3; a caller holding two numbers per sample cannot.
    #[must_use]
    pub fn what_each_fit_came_to(&self) -> Vec<Option<WhatTheFitCameTo>> {
        self.coverage_models
            .iter()
            .map(|model| {
                model.as_ref().map(|model| WhatTheFitCameTo {
                    one_copy_depth: model.single_copy_scale(),
                    single_copy_depth_sd: model.single_copy_depth_sd(),
                })
            })
            .collect()
    }
}

/// **One sample's fitted coverage model, in the two numbers a reader needs** — spec §3.1's fit,
/// as the run report says it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhatTheFitCameTo {
    /// **Reads a one-copy window carries in this sample**, before the GC multiplier — the
    /// depth-distribution mode the fit anchors on
    /// ([`single_copy_scale`](SingleCopyCoverageModel::single_copy_scale)). This is what a
    /// record's window depth is divided by to give its copy number, so it is the number a
    /// reader needs to turn "19 reads" into "three copies", and it is **not** the median depth
    /// of the records the run wrote.
    pub one_copy_depth: f64,
    /// **σ₀ — how much a one-copy window's relative depth scatters in this sample.** The score's
    /// coverage half measures a record's departure from one copy in these units, so it is what
    /// says whether a window at twice one copy's depth is remarkable or ordinary.
    pub single_copy_depth_sd: f64,
}

#[cfg(test)]
mod tests;
