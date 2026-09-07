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
//! bins at; or the histogram existed and the fit refused it — an unresolvable single-copy peak,
//! a mode on the wrong copy-number peak, too much of the sample in the overflow bin. The first
//! three come from [`SampleHistogram`](crate::ng::window_coverage::SampleHistogram) and the
//! fourth from the fit. The run report says how many samples fell to each, because "the
//! coverage model rested on 61 of 63 samples" is a different statement from "on 61 of 63, and
//! the two that dropped out were nearly uncovered".
//!
//! **A sample with no model is absent from every record's score**, not scored as neutral. That
//! is production's behaviour and the copied scorer's: `None` in the observation slice, and the
//! sample simply does not enter the sums.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md` §3.1, §3.2.

use crate::ng::paralog::{
    CoverageFitConfig, CoverageModelError, ParalogModelParams, ParalogScorePrecompute,
    SampleObservation, SingleCopyCoverageModel,
};
use crate::ng::window_coverage::SampleHistogram;

use super::{SpillEntry, SpilledSamples};

/// Why a sample has no fitted coverage model, and therefore enters no record's score.
///
/// **Kept per sample rather than counted**, because the run report groups them and because a
/// sample that was never covered and one whose fit was refused are different things to tell an
/// operator about.
/// **Not `Clone`**, because the fit's own error is not: it is a byte-for-byte copy of
/// production's type and gaining a derive would be an edit to it.
#[derive(Debug, PartialEq)]
#[non_exhaustive]
pub enum WhyNoCoverageModel {
    /// The calling pass reached no covered position for this sample, so no window was ever
    /// finalised. Ordinarily this is a sample with no reads in the run's intervals.
    NoWindowFinalised,
    /// Every window this sample finalised held too few covered positions to be trusted, so all
    /// of them were refused and none could train a yardstick.
    EveryWindowUnderTheFloor,
    /// The windows the depth axis would have been fitted from had no positive median depth, so
    /// there was no bin width to cut the histogram at.
    MedianDepthNotPositive,
    /// The histogram existed and the fit refused it. **Carries the fit's own reason**, which
    /// distinguishes an essentially uncovered sample from one whose depth ran off the top of
    /// the histogram — opposite problems that a bare count would merge.
    TheFitRefusedTheHistogram(CoverageModelError),
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
    /// # Panics
    ///
    /// If the two slices disagree in length. They come from the same run, indexed by the same
    /// sample order, so a mismatch is a wiring error rather than an input the caller can meet —
    /// and the alternative is a scorer that silently drops the tail of the cohort.
    #[must_use]
    pub fn new(
        histograms: Vec<SampleHistogram>,
        inbreeding: &[f64],
        params: &ParalogModelParams,
        fit: &CoverageFitConfig,
    ) -> Self {
        assert_eq!(
            histograms.len(),
            inbreeding.len(),
            "the run has one histogram and one inbreeding coefficient per sample, and these \
             disagree; a scorer built from them would drop the tail of the cohort"
        );

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

        Self {
            coverage_models,
            why_no_model,
            inbreeding: inbreeding.to_vec(),
            single_copy_depth_sd,
            precompute: ParalogScorePrecompute::new(params, inbreeding),
        }
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
    pub fn observation_of(&self, entry: &SpillEntry, sample: usize) -> Option<SampleObservation> {
        let model = self.coverage_models.get(sample)?.as_ref()?;
        let window = entry.samples.window(sample)?;

        // The `NaN` pair is *this sample has no usable window here*, and it is compared as a
        // number rather than by bits on purpose: any non-finite value is unusable, however it
        // arrived.
        if !window.gc_fraction.is_finite() || !window.mean_depth.is_finite() {
            return None;
        }
        let relative_copy_number =
            model.relative_copy_number(f64::from(window.gc_fraction), f64::from(window.mean_depth));

        let (alt_reads, total_reads) = match &entry.samples {
            SpilledSamples::GenericLocus(samples) => {
                let sample = samples.get(sample)?;
                let total = sample.ref_reads.checked_add(sample.alt_reads)?;
                if total == 0 {
                    return None;
                }
                (sample.alt_reads, total)
            }
            // No read counts to read, and no skip. The allele term is
            // `alt·ln(vaf) + (total − alt)·ln(1 − vaf)`, which is zero at zero reads under every
            // genotype and every carrier count, so it cancels and the ratio rests on coverage
            // alone (spec §3.2).
            SpilledSamples::RepeatTract(_) => (0, 0),
        };

        Some(SampleObservation {
            relative_copy_number,
            alt_reads,
            total_reads,
            inbreeding_coefficient: *self.inbreeding.get(sample)?,
        })
    }

    /// How many samples the run has — the cohort size every scored record must match.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.coverage_models.len()
    }

    /// How many samples have a fitted coverage model.
    #[must_use]
    pub fn samples_with_a_coverage_model(&self) -> usize {
        self.coverage_models.iter().flatten().count()
    }

    /// Why each sample without a model has none, in the run's sample order — `None` where it
    /// has one. What the run report groups (spec §3.5).
    #[must_use]
    pub fn why_no_model(&self) -> &[Option<WhyNoCoverageModel>] {
        &self.why_no_model
    }

    /// σ₀ per sample, the slice the scorer takes beside the observations.
    #[must_use]
    pub fn single_copy_depth_sd(&self) -> &[f64] {
        &self.single_copy_depth_sd
    }

    /// The per-pass tables the scorer reuses across every record.
    #[must_use]
    pub fn precompute(&self) -> &ParalogScorePrecompute {
        &self.precompute
    }
}

#[cfg(test)]
mod tests;
