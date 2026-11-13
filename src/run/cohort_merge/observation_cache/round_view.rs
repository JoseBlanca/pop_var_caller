//! **A round's records, sealed so they can be read while the cache draws on** — what lets the
//! next round's ground be read while this round's regions are still being built
//! (`spec/cohort_merge_overlapping_rounds.md`).
//!
//! A builder reads the cache through a shared borrow, and a cover draws it forward through an
//! exclusive one, so the two could not overlap (`cohort_merge.md` §6.2). [`ObservationCache::seal_round`]
//! ends that: once a round's ground is covered, the records each sample holds are **moved** into a
//! [`RoundView`] the builders read instead, and only the records the next round also needs — those
//! reaching past this round's end — are copied back into the cache. A source that keeps its
//! evidence seals the stretch of its arena the round's bodies lie in, so the view builds them
//! without the source ([`ObservationSource::seal`]).
//!
//! **What it costs is the copy back**: the records and windows reaching past the round's end,
//! which is the look-ahead half a window long and whatever chains across the edge — small beside
//! the round. The view itself is the round's records as the cache held them, not a copy.

use std::sync::Arc;

use super::{
    BuildsKeptEvidence, LocusSummary, ObservationCache, ObservationSource, WindowedCohort,
    first_reaching_summary,
};
use crate::locus_generation::SampleLocusObservations;
use crate::run::explain::DroppedLoci;
use crate::types::{GenomePosition, GenomeRegion};
use crate::window_coverage::WindowCoverage;

/// One sample's part of a sealed round.
struct SampleView<E> {
    summaries: Vec<LocusSummary>,
    bodies: Vec<core::ops::Range<usize>>,
    observations: Vec<SampleLocusObservations>,
    windows: Vec<(GenomePosition, WindowCoverage)>,
    evidence: Option<Arc<dyn BuildsKeptEvidence<E>>>,
}

/// **Every sample's records over one round's ground, read-only**, as the cache held them when the
/// round's cover ended ([`ObservationCache::seal_round`]). It answers a builder the questions the
/// cache answers — [`with_observations`](Self::with_observations),
/// [`window_starts_into`](Self::window_starts_into), [`build_at`](Self::build_at) — the same way.
pub struct RoundView<E> {
    samples: Vec<SampleView<E>>,
    keeps_evidence: bool,
    covered_to: Option<GenomePosition>,
    dropped_loci: Option<Arc<DroppedLoci>>,
}

impl<S, E> ObservationCache<S>
where
    S: ObservationSource<Error = E> + Send,
    E: 'static,
{
    /// **Seal what every sample holds into a [`RoundView`]**, keeping in the cache only what
    /// reaches `next_from` or past it — the first base the next round can start at.
    ///
    /// Called after a round's cover. What the view holds is what the round's builders would have
    /// read from the cache; what the cache keeps is what the next round's eviction would have
    /// kept, so its cover and its builders see what they saw before.
    pub(in crate::run::cohort_merge) fn seal_round(
        &mut self,
        next_from: GenomePosition,
    ) -> RoundView<E> {
        use rayon::prelude::*;

        let keeps_evidence = self.keeps_evidence;
        let samples = self
            .samples
            .par_iter_mut()
            .map(|sample| {
                let tail = first_reaching_summary(&sample.held_summaries, next_from);
                let summaries = std::mem::take(&mut sample.held_summaries);
                sample.held_summaries = summaries[tail..].to_vec();
                let bodies = std::mem::take(&mut sample.held_bodies);
                sample.held_bodies = bodies[tail.min(bodies.len())..].to_vec();
                let observations = std::mem::take(&mut sample.held_observations);
                sample.held_observations = observations[tail.min(observations.len())..].to_vec();
                let finalised = &mut sample.window_coverage.finalised;
                let windows = std::mem::take(finalised);
                let windows_tail = windows.partition_point(|(centre, _)| *centre < next_from);
                *finalised = windows[windows_tail..].to_vec();
                let evidence = if sample.keeps_evidence {
                    sample
                        .source
                        .seal(bodies.first().map_or(usize::MAX, |body| body.start))
                } else {
                    None
                };
                SampleView {
                    summaries,
                    bodies,
                    observations,
                    windows,
                    evidence,
                }
            })
            .collect();
        RoundView {
            samples,
            keeps_evidence,
            covered_to: self.covered_to,
            dropped_loci: self.dropped_loci.clone(),
        }
    }
}

impl<E> RoundView<E> {
    /// [`ObservationCache::with_observations`], over the sealed round.
    ///
    /// # Panics
    ///
    /// If `span` reaches past the ground the round's cover drew, as the cache's does.
    pub fn with_observations<R>(
        &self,
        span: GenomeRegion,
        f: impl FnOnce(&WindowedCohort<'_>) -> R,
    ) -> R {
        let left_edge = GenomePosition {
            contig: span.contig,
            position: span.start.min(span.end),
        };
        let last_base = GenomePosition {
            contig: span.contig,
            position: span.start.max(span.end),
        };
        assert!(
            self.covered_to
                .is_some_and(|covered_to| last_base <= covered_to),
            "the window was asked for {span} but the round's cover reached only {:?}",
            self.covered_to,
        );
        let mut observations: Vec<&[SampleLocusObservations]> =
            Vec::with_capacity(self.samples.len());
        let mut summaries: Vec<&[LocusSummary]> = Vec::with_capacity(self.samples.len());
        let mut windows: Vec<&[(GenomePosition, WindowCoverage)]> =
            Vec::with_capacity(self.samples.len());
        for sample in &self.samples {
            let from = first_reaching_summary(&sample.summaries, left_edge);
            summaries.push(&sample.summaries[from..]);
            windows.push(&sample.windows);
            if !self.keeps_evidence {
                observations.push(&sample.observations[from..]);
            }
        }
        let records = (!self.keeps_evidence).then_some(&observations[..]);
        f(&WindowedCohort {
            observations: records,
            summaries: Some(&summaries),
            finalised_windows: Some(&windows),
        })
    }

    /// [`ObservationCache::window_starts_into`], over the sealed round.
    pub(in crate::run::cohort_merge) fn window_starts_into(
        &self,
        span: GenomeRegion,
        into: &mut Vec<usize>,
    ) {
        let left_edge = GenomePosition {
            contig: span.contig,
            position: span.start.min(span.end),
        };
        into.clear();
        into.extend(
            self.samples
                .iter()
                .map(|sample| first_reaching_summary(&sample.summaries, left_edge)),
        );
    }

    /// [`ObservationCache::build_at`], over the sealed round.
    ///
    /// # Errors
    ///
    /// Whatever the sample's sealed evidence refuses.
    pub(in crate::run::cohort_merge) fn build_at(
        &self,
        sample: usize,
        index: usize,
    ) -> Result<SampleLocusObservations, E> {
        let sample = &self.samples[sample];
        sample
            .evidence
            .as_ref()
            .expect("a round that builds from kept evidence sealed its sample's arena")
            .build(sample.bodies[index].clone())
    }

    /// The `--explain-loci` log, shared with the cache.
    pub fn dropped_loci(&self) -> Option<&DroppedLoci> {
        self.dropped_loci.as_deref()
    }
}
