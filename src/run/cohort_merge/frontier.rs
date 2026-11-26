//! **Spec §6.1's frontier rule**: when the regions built in parallel are handed over in region
//! order, a locus that begins on ground an earlier locus already owns is dropped and counted.
//!
//! The streaming driver ([`super::stream`]) and the round driver it is checked against both hand
//! their regions over through [`Frontier`].

use crate::types::{GenomePosition, GenomeRegion};

/// What one region's worker produced: what the closure made of each locus, in genome order,
/// and the spans the width bound refused.
pub(super) struct CalledRegion<T> {
    /// Each locus's span beside what the closure made of it.
    pub(super) calls: Vec<(GenomeRegion, T)>,
    /// The spans this region's builder refused, in genome order.
    pub(super) failed: Vec<GenomeRegion>,
}

/// **Spec §6.1's frontier rule**, for the drivers that build regions in parallel: each region's
/// calls and refused spans, handed over in region order, are interleaved by first base, and a
/// span that begins on ground an earlier one owns is dropped and counted.
#[derive(Default)]
pub(super) struct Frontier {
    /// The last base claimed so far.
    owned_through: Option<GenomePosition>,
    /// How many spans were dropped.
    pub(super) displaced: u64,
}

impl Frontier {
    /// Hand one region's outcome over: each call that claims its ground to `take`, each refused
    /// span that does to `refused`.
    pub(super) fn hand_over<T>(
        &mut self,
        region: CalledRegion<T>,
        take: &mut impl FnMut(T),
        refused: &mut Vec<GenomeRegion>,
    ) {
        let CalledRegion { calls, failed } = region;
        let mut calls = calls.into_iter().peekable();
        for span in failed {
            let span_first = first_base_of(span);
            while let Some((locus, made)) =
                calls.next_if(|(locus, _)| first_base_of(*locus) <= span_first)
            {
                if self.claim(locus) {
                    take(made);
                }
            }
            if self.claim(span) {
                refused.push(span);
            }
        }
        for (locus, made) in calls {
            if self.claim(locus) {
                take(made);
            }
        }
    }

    /// Claim `span`'s ground, or refuse it to whatever got there first.
    fn claim(&mut self, span: GenomeRegion) -> bool {
        if self
            .owned_through
            .is_some_and(|owned| first_base_of(span) <= owned)
        {
            self.displaced = self.displaced.saturating_add(1);
            return false;
        }
        self.owned_through = Some(last_base_of(span));
        true
    }
}

/// A locus span's first base, genome-wide — what decides which of two loci is the earlier
/// owner.
///
/// **The ends are ordered rather than trusted**, as [`ObservationCache::with_observations`] and
/// [`ObservationCache::cover`] order theirs and for the same reason: `GenomeRegion` has public
/// fields and no constructor putting them in order, and an inverted span would otherwise walk
/// the frontier backwards — past its own first base, releasing the locus behind it.
///
/// [`ObservationCache::with_observations`]: super::observation_cache::ObservationCache::with_observations
/// [`ObservationCache::cover`]: super::observation_cache::ObservationCache::cover
pub(super) fn first_base_of(locus_span: GenomeRegion) -> GenomePosition {
    GenomePosition {
        contig: locus_span.contig,
        position: locus_span.start.min(locus_span.end),
    }
}

/// A locus span's last base, genome-wide — where the frontier stands once the span is claimed.
/// Ordered for [`first_base_of`]'s reason, and never before it.
pub(super) fn last_base_of(locus_span: GenomeRegion) -> GenomePosition {
    GenomePosition {
        contig: locus_span.contig,
        position: locus_span.end.max(locus_span.start),
    }
}
