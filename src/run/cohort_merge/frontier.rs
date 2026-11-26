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

#[cfg(test)]
mod tests {
    use super::super::fixtures::{region, region_on};
    use super::*;

    /// What a [`Frontier`] let through when the regions were handed over in the order given:
    /// the loci it handed to `take`, the refused spans it kept, and how many spans it dropped.
    #[derive(Debug, PartialEq, Eq)]
    struct HandedOver {
        taken: Vec<GenomeRegion>,
        refused: Vec<GenomeRegion>,
        displaced: u64,
    }

    /// Hand each region over in turn — its loci, then its refused spans — through one
    /// [`Frontier`], with each locus's span standing for what the closure made of it.
    ///
    /// **Within one region both lists must be ascending and must not overlap each other**,
    /// because that is what one walk over one region produces. Every overlap in these fixtures
    /// is between regions, which is the only place a builder could put one.
    fn hand_over(regions: &[(&[GenomeRegion], &[GenomeRegion])]) -> HandedOver {
        let mut frontier = Frontier::default();
        let mut taken = Vec::new();
        let mut refused = Vec::new();
        for (loci, failed) in regions {
            let region = CalledRegion {
                calls: loci.iter().map(|span| (*span, *span)).collect(),
                failed: failed.to_vec(),
            };
            frontier.hand_over(region, &mut |made| taken.push(made), &mut refused);
        }
        HandedOver {
            taken,
            refused,
            displaced: frontier.displaced,
        }
    }

    fn outcome(taken: &[GenomeRegion], refused: &[GenomeRegion], displaced: u64) -> HandedOver {
        HandedOver {
            taken: taken.to_vec(),
            refused: refused.to_vec(),
            displaced,
        }
    }

    // Every fixture below is fabricated: under `build_region`'s input contract no builder can
    // produce two loci that overlap, so nothing in a real merge reaches the drop
    // (`super::super::serial::tests::refuse_displaced_loci` checks that through real builders).
    // These state what the rule does if it is ever reached.

    /// Within one region the builder's order is the genome's, and the frontier keeps it.
    #[test]
    fn the_loci_of_one_region_come_out_in_the_order_the_builder_gave_them() {
        let loci = [region(1, 3), region(7, 7), region(11, 18)];
        assert_eq!(hand_over(&[(&loci, &[])]), outcome(&loci, &[], 0));
    }

    /// The rule itself: of two loci that overlap, the one that starts earlier stands.
    #[test]
    fn a_locus_starting_on_ground_an_earlier_locus_owns_is_dropped() {
        assert_eq!(
            hand_over(&[(&[region(10, 40)], &[]), (&[region(25, 60)], &[])]),
            outcome(&[region(10, 40)], &[], 1),
        );
    }

    /// **The test is the first base, not the overlap.** A locus that opens on the very last
    /// base an earlier one owns is dropped; one that opens on the next base is kept, even
    /// though the two are adjacent.
    #[test]
    fn the_boundary_is_the_earlier_locus_last_base() {
        assert_eq!(
            hand_over(&[(&[region(10, 40)], &[]), (&[region(40, 55)], &[])]),
            outcome(&[region(10, 40)], &[], 1),
            "a locus opening on base 40 starts inside ground owned through base 40",
        );
        assert_eq!(
            hand_over(&[(&[region(10, 40)], &[]), (&[region(41, 55)], &[])]),
            outcome(&[region(10, 40), region(41, 55)], &[], 0),
        );
    }

    /// **A refused span displaces exactly as an emitted locus does** (spec §3.2, §6.1). This is
    /// the case the rule exists for: a wide deletion refused by the width bound owns ground that
    /// a later builder, which never saw what opened there, would build from a partial picture.
    #[test]
    fn a_refused_span_displaces_what_starts_inside_it() {
        assert_eq!(
            hand_over(&[
                (&[], &[region(10, 200)]),
                (&[region(120, 130)], &[]),
                (&[region(210, 214)], &[]),
            ]),
            outcome(&[region(210, 214)], &[region(10, 200)], 1),
            "the locus built inside the refused span survived with nothing to displace it",
        );
    }

    /// And the other way round: a refused span starting inside an emitted locus's ground is
    /// dropped, and so is not reported as refused.
    #[test]
    fn a_refused_span_displaced_by_an_earlier_locus_is_not_reported() {
        assert_eq!(
            hand_over(&[(&[region(10, 90)], &[]), (&[], &[region(50, 300)])]),
            outcome(&[region(10, 90)], &[], 1),
            "a refusal the run never owned must not reach the run summary",
        );
    }

    /// **The two lists are resolved as one sequence, not one after the other.** The first
    /// region's refusal at 10-200 and its locus at 220-224 interleave with the second region's
    /// locus at 100-110; resolving all the loci first and then all the refusals would let the
    /// refusal be judged against a frontier that had already run to 224.
    #[test]
    fn loci_and_refused_spans_are_resolved_in_one_genome_order() {
        assert_eq!(
            hand_over(&[
                (&[region(220, 224)], &[region(10, 200)]),
                (&[region(100, 110)], &[]),
            ]),
            outcome(&[region(220, 224)], &[region(10, 200)], 1),
        );
    }

    /// Inside one region too: a refused span between two loci is claimed between them, so the
    /// frontier is at the refusal's end, not the second locus's, when the refusal is judged.
    #[test]
    fn a_refused_span_between_two_loci_of_one_region_keeps_its_place() {
        assert_eq!(
            hand_over(&[(&[region(1, 3), region(61, 64)], &[region(10, 60)])]),
            outcome(&[region(1, 3), region(61, 64)], &[region(10, 60)], 0),
        );
    }

    /// One wide locus can displace several later ones, and the frontier does not retreat to
    /// the last one kept.
    #[test]
    fn one_wide_locus_displaces_every_locus_that_opens_inside_it() {
        assert_eq!(
            hand_over(&[
                (&[region(10, 300)], &[]),
                (&[region(40, 44)], &[]),
                (&[region(120, 124)], &[]),
                (&[region(280, 284)], &[]),
                (&[region(301, 305)], &[]),
            ]),
            outcome(&[region(10, 300), region(301, 305)], &[], 3),
        );
    }

    /// **The frontier is the latest ground owned, not the first.** A locus kept after another
    /// owns what it covers just as the first did, so the one that opens inside *it* is
    /// displaced — which a frontier that stopped at the first locus would let through.
    #[test]
    fn a_locus_kept_after_another_owns_its_own_ground_too() {
        assert_eq!(
            hand_over(&[
                (&[region(10, 20)], &[]),
                (&[region(30, 100)], &[]),
                (&[region(50, 60)], &[]),
            ]),
            outcome(&[region(10, 20), region(30, 100)], &[], 1),
        );
    }

    /// **A displaced locus owns nothing** — the frontier stays on the standing owner's last
    /// base, not the losing locus's. The second region's locus reaches to 100 and loses to the
    /// first's, so the ground from 41 to 100 is nobody's and the third region's locus there
    /// stands.
    #[test]
    fn a_displaced_locus_does_not_own_the_ground_it_lost() {
        assert_eq!(
            hand_over(&[
                (&[region(10, 40)], &[]),
                (&[region(25, 100)], &[]),
                (&[region(50, 60)], &[]),
            ]),
            outcome(&[region(10, 40), region(50, 60)], &[], 1),
        );
    }

    /// **"Earlier" means handed over first, not starting first.** A later region delivering a
    /// locus that starts before the standing owner's cannot happen under `build_region`'s input
    /// contract; where it does, the standing owner keeps its ground.
    #[test]
    fn a_later_region_delivering_an_earlier_start_still_loses() {
        assert_eq!(
            hand_over(&[(&[region(30, 40)], &[]), (&[region(10, 35)], &[])]),
            outcome(&[region(30, 40)], &[], 1),
        );
    }

    /// **A displaced span is neither taken nor refused.** The refused list is what the width
    /// bound refused; displacement — which here drops one locus and one refusal — adds to the
    /// displaced count alone.
    #[test]
    fn a_displaced_span_is_counted_apart_from_the_refused() {
        assert_eq!(
            hand_over(&[
                (&[region(10, 300)], &[]),
                (&[region(50, 60)], &[region(70, 80)]),
            ]),
            outcome(&[region(10, 300)], &[], 2),
        );
    }

    /// A span whose ends are the wrong way round is read as the ground it names. Taken raw, its
    /// last base would put the frontier *behind* its own first base, and the locus that follows
    /// would be let through whether or not it overlapped.
    #[test]
    fn a_span_with_its_ends_the_wrong_way_round_still_owns_its_ground() {
        assert_eq!(
            hand_over(&[(&[region(50, 40)], &[]), (&[region(45, 60)], &[])]),
            outcome(&[region(50, 40)], &[], 1),
            "45 falls inside the ground 40-50 names, whichever way round it was written",
        );

        // And the other end: an inverted span *arriving* on owned ground must lose it. Read raw,
        // its first base would be 50 — past the frontier at 45 — and it would claim ground 40-45
        // that the first region already owns.
        assert_eq!(
            hand_over(&[(&[region(10, 45)], &[]), (&[region(50, 40)], &[])]),
            outcome(&[region(10, 45)], &[], 1),
        );
    }

    /// Ownership does not cross a contig: a locus at the start of the next contig is kept
    /// however far the last locus of the previous one reached.
    #[test]
    fn a_locus_on_the_next_contig_is_never_displaced() {
        let loci = [region_on(0, 10, 4_000), region_on(1, 5, 9)];
        assert_eq!(
            hand_over(&[(&loci[..1], &[]), (&loci[1..], &[])]),
            outcome(&loci, &[], 0),
        );
    }

    /// **The frontier is the run's, not one region's**, so it carries across the join between
    /// two regions and across a region that had nothing in it.
    #[test]
    fn the_frontier_carries_across_an_empty_region() {
        assert_eq!(
            hand_over(&[
                (&[region(10, 300)], &[]),
                (&[], &[]),
                (&[region(250, 260)], &[]),
            ]),
            outcome(&[region(10, 300)], &[], 1),
        );
    }

    /// Nothing is displaced in a run where nothing overlaps, which is every run the builders
    /// can actually produce.
    #[test]
    fn nothing_is_displaced_when_no_two_loci_overlap() {
        assert_eq!(
            hand_over(&[
                (&[region(1, 3), region(7, 9)], &[region(20, 40)]),
                (&[region(41, 44)], &[region(50, 70)]),
            ]),
            outcome(
                &[region(1, 3), region(7, 9), region(41, 44)],
                &[region(20, 40), region(50, 70)],
                0,
            ),
        );
    }
}
