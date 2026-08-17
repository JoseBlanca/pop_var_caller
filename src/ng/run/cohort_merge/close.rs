//! Closing cohort loci: the walk that groups the samples' observations by shared
//! reference bases.
//!
//! One pass over the samples' observations, merged by position, closing each locus as
//! its reach stops growing. This is the **grouping half** of production's group walk
//! (`derive_is_kept`, `var_calling/cohort_integration.rs`) with the columns removed, and
//! the design is `doc/devel/ng/spec/cohort_merge.md` §4.1 with the types in
//! `doc/devel/ng/arch/cohort_merge.md` §3.
//!
//! **Two things differ from production besides the columns, and both are deliberate.**
//! Production judges in the same loop — it decides which positions to keep as it groups
//! them — where this walk only groups and hands every closed locus out, judging being
//! the next step's (`doc/devel/ng/impl_plan/cohort_merge.md` A4). And the quantity
//! summed is not the same: production adds the per-position *maximum* over samples,
//! this walk adds every sample's own non-reference reads, which is the difference
//! [`super::MinAltObs`] sets out — a maximum never exceeds a sum, so at one threshold
//! the two rules keep different loci.

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomePosition, GenomeRegion};

/// One sample's observations inside one cohort locus — a borrowed run, never copied
/// out of what the walk was given (spec §6.4).
///
/// **A run, not a scatter.** Each sample's observations arrive in coordinate order and a
/// locus covers one stretch of it, so the ones inside a locus are contiguous and a slice
/// says so. The sample is named because the caller's evidence is indexed by the run's
/// sample order, and a locus most samples did not cover has no entry for them at all —
/// which is a different fact from an entry with reference-only support, and stays one
/// (spec §4.2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampleMembers<'a> {
    /// Index into the run's sample table — the order the walk was handed.
    pub sample: usize,
    /// This sample's observations inside the locus, in coordinate order. Never empty:
    /// a sample with nothing here has no `SampleMembers` at all.
    pub observations: &'a [SampleLocusObservations],
}

/// A cohort locus as it comes off the walk, before anything is judged or assembled.
#[derive(Debug, Clone, PartialEq)]
pub struct ClosedLocus<'a> {
    /// First position to furthest reach. Uncapped here — how wide a locus may be is a
    /// verdict passed on it afterwards, not a limit on closing it (spec §4.1).
    pub region: GenomeRegion,
    /// The samples that covered it, in ascending sample order. Only the samples that
    /// contributed appear.
    ///
    /// **Owned, where the architecture sketches a borrowed slice, and that is what
    /// [`Iterator`] costs.** An iterator's item cannot borrow a buffer the iterator
    /// owns, so a slice into closer-owned scratch is not expressible; what the
    /// architecture's promise is really about — the observations are borrowed, never
    /// copied out of the cache — holds either way, since each entry carries a slice. The
    /// price is one allocation per closed locus, sized exactly during the walk. If a
    /// caller ever wants that back, the shape to reach for is internal iteration, and
    /// A4 and B1 are where a real caller decides.
    pub members: Vec<SampleMembers<'a>>,
    /// Non-reference reads summed across every member — what the keep threshold is
    /// compared against (spec §4.3). Summed as the walk passes each observation, so
    /// nothing revisits the locus to compute it.
    ///
    /// Saturating, and it costs nothing: the only consumer compares the total against
    /// `min_alt_obs`, so a total pinned at `u32::MAX` gives the same verdict as the true
    /// one — and reaching the cap needs four billion non-reference reads in one locus.
    pub non_reference_reads: u32,
}

impl ClosedLocus<'_> {
    /// How many reference bases the locus covers.
    ///
    /// **Not `region.len()`**, which overflows its `+ 1` before subtracting when a locus
    /// reaches the top of the coordinate space — a panic in a debug build and a span of
    /// **0** in the release profile, where overflow checks are off (that method's own
    /// doc). Zero is the dangerous answer: this is what A4's width verdict compares
    /// against, so the widest locus expressible would be judged narrow enough to build.
    /// Subtracting first and adding after cannot do that.
    pub fn span(&self) -> u64 {
        self.region
            .end
            .get()
            .saturating_sub(self.region.start.get())
            .saturating_add(1)
    }
}

/// Walks the samples' observations merged by position and closes each locus as its
/// reach stops growing (spec §4.1).
///
/// **A locus closes when the next observation starts beyond the current reach.** Until
/// then each observation extends the reach if it goes further, and adds its
/// non-reference reads to the running total. Extending the reach *as the walk scans* is
/// what makes a second pass unnecessary however late a deletion widens the locus.
///
/// **One pass is enough only because the observations arrive merged by position.** A
/// deletion widens the locus while it is being closed, and the records the widening
/// pulls in sit at later positions in other samples — so the walk simply keeps taking
/// them while the next position is within the reach, each possibly pushing the reach
/// further again. Written as a loop over samples at each position instead, a widening
/// sends you back over samples already visited, repeatedly, until the span stops
/// growing: correct, but a fixpoint where this is one pass (spec §11's first trap).
///
/// A locus never crosses a contig, because no observation does.
///
/// **The four per-sample vectors are parallel on purpose**, and `over` is the one place
/// their lengths are set. Folding them into one vector of a per-sample struct would let
/// the compiler check that, at the cost of the thing the layout exists for: picking the
/// next observation scans `keys` and nothing else, and it is contiguous only while it is
/// its own array. The sibling merge makes the same trade for the same reason
/// (`MergedCursors`, `read/input/sample_cursor.rs`).
pub struct LocusCloser<'a> {
    /// One slice per sample, each in coordinate order.
    ///
    /// **The k slice references are copied in; the observations are not.** The cache
    /// hands its windows out for the length of a call (arch §4), so a closer that
    /// borrowed the outer array as well would tie every caller to keeping that array
    /// alive alongside the data. Copying k pointers once per walk costs nothing and
    /// leaves the borrow where it belongs — on the observations.
    observations_per_sample: Vec<&'a [SampleLocusObservations]>,
    /// How far each sample has been consumed — the cursor the argmin advances.
    cursors: Vec<usize>,
    /// `keys[sample]` is where `observations_per_sample[sample][cursors[sample]]`
    /// starts; `None` once that sample is spent. Refreshed only for the sample just
    /// advanced.
    ///
    /// **Keys beside the cursors, not read through them** — the layout `MergedCursors`
    /// settled on for the same argmin (`read/input/sample_cursor.rs`). Picking the next
    /// observation compares k keys, and holding them contiguously makes that a linear
    /// scan of one array in cache rather than k jumps into observations scattered across
    /// the heap. It matters more here than there, because k is the cohort size: this
    /// walk runs over every observation of every sample on every calling run.
    keys: Vec<Option<GenomePosition>>,
    /// Where each cursor stood when the open locus opened — scratch, refilled per locus
    /// and never reallocated, since its length is the sample count for the whole walk.
    cursors_at_open: Vec<usize>,
}

impl<'a> LocusCloser<'a> {
    /// Walk `observations_per_sample`, one slice per sample in the run's sample order,
    /// each in coordinate order.
    ///
    /// **Coordinate order within a sample is a precondition, and the walk panics rather
    /// than grouping ground it cannot group.** The generator mints in order and the cache
    /// preserves it, so this fires only on a defect upstream — but an unsorted sample
    /// does not fail loudly on its own: it closes a locus over ground it does not cover,
    /// and every verdict downstream would be passed on the wrong span. The check is one
    /// comparison per observation, inside the walk (see [`Iterator::next`]).
    ///
    /// Whatever it is handed, the walk terminates and consumes every observation exactly
    /// once.
    ///
    /// Zero samples yields nothing rather than failing: refusing a zero-sample *run* is
    /// the caller's job and happens at construction (spec §7.2), and a walk over an
    /// empty cohort is not the place to discover it.
    pub fn over(observations_per_sample: &[&'a [SampleLocusObservations]]) -> Self {
        let keys = observations_per_sample
            .iter()
            .map(|observations| observations.first().map(key_of))
            .collect();
        Self {
            cursors: vec![0; observations_per_sample.len()],
            cursors_at_open: vec![0; observations_per_sample.len()],
            observations_per_sample: observations_per_sample.to_vec(),
            keys,
        }
    }

    /// The head of sample `sample`, if it has one left.
    fn head(&self, sample: usize) -> Option<&'a SampleLocusObservations> {
        self.observations_per_sample[sample].get(self.cursors[sample])
    }

    /// Consume sample `sample`'s head and refresh just that sample's key.
    fn advance(&mut self, sample: usize) {
        self.cursors[sample] += 1;
        self.keys[sample] = self.head(sample).map(key_of);
    }

    /// The sample whose head starts earliest, **ties to the lowest sample index**.
    ///
    /// **Nothing in a closed locus depends on the tie-break today**, and saying so is
    /// the point: the members come out in ascending sample order because they are
    /// collected that way, the reach is a `max` and the total a sum, and *which*
    /// observations land inside a locus is fixed by the growing-reach condition rather
    /// than by the order they were visited in. A tie-break toward the highest sample
    /// index would give identical loci.
    ///
    /// The rule is fixed anyway, at no cost. It is the read layer's rule for the same
    /// argmin (`MergedCursors::argmin`, `read/input/sample_cursor.rs`), and a later step
    /// that emitted members in the order they were consumed — rather than collecting
    /// them by sample afterwards — would need it.
    fn sample_with_earliest_head(&self) -> Option<(usize, &'a SampleLocusObservations)> {
        let mut best: Option<(usize, GenomePosition)> = None;
        for (sample, key) in self.keys.iter().enumerate() {
            let Some(key) = *key else { continue };
            // Strictly less, so an equal key leaves the earlier sample in place — that
            // *is* the tie-break.
            if best.is_none_or(|(_, best_key)| key < best_key) {
                best = Some((sample, key));
            }
        }
        // The head comes back with the sample it was found in, so no caller looks it up a
        // second time or has to state an invariant the scan already knows.
        best.map(|(sample, _)| (sample, self.head(sample).expect("a key implies a head")))
    }
}

/// Where an observation starts, as a genome-wide key.
fn key_of(observation: &SampleLocusObservations) -> GenomePosition {
    GenomePosition {
        contig: observation.region.contig,
        position: observation.region.start,
    }
}

impl<'a> Iterator for LocusCloser<'a> {
    type Item = ClosedLocus<'a>;

    /// **The walk always advances.** The first observation it looks at is the one that
    /// opened the locus, whose start *is* the reach so far and whose contig *is* the
    /// locus's, so neither break condition can hold on the first pass: every call
    /// consumes at least one observation, and the iterator cannot spin. That is a
    /// property of the code rather than of the input, which is why the sortedness
    /// precondition on [`LocusCloser::over`] carries no progress obligation — and why
    /// production's release-level `StalledCut` guard has no counterpart here (spec §10,
    /// lesson 2).
    fn next(&mut self) -> Option<Self::Item> {
        let (_, opening) = self.sample_with_earliest_head()?;
        let contig = opening.region.contig;
        let start = opening.region.start;

        // Where each sample stands now. What it consumes before the locus closes is its
        // run of members, so no second scan is needed to find them.
        self.cursors_at_open.clear();
        self.cursors_at_open.extend_from_slice(&self.cursors);

        let mut reach = start;
        let mut non_reference_reads = 0u32;
        let mut covering_samples = 0usize;

        // Take heads while the next one begins inside the reach, extending the reach as
        // we go.
        while let Some((sample, head)) = self.sample_with_earliest_head() {
            // A locus never crosses a contig, and neither does an observation, so a
            // change of contig closes the locus whatever the positions say.
            if head.region.contig != contig || head.region.start > reach {
                break;
            }
            // **The argmin is a merge only if each sample's own slice ascends**, and
            // nothing upstream of here can be made to prove it. Handed a sample whose
            // observations arrive out of order, the walk still terminates and still
            // consumes each one exactly once — but it closes a locus over ground it does
            // not cover: `[50–50, 10–10]` in one sample comes out as one locus spanning
            // 50–50 holding an observation forty bases outside it. Every verdict
            // downstream would then be passed on the wrong span, and nothing could
            // detect it.
            //
            // A release check, not a `debug_assert!`: the release profile is the one
            // this repo runs, a trap it has recorded hitting twice
            // (`locus_generation/mod.rs`). The cost is one comparison against
            // per-observation work that already includes a base-by-base sequence
            // comparison in `non_reference_reads`.
            assert!(
                head.region.start >= start,
                "sample {sample}'s observations are not in coordinate order: {} starts \
                 before the locus that opened at {}",
                head.region,
                start.get(),
            );
            reach = reach.max(head.reach());
            non_reference_reads = non_reference_reads.saturating_add(head.non_reference_reads());
            // A sample's first observation inside this locus is exactly the one taken
            // while its cursor still stands where the locus opened — so the member count
            // is exact and free, and the vector below is allocated once.
            if self.cursors[sample] == self.cursors_at_open[sample] {
                covering_samples += 1;
            }
            self.advance(sample);
        }

        // **Progress, asserted rather than assumed.** Every call must consume at least
        // the observation that opened the locus, or the iterator yields empty loci
        // forever and the caller hangs. It holds because `reach` starts at that
        // observation's own start and the comparison is strict — three separately
        // editable facts, and the failure mode is the one class no test can catch by
        // failing: a mutation of the comparison during review did not fail an assertion,
        // it killed the test binary with SIGKILL. Production guards the same class at
        // release level (`StalledCut`); here the guarantee is structural (spec §10,
        // lesson 2), so a debug assert is the whole cost.
        debug_assert!(
            covering_samples > 0,
            "a locus closed at {contig:?}:{} without consuming anything",
            start.get()
        );

        let mut members = Vec::with_capacity(covering_samples);
        members.extend(
            self.cursors_at_open
                .iter()
                .enumerate()
                .filter(|&(sample, &from)| from < self.cursors[sample])
                .map(|(sample, &from)| SampleMembers {
                    sample,
                    observations: &self.observations_per_sample[sample][from..self.cursors[sample]],
                }),
        );

        Some(ClosedLocus {
            region: GenomeRegion {
                contig,
                start,
                end: reach,
            },
            members,
            non_reference_reads,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, Position, ReadGroupId};

    fn region_on(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        region_on(0, start, end)
    }

    fn sequence(bases: &[u8], num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// **One sample's observation over `region`** — the spec's sense of the word
    /// (§1.3): one sample's whole record over a stretch of genome, which is what a
    /// `SampleLocusObservations` is. Not a `SequenceObservation`, which is one distinct
    /// sequence inside it and which `sequence` above builds.
    ///
    /// It carries `non_reference_reads` reads of one non-reference sequence. The
    /// reference bases are `A` repeated and the alternative `C` repeated, so every
    /// observation here differs from the reference in the same uninteresting way — the
    /// walk under test reads positions and counts, never bases.
    fn observation(region: GenomeRegion, non_reference_reads: u32) -> SampleLocusObservations {
        let width = region.len().max(1) as usize;
        SampleLocusObservations {
            region,
            reference_bases: vec![b'A'; width].into_boxed_slice(),
            observations: vec![sequence(&vec![b'C'; width], non_reference_reads)],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// An observation whose reads all matched the reference — `reference_reads` of them,
    /// and no non-reference reads at all.
    fn all_reference_observation(
        region: GenomeRegion,
        reference_reads: u32,
    ) -> SampleLocusObservations {
        let width = region.len().max(1) as usize;
        // Every field spelled out rather than `..observation(…)`: a field added to
        // `SampleLocusObservations` should be a compile error here, not something filled
        // in silently from a base expression.
        SampleLocusObservations {
            region,
            reference_bases: vec![b'A'; width].into_boxed_slice(),
            observations: vec![sequence(&vec![b'A'; width], reference_reads)],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// The spans of every locus the walk closes, in order — what most of these tests
    /// assert on. **The contig is dropped**, so a test about contig boundaries must not
    /// use this.
    fn closed_spans(observations_per_sample: &[&[SampleLocusObservations]]) -> Vec<(u64, u64)> {
        LocusCloser::over(observations_per_sample)
            .map(|closed| (closed.region.start.get(), closed.region.end.get()))
            .collect()
    }

    /// **Neither observation covers the other's base, so they cannot be called apart** —
    /// spec §4.1's own example, and the case that separates the chaining rule from
    /// "anything nearby joins".
    #[test]
    fn two_adjacent_snps_are_two_loci() {
        let sample = [
            observation(region(10, 10), 4),
            observation(region(11, 11), 4),
        ];
        assert_eq!(closed_spans(&[&sample]), vec![(10, 10), (11, 11)]);
    }

    /// **A deletion covers the SNP's base, so they are one locus** — spec §4.1's other
    /// example. The two are in different samples, which is the case the merge exists for.
    #[test]
    fn a_deletion_and_a_snp_inside_it_are_one_locus() {
        let with_deletion = [observation(region(10, 14), 3)];
        let with_snp = [observation(region(12, 12), 2)];
        assert_eq!(
            closed_spans(&[&with_deletion, &with_snp]),
            vec![(10, 14)],
            "the deletion covers position 12, so nothing separates them"
        );
    }

    /// **The reach grows while the locus is being closed, and that is what buys the
    /// single pass.**
    ///
    /// Positions 10–12, then a deletion at 11 reaching 20, then a SNP at 15. The third
    /// observation is only inside the locus *because* the second widened it: at the
    /// moment the walk passed position 10 the reach was 12, and 15 is well beyond it. An
    /// implementation that fixed the reach when the locus opened would close 10–12 and
    /// start a new locus at 15; one that re-scanned from the beginning after each
    /// widening would give the same answer as this, but in as many passes as the locus
    /// has wideners.
    #[test]
    fn a_late_widening_pulls_in_a_following_observation() {
        let narrow = [observation(region(10, 12), 2)];
        let deletion = [observation(region(11, 20), 3)];
        let following = [observation(region(15, 15), 4)];

        assert_eq!(
            closed_spans(&[&narrow, &deletion, &following]),
            vec![(10, 20)],
            "15 is inside the locus only because the deletion widened it to 20"
        );
    }

    /// The members are the observations the locus actually covered, grouped by sample
    /// and in ascending sample order — and a sample that covered nothing has no entry,
    /// which is what keeps *no coverage* distinct from *reference-only coverage*.
    #[test]
    fn the_members_are_the_covering_samples_in_order() {
        let covering = [observation(region(10, 14), 3)];
        let inside = [observation(region(12, 12), 2)];
        let elsewhere = [observation(region(50, 50), 5)];

        let mut closer = LocusCloser::over(&[&covering, &inside, &elsewhere]);
        let first = closer.next().expect("a locus at 10");

        assert_eq!(first.region, region(10, 14));
        assert_eq!(
            first.members.iter().map(|m| m.sample).collect::<Vec<_>>(),
            vec![0, 1],
            "sample 2 covered nothing here and so has no entry at all"
        );
        assert_eq!(first.members[0].observations.len(), 1);
        assert_eq!(first.members[1].observations.len(), 1);

        let second = closer.next().expect("a locus at 50");
        assert_eq!(second.region, region(50, 50));
        assert_eq!(
            second.members.iter().map(|m| m.sample).collect::<Vec<_>>(),
            vec![2]
        );
        assert!(closer.next().is_none());
    }

    /// One sample's several observations inside one locus stay one run, in order — the
    /// shape `SampleMembers` promises, and the one assembly will project.
    #[test]
    fn a_samples_own_observations_inside_one_locus_are_one_run() {
        let wide = [observation(region(10, 20), 3)];
        let two_inside = [
            observation(region(12, 12), 2),
            observation(region(14, 14), 1),
        ];

        let mut closer = LocusCloser::over(&[&wide, &two_inside]);
        let closed = closer.next().expect("one locus");

        assert_eq!(closed.region, region(10, 20));
        assert_eq!(closed.members[1].sample, 1);
        assert_eq!(closed.members[1].observations.len(), 2);
        assert_eq!(closed.members[1].observations[0].region, region(12, 12));
        assert_eq!(closed.members[1].observations[1].region, region(14, 14));
        assert!(closer.next().is_none());
    }

    /// **A locus never crosses a contig**, whatever the positions say. Position 10 on
    /// contig 1 is not inside a locus reaching position 14 on contig 0, and comparing
    /// bare positions would say it is.
    #[test]
    fn a_locus_never_crosses_a_contig() {
        let sample = [
            observation(region_on(0, 10, 14), 3),
            observation(region_on(1, 10, 10), 2),
        ];
        let closed: Vec<_> = LocusCloser::over(&[&sample]).collect();

        assert_eq!(closed.len(), 2);
        assert_eq!(closed[0].region, region_on(0, 10, 14));
        assert_eq!(closed[1].region, region_on(1, 10, 10));
    }

    /// **Production's `over_approximation_is_max_not_sum`, inverted** — the divergence
    /// spec §4.3 chose deliberately and spec §15 requires pinned.
    ///
    /// Two samples, one non-reference read each, at the *same* position. Production sums
    /// the per-position maximum across samples and answers 1, so at the default
    /// `min_alt_obs` of 2 it drops the group; ng sums every sample's reads and answers 2,
    /// so the locus reaches the threshold. The shared position is what makes this the
    /// discriminating fixture: spread the same reads across different positions and the
    /// two rules agree.
    #[test]
    fn one_non_reference_read_in_each_of_two_samples_sums_to_two() {
        let first = [observation(region(20, 20), 1)];
        let second = [observation(region(20, 20), 1)];

        let closed: Vec<_> = LocusCloser::over(&[&first, &second]).collect();
        assert_eq!(closed.len(), 1);
        assert_eq!(
            closed[0].non_reference_reads, 2,
            "production's maximum over these same two samples is 1"
        );
    }

    /// The non-reference total is summed across every member of the locus, over samples
    /// and positions alike — the number A4's keep verdict is compared against.
    ///
    /// Three samples, 3 + 2 + 4 non-reference reads, and a fourth sample whose reads all
    /// matched the reference contributing none. Nine reference reads sit in that fourth
    /// sample, so an implementation summing `num_obs` rather than the non-reference part
    /// answers 18 here. **It does not separate ng's rule from production's** — the three
    /// non-reference observations sit at three different positions, where a per-position
    /// maximum also sums to 9; the test above is the one that does.
    #[test]
    fn the_non_reference_total_is_summed_across_the_locus() {
        let wide = [observation(region(10, 20), 3)];
        let inside = [observation(region(12, 12), 2)];
        let also_inside = [observation(region(14, 14), 4)];
        let all_reference = [all_reference_observation(region(16, 16), 9)];

        let closed: Vec<_> =
            LocusCloser::over(&[&wide, &inside, &also_inside, &all_reference]).collect();

        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].non_reference_reads, 9);
        assert_eq!(
            closed[0].members.len(),
            4,
            "the all-reference sample is a member: its depth counts even though it varied at nothing"
        );
    }

    /// **An observation starting exactly at the reach joins the locus** — the boundary
    /// where production's `positions[j] <= group_end` and a `<` would part company, and
    /// the one input that separates them.
    ///
    /// A deletion at 10 reaching 14, and a SNP at exactly 14. Under this rule they are
    /// one locus, because the deletion covers base 14 and nothing can call them apart;
    /// under a strict `<` the SNP would open a second locus over ground the first
    /// already owns. The SNP at 15, one base further, is what says the rule stops where
    /// the reach does.
    #[test]
    fn an_observation_starting_exactly_at_the_reach_joins() {
        let deletion = [observation(region(10, 14), 3)];
        let at_the_reach = [observation(region(14, 14), 2)];
        assert_eq!(closed_spans(&[&deletion, &at_the_reach]), vec![(10, 14)]);

        let one_past_the_reach = [observation(region(15, 15), 2)];
        assert_eq!(
            closed_spans(&[&deletion, &one_past_the_reach]),
            vec![(10, 14), (15, 15)]
        );
    }

    /// **A locus at the top of the coordinate space reports its true width**, where
    /// `region.len()` would panic in debug and answer 0 in release — and 0 is the answer
    /// that matters, since A4 compares this against `max_cohort_locus_span` and would
    /// judge the widest locus expressible narrow enough to build.
    #[test]
    fn the_span_is_the_base_count_at_the_coordinate_ceiling() {
        let at_the_ceiling = ClosedLocus {
            region: region(u64::MAX - 3, u64::MAX),
            members: Vec::new(),
            non_reference_reads: 0,
        };
        assert_eq!(at_the_ceiling.span(), 4);
    }

    /// **The members come back in sample order, not in the order the walk met them.**
    /// Sample 1 opens this locus at 10 and sample 0 joins inside it at 12, so a closer
    /// that emitted members as it consumed them would answer `[1, 0]` — and every other
    /// fixture in this file happens to make the two orders agree, so nothing else here
    /// could tell.
    #[test]
    fn the_members_are_in_sample_order_not_consumption_order() {
        let joins_later = [observation(region(12, 12), 2)];
        let opens_the_locus = [observation(region(10, 14), 3)];

        let closed: Vec<_> = LocusCloser::over(&[&joins_later, &opens_the_locus]).collect();
        assert_eq!(closed.len(), 1);
        assert_eq!(
            closed[0]
                .members
                .iter()
                .map(|m| m.sample)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
    }

    /// The total saturates rather than wrapping. Toward the top of the committed cohort
    /// range a locus can carry more non-reference reads than a `u32` holds, and in the
    /// release profile — overflow checks off — a wrapped total would fall *below*
    /// `min_alt_obs` and drop the deepest loci in the cohort, counting nothing.
    #[test]
    fn the_non_reference_total_saturates_at_the_u32_ceiling() {
        let first = [observation(region(10, 10), u32::MAX)];
        let second = [observation(region(10, 10), u32::MAX)];

        let closed: Vec<_> = LocusCloser::over(&[&first, &second]).collect();
        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].non_reference_reads, u32::MAX);
    }

    /// A sample that covered nothing here sits between two that did. It is skipped rather
    /// than treated as a head at position zero, and it appears in no locus's members — the
    /// ordinary case at a few reads a position, not an edge.
    #[test]
    fn an_empty_sample_among_covered_ones_is_skipped() {
        let before: [SampleLocusObservations; 0] = [];
        let covered = [observation(region(10, 14), 3)];
        let after: [SampleLocusObservations; 0] = [];

        let closed: Vec<_> = LocusCloser::over(&[&before, &covered, &after]).collect();

        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].region, region(10, 14));
        assert_eq!(
            closed[0]
                .members
                .iter()
                .map(|m| m.sample)
                .collect::<Vec<_>>(),
            vec![1],
            "the empty samples contribute no entry at all"
        );
    }

    /// The merge is keyed genome-wide, not by bare position: contig 1's position 5 comes
    /// *after* contig 0's position 10, and a walk comparing positions alone would open it
    /// first and emit the two loci in the wrong order.
    #[test]
    fn the_merge_orders_by_contig_before_position() {
        let on_contig_zero = [observation(region_on(0, 10, 14), 3)];
        let on_contig_one = [observation(region_on(1, 5, 5), 2)];

        let closed: Vec<_> = LocusCloser::over(&[&on_contig_zero, &on_contig_one]).collect();

        assert_eq!(closed.len(), 2);
        assert_eq!(closed[0].region, region_on(0, 10, 14));
        assert_eq!(closed[1].region, region_on(1, 5, 5));
    }

    /// **A sample whose observations are out of order is refused, not quietly grouped.**
    ///
    /// Without the check the walk still terminates and still consumes each observation
    /// once — it closes one locus spanning 50–50 holding an observation at 10–10, forty
    /// bases outside the region that claims to cover it. A4's width verdict would then be
    /// passed on the wrong span, and nothing downstream could tell.
    #[test]
    #[should_panic(expected = "not in coordinate order")]
    fn a_sample_out_of_coordinate_order_is_refused() {
        let unsorted = [
            observation(region(50, 50), 1),
            observation(region(10, 10), 1),
        ];
        let _ = LocusCloser::over(&[&unsorted]).count();
    }

    /// **The argmin scans every sample, not the first two.** The observation at 50 is
    /// carried only by the third sample, so a merge that compared two heads would never
    /// reach it — the sibling read-layer merge keeps a test for exactly this shape
    /// (`MergedCursors`, `read/input/sample_cursor.rs`).
    #[test]
    fn the_argmin_scans_every_sample() {
        let first = [observation(region(10, 10), 1)];
        let second = [observation(region(20, 20), 1)];
        let third = [observation(region(50, 50), 1)];

        assert_eq!(
            closed_spans(&[&first, &second, &third]),
            vec![(10, 10), (20, 20), (50, 50)]
        );
    }

    /// A cohort where nobody varied still closes its loci — the walk groups ground, and
    /// judging it empty is a verdict passed afterwards (spec §4.3), not something the
    /// closing hides.
    #[test]
    fn an_all_reference_locus_still_closes() {
        let sample = [all_reference_observation(region(10, 10), 4)];
        let closed: Vec<_> = LocusCloser::over(&[&sample]).collect();

        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].non_reference_reads, 0);
    }

    /// Nothing to walk is not a failure: an empty cohort and a cohort of empty samples
    /// both yield nothing. Refusing a zero-sample run happens at the caller's
    /// construction (spec §7.2), not here.
    #[test]
    fn nothing_to_walk_yields_nothing() {
        assert_eq!(LocusCloser::over(&[]).count(), 0);

        let empty: [SampleLocusObservations; 0] = [];
        assert_eq!(LocusCloser::over(&[&empty, &empty]).count(), 0);
    }

    /// **Where the samples' observations were split between slices changes nothing** —
    /// the closing rule reads observations, not which sample carried them (spec §9).
    /// Here the same five observations are dealt out three ways, and the loci come back
    /// identical.
    #[test]
    fn the_loci_do_not_depend_on_which_sample_carried_what() {
        let all = [
            observation(region(10, 12), 2),
            observation(region(11, 20), 3),
            observation(region(15, 15), 4),
            observation(region(40, 40), 1),
            observation(region(41, 41), 1),
        ];
        let expected = vec![(10, 20), (40, 40), (41, 41)];

        // One sample carrying everything.
        assert_eq!(closed_spans(&[&all]), expected);

        // Split so that the widener and what it pulls in sit in different samples.
        let a = [
            observation(region(10, 12), 2),
            observation(region(40, 40), 1),
        ];
        let b = [
            observation(region(11, 20), 3),
            observation(region(41, 41), 1),
        ];
        let c = [observation(region(15, 15), 4)];
        assert_eq!(closed_spans(&[&a, &b, &c]), expected);

        // And with the samples in a different order, which the tie-break makes
        // irrelevant to the spans.
        assert_eq!(closed_spans(&[&c, &b, &a]), expected);
    }

    /// Two samples with an observation at the same position is the ordinary case, and
    /// the members come back in sample order whichever slice is passed first.
    #[test]
    fn observations_at_the_same_position_close_as_one_locus() {
        let first = [observation(region(10, 10), 2)];
        let second = [observation(region(10, 10), 3)];

        let closed: Vec<_> = LocusCloser::over(&[&first, &second]).collect();
        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].region, region(10, 10));
        assert_eq!(
            closed[0]
                .members
                .iter()
                .map(|m| m.sample)
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
        assert_eq!(closed[0].non_reference_reads, 5);
    }

    /// An insertion's reference span is its anchor base, so it cannot widen a locus
    /// however long the inserted sequence is (spec §3.1). Two insertions at adjacent
    /// anchors are two loci, exactly as two SNPs are.
    #[test]
    fn an_insertion_does_not_widen_a_locus() {
        let mut inserted = observation(region(10, 10), 3);
        inserted.observations = vec![sequence(b"ACCCCCCCCCT", 3)];
        let sample = [inserted, observation(region(11, 11), 2)];

        assert_eq!(closed_spans(&[&sample]), vec![(10, 10), (11, 11)]);
    }

    /// Closing is uncapped: a chain of overlapping observations grows as long as the
    /// data keeps overlapping, and how wide that is allowed to be is A4's verdict rather
    /// than a limit here (spec §4.1).
    #[test]
    fn closing_is_uncapped() {
        let chained: Vec<_> = (0..40)
            .map(|step| observation(region(10 + step * 2, 12 + step * 2), 1))
            .collect();
        let closed: Vec<_> = LocusCloser::over(&[&chained]).collect();

        assert_eq!(closed.len(), 1);
        assert_eq!(closed[0].region, region(10, 90));
        assert_eq!(closed[0].span(), 81);
    }

    /// Consecutive loci never share a base, which is what makes every observation a
    /// member of exactly one locus (spec §4.2) and lets position order the loci
    /// totally.
    #[test]
    fn closed_loci_are_disjoint_and_ascending() {
        let a = [
            observation(region(10, 12), 1),
            observation(region(13, 13), 1),
            observation(region(30, 31), 1),
        ];
        let b = [
            observation(region(11, 11), 1),
            observation(region(14, 20), 1),
        ];

        let closed: Vec<_> = LocusCloser::over(&[&a, &b]).collect();
        for pair in closed.windows(2) {
            assert!(
                pair[0].region.end < pair[1].region.start,
                "{:?} and {:?} share ground",
                pair[0].region,
                pair[1].region
            );
        }

        let total_members: usize = closed
            .iter()
            .flat_map(|locus| locus.members.iter())
            .map(|members| members.observations.len())
            .sum();
        assert_eq!(
            total_members, 5,
            "every observation is a member exactly once"
        );
    }
}
