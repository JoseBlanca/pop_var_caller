//! One STR locus, reduced to the entry a stratum's table is keyed on.
//!
//! This is the only place that decides what a read's repeat count *is*, which the locus
//! type cannot answer: a locus carries its reference bases, its motif and the sequence each
//! read showed, and turning those into "this read sits two copies short of the reference"
//! is a division that has to be done somewhere and done once.
//!
//! Its own file because the shaping of data and the mathematics on it never live together
//! (`arch/parameter_prepass_ssr.md` §Module home) — [`super::slippage`] is the mathematics.
//!
//! Three reductions live here, in the order a locus meets them: which stratum it belongs to at
//! all ([`stratum_of`]), how one read group's reads at it fall across the offset buckets
//! ([`tally_of`]), and the thinning of that tally down to the read cap ([`shape_of`]) that turns
//! it into the entry a stratum's table is keyed on. The architecture sketches the last two as one
//! function; they are two because the caller needs what the tally counts *besides* the shape —
//! the reads it left out — and because the cap's draw needs the whole tally to draw from.

use crate::ng::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::ng::parameter_estimation::subsample::{SelectionWalk, seed_at};
use crate::ng::types::{Bp, GenomeRegion, ReadGroupId, SsrPeriod};

use super::stratum_table::LocusShape;
use super::{MAX_LOCUS_READS, OFFSET_BUCKETS, RepeatCount, Stratum, WholeRepeatOffset, bucket_of};

/// Which stratum a locus belongs to, or why none does.
///
/// **Three answers and not two**, where the architecture sketches an `Option`: the two ways a
/// locus can fail to have a stratum mean opposite things to the accumulator that asks. A locus
/// that is not one repeat tract is the SNP/indel path's business, so it is passed over in
/// silence; a repeat tract whose reference length is *not a whole number of motif copies* is a
/// delimiting fault in whatever admitted it, so it is counted and reported. The counter it
/// feeds — `SsrAccumulationCounts::loci_without_whole_repeat_reference`, which
/// `arch/parameter_prepass_ssr.md` §4 says should be near zero — is named after this variant
/// for that reason. An `Option` would collapse the two, and the counter would then have to
/// re-derive from the locus what this function had already decided.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum LocusStratum {
    /// The stratum this locus's **reference** tract puts it in.
    Stratified(Stratum),
    /// Not one delimited repeat tract, so this path does not fit it: a SNP/indel candidate, or
    /// a repeat cluster with no clean flanks. Nothing to record — the STR accumulator passes it
    /// over and the generic one has its own.
    NotOneRepeatTract,
    /// A repeat tract whose reference length is not a whole number of motif copies, so no
    /// stratum holds it: 13 bases of a 3-base motif is four copies and a base, and rounding
    /// either way would file it with tracts it is not comparable to.
    ///
    /// **Counted and skipped, never rounded**, and a run where this is not near zero is a bug
    /// report against the classification that admitted the locus rather than something this
    /// unit absorbs. Named to match the counter it feeds.
    WithoutWholeRepeatReference {
        /// How long the reference tract is.
        tract_len: Bp,
        /// The motif period it failed to divide by.
        period: SsrPeriod,
    },
}

/// Which stratum a locus belongs to — its reference tract's motif period and how many whole
/// copies of that motif the tract holds.
///
/// **A pure function of the reference**, which is the property the whole design rests on: every
/// sample files a locus under the same stratum whatever its own reads showed, so one sample's
/// stutter can be compared with another's and a cohort can pool them
/// (`spec/parameter_prepass_ssr.md` §4.1). A sample whose alleles differ from the reference does
/// not move between strata; its reads land at an offset instead.
///
/// Both halves come from the locus itself — `reference_bases` and the motif the STR generator
/// put in the locus's own detail — so nothing here re-reads the reference or the catalog.
///
/// **A tract of no bases is filed at zero copies rather than reported.** Zero divides by every
/// period, and zero copies is a stratum such tracts could be compared inside, where a
/// fractional count is comparable with nothing. Region typing emits none — the catalog's
/// loosest floor is three copies — so this is an answer to arithmetic rather than to data.
///
/// # Panics
///
/// If a tract holds more than `u32::MAX` whole copies. Unreachable twice over: the reference
/// bases are held in memory, so that would be a 4.29-billion-base allocation, and the catalog
/// refuses to serve any tract past `CATALOG_MAX_STR_LEN_BP` — 500 bases, against step 3's
/// calling default of 100 — so a tract that arrives here holds at most a few hundred copies.
/// It panics rather than saturating because saturating would file the locus under a stratum it
/// does not belong to and say nothing, which is the one thing this function is written not to
/// do.
#[must_use]
pub fn stratum_of(locus: &SampleLocusObservations) -> LocusStratum {
    // Matched kind by kind rather than through a wildcard, so that a locus kind added later
    // is a compile error here and its routing is a decision rather than a default: the one
    // outcome this function leaves no record of is the locus it passes over.
    let detail = match &locus.kind {
        LocusKind::Ssr(detail) => detail,
        LocusKind::Generic | LocusKind::SsrBundle => return LocusStratum::NotOneRepeatTract,
    };

    let period = detail.motif.ssr_period();
    let tract_bases = locus.reference_bases.len();
    let motif_bases = usize::from(period.get());

    if !tract_bases.is_multiple_of(motif_bases) {
        return LocusStratum::WithoutWholeRepeatReference {
            tract_len: Bp(tract_bases as u64),
            period,
        };
    }

    let repeats = u32::try_from(tract_bases / motif_bases)
        .expect("a tract holds at most a few hundred copies; see this function's panic note");
    LocusStratum::Stratified(Stratum::new(period, RepeatCount(repeats)))
}

/// One read group's reads at one locus, counted into the offset buckets — **before the read cap
/// thins them**.
///
/// It is what a [`LocusShape`](super::stratum_table::LocusShape) is made from, and it is a
/// separate object because a locus can carry hundreds of reads while a shape holds at most
/// `MAX_LOCUS_READS`: the thinning is a draw, and a draw needs the whole tally to draw from.
///
/// **The reads it leaves out are as much a decision as the ones it counts**, so it carries the
/// count of those too — see [`reads_with_partial_witness`](Self::reads_with_partial_witness).
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct OffsetTally {
    /// Reads at each bucket, in bucket order, with the ends absorbing everything beyond them.
    reads_by_bucket: [u32; OFFSET_BUCKETS],
    /// Reads whose length differs from the reference tract's by something other than a whole
    /// number of motif copies.
    reads_not_whole_repeat: u32,
    /// Reads of this group that witnessed only part of the tract, and so were left out.
    reads_with_partial_witness: u32,
}

impl OffsetTally {
    /// Reads at each bucket, in bucket order.
    #[inline]
    #[must_use]
    pub fn reads_by_bucket(self) -> [u32; OFFSET_BUCKETS] {
        self.reads_by_bucket
    }

    /// Reads whose length differs from the reference tract's by something other than a whole
    /// number of motif copies.
    #[inline]
    #[must_use]
    pub fn reads_not_whole_repeat(self) -> u32 {
        self.reads_not_whole_repeat
    }

    /// Reads of this group that saw only part of the tract, so their length is a **lower
    /// bound** — left out of the tally, and counted here because leaving them out is a
    /// decision a run has to be able to see the size of.
    ///
    /// A large share means the reads are short against these tracts, which is a fact about the
    /// library rather than about the repeats.
    #[inline]
    #[must_use]
    pub fn reads_with_partial_witness(self) -> u32 {
        self.reads_with_partial_witness
    }

    /// Every read this tally counted, guard included and partial witnesses excluded — **the
    /// depth a shape built from it would have before the cap thinned it**, which is the same
    /// quantity `LocusShape::depth` reports afterwards and is named to match.
    #[must_use]
    pub fn depth(self) -> u32 {
        self.reads_by_bucket
            .iter()
            .copied()
            .chain([self.reads_not_whole_repeat])
            .sum()
    }
}

/// Count one read group's reads at one locus into the offset buckets, by how far each read's
/// tract sits from the **reference** tract's length in whole motif copies.
///
/// **Complete witnesses only, and that is the whole reason this function is worth its own
/// step.** A read that reached both borders of the tract measured its length; a read that
/// reached only one border saw *part* of it, so what it shows is a **lower bound**. Scoring a
/// lower bound as a length reads as a read that lost repeats — a direct bias in the direction
/// split, which is the one parameter this path exists to protect and the one that inverts on
/// real data when an estimator gets it wrong (`spec/parameter_prepass_ssr.md` §3). The partial
/// witnesses are counted instead, in
/// [`OffsetTally::reads_with_partial_witness`](OffsetTally::reads_with_partial_witness).
///
/// **`reads_without_observation` never enters a depth** and is not read here: those reads
/// covered the tract and witnessed nothing at all, so they are neither a length nor a lower
/// bound. The accumulator counts them from the locus itself.
///
/// **A locus whose depth cap fired is tallied as it stands**, and its reads are not thinned
/// here: the generator's reservoir is a uniform subsample, so a capped locus is the same locus
/// observed at a lower depth rather than a locus to skip — and skipping it would drop deep loci,
/// which is depth-dependent selection, the bias this whole step exists to remove. What the cap
/// took is the accumulator's to report, not this function's to react to.
///
/// A read's offset is `(what the read showed − the reference tract) / period`, whole only when
/// that difference divides by the period; a read whose length differs by anything else goes to
/// the guard. Offsets past the recorded range land in the two saturating end buckets, so a read
/// forty copies short is an observation this keeps rather than one it drops.
///
/// **The period comes from the locus's own motif rather than from the caller**, so it cannot
/// disagree with the one [`stratum_of`] used: a read counted in a different unit from the
/// stratum it is filed under would be a plausible tally measured against the wrong origin, with
/// nothing to say so. A locus that is not one repeat tract has no motif and no offsets, so it
/// tallies as empty.
#[must_use]
pub fn tally_of(locus: &SampleLocusObservations, read_group: ReadGroupId) -> OffsetTally {
    let mut tally = OffsetTally::default();
    let LocusKind::Ssr(detail) = &locus.kind else {
        return tally;
    };
    // Both lengths are slices held in memory, so neither can exceed `i64::MAX` on any machine
    // that could hold them; an `expect` rather than a saturation because a saturated length
    // would make the subtraction below overflow, silently, in a release build.
    let reference_bases =
        i64::try_from(locus.reference_bases.len()).expect("a tract held in memory fits an i64");
    let motif_bases = i64::from(detail.motif.ssr_period().get());

    // The counters below add with a plain `+=`, and cannot wrap: every increment is one
    // observation's `num_obs`, and their sum over a locus is the reads the generator kept there
    // — capped at `DEFAULT_SSR_MAX_READS_PER_LOCUS`, a thousand, and bounded in any case by the
    // reads it holds in memory.
    for observation in &locus.observations {
        if observation.read_group != read_group {
            continue;
        }
        if observation.read_witness != ReadWitness::Complete {
            tally.reads_with_partial_witness += observation.num_obs;
            continue;
        }

        let shown_bases = i64::try_from(observation.bases.len())
            .expect("a read's bases are held in memory and fit an i64");
        let difference = shown_bases - reference_bases;

        if difference % motif_bases == 0 {
            // Saturating into the end buckets, which is what they are for: a read forty copies
            // short of the reference is a real observation — a locus carrying a long deletion —
            // and `bucket_of` absorbs anything past the recorded range anyway. What the clamp
            // adds is that the arithmetic cannot overflow the offset type on the way there.
            let offset = (difference / motif_bases).clamp(i64::from(i8::MIN), i64::from(i8::MAX));
            let bucket = bucket_of(WholeRepeatOffset(offset as i8));
            tally.reads_by_bucket[bucket.index()] += observation.num_obs;
        } else {
            tally.reads_not_whole_repeat += observation.num_obs;
        }
    }

    tally
}

/// The shape a locus enters a stratum's table as, and the depth it was thinned from if the
/// read cap fired.
///
/// **The depth it came from is carried rather than recomputed**, because it is the number the
/// accumulator reports: a run where most loci were thinned is a run whose depths are the cap's
/// and not the data's, and nothing downstream could tell from the shape alone.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct EnteredShape {
    shape: LocusShape,
    subsampled_from: Option<u32>,
}

impl EnteredShape {
    /// The shape itself — what the table is keyed on.
    #[inline]
    #[must_use]
    pub fn shape(self) -> LocusShape {
        self.shape
    }

    /// The depth this locus was thinned from, where the cap fired; `None` where every read it
    /// showed was entered.
    #[inline]
    #[must_use]
    pub fn subsampled_from(self) -> Option<u32> {
        self.subsampled_from
    }
}

/// Reduce one read group's tally at one locus to the shape its table is keyed on, thinning to
/// [`MAX_LOCUS_READS`] where the reads are deeper — `None` where the tally holds no reads at
/// all, since a locus with no evidence is one to leave out rather than to enter empty.
///
/// **It takes the tally rather than the locus**, where the architecture sketches one function
/// from the locus. Two reasons, and the second is the load-bearing one: the caller has the
/// tally already — it needs the partial-witness count the tally carries and this function does
/// not — and re-deriving it here would walk the locus's observations twice per locus on a walk
/// of millions. The position comes separately because it is the only other thing the draw may
/// depend on.
///
/// **The thinning is a uniform subsample, and a subsample is exact rather than approximate**:
/// it leaves the bucket counts distributed exactly as they would have been at the lower depth,
/// so what the cap costs is precision and never a bias. Dropping the deep loci instead would be
/// depth-dependent selection, which is the bias step 4 exists to remove.
///
/// **The draw is seeded from the locus's position and from nothing else, and that is this
/// step's whole difficulty.** Seeded so, a region-sharded walk and a single-threaded one keep
/// the same reads at every locus, and merging their tables stays an equality; seeded from a
/// counter, a thread or the clock, the two diverge by a few reads at each deep locus — a
/// difference no test that does not compare two whole walks would ever show, and one that would
/// make a fitted level move with the thread count. **That whole-walk comparison cannot be
/// written here**, because nothing at this step has a shard to lay out differently: it belongs
/// to the accumulator, and what this file can hold is that the draw is reproducible, that it
/// moves with the position, and that it carries no state between calls.
///
/// The population is the tally's own depth — the reads that measured a length. A read that saw
/// only part of the tract is a lower bound rather than a length, so it is not in the tally and
/// must not be in the population either: drawing twelve from a population the walk never walks
/// would enter the locus with fewer reads than the cap, silently, since a short shape is a legal
/// shape.
///
/// The buckets are drawn from in order and the guard last, through one walk over that whole
/// depth, so exactly `MAX_LOCUS_READS` reads survive and each is as likely to be one of them as
/// any other.
#[must_use]
pub fn shape_of(tally: OffsetTally, at: GenomeRegion) -> Option<EnteredShape> {
    let depth = tally.depth();
    if depth == 0 {
        return None;
    }

    let (reads_by_bucket, reads_not_whole_repeat, subsampled_from) = if depth <= MAX_LOCUS_READS {
        (
            tally.reads_by_bucket(),
            tally.reads_not_whole_repeat(),
            None,
        )
    } else {
        let mut walk = SelectionWalk::new(seed_at(at), depth, MAX_LOCUS_READS);
        let mut kept = [0; OFFSET_BUCKETS];
        for (bucket, reads) in tally.reads_by_bucket().into_iter().enumerate() {
            kept[bucket] = walk.keep_from(reads);
        }
        // The guard last, so the walk covers the locus's whole depth exactly once — the reads
        // that showed a non-whole-repeat length are as eligible for the draw as any other, and
        // leaving them out of it would make the guard share a property of the depth.
        let kept_guard = walk.keep_from(tally.reads_not_whole_repeat());
        (kept, kept_guard, Some(depth))
    };

    // PANIC-FREE: the counts are either the tally's, which the branch above holds at or below
    // the cap, or a draw of exactly `MAX_LOCUS_READS` from it; and the depth is at least one,
    // which the early return above establishes.
    let shape = LocusShape::try_new(reads_by_bucket, reads_not_whole_repeat)
        .expect("a capped tally holds between one read and the cap");
    Some(EnteredShape {
        shape,
        subsampled_from,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusLen, SequenceObservation, SsrDetail};
    use crate::ng::types::{ContigId, GenomeRegion, MAX_MOTIF_LEN, Motif, Position, ReadGroupId};

    /// A locus over `reference_bases`, tiled by `motif`, carrying whatever reads a test wants.
    fn ssr_locus(
        reference_bases: &[u8],
        motif: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1_000),
                end: Position(1_000 + reference_bases.len().max(1) as u64 - 1),
            },
            reference_bases: Box::from(reference_bases),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(motif).expect("a motif inside the STR period range"),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        }
    }

    /// One read's worth of observation, at whatever length it showed.
    fn observation(bases: &[u8], num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: num_obs / 2,
            q_sum: -10.0,
            mapq_sum: 60 * num_obs,
            mapq_sum_sq: 3_600 * u64::from(num_obs),
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn stratum(period: u8, repeats: u32) -> Stratum {
        Stratum::new(
            SsrPeriod::try_new(usize::from(period)).expect("a period inside the STR scope"),
            RepeatCount(repeats),
        )
    }

    /// The tract's length over its motif's period, which is the whole of the arithmetic: twenty
    /// bases of a dinucleotide is a stratum of ten repeats.
    #[test]
    fn a_tract_is_stratified_by_its_reference_length_over_its_motif_period() {
        let locus = ssr_locus(b"ATATATATATATATATATAT", b"AT", Vec::new());

        assert_eq!(stratum_of(&locus), LocusStratum::Stratified(stratum(2, 10)));
    }

    /// **The stratum is a property of the reference alone**, which is what lets a cohort compare
    /// one sample's stutter with another's. Four samples that saw completely different things at
    /// the same tract — every read at the reference length, every read two copies short, every
    /// read two copies long, and nothing at all — file it under the same stratum.
    ///
    /// **The expansion is the sample a hand-written fixture forgets**, and the one this path's
    /// own measurements make rarest: a read is about five times as likely to have lost a copy as
    /// gained one, so a length rule reading the *longest* thing it can see would pass every
    /// other fixture here. A sample carrying an expansion that migrated to another stratum would
    /// be compared against nothing.
    #[test]
    fn every_sample_files_the_same_tract_under_the_same_stratum_whatever_its_reads_showed() {
        let reference = b"ATATATATATATATATATAT";
        let at_reference = ssr_locus(reference, b"AT", vec![observation(reference, 9)]);
        let two_copies_short =
            ssr_locus(reference, b"AT", vec![observation(b"ATATATATATATATAT", 9)]);
        let two_copies_long = ssr_locus(
            reference,
            b"AT",
            vec![observation(b"ATATATATATATATATATATATAT", 9)],
        );
        let saw_nothing = ssr_locus(reference, b"AT", Vec::new());

        let expected = LocusStratum::Stratified(stratum(2, 10));
        assert_eq!(stratum_of(&at_reference), expected);
        assert_eq!(stratum_of(&two_copies_short), expected);
        assert_eq!(stratum_of(&two_copies_long), expected);
        assert_eq!(stratum_of(&saw_nothing), expected);
    }

    /// **The divisor is the motif's own period at every period the STR scope admits**, not only
    /// at the mononucleotides and dinucleotides the other fixtures use. Periods 1 to 6 are all
    /// ordinary loci, and a rule that clamped the divisor — at 2, say — would file a
    /// sixteen-base tetranucleotide tract at eight copies instead of four, pooling it with a
    /// stratum whose reads are measured from a different origin, with nothing to say so.
    #[test]
    fn a_tract_is_stratified_by_its_own_period_at_every_period_the_str_scope_admits() {
        for (bases, motif, period, repeats) in [
            (&b"AAAAAAAA"[..], &b"A"[..], 1u8, 8u32),
            (&b"ATATATATATATATATATAT"[..], &b"AT"[..], 2, 10),
            (&b"CAGCAGCAGCAGCAG"[..], &b"CAG"[..], 3, 5),
            (&b"AAAGAAAGAAAGAAAG"[..], &b"AAAG"[..], 4, 4),
            (&b"AAAAGAAAAGAAAAG"[..], &b"AAAAG"[..], 5, 3),
            (&b"AAAAAGAAAAAGAAAAAG"[..], &b"AAAAAG"[..], 6, 3),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(bases, motif, Vec::new())),
                LocusStratum::Stratified(stratum(period, repeats)),
                "{} bases of a {period}-base motif",
                bases.len()
            );
        }
    }

    /// A tract that is exactly its own motif is one copy — not zero, and not "not a repeat".
    /// The catalog's floors emit none this short, but the arithmetic answers for them, and the
    /// mononucleotide loop below reaches one copy only at period 1, where a tract's length and
    /// its repeat count are the same number and so cannot separate an off-by-one from a
    /// division.
    #[test]
    fn a_tract_of_exactly_one_motif_copy_is_a_stratum_of_one_repeat() {
        for (motif, period) in [
            (&b"AT"[..], 2u8),
            (&b"CAG"[..], 3),
            (&b"AAAG"[..], 4),
            (&b"AAAAG"[..], 5),
            (&b"AAAAAG"[..], 6),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(motif, motif, Vec::new())),
                LocusStratum::Stratified(stratum(period, 1)),
                "one copy of a {period}-base motif"
            );
        }
    }

    /// **A tract that is not a whole number of copies is counted and skipped, not rounded.**
    /// Thirteen bases of a three-base motif is four copies and a base; filing it at four would
    /// pool it with clean four-copy tracts, and at five with clean five-copy ones, and its
    /// reads' offsets would be measured from a length no allele has.
    #[test]
    fn a_tract_whose_length_is_not_a_whole_number_of_copies_is_reported_rather_than_rounded() {
        let locus = ssr_locus(b"CAGCAGCAGCAGC", b"CAG", Vec::new());

        assert_eq!(
            stratum_of(&locus),
            LocusStratum::WithoutWholeRepeatReference {
                tract_len: Bp(13),
                period: SsrPeriod::try_new(3).expect("a trinucleotide"),
            },
            "13 bases of a 3-base motif is four copies and a base"
        );
    }

    /// The report names the two numbers a reader needs to find the fault upstream: how long the
    /// tract was, and what it failed to divide by. A count alone would say a locus was skipped
    /// and nothing about which.
    ///
    /// **Every remainder, not only one base over.** A rule that reported a tract only when it
    /// sat one base past a whole number of copies would floor all the others silently — and the
    /// counter that is supposed to catch exactly that would read zero while it happened.
    #[test]
    fn the_report_names_the_tract_length_and_the_period_it_did_not_divide_by() {
        for (bases, motif, tract_len, period) in [
            (&b"CAGCAGCAGCAGC"[..], &b"CAG"[..], 13, 3),
            (&b"ATATATATA"[..], &b"AT"[..], 9, 2),
            (&b"AAAGAAAGAAAGA"[..], &b"AAAG"[..], 13, 4),
            (&b"CAGCAGCAGCAGCA"[..], &b"CAG"[..], 14, 3),
            (&b"AAAGAAAGAAAGAA"[..], &b"AAAG"[..], 14, 4),
            (&b"AAAAAGAAAAAGAAAA"[..], &b"AAAAAG"[..], 16, 6),
            (&b"AAAAAGAAAAAGAAAAA"[..], &b"AAAAAG"[..], 17, 6),
        ] {
            assert_eq!(
                stratum_of(&ssr_locus(bases, motif, Vec::new())),
                LocusStratum::WithoutWholeRepeatReference {
                    tract_len: Bp(tract_len),
                    period: SsrPeriod::try_new(period).expect("a period inside the STR scope"),
                }
            );
        }
    }

    /// **At period 1 every length divides**, so a mononucleotide tract can never be reported as
    /// a non-whole number of copies — which matters because mononucleotides are the bulk of what
    /// this path fits: the copy floors admit them at 8 copies where a dinucleotide needs 6.
    ///
    /// What this cannot check, and what the every-period test above exists for: at period 1 a
    /// tract's length and its repeat count are the same number, so nothing here separates
    /// dividing by the period from not dividing at all.
    #[test]
    fn a_mononucleotide_tract_is_always_a_whole_number_of_copies() {
        for length in 1..40usize {
            let locus = ssr_locus(&b"A".repeat(length), b"A", Vec::new());

            assert_eq!(
                stratum_of(&locus),
                LocusStratum::Stratified(stratum(1, u32::try_from(length).expect("small"))),
                "a run of {length} A's"
            );
        }
    }

    /// A locus of any other kind is passed over in silence rather than counted: the SNP/indel
    /// path has its own accumulator, and a repeat *cluster* is not one delimited tract, so
    /// neither is a fault to report.
    #[test]
    fn a_locus_that_is_not_one_repeat_tract_has_no_stratum_and_is_not_a_fault() {
        let mut generic = ssr_locus(b"ATATATATAT", b"AT", Vec::new());
        generic.kind = LocusKind::Generic;
        assert_eq!(stratum_of(&generic), LocusStratum::NotOneRepeatTract);

        let mut bundle = ssr_locus(b"ATATATATAT", b"AT", Vec::new());
        bundle.kind = LocusKind::SsrBundle;
        assert_eq!(stratum_of(&bundle), LocusStratum::NotOneRepeatTract);
    }

    /// A tract of no bases at all is not something the catalog emits — its floors admit nothing
    /// below three copies — but the arithmetic has to answer for it, and zero copies is what
    /// zero bases is. Pinned so that the answer is a decision rather than an accident.
    ///
    /// It also pins **where the length comes from**: zero is the only length at which this
    /// fixture's region and its reference bases disagree, the helper giving every locus a region
    /// at least one base long, so a rule reading `region.len()` instead of the tract answers one
    /// copy here and agrees everywhere else.
    #[test]
    fn a_reference_tract_of_no_bases_is_a_stratum_of_no_copies() {
        let locus = ssr_locus(b"", b"AT", Vec::new());

        assert_eq!(stratum_of(&locus), LocusStratum::Stratified(stratum(2, 0)));
    }

    // -----------------------------------------------------------------
    // One read group's reads, counted into the offset buckets.
    // -----------------------------------------------------------------

    /// One observation of `num_obs` reads that reached one border of the tract and not the
    /// other, so the length they show is a lower bound.
    fn partial_observation(bases: &[u8], num_obs: u32) -> SequenceObservation {
        let covered = u16::try_from(bases.len()).expect("a short fixture");
        SequenceObservation {
            read_witness: ReadWitness::from_left(covered, LocusLen::from_positions(200))
                .expect("a run of at least one position"),
            ..observation(bases, num_obs)
        }
    }

    /// The same observation, from another library.
    fn from_group(group: u32, observation: SequenceObservation) -> SequenceObservation {
        SequenceObservation {
            read_group: ReadGroupId(group),
            ..observation
        }
    }

    /// The tally of a locus with `reads` reads at `offset` whole copies from the reference,
    /// written the way a test reads it.
    fn tally_at(offsets: &[(i32, u32)]) -> [u32; OFFSET_BUCKETS] {
        let mut buckets = [0; OFFSET_BUCKETS];
        for &(offset, reads) in offsets {
            // Clamped rather than cast, so a fixture may name an offset far past the recorded
            // range — which is the point of two of the tests below — without the helper itself
            // deciding anything: every offset past the range lands in the end that absorbs it.
            let clamped = offset.clamp(i32::from(i8::MIN), i32::from(i8::MAX)) as i8;
            buckets[bucket_of(WholeRepeatOffset(clamped)).index()] += reads;
        }
        buckets
    }

    /// A read's offset is how far its own tract sits from the reference's, in whole copies:
    /// sixteen bases against a twenty-base dinucleotide reference is two copies short, and
    /// twenty-four is two copies long.
    #[test]
    fn a_read_is_counted_at_the_whole_copies_between_its_length_and_the_references() {
        let locus = ssr_locus(
            b"ATATATATATATATATATAT",
            b"AT",
            vec![
                observation(b"ATATATATATATATATATAT", 5),
                observation(b"ATATATATATATATAT", 3),
                observation(b"ATATATATATATATATATATATAT", 2),
            ],
        );

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(
            tally.reads_by_bucket(),
            tally_at(&[(0, 5), (-2, 3), (2, 2)])
        );
        assert_eq!(tally.reads_not_whole_repeat(), 0);
        assert_eq!(tally.depth(), 10);
    }

    /// **A read whose length is not a whole number of copies from the reference goes to the
    /// guard**, never to the nearest bucket: a read one base short of the reference has not
    /// lost a repeat, it has lost a base, and charging that to slippage is what the guard
    /// exists to prevent.
    #[test]
    fn a_read_at_a_length_that_is_not_whole_copies_from_the_reference_goes_to_the_guard() {
        let locus = ssr_locus(
            b"ATATATATATATATATATAT",
            b"AT",
            vec![
                observation(b"ATATATATATATATATATAT", 4),
                observation(b"ATATATATATATATATATA", 3),
                observation(b"ATATATATATATATATATATA", 2),
            ],
        );

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(tally.reads_by_bucket(), tally_at(&[(0, 4)]));
        assert_eq!(
            tally.reads_not_whole_repeat(),
            5,
            "one base short and one base long are both non-whole-repeat"
        );
        assert_eq!(tally.depth(), 9);
    }

    /// **The oracle for this step, and it is proven to bite.** A partial witness saw only part
    /// of the tract, so the length it shows is a lower bound; scoring it as a length reads as a
    /// read that lost repeats, which biases the direction split — the one parameter this path
    /// exists to protect.
    ///
    /// The fixture is a locus where every partial witness shows a *short* tract. Counting only
    /// the complete witnesses puts every read at the origin. Counting the partials too — the
    /// same locus, scored without the guard — puts eight reads two and three copies short, so
    /// the test would fail loudly rather than pass on a fixture that could not tell.
    #[test]
    fn a_partial_witness_is_left_out_rather_than_counted_as_a_short_tract() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 6),
                partial_observation(b"ATATATATATATATAT", 5),
                partial_observation(b"ATATATATATATAT", 3),
            ],
        );

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(
            tally.reads_by_bucket(),
            tally_at(&[(0, 6)]),
            "every counted read is at the reference length"
        );
        assert_eq!(tally.depth(), 6, "the partial witnesses are not a depth");
        assert_eq!(
            tally.reads_with_partial_witness(),
            8,
            "they are counted, because leaving eight reads out is a size a run must see"
        );

        // The same locus with those reads scored as lengths — what the guard prevents. Eight of
        // the fourteen reads would sit two and three copies short, so this test cannot pass by
        // accident on a fixture where the two rules agree.
        let scored_as_lengths = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 6),
                observation(b"ATATATATATATATAT", 5),
                observation(b"ATATATATATATAT", 3),
            ],
        );
        assert_eq!(
            tally_of(&scored_as_lengths, ReadGroupId(0)).reads_by_bucket(),
            tally_at(&[(0, 6), (-2, 5), (-3, 3)])
        );
    }

    /// **Reads that covered the tract and witnessed nothing are not a depth either**, and they
    /// are not this function's to count: they carry no length at all, not even a lower bound.
    #[test]
    fn reads_that_witnessed_nothing_do_not_enter_the_tally() {
        let reference = b"ATATATATATATATATATAT";
        let mut locus = ssr_locus(reference, b"AT", vec![observation(reference, 4)]);
        locus.reads_without_observation = 11;

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(tally.depth(), 4);
        assert_eq!(tally.reads_with_partial_witness(), 0);
    }

    /// **One read group at a time**, because slippage is a property of the chemistry: a locus
    /// covered by two libraries makes one tally in each, and a tally that pooled them would fit
    /// one stutter model to two.
    #[test]
    fn only_the_named_read_groups_reads_are_counted() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 4),
                from_group(1, observation(b"ATATATATATATATAT", 6)),
                from_group(1, partial_observation(reference, 2)),
            ],
        );

        let first = tally_of(&locus, ReadGroupId(0));
        assert_eq!(first.reads_by_bucket(), tally_at(&[(0, 4)]));
        assert_eq!(
            first.reads_with_partial_witness(),
            0,
            "another library's partial witnesses are not this library's"
        );

        let second = tally_of(&locus, ReadGroupId(1));
        assert_eq!(second.reads_by_bucket(), tally_at(&[(-2, 6)]));
        assert_eq!(second.reads_with_partial_witness(), 2);

        let absent = tally_of(&locus, ReadGroupId(7));
        assert_eq!(absent, OffsetTally::default());
    }

    /// **A read far past the recorded range is kept, in the end bucket that absorbs it** — a
    /// locus carrying a long deletion is a real observation, and dropping it would take the
    /// reads of exactly the loci whose alleles are furthest from the reference.
    #[test]
    fn reads_far_beyond_the_recorded_range_land_in_the_saturating_end_buckets() {
        let reference = b"A".repeat(200);
        let locus = ssr_locus(
            &reference,
            b"A",
            vec![
                observation(b"A", 3),
                observation(&b"A".repeat(199), 2),
                observation(&b"A".repeat(500), 4),
            ],
        );

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(
            tally.reads_by_bucket(),
            tally_at(&[(-199, 3), (-1, 2), (300, 4)]),
            "199 copies short saturates into the low end, 300 long into the high end"
        );
        assert_eq!(tally.depth(), 9);
    }

    /// **The offsets are measured in the locus's own motif period**, which is the same one
    /// `stratum_of` files it under — so a read cannot be counted in one unit and filed in
    /// another. The same four-base shortfall is two copies at a dinucleotide and four at a
    /// mononucleotide, so a rule reading the period from anywhere else answers differently.
    #[test]
    fn the_offsets_are_measured_in_the_locus_own_motif_period() {
        let dinucleotide = ssr_locus(
            b"ATATATATATATATATATAT",
            b"AT",
            vec![observation(b"ATATATATATATATAT", 5)],
        );
        let mononucleotide = ssr_locus(
            &b"A".repeat(20),
            b"A",
            vec![observation(&b"A".repeat(16), 5)],
        );

        assert_eq!(
            tally_of(&dinucleotide, ReadGroupId(0)).reads_by_bucket(),
            tally_at(&[(-2, 5)])
        );
        assert_eq!(
            tally_of(&mononucleotide, ReadGroupId(0)).reads_by_bucket(),
            tally_at(&[(-4, 5)]),
            "four bases short is four copies at period 1 and two at period 2"
        );
    }

    /// **A bucket holds every observation that lands in it, not the last one.** Two reads of the
    /// same length are two observations whenever their bases differ — an interruption, or a
    /// substitution inside the tract — and the two saturating ends pool many offsets by design,
    /// so adding rather than replacing is the whole of what a bucket does. Every other fixture
    /// here puts at most one observation in a bucket, which is why this one exists.
    #[test]
    fn a_bucket_holding_several_observations_counts_all_of_them() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                // Both one copy short, at eighteen bases, differing in their last base.
                observation(b"ATATATATATATATATAT", 4),
                observation(b"ATATATATATATATATAC", 3),
                // Both far past the low end, so both in the bucket that absorbs it.
                observation(b"ATATAT", 6),
                observation(b"AT", 2),
            ],
        );

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(
            tally.reads_by_bucket(),
            tally_at(&[(-1, 7), (-7, 6), (-9, 2)])
        );
        assert_eq!(tally.depth(), 15);
    }

    /// **A locus whose depth cap fired is still tallied**, at the depth it was observed at. The
    /// generator's reservoir is a uniform subsample, so a capped locus is the same locus seen
    /// less deeply — and skipping it would drop the deep loci, which is the depth-dependent
    /// selection this whole step exists to remove.
    #[test]
    fn a_locus_whose_depth_cap_fired_is_still_tallied() {
        let reference = b"ATATATATATATATATATAT";
        let mut locus = ssr_locus(reference, b"AT", vec![observation(reference, 4)]);
        locus.reads_discarded_by_cap = 37;

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(tally.reads_by_bucket(), tally_at(&[(0, 4)]));
        assert_eq!(tally.depth(), 4);
    }

    /// **A read that shows no tract bases at all is the tract deleted, not a read to drop.** It
    /// arrives as a complete witness — the generator slices the read's own tract span, which is
    /// empty when a read spans both borders across a full deletion — and it belongs in the low
    /// end bucket, which exists for alleles exactly this far from the reference.
    #[test]
    fn a_read_showing_no_tract_bases_is_counted_at_the_far_low_end() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(reference, b"AT", vec![observation(b"", 3)]);

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(tally.reads_by_bucket(), tally_at(&[(-10, 3)]));
        assert_eq!(tally.depth(), 3);
    }

    /// A reference tract of no bases is a stratum of no copies rather than a rejected locus, so
    /// this function is reachable with one and a read's offset is then simply its own length
    /// over the period. Pinned so that the two halves of this file answer for the same loci.
    #[test]
    fn reads_at_a_reference_tract_of_no_bases_are_measured_from_zero() {
        let locus = ssr_locus(b"", b"AT", vec![observation(b"ATAT", 3)]);

        let tally = tally_of(&locus, ReadGroupId(0));

        assert_eq!(tally.reads_by_bucket(), tally_at(&[(2, 3)]));
        assert_eq!(tally.reads_not_whole_repeat(), 0);
    }

    /// A locus that is not one repeat tract has no motif to measure offsets in, so it tallies
    /// as empty rather than as a locus whose reads all sat at the reference length.
    #[test]
    fn a_locus_that_is_not_one_repeat_tract_tallies_as_empty() {
        let mut generic = ssr_locus(
            b"ATATATATATATATATATAT",
            b"AT",
            vec![observation(b"ATATATATATATATAT", 5)],
        );
        generic.kind = LocusKind::Generic;

        assert_eq!(tally_of(&generic, ReadGroupId(0)), OffsetTally::default());
    }

    // -----------------------------------------------------------------
    // The read cap: a uniform subsample, seeded from the locus's position.
    // -----------------------------------------------------------------

    /// The same locus, moved to another position — the span goes with it, so the position is the
    /// only thing that changes. (Assigning the start first and deriving the end from it afterwards
    /// leaves the end where it was and varies the *span* instead, which is a different fixture and
    /// one a span-seeded draw would pass.)
    fn at_position(start: u64, mut locus: SampleLocusObservations) -> SampleLocusObservations {
        let span = locus.region.end.0 - locus.region.start.0;
        locus.region.start = Position(start);
        locus.region.end = Position(start + span);
        locus
    }

    /// The shape one read group's reads give a locus, cap and all.
    fn shape_at(locus: &SampleLocusObservations, read_group: ReadGroupId) -> Option<EnteredShape> {
        shape_of(tally_of(locus, read_group), locus.region)
    }

    /// A locus of `short` reads one copy short and `at_reference` reads at the reference length.
    fn deep_locus(at_reference: u32, short: u32) -> SampleLocusObservations {
        ssr_locus(
            b"ATATATATATATATATATAT",
            b"AT",
            vec![
                observation(b"ATATATATATATATATATAT", at_reference),
                observation(b"ATATATATATATATATAT", short),
            ],
        )
    }

    /// A locus no deeper than the cap enters every read it showed, and says it was not thinned.
    #[test]
    fn a_locus_no_deeper_than_the_cap_enters_every_read_it_showed() {
        let entered = shape_at(&deep_locus(8, 4), ReadGroupId(0)).expect("twelve reads");

        assert_eq!(entered.subsampled_from(), None);
        assert_eq!(entered.shape().depth(), MAX_LOCUS_READS);
        assert_eq!(entered.shape().reads_in(bucket_of(WholeRepeatOffset(0))), 8);
        assert_eq!(
            entered.shape().reads_in(bucket_of(WholeRepeatOffset(-1))),
            4
        );
    }

    /// **A deeper locus is thinned to the cap and says what it was thinned from** — the count a
    /// run reads to tell "the data was this shallow" from "the cap set the depth".
    #[test]
    fn a_locus_deeper_than_the_cap_is_thinned_to_it_and_says_so() {
        let entered = shape_at(&deep_locus(200, 100), ReadGroupId(0)).expect("three hundred reads");

        assert_eq!(entered.subsampled_from(), Some(300));
        assert_eq!(
            entered.shape().depth(),
            MAX_LOCUS_READS,
            "exactly the cap survives, never more and never fewer"
        );
    }

    /// **A locus with no reads of this group enters nothing**, rather than an empty shape: a
    /// locus with no evidence would still count towards `MIN_LOCI_TO_FIT` and towards every
    /// "loci behind this fit" number while saying nothing about any parameter.
    #[test]
    fn a_locus_this_read_group_did_not_witness_enters_nothing() {
        assert_eq!(shape_at(&deep_locus(8, 4), ReadGroupId(3)), None);

        let reference = b"ATATATATATATATATATAT";
        let only_partial = ssr_locus(
            reference,
            b"AT",
            vec![partial_observation(b"ATATATATATATATAT", 9)],
        );
        assert_eq!(
            shape_at(&only_partial, ReadGroupId(0)),
            None,
            "nine partial witnesses are nine lower bounds, which is no evidence about a length"
        );
    }

    /// **The draw is a function of the locus's position and of nothing else.** This is the
    /// property the whole sharded design rests on: a region-sharded walk and a single-threaded
    /// one meet the same locus at the same position, so they keep the same reads and their
    /// tables merge as an equality. Seeded from a counter or a thread, the two would differ by
    /// a few reads at each deep locus, which nothing short of comparing two whole walks would
    /// ever show.
    ///
    /// Checked over two hundred positions: every locus draws the same shape twice, and the
    /// shapes across positions are not all one shape — so the seed is genuinely being used.
    #[test]
    fn the_draw_is_the_same_every_time_and_depends_on_the_position() {
        let mut shapes = std::collections::BTreeSet::new();

        for start in 0..200u64 {
            let locus = at_position(start, deep_locus(200, 100));

            let once = shape_at(&locus, ReadGroupId(0)).expect("three hundred reads");
            let again = shape_at(&locus, ReadGroupId(0)).expect("three hundred reads");

            assert_eq!(
                once, again,
                "the locus at {start} drew two different shapes"
            );
            shapes.insert(once.shape());
        }

        assert!(
            shapes.len() >= 5,
            "200 positions drew only {} distinct shapes, so the position is barely reaching \
             the draw",
            shapes.len()
        );
    }

    /// **What the cap keeps is a uniform subsample, and this is what makes that measurable**:
    /// over many loci the reads kept from one bucket are hypergeometric — mean
    /// `cap · k / depth`, variance `cap · (k/depth) · (1 − k/depth) · (depth − cap) / (depth − 1)`.
    /// A draw that favoured the buckets it visited first, or that drew each bucket
    /// independently, would match the mean and miss the variance.
    #[test]
    fn the_kept_reads_are_hypergeometric_in_mean_and_variance() {
        let (at_reference, short) = (200u32, 100u32);
        let depth = f64::from(at_reference + short);
        let cap = f64::from(MAX_LOCUS_READS);
        let share = f64::from(short) / depth;

        let kept: Vec<f64> = (0..100_000u64)
            .map(|start| {
                let locus = at_position(start * 7, deep_locus(at_reference, short));
                let entered = shape_at(&locus, ReadGroupId(0)).expect("three hundred reads");
                f64::from(entered.shape().reads_in(bucket_of(WholeRepeatOffset(-1))))
            })
            .collect();

        let draws = kept.len() as f64;
        let mean = kept.iter().sum::<f64>() / draws;
        let variance = kept.iter().map(|k| (k - mean).powi(2)).sum::<f64>() / draws;

        let expected_mean = cap * share;
        let expected_variance = cap * share * (1.0 - share) * (depth - cap) / (depth - 1.0);

        assert!(
            (mean - expected_mean).abs() < 0.06,
            "mean {mean} against {expected_mean}"
        );
        assert!(
            (variance - expected_variance).abs() < 0.05,
            "variance {variance} against {expected_variance}; a draw made **with** replacement \
             would sit at {}, and the whole difference between the two models is {}",
            cap * share * (1.0 - share),
            cap * share * (1.0 - share) - expected_variance
        );
    }

    /// **The draw's population is the reads that measured a length, not every read at the
    /// locus.** A read that reached only one border of the tract is a lower bound rather than a
    /// length, so it is not in the tally and must not be in the population either: drawing
    /// twelve from a population the walk never walks would enter this locus at five reads
    /// instead of twelve — silently, because a short shape is a legal shape.
    #[test]
    fn a_capped_locus_enters_the_full_cap_however_many_partial_witnesses_it_carries() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 200),
                observation(b"ATATATATATATATATAT", 100),
                partial_observation(b"ATATATATATATATAT", 500),
            ],
        );

        let entered = shape_at(&locus, ReadGroupId(0)).expect("three hundred complete reads");

        assert_eq!(entered.shape().depth(), MAX_LOCUS_READS);
        assert_eq!(
            entered.subsampled_from(),
            Some(300),
            "the five hundred lower bounds are not part of the depth the cap thinned"
        );
    }

    /// **A locus no deeper than the cap keeps its guard reads too.** The shallow path copies the
    /// tally across rather than drawing, and dropping the guard there would bias the guard share
    /// towards zero on the majority of loci — the diagnostic that decides whether this noise
    /// model describes a stratum at all — without anything failing.
    #[test]
    fn a_locus_no_deeper_than_the_cap_keeps_its_guard_reads() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 6),
                // A base short of the reference: off the whole-repeat grid, so the guard's.
                observation(b"ATATATATATATATATATA", 4),
            ],
        );

        let entered = shape_at(&locus, ReadGroupId(0)).expect("ten reads");

        assert_eq!(entered.subsampled_from(), None);
        assert_eq!(entered.shape().reads_not_whole_repeat(), 4);
        assert_eq!(entered.shape().depth(), 10);
    }

    /// The two ends of the cap's range: one read deeper than the cap is the first depth at which
    /// the walk runs at all, and a single read is the shallowest locus there is.
    #[test]
    fn the_cap_thins_a_locus_one_read_deeper_than_it_and_leaves_a_single_read_alone() {
        let one_over = shape_at(&deep_locus(9, 4), ReadGroupId(0)).expect("thirteen reads");
        assert_eq!(one_over.subsampled_from(), Some(13));
        assert_eq!(one_over.shape().depth(), MAX_LOCUS_READS);

        let alone = shape_at(&deep_locus(1, 0), ReadGroupId(0)).expect("one read");
        assert_eq!(alone.subsampled_from(), None);
        assert_eq!(alone.shape().depth(), 1);
    }

    /// **The draw is a format, not only a relation.** `seed_at`'s own note says it may not
    /// change silently with a dependency bump or a compiler version, because it decides which
    /// reads a fit sees — and every other test here asserts a *relation* (reproducible,
    /// position-dependent, sums to the cap, hypergeometric in aggregate), all of which survive
    /// reordering the walk or moving where its state advances. These are the shapes the current
    /// algorithm draws at four positions; a diff to them is a change to that format, and a
    /// decision rather than a refactor.
    #[test]
    fn the_draw_at_a_known_position_is_the_recorded_one() {
        for (start, at_reference, short) in [
            (1_000u64, 7u32, 5u32),
            (1_001, 7, 5),
            (5_000, 11, 1),
            (900_000, 7, 5),
        ] {
            let entered = shape_at(&at_position(start, deep_locus(200, 100)), ReadGroupId(0))
                .expect("three hundred reads");

            assert_eq!(
                (
                    entered.shape().reads_in(bucket_of(WholeRepeatOffset(0))),
                    entered.shape().reads_in(bucket_of(WholeRepeatOffset(-1))),
                ),
                (at_reference, short),
                "the draw at {start} changed: that is a format change, not a refactor"
            );
        }
    }

    /// **The guard's reads are in the draw like any others.** Leaving them out of it — capping
    /// only the reads that landed in a bucket — would make the guard share rise with depth at
    /// every locus the cap fires on, and the guard share is the diagnostic that decides whether
    /// this noise model describes a stratum at all.
    #[test]
    fn the_guard_reads_are_thinned_with_the_rest() {
        let reference = b"ATATATATATATATATATAT";
        let locus = ssr_locus(
            reference,
            b"AT",
            vec![
                observation(reference, 150),
                // A base short of the reference: off the whole-repeat grid, so the guard's.
                observation(b"ATATATATATATATATATA", 150),
            ],
        );

        let guard_reads: u32 = (0..400u64)
            .map(|start| {
                let entered = shape_at(&at_position(start * 11, locus.clone()), ReadGroupId(0))
                    .expect("three hundred reads");
                assert_eq!(entered.shape().depth(), MAX_LOCUS_READS);
                entered.shape().reads_not_whole_repeat()
            })
            .sum();

        let mean = f64::from(guard_reads) / 400.0;
        assert!(
            (mean - 6.0).abs() < 0.2,
            "half the reads are the guard's, so half the cap should be too: {mean} against 6"
        );
    }

    proptest::proptest! {
        /// **Every complete read of the named group is somewhere in the tally, and nothing else
        /// is** — the conservation law the tables above can only sample. Whatever the lengths,
        /// the period and the mix of libraries, the depth is the sum of `num_obs` over that
        /// group's complete observations and the partial count is the sum over its partial ones.
        /// A bucket that replaced rather than added, or a length that fell through unscored,
        /// breaks the first equality without any table having to think of the case.
        #[test]
        fn every_complete_read_of_the_named_group_is_somewhere_in_the_tally(
            period_bases in 1usize..=MAX_MOTIF_LEN,
            reference_copies in 0usize..12,
            reads in proptest::collection::vec(
                (0usize..40, 0u32..3, 1u32..20, proptest::bool::ANY),
                0..12,
            ),
        ) {
            let motif = &b"ACGTAC"[..period_bases];
            let reference: Vec<u8> = motif
                .iter()
                .copied()
                .cycle()
                .take(period_bases * reference_copies)
                .collect();
            let observations: Vec<SequenceObservation> = reads
                .iter()
                .map(|&(shown_bases, group, num_obs, complete)| {
                    let bases: Vec<u8> =
                        motif.iter().copied().cycle().take(shown_bases).collect();
                    // A partial witness needs at least one witnessed position, so a read
                    // showing nothing goes down the complete arm — which is what it would be
                    // on real data anyway: a read spanning both borders across a full deletion.
                    let witnessed = if complete || bases.is_empty() {
                        observation(&bases, num_obs)
                    } else {
                        partial_observation(&bases, num_obs)
                    };
                    from_group(group, witnessed)
                })
                .collect();
            let locus = ssr_locus(&reference, motif, observations);

            let tally = tally_of(&locus, ReadGroupId(1));

            let counted: u32 = locus
                .observations
                .iter()
                .filter(|obs| {
                    obs.read_group == ReadGroupId(1) && obs.read_witness == ReadWitness::Complete
                })
                .map(|obs| obs.num_obs)
                .sum();
            let partial: u32 = locus
                .observations
                .iter()
                .filter(|obs| {
                    obs.read_group == ReadGroupId(1) && obs.read_witness != ReadWitness::Complete
                })
                .map(|obs| obs.num_obs)
                .sum();

            proptest::prop_assert_eq!(tally.depth(), counted);
            proptest::prop_assert_eq!(tally.reads_with_partial_witness(), partial);
        }

        /// **A tract is stratified exactly when its period divides its length, and its count is
        /// then the quotient** — stated as a law over every period the scope admits and every
        /// remainder each of them can leave, so no period and no remainder goes untried. The
        /// tables above are the readable pins; this is the one that cannot be satisfied by a
        /// rule that happens to fit them.
        #[test]
        fn a_tract_is_stratified_exactly_when_its_period_divides_its_length(
            period in 1usize..=MAX_MOTIF_LEN,
            tract_bases in 0usize..200,
        ) {
            let motif = &b"ACGTAC"[..period];
            let bases: Vec<u8> = motif.iter().copied().cycle().take(tract_bases).collect();

            let expected = if tract_bases.is_multiple_of(period) {
                LocusStratum::Stratified(stratum(
                    u8::try_from(period).expect("a period inside the STR scope"),
                    u32::try_from(tract_bases / period).expect("under 200 bases"),
                ))
            } else {
                LocusStratum::WithoutWholeRepeatReference {
                    tract_len: Bp(tract_bases as u64),
                    period: SsrPeriod::try_new(period).expect("a period inside the STR scope"),
                }
            };

            proptest::prop_assert_eq!(stratum_of(&ssr_locus(&bases, motif, Vec::new())), expected);
        }
    }
}
