//! One locus reduced to one cell key: how many reads covered the site, and how many
//! of those showed something other than the reference.
//!
//! **The only place that decides what counts as an alternative read**, which is why
//! it is its own file rather than a method on the locus type — a locus cannot know a
//! model's answer to that question
//! (`arch/parameter_prepass_generic.md` §2.3).
//!
//! The answer here is byte equality against the locus's own reference bases, and it
//! reduces to production's `allele_index == 0` rule exactly where the two can be
//! compared. Three things follow from it, and each is a decision rather than an
//! implementation detail — see [`count_whole_site`].

use smallvec::SmallVec;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::{GenomeRegion, ReadGroupId};

use super::depth_bins::DepthBinEdges;
use super::histogram::DepthAndAltReads;

/// Count one generic locus's reads, taking the site **whole** — the spec's word for the
/// un-split site (`spec/parameter_prepass_generic.md` §4). This is the pair the windowed
/// histogram enters, because a site's genotype is one thing however many read groups
/// covered it.
///
/// **What counts as an alternative read.** A complete witness whose bases differ from
/// the locus's reference bases, byte for byte. That is the per-read form of production's
/// `allele_index == 0` rule and it needs no special cases: `bases` is in *read*
/// coordinates, so an insertion makes it longer than the reference slice and a deletion
/// shorter, and either way it cannot compare equal — which is the right answer for both.
/// The reference bases arrive canonical `{A, C, G, T, N}` from
/// [`RefSeq::fetch_into`](crate::ng::ref_seq::RefSeq::fetch_into), so a soft-masked
/// lowercase reference cannot make every read at a repeat look alternative; the verbatim
/// reader that preserves soft-masking is a different capability, used by the typed-region
/// catalog and not on this path.
///
/// **Complete witnesses only.** A read that spanned part of the locus witnessed neither
/// the reference allele nor an alternative one at the positions it missed, and scoring it
/// as either is the mis-scoring `complete_observations` exists to guard. Loci one base
/// wide — the overwhelming majority on this path — lose nothing, since a read either
/// covers the base or does not.
///
/// **`reads_without_observation` does not enter the depth.** Those reads covered the
/// locus and witnessed nothing, so counting them would assert they showed the reference,
/// and they did not show anything. The depth here is reads whose evidence is present,
/// which is what every likelihood in this design conditions on. They are counted
/// separately by the accumulator instead: a locus where many reads say nothing is one
/// whose depth is not what it appears, and that is a fact about the mapping rather than
/// about the genotype. It also has to be this way for
/// `spec/parameter_prepass_generic.md` §12.6 to hold at all — the field is a bare scalar
/// with no read-group attribution, so including it would make the two histograms differ
/// by construction on a single-library sample.
///
/// **`reads_discarded_by_cap` does not skip the locus.** An earlier draft skipped any
/// locus whose reads an upstream cap had already subsampled. That is depth-dependent
/// selection of exactly the kind this step exists to remove: ng's generic pileup caps
/// indel-bearing columns at 250 reads against 8,000 for the rest, so at 300× every
/// indel-bearing column would be dropped and at 30× none would, and the coverage-invariance
/// anchor of §9 would be measuring the skip rule rather than the estimator. The locus is
/// entered at the depth observed and the accumulator counts it
/// (`AccumulationCounts::loci_with_upstream_subsample`).
///
/// **A site deeper than the ladder's cap is subsampled down to it here**, before the
/// pair is built, so the depth recorded and the depth the alternative count belongs to
/// are the same number. `edges` is passed rather than a bare depth so that the cap and
/// the ladder cannot disagree — an earlier draft of the design set the cap to 300 while
/// the ladder ended at 124, which left depths 125–300 with no bin at all. See
/// [`CountedSite`] for what the subsample is and why it is not a rescale.
///
/// # Panics
///
/// If the locus's supports sum past `u32`, which no pileup cap allows.
#[must_use]
pub fn count_whole_site(locus: &SampleLocusObservations, edges: &DepthBinEdges) -> CountedSite {
    let mut depth = 0u32;
    let mut alt_reads = 0u32;

    for observation in locus.complete_observations() {
        depth = add_support(depth, observation.num_obs, "depth");
        if *observation.bases != *locus.reference_bases {
            alt_reads = add_support(alt_reads, observation.num_obs, "alternative reads");
        }
    }

    CountedSite::capped(depth, alt_reads, edges.max_depth(), seed_at(locus.region))
}

/// The same count, split between the read groups that covered the site — one entry per
/// group, ascending by [`ReadGroupId`].
///
/// This is what the read-group histogram enters: a site with 20 reads from one library
/// and 10 from another becomes two entries, one at depth 20 and one at depth 10, because
/// an error rate describes the chemistry and two libraries of one sample can genuinely
/// differ. The split costs that fit nothing, since what it counts is reads.
///
/// **It is also where the windowed histogram's library attribution comes from**, whose
/// depth is the *sum* of these — exact rather than approximate, since these are raw
/// counts at one position (`spec/parameter_prepass_generic.md` §4).
///
/// `out` is scratch: cleared, then filled. It is an argument rather than a return value
/// because this runs once per covered position over hundreds of millions of them, and a
/// fresh `Vec` per locus would be an allocation on the hottest path in the pre-pass.
///
/// What counts as an alternative read, which witnesses are read, and what is left out of
/// the depth are all [`count_whole_site`]'s doc — the two functions differ only in grain.
///
/// **Each group's entry is capped on its own depth**, because each is a histogram entry
/// in its own right and has to have a bin. A group's reads cannot outnumber the site's,
/// so at a site inside the cap no group is touched.
///
/// # Panics
///
/// If any group's supports sum past `u32`, which no pileup cap allows.
pub fn count_by_read_group(
    locus: &SampleLocusObservations,
    edges: &DepthBinEdges,
    out: &mut Vec<(ReadGroupId, CountedSite)>,
) {
    let mut running = ReadGroupTotals::new();
    accumulate_by_read_group(locus, &mut running);

    let cap = edges.max_depth();
    let seed = seed_at(locus.region);
    out.clear();
    out.extend(running.into_iter().map(|(group, depth, alt_reads)| {
        (group, CountedSite::capped(depth, alt_reads, cap, seed))
    }));
}

/// The whole-site count **with the kept alternative reads split between the libraries
/// they came from** — what a multi-library sample's windowed entry is keyed on.
///
/// **This exists because the split has to be of the *same* subsample.** The windowed
/// table enters a site once, at its total depth, with the attribution of its alternative
/// reads; that attribution must sum to exactly the alternative count of the entry it
/// belongs to, and [`add_attributed_site`](super::histogram::DepthAltHistogram::add_attributed_site)
/// refuses it otherwise. Above
/// the cap the per-group counts of [`count_by_read_group`] cannot serve: each of those is
/// its own subsample of its own group's reads, drawn against its own depth, so their
/// alternative counts sum to something else entirely. The two tables ask different
/// questions of the same site — one about the chemistry of each library, one about the
/// genotype of the individual — and above the cap they are answered by two different
/// draws, which is right. It is only their *internal* consistency that matters, and this
/// is what gives the windowed entry its own.
///
/// `alt_by_group` is scratch: cleared, then filled in read-group order, and it holds only
/// groups that kept at least one alternative read. Below the cap it is the exact per-group
/// alternative counts, with no draw at all.
///
/// # Panics
///
/// If the locus's supports sum past `u32`, which no pileup cap allows.
pub fn count_whole_site_by_library(
    locus: &SampleLocusObservations,
    edges: &DepthBinEdges,
    alt_by_group: &mut Vec<(ReadGroupId, u32)>,
) -> CountedSite {
    let mut running = ReadGroupTotals::new();
    accumulate_by_read_group(locus, &mut running);

    let depth = running.iter().fold(0u32, |total, &(_, group_depth, _)| {
        add_support(total, group_depth, "depth")
    });
    let alt_reads = running.iter().fold(0u32, |total, &(_, _, group_alt)| {
        add_support(total, group_alt, "alternative reads")
    });

    alt_by_group.clear();
    let cap = edges.max_depth();
    if depth <= cap {
        alt_by_group.extend(
            running
                .iter()
                .filter(|&&(_, _, group_alt)| group_alt > 0)
                .map(|&(group, _, group_alt)| (group, group_alt)),
        );
        return CountedSite {
            counts: DepthAndAltReads::new(depth, alt_reads),
            subsampled_from: None,
        };
    }

    // One walk over the site's alternative reads, group by group, sharing a single
    // selection state — so the totals it yields are the *same* draw
    // `CountedSite::capped` would have made, split. That is what keeps the entry's
    // attribution summing to its own alternative count.
    let mut walk = SelectionWalk::new(seed_at(locus.region), depth, cap);
    let mut kept_total = 0u32;
    for &(group, _, group_alt) in &running {
        let kept = walk.keep_from(group_alt);
        kept_total += kept;
        if kept > 0 {
            alt_by_group.push((group, kept));
        }
    }

    CountedSite {
        counts: DepthAndAltReads::new(cap, kept_total),
        subsampled_from: Some(depth),
    }
}

/// One site's `(read group, depth, alternative reads)` totals, in read-group order.
///
/// **Eight inline, not two, and the figure is measured.** The tomato archive survey has
/// 133 samples carrying two libraries, 20 carrying three, and four carrying 7, 16, 16
/// and 42 (`spec/read_groups.md`) — so "two is every sample with more than one library"
/// is false. At two inline, `add_locus` allocated twice per locus from three groups
/// upward (200 allocations per 100 steady-state loci, measured) and cost about 47 ns a
/// locus more, breaking the architecture's contract that it never allocates after a
/// key's first locus — on exactly the multi-library samples the machinery exists for.
/// Eight takes 1,703 of the 1,707 surveyed samples to zero; the four above it spill once
/// per locus rather than always.
type ReadGroupTotals = SmallVec<[(ReadGroupId, u32, u32); 8]>;

/// Total each read group's complete witnesses into `(group, depth, alternative reads)`,
/// in read-group order. Shared so that the two grains cannot come to disagree about what
/// they are counting.
fn accumulate_by_read_group(locus: &SampleLocusObservations, running: &mut ReadGroupTotals) {
    running.clear();
    for observation in locus.complete_observations() {
        let group = observation.read_group;
        let at = match running.iter().position(|&(seen, _, _)| seen == group) {
            Some(at) => at,
            None => {
                running.push((group, 0, 0));
                running.len() - 1
            }
        };
        let entry = &mut running[at];
        entry.1 = add_support(entry.1, observation.num_obs, "depth");
        if *observation.bases != *locus.reference_bases {
            entry.2 = add_support(entry.2, observation.num_obs, "alternative reads");
        }
    }
    running.sort_unstable_by_key(|&(group, _, _)| group);
}

/// One site counted, and whether the ladder's cap had to be applied to it.
///
/// **A site deeper than the cap keeps `cap` of its reads and counts the alternative ones
/// among them.** Losing the rest costs nothing: at 124 reads a heterozygote shows about
/// 62 alternative reads and a homozygous-reference site about 0.12, so the genotype is
/// already certain and more depth cannot make it so.
///
/// **Draw the kept alternative count; do not rescale it.** The two are not the same
/// thing, and the difference lands on the `k = 0, 1, 2` cells that carry every bit of
/// the error-rate evidence:
///
/// - *Rescaling and rounding to nearest* is not a subsample at all. A 200-read site with
///   one alternative read becomes `(124, 1)` — an alternative fraction of 1/124 against
///   a true 1/200, nearly double — while the same lone read at 500 reads becomes
///   `(124, 0)` and vanishes. **The bias reverses sign at depth 248**, where
///   `1 × 124/depth` stops rounding up to one, so the rule inflates rare alleles below
///   that depth and erases them above it.
/// - *Rescaling with a stochastic round* — take the floor, add one with probability
///   equal to the fraction — fixes the **mean** and breaks the **spread**. Thinning `k`
///   that way gives a variance of about `r²·Var[k]` where a real subsample of the same
///   reads gives `r·Var[k]`; at the 124-from-500 ratio that is four times too narrow,
///   and the fit reads an under-dispersed alternative-count distribution as a cleaner
///   read group than the data supports.
/// - **Subsampling the reads is exact**: keep `cap` of the site's reads and count the
///   alternative ones among them, which is a hypergeometric draw and marginally the
///   binomial the likelihood assumes. No rescaling, no correction, nothing to argue
///   about.
///
/// **This is the one depth mechanism shipping without a measurement behind it.** No
/// world in either research harness reaches a cap — at mean depth 60 the deepest site
/// the Poisson truncation keeps is 125, against a cap of 124 — so what a deeper site
/// costs is argued rather than measured (research note §4.6). It fires on samples above
/// about 124×, which is HG002 and never tomato.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct CountedSite {
    counts: DepthAndAltReads,
    subsampled_from: Option<u32>,
}

impl CountedSite {
    /// The pair the histogram is keyed on — after the cap, where it fired.
    #[inline]
    #[must_use]
    pub fn counts(self) -> DepthAndAltReads {
        self.counts
    }

    /// The depth this site was observed at, when that was deeper than the cap and the
    /// reads were subsampled down to it. `None` when the site was entered whole.
    ///
    /// Reported rather than swallowed because a run where most sites appear here is one
    /// where the ladder, not the data, is setting the depth — the accumulator tallies it
    /// as `sites_subsampled_to_cap` (Milestone C3).
    #[inline]
    #[must_use]
    pub fn subsampled_from(self) -> Option<u32> {
        self.subsampled_from
    }

    /// Count `depth` reads of which `alt_reads` showed something other than the
    /// reference, keeping at most `cap` of them.
    fn capped(depth: u32, alt_reads: u32, cap: u32, seed: u64) -> Self {
        if depth <= cap {
            return Self {
                counts: DepthAndAltReads::new(depth, alt_reads),
                subsampled_from: None,
            };
        }
        let kept_alt = hypergeometric_draw(seed, depth, alt_reads, cap);
        Self {
            counts: DepthAndAltReads::new(cap, kept_alt),
            subsampled_from: Some(depth),
        }
    }
}

/// The random stream a site's subsample is drawn from: **a function of where the site
/// is, and of nothing else**.
///
/// That is what makes a region-sharded walk and a single-threaded one keep the same
/// reads, so `merge` stays exact and a fitted rate does not move with the thread count.
/// Seeding from a shared stream, or from anything carrying the walk's history, would
/// make the evidence at a capped site a fact about the scheduling rather than about the
/// data — the same fault the pileup's own read sampler was rewritten to remove.
///
/// **Spelled out here rather than taken from a hashing crate**, for that sampler's
/// reason: this decides which reads a fit sees, so it is a format, and it may not change
/// silently with a dependency bump or a compiler version. `ahash` and `DefaultHasher`
/// are both free to seed themselves per process.
fn seed_at(region: GenomeRegion) -> u64 {
    // splitmix64's mixing constants. The contig and the start are folded separately so
    // that position 1 of contig 2 and position 2 of contig 1 are different streams.
    let mut seed = region.contig.0 as u64;
    seed = splitmix64(seed ^ 0x9e37_79b9_7f4a_7c15);
    seed = splitmix64(seed ^ region.start.0);
    seed
}

/// splitmix64's finaliser — a bijection, so distinct inputs stay distinct.
fn splitmix64(mut state: u64) -> u64 {
    state = (state ^ (state >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    state = (state ^ (state >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    state ^ (state >> 31)
}

/// How many of `successes` survive when `draws` of `population` are kept — a
/// hypergeometric draw, exactly.
///
/// **Selection sampling, walked over the successes only.** Taking a uniformly random
/// subset of size `draws` is the same as walking the population and keeping each item
/// with probability `remaining_draws / remaining_population`, decrementing both; the
/// items are exchangeable, so walking the successes first and stopping there gives their
/// kept count without touching the rest. That is what makes this cheap where it runs:
/// at a homozygous-reference site of depth 300 the expected number of alternative reads
/// is 0.3, so the loop almost always draws nothing at all, and the expensive case — a
/// homozygous non-reference site, where every read is a success — is one site in a
/// thousand.
///
/// **The draw is a pure function of `(seed, population, successes, draws)`, and that is
/// deliberate.** It means the whole-site count and a lone read group's count at the same
/// site are the *same* draw rather than two, which is what
/// `spec/parameter_prepass_generic.md` §12.6 needs: on a single-library sample the
/// read-group histogram must equal the windowed one cell for cell, and above the cap it
/// could not if the two grains drew independently. Two read groups that happen to arrive
/// at a site with identical depths and identical alternative counts do then get identical
/// draws; that correlates their entries at that one site and leaves each one's marginal
/// exactly hypergeometric, which is all any fit here reads.
fn hypergeometric_draw(seed: u64, population: u32, successes: u32, draws: u32) -> u32 {
    SelectionWalk::new(seed, population, draws).keep_from(successes)
}

/// A partly-walked selection sample: how much of the population is still to be passed
/// over, and how many places in the kept set are still open.
///
/// **Resumable, which is the point.** One site's alternative reads are walked group by
/// group through a single walk, so the per-group kept counts sum to exactly what one
/// call over all of them would have returned — that is what
/// [`count_whole_site_by_library`] rests on, and it holds by construction rather than by
/// two implementations agreeing.
struct SelectionWalk {
    state: u64,
    remaining_population: u64,
    remaining_draws: u64,
}

impl SelectionWalk {
    fn new(seed: u64, population: u32, draws: u32) -> Self {
        Self {
            state: seed,
            remaining_population: u64::from(population),
            remaining_draws: u64::from(draws),
        }
    }

    /// Walk the next `successes` items of the population and return how many were kept.
    fn keep_from(&mut self, successes: u32) -> u32 {
        let mut kept = 0u32;
        for done in 0..successes {
            if self.remaining_draws == 0 {
                // The kept set is full; the rest of the population is passed over.
                self.remaining_population -= u64::from(successes - done);
                return kept;
            }
            if self.remaining_draws == self.remaining_population {
                // Every item still in the running is kept, so every success still to be
                // walked is kept. **The count is the remaining successes, not the
                // remaining draws** — those are equal only when every item left is a
                // success, and taking the draws would have inflated the alternative
                // count at exactly the deep, allele-rich sites where the cap fires.
                let left = successes - done;
                self.remaining_draws -= u64::from(left);
                self.remaining_population -= u64::from(left);
                return kept + left;
            }
            self.state = splitmix64(self.state);
            // The modulo's bias is `2^64 mod remaining_population` out of `2^64` — at
            // most one part in 10^15 for any depth a pileup can produce, which is far
            // below the sampling noise the draw is made of.
            if self.state % self.remaining_population < self.remaining_draws {
                kept += 1;
                self.remaining_draws -= 1;
            }
            self.remaining_population -= 1;
        }
        kept
    }
}

/// Add one observation's support, refusing to wrap.
///
/// The release profile leaves `overflow-checks` off, so a wrapped depth would come back
/// as a small number and be scored as a shallow site. It cannot happen through the
/// generic pileup, which caps a column at 8,000 reads — this is the module's rule that a
/// counter says so rather than wrapping quietly.
fn add_support(running: u32, support: u32, what: &str) -> u32 {
    running.checked_add(support).unwrap_or_else(|| {
        panic!("a locus's {what} passed u32 at {running}, with {support} more to add")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, LocusLen, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, Position};

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    /// One observation, with only the fields this step reads set to anything meaningful.
    fn observation(
        bases: &[u8],
        witness: ReadWitness,
        read_group: ReadGroupId,
        num_obs: u32,
    ) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: witness,
            read_group,
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus(
        reference_bases: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        let span = reference_bases.len().max(1) as u64;
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(span),
            },
            reference_bases: reference_bases.into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    fn ladder() -> DepthBinEdges {
        DepthBinEdges::new()
    }

    fn by_group(locus: &SampleLocusObservations) -> Vec<(ReadGroupId, CountedSite)> {
        let mut out = Vec::new();
        count_by_read_group(locus, &ladder(), &mut out);
        out
    }

    /// The pair a site is entered at, for the tests that predate the cap and are not
    /// about it.
    fn whole(locus: &SampleLocusObservations) -> DepthAndAltReads {
        count_whole_site(locus, &ladder()).counts()
    }

    /// **The one-base site, which is the overwhelming majority of this path.** Eleven
    /// reads showed the reference base and three showed something else, so the site is
    /// entered at depth 14 with 3 alternative reads.
    #[test]
    fn a_one_base_site_counts_the_reads_that_disagreed_with_the_reference() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 11),
                observation(b"C", ReadWitness::Complete, group(0), 3),
            ],
        );

        let counted = whole(&site);
        assert_eq!(counted.depth(), 14);
        assert_eq!(counted.alt_reads(), 3);
    }

    /// **An insertion and a deletion are alternative reads without either being a
    /// special case.** `bases` is in read coordinates, so an inserted base makes the
    /// observation longer than the reference slice and a deletion shorter; neither can
    /// compare equal, which is the right answer for both. A rule that compared only
    /// same-length sequences, or that compared position by position, would need to name
    /// them — and would then have to decide what a length mismatch means.
    #[test]
    fn an_insertion_and_a_deletion_are_both_alternative_reads() {
        let widened = locus(
            b"ACGT",
            vec![
                observation(b"ACGT", ReadWitness::Complete, group(0), 20),
                observation(b"ACCGT", ReadWitness::Complete, group(0), 4), // an insertion
                observation(b"AGT", ReadWitness::Complete, group(0), 6),   // a deletion
            ],
        );

        let counted = whole(&widened);
        assert_eq!(counted.depth(), 30);
        assert_eq!(counted.alt_reads(), 10);
    }

    /// **A read that witnessed part of the locus witnessed neither allele at the
    /// positions it missed**, so it enters neither the depth nor the alternative count.
    /// Scoring it as the reference would assert it agreed where it saw nothing; scoring
    /// it as alternative would assert the opposite. Both are the mis-scoring
    /// `complete_observations` exists to guard.
    #[test]
    fn a_partial_witness_enters_neither_the_depth_nor_the_alternative_count() {
        // Two of the four locus positions, witnessed from the left — a read that ran out
        // before the locus's right border.
        let saw_two_of_four =
            ReadWitness::from_left(2, LocusLen::from_positions(4)).expect("a non-empty run");
        let widened = locus(
            b"ACGT",
            vec![
                observation(b"ACGT", ReadWitness::Complete, group(0), 9),
                observation(b"AC", saw_two_of_four.clone(), group(0), 5),
                observation(b"AG", saw_two_of_four, group(0), 7),
            ],
        );

        let counted = whole(&widened);
        assert_eq!(counted.depth(), 9, "only the complete witness");
        assert_eq!(counted.alt_reads(), 0);
    }

    /// A locus every read covered without witnessing anything is entered at depth zero,
    /// not skipped and not counted as reference. `reads_without_observation` is a bare
    /// scalar with no read-group attribution, so it could not be split between the two
    /// histograms even if it were wanted — including it would break the cell-for-cell
    /// equality of `spec/parameter_prepass_generic.md` §12.6 by construction.
    #[test]
    fn reads_that_witnessed_nothing_do_not_enter_the_depth() {
        let mut silent = locus(b"A", Vec::new());
        silent.reads_without_observation = 30;

        let counted = whole(&silent);
        assert_eq!(counted.depth(), 0);
        assert_eq!(counted.alt_reads(), 0);
        assert!(by_group(&silent).is_empty(), "no group covered it");
    }

    /// A locus whose reads an upstream cap already subsampled is entered **at the depth
    /// observed**, not skipped. Skipping it would be depth-dependent selection: ng's
    /// pileup caps indel-bearing columns at 250 reads and everything else at 8,000, so
    /// at 300× every indel-bearing column would go and at 30× none would.
    #[test]
    fn a_locus_an_upstream_cap_subsampled_is_entered_at_the_depth_observed() {
        // Kept below this step's own cap of 124, so that the two subsamplings do not
        // tangle: what is being tested is that an *upstream* cap does not make the
        // locus disappear, and `CountedSite` owns what our own cap does.
        let mut capped = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 96),
                observation(b"T", ReadWitness::Complete, group(0), 4),
            ],
        );
        capped.reads_discarded_by_cap = 1_400;

        let counted = count_whole_site(&capped, &ladder());
        assert_eq!(
            counted.counts().depth(),
            100,
            "the reads that are here, not the 1,500"
        );
        assert_eq!(counted.counts().alt_reads(), 4);
        assert_eq!(
            counted.subsampled_from(),
            None,
            "this step did not subsample it — an upstream cap already had"
        );
    }

    /// **The read-group split, and the entries come back in read-group order** whatever
    /// order the observations arrived in — the histogram keys on this, and an unsorted
    /// listing would split one site's cell in two.
    #[test]
    fn the_read_group_split_comes_back_in_read_group_order() {
        let site = locus(
            b"A",
            vec![
                observation(b"T", ReadWitness::Complete, group(7), 2),
                observation(b"A", ReadWitness::Complete, group(3), 18),
                observation(b"A", ReadWitness::Complete, group(7), 8),
                observation(b"T", ReadWitness::Complete, group(3), 2),
            ],
        );

        let split = by_group(&site);
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].0, group(3));
        assert_eq!(
            (split[0].1.counts().depth(), split[0].1.counts().alt_reads()),
            (20, 2)
        );
        assert_eq!(split[1].0, group(7));
        assert_eq!(
            (split[1].1.counts().depth(), split[1].1.counts().alt_reads()),
            (10, 2)
        );
    }

    /// **The two grains agree, which is the property the two histograms rest on.** The
    /// windowed table enters a site once at its total depth and the read-group table
    /// once per group; if the split did not sum to the whole, the two tables would
    /// describe different data and §12.6's cell-for-cell equality on a single-library
    /// sample could not hold.
    ///
    /// **Every rule the two functions share has to appear in at least one of these
    /// loci, or this cannot see them drifting apart.** Making `count_whole_site` alone
    /// read the partial witnesses left it green while the two grains disagreed by
    /// twelve reads, because no fixture had one — so the last locus carries a partial
    /// witness, a widened span and two libraries at once.
    #[test]
    fn the_read_group_split_sums_to_the_whole_site() {
        let saw_two_of_four =
            ReadWitness::from_left(2, LocusLen::from_positions(4)).expect("a non-empty run");
        for site in [
            locus(
                b"A",
                vec![
                    observation(b"A", ReadWitness::Complete, group(0), 11),
                    observation(b"C", ReadWitness::Complete, group(1), 3),
                    observation(b"A", ReadWitness::Complete, group(1), 7),
                    observation(b"G", ReadWitness::Complete, group(4), 1),
                ],
            ),
            locus(
                b"A",
                vec![observation(b"A", ReadWitness::Complete, group(2), 5)],
            ),
            locus(b"ACGT", Vec::new()),
            locus(
                b"ACGT",
                vec![
                    observation(b"ACGT", ReadWitness::Complete, group(0), 9),
                    observation(b"ACCGT", ReadWitness::Complete, group(1), 4),
                    observation(b"AC", saw_two_of_four.clone(), group(0), 5),
                    observation(b"AG", saw_two_of_four.clone(), group(1), 7),
                ],
            ),
        ] {
            let whole_site = whole(&site);
            let split = by_group(&site);

            let depth: u32 = split
                .iter()
                .map(|(_, counted)| counted.counts().depth())
                .sum();
            let alt_reads: u32 = split
                .iter()
                .map(|(_, counted)| counted.counts().alt_reads())
                .sum();
            assert_eq!(depth, whole_site.depth());
            assert_eq!(alt_reads, whole_site.alt_reads());
        }
    }

    /// A single-library site produces exactly one entry, which is the case that matters
    /// most: 1,550 of the 1,707 samples in the tomato archive survey carry one library,
    /// and for them the read-group table and the windowed one hold the same numbers.
    #[test]
    fn a_single_library_site_produces_one_entry_equal_to_the_whole() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(6), 11),
                observation(b"C", ReadWitness::Complete, group(6), 3),
            ],
        );

        let split = by_group(&site);
        assert_eq!(split.len(), 1);
        assert_eq!(split[0].0, group(6));
        assert_eq!(split[0].1.counts(), whole(&site));
    }

    /// The scratch vector is cleared before it is filled, so a caller reusing one across
    /// loci does not accumulate the previous locus's groups.
    #[test]
    fn the_scratch_vector_is_cleared_between_loci() {
        let mut out = Vec::new();
        let first = locus(
            b"A",
            vec![observation(b"C", ReadWitness::Complete, group(9), 4)],
        );
        count_by_read_group(&first, &ladder(), &mut out);
        assert_eq!(out.len(), 1);

        let second = locus(
            b"A",
            vec![observation(b"A", ReadWitness::Complete, group(1), 6)],
        );
        count_by_read_group(&second, &ladder(), &mut out);
        assert_eq!(out.len(), 1, "the previous locus's group is gone");
        assert_eq!(out[0].0, group(1));
        assert_eq!(out[0].1.counts().alt_reads(), 0);
    }

    // ----------------------------------------------------------------------
    // The depth cap (C2). A site deeper than the ladder's top keeps 124 of its
    // reads and counts the alternative ones among them.
    // ----------------------------------------------------------------------

    /// A site inside the cap is entered whole and says so, which is every site of every
    /// sample below about 124× — tomato's whole cohort.
    #[test]
    fn a_site_inside_the_cap_is_entered_whole() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 120),
                observation(b"T", ReadWitness::Complete, group(0), 4),
            ],
        );

        let counted = count_whole_site(&site, &ladder());
        assert_eq!(counted.counts().depth(), 124, "the ladder's top, exactly");
        assert_eq!(counted.counts().alt_reads(), 4);
        assert_eq!(counted.subsampled_from(), None);
    }

    /// A site deeper than the cap is entered **at the cap**, and reports the depth it
    /// was observed at so the accumulator can tally how much of a run rests on it.
    #[test]
    fn a_site_deeper_than_the_cap_is_entered_at_the_cap_and_reports_its_true_depth() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 400),
                observation(b"T", ReadWitness::Complete, group(0), 100),
            ],
        );

        let counted = count_whole_site(&site, &ladder());
        assert_eq!(counted.counts().depth(), 124);
        assert_eq!(counted.subsampled_from(), Some(500));
        assert!(counted.counts().alt_reads() <= 124);
    }

    /// **The same locus position gives the same draw on every run**, which is what makes
    /// a region-sharded walk and a single-threaded one keep the same reads — and so what
    /// makes `merge` exact and a fitted rate independent of the thread count. A
    /// *different* position draws differently, or the stream would not be a stream.
    #[test]
    fn the_draw_is_a_function_of_the_locus_position_and_nothing_else() {
        let deep = |contig: u32, start: u64| {
            let mut site = locus(
                b"A",
                vec![
                    observation(b"A", ReadWitness::Complete, group(0), 250),
                    observation(b"T", ReadWitness::Complete, group(0), 250),
                ],
            );
            site.region.contig = ContigId(contig);
            site.region.start = Position(start);
            site.region.end = Position(start);
            count_whole_site(&site, &ladder()).counts().alt_reads()
        };

        assert_eq!(deep(0, 1_000), deep(0, 1_000), "same position, same draw");
        // Contig and start are folded separately, so position 1 of contig 2 and
        // position 2 of contig 1 are different streams.
        let across_the_genome: Vec<u32> = (0..8).map(|at| deep(0, 1_000 + at)).collect();
        assert!(
            across_the_genome.windows(2).any(|pair| pair[0] != pair[1]),
            "neighbouring positions draw independently: {across_the_genome:?}"
        );
        assert_ne!(
            (1..40).map(|c| deep(c, 7)).collect::<Vec<_>>(),
            vec![deep(0, 7); 39],
            "the contig enters the seed"
        );
    }

    /// **The oracle: the kept alternative count is hypergeometric.** Over many positions
    /// — each an independent draw, since the seed is the position — the mean and
    /// variance match the closed form for drawing `draws` of `population` without
    /// replacement from `successes`:
    ///
    /// ```text
    /// E[k]   = draws · successes / population
    /// Var[k] = draws · p · (1 − p) · (population − draws) / (population − 1),  p = successes/population
    /// ```
    ///
    /// **The variance is what separates this from the two shortcuts, and it separates
    /// them by two orders of magnitude.** Both were written and run against this test:
    /// rescaling and rounding to nearest gives every site the same answer, so variance
    /// **0**; rescaling with a stochastic round gives variance `f·(1−f)` where `f` is
    /// the fraction dropped — here about **0.19** — against a true **23.4**. The mean is
    /// right in all three, to well inside this test's tolerance, so asserting the mean
    /// alone would have accepted either shortcut. A fit reads that missing spread as a
    /// cleaner read group than the data supports.
    ///
    /// The population is deliberately **not** a whole multiple of the draw ratio —
    /// 251 of 500 kept 124 gives 62.248 — so that the stochastic round is genuinely
    /// stochastic here rather than degenerating to the deterministic one.
    #[test]
    fn the_kept_alternative_count_is_hypergeometric_in_mean_and_variance() {
        let (population, successes, draws) = (500u32, 251u32, 124u32);
        let positions = 4_000u32;

        let kept: Vec<f64> = (0..positions)
            .map(|at| {
                let seed = seed_at(GenomeRegion {
                    contig: ContigId(0),
                    start: Position(u64::from(at) * 7 + 1),
                    end: Position(u64::from(at) * 7 + 1),
                });
                f64::from(hypergeometric_draw(seed, population, successes, draws))
            })
            .collect();

        let n = f64::from(positions);
        let mean = kept.iter().sum::<f64>() / n;
        let variance = kept.iter().map(|k| (k - mean).powi(2)).sum::<f64>() / (n - 1.0);

        let p = f64::from(successes) / f64::from(population);
        let expected_mean = f64::from(draws) * p;
        let expected_variance = f64::from(draws) * p * (1.0 - p) * f64::from(population - draws)
            / f64::from(population - 1);

        // The standard error of the mean over 4,000 draws is about 0.07, and of the
        // variance about 0.5; three of each is a tolerance that a real fault clears by
        // far and sampling noise does not.
        assert!(
            (mean - expected_mean).abs() < 0.21,
            "mean {mean} against {expected_mean}"
        );
        assert!(
            (variance - expected_variance).abs() < 1.6,
            "variance {variance} against {expected_variance}"
        );
    }

    /// **`SelectionWalk`'s one contract, over the whole domain rather than at a
    /// fixture.** Walking the successes in pieces must give exactly what walking them at
    /// once gives — that is what `count_whole_site_by_library` rests on when it splits a
    /// site's alternative reads between libraries, and the two fast paths are where a
    /// resumed walk's state can be left wrong without any *single* call looking odd.
    ///
    /// It is here because the fixtures could not reach it: deleting the population
    /// decrement in the all-kept arm — leaving every later group drawing against a
    /// population one too high — left the whole suite green. That arm fires exactly
    /// where the design says the damage lands, at the deep allele-rich sites where the
    /// cap bites, and nothing would have panicked: the entry's attribution still sums to
    /// its own total, so `add_attributed_site`'s cross-check passes.
    #[test]
    fn a_resumed_walk_sums_to_a_single_call_over_every_split() {
        for population in [1u32, 2, 3, 5, 8, 124, 125, 200, 500] {
            for draws in [0u32, 1, 2, 3, population / 2, population - 1, population] {
                if draws > population {
                    continue;
                }
                for successes in 0..=population.min(40) {
                    for seed in [1u64, 7, 0xdead_beef, u64::MAX] {
                        let whole =
                            SelectionWalk::new(seed, population, draws).keep_from(successes);
                        assert!(whole <= draws.min(successes));
                        for cut in 0..=successes {
                            let mut walk = SelectionWalk::new(seed, population, draws);
                            let first = walk.keep_from(cut);
                            let rest = walk.keep_from(successes - cut);
                            assert_eq!(
                                first + rest,
                                whole,
                                "population {population}, draws {draws}, \
                                 successes {successes}, cut after {cut}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// The two degenerate ends are exact rather than drawn: a site with no alternative
    /// read keeps none however deep it is, and one where every read is alternative keeps
    /// the cap. The first is the overwhelming majority of a genome, and it costs no
    /// random number at all.
    #[test]
    fn a_site_with_no_alternative_reads_or_only_alternative_reads_needs_no_draw() {
        assert_eq!(hypergeometric_draw(1, 500, 0, 124), 0);
        assert_eq!(hypergeometric_draw(1, 500, 500, 124), 124);
        assert_eq!(hypergeometric_draw(1, 124, 60, 124), 60, "nothing dropped");
    }

    /// **A lone read group's entry is the same draw as the whole site's**, which is what
    /// `spec/parameter_prepass_generic.md` §12.6 needs above the cap: on a
    /// single-library sample the read-group histogram must equal the windowed one cell
    /// for cell, and it could not if the two grains drew independently. Below the cap
    /// this is trivially true; the test is here for the deep samples where it is not.
    #[test]
    fn a_single_library_deep_site_draws_once_for_both_grains() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(3), 380),
                observation(b"T", ReadWitness::Complete, group(3), 120),
            ],
        );

        let whole_site = count_whole_site(&site, &ladder());
        let split = by_group(&site);

        assert_eq!(split.len(), 1);
        assert_eq!(split[0].1, whole_site, "one site, one draw");
        assert_eq!(whole_site.subsampled_from(), Some(500));
    }

    /// **The windowed entry's attribution sums to its own alternative count — at every
    /// depth, capped or not.** This is the property `add_attributed_site` asserts, and
    /// the reason the split cannot be taken from `count_by_read_group`: above the cap
    /// those are each a subsample of their own group's reads against their own depth, so
    /// their alternative counts sum to something else. Checked here across the cap, at
    /// 40, 124, 125 and 500 reads.
    #[test]
    fn the_windowed_entrys_attribution_sums_to_its_own_alternative_count() {
        let mut alt_by_group = Vec::new();
        for depth in [40u32, 124, 125, 500] {
            let from_first = depth / 2;
            let site = locus(
                b"A",
                vec![
                    observation(b"A", ReadWitness::Complete, group(0), from_first - 3),
                    observation(b"T", ReadWitness::Complete, group(0), 3),
                    observation(
                        b"A",
                        ReadWitness::Complete,
                        group(5),
                        depth - from_first - 7,
                    ),
                    observation(b"T", ReadWitness::Complete, group(5), 7),
                ],
            );

            let entry = count_whole_site_by_library(&site, &ladder(), &mut alt_by_group);
            let attributed: u32 = alt_by_group.iter().map(|&(_, kept)| kept).sum();
            assert_eq!(
                attributed,
                entry.counts().alt_reads(),
                "at depth {depth}: {alt_by_group:?}"
            );
            assert!(alt_by_group.iter().all(|&(_, kept)| kept > 0));
            assert!(
                alt_by_group.windows(2).all(|pair| pair[0].0 < pair[1].0),
                "read-group order"
            );
        }
    }

    /// And the whole-site pair it reports is the very pair [`count_whole_site`] reports,
    /// so the two ways in cannot key one site to two cells.
    #[test]
    fn splitting_the_attribution_does_not_move_the_whole_site_pair() {
        let mut alt_by_group = Vec::new();
        for depth in [40u32, 500] {
            let site = locus(
                b"A",
                vec![
                    observation(b"A", ReadWitness::Complete, group(0), depth - 30),
                    observation(b"T", ReadWitness::Complete, group(0), 20),
                    observation(b"T", ReadWitness::Complete, group(5), 10),
                ],
            );

            assert_eq!(
                count_whole_site_by_library(&site, &ladder(), &mut alt_by_group),
                count_whole_site(&site, &ladder()),
                "at depth {depth}"
            );
        }
    }

    /// **The reference bases arrive canonical, and this is what would break if they did
    /// not.** ng has two reference readers: `RefSeq::fetch_into` returns canonical
    /// `{A, C, G, T, N}`, and the raw reader preserves soft-masked lowercase for the
    /// typed-region catalog's byte oracle. A locus built from the second would call
    /// every read at a soft-masked position alternative — about half the human genome,
    /// and the repeat-rich half. This pins the comparison as byte equality so the
    /// dependency is visible rather than implied.
    #[test]
    fn byte_equality_is_the_rule_and_a_lowercase_reference_would_break_it() {
        let canonical = locus(
            b"A",
            vec![observation(b"A", ReadWitness::Complete, group(0), 12)],
        );
        assert_eq!(whole(&canonical).alt_reads(), 0);

        let soft_masked = locus(
            b"a",
            vec![observation(b"A", ReadWitness::Complete, group(0), 12)],
        );
        assert_eq!(
            whole(&soft_masked).alt_reads(),
            12,
            "every read would read as alternative — which is why the generic locus \
             generator fetches through RefSeq and not the raw reader"
        );
    }
}
