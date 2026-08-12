//! Thinning a locus's reads down to a cap, the same way twice.
//!
//! Both paths of step 4 cap how many of a locus's reads enter their evidence, and both must
//! thin the same locus to the same reads however the genome was cut into shards — otherwise
//! two runs of the same sample disagree, and merging a sharded walk stops being an equality.
//! What makes that hold is where the randomness comes from: **the locus's position, and
//! nothing else**.
//!
//! **A uniform subsample is exact rather than approximate.** Thinning a locus's reads
//! uniformly leaves the counts distributed exactly as they would have been at the lower
//! depth, so what a cap costs is precision and never a bias — which is the whole reason the
//! design caps by subsampling instead of by dropping deep loci, that being depth-dependent
//! selection.
//!
//! This module was lifted out of the SNP/indel path when the STR path's read cap needed the
//! same draw. Nothing here knows what a read is: it walks a population, keeps a fixed-size
//! uniformly random subset of it, and can be **resumed**, which is what lets one draw be
//! split across the categories a caller counts its population in.

use crate::ng::types::GenomeRegion;

/// The random stream a locus's subsample is drawn from: **a function of where the locus is,
/// and of nothing else**.
///
/// That is what makes a region-sharded walk and a single-threaded one keep the same reads, so
/// merging stays exact and a fitted rate does not move with the thread count. Seeding from a
/// shared stream, or from anything carrying the walk's history, would make the evidence at a
/// capped locus a fact about the scheduling rather than about the data — the same fault the
/// pileup's own read sampler was rewritten to remove.
///
/// **Spelled out here rather than taken from a hashing crate**, for that sampler's reason:
/// this decides which reads a fit sees, so it is a format, and it may not change silently with
/// a dependency bump or a compiler version. `ahash` and `DefaultHasher` are both free to seed
/// themselves per process.
pub(crate) fn seed_at(region: GenomeRegion) -> u64 {
    // splitmix64's mixing constants. The contig and the start are folded separately so
    // that position 1 of contig 2 and position 2 of contig 1 are different streams.
    let mut seed = u64::from(region.contig.0);
    seed = splitmix64(seed ^ 0x9e37_79b9_7f4a_7c15);
    seed = splitmix64(seed ^ region.start.0);
    seed
}

/// splitmix64's finaliser — a bijection, so distinct inputs stay distinct.
pub(crate) fn splitmix64(mut state: u64) -> u64 {
    state = (state ^ (state >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    state = (state ^ (state >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    state ^ (state >> 31)
}

/// A partly-walked selection sample: how much of the population is still to be passed over,
/// and how many places in the kept set are still open.
///
/// **Resumable, which is the point.** A caller counts its population in categories — one
/// site's alternative reads by library, one locus's reads by offset bucket — and walks the
/// categories one at a time through a single walk, so the per-category kept counts sum to
/// exactly what one call over all of them would have returned. That holds by construction
/// rather than by two implementations agreeing.
pub(crate) struct SelectionWalk {
    state: u64,
    remaining_population: u64,
    remaining_draws: u64,
}

impl SelectionWalk {
    pub(crate) fn new(seed: u64, population: u32, draws: u32) -> Self {
        Self {
            state: seed,
            remaining_population: u64::from(population),
            remaining_draws: u64::from(draws),
        }
    }

    /// Walk the next `members` items of the population and return how many were kept.
    ///
    /// **Selection sampling.** Taking a uniformly random subset of size `draws` is the same
    /// as walking the population and keeping each item with probability
    /// `remaining_draws / remaining_population`, decrementing both; the items are
    /// exchangeable, so walking one category's members and stopping there gives that
    /// category's kept count without touching the rest.
    ///
    /// **The caller's obligation, and the compiler cannot state it:** the members summed over
    /// every call to one walk may not exceed the population it was built with. Walk further and
    /// `remaining_population` wraps — the release profile leaves `overflow-checks` off — after
    /// which every later category keeps nothing, silently. The debug assertion below is what
    /// makes a test catch it; a release build cannot afford the check per category and does not
    /// need it, since the callers walk a population they counted themselves.
    pub(crate) fn keep_from(&mut self, members: u32) -> u32 {
        debug_assert!(
            u64::from(members) <= self.remaining_population,
            "a selection walk was asked for {members} more members than the {} left in its \
             population",
            self.remaining_population
        );
        let mut kept = 0u32;
        for done in 0..members {
            if self.remaining_draws == 0 {
                // The kept set is full; the rest of the population is passed over.
                self.remaining_population -= u64::from(members - done);
                return kept;
            }
            if self.remaining_draws == self.remaining_population {
                // Every item still in the running is kept, so every member still to be
                // walked is kept. **The count is the remaining members, not the remaining
                // draws** — those are equal only when every item left is a member of this
                // category, and taking the draws would have inflated the category's count at
                // exactly the deep, category-rich loci where the cap fires.
                let left = members - done;
                self.remaining_draws -= u64::from(left);
                self.remaining_population -= u64::from(left);
                return kept + left;
            }
            self.state = splitmix64(self.state);
            // The modulo's bias is `2^64 mod remaining_population` out of `2^64` — at most
            // one part in 10^15 for any depth a pileup can produce, which is far below the
            // sampling noise the draw is made of.
            if self.state % self.remaining_population < self.remaining_draws {
                kept += 1;
                self.remaining_draws -= 1;
            }
            self.remaining_population -= 1;
        }
        kept
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    fn region(contig: u32, start: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(start),
        }
    }

    /// **`SelectionWalk`'s one contract, over the whole domain rather than at a fixture.**
    /// Walking the population in pieces must give exactly what walking it at once gives —
    /// that is what a caller rests on when it splits a locus's reads between categories, and
    /// the two fast paths are where a resumed walk's state can be left wrong without any
    /// *single* call looking odd.
    ///
    /// It is here because fixtures could not reach it: deleting the population decrement in
    /// the all-kept arm — leaving every later category drawing against a population one too
    /// high — left the whole suite green. That arm fires exactly where the design says the
    /// damage lands, at the deep loci where the cap bites, and nothing would have panicked.
    #[test]
    fn a_resumed_walk_sums_to_a_single_call_over_every_split() {
        for population in [1u32, 2, 3, 5, 8, 124, 125, 200, 500] {
            for draws in [0u32, 1, 2, 3, population / 2, population - 1, population] {
                if draws > population {
                    continue;
                }
                for members in 0..=population.min(40) {
                    for seed in [1u64, 7, 0xdead_beef, u64::MAX] {
                        let whole = SelectionWalk::new(seed, population, draws).keep_from(members);
                        assert!(whole <= draws.min(members));
                        for cut in 0..=members {
                            let mut walk = SelectionWalk::new(seed, population, draws);
                            let first = walk.keep_from(cut);
                            let rest = walk.keep_from(members - cut);
                            assert_eq!(
                                first + rest,
                                whole,
                                "population {population}, draws {draws}, \
                                 members {members}, cut after {cut}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// **The seed is a function of the position and of nothing else**, so the same locus is
    /// thinned to the same reads in every run and in every shard layout — and two different
    /// loci are different streams, including the pair a fold in one order would collide:
    /// position 1 of contig 2 against position 2 of contig 1.
    #[test]
    fn the_seed_is_the_position_and_neighbouring_positions_are_different_streams() {
        assert_eq!(seed_at(region(3, 1_000)), seed_at(region(3, 1_000)));

        let streams = [
            seed_at(region(1, 2)),
            seed_at(region(2, 1)),
            seed_at(region(1, 1)),
            seed_at(region(1, 3)),
            seed_at(region(0, 0)),
        ];
        for (at, seed) in streams.iter().enumerate() {
            for other in &streams[at + 1..] {
                assert_ne!(seed, other, "two loci share a random stream");
            }
        }
    }

    /// The two degenerate ends need no random number at all: a category holding none of the
    /// population keeps none, and a walk whose draws equal its population keeps everything.
    #[test]
    fn the_degenerate_ends_are_exact_rather_than_drawn() {
        assert_eq!(SelectionWalk::new(1, 500, 124).keep_from(0), 0);
        assert_eq!(SelectionWalk::new(1, 500, 500).keep_from(500), 500);
        assert_eq!(
            SelectionWalk::new(1, 124, 124).keep_from(60),
            60,
            "nothing dropped"
        );
    }
}
