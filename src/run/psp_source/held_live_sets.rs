//! **The reads live at every record a sample's source still holds, stored as what changed.**
//!
//! A held record's body cannot be decoded without the set of reads live at its position: one
//! observation's read list is stored as that set minus every other observation's
//! ([`KeptRecord`](super::KeptRecord)). The source used to keep a full copy of the set for every
//! record it held, which costs *records held × reads spanning each record*. That is nothing at
//! three reads a position and it is what killed the 2,169-sample tomato run at a pile-up of
//! 5,669 reads: one sample's source held 20 million identifiers for 8,816 records
//! (`doc/devel/reports/psp_pileup_memory_2026-09-28/report.md`).
//!
//! # What is kept instead
//!
//! **Per record, only the reads that arrived and departed there** — which the walk has already
//! decoded from the record's head, because that is how the file itself stores them. **Every so
//! often, a whole copy of the set**, called a restart point here. A set is rebuilt by starting
//! from the nearest restart point at or before the record and applying each later record's
//! changes in turn.
//!
//! **A restart point is written as the changes from an empty set**: nothing departed, and the
//! whole set arrived. So a replay has one rule and no special case — empty the set, then apply
//! every entry from the restart point to the record — and a psp block's own restart, where the
//! file's first record of a block restates the whole set as arrivals, is simply a restart point
//! the file chose.
//!
//! # When a restart point is taken
//!
//! **When the changes stored since the last one add up to as many identifiers as the set now
//! holds.** That puts one copy of the set per set's worth of changes, so memory is at most about
//! twice the changes, and a rebuild reads at most about twice the set — whatever the depth. At
//! three reads a position with 150-base reads a restart point lands every hundred records or so;
//! at 5,669 reads with about 76 changes a position, every 75. Nothing is tuned to either.
//! [`LEAST_CHANGES_BETWEEN_RESTARTS`] keeps a restart point from being taken at every record
//! where almost nothing is live, where it would cost more in bookkeeping than it saves.

use crate::psp::{LiveSet, LiveSetChanges};
use crate::types::ChainId;

/// Fewest changes stored between two restart points, whatever the set's size. At a gap in
/// coverage the set is empty and any change would otherwise trigger a new restart point.
const LEAST_CHANGES_BETWEEN_RESTARTS: usize = 16;

/// One held record's entry: where its identifiers sit in the arena, as the reads that departed
/// followed by the reads that arrived.
#[derive(Debug, Clone, Copy)]
struct Entry {
    /// First identifier, counted from the file's first record — the same file-wide addressing
    /// the source's body arena uses, so an entry outlives the release that moves the arena.
    start: u64,
    departed: u32,
    arrived: u32,
}

/// The live sets of one sample's held records. See the module documentation.
#[derive(Debug, Default)]
pub(super) struct HeldLiveSets {
    /// One per record from the oldest restart point still needed, in draw order.
    entries: Vec<Entry>,
    /// The ordinal (counted from the file's first record) of `entries[0]`.
    entries_released: u64,
    /// Every entry's identifiers, back to back.
    ids: Vec<ChainId>,
    /// How many identifiers have been dropped from the front of `ids`.
    ids_released: u64,
    /// Ordinals of the records that are restart points, ascending.
    restarts: Vec<u64>,
    /// Identifiers stored as changes since the last restart point.
    changes_since_restart: usize,
    /// How many blocks the walk had begun at the record last pushed, to see a block restart.
    blocks_begun: u64,
}

impl HeldLiveSets {
    /// The ordinal the next pushed record will get.
    fn next_ordinal(&self) -> u64 {
        self.entries_released + self.entries.len() as u64
    }

    /// Record one drawn record's live set and return its ordinal, which is what
    /// [`rebuild`](Self::rebuild) takes.
    ///
    /// `live` is the set as of this record, `changes` the step to it from the record before, and
    /// `blocks_begun` the walk's count of blocks opened — a change in it means `changes` is a
    /// restatement from an empty set rather than a step.
    pub(super) fn push(
        &mut self,
        live: &LiveSet,
        changes: &LiveSetChanges,
        blocks_begun: u64,
    ) -> u64 {
        let ordinal = self.next_ordinal();
        let block_restarted = blocks_begun != self.blocks_begun;
        self.blocks_begun = blocks_begun;
        let changed = changes.departed().len() + changes.arrived().len();
        let restart = block_restarted
            || self.entries.is_empty()
            || self.changes_since_restart + changed
                >= live.len().max(LEAST_CHANGES_BETWEEN_RESTARTS);
        let start = self.ids_released + self.ids.len() as u64;
        let entry = if restart {
            self.ids.extend_from_slice(live.ids());
            self.restarts.push(ordinal);
            self.changes_since_restart = 0;
            Entry {
                start,
                departed: 0,
                arrived: u32::try_from(live.len()).expect("fewer than 2³² reads live at once"),
            }
        } else {
            self.ids.extend_from_slice(changes.departed());
            self.ids.extend_from_slice(changes.arrived());
            self.changes_since_restart += changed;
            Entry {
                start,
                departed: changes.departed().len() as u32,
                arrived: changes.arrived().len() as u32,
            }
        };
        self.entries.push(entry);
        ordinal
    }

    /// Rebuild into `set` the reads live at the record `ordinal`, which must still be held.
    ///
    /// # Panics
    ///
    /// If `ordinal` was released or never pushed — the merge only builds inside its own window,
    /// so this is a defect in the release bookkeeping, not anything about a file.
    pub(super) fn rebuild(&self, ordinal: u64, set: &mut LiveSet) {
        let restart = self.restarts[self.restarts.partition_point(|at| *at <= ordinal) - 1];
        let first = usize::try_from(restart - self.entries_released).expect("an index");
        let last = usize::try_from(ordinal - self.entries_released).expect("an index");
        set.clear();
        for entry in &self.entries[first..=last] {
            let at = usize::try_from(entry.start - self.ids_released).expect("an index");
            let departed_end = at + entry.departed as usize;
            let arrived_end = departed_end + entry.arrived as usize;
            set.apply(
                &self.ids[at..departed_end],
                &self.ids[departed_end..arrived_end],
            );
        }
    }

    /// Drop everything no record at or after `first_held` needs: entries and identifiers before
    /// the last restart point at or before it. With `first_held` past every pushed record, drop
    /// everything.
    pub(super) fn release_before(&mut self, first_held: u64) {
        let keep_from = if first_held >= self.next_ordinal() {
            self.next_ordinal()
        } else {
            let restarts_before = self.restarts.partition_point(|at| *at <= first_held);
            self.restarts[restarts_before - 1]
        };
        if keep_from <= self.entries_released {
            return;
        }
        let entries = usize::try_from(keep_from - self.entries_released).expect("an index");
        let ids_from = self
            .entries
            .get(entries)
            .map_or(self.ids_released + self.ids.len() as u64, |entry| {
                entry.start
            });
        self.entries.drain(..entries);
        self.entries_released = keep_from;
        self.ids
            .drain(..usize::try_from(ids_from - self.ids_released).expect("an index"));
        self.ids_released = ids_from;
        let restarts = self.restarts.partition_point(|at| *at < keep_from);
        self.restarts.drain(..restarts);
    }

    /// How many read identifiers are held, changes and restart points together.
    pub(super) fn held_ids(&self) -> usize {
        self.ids.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A walk's view of a sequence of live sets: each set, and the changes that led to it.
    fn changes_between(before: &[ChainId], after: &[ChainId]) -> (Vec<ChainId>, Vec<ChainId>) {
        let departed = before
            .iter()
            .copied()
            .filter(|id| after.binary_search(id).is_err())
            .collect();
        let arrived = after
            .iter()
            .copied()
            .filter(|id| before.binary_search(id).is_err())
            .collect();
        (departed, arrived)
    }

    /// Sets that slide along like reads over a genome: at record `r`, reads `r .. r + depth`
    /// are live.
    fn sliding_sets(records: u64, depth: u64) -> Vec<Vec<ChainId>> {
        (0..records).map(|r| (r..r + depth).collect()).collect()
    }

    /// **Every held record rebuilds to exactly the set it was pushed with**, across block
    /// restarts, at a shallow and a deep corner, and after releases that leave the oldest
    /// held record between two restart points.
    #[test]
    fn every_held_record_rebuilds_to_the_set_it_was_pushed_with() {
        for depth in [0u64, 3, 40, 3_000] {
            let sets = sliding_sets(2_000, depth);
            let mut held = HeldLiveSets::default();
            let mut blocks = 1;
            let mut previous: Vec<ChainId> = Vec::new();
            let mut first_held = 0u64;
            for (r, now) in sets.iter().enumerate() {
                // A block restarts every 700 records: the file restates the whole set.
                if r > 0 && r % 700 == 0 {
                    blocks += 1;
                    previous.clear();
                }
                let (departed, arrived) = changes_between(&previous, now);
                let changes = LiveSetChanges::for_tests(departed, arrived);
                let live = LiveSet::from_sorted_slice(now);
                assert_eq!(held.push(&live, &changes, blocks), r as u64);
                previous = now.clone();
                // Release in uneven steps, as the merge's rounds do.
                if r % 333 == 332 {
                    first_held = r as u64 - 50;
                    held.release_before(first_held);
                }
            }
            let mut rebuilt = LiveSet::new();
            for ordinal in first_held..sets.len() as u64 {
                held.rebuild(ordinal, &mut rebuilt);
                assert_eq!(
                    rebuilt.ids(),
                    &sets[ordinal as usize][..],
                    "depth {depth}, record {ordinal}"
                );
            }
            held.release_before(sets.len() as u64);
            assert_eq!(held.held_ids(), 0, "a full release holds nothing");
        }
    }

    /// **At depth the arena holds a small fraction of a full copy per record.** 3,000 reads
    /// live and one arriving and one departing per record: 3,000 identifiers a record the old
    /// way, and about four here — two changes, plus a restart point's 3,000 spread over the
    /// 1,500 records between restart points.
    #[test]
    fn at_depth_the_arena_holds_the_changes_not_a_copy_per_record() {
        let depth = 3_000;
        let records = 6_000u64;
        let sets = sliding_sets(records, depth);
        let mut held = HeldLiveSets::default();
        let mut previous: Vec<ChainId> = Vec::new();
        for now in &sets {
            let (departed, arrived) = changes_between(&previous, now);
            held.push(
                &LiveSet::from_sorted_slice(now),
                &LiveSetChanges::for_tests(departed, arrived),
                1,
            );
            previous = now.clone();
        }
        let per_record = held.held_ids() as f64 / records as f64;
        assert!(per_record < 6.0, "{per_record} identifiers a record");
    }
}
