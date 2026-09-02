//! Retired locus records, kept so the next one can be filled into them — **G1.**
//!
//! # What it is for
//!
//! The merge is the last owner of every record a walk emits: it draws one, walks it,
//! evicts it, and lets it go. Before this pool existed the walk then allocated the next
//! record's buffers from scratch, so a run over 8 Mb of reference at 63 accessions
//! allocated and freed the same four heap blocks a million times a sample — the
//! reference bases, an observation's bases, its chain-id list, and the observation list
//! itself — to hold, at any instant, one record.
//!
//! The hook for handing a record back has been in
//! [`ObservationSource::next_observation`](crate::ng::run::ObservationSource::next_observation)
//! since Milestone B; nothing filled it. This is what fills it.
//!
//! # Why it is bounded, and why the bound is asserted rather than documented
//!
//! **A pool that never refuses a record is a leak wearing a pool's clothes.** B1's suite
//! pins that a returned record does not come back out as an observation; it cannot pin
//! that the record is *released*, because a walker that simply drops it has nothing to
//! count. Measured at the time: a walker that stashed every offered record for ever
//! passed all fourteen of those tests — the one survivor of a twenty-one-mutation pass
//! that killed the other twenty. At 63 samples that is unbounded growth of exactly the
//! records this pool exists to stop allocating.
//!
//! So the pool takes a hard bound, refuses past it, and
//! [`RecordPool::kept`](RecordPool::kept) makes the count observable to a test.
//!
//! # Why the bound is two
//!
//! One walk holds at most one record in flight at a time: the walk mints a record, the
//! merge draws it, walks it, evicts it and offers it back, and only then does the walk
//! mint another. A second slot absorbs the one round of overlap the merge's window
//! allows without ever growing. Anything above that would be holding buffers against a
//! demand that does not exist.

use crate::ng::locus_generation::SampleLocusObservations;

/// How many retired records one walk keeps. See the module note.
pub(super) const RECORDS_KEPT: usize = 2;

/// Retired locus records, waiting to be filled again.
///
/// **Not `Clone`, and deliberately.** Two pools holding clones of one record would
/// double the buffers this exists to keep singular, and nothing needs to copy a pool.
#[derive(Debug, Default)]
pub(super) struct RecordPool {
    kept: Vec<SampleLocusObservations>,
}

impl RecordPool {
    /// An empty pool.
    pub(super) fn new() -> Self {
        Self {
            kept: Vec::with_capacity(RECORDS_KEPT),
        }
    }

    /// Offer a finished record back. **Refused, and dropped, once the pool is full.**
    ///
    /// Refusing rather than growing is the whole safety property: see the module note on
    /// what an unbounded pool passed.
    pub(super) fn put(&mut self, record: SampleLocusObservations) {
        if self.kept.len() < RECORDS_KEPT {
            self.kept.push(record);
        }
    }

    /// A record to fill. One that was handed back if there is one, a fresh one otherwise.
    ///
    /// **What comes back is not empty** — it holds the previous locus's values in buffers
    /// of the previous locus's size, and the caller overwrites every field. That is the
    /// point, and it is why the fill sites build an exhaustive struct literal out of the
    /// taken buffers rather than assigning field by field: a field added to
    /// [`SequenceObservation`](crate::ng::locus_generation::SequenceObservation) then
    /// fails to compile at the fill site instead of silently arriving stale.
    pub(super) fn take(&mut self) -> SampleLocusObservations {
        self.kept
            .pop()
            .unwrap_or_else(SampleLocusObservations::empty_shell)
    }

    /// How many records the pool is holding — **the bound, made observable.**
    ///
    /// Exists for the test the milestone owes: without a count, a pool that keeps
    /// everything and a pool that keeps two are indistinguishable from the outside.
    #[cfg(test)]
    pub(super) fn kept(&self) -> usize {
        self.kept.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The bound holds against a source that never stops offering.** This is the
    /// property B1 could not state, and the mutation that survived its whole suite was a
    /// walker that kept every record: here that shows up as a count that climbs.
    #[test]
    fn the_pool_refuses_records_past_its_bound() {
        let mut pool = RecordPool::new();
        for _ in 0..1_000 {
            pool.put(SampleLocusObservations::empty_shell());
        }
        assert_eq!(
            pool.kept(),
            RECORDS_KEPT,
            "the pool grew past its bound, which is the leak this bound exists to stop",
        );
    }

    /// **A record handed back is the record handed out**, buffers and all — otherwise the
    /// pool is an allocation with extra steps.
    #[test]
    fn a_record_comes_back_with_the_buffers_it_went_in_with() {
        let mut pool = RecordPool::new();
        let mut record = SampleLocusObservations::empty_shell();
        record.reference_bases.reserve(64);
        let capacity_lent = record.reference_bases.capacity();
        pool.put(record);
        let back = pool.take();
        assert_eq!(
            back.reference_bases.capacity(),
            capacity_lent,
            "the pool handed back a record whose buffers had been reallocated",
        );
    }

    /// **An empty pool still answers**, so no fill site needs a branch for the first
    /// locus of a walk or for a walk the merge never draws from.
    #[test]
    fn an_empty_pool_hands_out_a_fresh_record() {
        let mut pool = RecordPool::new();
        assert!(pool.take().observations.is_empty());
    }
}
