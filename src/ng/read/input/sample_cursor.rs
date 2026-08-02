//! **One sample's reads over one chromosome** — the layer a caller actually holds.
//!
//! A sample is *k* files: one for most, several when it was sequenced in more than one run.
//! Both generators ask a *sample* for reads, not a file, so the cursor a caller keeps has to
//! be a sample-level one that holds a file cursor each and merges what they yield
//! (`arch/alignment_cursor.md` §2.4).

#![allow(
    dead_code,
    reason = "the generators move onto this at Milestone D; until then it is reached from its \
              own tests. Remove at D2 — if it is still needed then, the generators are not \
              reading through the cursor that was built for them."
)]

use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::types::{ContigId, GenomePosition, Position};

use super::cursor::{AlignmentCursor, CursorError};

/// One sample's reads over one chromosome.
///
/// **An enum, and the single-file arm is the *absence* of a merge rather than a merge with
/// k = 1.** A one-file sample — which is most of them — pays no merge machinery at all. The
/// two arms are unified by an enum and not by `Box<dyn …>` because dynamic dispatch is opaque
/// to the optimiser: the whole per-read chain would stop being inlined at the boundary, on
/// the hottest loop in the module. `SampleRegionReads`, which this replaces, is an enum for
/// exactly that reason and this keeps the property.
/// **On the size difference clippy reports here.** `Single` holds a cursor inline (about a
/// kilobyte, most of it the read filter's reused buffers) while `Merged` holds a `Vec` of
/// them behind a pointer, so the two arms differ by roughly 900 bytes. Boxing `Single` would
/// even them up and put a pointer chase on every read of the common case — which is the exact
/// cost this enum exists to avoid, and a worse trade than one under-filled stack slot per
/// sample. There is one of these per sample per worker, not one per region.
#[allow(clippy::large_enum_variant, reason = "see the note above the enum")]
pub enum SampleCursor<R: RawRefSeq> {
    Single(AlignmentCursor<R>),
    Merged(MergedCursors<R>),
}

impl<R: RawRefSeq> SampleCursor<R> {
    /// Build from one cursor per file. Panics on an empty list, which is not a runtime
    /// condition: a sample with no files cannot be opened.
    pub fn new(mut files: Vec<AlignmentCursor<R>>) -> Self {
        assert!(
            !files.is_empty(),
            "a sample cursor needs at least one file cursor; a sample with no files cannot \
             be opened",
        );
        if files.len() == 1 {
            Self::Single(files.pop().expect("length checked as 1"))
        } else {
            Self::Merged(MergedCursors::new(files))
        }
    }

    /// The chromosome every file cursor here covers.
    pub fn contig(&self) -> ContigId {
        match self {
            Self::Single(cursor) => cursor.contig(),
            Self::Merged(merged) => merged.cursors[0].contig(),
        }
    }

    /// Point every file cursor at `region`.
    ///
    /// **All or nothing.** If any file refuses, the whole move fails — a sample answered from
    /// some of its files is worse than one that answers with an error, because the reads it
    /// returns are real and the ones it drops are invisible.
    pub fn move_to_region(
        &mut self,
        region: crate::ng::types::GenomeRegion,
    ) -> Result<(), CursorError> {
        match self {
            Self::Single(cursor) => cursor.move_to_region(region),
            Self::Merged(merged) => merged.move_to_region(region),
        }
    }

    /// The next read of the current region, in position order across every file.
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>> {
        match self {
            Self::Single(cursor) => cursor.next_read(),
            Self::Merged(merged) => merged.next_read(),
        }
    }
}

/// Argmin k-way merge over a sample's per-file cursors.
///
/// **Keys are held beside the heads, not read through them.** `keys[i]` mirrors `heads[i]` and
/// is refreshed only when that head is refilled, so each step scans a small contiguous array
/// of positions in cache rather than dereferencing k reads — the same layout
/// `MergedRegionReads` uses, and for the same reason.
pub struct MergedCursors<R: RawRefSeq> {
    cursors: Vec<AlignmentCursor<R>>,
    /// `None` = that cursor has no more reads for this region.
    heads: Vec<Option<AlignedRead>>,
    /// `keys[i]` is `heads[i]`'s key; `None` in lockstep.
    keys: Vec<Option<GenomePosition>>,
    /// Cleared by `move_to_region`, because a primed head belongs to a region.
    filled: bool,
    /// An error met while refilling *after* a good read was already in hand. The read is
    /// handed over first and this is yielded next — without it the read would be dropped on
    /// the floor, which is the one outcome this layer exists to prevent.
    pending_error: Option<CursorError>,
    /// Set once every cursor is exhausted or an error has been yielded.
    done: bool,
}

impl<R: RawRefSeq> MergedCursors<R> {
    fn new(cursors: Vec<AlignmentCursor<R>>) -> Self {
        let k = cursors.len();
        Self {
            cursors,
            heads: (0..k).map(|_| None).collect(),
            keys: vec![None; k],
            filled: false,
            pending_error: None,
            done: false,
        }
    }

    fn move_to_region(
        &mut self,
        region: crate::ng::types::GenomeRegion,
    ) -> Result<(), CursorError> {
        for cursor in &mut self.cursors {
            cursor.move_to_region(region)?;
        }
        // The heads belong to the region that primed them.
        self.heads.iter_mut().for_each(|head| *head = None);
        self.keys.iter_mut().for_each(|key| *key = None);
        self.filled = false;
        self.pending_error = None;
        self.done = false;
        Ok(())
    }

    /// Pull one read into slot `i`, recomputing just that key.
    fn refill(&mut self, i: usize) -> Result<(), CursorError> {
        match self.cursors[i].next_read() {
            None => {
                self.heads[i] = None;
                self.keys[i] = None;
                Ok(())
            }
            Some(Ok(read)) => {
                self.keys[i] = Some(key_of(&read));
                self.heads[i] = Some(read);
                Ok(())
            }
            Some(Err(error)) => {
                self.heads[i] = None;
                self.keys[i] = None;
                Err(error)
            }
        }
    }

    /// The lowest key, **ties to the lowest file index**.
    ///
    /// The tie-break is not incidental: a sample's files usually cover the same coordinate
    /// range, so ties are routine rather than rare, and this is what makes the output order
    /// reproducible run to run.
    fn argmin(&self) -> Option<usize> {
        let mut best: Option<usize> = None;
        for (i, key) in self.keys.iter().enumerate() {
            let Some(key) = key else { continue };
            match best {
                // Strictly less, so an equal key leaves the earlier index in place — that
                // *is* the tie-break.
                Some(b) if key >= &self.keys[b].expect("best always has a key") => {}
                _ => best = Some(i),
            }
        }
        best
    }

    fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>> {
        if self.done {
            return None;
        }
        if let Some(error) = self.pending_error.take() {
            self.done = true;
            return Some(Err(error));
        }

        if !self.filled {
            self.filled = true;
            for i in 0..self.cursors.len() {
                if let Err(error) = self.refill(i) {
                    self.done = true;
                    return Some(Err(error));
                }
            }
        }

        let winner = self.argmin()?;

        // One move per read: `take` hands the read out and empties the slot, which is then
        // refilled from its own cursor. No clone.
        let read = self.heads[winner].take().expect("the argmin found a head");
        self.keys[winner] = None;
        if let Err(error) = self.refill(winner) {
            // The read in hand is good and goes out first; the failure is yielded next.
            self.pending_error = Some(error);
        }
        Some(Ok(read))
    }
}

/// Where a read sits in the genome. Sound as a cross-file comparison key only because the
/// open gate proved every file's `ref_id`s are the reference's `ContigId`s.
fn key_of(read: &AlignedRead) -> GenomePosition {
    GenomePosition {
        // PANIC-FREE: `ref_id` comes from a 32-bit field in both formats, so it fits by
        // construction. Checked rather than `as`-cast because a wrapped value would collapse
        // two contigs onto one key and silently disarm the ordering this merge rests on.
        contig: ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32")),
        position: Position(read.pos),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::ReadFilterConfig;
    use crate::ng::read::input::record_reader::{InMemoryRecordReader, RecordReader};
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, bam_header, fixture_read_group, matching_contigs, read_named_with_length,
    };
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::GenomeRegion;
    use noodles_sam::alignment::RecordBuf;
    use std::path::Path;
    use std::sync::Arc;

    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn read_at(qname: &str, start: usize) -> RecordBuf {
        read_named_with_length(qname, 0, start, 30)
    }

    fn cursor_over(name: &str, records: Vec<RecordBuf>) -> AlignmentCursor<InMemoryRefSeq> {
        AlignmentCursor::over_records(
            RecordReader::InMemory(InMemoryRecordReader::new(
                bam_header(&matching_contigs()),
                records,
            )),
            ContigId(0),
            fixture_read_group(),
            reference_bases(),
            ReadFilterConfig::default(),
            Arc::from(Path::new(name)),
        )
        .expect("the fixture header resolves")
    }

    fn names_of(sample: &mut SampleCursor<InMemoryRefSeq>, asked: GenomeRegion) -> Vec<String> {
        sample.move_to_region(asked).expect("on this chromosome");
        let mut names = Vec::new();
        while let Some(read) = sample.next_read() {
            let read = read.expect("the scripted reads decode");
            names.push(String::from_utf8_lossy(&read.qname).into_owned());
        }
        names
    }

    /// **A one-file sample is the `Single` arm**, which is the absence of a merge rather than
    /// a merge with k = 1 — the property `SampleRegionReads` exists as an enum to keep, on
    /// the hottest loop in the module.
    #[test]
    fn a_one_file_sample_takes_the_arm_with_no_merge_in_it() {
        let sample = SampleCursor::new(vec![cursor_over("/a.bam", vec![read_at("r0", 1)])]);

        assert!(matches!(sample, SampleCursor::Single(_)));
    }

    /// Two files' reads come out **in position order across both**, which is what the layer
    /// above assumes and what nothing else would notice going wrong: an out-of-order read
    /// reaches the walker, which is where it becomes a hard error rather than a wrong answer.
    #[test]
    fn two_files_are_merged_in_position_order() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("a1", 1), read_at("a3", 31)]),
            cursor_over("/b.bam", vec![read_at("b2", 16), read_at("b4", 46)]),
        ]);

        assert!(matches!(sample, SampleCursor::Merged(_)));
        assert_eq!(
            names_of(&mut sample, region(1, 100)),
            ["a1", "b2", "a3", "b4"],
        );
    }

    /// **Ties break to the lowest file index**, and that is what makes a run reproducible: a
    /// sample's files usually cover the same coordinate range, so equal positions are routine
    /// rather than an edge case.
    #[test]
    fn reads_at_the_same_position_break_to_the_first_file() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("from-a", 16)]),
            cursor_over("/b.bam", vec![read_at("from-b", 16)]),
        ]);

        assert_eq!(names_of(&mut sample, region(1, 100)), ["from-a", "from-b"]);
    }

    /// A file with nothing in the region does not stall the others — the merge has to skip an
    /// exhausted cursor rather than wait on it.
    #[test]
    fn a_file_with_no_reads_in_the_region_does_not_stall_the_others() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("a1", 1), read_at("a2", 16)]),
            cursor_over("/b.bam", vec![read_at("b-far", 61)]),
        ]);

        assert_eq!(names_of(&mut sample, region(1, 40)), ["a1", "a2"]);
    }

    /// **Every file is repositioned, and the heads are dropped with the region that primed
    /// them.** A head left over from the last region would be handed to this one as though it
    /// belonged — a read in the wrong region, from a file the caller cannot see.
    #[test]
    fn moving_to_a_region_repositions_every_file_and_drops_the_old_heads() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("a1", 1), read_at("a3", 61)]),
            cursor_over("/b.bam", vec![read_at("b2", 16), read_at("b4", 66)]),
        ]);

        // Prime both heads, then abandon the region after one read.
        sample
            .move_to_region(region(1, 40))
            .expect("on this chromosome");
        assert!(sample.next_read().is_some());

        assert_eq!(names_of(&mut sample, region(55, 100)), ["a3", "b4"]);
    }

    /// The two arms must answer the same, or a sample's reads would depend on how many files
    /// it happened to be sequenced across.
    #[test]
    fn one_file_through_the_merge_answers_as_the_single_arm_does() {
        let script = vec![read_at("r0", 1), read_at("r1", 16), read_at("r2", 31)];

        let mut single = SampleCursor::new(vec![cursor_over("/a.bam", script.clone())]);
        let mut merged = SampleCursor::Merged(MergedCursors::new(vec![cursor_over(
            "/a.bam",
            script.clone(),
        )]));

        for asked in [region(1, 100), region(16, 40), region(1, 5)] {
            assert_eq!(names_of(&mut single, asked), names_of(&mut merged, asked));
        }
    }
}
