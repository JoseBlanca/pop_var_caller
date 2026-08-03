//! **One sample's reads over one chromosome** — the layer a caller actually holds.
//!
//! A sample is *k* files: one for most, several when it was sequenced in more than one run.
//! Both generators ask a *sample* for reads, not a file, so the cursor a caller keeps has to
//! be a sample-level one that holds a file cursor each and merges what they yield
//! (`arch/alignment_cursor.md` §2.4).

use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::read::filtering::ReadGroupCounts;
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::types::{ContigId, GenomePosition, Position};

use super::cursor::{AlignmentCursor, CursorError};

/// One sample's reads over one chromosome.
///
/// **An enum, and the single-file arm is the *absence* of a merge rather than a merge with
/// k = 1.** A one-file sample — which is most of them — pays no merge machinery at all. The
/// two arms are unified by an enum and not by `Box<dyn …>` because dynamic dispatch is opaque
/// to the optimiser: the whole per-read chain would stop being inlined at the boundary, on
/// the hottest loop in the module. The per-region sample stream this replaced was an enum for
/// exactly that reason, and this keeps the property.
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

    /// What this sample's cursors did, summed across its files — see
    /// [`CursorCounts`](crate::ng::read::input::cursor::CursorCounts). Summed rather than
    /// per file because the question it answers is about the sample's walk.
    pub fn counts(&self) -> crate::ng::read::input::cursor::CursorCounts {
        match self {
            Self::Single(cursor) => cursor.counts(),
            Self::Merged(merged) => merged.cursors.iter().fold(
                crate::ng::read::input::cursor::CursorCounts::default(),
                |mut total, cursor| {
                    total += cursor.counts();
                    total
                },
            ),
        }
    }

    /// Step-1's tally for this sample, **one entry per read group met, never summed across
    /// them**.
    ///
    /// A drop rate is a read group's property: a run that went badly shows up as one library
    /// with an anomalous rate, and adding them together erases exactly that signal. So the k
    /// files' tallies are folded together *by read group* rather than accumulated into one —
    /// which is also why a sample spanning two files with one library reports one entry, not
    /// two.
    ///
    /// A read group appears only once some read has met it, so a cursor that has walked
    /// nothing reports nothing rather than a row of zeroes.
    ///
    /// **This is what `SampleReads::counts` used to answer**, and the difference is where the
    /// number lives: the per-region streams folded their tallies back into the file as each
    /// one ended, so a stream still running was not yet reflected. The filter here lives as
    /// long as the cursor, so this is current the moment it is asked.
    pub fn read_group_counts(&self) -> Vec<ReadGroupCounts> {
        match self {
            Self::Single(cursor) => cursor.read_group_counts(),
            Self::Merged(merged) => {
                let mut total: Vec<ReadGroupCounts> = Vec::new();
                for cursor in &merged.cursors {
                    for (read_group, counts) in cursor.read_group_counts() {
                        match total.iter_mut().find(|(id, _)| *id == read_group) {
                            Some((_, running)) => running.add(&counts),
                            // First **seen**, appended — so the order is the order reads
                            // met these groups, file by file, not the order the headers
                            // declared them. A file whose first read is its second `@RG`
                            // reports that one first. Same rule the per-file fold this
                            // replaces had.
                            None => total.push((read_group, counts)),
                        }
                    }
                }
                total
            }
        }
    }

    /// Release the reference bases every file's read filter has gone past — see
    /// [`AlignmentCursor::evict_reference_before`].
    ///
    /// Each file's cursor holds its own reference reader, because they are stateful readers
    /// and sharing one between k files would give them one file position between them. So the
    /// release goes to all of them.
    pub fn evict_reference_before(&self, pos: u64)
    where
        R: crate::ng::ref_seq::EvictableRefSeq,
    {
        match self {
            Self::Single(cursor) => cursor.evict_reference_before(pos),
            Self::Merged(merged) => merged
                .cursors
                .iter()
                .for_each(|cursor| cursor.evict_reference_before(pos)),
        }
    }

    /// How many reference bases this sample's readers hold between them.
    pub fn resident_reference_bases(&self) -> usize
    where
        R: crate::ng::ref_seq::EvictableRefSeq,
    {
        match self {
            Self::Single(cursor) => cursor.resident_reference_bases(),
            Self::Merged(merged) => merged
                .cursors
                .iter()
                .map(|cursor| cursor.resident_reference_bases())
                .sum(),
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
/// the per-region merge this replaced used, and for the same reason.
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

    /// **The same-file-twice check, run only on a tie.**
    ///
    /// Two reads at different positions cannot be the same read, so this has nothing to do
    /// unless the argmin already found two heads sharing a position — which it learned as a
    /// by-product of the comparison it had to make anyway. Within a tie, `flag` (a `u16`) is
    /// compared before `qname` (a short string), so the string compare happens only for reads
    /// that already agree on position *and* flags.
    ///
    /// **Why it must exist here and not only at open.** `SampleReads::open` rejects a
    /// duplicate *path*, which a copy or a symlink walks straight past. Two copies of one read
    /// are two votes at a locus: the depth doubles and so does every allele count derived from
    /// it, which looks like real evidence rather than like a fault. The per-region merge this
    /// replaced always had this guard; the first version of this type dropped it, and two
    /// cursors over one file yielded every read twice with nothing failing.
    fn same_read_across_files(&self, winner: usize) -> Option<CursorError> {
        let winning_key = self.keys[winner]?;
        let winning_read = self.heads[winner].as_ref()?;

        for (other, key) in self.keys.iter().enumerate().skip(winner + 1) {
            if *key != Some(winning_key) {
                continue;
            }
            let Some(other_read) = self.heads[other].as_ref() else {
                continue;
            };
            if other_read.flag == winning_read.flag && other_read.qname == winning_read.qname {
                return Some(CursorError::DuplicateReadAcrossFiles {
                    qname: winning_read.qname.clone(),
                    contig: winning_key.contig,
                    position: winning_key.position.get(),
                    first_file: winner,
                    second_file: other,
                });
            }
        }
        None
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

        if let Some(error) = self.same_read_across_files(winner) {
            self.done = true;
            return Some(Err(error));
        }

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
    /// a merge with k = 1 — the property the per-region enum this replaced existed to keep, on
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

    /// **Three files, to prove the argmin is not accidentally a two-way compare.**
    ///
    /// Every other test of the merge uses two files, and with two files a scan that examines
    /// only the first two slots is indistinguishable from one that examines all of them. So a
    /// sample of three or more would silently lose every read of its third file and later —
    /// which reads as a sample sequenced less deeply, not as a fault.
    ///
    /// Inherited from the per-region merge, whose own version of this test carried the same
    /// sentence. Milestone F replaced that merge and this was the one rule of it left without
    /// a successor: mutating [`MergedCursors::argmin`]'s scan to `.take(2)` left the whole
    /// suite green.
    ///
    /// The files are given **out of coordinate order** — the third holds the earliest read —
    /// so a merge that simply concatenated in file order would fail too.
    #[test]
    fn three_files_merge_in_coordinate_order() {
        let mut sample = SampleCursor::new(vec![
            // 30-base reads on a 100-base fixture contig, so the last starts at 51.
            cursor_over("/a.bam", vec![read_at("a2", 11), read_at("a5", 41)]),
            cursor_over("/b.bam", vec![read_at("b3", 21), read_at("b6", 51)]),
            cursor_over("/c.bam", vec![read_at("c1", 1), read_at("c4", 31)]),
        ]);

        assert_eq!(
            names_of(&mut sample, region(1, 100)),
            ["c1", "a2", "b3", "c4", "a5", "b6"],
            "every file's reads must reach the merge, interleaved by position",
        );
    }

    /// **"All or nothing" has to be a fact, not a claim.** A move that fails part-way leaves
    /// the earlier files pointed at the new region and the later ones where they were, with
    /// heads still primed — and reads from a subset of a sample's files are real reads at real
    /// positions, so nothing downstream can tell they are half an answer. Measured before the
    /// latch: 13 reads served after `move_to_region` returned `Err`.
    #[test]
    fn a_failed_move_serves_nothing_afterwards() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("a0", 1), read_at("a1", 31)]),
            cursor_over("/b.bam", vec![read_at("b0", 1), read_at("b1", 31)]),
        ]);
        // Prime both files on a good region first.
        assert!(!names_of(&mut sample, region(1, 100)).is_empty());

        // A region on another chromosome: the first file refuses, and the second is never
        // reached.
        let foreign = GenomeRegion {
            contig: ContigId(1),
            start: Position(1),
            end: Position(50),
        };
        assert!(sample.move_to_region(foreign).is_err());

        assert!(
            sample.next_read().is_none(),
            "a sample whose move failed answered anyway, from the files that had moved",
        );
    }

    /// The `Merged` arm's accessors must answer for the sample, not for whichever file
    /// happens to be first — `counts` sums, `contig` is the one they all cover.
    #[test]
    fn the_merged_arm_reports_the_whole_samples_work() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("a0", 1)]),
            cursor_over("/b.bam", vec![read_at("b0", 31)]),
        ]);

        assert_eq!(sample.contig(), ContigId(0));
        assert_eq!(sample.counts(), Default::default());
        let _ = names_of(&mut sample, region(1, 100));

        let counts = sample.counts();
        assert_eq!(
            counts.reads_decoded, 2,
            "the sample decoded one read from each of its two files",
        );
        assert_eq!(
            counts.regions_jumping, 2,
            "one first region per file, and both had to jump",
        );
    }

    /// **The guard the first version of this type dropped.** Two cursors over the same file
    /// are two votes at every locus: the depth doubles and so does every allele count derived
    /// from it, which looks like evidence rather than like a fault.
    ///
    /// `SampleReads::open` rejects a duplicate *path*; a copy or a symlink walks straight
    /// past it, which is why the check has to be here too. Measured before the fix: the same
    /// script through two cursors yielded `["r0", "r0", "r1", "r1"]` with no error.
    #[test]
    fn the_same_read_from_two_files_is_a_hard_error() {
        let script = vec![read_at("r0", 1), read_at("r1", 31)];
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", script.clone()),
            cursor_over("/a-copy.bam", script),
        ]);

        sample
            .move_to_region(region(1, 100))
            .expect("on this chromosome");
        let first = sample.next_read().expect("a read");
        assert!(
            matches!(first, Err(CursorError::DuplicateReadAcrossFiles { .. })),
            "two files holding the same read must be refused, not counted twice: {first:?}",
        );
    }

    /// And the guard must not fire on reads that merely *share a position* — which is
    /// routine, because a sample's files usually cover the same coordinate range.
    #[test]
    fn different_reads_at_the_same_position_are_not_a_duplicate() {
        let mut sample = SampleCursor::new(vec![
            cursor_over("/a.bam", vec![read_at("from-a", 16)]),
            cursor_over("/b.bam", vec![read_at("from-b", 16)]),
        ]);

        assert_eq!(names_of(&mut sample, region(1, 100)), ["from-a", "from-b"]);
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
