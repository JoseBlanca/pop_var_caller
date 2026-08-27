//! The chain ids: which reads are live, and what changes between one record and the next.
//!
//! **A chain id names the read that produced a piece of evidence** — one identifier per read
//! pair, with the mates collapsed onto it, allocated in order and never reused. Since the
//! owner's ruling of 2026-08-17 ng names *every* read it folds, including the ones that agreed
//! with the reference, because the cohort merge has to tell a read that covered a position and
//! agreed from a read that never reached it. Unnamed, those are the same absence.
//!
//! **That makes this the field that decides the file's size at depth.** A read of length `L` is
//! named at every one of the `L` positions it covers, so the naive column stores each read about
//! `L` times: measured, the ids are 16 % of the file at eleven reads a position and **89 % of it
//! at three hundred** (spec [`psp_record_encoding.md`] §6).
//!
//! # What this module stores instead
//!
//! **Only the changes.** At each record, which ids started covering it and which stopped; a
//! reader carries the set forward. A read covering 150 positions then costs about two entries in
//! total rather than appearing in 150 lists. Measured on the same two corners: 0.432 bytes a
//! position at 11.4 reads and 6.42 at 293, against 43.78 for the whole list as raw identifiers.
//!
//! **The set restarts at every block** (spec [`psp_file_format.md`] §3.2), which is what lets a
//! reader begin at any block: [`start_block`](LiveSetWriter::start_block) empties it, so a
//! block's first record has no departures and its arrivals *are* the whole live set. The
//! restatement is not a separate field — it falls out of the reset, and costs 12 % of this
//! form's own bytes at the settled block size.
//!
//! # What is not here
//!
//! - **The exception lists** — the ids of every observation except the residual one — stay in
//!   the record's skippable body, because they carry no state. This module holds the half that
//!   does. ([`record`](super::record), Milestone E3.)
//! - **The residual observation's ids**, which are the live set minus every other
//!   observation's. That is where the scheme is cheap and where it fails silently, and it is
//!   Milestone E4's, with the count check spec [`psp_chain_id_encoding.md`] §5 names.
//! - **Re-entry** is Milestone E2's oracle rather than a change to this encoding: an id may go
//!   live, stop, and go live again, because a pair's mates rarely overlap — 83 % of ids on the
//!   human sample and 91 % on tomato cover two stretches. Arrivals and departures already
//!   express that; what E2 owes is a fixture that *contains* it and a test that counts it, and
//!   nothing here may grow an assumption of one stretch per id in the meantime.
//!
//! [`psp_file_format.md`]: ../../../../doc/devel/ng/spec/psp_file_format.md
//! [`psp_record_encoding.md`]: ../../../../doc/devel/ng/spec/psp_record_encoding.md
//! [`psp_chain_id_encoding.md`]: ../../../../doc/devel/ng/spec/psp_chain_id_encoding.md

use crate::pileup_record::ChainId;

use super::record::{FieldReader, RecordDecodeError, put_varint};

/// The name this module's faults are reported under, so a message says which field it was.
const DEPARTURE_COUNT: &str = "chain-id departure count";
const DEPARTURE_POSITION: &str = "chain-id departure position";
const ARRIVAL_COUNT: &str = "chain-id arrival count";
const ARRIVAL_ID: &str = "chain-id arrival";

/// The fewest bytes one departure or one arrival can take, for the count bound.
///
/// **A count is refused when no record could hold that many entries**, the same guard
/// [`FieldReader::read_count`] applies to observations: a declared 2⁴⁰ arrivals is not a buffer
/// that stopped early, and reporting it as one would ask Milestone D's reader to grow its buffer
/// to a terabyte instead of reporting damage.
const LEAST_BYTES_PER_ENTRY: usize = 1;

/// The reads live at one record — **sorted ascending and without duplicates**, which is what
/// makes a departure spellable as a position rather than as an identifier.
///
/// **Its own type rather than a bare `Vec<ChainId>`**, because this module deals in two kinds of
/// integer that must never be transposed — identifiers, which run into the millions, and
/// positions in this set, which are small (arch §7 names that hazard). A `Vec<u64>` for both
/// would make the transposition a silent wrong answer; here only the codec sees a position at
/// all, and it converts back to an identifier before anything else does.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LiveSet {
    ids: Vec<ChainId>,
}

impl LiveSet {
    /// The set a block begins with: empty.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The live ids, ascending.
    #[must_use]
    pub fn ids(&self) -> &[ChainId] {
        &self.ids
    }

    /// How many reads are live.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Whether no read is live — which is a block's opening state, and an ordinary one at a gap
    /// in coverage.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Whether `id` is live. A binary search, because the set is sorted.
    #[must_use]
    pub fn holds(&self, id: ChainId) -> bool {
        self.ids.binary_search(&id).is_ok()
    }
}

/// What changed between the previous record and this one.
///
/// **Both halves are identifiers**, not positions. The positions exist only inside the encoding,
/// where a departure is spelled as an index into the live set because that is a small number
/// where an identifier is a large one.
///
/// **Constructed only by deriving it from two sets or by decoding it**, so applying one needs no
/// check of its own — the checks happen where the bytes are, which is where there is an offset to
/// report.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LiveSetChanges {
    departed: Vec<ChainId>,
    arrived: Vec<ChainId>,
}

impl LiveSetChanges {
    /// The ids that stopped covering, ascending.
    #[must_use]
    pub fn departed(&self) -> &[ChainId] {
        &self.departed
    }

    /// The ids that started covering, ascending.
    #[must_use]
    pub fn arrived(&self) -> &[ChainId] {
        &self.arrived
    }

    /// Whether nothing changed — the common case at depth between two adjacent positions, and
    /// the reason this form is cheap.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.departed.is_empty() && self.arrived.is_empty()
    }

    fn clear(&mut self) {
        self.departed.clear();
        self.arrived.clear();
    }
}

/// Sort `ids` ascending and drop duplicates, in place.
///
/// **A record's ids arrive as the union of its observations' lists**, which is neither sorted
/// nor distinct: one read pair that showed the same sequence at two observations of one record
/// is named twice, and observations are stored in no particular id order.
fn sort_and_dedup(ids: &mut Vec<ChainId>) {
    ids.sort_unstable();
    ids.dedup();
}

/// Fill `changes` with what it takes to get from `live` to `wanted`.
///
/// Both inputs are sorted and distinct, so this is one merge pass rather than two set
/// differences: an id in `live` and not in `wanted` departed, one in `wanted` and not in `live`
/// arrived, and one in both did not change — which is nearly all of them.
fn derive_changes(live: &LiveSet, wanted: &[ChainId], changes: &mut LiveSetChanges) {
    changes.clear();
    let (mut in_live, mut in_wanted) = (0usize, 0usize);
    while in_live < live.ids.len() && in_wanted < wanted.len() {
        match live.ids[in_live].cmp(&wanted[in_wanted]) {
            std::cmp::Ordering::Equal => {
                in_live += 1;
                in_wanted += 1;
            }
            std::cmp::Ordering::Less => {
                changes.departed.push(live.ids[in_live]);
                in_live += 1;
            }
            std::cmp::Ordering::Greater => {
                changes.arrived.push(wanted[in_wanted]);
                in_wanted += 1;
            }
        }
    }
    changes.departed.extend_from_slice(&live.ids[in_live..]);
    changes.arrived.extend_from_slice(&wanted[in_wanted..]);
}

/// Take `departed` out of `live`. Both are ascending, so this is one pass.
///
/// **Infallible, and that is a property of how the two lists can be made**: deriving them gives
/// departures that are all live and arrivals that are none of them, and decoding them refuses
/// anything else before it returns. Nothing else produces either.
fn apply_departures(live: &mut LiveSet, departed: &[ChainId]) {
    if departed.is_empty() {
        return;
    }
    let mut departing = departed.iter().peekable();
    live.ids.retain(|id| match departing.peek() {
        Some(next) if *next == id => {
            departing.next();
            false
        }
        _ => true,
    });
}

/// Put `arrived` into `live`, keeping it ascending.
///
/// **A merge rather than an append and a sort**: both sides are already ascending, and at depth
/// the live set is hundreds of ids where the arrivals are a handful, so re-sorting the whole set
/// once a record would be the expensive part of the walk.
fn apply_arrivals(live: &mut LiveSet, arrived: &[ChainId], scratch: &mut Vec<ChainId>) {
    if arrived.is_empty() {
        return;
    }
    scratch.clear();
    scratch.reserve(live.ids.len() + arrived.len());
    let (mut in_live, mut in_arrived) = (0usize, 0usize);
    while in_live < live.ids.len() && in_arrived < arrived.len() {
        if live.ids[in_live] < arrived[in_arrived] {
            scratch.push(live.ids[in_live]);
            in_live += 1;
        } else {
            scratch.push(arrived[in_arrived]);
            in_arrived += 1;
        }
    }
    scratch.extend_from_slice(&live.ids[in_live..]);
    scratch.extend_from_slice(&arrived[in_arrived..]);
    std::mem::swap(&mut live.ids, scratch);
}

/// Write a record's live-set changes, and carry the set forward.
///
/// **Nothing calls this yet.** Milestone E3 puts the bytes it writes into the record head, which
/// is where they have to be: they carry state across records, so a reader that skipped a body
/// holding them would have a stale set from that point on. E1 is the stream and its oracle; E3
/// is where it sits in a record.
///
/// One of these per sample being written. **The live set is the per-block state**, and
/// [`start_block`](Self::start_block) is its whole reset — the same shape `RecordEncoder` and
/// `BlockCursor` use for theirs, for the reason spec §3.2 gives: a difference that survives a
/// block boundary reads back wrong from that block's first record, and plausibly wrong.
#[derive(Debug, Default)]
pub struct LiveSetWriter {
    // ---- what lives for the whole file: scratch, reused so a record does not allocate ----
    /// The ids this record wants, gathered, sorted and deduplicated.
    wanted: Vec<ChainId>,
    /// This record's changes.
    changes: LiveSetChanges,
    /// Where the live set is rebuilt when arrivals are merged into it. Swapped with the set's own
    /// vector rather than copied back, so the two trade places and neither reallocates.
    merged: Vec<ChainId>,

    // ---- what lives for one block, and is replaced whole at every boundary ----
    //
    // **Nothing per-block goes beside `block`.** A field here is initialised once and never
    // reset; a field inside `PerBlockState` cannot be, because `at_block_start` is its only
    // constructor and its literal is exhaustive.
    block: PerBlockState,
}

/// Everything the chain-id stream **restarts when a block does** (spec §3.2).
///
/// **A struct for one field, and deliberately.** It is the same shape `RecordEncoder` and
/// `BlockCursor` carry theirs in, for the same reason: Milestone E2's re-entry bookkeeping and
/// E4's residual counts are the next things that might live for a block, and a field added beside
/// this one and initialised once is one that silently never resets — a file that then reads back
/// wrong from each block's first record, and plausibly wrong.
///
/// **Not `Copy`**, so a collection added here gives the `error[E0063]` that points at the reset
/// rather than an `error[E0204]` whose cheapest answer is to put the field somewhere else.
#[derive(Debug, Default, Clone)]
struct PerBlockState {
    /// The reads live at the record last handled.
    live: LiveSet,
}

impl PerBlockState {
    /// A block's opening state: nothing is live, so the first record restates everything as
    /// arrivals and there is no separate field to carry the restatement.
    fn at_block_start() -> Self {
        Self {
            live: LiveSet::new(),
        }
    }
}

impl LiveSetWriter {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Start a block: nothing is live, so the next record restates everything as arrivals.
    ///
    /// **⚠ A writer that forgets to call this at a block boundary produces a file that is wrong
    /// from that block's first record and parses perfectly** — the reader restarts whether the
    /// writer did or not, so it reconstructs a set the writer never meant. The call site is
    /// Milestone E3's, and so is the test that holds it.
    pub fn start_block(&mut self) {
        self.block = PerBlockState::at_block_start();
    }

    /// The reads live at the record last written. **After** `write_changes`, not before.
    #[must_use]
    pub fn live(&self) -> &LiveSet {
        &self.block.live
    }

    /// Append this record's changes to `out`, and move the live set to what `live_now` names.
    ///
    /// `live_now` is every id the record names, in any order and with any repeats — the union of
    /// its observations' lists.
    pub fn write_changes(
        &mut self,
        live_now: impl IntoIterator<Item = ChainId>,
        out: &mut Vec<u8>,
    ) {
        self.wanted.clear();
        self.wanted.extend(live_now);
        sort_and_dedup(&mut self.wanted);
        derive_changes(&self.block.live, &self.wanted, &mut self.changes);
        encode_changes(&self.changes, &self.block.live, out);
        apply_departures(&mut self.block.live, &self.changes.departed);
        apply_arrivals(
            &mut self.block.live,
            &self.changes.arrived,
            &mut self.merged,
        );
    }
}

/// **Departures first, then arrivals**, and the order is load-bearing rather than a convention.
///
/// A departure is written as its *position* in the live set, which is small where an identifier
/// is large. Positions only mean anything against a particular set, so writing the departures
/// first — and applying them before the arrivals are read — means the set a position indexes is
/// always the set as it stands at that moment. Written the other way round, both sides would
/// have to agree to resolve positions against a set that no longer existed, which is the kind of
/// agreement that holds until someone reorders two lines.
fn encode_changes(changes: &LiveSetChanges, live: &LiveSet, out: &mut Vec<u8>) {
    put_varint(out, changes.departed.len() as u64);
    let mut previous_position: Option<u64> = None;
    for id in &changes.departed {
        let position = live
            .ids
            .binary_search(id)
            .expect("a derived departure is live") as u64;
        put_varint(out, gap_from(previous_position, position));
        previous_position = Some(position);
    }

    put_varint(out, changes.arrived.len() as u64);
    let mut previous_id: Option<u64> = None;
    for id in &changes.arrived {
        put_varint(out, gap_from(previous_id, *id));
        previous_id = Some(*id);
    }
}

/// What to write for a value in a strictly ascending run: the first absolutely, each later one
/// as the gap since the one before.
///
/// **Biased by one**, because the run is strictly ascending: consecutive values differ by at
/// least one, so storing `next - previous - 1` puts a dense run — which is what monotonically
/// allocated identifiers give — in a single byte each.
fn gap_from(previous: Option<u64>, value: u64) -> u64 {
    match previous {
        None => value,
        Some(previous) => value - previous - 1,
    }
}

/// Read a record's live-set changes, and carry the set forward.
///
/// The reading half of [`LiveSetWriter`], with the same per-block reset.
#[derive(Debug, Default)]
pub struct LiveSetReader {
    // ---- what lives for the whole file ----
    changes: LiveSetChanges,
    /// Scratch for the arrival merge — see [`LiveSetWriter`]'s field of the same name.
    merged: Vec<ChainId>,

    // ---- what lives for one block, and is replaced whole at every boundary ----
    block: PerBlockState,
}

impl LiveSetReader {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Start a block: nothing is live, and the block's first record restates everything.
    pub fn start_block(&mut self) {
        self.block = PerBlockState::at_block_start();
    }

    /// The reads live at the record last read.
    #[must_use]
    pub fn live(&self) -> &LiveSet {
        &self.block.live
    }

    /// What changed at the record last read. Handed back so a caller can see the shape of the
    /// stream without diffing two sets.
    #[must_use]
    pub fn changes(&self) -> &LiveSetChanges {
        &self.changes
    }

    /// Read one record's changes off the front of `bytes` and apply them, returning how many
    /// bytes they took.
    ///
    /// **A stream that is internally inconsistent is corrupt input, not a panic** — a departure
    /// for a position past the end of the set, or an arrival for an id already live. Both are
    /// `Malformed` rather than `Truncated`: the bytes were there and cannot mean what they say,
    /// so no quantity of further bytes changes the answer and Milestone D's reader must not
    /// fetch more and retry.
    ///
    /// **⚠ A fault's `bytes_in` counts from `bytes`'s own first byte**, not from the record's.
    /// The caller that puts this stream inside a record re-bases it, the way `record.rs` already
    /// re-bases a body's faults into the record that holds it.
    ///
    /// **A slice rather than this module's byte cursor**, so that where the stream sits inside a
    /// record stays E3's question and not a constraint the type imposes now.
    pub fn read_changes(&mut self, bytes: &[u8]) -> Result<usize, RecordDecodeError> {
        let reader = &mut FieldReader::new(bytes);
        self.changes.clear();

        let departures = reader.read_count(DEPARTURE_COUNT, LEAST_BYTES_PER_ENTRY)?;
        let mut previous_position: Option<u64> = None;
        for _ in 0..departures {
            let gap = reader.read_varint(DEPARTURE_POSITION)?;
            let position = value_after(previous_position, gap).ok_or_else(|| {
                reader.malformed(
                    DEPARTURE_POSITION,
                    format!(
                        "a gap of {gap} past position {}, which is past every position there is",
                        previous_position.unwrap_or(0)
                    ),
                )
            })?;
            let id = *usize::try_from(position)
                .ok()
                .and_then(|position| self.block.live.ids.get(position))
                .ok_or_else(|| {
                    reader.malformed(
                        DEPARTURE_POSITION,
                        format!(
                            "position {position} of a live set holding {} reads",
                            self.block.live.len()
                        ),
                    )
                })?;
            self.changes.departed.push(id);
            previous_position = Some(position);
        }
        // **Applied before the arrivals are read**, so the positions above indexed the set they
        // were written against — see `encode_changes`.
        apply_departures(&mut self.block.live, &self.changes.departed);

        let arrivals = reader.read_count(ARRIVAL_COUNT, LEAST_BYTES_PER_ENTRY)?;
        let mut previous_id: Option<u64> = None;
        for _ in 0..arrivals {
            let gap = reader.read_varint(ARRIVAL_ID)?;
            let id = value_after(previous_id, gap).ok_or_else(|| {
                reader.malformed(
                    ARRIVAL_ID,
                    format!(
                        "a gap of {gap} past id {}, which is past every id there is",
                        previous_id.unwrap_or(0)
                    ),
                )
            })?;
            if self.block.live.holds(id) {
                return Err(reader.malformed(
                    ARRIVAL_ID,
                    format!(
                        "id {id}, which is already live — an arrival names a read that was not"
                    ),
                ));
            }
            self.changes.arrived.push(id);
            previous_id = Some(id);
        }
        apply_arrivals(
            &mut self.block.live,
            &self.changes.arrived,
            &mut self.merged,
        );

        Ok(reader.bytes_read())
    }
}

/// The value a gap names, or `None` when it names one past `u64`.
///
/// **The inverse of [`gap_from`]**, and it is fallible where that one is not: a gap comes off
/// the wire and may say anything, so a run walking past `u64::MAX` is damage rather than a wrap.
fn value_after(previous: Option<u64>, gap: u64) -> Option<u64> {
    match previous {
        None => Some(gap),
        Some(previous) => previous.checked_add(gap)?.checked_add(1),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write a whole file's worth of records and read them back, restarting at every block.
    ///
    /// `blocks` is one entry per block, each a list of records, each record the ids it names —
    /// in whatever order and with whatever repeats the caller gives, because that is the shape a
    /// record's observations hand over.
    fn round_trip(blocks: &[Vec<Vec<ChainId>>]) -> (Vec<Vec<ChainId>>, usize) {
        let written = write_a_file(blocks);

        let mut reader = LiveSetReader::new();
        let mut back = Vec::new();
        let mut at = 0usize;
        for block in blocks {
            reader.start_block();
            for _ in block {
                let used = reader
                    .read_changes(&written[at])
                    .unwrap_or_else(|refused| panic!("record {at}: {refused}"));
                assert_eq!(
                    used,
                    written[at].len(),
                    "record {at} left bytes of its own stream unread"
                );
                back.push(reader.live().ids().to_vec());
                at += 1;
            }
        }
        (back, written.iter().map(Vec::len).sum())
    }

    /// One block at a time, one record at a time: the bytes a writer lays down for each record.
    fn write_a_file(blocks: &[Vec<Vec<ChainId>>]) -> Vec<Vec<u8>> {
        let mut writer = LiveSetWriter::new();
        let mut written = Vec::new();
        for block in blocks {
            writer.start_block();
            for record in block {
                let mut bytes = Vec::new();
                writer.write_changes(record.iter().copied(), &mut bytes);
                written.push(bytes);
            }
        }
        written
    }

    /// What the ids of one record ought to read back as: sorted, without duplicates.
    fn as_a_set(ids: &[ChainId]) -> Vec<ChainId> {
        let mut set = ids.to_vec();
        sort_and_dedup(&mut set);
        set
    }

    /// Coverage that looks like reads: `reads_starting` reads of `read_length` bases begin at
    /// every position, ids allocated in order, and a read is live over the positions it covers.
    ///
    /// One record a position, which is the generic path's ordinary case.
    fn a_sliding_coverage(
        positions: u64,
        reads_starting: u64,
        read_length: u64,
    ) -> Vec<Vec<ChainId>> {
        (0..positions)
            .map(|at| {
                let earliest = at.saturating_sub(read_length - 1);
                (earliest..=at)
                    .flat_map(|start| {
                        (0..reads_starting).map(move |which| start * reads_starting + which)
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn a_run_of_records_reads_back_the_set_each_one_named() {
        let records = a_sliding_coverage(200, 3, 40);
        let wanted: Vec<_> = records.iter().map(|ids| as_a_set(ids)).collect();
        let (back, _) = round_trip(&[records]);
        assert_eq!(back, wanted);
    }

    /// **The ids a record names arrive unsorted and with repeats**, because they are the union of
    /// its observations' lists and one read pair can be named by two of them. The set is the
    /// writer's to build; a caller sorting first would be a precondition nothing enforces.
    #[test]
    fn unsorted_and_repeated_ids_name_the_same_set() {
        let records = vec![vec![900, 4, 4, 71, 4, 900], vec![71, 900, 71]];
        let (back, _) = round_trip(&[records]);
        assert_eq!(back, vec![vec![4, 71, 900], vec![71, 900]]);
    }

    /// **A block restates the whole live set, and that is what lets a reader start at one.**
    ///
    /// The same records read from the file's beginning and read from the second block's first
    /// record give the same sets from that point on — which they cannot if any part of the set
    /// was carried across the boundary. This is spec §3.2's property for this field.
    #[test]
    fn a_block_read_alone_gives_the_sets_it_gives_in_the_middle_of_a_file() {
        let first = a_sliding_coverage(60, 2, 25);
        let second: Vec<Vec<ChainId>> = a_sliding_coverage(60, 2, 25)
            .into_iter()
            .map(|ids| ids.into_iter().map(|id| id + 10_000).collect())
            .collect();
        let blocks = [first.clone(), second.clone()];

        // **One file, read twice — and reading the file's own bytes is the point.** A test that
        // re-wrote the second block on its own would pass even if the writer *and* the reader
        // both carried the set across the boundary, because the two would agree with each other.
        // That is the failure spec §3.2 names: self-consistent, and only a reader starting
        // mid-file sees it.
        let written = write_a_file(&blocks);

        let (whole, _) = round_trip(&blocks);

        let mut reader = LiveSetReader::new();
        reader.start_block();
        let mut alone = Vec::new();
        for bytes in &written[first.len()..] {
            reader.read_changes(bytes).expect("the second block reads");
            alone.push(reader.live().ids().to_vec());
        }

        assert_eq!(
            &whole[first.len()..],
            alone.as_slice(),
            "a reader starting at the second block must get what a reader that read the first \
             one gets from there"
        );
        assert_eq!(
            alone,
            second.iter().map(|ids| as_a_set(ids)).collect::<Vec<_>>(),
            "and both must be the sets the writer was given"
        );
    }

    /// **The restatement is the first record's arrivals**, which is why no field carries it.
    #[test]
    fn a_blocks_first_record_arrives_with_everything_and_departs_nothing() {
        let mut writer = LiveSetWriter::new();
        let mut reader = LiveSetReader::new();
        let mut bytes = Vec::new();

        writer.start_block();
        writer.write_changes([5u64, 9, 12], &mut bytes);
        writer.write_changes([9u64, 12, 40], &mut Vec::new());

        // A second block whose first record is that same set: it must restate all three.
        let mut restated = Vec::new();
        writer.start_block();
        writer.write_changes([9u64, 12, 40], &mut restated);

        reader.start_block();
        reader.read_changes(&restated).expect("it reads");
        assert_eq!(reader.changes().arrived(), [9, 12, 40]);
        assert!(reader.changes().departed().is_empty());
        assert_eq!(reader.live().ids(), [9, 12, 40]);

        // And the same set written *inside* a block costs almost nothing, which is the contrast.
        reader.start_block();
        reader.read_changes(&bytes).expect("it reads");
        assert_eq!(reader.live().ids(), [5, 9, 12]);
    }

    /// **A departure is a position in the live set, not an identifier**, and that is worth a test
    /// because both are `u64` and the saving is the whole reason the departures half is cheap.
    #[test]
    fn a_departure_costs_a_position_and_not_an_identifier() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();
        // Identifiers past 2²¹, so each needs four bytes as a varint.
        let live = [4_000_000u64, 4_000_001, 4_000_002];
        writer.write_changes(live, &mut Vec::new());

        let mut bytes = Vec::new();
        writer.write_changes([4_000_001u64, 4_000_002], &mut bytes);

        // one departure, position 0, no arrivals — and nothing else.
        assert_eq!(
            bytes,
            [1, 0, 0],
            "a departure spelled as its identifier would cost four bytes more"
        );
    }

    /// **Nothing changing costs two bytes**, which is the common case at depth between two
    /// adjacent positions and the reason the whole form is cheap.
    #[test]
    fn a_record_that_changes_nothing_costs_two_bytes() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();
        writer.write_changes([7u64, 8, 9], &mut Vec::new());

        let mut bytes = Vec::new();
        writer.write_changes([7u64, 8, 9], &mut bytes);
        assert_eq!(bytes, [0, 0], "no departures and no arrivals");
    }

    /// **The changes cost a fraction of the lists**, which is what the whole encoding is for.
    ///
    /// Measured against the cheapest honest alternative rather than the naive one: each record's
    /// ids written as ascending varint gaps — spec `psp_record_encoding.md` §6's middle arm,
    /// which is the one that captures 60 % to 86 % of the available saving on its own. Beating
    /// the raw-identifier arm would be a weaker claim.
    ///
    /// Measured on this fixture — 400 positions, three reads starting at each, reads 100 bases
    /// long, so about 300 live at once: **3,257 bytes against 106,166**, 32.6 times smaller. The
    /// assertion below asks for 16, which leaves the fixture room to change shape without
    /// becoming a number nobody re-measures.
    #[test]
    fn the_changes_cost_a_fraction_of_writing_every_record_s_list() {
        let records = a_sliding_coverage(400, 3, 100);
        let (_, changes_bytes) = round_trip(std::slice::from_ref(&records));

        let mut list_bytes = 0usize;
        for record in &records {
            let ids = as_a_set(record);
            let mut out = Vec::new();
            put_varint(&mut out, ids.len() as u64);
            let mut previous: Option<u64> = None;
            for id in ids {
                put_varint(&mut out, gap_from(previous, id));
                previous = Some(id);
            }
            list_bytes += out.len();
        }

        assert!(
            changes_bytes * 16 < list_bytes,
            "the changes must be a small fraction of the lists to be worth the state they \
             carry: {changes_bytes} bytes against {list_bytes}"
        );
    }

    /// **A departure naming a position the live set does not have is damage**, not a short read:
    /// the bytes were there and cannot mean what they say, so Milestone D's reader must refuse
    /// rather than fetch more bytes and try again.
    #[test]
    fn a_departure_past_the_live_set_is_damage() {
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&[0, 2, 4, 0]).expect("two arrive");
        assert_eq!(reader.live().ids(), [4, 5]);

        // one departure at position 7, of a set holding two.
        let refused = reader
            .read_changes(&[1, 7, 0])
            .expect_err("a position past the end is damage");
        assert!(
            matches!(refused, RecordDecodeError::Malformed { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("position 7"),
            "the message must say which position: {refused}"
        );
    }

    /// **An arrival for a read already live is damage too.** It is the other way a differential
    /// stream can be internally inconsistent in a way a list cannot, and left unchecked it would
    /// put an id in the set twice — which the residual arithmetic of Milestone E4 then subtracts
    /// once, silently.
    #[test]
    fn an_arrival_for_a_read_already_live_is_damage() {
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&[0, 1, 12]).expect("one arrives");
        assert_eq!(reader.live().ids(), [12]);

        let refused = reader
            .read_changes(&[0, 1, 12])
            .expect_err("12 is already live");
        assert!(
            matches!(refused, RecordDecodeError::Malformed { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("already live"),
            "the message must say why: {refused}"
        );
    }

    /// **A stream cut short anywhere is `Truncated` and never `Malformed`** — at every cut, one
    /// at a time. This is the line Milestone D's restartable reader branches on: a fault put in
    /// the wrong class here makes a streaming reader either reject a good record or retry for
    /// ever on a bad one.
    #[test]
    fn a_stream_cut_short_anywhere_is_a_short_read_and_never_damage() {
        let opening = [3u64, 400, 90_000, 12_000_000];
        let then = [400u64, 12_000_000, 12_000_001];

        let mut writer = LiveSetWriter::new();
        writer.start_block();
        let mut opening_bytes = Vec::new();
        writer.write_changes(opening, &mut opening_bytes);
        let mut bytes = Vec::new();
        writer.write_changes(then, &mut bytes);
        assert!(
            bytes.len() > 4,
            "the stream must have room to be cut in several places: {} bytes",
            bytes.len()
        );

        for cut in 0..bytes.len() {
            let mut reader = LiveSetReader::new();
            reader.start_block();
            reader
                .read_changes(&opening_bytes)
                .expect("the opening record reads");
            let refused = reader
                .read_changes(&bytes[..cut])
                .expect_err("a strict prefix of a record's stream cannot be complete");
            assert!(
                matches!(refused, RecordDecodeError::Truncated { .. }),
                "cutting at {cut} of {} gave {refused}",
                bytes.len()
            );
        }
    }

    /// **A count no stream could hold is damage rather than a short read**, the same bound
    /// `record.rs` puts on an observation count: a declared 2⁴⁰ arrivals would otherwise ask a
    /// streaming reader to grow its buffer to a terabyte instead of reporting a corrupt block.
    ///
    /// **Both counts, because there are two.** Testing only the arrivals left the departures'
    /// bound free: replacing its `read_count` with a plain varint read passed every other test
    /// in this module.
    #[test]
    fn a_count_larger_than_any_record_is_damage() {
        let mut absurd_arrivals = Vec::new();
        put_varint(&mut absurd_arrivals, 0);
        put_varint(&mut absurd_arrivals, u64::MAX / 2);

        let mut absurd_departures = Vec::new();
        put_varint(&mut absurd_departures, u64::MAX / 2);

        for (what, bytes) in [
            ("arrivals", absurd_arrivals),
            ("departures", absurd_departures),
        ] {
            let mut reader = LiveSetReader::new();
            reader.start_block();
            let refused = reader
                .read_changes(&bytes)
                .expect_err("no record holds that many entries");
            assert!(
                matches!(refused, RecordDecodeError::Malformed { .. }),
                "a count of {what} no record could hold gave {refused}"
            );
        }
    }

    /// **A gap that walks past the largest identifier there is, is damage** — not a wrap to a
    /// small id, which would put a read in the set that the writer never named.
    #[test]
    fn a_gap_past_the_largest_identifier_is_damage() {
        let mut past_the_end = Vec::new();
        put_varint(&mut past_the_end, 0);
        put_varint(&mut past_the_end, 2);
        put_varint(&mut past_the_end, u64::MAX);
        put_varint(&mut past_the_end, 1);

        let mut reader = LiveSetReader::new();
        reader.start_block();
        let refused = reader
            .read_changes(&past_the_end)
            .expect_err("the second arrival is past u64::MAX");
        assert!(
            matches!(refused, RecordDecodeError::Malformed { .. }),
            "got {refused}"
        );
    }

    /// A gap of nothing and a set that empties: coverage stops, and the next record starts again.
    #[test]
    fn a_set_that_empties_and_fills_again_reads_back() {
        let records = vec![
            vec![1u64, 2, 3],
            vec![],
            vec![],
            vec![50u64, 51],
            vec![51u64],
            vec![],
        ];
        let (back, _) = round_trip(&[records]);
        assert_eq!(
            back,
            vec![
                vec![1, 2, 3],
                vec![],
                vec![],
                vec![50, 51],
                vec![51],
                vec![],
            ]
        );
    }
}
