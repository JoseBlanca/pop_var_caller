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
//! total rather than appearing in 150 lists. Measured on the same two corners: **0.432 bytes a
//! position at 11.4 reads against 1.020 for the whole list as raw identifiers, and 6.42 at 293
//! reads against 43.78** — so 2.4 times smaller at the shallow corner and 6.8 times at the deep
//! one. (An earlier version of this paragraph gave 43.78 as *the* baseline for both, which reads
//! as a hundredfold saving at 11.4 reads.)
//!
//! **The set restarts at every block** (spec [`psp_file_format.md`] §3.2), which is what lets a
//! reader begin at any block: [`start_block`](LiveSetWriter::start_block) empties it, so a
//! block's first record has no departures and its arrivals *are* the whole live set. The
//! restatement is not a separate field — it falls out of the reset. Spec
//! [`psp_record_encoding.md`] §6 measured that reset at 12 % of this form's own bytes **with
//! blocks cut every 1,500 positions** (0.385 → 0.432 on tomato); the settled genomic block size
//! is 100 kb, over which the same per-block restatement amortises to a fraction of a percent.
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

use std::cmp::Ordering;

use crate::ng::psp::record::{FieldReader, RecordDecodeError, entries_to_reserve, put_varint};
use crate::pileup_record::ChainId;

// The names this module's faults are reported under, so a message says which field it was.
const DEPARTURE_COUNT: &str = "chain-id departure count";
const DEPARTURE_POSITION: &str = "chain-id departure position";
const ARRIVAL_COUNT: &str = "chain-id arrival count";
const ARRIVAL_ID: &str = "chain-id arrival";

/// The fewest bytes one departure or one arrival can take, for the count bound.
///
/// **A count is refused when no record could hold that many entries**, the same guard
/// [`FieldReader::read_count`] applies to observations: a declared 2⁶³ arrivals is not a buffer
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

    /// Build a set straight from ids that are already ascending and distinct.
    ///
    /// **Tests only, and it used to have a caller.** The writer deciding which observation it can
    /// derive built one of these out of the record's own identifiers, so that
    /// [`residual_reads`] could take a `&LiveSet`. That function takes a slice now — the writer
    /// holds those identifiers in a buffer it reuses, and handing the buffer away to build a set
    /// out of it was an allocation a record for a type check the sortedness already carries. What
    /// is left is the tests, which build a set to derive against directly.
    #[cfg(test)]
    pub(super) fn from_sorted_ids(ids: Vec<ChainId>) -> Self {
        debug_assert!(
            ids.windows(2).all(|pair| pair[0] < pair[1]),
            "a live set is ascending and without duplicates"
        );
        Self { ids }
    }

    /// Whether `id` is live. A binary search, because the set is sorted.
    #[must_use]
    pub fn contains(&self, id: ChainId) -> bool {
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

    /// Whether nothing changed — the common case at depth between two adjacent positions, and the
    /// reason this form is cheap.
    ///
    /// **Kept because a test exercises the claim.** A review found it uncalled and unchecked,
    /// which is where an inverted condition survives to its first caller; it is
    /// `an_empty_set_and_an_empty_change_say_so` that makes it a fact rather than a comment.
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

/// Fill `changes` with what it takes to get from `previously_live` to `now_live`.
///
/// Both inputs are sorted and distinct, so this is one merge pass rather than two set
/// differences: an id in `previously_live` and not in `now_live` departed, one in `now_live` and
/// not in `previously_live` arrived, and one in both did not change — which is nearly all of them.
fn derive_changes(previously_live: &LiveSet, now_live: &[ChainId], changes: &mut LiveSetChanges) {
    changes.clear();
    let (mut was, mut now) = (0usize, 0usize);
    while was < previously_live.ids.len() && now < now_live.len() {
        match previously_live.ids[was].cmp(&now_live[now]) {
            Ordering::Equal => {
                was += 1;
                now += 1;
            }
            Ordering::Less => {
                changes.departed.push(previously_live.ids[was]);
                was += 1;
            }
            // **An arrival below an id already live**, which is what a returning read looks like:
            // it sorts under everything allocated since it left. Spec `psp_record_encoding.md` §6
            // measures that at 83 % of ids on the human sample and 91 % on tomato, so this arm is
            // the ordinary case rather than the exotic one — see
            // `an_id_that_goes_live_again_below_the_ids_allocated_since_reads_back`.
            Ordering::Greater => {
                changes.arrived.push(now_live[now]);
                now += 1;
            }
        }
    }
    changes
        .departed
        .extend_from_slice(&previously_live.ids[was..]);
    changes.arrived.extend_from_slice(&now_live[now..]);
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
    // **Only the stretch from the first departure onward moves.** Everything below it keeps both
    // its value and its index, so the compaction starts there rather than at zero; and once the
    // last departure has been passed, the rest moves in one block copy rather than one predicate
    // call each. This was `Vec::retain`, which evaluates a closure over every identifier live —
    // 280 of them at 280 reads a position, to remove about two.
    let ids = &mut live.ids;
    let from = ids.partition_point(|id| *id < departed[0]);
    let len = ids.len();
    let (mut write, mut read, mut leaving) = (from, from, 0usize);
    while read < len {
        if leaving == departed.len() {
            ids.copy_within(read..len, write);
            write += len - read;
            break;
        }
        let id = ids[read];
        if departed[leaving] == id {
            leaving += 1;
        } else {
            ids[write] = id;
            write += 1;
        }
        read += 1;
    }
    ids.truncate(write);
}

/// Put `arrived` into `live`, keeping it ascending.
///
/// **In place, and it never reads what sits below the lowest arrival**: those identifiers keep
/// their value and their index whatever arrives above them. When every arrival is newer than the
/// whole set — which is what a read not seen before gives — that is the entire set, and the merge
/// becomes an append that does not read it at all.
///
/// This rebuilt the whole set into a scratch buffer and swapped it in. At 280 reads a position
/// about two identifiers arrive at every record, so 280 were read and 280 written to insert two;
/// the scratch buffer that made that possible is gone, and with it 2,264 bytes an open sample.
fn apply_arrivals(live: &mut LiveSet, arrived: &[ChainId]) {
    if arrived.is_empty() {
        return;
    }
    let ids = &mut live.ids;
    let from = ids.partition_point(|id| *id < arrived[0]);
    if from == ids.len() {
        ids.extend_from_slice(arrived);
        return;
    }
    // **The interleaving case**, reached when an arrival sorts below an identifier already live —
    // a returning read, which spec `psp_record_encoding.md` §6 measures at 83 % of identifiers on
    // the human sample and 91 % on tomato, because a chain id names a read *pair* and a pair
    // covers the reference twice.
    //
    // Merged from the back, so a value is only ever written to a slot whose own value has already
    // moved: `slot` is always `unmoved + still_to_place - 1`, one past the highest index either
    // source will still be read from, which is why the two never collide.
    let was = ids.len();
    ids.resize(was + arrived.len(), 0);
    let (mut unmoved, mut still_to_place) = (was, arrived.len());
    while still_to_place > 0 {
        let slot = unmoved + still_to_place - 1;
        if unmoved > from && ids[unmoved - 1] > arrived[still_to_place - 1] {
            ids[slot] = ids[unmoved - 1];
            unmoved -= 1;
        } else {
            ids[slot] = arrived[still_to_place - 1];
            still_to_place -= 1;
        }
    }
}

/// Write one observation's own list of reads: a count, then the identifiers as ascending gaps.
///
/// **The same shape as the arrivals half of a record's changes**, for the same reason — the ids
/// are ascending and distinct, so consecutive ones differ by at least one and a dense run costs a
/// byte each.
///
/// **This is the half of the column that carries no state**, so it lives in a record's skippable
/// body while the changes live in its head (spec `psp_record_encoding.md` §6). It is the ~3.4 % of
/// ids production stored before the 2026-08-17 ruling.
/// **`ids` need not arrive ascending or distinct, and that is the point of the check below.**
/// The bytes are gaps, so a list that is not a set cannot be written as it stands — and the
/// alternative to normalising here was a precondition on a `pub` path whose only enforcement was
/// a `debug_assert` the shipped profile removes. ⚠ Measured on that version: an observation handed
/// `[3, 3]` wrote bytes that read back as the reads `[3, 4]`, gaining a read nothing folded.
///
/// **It costs one pass and no allocation in the case that happens.** Both of ng's pileup paths
/// sort and deduplicate before handing an observation over, so the copy below is the branch
/// nothing takes.
pub(super) fn encode_read_list(ids: &[ChainId], out: &mut Vec<u8>) {
    if ids.windows(2).all(|pair| pair[0] < pair[1]) {
        write_a_read_list(ids, out);
    } else {
        let mut a_set = ids.to_vec();
        as_a_read_set(&mut a_set);
        write_a_read_list(&a_set, out);
    }
}

fn write_a_read_list(ids: &[ChainId], out: &mut Vec<u8>) {
    put_varint(out, ids.len() as u64);
    let mut previous: Option<u64> = None;
    for id in ids {
        put_varint(out, gap_from(previous, *id));
        previous = Some(*id);
    }
}

/// Read one observation's own list of reads into `into`, which is cleared first.
///
/// # Errors
///
/// [`RecordDecodeError::Truncated`] when the bytes run out, and
/// [`RecordDecodeError::Malformed`] for a count no record could hold or a gap walking past the
/// largest identifier there is.
pub(super) fn decode_read_list(
    reader: &mut FieldReader<'_>,
    field: &'static str,
    into: &mut Vec<ChainId>,
) -> Result<(), RecordDecodeError> {
    into.clear();
    let count = reader.read_count(field, LEAST_BYTES_PER_ENTRY)?;
    // **Bounded by what the remaining bytes could hold, never the declared count alone** — the
    // same guard the observation count already gets, and for the same reason: a hostile body says
    // a million reads in eleven bytes. Without the reservation the list grows from nothing, which
    // is five allocations for an observation naming seventy reads.
    into.reserve(entries_to_reserve(
        count,
        LEAST_BYTES_PER_ENTRY,
        reader.bytes_left(),
    ));
    let mut previous: Option<u64> = None;
    for _ in 0..count {
        let id = read_ascending(reader, field, "id", previous)?;
        into.push(id);
        previous = Some(id);
    }
    Ok(())
}

/// The reads the residual observation names: **the live set minus every other observation's**.
///
/// `live` is the live set's own identifiers, ascending, and `named_elsewhere` the ascending,
/// deduplicated union of the lists the record stores explicitly. Both inputs are sorted, so this
/// is one pass.
///
/// **A slice rather than a [`LiveSet`]**, because the writer's caller holds the union of a
/// record's own lists in a reused buffer and would otherwise have to give that buffer away to
/// build a set out of it — for a type check that the sortedness already carries.
///
/// **This is where most of the column's saving is and where it fails silently** (spec
/// `psp_chain_id_encoding.md` §5): derive one id too many and the reference allele gains a read
/// that does not exist, which the cohort merge composes an allele for without complaint. The
/// guard is the caller's — an observation's derived count against its own read count — and
/// `record.rs` applies it.
pub(super) fn residual_reads(
    live: &[ChainId],
    named_elsewhere: &[ChainId],
    into: &mut Vec<ChainId>,
) {
    into.clear();
    // **The exact length when the record is sound, and one allocation either way.** Every
    // identifier in `named_elsewhere` is live, so what this builds is exactly the difference of
    // the two lengths. A record where that does not hold is one the caller's guard refuses, and
    // the worst it costs here is a single growth. Without this the list grows from nothing:
    // measured at 280 reads a position, eight allocations a record where one will do.
    into.reserve(live.len().saturating_sub(named_elsewhere.len()));
    let mut elsewhere = 0usize;
    for id in live {
        while elsewhere < named_elsewhere.len() && named_elsewhere[elsewhere] < *id {
            elsewhere += 1;
        }
        if elsewhere < named_elsewhere.len() && named_elsewhere[elsewhere] == *id {
            continue;
        }
        into.push(*id);
    }
}

/// Sort `ids` ascending and drop duplicates — the shape both halves of this column want.
pub(super) fn as_a_read_set(ids: &mut Vec<ChainId>) {
    sort_and_dedup(ids);
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
    /// The ids this record names, gathered, sorted and deduplicated.
    now_live: Vec<ChainId>,
    /// This record's changes.
    changes: LiveSetChanges,
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
/// **`Default` routes through [`at_block_start`](Self::at_block_start) rather than being
/// derived**, so there is genuinely one constructor. A derived `Default` would let
/// `LiveSetWriter::new` build this state without touching the literal the `error[E0063]` above
/// depends on, and a field added here would be default-initialised on that path in silence.
#[derive(Debug)]
struct PerBlockState {
    /// The reads live at the record last handled.
    live: LiveSet,
}

impl Default for PerBlockState {
    fn default() -> Self {
        Self::at_block_start()
    }
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
    /// A writer at the start of a file's first block: nothing is live.
    ///
    /// **Equivalent to a [`start_block`](Self::start_block) that has already happened**, so the
    /// first record written restates its whole set as arrivals. Every *later* block boundary
    /// still needs the call — that is the one whose absence produces a file that parses perfectly
    /// and is wrong from that block's first record.
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
    /// `now_live` is every id the record names, in any order and with any repeats — the union of
    /// its observations' lists.
    pub fn write_changes(
        &mut self,
        now_live: impl IntoIterator<Item = ChainId>,
        out: &mut Vec<u8>,
    ) {
        self.now_live.clear();
        self.now_live.extend(now_live);
        sort_and_dedup(&mut self.now_live);
        derive_changes(&self.block.live, &self.now_live, &mut self.changes);
        encode_changes(&self.changes, &self.block.live, out);
        apply_departures(&mut self.block.live, &self.changes.departed);
        apply_arrivals(&mut self.block.live, &self.changes.arrived);
    }
}

/// **Departures first, then arrivals**, and the order is load-bearing rather than a convention.
///
/// A departure is written as its *position* in the live set, which is small where an identifier
/// is large. **Both sides resolve a position against the set as the previous record left it** —
/// this function binary-searches `live` before `write_changes` applies anything, and the reader
/// decodes the whole record before it applies anything either. So there is one set a position can
/// mean, and neither side has to remember a set that has already moved.
///
/// The order is still load-bearing for a different reason: an arrival is checked against the set
/// minus this record's departures, which the reader can only do once it has read them.
///
/// ⚠ *An earlier version of this comment said the departures were applied before the arrivals
/// were read, and that this was what made positions resolvable. Neither half was true, and the
/// eager apply was a Blocker — see `read_changes`.*
fn encode_changes(changes: &LiveSetChanges, live: &LiveSet, out: &mut Vec<u8>) {
    put_varint(out, changes.departed.len() as u64);
    let mut previous_position: Option<u64> = None;
    for id in &changes.departed {
        // PANIC-FREE: `changes.departed` comes from `derive_changes`, which pushes only ids it
        // read out of `live.ids`, and `LiveSetChanges`'s fields are private, so nothing else can
        // build one. `apply_departures` has not run yet — see `write_changes` — so `live` still
        // holds every one of them.
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
        // ⚠ **Wrapping, and deliberately not saturating — and today unreachable.** Every run
        // written through here is ascending and distinct: the arrivals and the departure
        // positions by construction, an observation's reads because [`encode_read_list`] makes
        // them so before it gets here. So no test can reach this arm, and swapping wrapping for
        // saturating fails nothing — which is reported rather than left implying otherwise.
        //
        // It is wrapping anyway, because neither a panic nor a plausible wrong byte is the right
        // answer if that ever stops being true.
        //
        // A wrapped gap is enormous, so [`value_after`]'s `checked_add` refuses it and the reader
        // names the file as damaged. **Saturating instead, which this line briefly was, turns
        // that refusal into silence**: measured, an observation handed `[3, 3]` wrote bytes that
        // read back as the reads `[3, 4]`, so it gained read 4 — a read nothing folded, which is
        // spec §5's failure reached from the writer. A reader that cannot trust the bytes must at
        // least be able to refuse them.
        Some(previous) => value.wrapping_sub(previous).wrapping_sub(1),
    }
}

/// Read a record's live-set changes, and carry the set forward.
///
/// The reading half of [`LiveSetWriter`], with the same per-block reset.
#[derive(Debug, Default)]
pub struct LiveSetReader {
    // ---- what lives for the whole file ----
    changes: LiveSetChanges,
    /// Whether `changes` holds a record that has been parsed and not yet applied.
    ///
    /// **The split exists because the caller can still refuse the record.**
    /// `read_record_head` bounds the record's body *after* the changes, and a body that stops
    /// early is a `Truncated` whose contract is *fetch more bytes and re-parse this record from
    /// its first byte* — so a set already moved would meet those bytes a second time.
    parsed_but_not_applied: bool,
    // ---- what lives for one block, and is replaced whole at every boundary ----
    block: PerBlockState,
}

impl LiveSetReader {
    /// A reader at the start of a block: nothing is live.
    ///
    /// **Equivalent to a [`start_block`](Self::start_block) that has already happened**, which is
    /// what lets a caller that seeked to an arbitrary block just build one of these and read.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Start a block: nothing is live, and the block's first record restates everything.
    pub fn start_block(&mut self) {
        self.block = PerBlockState::at_block_start();
        self.parsed_but_not_applied = false;
    }

    /// The reads live at the record last read.
    ///
    /// **Exact after a fault too**: `read_changes` moves the set only once both halves of a
    /// record's stream have been read, so a refused record leaves this as the record before it
    /// left it — which is what makes a retry from the record's first byte safe.
    #[must_use]
    pub fn live(&self) -> &LiveSet {
        &self.block.live
    }

    /// What changed at the record last read. Handed back so a caller can see the shape of the
    /// stream without diffing two sets.
    ///
    /// **Only after a [`read_changes`](Self::read_changes) that returned `Ok`.** After a fault it
    /// holds however much of the record was decoded before the fault, which is not a record's
    /// changes — spec §6.7's *"none of these may reach a caller as a half-built record"* is about
    /// exactly this.
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
    ///
    /// # Errors
    ///
    /// - [`RecordDecodeError::Truncated`] when the bytes run out mid-stream — **fetch more and
    ///   call this again with the record's whole stream.** The set is untouched, so the retry is
    ///   a retry.
    /// - [`RecordDecodeError::Malformed`] when the stream cannot mean what it says: a departure
    ///   count larger than the live set, a departure position past its end, an arrival for a read
    ///   already live, or a gap walking past the largest identifier there is. **The block is
    ///   damaged**; no quantity of further bytes changes the answer.
    ///
    /// [`encode_changes`] is the mirror of this function, and the field order the two agree on is
    /// written there.
    pub fn read_changes(&mut self, bytes: &[u8]) -> Result<usize, RecordDecodeError> {
        let used = self.parse_changes(bytes)?;
        self.apply_the_changes_just_parsed();
        Ok(used)
    }

    /// Read one record's changes off the front of `bytes` **without moving the live set**,
    /// returning how many bytes they took.
    ///
    /// **For a caller that can still refuse the record after this point**, which is
    /// `read_record_head`: it bounds the record's body *after* the changes, and a body that stops
    /// early is a `Truncated` answered by fetching more bytes and re-parsing the record from its
    /// first byte. A set already moved would meet those bytes a second time — an arrival into a
    /// set it is already in, refused as damage on a perfectly good file, or worse, a departure
    /// position resolved against an already-shrunk set, which removes a different read and says
    /// nothing.
    pub fn parse_changes(&mut self, bytes: &[u8]) -> Result<usize, RecordDecodeError> {
        let reader = &mut FieldReader::new(bytes);
        self.changes.clear();
        self.parsed_but_not_applied = false;

        self.read_departures(reader)?;
        self.read_arrivals(reader)?;
        self.parsed_but_not_applied = true;
        Ok(reader.bytes_read())
    }

    /// Move the live set by the changes [`parse_changes`](Self::parse_changes) last read.
    ///
    /// **A no-op unless a parse is waiting**, so a caller that applies twice moves the set once,
    /// and a caller that refused the record between the two applies nothing.
    ///
    /// **Defensive, and measured to be so**: removing the guard passes every test in the module,
    /// because no caller applies twice today. It is here because `read_changes` parses *and*
    /// applies while `read_record_head` calls the two separately, so the two spellings sit beside
    /// each other and the mistake is one line away — the same reason `begin_next_block` resets
    /// zstd's decoder.
    pub fn apply_the_changes_just_parsed(&mut self) {
        if !self.parsed_but_not_applied {
            return;
        }
        self.parsed_but_not_applied = false;

        // **The set moves here and nowhere else.** Spec §8: *"a parse that half-advances that
        // state before failing corrupts every record after it, plausibly."*
        //
        // ⚠ **Twice now that has been got wrong, one level apart.** `read_changes` used to apply
        // the departures between reading them and reading the arrivals, so a stream cut in the
        // arrivals half left the set shrunk — five of six cut points retried to `Ok` with a read
        // silently gone. And then `read_record_head` applied the whole thing before it bounded
        // the record's body, so a body that stopped early did the same: measured end to end, a
        // well-formed file of 1,999 records was refused at record 149 as *"already live"*, and a
        // record that only departs reads retried to `Ok` with `[1, 2]` where the truth was
        // `[1, 2, 4]`.
        apply_departures(&mut self.block.live, &self.changes.departed);
        apply_arrivals(&mut self.block.live, &self.changes.arrived);
    }

    /// The departures, each a position in the live set as it stands, resolved back to the
    /// identifier it names.
    ///
    /// **A position past the end of the set is damage.** It is the one place a `usize` from the
    /// wire indexes something, and an unchecked index here would be a panic on a corrupt file
    /// where spec §6.7 requires a refusal.
    fn read_departures(&mut self, reader: &mut FieldReader<'_>) -> Result<(), RecordDecodeError> {
        let departures = reader.read_count(DEPARTURE_COUNT, LEAST_BYTES_PER_ENTRY)?;
        // **A departure is a position in the live set**, and the positions are strictly ascending,
        // so there cannot be more of them than the set holds. That is exact where the count bound
        // above is the record body's byte ceiling borrowed from a different container — and it
        // names the fault at the count rather than at whichever position first runs off the end.
        if departures > self.block.live.len() as u64 {
            return Err(reader.malformed(
                DEPARTURE_COUNT,
                format!(
                    "{departures} departures from a live set holding {} reads",
                    self.block.live.len()
                ),
            ));
        }
        let mut previous: Option<u64> = None;
        for _ in 0..departures {
            let position = read_ascending(reader, DEPARTURE_POSITION, "position", previous)?;
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
            previous = Some(position);
        }
        Ok(())
    }

    /// The arrivals, each an identifier.
    ///
    /// **An identifier already live is damage**: accepted, it would put a read in the set twice,
    /// which Milestone E4's residual arithmetic then subtracts once — silently.
    fn read_arrivals(&mut self, reader: &mut FieldReader<'_>) -> Result<(), RecordDecodeError> {
        let arrivals = reader.read_count(ARRIVAL_COUNT, LEAST_BYTES_PER_ENTRY)?;
        let mut previous: Option<u64> = None;
        for _ in 0..arrivals {
            let id = read_ascending(reader, ARRIVAL_ID, "id", previous)?;
            // **The set has not moved yet**, so this asks the set as the *previous* record left
            // it — which also refuses a record that departs an id and arrives the same id.
            // `derive_changes` puts an id in at most one of the two lists, so no writer produces
            // that; a stream that does is internally inconsistent, and honouring it would report
            // one departure and one arrival for a record where nothing changed, which is what
            // Milestone E4's residual arithmetic counts.
            if self.block.live.contains(id) {
                return Err(reader.malformed(
                    ARRIVAL_ID,
                    format!(
                        "id {id}, which is already live — an arrival names a read that was not"
                    ),
                ));
            }
            self.changes.arrived.push(id);
            previous = Some(id);
        }
        Ok(())
    }
}

/// The next value of a strictly ascending run, read as its gap from the one before.
///
/// **One helper for both runs**, because the departures and the arrivals differ only in what the
/// numbers mean — `noun` is what the message calls them. Written twice, the two would be free to
/// drift in exactly the place where identifiers and positions must not be confused.
fn read_ascending(
    reader: &mut FieldReader<'_>,
    field: &'static str,
    noun: &str,
    previous: Option<u64>,
) -> Result<u64, RecordDecodeError> {
    let gap = reader.read_varint(field)?;
    value_after(previous, gap).ok_or_else(|| {
        reader.malformed(
            field,
            format!(
                "a gap of {gap} past {noun} {}, which is past every {noun} there is",
                previous.unwrap_or(0)
            ),
        )
    })
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

    use std::collections::{BTreeMap, BTreeSet};

    use proptest::prelude::*;

    /// Write a whole file's worth of records and read them back, restarting at every block.
    ///
    /// `blocks` is one entry per block, each a list of records, each record the ids it names —
    /// in whatever order and with whatever repeats the caller gives, because that is the shape a
    /// record's observations hand over.
    /// **What one write-and-read gives back**, so a test that needs the arrival count does not
    /// write the file a second time to get it. Four tests did, which is four places the block
    /// restart could be forgotten independently — the one mistake this module's tests exist to
    /// catch.
    struct WhatCameBack {
        /// The live set after each record, in file order.
        sets: Vec<Vec<ChainId>>,
        /// How many arrivals the reader met over the whole file. **One a stretch, not one an
        /// identifier** — which is what a stream assuming one stretch per read gets wrong.
        arrivals: usize,
        /// What the whole file's changes cost on the wire.
        bytes: usize,
    }

    fn round_trip(blocks: &[Vec<Vec<ChainId>>]) -> WhatCameBack {
        let written = write_a_file(blocks);

        let mut reader = LiveSetReader::new();
        let mut came_back = WhatCameBack {
            sets: Vec::new(),
            arrivals: 0,
            bytes: written.iter().map(Vec::len).sum(),
        };
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
                came_back.arrivals += reader.changes().arrived().len();
                came_back.sets.push(reader.live().ids().to_vec());
                at += 1;
            }
        }
        came_back
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

    /// Coverage that looks like **read pairs**: pairs begin at every position, each mate a given
    /// number of bases long with an unsequenced gap between them, and one identifier names the
    /// whole pair.
    ///
    /// **The shape is a struct rather than four `u64`s**, because the mate length and the gap
    /// between the mates are different biology and the fixture cannot tell them apart: measured,
    /// transposing them gives the same 800 identifiers, the same 660 covering two stretches and
    /// the same 1,460 stretches, because a mate starts at the sum of the two either way. Named at
    /// the call site, a later edit that swaps them is visible.
    ///
    /// **So an identifier covers two stretches with a hole between them**, which is what makes a
    /// read go live, stop, and go live again. Spec `psp_record_encoding.md` §6 measures that on
    /// real alignments: 83 % of identifiers on the human sample and 91 % on tomato.
    struct PairedEndFixture {
        positions: u64,
        pairs_starting_at_each_position: u64,
        mate_length: u64,
        unsequenced_gap_between_mates: u64,
    }

    fn a_paired_end_coverage(shape: &PairedEndFixture) -> Vec<Vec<ChainId>> {
        let PairedEndFixture {
            positions,
            pairs_starting_at_each_position,
            mate_length,
            unsequenced_gap_between_mates,
        } = *shape;
        let pair_span = mate_length * 2 + unsequenced_gap_between_mates;
        (0..positions)
            .map(|at| {
                let earliest = at.saturating_sub(pair_span.saturating_sub(1));
                (earliest..=at)
                    .filter(|start| {
                        let into_the_pair = at - start;
                        into_the_pair < mate_length
                            || (mate_length + unsequenced_gap_between_mates..pair_span)
                                .contains(&into_the_pair)
                    })
                    .flat_map(|start| {
                        (0..pairs_starting_at_each_position)
                            .map(move |which| start * pairs_starting_at_each_position + which)
                    })
                    .collect()
            })
            .collect()
    }

    /// How many separate stretches of records each identifier covers — one for a read that is
    /// live over a run and then gone, two for a pair whose mates do not overlap.
    ///
    /// **The oracle's own measure**, so a fixture that quietly stopped containing re-entry shows
    /// up as a number rather than as a test that still passes.
    fn stretches_of(records: &[Vec<ChainId>]) -> BTreeMap<ChainId, usize> {
        let mut stretches: BTreeMap<ChainId, usize> = BTreeMap::new();
        let mut previously_live: BTreeSet<ChainId> = BTreeSet::new();
        for record in records {
            let now_live: BTreeSet<ChainId> = record.iter().copied().collect();
            for id in now_live.difference(&previously_live) {
                *stretches.entry(*id).or_default() += 1;
            }
            previously_live = now_live;
        }
        stretches
    }

    #[test]
    fn a_run_of_records_reads_back_the_set_each_one_named() {
        let records = a_sliding_coverage(200, 3, 40);
        let wanted: Vec<_> = records.iter().map(|ids| as_a_set(ids)).collect();
        let back = round_trip(&[records]).sets;
        assert_eq!(back, wanted);
    }

    /// **The ids a record names arrive unsorted and with repeats**, because they are the union of
    /// its observations' lists and one read pair can be named by two of them. The set is the
    /// writer's to build; a caller sorting first would be a precondition nothing enforces.
    #[test]
    fn unsorted_and_repeated_ids_name_the_same_set() {
        let records = vec![vec![900, 4, 4, 71, 4, 900], vec![71, 900, 71]];
        let back = round_trip(&[records]).sets;
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

        let whole = round_trip(&blocks).sets;

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
        // The same set, reached from inside a block. ⚠ This was written into a throwaway `Vec`,
        // so the contrast the comment below claims was never made.
        let mut inside_a_block = Vec::new();
        writer.write_changes([9u64, 12, 40], &mut inside_a_block);

        // A second block whose first record is that same set: it must restate all three.
        let mut restated = Vec::new();
        writer.start_block();
        writer.write_changes([9u64, 12, 40], &mut restated);

        reader.start_block();
        reader.read_changes(&restated).expect("it reads");
        assert_eq!(reader.changes().arrived(), [9, 12, 40]);
        assert!(reader.changes().departed().is_empty());
        assert_eq!(reader.live().ids(), [9, 12, 40]);

        // **And the contrast, as bytes.** Restating {9, 12, 40} is a count of three and three
        // gaps; stepping to it from {5, 9, 12} is one departure at position 0 and one arrival,
        // and it never names 9 or 12 at all.
        assert_eq!(restated, [0, 3, 9, 2, 27], "no departures, three arrivals");
        assert_eq!(
            inside_a_block,
            [1, 0, 1, 40],
            "one departure at position 0, and 40 arriving"
        );
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
            "a departure spelled as its identifier would cost four bytes where the position \
             costs one"
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
        let changes_bytes = round_trip(std::slice::from_ref(&records)).bytes;

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

        // What an uninterrupted read gives, to compare every interrupted one against.
        let mut uninterrupted = LiveSetReader::new();
        uninterrupted.start_block();
        uninterrupted
            .read_changes(&opening_bytes)
            .expect("the opening record reads");
        uninterrupted.read_changes(&bytes).expect("and the next");
        let uninterrupted = uninterrupted.live().ids().to_vec();

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

            // **And the retry the class instructs gives the uninterrupted answer**, on the same
            // reader. ⚠ This half was missing, and its absence is why a Blocker survived: with a
            // fresh reader per cut the sweep checks the fault's *class* and never the state the
            // fault leaves behind, which is the only thing the class split exists to protect.
            reader
                .read_changes(&bytes)
                .unwrap_or_else(|refused| panic!("the retry after a cut at {cut}: {refused}"));
            assert_eq!(
                reader.live().ids(),
                uninterrupted.as_slice(),
                "retrying after a cut at {cut} of {} gave a different set",
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

    // -----------------------------------------------------------------
    // Re-entry: a read goes live, stops, and goes live again
    // -----------------------------------------------------------------

    /// The share of identifiers that must cover two stretches for this fixture to be saying
    /// anything about re-entry at all. **The fixture measures 82.5 %**; the corpora spec
    /// `psp_record_encoding.md` §6 gives are 83 % on the human sample and 91 % on tomato. Written
    /// out because the number is what a later shortening of the fixture would erode, and an
    /// expression like `× 10 > × 8` does not say how much headroom was meant.
    const LEAST_TWO_STRETCH_SHARE: f64 = 0.80;

    /// **Most reads go live twice, and every one of them comes back.**
    ///
    /// A chain id names a read *pair*, with the mates collapsed onto one identifier, and a pair's
    /// mates rarely overlap — so the id covers two stretches with an unsequenced hole between
    /// them. Spec `psp_record_encoding.md` §6 measures that at **83 % of identifiers on the human
    /// sample and 91 % on tomato**, and says what a stream that assumed one stretch per id would
    /// do: *"loses the second mate of nine reads in ten — silently, because the merge would simply
    /// see a read that was not there."*
    ///
    /// So the assertion is not that the walk came out right. It is that **the fixture contains the
    /// thing**, counted, and that the thing survives:
    ///
    /// - how many identifiers cover two stretches rather than one, as a fraction, so a fixture
    ///   that quietly stopped containing re-entry fails here rather than passing quietly;
    /// - that the writer emitted **one arrival per stretch**, not one per identifier — which is
    ///   exactly what a writer that assumed one stretch per id would get wrong, and what nothing
    ///   downstream could then recover;
    /// - and that every record's set reads back, so a second stretch that was written is also
    ///   read.
    #[test]
    fn most_reads_go_live_twice_and_every_one_of_them_comes_back() {
        // Measured on this fixture: 800 identifiers over 400 records, 660 of them (82.5 %)
        // covering two stretches and 1,460 stretches in all — so a stream naming each identifier
        // once would carry 800 arrivals and lose 660 second mates. The 140 that cover one stretch
        // are the pairs whose second mate starts past the fixture's last record.
        let records = a_paired_end_coverage(&PairedEndFixture {
            positions: 400,
            pairs_starting_at_each_position: 2,
            mate_length: 30,
            unsequenced_gap_between_mates: 40,
        });
        let stretches = stretches_of(&records);

        // **The mate length is measured, not only stated.** Transposing it with the gap gives the
        // same identifier and stretch counts, so nothing else in this test would notice.
        assert!(
            records[29].contains(&0) && !records[30].contains(&0),
            "a mate is 30 bases, so the pair that started at record 0 covers records 0 to 29"
        );

        let went_live_twice = stretches.values().filter(|count| **count == 2).count();
        let identifiers = stretches.len();
        // **The three numbers the comment above states, asserted.** The fraction below is the
        // statement of intent and catches a gross change; these catch drift, and they close a
        // circularity — `stretches_of` and `derive_changes` are two implementations of the same
        // set difference, so the arrival equality further down is a differential between them
        // rather than a check against a number derived without the code.
        assert_eq!(
            identifiers, 800,
            "two pairs starting at each of 400 positions"
        );
        assert_eq!(
            went_live_twice, 660,
            "the pairs whose second mate fits inside the fixture"
        );
        assert!(
            went_live_twice as f64 > identifiers as f64 * LEAST_TWO_STRETCH_SHARE,
            "the fixture has to contain re-entry to say anything about it: {went_live_twice} of \
             {identifiers} identifiers cover two stretches, and this asks for more than {}%",
            LEAST_TWO_STRETCH_SHARE * 100.0
        );
        assert!(
            stretches.values().all(|count| *count <= 2),
            "a pair has two mates, so no identifier should reach three stretches here"
        );

        // **One arrival per stretch, not one per identifier.** This is the number a writer that
        // assumed one stretch per id would get wrong, and it is measured on the bytes rather than
        // inferred: the reader reports what it read.
        let stretches_in_all: usize = stretches.values().sum();
        assert_eq!(
            stretches_in_all, 1_460,
            "800 first mates and 660 second ones"
        );
        let came_back = round_trip(std::slice::from_ref(&records));
        assert_eq!(
            came_back.arrivals, stretches_in_all,
            "one arrival a stretch: {identifiers} identifiers cover {stretches_in_all} stretches, \
             and a stream that named each identifier once would carry {identifiers}"
        );

        // And every record's set, so a second stretch that was written is also read.
        let wanted: Vec<_> = records.iter().map(|ids| as_a_set(ids)).collect();
        assert_eq!(came_back.sets, wanted);
    }

    /// **A read that comes back after a block boundary is restated, not remembered.**
    ///
    /// The two rules meet here: the live set restarts at every block (spec §3.2), and a read may
    /// go live twice. **Two shapes cross the boundary and both come out as ordinary arrivals** —
    /// a read that never stopped covering, which is *restated*, and a read whose second mate
    /// begins after the boundary, which arrives there with nothing saying it is the same read.
    /// The reader carries nothing across, so it cannot tell them apart, and it does not need to.
    ///
    /// ⚠ **The fixture had no read live across the boundary when this test was written**, so the
    /// restatement half of that claim was not exercised — and making `start_block` a no-op on the
    /// writing side left this test green while two older ones caught it. Id 7 below is the read
    /// the docstring was talking about.
    #[test]
    fn a_read_that_comes_back_across_a_block_boundary_arrives_again() {
        // Id 4 is a pair whose mates cover records 0-1 and 4-5, so its second mate begins after
        // the cut. Id 7 is a read still covering *across* the cut: live from record 0 to 6.
        let first_block = vec![vec![4u64, 7], vec![4u64, 7], vec![7u64], vec![7u64]];
        let second_block = vec![vec![4u64, 7], vec![4u64, 7], vec![7u64]];

        let written = write_a_file(&[first_block.clone(), second_block.clone()]);
        let mut reader = LiveSetReader::new();

        reader.start_block();
        let mut arrivals_in_the_first_block = 0usize;
        for bytes in &written[..first_block.len()] {
            reader.read_changes(bytes).expect("it reads");
            arrivals_in_the_first_block += reader.changes().arrived().len();
        }
        assert_eq!(
            arrivals_in_the_first_block, 2,
            "two reads, one arrival each"
        );
        assert_eq!(
            reader.live().ids(),
            [7],
            "the pair's first mate has gone by the block's end, and the other read has not"
        );

        // **The same reader, carried across the boundary.** ⚠ This half is what makes the test
        // able to see a *reader* that forgot to reset: a fresh reader has no state to carry, so
        // the alternative below cannot fail on that at all.
        reader.start_block();
        let mut carried_on = Vec::new();
        let mut arrivals_at_the_new_block = Vec::new();
        for bytes in &written[first_block.len()..] {
            reader.read_changes(bytes).expect("it reads");
            arrivals_at_the_new_block.push(reader.changes().arrived().to_vec());
            carried_on.push(reader.live().ids().to_vec());
        }
        assert_eq!(
            arrivals_at_the_new_block[0],
            vec![4, 7],
            "id 7 is restated even though it never stopped covering, and id 4's second mate \
             arrives beside it with nothing saying it is the same read"
        );

        // And a reader that starts *here*, with no history at all — Milestone F's index will hand
        // one exactly this. It must agree with the reader that came through the first block.
        let mut reader_starting_fresh = LiveSetReader::new();
        reader_starting_fresh.start_block();
        let mut read_alone = Vec::new();
        for bytes in &written[first_block.len()..] {
            reader_starting_fresh.read_changes(bytes).expect("it reads");
            read_alone.push(reader_starting_fresh.live().ids().to_vec());
        }
        assert_eq!(
            carried_on, read_alone,
            "the block reads the same either way"
        );
        assert_eq!(carried_on, vec![vec![4, 7], vec![4, 7], vec![7]]);
    }

    /// **An identifier may go live more than twice**, and nothing here counts.
    ///
    /// Two mates is what a pair gives, but **a spliced alignment can leave more than one hole**,
    /// and nothing in the codec caps the count: an identifier that departed and arrives again is
    /// an ordinary arrival however many times it has done so. This is the test that says so.
    ///
    /// *⚠ An earlier version of this sentence also offered "a read straddling two walked regions
    /// is named twice". That argues the other way: under an allocator that never reuses an id,
    /// being named twice gives **two identifiers of one stretch each**, not one identifier of
    /// two.*
    #[test]
    fn an_identifier_that_goes_live_five_times_reads_back() {
        let records: Vec<Vec<ChainId>> = (0..20u64)
            .map(|at| if at % 4 < 2 { vec![9u64] } else { Vec::new() })
            .collect();
        let stretches = stretches_of(&records);
        assert_eq!(
            stretches.get(&9),
            Some(&5),
            "the fixture must give this identifier five stretches"
        );

        let back = round_trip(std::slice::from_ref(&records)).sets;
        assert_eq!(back, records);
        // ⚠ **What carries this test is the fixture assertion above, not the round trip.** Only
        // one read is ever live, so its departure is always position 0 and its arrival always id
        // 9 — a codec that confused the two would pass here. The five-stretch shape is the point.
    }

    /// **An id that goes live, stops, and goes live again while larger ids are live.**
    ///
    /// A pair's mates rarely overlap, so most identifiers cover two stretches — spec
    /// `psp_record_encoding.md` §6 measures 83 % of them on the human sample and 91 % on tomato.
    /// **The returning id sorts below every id allocated since it left**, and that is the only
    /// shape that makes either merge in this module interleave: `derive_changes`'s `Greater` arm
    /// and `apply_arrivals`'s `else`.
    ///
    /// ⚠ **Both arms were unreached by all 214 `ng::psp` tests** before this one — measured by
    /// replacing each with `panic!` and watching the suite stay green. Every other fixture
    /// allocates ids in walk order, so an arrival was always the largest. Counting two-stretch
    /// ids on a real corpus is Milestone E2's; making the arms reachable at all is this step's,
    /// because they are E1's code.
    #[test]
    fn an_id_that_goes_live_again_below_the_ids_allocated_since_reads_back() {
        let records = vec![
            vec![10u64, 20, 30],
            vec![20u64, 30],
            vec![20u64, 30, 40, 50],
            // 10 comes back, below 20, 30, 40 and 50.
            vec![10u64, 20, 30, 40, 50],
            // and 15 arrives between two live ids while three others leave.
            vec![10u64, 15, 50],
        ];
        let back = round_trip(std::slice::from_ref(&records)).sets;
        assert_eq!(
            back,
            vec![
                vec![10, 20, 30],
                vec![20, 30],
                vec![20, 30, 40, 50],
                vec![10, 20, 30, 40, 50],
                vec![10, 15, 50],
            ]
        );
    }

    /// **The stream reports the bytes it used, not the buffer's length**, because Milestone E3
    /// hands it the rest of the record head and the count is how the next field is found.
    ///
    /// ⚠ Nothing could see this before: every fixture handed `read_changes` a `Vec` holding
    /// exactly one record's stream, so returning `bytes.len()` instead was the same number.
    #[test]
    fn read_changes_reports_only_the_bytes_it_used() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();
        let mut bytes = Vec::new();
        writer.write_changes([1u64, 2, 3], &mut bytes);
        let its_own = bytes.len();
        bytes.extend_from_slice(&[0xAA; 8]);

        let mut reader = LiveSetReader::new();
        reader.start_block();
        assert_eq!(reader.read_changes(&bytes).expect("it reads"), its_own);
        assert_eq!(
            reader.live().ids(),
            [1, 2, 3],
            "and whatever follows the stream is not read as arrivals"
        );
    }

    /// **What the two halves say when nothing is happening.** A gap in coverage empties the set; a
    /// record between two positions of the same reads changes nothing, which is the common case at
    /// depth and the reason this form is cheap.
    #[test]
    fn an_empty_set_and_an_empty_change_say_so() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();
        assert!(writer.live().is_empty(), "a block opens with nothing live");

        let mut bytes = Vec::new();
        writer.write_changes([7u64, 8], &mut bytes);
        assert!(!writer.live().is_empty());
        assert_eq!(
            writer.live().ids(),
            [7, 8],
            "the writer's set is the record just written, not the one before it"
        );

        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&bytes).expect("it reads");
        assert!(!reader.changes().is_empty(), "two reads arrived");

        let mut unchanged = Vec::new();
        writer.write_changes([7u64, 8], &mut unchanged);
        reader.read_changes(&unchanged).expect("it reads");
        assert!(
            reader.changes().is_empty(),
            "the cheap case: nothing arrived, nothing left"
        );

        let mut gone = Vec::new();
        writer.write_changes([], &mut gone);
        reader.read_changes(&gone).expect("it reads");
        assert!(reader.live().is_empty(), "coverage stopped");
    }

    /// **`departed()` names the reads that stopped covering**, and until this test the accessor
    /// could have returned anything.
    ///
    /// ⚠ Measured: replacing its body with `&[]` left all 228 `ng::psp` tests green, because the
    /// only test that read it asserted the list was *empty*. Milestone E4's residual arithmetic is
    /// specified against these two counts.
    #[test]
    fn changes_departed_names_the_reads_that_stopped_covering() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();
        writer.write_changes([7u64, 8, 9], &mut Vec::new());
        let mut bytes = Vec::new();
        writer.write_changes([9u64, 40], &mut bytes);

        let written = write_a_file(&[vec![vec![7u64, 8, 9], vec![9u64, 40]]]);
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&written[0]).expect("it reads");
        assert_eq!(written[1], bytes, "the same record, written the same way");
        reader.read_changes(&written[1]).expect("it reads");
        assert_eq!(
            reader.changes().departed(),
            [7, 8],
            "7 and 8 stopped covering"
        );
        assert_eq!(reader.changes().arrived(), [40], "and 40 started");
        assert!(!reader.changes().is_empty());
    }

    /// **The departure count's own boundary.** One departure per live read is the most a record
    /// can have, and one more than that is damage — named at the count, where the offset points at
    /// the right byte, rather than at whichever position first runs off the end.
    ///
    /// ⚠ The guard test beside this one uses nine departures against a set of two, which is far
    /// from the boundary: relaxing the bound by one left all 228 tests green.
    #[test]
    fn the_departure_count_may_equal_the_live_set_and_no_more() {
        // Two reads live, then a record that departs both — the most a record can depart.
        let mut all_of_them = LiveSetReader::new();
        all_of_them.start_block();
        all_of_them.read_changes(&[0, 2, 4, 0]).expect("two arrive");
        assert_eq!(all_of_them.live().ids(), [4, 5]);
        all_of_them
            .read_changes(&[2, 0, 0, 0])
            .expect("a whole set may depart at once");
        assert!(all_of_them.live().is_empty());

        // Three, one past the set: damage, and named at the count.
        let mut one_too_many = LiveSetReader::new();
        one_too_many.start_block();
        one_too_many
            .read_changes(&[0, 2, 4, 0])
            .expect("two arrive");
        let refused = one_too_many
            .read_changes(&[3, 0, 0, 0, 0])
            .expect_err("three departures from a set of two");
        assert!(
            refused.to_string().contains("3 departures"),
            "the message must name the count rather than a position: {refused}"
        );
    }

    /// **A departure count larger than the live set is damage at the count**, not at whichever
    /// position first happens to run off the end — which is where the record-body ceiling alone
    /// would have reported it.
    #[test]
    fn more_departures_than_the_live_set_holds_is_damage() {
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&[0, 2, 4, 0]).expect("two arrive");

        let refused = reader
            .read_changes(&[9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            .expect_err("nine departures from a set of two");
        assert!(
            matches!(refused, RecordDecodeError::Malformed { .. }),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("9 departures"),
            "the message must name the count: {refused}"
        );
    }

    /// **The gap is biased by one on both halves**, and only a run of more than one entry shows
    /// it: with a run of one, the first value is written absolutely and the bias never applies.
    /// Every other byte-exact test in this module holds a run of one.
    #[test]
    fn a_run_of_arrivals_costs_its_gaps_biased_by_one() {
        let mut writer = LiveSetWriter::new();
        writer.start_block();

        // Three arrivals: 4, then 5 and 6 adjacent — gaps of nothing if the bias is there, of one
        // if it is not.
        let mut opening = Vec::new();
        writer.write_changes([4u64, 5, 6], &mut opening);
        assert_eq!(
            opening,
            [0, 3, 4, 0, 0],
            "no departures, three arrivals: 4 absolutely, then two gaps of nothing"
        );

        // Two departures at adjacent positions 0 and 1 of [4, 5, 6].
        let mut departing = Vec::new();
        writer.write_changes([6u64], &mut departing);
        assert_eq!(
            departing,
            [2, 0, 0, 0],
            "two departures at positions 0 and 1, then no arrivals"
        );
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(600))]

        /// **A `Truncated` fault must leave the reader exactly where it was**, because
        /// Milestone D's reader answers it by re-parsing the record from its first byte (spec
        /// `psp_file_format.md` §8).
        ///
        /// ⚠ **This found the Blocker on its first generated case** — `records = [[0], [1, 2]]` —
        /// where thirteen hand-written tests could not, because the sweep over cuts built a fresh
        /// reader each time and so never retried.
        #[test]
        fn a_truncated_read_retried_with_the_whole_record_gives_the_uninterrupted_answer(
            records in prop::collection::vec(prop::collection::vec(0u64..3_000, 0..8), 2..8),
            cut_seed in any::<usize>(),
        ) {
            let mut writer = LiveSetWriter::new();
            writer.start_block();
            let mut per_record = Vec::new();
            for record in &records {
                let mut bytes = Vec::new();
                writer.write_changes(record.iter().copied(), &mut bytes);
                per_record.push(bytes);
            }
            let last = per_record.len() - 1;
            if per_record[last].len() >= 2 {
                let cut = 1 + cut_seed % (per_record[last].len() - 1);

                let mut uninterrupted = LiveSetReader::new();
                uninterrupted.start_block();
                for bytes in &per_record {
                    uninterrupted
                        .read_changes(bytes)
                        .map_err(|refused| TestCaseError::fail(refused.to_string()))?;
                }
                let uninterrupted = uninterrupted.live().ids().to_vec();

                let mut retried = LiveSetReader::new();
                retried.start_block();
                for bytes in &per_record[..last] {
                    retried
                        .read_changes(bytes)
                        .map_err(|refused| TestCaseError::fail(refused.to_string()))?;
                }
                let fault = retried.read_changes(&per_record[last][..cut]);
                prop_assert!(
                    matches!(fault, Err(RecordDecodeError::Truncated { .. })),
                    "a strict prefix gave {fault:?}"
                );
                let after = retried.read_changes(&per_record[last]);
                prop_assert!(after.is_ok(), "the retry refused a good record: {after:?}");
                prop_assert_eq!(retried.live().ids().to_vec(), uninterrupted);
            }
        }

        /// **Any run of records round trips through one buffer**, over the whole `u64` range
        /// rather than identifiers allocated in order — and read back by advancing on the byte
        /// count the reader returns, which is what Milestone E3 will do.
        #[test]
        fn any_run_of_records_round_trips_through_one_buffer(
            records in prop::collection::vec(prop::collection::vec(any::<u64>(), 0..8), 1..25),
        ) {
            let mut writer = LiveSetWriter::new();
            writer.start_block();
            let mut stream = Vec::new();
            let mut wanted = Vec::new();
            for record in &records {
                writer.write_changes(record.iter().copied(), &mut stream);
                wanted.push(as_a_set(record));
            }

            let mut reader = LiveSetReader::new();
            reader.start_block();
            let mut at = 0usize;
            let mut back = Vec::new();
            for _ in &records {
                let used = reader
                    .read_changes(&stream[at..])
                    .map_err(|refused| TestCaseError::fail(format!("record at {at}: {refused}")))?;
                at += used;
                back.push(reader.live().ids().to_vec());
            }
            prop_assert_eq!(back, wanted);
            prop_assert_eq!(at, stream.len());
        }

        /// **Any pattern of comings and goings round trips**, however many times a read leaves
        /// and returns.
        ///
        /// Each generated read is a mask over the same records — live where the mask says so — so
        /// the generator reaches re-entry shapes a pair does not make: three stretches, four, a
        /// read live only at the last record, a read that never goes live at all.
        ///
        /// **`prop_assume!` throws away the draws with no re-entry in them.** Without it the
        /// count assertion below is self-consistent — both sides are computed from the same
        /// records, so it says the reader met one arrival per stretch and nothing about how many
        /// stretches there are, and a draw with no re-entry at all satisfies it exactly. ⚠ The
        /// first version of this test claimed the opposite in its own docstring.
        #[test]
        fn any_pattern_of_comings_and_goings_round_trips(
            (records_in_the_run, masks) in (2usize..30).prop_flat_map(|records_in_the_run| {
                (
                    Just(records_in_the_run),
                    prop::collection::vec(
                        prop::collection::vec(any::<bool>(), records_in_the_run),
                        1..12,
                    ),
                )
            }),
            records_in_the_first_block in 1usize..29,
        ) {
            let records: Vec<Vec<ChainId>> = (0..records_in_the_run)
                .map(|at| {
                    masks
                        .iter()
                        .enumerate()
                        .filter(|(_, mask)| mask[at])
                        .map(|(read, _)| read as ChainId)
                        .collect()
                })
                .collect();

            let stretches = stretches_of(&records);
            prop_assume!(stretches.values().any(|count| *count > 1));

            // **Cut into two blocks**, at a point the generator chose, so a read that is still
            // covering at the cut is restated and one whose next stretch begins after it arrives
            // fresh. Milestone E3 stresses exactly this seam, and two of the three hand-written
            // re-entry tests use one block.
            let cut = records_in_the_first_block.min(records.len() - 1);
            let blocks = vec![records[..cut].to_vec(), records[cut..].to_vec()];

            let came_back = round_trip(&blocks);
            prop_assert_eq!(&came_back.sets, &records);

            // The reader met every stretch of every read, and no more — **plus one restatement
            // for every read still covering at the cut**, which is what a block costs.
            let still_covering_at_the_cut = records[cut - 1]
                .iter()
                .filter(|id| records[cut].contains(id))
                .count();
            prop_assert_eq!(
                came_back.arrivals,
                stretches.values().sum::<usize>() + still_covering_at_the_cut
            );
        }

        /// **Damaged streams are refused, or leave the live set a set** — ascending and without
        /// duplicates. A decode that returns `Ok` having put an id in twice or out of order is
        /// what silently breaks Milestone E4's residual subtraction.
        ///
        /// **The bytes are a real record's, damaged, rather than uniform noise**, and the reader
        /// is seeded with a live set first. Both matter, and the difference is measured: over 600
        /// cases **113 of the damaged streams are accepted** and have their set checked, where
        /// uniform bytes are refused at the first count almost every time.
        ///
        /// ⚠ **What it still does not reach is the interleaving arm** of `apply_arrivals` — the
        /// shape that can actually break the ordering, where an arrival sorts below an id already
        /// live. Measured on both versions of this test: zero entries over 600 cases, against
        /// 1,320 by `most_reads_go_live_twice_and_every_one_of_them_comes_back`. **That arm's
        /// guards are that oracle and `any_pattern_of_comings_and_goings_round_trips`**, not this
        /// test; what this one holds is that damage is refused or harmless.
        #[test]
        fn damaged_streams_are_refused_or_keep_the_live_set_a_sorted_set(
            already_live in prop::collection::vec(0u64..300, 0..40),
            then_live in prop::collection::vec(0u64..300, 0..40),
            damage_at in any::<usize>(),
            damage_to in any::<u8>(),
        ) {
            let mut writer = LiveSetWriter::new();
            writer.start_block();
            let mut opening = Vec::new();
            writer.write_changes(already_live.iter().copied(), &mut opening);
            let mut next = Vec::new();
            writer.write_changes(then_live.iter().copied(), &mut next);
            prop_assume!(!next.is_empty());

            let mut reader = LiveSetReader::new();
            reader.start_block();
            reader
                .read_changes(&opening)
                .map_err(|refused| TestCaseError::fail(refused.to_string()))?;
            prop_assert!(
                reader.live().ids().windows(2).all(|pair| pair[0] < pair[1]),
                "the seeded set is a set"
            );

            let at = damage_at % next.len();
            next[at] = damage_to;
            if reader.read_changes(&next).is_ok() {
                let ids = reader.live().ids();
                prop_assert!(
                    ids.windows(2).all(|pair| pair[0] < pair[1]),
                    "a damaged stream was accepted and left the live set unsorted: {ids:?}"
                );
            }
        }
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
        let back = round_trip(&[records]).sets;
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

    /// **A record that departs an id and arrives the same id is damage.**
    ///
    /// No writer produces it — `derive_changes` puts an id in at most one of the two lists — so a
    /// stream that does is internally inconsistent. Honouring it would report one departure and
    /// one arrival for a record where nothing changed, and those two counts are what Milestone
    /// E4's residual arithmetic is specified against.
    #[test]
    fn an_id_that_departs_and_arrives_in_one_record_is_damage() {
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader.read_changes(&[0, 2, 7, 1]).expect("7 and 9 arrive");
        assert_eq!(reader.live().ids(), [7, 9]);

        // one departure at position 0 — id 7 — then one arrival of id 7.
        let refused = reader
            .read_changes(&[1, 0, 1, 7])
            .expect_err("departing and arriving one id in one record means nothing");
        assert!(
            matches!(refused, RecordDecodeError::Malformed { .. }),
            "got {refused}"
        );
    }

    /// **The edges of the range**: nothing at all, the largest identifier there is, and a live set
    /// far larger than any fixture above.
    #[test]
    fn the_edges_of_the_range_read_back() {
        let mut nothing = LiveSetReader::new();
        nothing.start_block();
        let refused = nothing
            .read_changes(&[])
            .expect_err("an empty slice is a stream that stopped before it started");
        assert!(
            matches!(refused, RecordDecodeError::Truncated { .. }),
            "got {refused}"
        );

        let a_thousand: Vec<ChainId> = (0..1_000).collect();
        let records = vec![
            a_thousand.clone(),
            // the last read of the thousand leaves, and the largest identifier there is arrives
            a_thousand[..999].to_vec(),
            {
                let mut with_the_largest = a_thousand[..999].to_vec();
                with_the_largest.push(u64::MAX);
                with_the_largest
            },
        ];
        let back = round_trip(std::slice::from_ref(&records)).sets;
        assert_eq!(back, records);
    }
}
