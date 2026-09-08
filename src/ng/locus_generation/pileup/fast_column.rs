//! The **ordinary column** — one covered reference base where every read that covers it
//! simply matches, and there is nothing to reconcile.
//!
//! # What "ordinary" means, and why it gets its own path
//!
//! The walk runs one fully general machine at every covered base: each read's CIGAR is
//! re-walked into a list of [`ReadEvent`](super::decompose::ReadEvent)s, the events open or
//! widen a record, the record's reference bytes are re-derived per read into an allele
//! string, the positions the read witnessed are accumulated as runs and canonicalised, the
//! read's state goes into a hash table keyed by its id, and at close the table is collected,
//! sorted, and re-read to produce the observations. Every one of those steps exists for a
//! case that is rare: an indel, a record widened by one, a hole in a witness, a mate-overlap
//! replay, a read folded twice.
//!
//! **Measured, on the two review fixtures**, the column where none of that applies is
//! 7,898 in 10,000 (tomato ~130×) and 7,789 in 10,000 (HG002 chr1 30×). For it the answer
//! is a handful of scalars per read: a base, a quality, a strand bit, two mapping-quality
//! moments, a `q_sum` term, and a witness that is by construction complete.
//!
//! # The predicate is decided **before** any per-read work
//!
//! Three tests, all cheap, and each one is a reason the general path would have had
//! something to do:
//!
//! 1. **No record is already open over this base.** Then the record this column opens is
//!    one base wide, anchored at the walker, and can never be found again — the next
//!    position's events start at `pos + 1` and half-open intervals that touch do not
//!    overlap. So no widen, no re-fold, and no read folded into it twice.
//! 2. **This base is not a region's first base beside repeat ground.** There, and nowhere
//!    else in the walk, the general path builds its record from an event window that starts
//!    one base earlier, so that an insertion anchored on the repeat's last base can be
//!    claimed by this region ([`window_of_read`](super::open_record)). See the test's own
//!    note for what that window returns that this lane cannot see.
//! 3. **No active read shows anything but one letter at this base.**
//!    Until 2026-09-08 this was a question about the whole read — every active read's CIGAR
//!    free of `I` and `D` ops — which bought two things at once: no read answers with
//!    anything but a `Match` here, *and* no event reaches in from an earlier anchor. A read
//!    with an indel at its far end failed that test at every one of the hundred bases it
//!    covers, and **that alone sent 19.1% of the 10.6 M bases of a 10 Mb walk down the
//!    general path**. It is now asked at this base
//!    ([`plain_match_at`](super::cigar_cursor::CigarCursor::plain_match_at)); the second
//!    half, which is only ever about a deletion, is what test 2 above closes.
//! 4. **No two contributors share a chain id**, i.e. no mate-overlap reconciliation fires
//!    here. **The active set answers this one before the pass, in O(1)**: a shared chain id
//!    is a mate pair, and
//!    [`may_have_mate_overlap_at`](super::active_read_set::ActiveReads::may_have_mate_overlap_at)
//!    says whether any pair the set holds still has both alignments on the reference at this
//!    column. Its `false` is exact. Only when it says a pair *could* be here — 1,664 columns
//!    in 10,000 at 130× — is the depth-sized sort over the contributors' chain ids run, and
//!    a column that fails it is handed back with nothing lost but the pass itself.
//!
//! A fifth, the per-column depth cap, is checked against the active-set size, which bounds
//! the contributor count from above.
//!
//! **Handing back is always safe.** The pass writes only into scratch buffers this module
//! owns and sets [`ever_contributed`](super::active_read_set::ActiveRead::ever_contributed),
//! which the general path sets for exactly the same reads a moment later. Nothing else is
//! touched, so a fallback re-runs the general path against untouched state.
//!
//! # Byte-identity, and the one place the risk lives
//!
//! `q_sum` is an `f64` sum, so it is not the value that has to match but the **order the
//! terms are added in**. The general path adds one term per read per observation in
//! ascending `read_id`, established by
//! [`keyed_observations_counting`](super::open_record::OpenPileupRecord)'s sort; this path
//! sorts its own compact per-read buffer by `read_id` and accumulates in that order, which
//! is the same sequence of additions. Everything else the locus carries is integer or
//! byte-valued.
//!
//! The observations themselves come out sorted by `(bases, read_witness, read_group)` as
//! `finalise` sorts them. Here every witness is `Complete` and every `bases` is one byte, so
//! that reduces to `(base, read_group)`.

use crate::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use crate::ng::ref_seq::RefSeq;
use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};
use crate::pileup_record::ChainId;

use super::active_read_set::ActiveReads;
use super::cigar_cursor::BaseShown;
use super::errors::WalkerError;
use super::open_record::{JunctionAtRegionStart, OpenPileupRecordTable, minted_ln_read_error};

/// What one read contributes to an ordinary column: everything the locus needs from it,
/// with no reference to the read left behind.
///
/// 40 bytes, so a 130-deep column's whole buffer is ~5 kB and stays in L1 across the sort
/// and the grouping pass that follow.
#[derive(Clone, Copy, Debug)]
struct PlainContribution {
    /// Sorted on — see the module note on `q_sum`.
    read_id: u32,
    chain_id: ChainId,
    read_group: ReadGroupId,
    /// The read's base at this position. Its allele, and the whole of it.
    base: u8,
    /// `ln(P_err)` for this read at this position — the BAQ-capped base quality floored by the
    /// read's mapping-quality log-error, minted by
    /// [`minted_ln_read_error`](super::open_record::minted_ln_read_error), which is the
    /// function the general fold mints with too.
    ///
    /// **This path reaches the same number by a shorter route, not over a narrower window.**
    /// It hands the mint one base's quality directly; the general fold reaches it through a
    /// window that, in every column this path accepts, holds exactly that one event — so
    /// `min_bq_for_read` returns this base's own quality and the two mints agree. The three
    /// column tests are what make that true: no record already open here, so the record the
    /// general path builds spans one base; no indel anchored here, so no second event shares
    /// the anchor; and not a region's first base beside repeat ground, so the window is
    /// `[p, p + 1)` and a deletion ending on the base before is outside it. The
    /// `debug_assert_eq!` below re-checks the event list at every column in a debug walk.
    ln_q: f64,
    /// Forward strand.
    fwd: bool,
    /// The read started strictly left of this position.
    placed_left: bool,
    mapq: u8,
}

/// One accumulating observation — an allele × read group, which at a one-base record with
/// every witness complete is the whole of an observation's identity.
#[derive(Debug)]
struct PlainObservation {
    base: u8,
    read_group: ReadGroupId,
    num_obs: u32,
    q_sum: f64,
    fwd: u32,
    placed_left: u32,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    /// **The ids of every read folded here, whether it departed from the reference or
    /// agreed with it** — the owner's ruling of 2026-08-17, whose reasons are on the field
    /// this one fills
    /// ([`SequenceObservation::chain_ids`](crate::ng::locus_generation::SequenceObservation::chain_ids)).
    /// This lane used to fill it only for a read whose base differed, matching the general
    /// path's old rule.
    chain_ids: Vec<ChainId>,
}

/// The buffers the fast lane reuses column to column, so an ordinary column allocates only
/// what it emits.
#[derive(Debug, Default)]
pub(super) struct FastColumnScratch {
    reads: Vec<PlainContribution>,
    /// `(chain_id, ())` sorted to find a mate overlap — the same shape
    /// `resolve_mate_overlap_at_pos` uses, and the same cost.
    chains: Vec<ChainId>,
    observations: Vec<PlainObservation>,
    ref_base: Vec<u8>,
    /// Which observation each contributing read landed in, indexed by that read's `read_id`
    /// less the smallest `read_id` contributing to this column — see
    /// [`OBSERVATION_OF_NO_CONTRIBUTOR`] and the chain-ordered pass that reads it.
    observation_of_read: Vec<u16>,
}

/// The entry [`FastColumnScratch::observation_of_read`] holds for a `read_id` that
/// contributed nothing to this column — a read masked to `N` here, or one inside a deleted
/// run, or simply an id no live read carries.
///
/// It doubles as the bound on how many observations a slot can name: the table is used only
/// where the live set is smaller than this, so no real observation index can reach it.
const OBSERVATION_OF_NO_CONTRIBUTOR: u16 = u16::MAX;

/// How far apart the smallest and largest contributing `read_id` may be before the column
/// gives up on the dense table above and sorts its chain ids the way the general path does.
///
/// **It is a guard and not a tuning knob.** Reads are admitted in ascending `read_id` and
/// leave in very nearly the same order, so the live set's ids are all but contiguous: on 10
/// Mb of tomato chromosome 1 at 103.5× the span never reached twice the column's depth. What
/// this bounds is the pathological set — one read held open across a long stretch while
/// thousands are admitted past it — where a dense table indexed by id offset would be a
/// megabyte of memset per column.
const MAX_READ_ID_SPAN: u32 = 4_096;

/// What [`try_ordinary_column`] decided.
#[derive(Debug)]
pub(super) enum FastColumn {
    /// The column was ordinary and this is its locus, ready to be emitted at the position
    /// the general path would have emitted it — see `WalkerState::close_aged_records_into`.
    Emitted(SampleLocusObservations),
    /// Not ordinary. Nothing was changed; run the general path.
    Fallback,
}

/// Try to answer this column with scalars. See the module note for the predicate and for
/// why handing back is free.
///
/// **Eight arguments, and each is a distinct fact the walker already holds** — where it is,
/// what is live, where the reference is, where to put the scratch, and the two limits the
/// column is judged against. Grouping them into a context struct would name the caller's own
/// fields a second time without hiding anything; the alternative that would genuinely shorten
/// the list is moving the predicate out to `genome_walk.rs`, which is the one thing this
/// module exists to keep in one place.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_ordinary_column(
    walker_pos: u32,
    chrom_id: u32,
    active_reads: &ActiveReads,
    open_records: &OpenPileupRecordTable,
    reference: &dyn RefSeq,
    scratch: &mut FastColumnScratch,
    max_snp_column_depth: usize,
    may_have_mate_overlap: bool,
    junction: Option<JunctionAtRegionStart>,
) -> Result<FastColumn, WalkerError> {
    // The contributor count is at most the active-read count, so this bounds the column
    // cap from above without knowing which reads contribute.
    if active_reads.len() > max_snp_column_depth || active_reads.is_empty() {
        return Ok(FastColumn::Fallback);
    }
    // A record already covering this base is a widen, a re-fold, or a footprint wider than
    // one — all three of which the scalars below cannot express.
    if open_records
        .find_overlapping(walker_pos, walker_pos.saturating_add(1))
        .is_some()
    {
        return Ok(FastColumn::Fallback);
    }
    // **The one base of the walk whose record does not ask about its own base.** Where the
    // region begins immediately after repeat ground, the general path builds every record on
    // the region's first base from a window that starts one base *earlier*
    // ([`window_of_read`](super::open_record)), so that an insertion anchored on the repeat's
    // last base can be claimed by this region. That widened window returns two things this
    // lane cannot see, and both change the emitted bytes:
    //
    // - **a deletion that ends on the repeat's last base.** `window_of_read` keeps every
    //   deletion whatever its anchor, so a read whose deleted run stops one base short of
    //   here still contributes its deletion's quality proxy to `min_bq_for_read`, and the
    //   read's `ln ε` at this base is that proxy rather than the base's own quality.
    // - **the claimed junction insertion itself**, which the fold spells before this base
    //   and this lane has no allele for.
    //
    // Refusing the column costs one base per region and removes both.
    if junction.is_some_and(|junction| junction.first_base == walker_pos) {
        return Ok(FastColumn::Fallback);
    }

    scratch.reads.clear();
    for active in active_reads.iter() {
        // **Asked at this base, not of the whole read.** Until 2026-09-08 a read carrying
        // any `I` or `D` op anywhere in its alignment was refused at every one of the
        // hundred bases it covers, and that alone sent 19.1% of the 10.6 M bases of a 10 Mb
        // walk down the general path. `plain_match_at` asks the question the lane actually
        // needs — does this read show one letter *here* — and answers `NotOneLetter` for the
        // single base an indel is anchored at.
        //
        // # Why nothing has to be asked of the read's history any more
        //
        // The whole-read test also bought a second thing: that no event anchored at an
        // earlier base is still reaching in here. A deletion's footprint is the only one
        // that does, and every way it can reach a record on this base is now closed by a
        // column test above rather than by a read test here:
        //
        // - **at the deletion's own anchor**, `plain_match_at` returns `NotOneLetter` — the
        //   anchor is the last reference base of the match run before the `D` op;
        // - **inside the deleted run**, the read shows no letter at all, and the record the
        //   deletion opened covers this base, so `find_overlapping` above has already handed
        //   the column back;
        // - **one base past the deleted run**, the ordinary window `[p, p + 1)` excludes the
        //   deletion — its footprint ends exactly at `p` — *except* on a region's first base
        //   beside repeat ground, where the window starts a base early. That column is
        //   refused above.
        //
        // # What refusing costs
        //
        // Nothing but this pass. Handing back re-runs the general path against untouched
        // state, which is what makes a conservative guard here free to be conservative.
        let shown = active.cursor.plain_match_at(walker_pos, &active.read);
        if shown == BaseShown::NotOneLetter {
            return Ok(FastColumn::Fallback);
        }
        // **The two paths cross-checked where they genuinely duplicate a computation**, and
        // only there: everything after this loop is a sum over what it produced. Armed in
        // every debug walk in the suite — the parity census alone runs it over ~257,000 loci
        // — which is what makes a second path through the walk's hottest code defensible.
        //
        // **It checks the whole event list and not only its head, and that is what carries the
        // per-base test.** Until 2026-09-08 the lane refused any read whose CIGAR held an indel
        // *anywhere*, which made "one Match here" true by construction; now it refuses only a
        // read with an indel anchored *here*, so the property that has to hold — this read
        // raises exactly one event at this base and it is the Match this lane is about to
        // record — is a claim rather than a consequence. This is where it is claimed.
        debug_assert_eq!(
            match shown {
                BaseShown::Letter(base, bq) => Some((base, bq)),
                BaseShown::Nothing | BaseShown::NotOneLetter => None,
            },
            match active.cursor.events_at(walker_pos, &active.read).as_slice() {
                [super::decompose::ReadEvent::Match { base, bq_baq, .. }] => {
                    Some((*base, *bq_baq))
                }
                [] => None,
                _ => Some((b'?', u8::MAX)),
            },
            "plain_match_at and events_at disagree at {walker_pos} for read {}",
            active.read_id,
        );
        let BaseShown::Letter(base, bq) = shown else {
            continue;
        };
        // Set here rather than in the loop below, and for the same reason the general path
        // sets it before the mate collapse and the depth cap: it exists to find the read
        // that reaches *no* contributor list anywhere. A read counted here and then handed
        // to the general path is counted there too, idempotently.
        active.ever_contributed.set(true);
        let mapq = active.read.mapq;
        scratch.reads.push(PlainContribution {
            read_id: active.read_id,
            chain_id: active.chain_id,
            read_group: active.read.read_group,
            base,
            ln_q: minted_ln_read_error(bq, active.read.mq_log_err),
            fwd: !active.read.is_reverse_strand,
            placed_left: active.read.alignment_start < walker_pos,
            mapq,
        });
    }
    if scratch.reads.is_empty() {
        return Ok(FastColumn::Fallback);
    }

    // Mate overlap: two contributors from one pair. The general path reconciles their
    // qualities against each other, which is a rewrite this path has no term for.
    //
    // **The sort only runs where a pair could be.** `may_have_mate_overlap` is the active
    // set's O(1) answer to exactly the question this sort asks; a `false` is exact, so the
    // set of columns this function accepts is the same either way — `fast_columns` is
    // unchanged by the shortcut. A `true` is an over-approximation (the pair may be silent
    // here, or one mate may be `N`-masked), so it still has to be settled by the sort.
    if may_have_mate_overlap {
        scratch.chains.clear();
        scratch
            .chains
            .extend(scratch.reads.iter().map(|r| r.chain_id));
        scratch.chains.sort_unstable();
        if scratch.chains.windows(2).any(|w| w[0] == w[1]) {
            return Ok(FastColumn::Fallback);
        }
    } else {
        // The pin, in debug builds only: the sort the shortcut skipped, run in full. A
        // wrong shortcut here would keep both mates' quality instead of zeroing the lower
        // one, change the emitted bytes, and show up nowhere else.
        debug_assert!(
            {
                let mut chains: Vec<ChainId> = scratch.reads.iter().map(|r| r.chain_id).collect();
                chains.sort_unstable();
                chains.windows(2).all(|w| w[0] != w[1])
            },
            "the ordinary-column path skipped its chain-id sort at {walker_pos} on a column \
             where two contributors share a chain id — the mate reconciliation was lost",
        );
    }

    // **The determinism guarantee**, and the only ordering in this function that is not
    // free to change: see the module note.
    //
    // **Kept even though the ordered active set makes it redundant.** With `ActiveReads` a
    // queue in ascending `read_id`, this buffer arrives sorted and the sort is a scan.
    // Deleting it was measured — −0.44 % at 130×, −0.18 % at 30×, −0.08 % at 300× — and
    // declined: it would make this module's `q_sum` summation order, and so the emitted
    // bytes, depend on a container choice in `active_read_set.rs`. Half a percent is not
    // worth that coupling. See `composed_full.md` §7.
    scratch.reads.sort_unstable_by_key(|r| r.read_id);

    scratch.ref_base.clear();
    reference
        .fetch_into(
            ContigId(chrom_id),
            u64::from(walker_pos),
            1,
            &mut scratch.ref_base,
        )
        .map_err(|source| WalkerError::Fasta {
            chrom_id,
            start: walker_pos,
            start_plus_len: walker_pos.saturating_add(1),
            source,
        })?;
    // PANIC-FREE: `fetch_into` was asked for one base and returned `Ok`.
    let ref_base = *scratch
        .ref_base
        .first()
        .expect("one base was fetched and the fetch succeeded");

    scratch.observations.clear();
    // **The dense table the chain-ordered pass below reads.** One slot per `read_id` between
    // the smallest and largest contributing to this column, which is the column's depth plus
    // whatever ids expired inside it.
    //
    // PANIC-FREE: `scratch.reads` was found non-empty above and is ascending by `read_id`.
    let first_read_id = scratch.reads.first().expect("checked non-empty").read_id;
    let last_read_id = scratch.reads.last().expect("checked non-empty").read_id;
    let read_id_span = last_read_id - first_read_id + 1;
    // The second condition is what makes the `u16` slot below exact rather than lossy: a
    // column has at most one observation per contributor, and at most one contributor per
    // live read, so a set smaller than the sentinel cannot name a slot the sentinel collides
    // with. The walk's own admission ceiling is far below it — 32,768 — so this never binds.
    let dense = read_id_span <= MAX_READ_ID_SPAN
        && active_reads.len() < OBSERVATION_OF_NO_CONTRIBUTOR as usize;
    if dense {
        scratch.observation_of_read.clear();
        scratch
            .observation_of_read
            .resize(read_id_span as usize, OBSERVATION_OF_NO_CONTRIBUTOR);
    }
    // **Measurement scaffolding, hoisted out of the loop** — this lane emits observations
    // without ever building an open record, so a census hooked only into `finalise` would
    // miss every read that came through here. See
    // [`minted_error_census`](super::minted_error_census); off unless armed.
    let census_armed = super::minted_error_census::enabled();
    for read in &scratch.reads {
        // The read's own `ln ε`, before the loop below pools it into an observation's
        // `q_sum` and the read is gone. Every witness in this lane is complete, and the
        // lane refuses any column an open record overlaps — so these are the same reads the
        // pre-pass's calibration accumulator sums over.
        if census_armed {
            super::minted_error_census::record_read(
                u64::from(walker_pos),
                read.read_group,
                read.ln_q,
            );
        }
        let at = match scratch
            .observations
            .iter()
            .position(|o| o.base == read.base && o.read_group == read.read_group)
        {
            Some(at) => at,
            None => {
                scratch.observations.push(PlainObservation {
                    base: read.base,
                    read_group: read.read_group,
                    num_obs: 0,
                    q_sum: 0.0,
                    fwd: 0,
                    placed_left: 0,
                    mapq_sum: 0,
                    mapq_sum_sq: 0,
                    chain_ids: Vec::new(),
                });
                scratch.observations.len() - 1
            }
        };
        let observation = &mut scratch.observations[at];
        let mapq = u32::from(read.mapq);
        observation.num_obs += 1;
        observation.q_sum += read.ln_q;
        observation.fwd += u32::from(read.fwd);
        observation.placed_left += u32::from(read.placed_left);
        observation.mapq_sum += mapq;
        observation.mapq_sum_sq += u64::from(mapq) * u64::from(mapq);
        if dense {
            // PANIC-FREE: `read_id` is between the first and last of an ascending buffer, and
            // the table was resized to that span.
            // PANIC-FREE: `at` indexes a list with one entry per contributor at most, and
            // `dense` is false unless the live set is smaller than the sentinel.
            scratch.observation_of_read[(read.read_id - first_read_id) as usize] =
                u16::try_from(at).expect("fewer observations than the sentinel");
        } else {
            observation.chain_ids.push(read.chain_id);
        }
    }

    // **The chain ids, in ascending order, without ordering them.** They have to come out
    // ascending — the general path reaches that by sorting each observation's list once the
    // reads are folded, and this lane did the same until the active set began keeping the
    // order (`ActiveReads::chain_order`). Walking that order and pushing as we go leaves each
    // observation's list ascending by construction, because a subsequence of an ascending
    // sequence is ascending.
    //
    // **No `dedup`, and it is not an omission.** Two entries of one observation could only
    // repeat a chain id if two contributors shared one, which is a mate overlap — and this
    // lane refuses those columns outright, exactly so that it never has to reconcile a pair.
    // The `debug_assert` states it.
    //
    // The fallback arm is the old shape and is reached only when the live read ids are spread
    // wider than `MAX_READ_ID_SPAN`.
    if dense {
        for &(chain_id, read_id) in active_reads.chain_order() {
            let Some(offset) = read_id.checked_sub(first_read_id) else {
                continue;
            };
            let Some(&at) = scratch.observation_of_read.get(offset as usize) else {
                continue;
            };
            if at != OBSERVATION_OF_NO_CONTRIBUTOR {
                scratch.observations[at as usize].chain_ids.push(chain_id);
            }
        }
        debug_assert!(
            scratch
                .observations
                .iter()
                .all(|o| o.chain_ids.windows(2).all(|w| w[0] < w[1])),
            "the chain-ordered pass left an observation's chain ids unordered or repeated at \
             {walker_pos}",
        );
    } else {
        for o in &mut scratch.observations {
            o.chain_ids.sort_unstable();
            o.chain_ids.dedup();
        }
    }

    // `finalise` sorts on `(bases, read_witness, read_group)`; one byte of bases and one
    // witness class collapse that to this.
    scratch
        .observations
        .sort_unstable_by_key(|o| (o.base, o.read_group.0));

    let observations = scratch
        .observations
        .iter_mut()
        .map(|o| {
            SequenceObservation {
                bases: vec![o.base].into_boxed_slice(),
                read_witness: ReadWitness::Complete,
                read_group: o.read_group,
                num_obs: o.num_obs,
                num_fwd: o.fwd,
                // Rounded once, here, where the sum is finished — not per read. See
                // `SummedLogError::from_nats`.
                q_sum: SummedLogError::from_nats(o.q_sum),
                mapq_sum: o.mapq_sum,
                mapq_sum_sq: o.mapq_sum_sq,
                placed_left: o.placed_left,
                chain_ids: std::mem::take(&mut o.chain_ids),
            }
        })
        .collect();

    super::column_census::FAST_COLUMNS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    Ok(FastColumn::Emitted(SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(chrom_id),
            start: Position(u64::from(walker_pos)),
            end: Position(u64::from(walker_pos)),
        },
        reference_bases: vec![ref_base].into_boxed_slice(),
        observations,
        // A read that produced no observation did not reach this path — every entry in
        // `reads` carries a `Match` — and the depth cap is gated out above.
        reads_without_observation: 0,
        reads_discarded_by_cap: 0,
        kind: LocusKind::Generic,
    }))
}
