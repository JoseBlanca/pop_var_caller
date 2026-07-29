//! Open-record formation, merging, and widening — the
//! walker's central operation. See `ia/specs/pileup_walker.md`
//! §"Open-record formation and merging".
//!
//! An open record is the in-flight version of a `PileupRecord`:
//! same shape, but its REF span and allele list grow as the walker
//! folds events into it. When the walker confirms no future event
//! can touch the record (its footprint is fully behind the walker
//! per the closure rule), it converts the open record into a
//! finalised `PileupRecord` and pushes it through the channel.
//!
//! **No longer a verbatim copy — A0 (plan 3).** Copied from
//! `src/pileup/walker/open_record.rs`, then changed: the reference is reached
//! through ng's [`RefSeq`] rather than production's `MultiChromRefFetcher`, and
//! through [`RefSeq::fetch_into`] rather than `fetch`, so `widen` writes into a
//! buffer this table owns instead of allocating a `Vec<u8>` per call.
//! `copy_fidelity.rs` released this file in that commit.

use std::collections::BTreeMap;

use ahash::AHashMap;

use crate::ng::ref_seq::RefSeq;
use crate::ng::types::{ContigId, ReadGroupId};
use crate::pileup_record::{
    AlleleObservation, AlleleSupportStats as ProductionAlleleSupportStats, ChainId, PileupRecord,
};

use super::DEFAULT_MAX_RECORD_SPAN;
use super::active_read_set::ActiveReads;
use super::decompose::ReadEvent;
use super::errors::WalkerError;

/// Pre-allocated capacity for `OpenPileupRecord::folded_reads` —
/// sized for typical WGS coverage so the per-record fold doesn't
/// pay 4–5 grow reallocations as contributors accumulate.
/// Previously the largest remaining alloc site by bytes (228 MB
/// on `pileup_walker_multi_op/L=5000`); see L6 in
/// `ia/reviews/perf_pileup_2026-05-10.md`.
///
/// H1 in `ia/reviews/perf_pileup_2026-05-12.md` tried swapping
/// the `AHashMap` for a sorted `Vec<(u32, FoldedReadState)>` —
/// the cumulative bench regressed 1.3 % on the mean with four of
/// eight fixtures regressing 3–12 % (worst on `multi_op/5000`,
/// +11.8 %). The `Vec` doubled-on-grow past cap 32 inflated bytes
/// (131 MB AHashMap → 166 MB Vec) and the `Vec::remove` shift on
/// re-fold was more expensive than the AHashMap probe. Reverted;
/// keep the `AHashMap` shape until a smaller-than-AHashMap
/// container (e.g. an arena-pooled map or a perfect hash on
/// dense `read_id` ranges) can be evaluated.
const RECORD_FOLDED_READS_INITIAL_CAPACITY: usize = 32;

/// A stretch of reference positions, **1-based and inclusive of both ends** —
/// the convention [`PreparedRead::alignment_start`](super::PreparedRead) and
/// `alignment_end` already use, so a span built from a read's own coordinates
/// needs no adjustment.
///
/// Deliberately **not** [`GenomeRegion`](crate::ng::types::GenomeRegion), which
/// carries a `ContigId` and a `u64` this fold has no use for: an open record is
/// already pinned to one contig, and every coordinate inside the walk is `u32`.
///
/// *Invariant:* `start <= end`. There is no empty span — a read that witnessed
/// nothing inside a record does not fold into it at all, so "no positions" is
/// the absence of a `RefSpan`, never a degenerate one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct RefSpan {
    pub start: u32,
    pub end: u32,
}

/// The support one allele bucket has accumulated — **ng's own**, and production's
/// [`AlleleSupportStats`](crate::pileup_record::AlleleSupportStats) **minus
/// `placed_start`**.
///
/// `placed_left` stays because something computes on it:
/// [`vcf::qual_refine`](crate::vcf) turns it into the read-position-bias term
/// subtracted from QUAL, live through `final_qual` into the cohort VCF writer
/// and the `--min-qual` gate, so dropping it would forfeit the ability to
/// reproduce production's QUAL. `placed_start` is merged, serialised and printed
/// by the `psp-to-pileup` dump and read by **nothing that computes**, and it is
/// cheap to reverse: both counters are pure functions of the read's start against
/// the record's anchor (spec §6, arch §3).
///
/// # It is still emitted, for now, and that is deliberate
///
/// [`OpenPileupRecord::finalise`] still returns production's [`PileupRecord`],
/// which has the field — so until Milestone B changes the output type, the value
/// is **reconstructed at that boundary** from a per-read flag on
/// [`FoldedReadState`], exactly and by definition. That is what keeps the stage-1
/// differential comparing every field of every record through Milestone A, which
/// is the oracle A2 and A3 are checked against; dropping the number here instead
/// would have blinded that oracle, and broken one inherited test inside the still
/// verbatim `tests.rs`, four steps before either needed to move.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct AlleleSupportStats {
    /// Number of supporting reads for this allele in this record.
    pub num_obs: u32,
    /// `Σ max(ln(P_err_BQ_BAQ), ln(P_err_MQ))` over supporting reads.
    pub q_sum: f64,
    /// Reads on the forward strand among `num_obs`.
    pub fwd: u32,
    /// Reads whose mapped 5′ end is strictly to the left of this record's anchor
    /// position (freebayes' `placedLeft`).
    pub placed_left: u32,
    /// Σ mapping quality over supporting reads.
    pub mapq_sum: u32,
    /// Σ mapq² over supporting reads.
    pub mapq_sum_sq: u64,
}

/// One in-flight allele bucket inside an `OpenPileupRecord`.
///
/// The bucket's `chain_ids` are *not* tracked here — they are
/// derived from `OpenPileupRecord::folded_reads` at `finalise()`
/// time. Storing them on the bucket and the fold state in parallel
/// invites the two to drift on re-fold (a read leaving a bucket
/// needs both its scalar contribution and its chain id removed —
/// historically only the scalars were subtracted, leaving stale
/// chain ids behind on `num_obs == 0` buckets). Projection at
/// emit time is impossible to drift by construction.
#[derive(Debug, Clone)]
pub(super) struct OpenAllele {
    pub seq: Vec<u8>,
    pub support: AlleleSupportStats,
}

impl OpenAllele {
    fn new(seq: Vec<u8>) -> Self {
        Self {
            seq,
            support: AlleleSupportStats::default(),
        }
    }
}

/// One in-flight per-position record. Its REF span is
/// `alleles[0].seq.len()`; no separate `ref_span` field is stored.
/// `alleles[0]` is always REF — its `seq` field is the canonical
/// reference sequence under the record's footprint and is the only
/// place those bytes are stored.
#[derive(Debug, Clone)]
pub(super) struct OpenPileupRecord {
    pub chrom_id: u32,
    /// 1-based anchor position.
    pub pos: u32,
    pub alleles: Vec<OpenAllele>,
    /// Per-read fold state — the contribution this record currently
    /// holds for each read that has folded into it. Used to enforce
    /// "fold each (record, read) pair exactly once" across walker
    /// steps: at re-fold time the previous contribution is
    /// subtracted from its old bucket before the new contribution is
    /// added to the new bucket. Without this state, a read with
    /// Match events at every position inside an open record's
    /// footprint would be re-folded once per walker step inside
    /// that footprint, multiplying every five-scalar value by
    /// `ref_span` (B1 in `ia/reviews/pileup_2026-05-06.md`).
    ///
    folded_reads: AHashMap<u32, FoldedReadState>,
}

/// What a single read currently contributes to one bucket of one
/// open record. Carries enough state to subtract the contribution
/// cleanly when the read re-folds (e.g. on widening that grows the
/// haplotype seq under the record). `chain_id` is also stored here
/// so the per-bucket chain id set can be projected at `finalise()`
/// time straight from the current fold state — see [`OpenAllele`].
#[derive(Debug, Clone, Copy)]
struct FoldedReadState {
    allele_index: usize,
    contribution: AlleleSupportStats,
    chain_id: ChainId,
    /// The read group this read was prepared in — part of the observation's
    /// identity from B1 on, so a per-library model gets the allele × group cross
    /// *with its quality moments* rather than a count beside a merged
    /// observation (spec §6). Copied off the contributor like the fields around
    /// it; ng's [`PreparedRead`](super::PreparedRead) carries it, so nothing here
    /// reconstructs anything.
    ///
    /// Unread until B1 puts it in the bucket key. The `expect` is a compiler
    /// backstop that forces its own removal the moment that lands — an `allow`
    /// would simply go stale.
    #[expect(
        dead_code,
        reason = "B1 makes the read group part of the bucket key; this expect must \
                  be deleted when it lands"
    )]
    read_group: ReadGroupId,
    /// The positions this read actually witnessed inside the record — the union
    /// of its event footprints, **not** its alignment span. The span is blind to
    /// `N`, adaptor-masked, ref-skipped and dropped-indel positions, all of which
    /// production fills from the reference (spec §6).
    ///
    /// Held in **absolute reference coordinates**, not relative to the footprint:
    /// the bucket a read folds into is chosen from its bases, which are fixed at
    /// fold time, while its coverage is relative to a footprint that is not.
    /// Coverage is therefore resolved once, at `finalise()`, against the record's
    /// **final** footprint — by which time the read may be long gone, since
    /// `expire_passed` touches no open record (spec §4).
    ///
    /// **A1 fills this with the record's whole footprint**, which is precisely
    /// what production's fill assumes every folded read witnessed — so this step
    /// changes no answer. A2 replaces it with the extent the events actually
    /// cover, and that is where the answers move.
    ///
    /// Unread until A4 resolves coverage from it. Same backstop, same reason as
    /// `read_group` above.
    #[expect(
        dead_code,
        reason = "A4 resolves ReadCoverage from this at finalise; this expect must \
                  be deleted when it lands"
    )]
    witnessed: RefSpan,
    /// Whether this read's mapped 5′ end *is* the record's anchor position
    /// (freebayes' `placedStart`).
    ///
    /// **Transitional, and it disappears with the `PileupRecord` boundary.** ng's
    /// [`AlleleSupportStats`] does not carry `placed_start`, but `finalise` still
    /// has to emit one; a per-read flag reconstructs the per-bucket count exactly,
    /// because that count *is* "how many currently-folded reads in this bucket
    /// started on the anchor". Milestone B deletes both this field and the
    /// conversion that reads it.
    placed_start: bool,
}

impl OpenPileupRecord {
    /// Open a fresh record at `pos` with the given initial REF
    /// sequence. The REF allele bucket is created up front (with
    /// zero observations) so the `alleles[0] == REF` invariant
    /// holds from the very start. The REF bytes are moved into the
    /// REF bucket — there is no separate copy on the record.
    fn new(chrom_id: u32, pos: u32, ref_seq: Vec<u8>) -> Self {
        Self {
            chrom_id,
            pos,
            alleles: vec![OpenAllele::new(ref_seq)],
            folded_reads: AHashMap::with_capacity(RECORD_FOLDED_READS_INITIAL_CAPACITY),
        }
    }

    pub fn ref_span(&self) -> u32 {
        self.alleles[0].seq.len() as u32
    }

    /// Footprint end (exclusive), in 1-based coordinates: the
    /// position one past the last reference base this record
    /// covers. `saturating_add` defends against `pos` near
    /// `u32::MAX` on multi-Gbp chromosomes — without it the
    /// returned end wraps and `drain_aged` misreads it. Mi8 in
    /// `ia/reviews/pileup_2026-05-11.md`.
    pub fn footprint_end_exclusive(&self) -> u32 {
        self.pos.saturating_add(self.ref_span())
    }

    /// Convert into a finalised `PileupRecord`. Per-bucket
    /// `chain_ids` are projected from `folded_reads` here rather
    /// than carried on each `OpenAllele`: a re-fold subtracts the
    /// read's scalar contribution from the old bucket and adds it
    /// to the new one, but the old bucket would retain the read's
    /// chain id unless the walker explicitly tracked it. Deriving
    /// the chain id set from `folded_reads` at finalise removes
    /// that drift surface entirely.
    ///
    /// REF (allele 0) chain ids are intentionally **not** emitted —
    /// see the skip in the projection loop below.
    ///
    /// **`placed_start` is projected here too, for the same reason and by the
    /// same route.** ng's [`AlleleSupportStats`] does not carry it (spec §6), but
    /// this function still returns production's `PileupRecord`, which has the
    /// field. Per bucket it is "how many currently-folded reads in that bucket
    /// started on the record's anchor", which is exactly what a per-read flag in
    /// `folded_reads` counts — so the reconstruction is exact rather than
    /// approximate, and it disappears with this whole conversion at Milestone B.
    pub fn finalise(self) -> PileupRecord {
        let mut per_bucket: Vec<Vec<ChainId>> =
            (0..self.alleles.len()).map(|_| Vec::new()).collect();
        let mut placed_start_per_bucket: Vec<u32> = vec![0; self.alleles.len()];
        for state in self.folded_reads.values() {
            // Counted before the REF skip below: `placed_start` is a property of
            // every bucket, REF included, where chain ids are not.
            if state.placed_start {
                placed_start_per_bucket[state.allele_index] += 1;
            }
            // Drop REF (allele 0) chain ids. The only downstream reader
            // of `chain_ids` is Stage 5's compound-haplotype check
            // (`per_group_merger::build_chain_proposals`), which links
            // *non-reference* alleles only and explicitly skips allele 0,
            // so REF chain ids are provably never read. They are the vast
            // majority of all chain ids (~96.6% on real cohorts), so not
            // materialising them keeps them out of the `.psp` entirely.
            // See doc/devel/specs/phase_chain.md §8.
            if state.allele_index == 0 {
                continue;
            }
            per_bucket[state.allele_index].push(state.chain_id);
        }
        for chain_ids in &mut per_bucket {
            chain_ids.sort_unstable();
            chain_ids.dedup();
        }
        let alleles = self
            .alleles
            .into_iter()
            .zip(per_bucket)
            .zip(placed_start_per_bucket)
            .map(|((a, chain_ids), placed_start)| AlleleObservation {
                seq: a.seq,
                support: to_production_support(a.support, placed_start),
                chain_ids,
            })
            .collect();
        PileupRecord {
            chrom_id: self.chrom_id,
            pos: self.pos,
            alleles,
            // The walker cannot compute the centred window (it needs
            // look-ahead); the pileup→psp seam fills these from the
            // sliding-window accumulator before writing.
            windowed_gc: f32::NAN,
            windowed_coverage: f32::NAN,
        }
    }
}

/// Widen ng's support stats back to production's, re-attaching the
/// `placed_start` this fold no longer tracks.
///
/// **Transitional — it dies with the `PileupRecord` return at Milestone B.**
/// Written as an exhaustive destructure of ng's type rather than a
/// field-by-field read, so a field added to ng's stats stops this compiling
/// instead of being silently dropped on the way out; that is the direction that
/// can lose information, since production's type is the frozen one.
fn to_production_support(
    stats: AlleleSupportStats,
    placed_start: u32,
) -> ProductionAlleleSupportStats {
    let AlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        mapq_sum,
        mapq_sum_sq,
    } = stats;
    ProductionAlleleSupportStats {
        num_obs,
        q_sum,
        fwd,
        placed_left,
        placed_start,
        mapq_sum,
        mapq_sum_sq,
    }
}

/// The set of currently-open records, keyed by anchor position.
/// Range queries (find records overlapping a given event span)
/// use the BTreeMap's ordered structure.
#[derive(Debug)]
pub(super) struct OpenPileupRecordTable {
    /// 1-based anchor position → record.
    records: BTreeMap<u32, OpenPileupRecord>,
    /// Reusable scratch buffer for the per-(record, contributor)
    /// haplotype string built by `apply_events_to_ref_into`. Hoisted
    /// here so the inner-fold's allele-equality check can run
    /// against a borrowed `&[u8]` and only `clone()` when adding a
    /// genuinely new allele. See L10/L11 in
    /// `ia/reviews/perf_pileup_2026-05-10.md`.
    allele_seq_buf: Vec<u8>,
    /// Reusable scratch buffer for `drain_aged_into`'s "keys to
    /// remove" list. Sized to ~1 entry per call at steady state;
    /// cleared between calls so a per-call `Vec` allocation is not
    /// paid each walker step. H2 in
    /// `ia/reviews/perf_pileup_2026-05-12.md`.
    closing_keys_buf: Vec<u32>,
    /// Reusable destination for [`RefSeq::fetch_into`] in `widen` —
    /// the buffer that came free with A0's move off
    /// `MultiChromRefFetcher::fetch`, which allocated a fresh
    /// `Vec<u8>` on every widen. `open_new` does **not** use it: the
    /// bytes it fetches are moved into the record and kept, so a
    /// scratch buffer there would only add a copy.
    widen_bases_buf: Vec<u8>,
    /// Per-instance cap on per-record reference span, mirrors
    /// `WalkerConfig::max_record_span`. M11 in
    /// `ia/reviews/pileup_2026-05-11.md`.
    max_record_span: u32,
}

impl Default for OpenPileupRecordTable {
    fn default() -> Self {
        Self::new()
    }
}

impl OpenPileupRecordTable {
    /// Construct with the default `max_record_span`. Used by
    /// tests; production code calls
    /// [`OpenPileupRecordTable::with_cap`] from `walker::run`.
    pub fn new() -> Self {
        Self::with_cap(DEFAULT_MAX_RECORD_SPAN)
    }

    /// Construct with an explicit per-record reference span cap.
    pub fn with_cap(max_record_span: u32) -> Self {
        Self {
            records: BTreeMap::new(),
            allele_seq_buf: Vec::new(),
            closing_keys_buf: Vec::new(),
            widen_bases_buf: Vec::new(),
            max_record_span,
        }
    }

    /// Reset chromosome-scoped state in place. Used at chromosome
    /// boundaries instead of replacing `self` with a fresh
    /// `OpenPileupRecordTable::new()`, which would discard the
    /// perf-hoisted `allele_seq_buf`. The caller is expected to
    /// have drained `records` already (via `drain_all`); the
    /// debug-assert pins that contract. Mi11 in
    /// `ia/reviews/pileup_2026-05-11.md`.
    pub fn reset(&mut self) {
        debug_assert!(
            self.records.is_empty(),
            "OpenPileupRecordTable::reset called with {} records still open",
            self.records.len(),
        );
        self.records.clear();
        self.allele_seq_buf.clear();
        self.closing_keys_buf.clear();
        self.widen_bases_buf.clear();
    }

    /// Drain every record whose footprint is fully behind the
    /// walker (`pos + ref_span ≤ walker_pos`), in coordinate
    /// order. Used by the walker's `close_aged_records` step.
    ///
    /// **No early break on first not-yet-aged key.** A record at
    /// key `q < walker_pos` whose footprint extends past
    /// `walker_pos` (a wide deletion still in flight) does *not*
    /// imply later keys are also not-aged: a narrower record
    /// opened after it may have closed already. See finding Mi6
    /// in `ia/reviews/pileup_2026-05-11.md`. Iteration is bounded
    /// to keys strictly less than `walker_pos` — a record at
    /// `pos ≥ walker_pos` has `footprint_end_exclusive ≥ pos + 1
    /// > walker_pos`, so it cannot be aged this step.
    pub fn drain_aged_into(&mut self, walker_pos: u32, out: &mut Vec<OpenPileupRecord>) {
        out.clear();
        self.closing_keys_buf.clear();
        for (&pos, rec) in self.records.range(..walker_pos) {
            if rec.footprint_end_exclusive() <= walker_pos {
                self.closing_keys_buf.push(pos);
            }
        }
        out.reserve(self.closing_keys_buf.len());
        for pos in self.closing_keys_buf.drain(..) {
            if let Some(rec) = self.records.remove(&pos) {
                out.push(rec);
            }
        }
    }

    /// Drain everything unconditionally (chromosome boundary or
    /// end-of-input). Records come out in coordinate order
    /// because `BTreeMap::into_values` iterates by key order.
    pub fn drain_all(&mut self) -> Vec<OpenPileupRecord> {
        std::mem::take(&mut self.records).into_values().collect()
    }

    /// Find the open record (if any) whose footprint overlaps the
    /// half-open interval `[event_start, event_end)`. "Overlap"
    /// here is non-empty interval intersection — touching
    /// intervals are not overlapping. Returns the anchor position
    /// of the matched record.
    ///
    /// Precondition: `event_start < event_end`. Empty events
    /// (`event_start == event_end`) are not produced by any caller —
    /// every `ReadEvent::footprint_span()` returns ≥ 1.
    pub fn find_overlapping(&self, event_start: u32, event_end: u32) -> Option<u32> {
        debug_assert!(
            event_start < event_end,
            "find_overlapping called with empty event [{event_start}, {event_end})",
        );
        // Candidates are records whose anchor `Q ≤ event_start`
        // (any record opened to the right of the event's start
        // would have its footprint start ≥ event_end > event_start
        // — they can't overlap). The search range is bounded by
        // `max_record_span`: a record's footprint can extend at
        // most `max_record_span` past its anchor, so any record
        // anchored before `event_start - max_record_span` cannot
        // reach `event_start`. Heterogeneous spans coexist (a
        // wide deletion record may stay open while shorter records
        // open and close around it), so an early break at the
        // first record whose footprint ends at or before
        // `event_start` would miss a wide earlier record sitting
        // behind a narrow intermediate one.
        //
        // The range bound `lo..=event_start` already guarantees
        // `q ≤ event_start`, and the precondition gives
        // `event_start < event_end`, so `q < event_end` is
        // implied. Mi8 in `ia/reviews/pileup_2026-05-09.md`.
        let lo = event_start.saturating_sub(self.max_record_span);
        for (&q, rec) in self.records.range(lo..=event_start).rev() {
            if rec.footprint_end_exclusive() > event_start {
                return Some(q);
            }
        }
        None
    }

    /// Widen the record at `key` so its REF span covers up to
    /// `new_end_exclusive`, fetching the additional reference
    /// bases. Existing alleles are rewritten by appending the
    /// new reference bases (alleles whose previous coverage
    /// already ran to the end of the old span); deletion alleles
    /// keep their existing length, expressing "more bases deleted
    /// relative to the wider REF" — see
    /// `ia/specs/pileup_walker.md` §"Step 4.2".
    ///
    /// Returns `true` when the record actually widened; `false`
    /// when `new_end_exclusive` was already covered (no-op).
    /// Callers use the bool to count real widen events without
    /// conflating them with fresh `open_new` calls.
    fn widen(
        &mut self,
        key: u32,
        new_end_exclusive: u32,
        reference: &dyn RefSeq,
    ) -> Result<bool, WalkerError> {
        // Split-borrow so the fetch can write into `widen_bases_buf` while
        // `rec` holds a mutable borrow of `records`. `max_record_span` is
        // read through the same destructure rather than through `self`,
        // which the outstanding `records` borrow would otherwise block.
        let Self {
            records,
            widen_bases_buf,
            max_record_span,
            ..
        } = self;
        let max_record_span = *max_record_span;
        // PANIC-FREE: `widen` is only called from
        // `process_position` immediately after `find_overlapping`
        // returned `Some(key)` for this key, and we hold an
        // exclusive borrow on `self.records` from that call site
        // through here. No path between the find and this lookup
        // removes the entry.
        let rec = records
            .get_mut(&key)
            .expect("widen called on absent record");
        let old_end = rec.footprint_end_exclusive();
        if new_end_exclusive <= old_end {
            return Ok(false);
        }
        let extra_len = new_end_exclusive - old_end;
        if (new_end_exclusive - rec.pos) > max_record_span {
            return Err(WalkerError::RecordTooWide {
                chrom_id: rec.chrom_id,
                pos: rec.pos,
                span: new_end_exclusive - rec.pos,
                cap: max_record_span,
            });
        }
        reference
            .fetch_into(
                ContigId(rec.chrom_id),
                u64::from(old_end),
                u64::from(extra_len),
                widen_bases_buf,
            )
            .map_err(|source| WalkerError::Fasta {
                chrom_id: rec.chrom_id,
                start: old_end,
                start_plus_len: new_end_exclusive,
                source,
            })?;
        let extra_bases = &*widen_bases_buf;

        // Rewrite each existing allele by appending the new REF
        // bases to its `seq`. `alleles[0]` is REF, so this loop is
        // also where the canonical REF sequence widens.
        //
        // **Invariant preserved by `apply_events_to_ref`:** every
        // allele's `seq` ends with a ref-aligned suffix — concretely,
        // the bytes `ref_seq[k..ref_seq.len()]` for some `k` ≥ the
        // offset of the last event the read had inside this record.
        // The function is by construction: after processing each
        // event it advances `ref_cursor` past the event's footprint,
        // and the post-event tail loop emits `ref_seq[ref_cursor..]`
        // verbatim. So whatever modifications (SNP, INS, DEL) sit in
        // the allele, the *trailing* bytes are always ref bases.
        //
        // Widening extends `ref_seq` with `extra_bases`. The
        // ref-aligned tail of every allele then becomes
        // `ref_seq_old[k..] ++ extra_bases`, which is still
        // ref-aligned (now to the wider `ref_seq`). Appending
        // `extra_bases` to each allele's `seq` therefore reproduces
        // exactly what `apply_events_to_ref` would emit if we
        // re-folded the read against the wider `ref_seq`, modulo
        // the new bytes never being event-modified by this read
        // (events at the new positions, if any, haven't been folded
        // yet — they are processed at a later walker step which
        // re-folds the read into the right bucket).
        //
        // So this loop is correct for all allele kinds: SNP/MNP
        // (`seq.len() == old_ref_seq.len()`), DEL (shorter), INS
        // (longer than old span by the inserted run). Each gets
        // exactly the same `extra_bases` appended.
        for allele in &mut rec.alleles {
            allele.seq.extend_from_slice(extra_bases);
        }
        Ok(true)
    }

    /// Open a fresh record at `pos` with REF span `span`, fetching
    /// the reference bases from `reference`.
    fn open_new(
        &mut self,
        chrom_id: u32,
        pos: u32,
        span: u32,
        reference: &dyn RefSeq,
    ) -> Result<&mut OpenPileupRecord, WalkerError> {
        if span > self.max_record_span {
            return Err(WalkerError::RecordTooWide {
                chrom_id,
                pos,
                span,
                cap: self.max_record_span,
            });
        }
        // A fresh `Vec` rather than the table's scratch: these bytes become the
        // record's REF allele and are kept for the record's whole lifetime, so
        // fetching into a reused buffer would only add a copy out of it.
        let mut ref_seq = Vec::new();
        reference
            .fetch_into(
                ContigId(chrom_id),
                u64::from(pos),
                u64::from(span),
                &mut ref_seq,
            )
            .map_err(|source| WalkerError::Fasta {
                chrom_id,
                start: pos,
                start_plus_len: pos + span,
                source,
            })?;
        let rec = OpenPileupRecord::new(chrom_id, pos, ref_seq);
        self.records.insert(pos, rec);
        // PANIC-FREE: the entry was just inserted on the previous
        // line; no concurrent mutation is possible because `&mut self`.
        Ok(self.records.get_mut(&pos).expect("just inserted"))
    }
}

/// Build the haplotype one read presents under an open record, **emitting only
/// what its events cover**, and return the extent they covered.
///
/// # What changed, and why it is the reason this whole step exists
///
/// Production's `apply_events_to_ref_into` emits a reference byte for every
/// offset no event covered — between events, and past the last one
/// ([open_record.rs:522-531](../../../../src/pileup/walker/open_record.rs#L522),
/// [:568-573](../../../../src/pileup/walker/open_record.rs#L568)). At a six-base
/// deletion locus a read that saw only the first two bases is folded as a full
/// witness of a six-base reference haplotype it never saw. **ng emits nothing for
/// a position no event covered** (spec §4, §6).
///
/// **Reads whose events tile the footprint come out byte for byte as before** —
/// there are no gaps to fill — which is what keeps the complete class
/// parity-comparable to production and makes every other divergence a measured
/// one rather than an accident.
///
/// # `ref_seq` is still needed, and not for filling
///
/// An `Insertion`/`Deletion` arm emits the **anchor base** from the reference
/// when no `Match` already emitted it. Normally the `Match` is there and the
/// read's own base wins, so nothing is borrowed; the corner is a read whose base
/// at the anchor was dropped — `N` or adaptor-masked — while its indel at that
/// position was still emitted. It witnessed *the indel* but not *the anchor
/// base's identity*. **Recorded as a known residual, not fixed** (spec §4): it is
/// one base inside an event the read genuinely witnessed, and discarding an
/// observed indel over a masked anchor loses more than it saves.
///
/// # The extent, and the trap in computing it
///
/// The returned span is the union of the events' **reference footprints**, in
/// absolute coordinates, **intersected with `[record_pos, record_end)`**.
/// `events_overlapping` clips a `Match` to the window but returns a `Deletion`
/// **whole** whenever its footprint intersects — so a deletion anchored before
/// the record comes back with its full run, and an unclipped union would put the
/// extent's start below `record_pos` at exactly the long-deletion loci this
/// change exists to fix (spec §8).
///
/// A `Deletion` witnesses every position it deletes, not just its anchor: the
/// read is evidence that those bases are absent. That is why `footprint_span()`
/// is `deleted_len + 1` and why `bases.len()` is **not** the number of positions
/// covered — an insertion adds bases without positions, a deletion positions
/// without bases (spec §8, §13).
///
/// # `None` means the witness is non-contiguous
///
/// An interior `N`, or a ref-skip, leaves a hole in the middle of the run. One
/// `Observed { offset_in_locus, positions_covered }` cannot describe two runs
/// honestly, so the read yields **no observation** and is counted in
/// `reads_without_observation` (spec §6). Rare by construction — adaptor masking
/// and the dropped-indel rule always truncate from one side, so they stay
/// expressible.
///
/// **Preconditions on `events`** (every caller in this module
/// satisfies them; new callers must too):
/// 1. Every event's anchor lies inside the record's footprint
///    (`event.anchor_pos() >= record_pos`), **except** a `Deletion`, which
///    `events_overlapping` may return anchored before the record.
/// 2. Events are sorted by anchor position, non-decreasing.
/// 3. At a tied anchor, events appear in the order
///    `Match → Insertion → Deletion`. The Match must come first
///    so its read base is emitted at the anchor offset before the
///    Insertion/Deletion's `consumed_until` guard suppresses the
///    REF base. The cursor's CIGAR walk satisfies this by
///    construction (the M op preceding an I/D op emits its last
///    Match before the I/D's anchored event).
pub(super) fn apply_events_into(
    allele_seq: &mut Vec<u8>,
    record_pos: u32,
    ref_seq: &[u8],
    events: &[ReadEvent],
) -> Option<RefSpan> {
    allele_seq.clear();
    allele_seq.reserve(ref_seq.len() + 8);

    debug_assert!(
        events.windows(2).all(|w| {
            let a = (w[0].anchor_pos(), event_kind_rank(&w[0]));
            let b = (w[1].anchor_pos(), event_kind_rank(&w[1]));
            a <= b
        }),
        "apply_events_into: events must be sorted by (anchor, Match<Insertion<Deletion); got {events:?}",
    );

    // Skip indices already consumed by an event so an indel does not re-emit an
    // anchor base a Match has already put down.
    let mut consumed_until: u32 = 0; // ref offset (exclusive) consumed by the last event
    let ref_len = ref_seq.len() as u32;
    let record_end = record_pos.saturating_add(ref_len); // exclusive

    // The witnessed run so far, in absolute reference coordinates, half-open.
    // `None` until the first event contributes a position.
    let mut witnessed: Option<(u32, u32)> = None;

    for ev in events {
        debug_assert!(
            ev.anchor_pos() >= record_pos || matches!(ev, ReadEvent::Deletion { .. }),
            "apply_events_into: event anchor {} below record_pos {}",
            ev.anchor_pos(),
            record_pos,
        );
        let offset = ev.anchor_pos().saturating_sub(record_pos);

        // The positions this event witnesses, clipped to the record. Clipping is
        // the `Deletion` trap above; without it a deletion anchored before the
        // record drags the extent's start below `record_pos`.
        let event_start = ev.anchor_pos().max(record_pos);
        let event_end = ev
            .anchor_pos()
            .saturating_add(ev.footprint_span())
            .min(record_end);
        if event_start < event_end {
            match witnessed {
                None => witnessed = Some((event_start, event_end)),
                Some((run_start, run_end)) => {
                    if event_start > run_end {
                        // A hole: the read said nothing about the positions in
                        // between. One run cannot describe that honestly.
                        allele_seq.clear();
                        return None;
                    }
                    witnessed = Some((run_start, run_end.max(event_end)));
                }
            }
        }

        match ev {
            ReadEvent::Match { base, .. } => {
                if offset < ref_len {
                    allele_seq.push(*base);
                }
                consumed_until = consumed_until.max(offset + 1);
            }
            ReadEvent::Insertion { seq, .. } => {
                // Insertion sits AFTER the anchor base; the anchor
                // base itself is unchanged (in the ref-matching
                // sense). Emit the anchor base if not already
                // emitted, then append the inserted bases.
                if offset < ref_len && offset >= consumed_until {
                    allele_seq.push(ref_seq[offset as usize]);
                }
                allele_seq.extend_from_slice(seq);
                consumed_until = consumed_until.max(offset + 1);
            }
            ReadEvent::Deletion {
                anchor_ref_pos,
                deleted_len,
                ..
            } => {
                // DEL: keep the anchor base, drop the next
                // `deleted_len` reference bases.
                let anchor_inside = *anchor_ref_pos >= record_pos;
                if anchor_inside && offset < ref_len && offset >= consumed_until {
                    allele_seq.push(ref_seq[offset as usize]);
                }
                let skip_until = anchor_ref_pos
                    .saturating_add(1)
                    .saturating_add(*deleted_len)
                    .saturating_sub(record_pos);
                consumed_until = consumed_until.max(skip_until);
            }
        }
    }

    // Half-open on the way in, inclusive on the way out — `RefSpan`'s convention.
    // `end > start` holds for every run that got here, so the `- 1` cannot wrap.
    witnessed.map(|(start, end)| RefSpan {
        start,
        end: end - 1,
    })
}

/// Owning wrapper around [`apply_events_into`]. Kept so unit tests (and any
/// caller that isn't on the hot path) don't have to hoist a buffer manually. The
/// hot-path call in `process_position` goes through the `_into` variant against a
/// buffer on `OpenPileupRecordTable`.
///
/// Returns the bases **and** the extent, or `None` for a non-contiguous witness —
/// a caller that wanted only the bases would be able to ignore the very
/// distinction this step introduces.
#[cfg(test)]
pub(super) fn apply_events(
    record_pos: u32,
    ref_seq: &[u8],
    events: &[ReadEvent],
) -> Option<(Vec<u8>, RefSpan)> {
    let mut out = Vec::new();
    let witnessed = apply_events_into(&mut out, record_pos, ref_seq, events)?;
    Some((out, witnessed))
}

/// Order rank for the tied-anchor precondition on
/// `apply_events_to_ref`: at the same anchor, Match must come
/// before Insertion which must come before Deletion. The cursor's
/// CIGAR walk produces this order naturally; the rank is just a
/// total order for the debug-assert.
fn event_kind_rank(ev: &ReadEvent) -> u8 {
    match ev {
        ReadEvent::Match { .. } => 0,
        ReadEvent::Insertion { .. } => 1,
        ReadEvent::Deletion { .. } => 2,
    }
}

/// Find or create the allele bucket inside a record matching `seq`.
/// Returns the bucket index. Linear scan is fine — records
/// typically carry ≤ a few alleles.
///
/// Takes `&mut Vec<OpenAllele>` rather than `&mut OpenPileupRecord`
/// so the caller can hold an independent borrow on the record's
/// `ref_seq` (e.g. for `apply_events_to_ref`) at the same time.
/// This is what lets `process_position` avoid cloning `ref_seq` per
/// affected record. Mi6 in `ia/reviews/pileup_2026-05-09.md`.
#[cfg(test)]
pub(super) fn find_or_create_allele_index(alleles: &mut Vec<OpenAllele>, seq: Vec<u8>) -> usize {
    if let Some(idx) = alleles.iter().position(|a| a.seq == seq) {
        idx
    } else {
        alleles.push(OpenAllele::new(seq));
        alleles.len() - 1
    }
}

/// Borrowed lookup variant: matches `seq` against existing alleles
/// without taking ownership. Used by the hot-path fold so the
/// caller can run an equality check against a reusable buffer and
/// only `clone()` the bytes when a genuinely new allele is added.
/// L10 in `ia/reviews/perf_pileup_2026-05-10.md`.
pub(super) fn find_allele_index(alleles: &[OpenAllele], seq: &[u8]) -> Option<usize> {
    alleles.iter().position(|a| a.seq.as_slice() == seq)
}

/// Process all events at `walker_pos` from the given list of
/// per-read contributions. Each contributor MUST have at least
/// one event at `walker_pos` (anchor == walker_pos); silent reads
/// (e.g. inside their own deletion or an N-skip at this pos) are
/// already filtered out by the caller.
///
/// 1. Step 3 — identify candidate records per event. Each event
///    either merges into an existing overlapping record or opens
///    a fresh one. Record keys touched here go into `affected`.
///
/// 2. Step 4-6 — for each affected record, fold every contributor
///    whose events overlap the record's footprint exactly once
///    into the matching allele bucket. The fold uses the
///    contributor's full event list (not just events_at_pos) so
///    that compound alleles spanning the whole record's footprint
///    collapse into a single allele bucket.
///
/// Records that already existed at walker_pos but were *not*
/// affected at this step are not re-folded — those reads were
/// folded at the walker step where the record was created or
/// widened, and folding again here would double-count.
pub(super) fn process_position(
    open: &mut OpenPileupRecordTable,
    walker_pos: u32,
    chrom_id: u32,
    contributors: &[ReadContribution],
    active_reads: &ActiveReads,
    reference: &dyn RefSeq,
) -> Result<ProcessOutcome, WalkerError> {
    let mut affected: Vec<u32> = Vec::new();
    let mut widen_count: u64 = 0;

    // Step 3: each event either lands in an existing record (and
    // possibly widens it) or opens a fresh one.
    for contrib in contributors {
        for ev in &contrib.events_at_pos {
            let event_start = ev.anchor_pos();
            // `saturating_add` for `event_end` per Mi8: on
            // multi-Gbp chromosomes a raw `+` would wrap and the
            // `find_overlapping` range lookup would search the
            // wrong region.
            let event_end = event_start.saturating_add(ev.footprint_span());

            let key = if let Some(k) = open.find_overlapping(event_start, event_end) {
                // PANIC-FREE: `find_overlapping` returned `Some(k)`
                // on the previous line; no mutation between then
                // and this lookup.
                let cur_end = open
                    .records
                    .get(&k)
                    .expect("just located")
                    .footprint_end_exclusive();
                if event_end > cur_end && open.widen(k, event_end, reference)? {
                    widen_count += 1;
                }
                k
            } else {
                let new = open.open_new(chrom_id, event_start, ev.footprint_span(), reference)?;
                new.pos
            };
            if !affected.contains(&key) {
                affected.push(key);
            }
        }
    }

    // Step 4-6: for each affected record (in coordinate order),
    // fold each contributor that has events overlapping the
    // record's footprint. Each (record, contributor) pair folds
    // exactly once across the record's lifetime: re-folds at
    // later walker steps subtract the prior contribution from
    // the old bucket before adding the new one.
    affected.sort_unstable();
    // Split-borrow the table fields so the inner fold can hold a
    // mutable borrow on `allele_seq_buf` simultaneously with a
    // mutable borrow on a record inside `records`. Equivalent to
    // the `OpenPileupRecord { alleles, folded_reads, .. }` split
    // below, just one level up.
    let OpenPileupRecordTable {
        records,
        allele_seq_buf,
        closing_keys_buf: _,
        widen_bases_buf: _,
        max_record_span: _,
    } = open;
    for key in affected {
        // PANIC-FREE: every key in `affected` was either inserted
        // by `open_new` or returned by `find_overlapping` in the
        // step-3 loop above; no path between that loop and here
        // removes records.
        let rec = records.get_mut(&key).expect("affected key must exist");
        let rec_pos = rec.pos;
        // Destructure through `&mut` so `alleles` and `folded_reads`
        // are independent mutable borrows of disjoint fields. The
        // REF sequence is read from `alleles[0].seq` — that borrow
        // ends inside each `apply_events_to_ref` call (which returns
        // an owned `Vec<u8>`), leaving the rest of the inner loop
        // free to mutate `alleles[other_idx]`. This is the same
        // disjoint-borrow shape Mi6 (`pileup_2026-05-09.md`)
        // introduced, with `ref_seq` collapsed onto `alleles[0]` —
        // the bytes are no longer duplicated on the record.
        let OpenPileupRecord {
            alleles,
            folded_reads,
            ..
        } = rec;
        let rec_end = rec_pos + alleles[0].seq.len() as u32;

        for contrib in contributors {
            // PANIC-FREE: every `ReadContribution` is built from
            // an entry in the active set in `walker::process_position`
            // (the same step that produces `contributors`), and
            // `process_position` runs before `expire_passed_reads`,
            // so the read_id is still in the active set here.
            let active_read = active_reads
                .get_by_read_id(contrib.read_id)
                .expect("contributor's read_id must still be in the active set");
            let mut window_events =
                active_read
                    .cursor
                    .events_overlapping(rec_pos, rec_end, &active_read.read);
            if contrib.bq_zero_in_window {
                for ev in &mut window_events {
                    zero_event_bq(ev);
                }
            }
            if let Some(override_bq) = contrib.bq_override_at_walker_pos {
                // S7 mate-overlap math: the agree-case keeper's
                // walker_pos event carries the summed BQ; the
                // disagree-case winner's carries the 0.8-scaled
                // BQ. Apply only at the walker_pos anchor —
                // window events at other positions keep cursor-
                // original BQ.
                for ev in &mut window_events {
                    if ev.anchor_pos() == walker_pos
                        && let ReadEvent::Match { bq_baq, .. } = ev
                    {
                        *bq_baq = override_bq;
                    }
                }
            }

            // A contributor only folds into a record if it has
            // events overlapping the record's footprint. (No
            // events overlapping = the read doesn't observe this
            // record's REF stretch at all — it shouldn't fold.)
            if window_events.is_empty() {
                continue;
            }

            let Some(witnessed) =
                apply_events_into(allele_seq_buf, rec_pos, &alleles[0].seq, &window_events)
            else {
                // A non-contiguous witness — an interior `N`, or a ref-skip —
                // yields no observation (spec §6). **The prior contribution has
                // to come off first**, and this is not a corner: a read that
                // folded contiguously *becomes* non-contiguous when the window
                // widens right across an interior gap, so a bare `continue` here
                // would strand a live contribution in a bucket for a read that
                // now has no row, breaking `chain_ids.len() <= num_obs` silently
                // and only on multi-base records (spec §4).
                //
                // A5 adds the other half — recording *which* reads these were,
                // as a per-record set of read ids rather than a counter, since
                // this path is reached at every position the record is affected
                // at and a counter would multiply by the footprint length.
                if let Some(prev) = folded_reads.remove(&contrib.read_id) {
                    subtract_contribution(
                        &mut alleles[prev.allele_index].support,
                        &prev.contribution,
                    );
                }
                continue;
            };
            let bq = ln_bq_for_read(&window_events, contrib.bq_baq_at_walker_pos);
            let ln_q = bq.max(contrib.mq_log_err);
            let mapq = u32::from(contrib.mapq);
            let new_contribution = AlleleSupportStats {
                num_obs: 1,
                q_sum: ln_q,
                fwd: u32::from(!contrib.is_reverse_strand),
                placed_left: u32::from(contrib.alignment_start < rec_pos),
                mapq_sum: mapq,
                mapq_sum_sq: (mapq as u64) * (mapq as u64),
            };

            // Subtract any prior contribution this read had to
            // this record. Re-folds happen when a record widens
            // and the read's haplotype under the wider span
            // either changes bucket or stays in the same bucket
            // with a possibly-different ln_q.
            if let Some(prev) = folded_reads.remove(&contrib.read_id) {
                subtract_contribution(&mut alleles[prev.allele_index].support, &prev.contribution);
            }
            // Borrowed lookup; only `clone()` the bytes when adding
            // a genuinely new allele. In SNP/REF steady state every
            // contributor lands in the same existing bucket, so
            // the clone path never fires.
            let new_index = match find_allele_index(alleles, allele_seq_buf) {
                Some(idx) => idx,
                None => {
                    alleles.push(OpenAllele::new(allele_seq_buf.clone()));
                    alleles.len() - 1
                }
            };
            add_contribution(&mut alleles[new_index].support, &new_contribution);
            folded_reads.insert(
                contrib.read_id,
                FoldedReadState {
                    allele_index: new_index,
                    contribution: new_contribution,
                    chain_id: contrib.chain_id,
                    read_group: contrib.read_group,
                    witnessed,
                    placed_start: contrib.alignment_start == rec_pos,
                },
            );
        }
    }

    Ok(ProcessOutcome { widen_count })
}

/// Side-effect counters returned by `process_position` so the
/// walker can update its run-summary fields without having to
/// inspect open-record-table state externally.
#[derive(Debug, Default, Clone, Copy)]
pub(super) struct ProcessOutcome {
    /// Number of records that actually widened during this call.
    /// Excludes fresh `open_new` calls and re-finds where the
    /// event already fits inside the record's footprint.
    pub widen_count: u64,
}

/// Zero a `Match` event's `bq_baq` in place; leave indel
/// `bq_proxy` untouched. Used by the open-record fold when a
/// contributor is flagged as the match-only mate-overlap loser:
/// the rule (`pileup_walker.md` §"Mate overlap") zeros the
/// Match at the overlap *position*, not all events in the
/// loser's haplotype window. Indels in the same window sit at
/// different anchors and carry independent evidence.
///
/// In today's pipeline the narrowing is invisible — the
/// haplotype-level `min` in `ln_bq_for_read` collapses to 0 via
/// the zeroed Match regardless. The narrow version is correct
/// against a wider class of future reductions (median, weighted)
/// and matches the spec wording. Mi3 in
/// `ia/reviews/pileup_2026-05-09.md`.
fn zero_event_bq(ev: &mut ReadEvent) {
    if let ReadEvent::Match { bq_baq, .. } = ev {
        *bq_baq = 0;
    }
}

fn add_contribution(support: &mut AlleleSupportStats, c: &AlleleSupportStats) {
    support.num_obs += c.num_obs;
    support.q_sum += c.q_sum;
    support.fwd += c.fwd;
    support.placed_left += c.placed_left;
    support.mapq_sum += c.mapq_sum;
    support.mapq_sum_sq += c.mapq_sum_sq;
}

fn subtract_contribution(support: &mut AlleleSupportStats, c: &AlleleSupportStats) {
    // Saturating-style: an internal-bookkeeping bug that produced
    // a negative would otherwise wrap silently. Saturate to zero
    // in release builds and rely on the upper invariant (num_obs
    // over the record total still adds up) to surface mistakes
    // through tests.
    //
    // `debug_assert!` peer per Mi7 in
    // `ia/reviews/pileup_2026-05-11.md`: in debug builds (and
    // tests) we trip loudly on an underflow rather than letting
    // it saturate silently. `q_sum` is signed `f64` and is left
    // as a raw subtract — negatives are legal there by design.
    debug_assert!(
        support.num_obs >= c.num_obs,
        "subtract_contribution underflow on num_obs: {} -= {}",
        support.num_obs,
        c.num_obs,
    );
    debug_assert!(
        support.fwd >= c.fwd,
        "subtract_contribution underflow on fwd: {} -= {}",
        support.fwd,
        c.fwd,
    );
    debug_assert!(
        support.placed_left >= c.placed_left,
        "subtract_contribution underflow on placed_left: {} -= {}",
        support.placed_left,
        c.placed_left,
    );
    debug_assert!(
        support.mapq_sum >= c.mapq_sum,
        "subtract_contribution underflow on mapq_sum: {} -= {}",
        support.mapq_sum,
        c.mapq_sum,
    );
    debug_assert!(
        support.mapq_sum_sq >= c.mapq_sum_sq,
        "subtract_contribution underflow on mapq_sum_sq: {} -= {}",
        support.mapq_sum_sq,
        c.mapq_sum_sq,
    );
    support.num_obs = support.num_obs.saturating_sub(c.num_obs);
    support.q_sum -= c.q_sum;
    support.fwd = support.fwd.saturating_sub(c.fwd);
    support.placed_left = support.placed_left.saturating_sub(c.placed_left);
    support.mapq_sum = support.mapq_sum.saturating_sub(c.mapq_sum);
    support.mapq_sum_sq = support.mapq_sum_sq.saturating_sub(c.mapq_sum_sq);
}

/// Per-read BQ for an allele's quality contribution: min over
/// the read's events in the record's footprint, converted to
/// `ln(P_err)`. Confirmed against freebayes
/// `AlleleParser.cpp:3151-3155` (haplotype-allele construction
/// uses `min(quality)`). For a contributor with no events in the
/// window (a clean REF read), we fall back to the read's BQ at
/// the walker_pos, which is a `Match` event's `bq_baq` already
/// stamped on the contribution.
fn ln_bq_for_read(window_events: &[ReadEvent], fallback_bq: u8) -> f64 {
    let min_bq = window_events
        .iter()
        .map(|e| match e {
            ReadEvent::Match { bq_baq, .. } => *bq_baq,
            ReadEvent::Insertion { bq_proxy, .. } => *bq_proxy,
            ReadEvent::Deletion { bq_proxy, .. } => *bq_proxy,
        })
        .min()
        .unwrap_or(fallback_bq);
    phred_to_ln_perr(min_bq)
}

/// `Q -> ln(P_err)` where `P_err = 10^(-Q/10)`.
///
/// `q` is a Phred score constrained to `0..=93` by the
/// [`PreparedRead`] spec, but the table is sized at 256 so the
/// `q as usize` index covers `u8`'s full domain — bounds-check
/// elision is unconditional. Built at compile time; the 2 KB
/// table sits in L1 and the function compiles to one load.
/// `q == 0` is pinned to `+0.0` to match the prior branch's
/// `return 0.0;` exactly.
///
/// H3 in `ia/reviews/perf_pileup_2026-05-12.md`. Replaces the
/// per-call FP multiply that round-2 perf attributed at 2.94 %
/// of walker self-time on `pileup_walker_multi_op/5000`.
fn phred_to_ln_perr(q: u8) -> f64 {
    static LN_PERR_TABLE: [f64; 256] = {
        let mut t = [0.0_f64; 256];
        let mut q = 1usize;
        while q < 256 {
            t[q] = -(q as f64) * std::f64::consts::LN_10 / 10.0;
            q += 1;
        }
        t
    };
    LN_PERR_TABLE[q as usize]
}

/// One read's contribution to the current walker position. The
/// walker assembles these from its active set before calling
/// `process_position`. The contribution carries no event window —
/// the fold pulls window events lazily from the read's
/// `CigarCursor` via `read_id` lookup against the `ActiveReads`.
#[derive(Debug, Clone)]
pub(super) struct ReadContribution {
    /// Active-set local id of the contributing read. Keys the
    /// per-record `folded_reads` map (so re-folds subtract the
    /// prior contribution) and is also how the fold looks the
    /// read up against the `ActiveReads` to query window events.
    pub read_id: u32,
    pub chain_id: ChainId,
    /// The read group this read was prepared in, copied off the active read like
    /// the fields around it — ng's [`PreparedRead`](super::PreparedRead) carries
    /// it, so nothing here reconstructs it from a side channel. Unread until B1
    /// makes it part of the observation's identity; carried to
    /// [`FoldedReadState`] from A1 so the fold has it when that lands.
    pub read_group: ReadGroupId,
    /// Events whose anchor *is* this walker_pos (used by step 3
    /// to identify candidate records). At most 2 events anchor at
    /// any walker_pos (one Match plus at most one indel), so the
    /// SmallVec keeps this list off the heap on every step.
    pub events_at_pos: super::cigar_cursor::EventsAt,
    /// BAQ-capped BQ at this walker position (the Match-event's
    /// quality, used as the fallback when the cursor's window
    /// returns no events — a clean REF read).
    pub bq_baq_at_walker_pos: u8,
    pub mq_log_err: f64,
    /// Raw BAM MAPQ for this contributor, forwarded from
    /// [`PreparedRead::mapq`] so the fold can accumulate
    /// `mapq_sum` / `mapq_sum_sq` on the chosen allele bucket.
    pub mapq: u8,
    pub is_reverse_strand: bool,
    pub alignment_start: u32,
    /// SAM-flag-derived mate role. Carries through from
    /// [`PreparedRead::mate_role`]; only
    /// [`MateRole::is_first_of_pair`] is read by the tie-break
    /// helpers (a deterministic tertiary key on equal-BQ
    /// mate-overlap positions per the spec — not a faithful proxy
    /// for "earlier alignment_start").
    pub mate_role: super::MateRole,
    /// Set by `resolve_mate_overlap_at_pos` when this contributor
    /// is the loser of a Match-only mate-overlap. The fold zeroes
    /// the BQ on every event it pulls from this contributor's
    /// cursor — equivalent to the eager design's in-place mutation
    /// of the cloned full-window event list.
    pub bq_zero_in_window: bool,
    /// Set by `resolve_mate_overlap_at_pos` (S7) on the surviving
    /// side of a Match-only mate-overlap pair: the agree-case
    /// keeper carries the *summed* BQ (capped at 200), and the
    /// disagree-case winner carries its BQ scaled by 0.8. The fold
    /// rewrites the BQ on any window event whose anchor is
    /// walker_pos to this value. Other window-event positions
    /// keep their cursor-original BQ — the override is per-position,
    /// matching our per-walker_pos resolution model rather than
    /// samtools' admission-time mutation of the entire overlap
    /// region.
    pub bq_override_at_walker_pos: Option<u8>,
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::pileup::tests::MockFasta;

    fn fa(s: &str) -> MockFasta {
        MockFasta::new(s)
    }

    #[test]
    fn open_new_creates_record_with_ref_allele_zero_obs() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTAC");
        let rec = t.open_new(0, 1, 1, &f).unwrap();
        assert_eq!(rec.pos, 1);
        assert_eq!(rec.alleles.len(), 1);
        assert_eq!(rec.alleles[0].seq, b"A");
        assert_eq!(rec.alleles[0].support.num_obs, 0);
    }

    #[test]
    fn widen_extends_ref_seq_and_existing_alleles() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTAC");
        // Open at pos 1 with span 1 ("A").
        t.open_new(0, 1, 1, &f).unwrap();
        // Now widen to span 3 ("ACG").
        t.widen(1, 4, &f).unwrap();
        let rec = t.records.get(&1).unwrap();
        assert_eq!(rec.alleles[0].seq, b"ACG");
    }

    #[test]
    fn drain_aged_emits_in_coordinate_order() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTACGTAC");
        t.open_new(0, 1, 1, &f).unwrap();
        t.open_new(0, 5, 1, &f).unwrap();
        t.open_new(0, 8, 1, &f).unwrap();
        // Walker at pos 6 — record at 1 (ends 2) and at 5 (ends 6) are aged out.
        let mut drained = Vec::new();
        t.drain_aged_into(6, &mut drained);
        assert_eq!(drained.len(), 2);
        assert_eq!(drained[0].pos, 1);
        assert_eq!(drained[1].pos, 5);
        assert_eq!(t.records.len(), 1);
    }

    #[test]
    fn drain_aged_does_not_break_early_when_a_wide_record_blocks_a_narrow_one() {
        // Heterogeneous spans: a wide deletion record sits at an
        // earlier key than a narrow record. At a walker position
        // where the narrow record is past closure but the wide one
        // is not, `drain_aged` must still return the narrow record.
        // An early break at the first "not yet aged" key would
        // miss it.
        let mut t = OpenPileupRecordTable::new();
        let f = fa(&"A".repeat(60));
        // Wide record at pos 5, span 50 → footprint [5, 55).
        t.open_new(0, 5, 50, &f).unwrap();
        // Narrow record at pos 10, span 1 → footprint [10, 11).
        t.open_new(0, 10, 1, &f).unwrap();
        // Walker at 11: narrow record should be aged out (11 >= 11);
        // wide record is still open (55 > 11).
        let mut drained = Vec::new();
        t.drain_aged_into(11, &mut drained);
        assert_eq!(
            drained.len(),
            1,
            "narrow record must drain even though the earlier wide one is still open",
        );
        assert_eq!(drained[0].pos, 10);
        assert_eq!(t.records.len(), 1, "wide record stays open");
    }

    #[test]
    fn find_overlapping_returns_record_when_event_falls_inside_footprint() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("AAAACCCCGGGG");
        // Open a deletion-shaped record at pos 1 with span 5
        // (footprint [1, 6)).
        t.open_new(0, 1, 5, &f).unwrap();
        // Event at pos 3 with span 1 (a SNP) — overlaps [1, 6).
        let key = t.find_overlapping(3, 4);
        assert_eq!(key, Some(1));
        // Event at pos 6 with span 1 — does NOT overlap (touching, not overlapping).
        let key = t.find_overlapping(6, 7);
        assert_eq!(key, None);
    }

    #[test]
    fn find_overlapping_walks_past_intermediate_narrow_record_to_wide_one() {
        // Two open records:
        //   - wide deletion at pos 5, span 50 → footprint [5, 55)
        //   - narrow SNP at pos 40, span 1 → footprint [40, 41)
        // Query for an event at pos 41 with span 1: the wide
        // record's footprint reaches into 41, so the answer must
        // be Some(5). The previous early-break terminated at the
        // narrow record because its footprint ends at 41 (≤ 41),
        // missing the wide record entirely.
        let mut t = OpenPileupRecordTable::new();
        let chrom = "A".repeat(60);
        let f = fa(&chrom);
        t.open_new(0, 5, 50, &f).unwrap();
        t.open_new(0, 40, 1, &f).unwrap();
        let key = t.find_overlapping(41, 42);
        assert_eq!(key, Some(5), "wide record at 5 must be found");
    }

    /// **The builder no longer fills.** No events means the read witnessed nothing
    /// inside this record — production emitted the whole reference sequence here, as
    /// if the read had seen every base of it.
    ///
    /// The caller never reaches this: `process_position` skips a contributor whose
    /// window is empty. It is pinned anyway because it is the fabrication in its
    /// purest form, and a re-introduced fill would show up here first.
    #[test]
    fn apply_events_with_no_events_witnesses_nothing() {
        let ref_seq = b"ACGTA";
        assert!(apply_events(100, ref_seq, &[]).is_none());
    }

    /// One `Match` inside a five-base record: the read witnessed **one** position and
    /// contributes one base — where production emitted `ACXTA`, four bases of which
    /// it invented.
    #[test]
    fn apply_events_a_lone_match_contributes_only_the_base_it_saw() {
        let ref_seq = b"ACGTA";
        let snp = ReadEvent::Match {
            ref_pos: 102,
            base: b'X',
            bq_baq: 30,
        };
        let (bases, witnessed) =
            apply_events(100, ref_seq, std::slice::from_ref(&snp)).expect("one run");
        assert_eq!(bases, b"X");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 102,
                end: 102
            }
        );
    }

    /// Events tiling the whole footprint are **byte-identical to production** — there
    /// are no gaps to fill. This is the class the permanent differential anchors on,
    /// and the reason the change can be measured rather than merely asserted.
    #[test]
    fn apply_events_tiling_the_footprint_reproduces_productions_bytes() {
        let ref_seq = b"ACGTA";
        let events: Vec<ReadEvent> = (0..5)
            .map(|k| ReadEvent::Match {
                ref_pos: 100 + k,
                base: if k == 2 { b'X' } else { ref_seq[k as usize] },
                bq_baq: 30,
            })
            .collect();
        let (bases, witnessed) = apply_events(100, ref_seq, &events).expect("one run");
        assert_eq!(
            bases, b"ACXTA",
            "exactly what production emits for this read"
        );
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 104
            }
        );
    }

    /// **A deletion witnesses every position it deletes**, not just its anchor: the
    /// read is evidence those bases are absent. So `bases.len()` (2) is not
    /// `positions_covered` (3) — the §13 consistency check is an inequality, and an
    /// implementer who derives the extent from the byte length has rebuilt the
    /// span-versus-events confusion this step exists to remove.
    #[test]
    fn apply_events_a_deletion_witnesses_the_positions_it_deleted() {
        let ref_seq = b"ACGTA";
        let del = ReadEvent::Deletion {
            anchor_ref_pos: 100,
            deleted_len: 2,
            bq_proxy: 30,
        };
        let (bases, witnessed) =
            apply_events(100, ref_seq, std::slice::from_ref(&del)).expect("one run");
        // The anchor base only — 101 and 102 are deleted, and 103/104 were never
        // witnessed, where production appended them from the reference ("ATA").
        assert_eq!(bases, b"A");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 102
            },
            "anchor plus the two deleted positions"
        );
    }

    /// A deletion **anchored before the record** contributes none of the bases it deleted,
    /// and no anchor base — the anchor belongs to another record. **And the extent is
    /// clipped to the record**: `events_overlapping` returns such a deletion *whole*, so
    /// an unclipped union would report a start below `record_pos` at exactly the
    /// long-deletion loci this change exists to fix (spec §8).
    ///
    /// The bases half is production's own defect, fixed in `5f32a62` and inherited here:
    /// the offset was computed with `saturating_sub`, which emitted `ref_seq[0]` — a base
    /// the read had *deleted* — and skipped one position too many. Reachable in
    /// production: a mate-overlap collapse in the indel regime removes the indel-carrying
    /// contributor, so no record opens at the anchor, and a later record opens inside the
    /// deleted run. Seen on real GIAB HG002 data at chr1:106,324,863.
    ///
    /// Here the record starts at 100 and the deletion is anchored at 98, removing 99–101.
    #[test]
    fn apply_events_deletion_anchored_before_the_record_emits_no_anchor_base() {
        let ref_seq = b"ACGTA"; // positions 100..=104
        let del = ReadEvent::Deletion {
            anchor_ref_pos: 98,
            deleted_len: 3,
            bq_proxy: 30,
        };
        let (bases, witnessed) =
            apply_events(100, ref_seq, std::slice::from_ref(&del)).expect("one run");
        // 99, 100, 101 deleted. 100 and 101 lie inside the record and are witnessed as
        // absent; 102–104 the read said nothing about, where production emitted "GTA".
        assert_eq!(bases, b"");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 101
            },
            "clipped to the record: the anchor at 98 is not inside it"
        );
    }

    /// The whole deleted run can also fall *past* the record, so every position of it is
    /// witnessed as absent and no base is contributed.
    #[test]
    fn apply_events_deletion_anchored_before_the_record_can_consume_all_of_it() {
        let ref_seq = b"ACGTA";
        let del = ReadEvent::Deletion {
            anchor_ref_pos: 98,
            deleted_len: 20,
            bq_proxy: 30,
        };
        let (bases, witnessed) =
            apply_events(100, ref_seq, std::slice::from_ref(&del)).expect("one run");
        assert_eq!(bases, b"");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 104
            },
            "clipped at the record's end, not at the deletion's"
        );
    }

    #[test]
    fn apply_events_match_before_insertion_at_same_anchor_emits_read_base_then_inserted() {
        // Tied-anchor case: Match and Insertion both anchored at
        // record_pos. The Match must come first (per the function's
        // precondition) so its read base is emitted; the Insertion
        // then appends its inserted run, with `consumed_until`
        // suppressing the would-be REF anchor push. The cursor's
        // CIGAR walk satisfies this order by construction (M op
        // before I op at the same anchor). Mi9 in
        // `ia/reviews/pileup_2026-05-09.md`.
        let ref_seq = b"ACGTA";
        let m = ReadEvent::Match {
            ref_pos: 100,
            base: b'X',
            bq_baq: 30,
        };
        let i = ReadEvent::Insertion {
            anchor_ref_pos: 100,
            seq: b"YY".to_vec(),
            bq_proxy: 30,
        };
        let (bases, witnessed) = apply_events(100, ref_seq, &[m, i]).expect("one run");
        // Anchor: read base 'X' (not REF 'A'), then inserted "YY". The REF tail
        // production appended is gone — the read never witnessed 101–104.
        assert_eq!(bases, b"XYY");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 100
            },
            "an insertion's footprint is one reference position, however many bases it adds"
        );
    }

    /// **The one residual: a borrowed reference base.** The insertion's anchor base was
    /// not emitted by any `Match` — the read's own base there was dropped as `N` or
    /// adaptor-masked — so the builder supplies the reference base for that one position.
    /// Recorded, not fixed: it is one base inside an event the read genuinely witnessed,
    /// and discarding an observed indel over a masked anchor loses more (spec §4).
    ///
    /// **Pinned at one base and no more:** the tail production appended is gone.
    #[test]
    fn apply_events_an_indel_over_a_masked_anchor_borrows_exactly_one_reference_base() {
        let ref_seq = b"ACGTA";
        let ins = ReadEvent::Insertion {
            anchor_ref_pos: 100,
            seq: b"XX".to_vec(),
            bq_proxy: 30,
        };
        let (bases, witnessed) =
            apply_events(100, ref_seq, std::slice::from_ref(&ins)).expect("one run");
        assert_eq!(bases, b"AXX", "the borrowed 'A', then the inserted run");
        assert_eq!(
            bases.iter().filter(|b| **b == b'A').count(),
            1,
            "exactly one reference base is borrowed — production emitted four more"
        );
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 100
            }
        );
    }

    /// **A hole in the middle yields no observation at all.** An interior `N` or a
    /// ref-skip leaves the read silent about positions inside a run it otherwise
    /// witnessed, and one `Observed` run cannot describe two runs honestly (spec §6).
    ///
    /// Production filled the hole from the reference and folded the read as a complete
    /// witness.
    #[test]
    fn apply_events_a_hole_in_the_middle_yields_no_observation() {
        let ref_seq = b"ACGTA";
        let events = vec![
            ReadEvent::Match {
                ref_pos: 100,
                base: b'A',
                bq_baq: 30,
            },
            // 101 and 102 witnessed by nothing.
            ReadEvent::Match {
                ref_pos: 103,
                base: b'T',
                bq_baq: 30,
            },
        ];
        assert!(apply_events(100, ref_seq, &events).is_none());
    }

    /// **Adjacent is not a hole.** The gap check is `event_start > run_end` on a
    /// half-open run, so two neighbouring positions stay one run — an off-by-one there
    /// would send every ordinary read to the no-observation path.
    #[test]
    fn apply_events_adjacent_events_stay_one_run() {
        let ref_seq = b"ACGTA";
        let events = vec![
            ReadEvent::Match {
                ref_pos: 101,
                base: b'C',
                bq_baq: 30,
            },
            ReadEvent::Match {
                ref_pos: 102,
                base: b'G',
                bq_baq: 30,
            },
        ];
        let (bases, witnessed) = apply_events(100, ref_seq, &events).expect("one run");
        assert_eq!(bases, b"CG");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 101,
                end: 102
            }
        );
    }

    /// A deletion's run closes what would otherwise be a hole: the positions it deletes
    /// are witnessed, so a `Match` on its far side continues the same run.
    #[test]
    fn apply_events_a_deletions_run_bridges_to_the_match_beyond_it() {
        let ref_seq = b"ACGTA";
        let events = vec![
            ReadEvent::Deletion {
                anchor_ref_pos: 100,
                deleted_len: 2,
                bq_proxy: 30,
            },
            ReadEvent::Match {
                ref_pos: 103,
                base: b'T',
                bq_baq: 30,
            },
        ];
        let (bases, witnessed) = apply_events(100, ref_seq, &events).expect("one run");
        assert_eq!(bases, b"AT");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 103
            },
            "the deleted positions are witnessed, so 103 continues the run"
        );
    }

    /// **A read that folded contiguously and then became non-contiguous is taken back
    /// out of the record it had folded into.**
    ///
    /// The `None` path is not reached only on first contact. A record widens, and a read
    /// whose witness *was* one run across the old footprint now has a hole inside the new
    /// one — so it stops having an observation, and the contribution it already made has to
    /// come off the bucket. A bare `continue` there leaves a live contribution behind for a
    /// read that has no row, which is silent, only happens on multi-base records, and is
    /// the failure spec §4 names.
    ///
    /// Nothing else in the suite catches it: the differential's census tolerates any
    /// difference in the allele lists by design, and the inherited tests never widen a
    /// record across a read's interior gap. So it is built by hand.
    ///
    /// ```text
    /// ref     1        5   7    ...              20      25
    /// wide    MMMMMMMMMMNMMMMMMMMMMMMMM     the N at 11 is silent
    /// opener      M--MMMMMMMMMMMMMMMMMM     D2 opens the record at 5, span 5..=7
    /// widener       M-------------MMMM      D13 widens it to 5..=20
    /// ```
    ///
    /// `wide` folds at position 5 over the three-base record (5, 6, 7 — one run). At
    /// position 7 `widener`'s deletion grows that record to sixteen bases, and `wide`'s
    /// witness over it is now `5..=10` and `12..=20`: two runs, so no observation.
    #[test]
    fn a_read_that_becomes_non_contiguous_when_the_record_widens_leaves_its_bucket() {
        use super::super::{CigarOp, MateRole, PreparedRead, WalkerConfig, run};
        use crate::ng::types::ReadGroupId;
        use std::sync::Arc;

        let contig = "ACGTACGTACGTACGTACGTACGTA"; // 25 bases
        let read =
            |qname: &str, start: u32, end: u32, cigar: Vec<CigarOp>, seq: Vec<u8>| PreparedRead {
                chrom_id: 0,
                alignment_start: start,
                alignment_end: end,
                cigar,
                bq_baq: vec![30; seq.len()],
                seq,
                mq_log_err: -3.0,
                mapq: 60,
                is_reverse_strand: false,
                qname: Arc::from(qname),
                mate_role: MateRole::Solo,
                adaptor_boundary: None,
                read_group: ReadGroupId(0),
            };

        // `wide` matches every position 1..=25 except 11, where its base is `N` and the
        // cursor emits nothing — the interior hole, and a hole its *alignment span* is
        // blind to, which is why coverage has to come from the events.
        let mut wide_seq = contig.as_bytes().to_vec();
        wide_seq[10] = b'N';
        let reads = vec![
            read("wide", 1, 25, vec![CigarOp::Match(25)], wide_seq),
            read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(18)],
                vec![b'A'; 19],
            ),
            read(
                "widener",
                7,
                24,
                vec![CigarOp::Match(1), CigarOp::Deletion(13), CigarOp::Match(4)],
                vec![b'A'; 5],
            ),
        ];

        let records: Vec<_> = run(reads, MockFasta::new(contig), &WalkerConfig::default())
            .map(|record| record.expect("the walk succeeds"))
            .collect();
        let widened = records
            .iter()
            .find(|record| record.pos == 5)
            .expect("the opener's deletion opens a record at 5");
        assert_eq!(
            widened.alleles[0].seq.len(),
            16,
            "the widener's deletion must grow the record to 5..=20, or this fixture does \
             not reach the path it exists for"
        );

        let folded: u32 = widened
            .alleles
            .iter()
            .map(|allele| allele.support.num_obs)
            .sum();
        assert_eq!(
            folded, 2,
            "only `opener` and `widener` witnessed this record as one run; `wide` folded \
             into it before the widen and must have been subtracted back out when its \
             witness split in two. Record: {widened:?}"
        );
    }

    #[test]
    fn find_or_create_allele_returns_same_bucket_on_match() {
        let mut rec = OpenPileupRecord::new(0, 100, b"ACG".to_vec());
        let idx1 = find_or_create_allele_index(&mut rec.alleles, b"ACT".to_vec());
        rec.alleles[idx1].support.num_obs = 1;
        let idx2 = find_or_create_allele_index(&mut rec.alleles, b"ACT".to_vec());
        assert_eq!(idx1, idx2);
        assert_eq!(rec.alleles[idx2].support.num_obs, 1);
        // REF + ACT = 2 buckets total.
        assert_eq!(rec.alleles.len(), 2);
    }

    #[test]
    fn zero_event_bq_zeros_match_only_preserving_indel_proxies() {
        // Mate-overlap loser's BQ is zeroed across its haplotype
        // window. Per `pileup_walker.md` §"Mate overlap", the rule
        // applies to the Match event at the overlap position — the
        // loser's *indel* events (which can sit at other positions
        // in the window) carry independent evidence and must keep
        // their `bq_proxy`. Today the haplotype-level `min`
        // collapses to 0 via the zeroed Match regardless, so this is
        // invisible at output level. Pin the helper's contract here
        // so a future change to `ln_bq_for_read`'s reduction (e.g.
        // median or weighted) cannot silently corrupt indel BQ
        // proxies. Mi3 in `ia/reviews/pileup_2026-05-09.md`.
        let mut m = ReadEvent::Match {
            ref_pos: 100,
            base: b'A',
            bq_baq: 25,
        };
        zero_event_bq(&mut m);
        match m {
            ReadEvent::Match { bq_baq, .. } => {
                assert_eq!(bq_baq, 0, "Match BQ must be zeroed on the overlap loser");
            }
            _ => panic!("Match shape must survive zeroing"),
        }

        let mut ins = ReadEvent::Insertion {
            anchor_ref_pos: 130,
            seq: vec![b'A', b'C'],
            bq_proxy: 25,
        };
        zero_event_bq(&mut ins);
        match ins {
            ReadEvent::Insertion { bq_proxy, .. } => {
                assert_eq!(
                    bq_proxy, 25,
                    "Insertion bq_proxy must survive — the indel sits at a different anchor than the overlap",
                );
            }
            _ => panic!("Insertion shape must survive zeroing"),
        }

        let mut del = ReadEvent::Deletion {
            anchor_ref_pos: 130,
            deleted_len: 3,
            bq_proxy: 25,
        };
        zero_event_bq(&mut del);
        match del {
            ReadEvent::Deletion { bq_proxy, .. } => {
                assert_eq!(
                    bq_proxy, 25,
                    "Deletion bq_proxy must survive — the indel sits at a different anchor than the overlap",
                );
            }
            _ => panic!("Deletion shape must survive zeroing"),
        }
    }

    // --- M18: subtract_contribution direct unit tests -----------
    //
    // The function is the most subtle correctness primitive in the
    // module: every record widen re-folds via subtract-then-add.
    // Pin the saturate-to-zero contract on the four `u32` fields
    // and the straight `f64` subtract on `q_sum`.

    // Mi7: in debug builds an underflow trips the `debug_assert!`
    // peer; in release builds the `saturating_sub` zeroes the
    // field. We pin both ends of the contract.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "subtract_contribution underflow")]
    fn subtract_contribution_panics_on_underflow_in_debug() {
        let mut s = AlleleSupportStats {
            num_obs: 1,
            q_sum: 0.0,
            fwd: 1,
            placed_left: 0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        let c = AlleleSupportStats {
            num_obs: 5,
            q_sum: 0.0,
            fwd: 0,
            placed_left: 0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        subtract_contribution(&mut s, &c);
    }

    #[cfg(not(debug_assertions))]
    #[test]
    fn subtract_contribution_saturates_u32_fields_to_zero_in_release() {
        let mut s = AlleleSupportStats {
            num_obs: 1,
            q_sum: -1.0,
            fwd: 1,
            placed_left: 0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        let c = AlleleSupportStats {
            num_obs: 5,
            q_sum: -1.0,
            fwd: 5,
            placed_left: 5,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        subtract_contribution(&mut s, &c);
        assert_eq!(s.num_obs, 0);
        assert_eq!(s.fwd, 0);
        assert_eq!(s.placed_left, 0);
        // q_sum is straight f64 subtract by design (signed):
        // (-1.0) - (-1.0) = 0.0.
        assert!((s.q_sum - 0.0).abs() < 1e-12);
    }

    #[test]
    fn add_then_subtract_contribution_round_trips_for_u32_fields() {
        // The widen path subtracts the prior contribution and adds
        // the new one. When old == new, the bucket must end up
        // unchanged.
        let mut bucket = AlleleSupportStats {
            num_obs: 7,
            q_sum: -3.0,
            fwd: 4,
            placed_left: 2,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        let c = AlleleSupportStats {
            num_obs: 1,
            q_sum: -2.0,
            fwd: 1,
            placed_left: 0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };
        let snapshot = bucket;
        add_contribution(&mut bucket, &c);
        subtract_contribution(&mut bucket, &c);
        assert_eq!(bucket.num_obs, snapshot.num_obs);
        assert_eq!(bucket.fwd, snapshot.fwd);
        assert_eq!(bucket.placed_left, snapshot.placed_left);
        assert!((bucket.q_sum - snapshot.q_sum).abs() < 1e-12);
    }
}
