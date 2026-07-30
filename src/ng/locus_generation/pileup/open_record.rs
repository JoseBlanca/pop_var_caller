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
use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};
use crate::pileup_record::ChainId;

use super::super::{
    LocusKind, LocusLen, ReadWitness, SampleLocusObservations, SequenceObservation,
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
/// # `placed_start` is gone, and this is where it went
///
/// It was reconstructed at the `PileupRecord` boundary through Milestone A, from a per-read
/// flag on [`FoldedReadState`], so the stage-1 differential could keep comparing every field
/// of every record while A2–A5 changed the fold underneath it. **B2 removed that boundary**
/// and both went with it. `parity.rs` now zeroes the field on *both* sides and names the
/// removal where it does so, which is what keeps a deliberate absence distinguishable from
/// the oversight the Milestone A review found.
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

/// How the reads folded into one finished record witnessed it — the tally
/// [`OpenPileupRecord::finalise`] produces on its way past every
/// [`FoldedReadState`], resolving [`coverage_of`] once against the **final**
/// footprint.
///
/// **Two counts rather than the per-read runs themselves, and only until B2.**
/// `finalise` still returns production's [`PileupRecord`], which has nowhere to
/// put a [`ReadWitness`]; B2 replaces the return with `SampleLocusObservations`,
/// where each row carries its own. What has to be true *before* that is the
/// resolution point — coverage read at fold time is measured against a footprint
/// the record may still outgrow — so A4 resolves it here and reports what it
/// found, and B2 turns the same loop into rows.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct RecordWitness {
    /// Folded reads that witnessed every position of the final footprint.
    pub reads_complete: u32,
    /// Folded reads that witnessed one contiguous run short of it.
    pub reads_partially_observed: u32,
    /// Reads that covered this record and yielded **no observation at all** — their
    /// witnessed positions inside the footprint were non-contiguous, which one
    /// `Observed` run cannot describe honestly (spec §6). The size of
    /// [`OpenPileupRecord::reads_without_observation`], which is a set of read ids
    /// precisely so this number cannot be inflated by the footprint's length.
    pub reads_without_observation: u32,
    /// Reads a depth cap truncated at **every** position of this record's footprint where
    /// they had events, so they folded nowhere — see
    /// [`OpenPileupRecord::reads_discarded_by_cap`] for why that is the quantity rather
    /// than "how often this record's columns were truncated".
    pub reads_discarded_by_cap: u32,
}

/// How much of the finished record a read witnessed, in **locus positions**.
///
/// Resolved once, at [`finalise`](OpenPileupRecord::finalise), from the read's
/// `witnessed` extent against the record's **final** footprint — never during the
/// fold, when that footprint is not yet known. A read that was a complete witness
/// when it folded becomes `Observed` after a later widen **with nothing about the
/// read having changed**, and there is no re-fold that would notice: the read may
/// have expired long before the record closed, since `expire_passed` touches no
/// open record (spec §4, §6).
///
/// The extent is **clamped into the footprint at both ends** before it is measured. That is
/// not belt-and-braces: `events_overlapping` does not clip a deletion, so a deletion
/// anchored before the record comes back whole and its run can reach past the record's end
/// (spec §8). `apply_events_into` already intersects for the extent it stores, and this
/// repeats it so a wrong `offset_in_locus` cannot arrive by another route — the two numbers
/// are only ever wrong quietly.
///
/// The clamp was **one-sided** as first written: the left edge was pulled up to
/// `record_pos` but never pushed down to `record_end_exclusive`, so an extent lying entirely
/// right of the footprint yielded an `offset_in_locus` past the end of the locus paired with
/// `positions_covered: 0` — a run of no positions, which both [`RefSpan`] and
/// [`ReadWitness`] document as not existing. Unreachable from the fold today, because a
/// read that witnessed nothing inside a record does not fold into it; the `debug_assert`
/// says so, and the clamp means a future caller gets a truthful answer rather than that one.
///
/// # The `u16` narrowing is bounded by config, not by a constant
///
/// Runs are expressed in `u16`, and the narrowing goes through [`LocusLen`], which owns the
/// saturating cast. **An earlier version of this comment claimed saturation was unreachable
/// because "production's `max_record_span` is 5000". That is the wrong bound**: the cap is
/// `--max-record-span`, an unbounded `u32` CLI flag, so a caller could configure a footprint
/// wider than `u16::MAX` and a partial witness inside it would report a truncated
/// `positions_covered` with no error.
///
/// **Settled: `PileupGeneratorConfig` caps `max_record_span` at `u16::MAX` and rejects more
/// (owner, 2026-07-29).** The cap costs nothing real — a locus is at most ~100 bp, and a
/// 5,000 bp record is already unreachable with Illumina reads, so the existing default is
/// generous by fifty-fold and this ceiling by six hundred. Widening the run to `u32` would
/// touch the shared locus type and the STR generator that also mints coverage, to buy a range
/// no data can occupy. **This is the one knob where ng's constant is not simply production's**
/// — inheriting it "by name" would inherit the hazard.
///
/// So the `debug_assert` below is the invariant's statement and **C1 is its enforcement**;
/// until C1 lands, ng's walker is reachable only from tests and the default 5,000 leaves
/// thirteen-fold headroom.
pub(super) fn coverage_of(
    witnessed: RefSpan,
    record_pos: u32,
    record_end_exclusive: u32,
) -> ReadWitness {
    debug_assert!(
        record_end_exclusive.saturating_sub(record_pos) <= u32::from(u16::MAX),
        "record footprint {record_pos}..{record_end_exclusive} is wider than a `u16` run can \
         describe; `max_record_span` needs a ceiling (C1)",
    );
    debug_assert!(
        witnessed.start < record_end_exclusive && witnessed.end >= record_pos,
        "witnessed extent {witnessed:?} does not intersect the footprint \
         {record_pos}..{record_end_exclusive}; a read that witnessed nothing inside a record \
         does not fold into it",
    );
    let first = witnessed.start.clamp(record_pos, record_end_exclusive);
    let past_last = witnessed
        .end
        .saturating_add(1)
        .clamp(first, record_end_exclusive);
    if first <= record_pos && past_last >= record_end_exclusive {
        return ReadWitness::Complete;
    }
    ReadWitness::Observed {
        offset_in_locus: LocusLen::from_positions(u64::from(first - record_pos)).get(),
        positions_covered: LocusLen::from_positions(u64::from(past_last - first)).get(),
    }
}

/// **The identity of one emitted observation** — what makes two reads the same row.
///
/// Three parts, and the two ng adds are the point (spec §6): the bases a read showed;
/// **how much of the locus it witnessed**, because a complete witness and a partial one
/// of the same bases are different evidence; and **which read group it came from**,
/// because a per-chemistry model needs the allele × group cross *with its quality
/// moments*, which a per-group count beside one merged observation cannot give.
///
/// **Only the bases are decidable while the record is open.** Coverage is relative to a
/// footprint that grows until the record closes (A4), so the fold keys its buckets on
/// bases alone and the full identity is realised at `finalise` — which is where arch §1.2
/// puts it. That is why rows are re-derived *per read* rather than read off the per-bucket
/// totals: coverage and group are facts about a read, not about a bucket.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ObservationKey {
    pub bases: Vec<u8>,
    pub read_witness: ReadWitness,
    pub read_group: ReadGroupId,
}

/// One row of a finished record, accumulated across the reads that share its
/// [`ObservationKey`].
#[derive(Debug, Clone)]
pub(super) struct ObservationRow {
    pub key: ObservationKey,
    pub support: AlleleSupportStats,
    /// The chain ids of the reads in this row — **absent for a read that agreed with the
    /// reference across everything it witnessed** (spec §6).
    ///
    /// Production's rule is positional: `allele_index == 0`, the REF bucket. That named a
    /// unique row while there was one row per allele, and it no longer does — rows split by
    /// coverage and by read group, so a reference-matching read can sit in a *partial* row
    /// whose bases are a prefix of the reference bytes and never compare equal to them.
    ///
    /// So the rule is stated **per read** instead: it is decidable at fold time from what
    /// the read is, it survives every split of the rows, and it reduces to production's
    /// exactly when the rows are one-per-allele. A chain id marks which haplotype a read
    /// came from, and the reference is the default — a default needs no tag.
    pub chain_ids: Vec<ChainId>,
}

/// A total order over [`ReadWitness`], which is `Eq` but not `Ord` — the shared type has
/// no natural ordering, and inventing an `Ord` impl on it would export this file's sorting
/// convention to every other consumer.
/// `pub(super)` for the differential: D1's projection has to lay production's alleles out in
/// **ng's** emission order, and a second spelling of this comparator in `parity.rs` is a
/// spelling that can drift from the one the walk actually uses.
pub(super) fn coverage_order(coverage: ReadWitness) -> (u8, u16, u16) {
    match coverage {
        ReadWitness::Complete => (0, 0, 0),
        ReadWitness::Observed {
            offset_in_locus,
            positions_covered,
        } => (1, offset_in_locus, positions_covered),
    }
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
    /// The reads that covered this record and produced no observation: their
    /// witnessed positions inside the footprint were non-contiguous — an interior
    /// `N`, a ref-skip — and `Observed` describes one run, so there is nothing
    /// honest to emit (spec §6).
    ///
    /// **A set of read ids, and that is the whole point of the field (spec §4).**
    /// The no-observation path is reached at *every* position the record is affected
    /// at, so a counter incremented there would multiply by the footprint length —
    /// the same once-per-record-not-once-per-position bug the subtract-then-add
    /// mechanism exists to prevent, on the one path with no inherited test to catch
    /// it. Membership is idempotent; the count is taken at `finalise`.
    ///
    /// A `Vec` rather than a hash set because the population is tiny by construction
    /// — adaptor masking and the dropped-indel rule always truncate from one side,
    /// so they stay expressible as a run, and only an interior hole lands here — and
    /// an empty `Vec` costs no allocation on a path that runs once per covered base.
    ///
    /// *Invariant: membership is monotone.* The window's left edge is the record's
    /// anchor and never moves; widening only extends the right. A hole inside the
    /// footprint therefore stays a hole, so a read recorded here can never fold
    /// successfully later — asserted in `fold_read_into_record` rather than handled,
    /// because a removal path would be code no input can reach.
    reads_without_observation: Vec<u32>,
    /// The reads a depth cap truncated at some position of this record's footprint —
    /// **candidates for `reads_discarded_by_cap`, not the count itself** (spec §6).
    ///
    /// # Why the obvious per-record count is wrong, and this is the correction
    ///
    /// Production counts `column_depth_truncations` on `RunSummary`: *positions* truncated,
    /// run-wide. A read can be truncated at one position of a footprint and survive at
    /// another, and **if it folds at all it folds with its whole window**, so its evidence
    /// is not subsampled. Counting truncation events per record would therefore flag
    /// records whose support is complete.
    ///
    /// What the locus type wants is "*the support counts are a subsample, not the depth*":
    /// **reads that had events inside this footprint and were truncated at every position
    /// where they did, so folded nowhere.** That is why this is a membership list resolved
    /// at `finalise` against `folded_reads` rather than a counter — a read here that folded
    /// later is not discarded, it is present.
    ///
    /// **The cap truncates in the walk, before any record is identified**, so these ids are
    /// plumbed in from `genome_walk` rather than discovered here. Two cases have no clean
    /// answer and are recorded rather than solved: a read truncated where no record is open
    /// is unattributable to any locus, and a truncated read carrying a *deletion* would have
    /// widened a record — so dropping it changes the footprint, and with it every other
    /// read's coverage.
    reads_discarded_by_cap: Vec<u32>,
}

/// Record `read_id` as having yielded no observation in this record, once.
///
/// Linear membership over a `Vec`: see [`OpenPileupRecord::reads_without_observation`]
/// for why the population is small enough that this beats a set that allocates.
fn note_no_observation(reads_without_observation: &mut Vec<u32>, read_id: u32) {
    if !reads_without_observation.contains(&read_id) {
        reads_without_observation.push(read_id);
    }
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
    /// It carried an `#[expect(dead_code)]` until the widen re-place began rebuilding the
    /// whole state from an exhaustive destructure, which reads it — so the backstop
    /// **fired**, as designed, and was removed rather than downgraded to an `allow`. Being
    /// read is not the same as being *used*: **B1 still owes making it part of the bucket
    /// key**, and until then it is carried through the fold without changing any answer.
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
    /// **A3 keeps it current across a widen** (owner, 2026-07-29): `refold_live_reads`
    /// re-places every live folded read and rewrites this field, so a read sitting inside
    /// its own deletion at the widening position no longer carries an extent measured
    /// against a footprint the record has outgrown. A4 resolves `ReadWitness` from it,
    /// and would otherwise have reported a wrong depth with no error.
    witnessed: RefSpan,
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
            // Unallocated until the first non-contiguous witness, which most records
            // never see — see the field's own note. Same for the cap list: a truncated
            // column is rare and most runs never see one at all.
            reads_without_observation: Vec::new(),
            reads_discarded_by_cap: Vec::new(),
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

    /// Whether this read **agreed with the reference across everything it witnessed** — the
    /// per-read form of production's `allele_index == 0` chain-id rule (spec §6).
    ///
    /// Production's rule is positional and unportable once rows split: a reference-matching
    /// read that witnessed only part of the footprint sits in a *partial* row whose bases
    /// are a prefix of the reference bytes, so "the row whose bases equal
    /// `reference_bases`" no longer names it. Asking the question of the **read** instead
    /// works at every split, and reduces to production's rule exactly when the rows are
    /// one-per-allele: a complete witness that agreed everywhere *is* the REF bucket.
    ///
    /// The comparison is against the reference **over the positions the read witnessed**,
    /// not the whole footprint — that is the entire difference. Note the two axes do not
    /// have to line up: `bases` is in read coordinates and an insertion makes it longer
    /// than the run it covers, which simply means it cannot equal the reference slice, and
    /// the read keeps its id. Correct, and for the right reason.
    fn read_agreed_with_reference(&self, state: &FoldedReadState) -> bool {
        let reference = &self.alleles[0].seq;
        let first = state.witnessed.start.saturating_sub(self.pos) as usize;
        let past_last = (state.witnessed.end.saturating_sub(self.pos) as usize)
            .saturating_add(1)
            .min(reference.len());
        if first >= past_last {
            return false;
        }
        self.alleles[state.allele_index].seq == reference[first..past_last]
    }

    /// Re-derive this record's rows **per read**, keyed on the full [`ObservationKey`].
    ///
    /// # Why per read, when the buckets already hold the totals
    ///
    /// Two of the three parts of a row's identity are facts about a *read*, not a bucket. A
    /// bucket knows its bases; it does not know that one of its reads witnessed the whole
    /// footprint while another saw one position of fourteen, nor that they came from
    /// different lanes. Reading rows off the bucket totals can only ever produce the merged
    /// answer, which spec §6 says is not good enough.
    ///
    /// With one read group and every read a complete witness the result **is** the bucket
    /// totals, row for row and sum for sum, because every read in a bucket then shares one
    /// key. That is what "free at one read group" means, and it is what keeps the stage-1
    /// differential green across B1.
    ///
    /// # Reads are taken in `read_id` order, and **that line is the determinism guarantee**
    ///
    /// `folded_reads` is an `AHashMap` whose iteration order is seeded per process, and
    /// `q_sum` is an `f64` sum, so accumulating in hash order would make the last bits of
    /// every quality total depend on the seed — the same input emitting different bytes run
    /// to run, which spec §7 forbids. The bucket totals this replaces had the property for
    /// free, from being accumulated in walk order; re-deriving is what puts it at risk, and
    /// `ids.sort_unstable()` below is what restores it.
    ///
    /// Pinned by `parity::ng_emits_the_same_bytes_in_a_second_process`, which walks the same
    /// input in two processes: delete that sort and it fails; delete `finalise`'s *row* sort
    /// and it stays green. No test inside one process can see any of this, because `ahash`
    /// seeds once per process.
    fn observation_rows(&self, record_end_exclusive: u32) -> Vec<ObservationRow> {
        let mut ids: Vec<u32> = self.folded_reads.keys().copied().collect();
        ids.sort_unstable();

        let mut rows: Vec<ObservationRow> = Vec::new();
        for read_id in ids {
            // PANIC-FREE: the id came from this map's own keys a few lines above, and
            // nothing between then and here mutates the map.
            let state = self
                .folded_reads
                .get(&read_id)
                .expect("the id came from this record's own fold state");
            // The identity is compared **borrowed** and the bases cloned only when a row is
            // genuinely new. Cloning into the key first and letting `find` compare owned
            // values costs one allocation *per read* where the rows need one *per row* — the
            // only cost in this function that scales with depth, and measured at 2.2 % of
            // the milestone's 15.1 %. (The `Vec<u32>` and its sort are another 2.2 % and
            // stay: they are the determinism guarantee. The linear `find` itself measured
            // 0 %, and hash-keying the rows measured *worse*.)
            let bases = self.alleles[state.allele_index].seq.as_slice();
            let read_witness = coverage_of(state.witnessed, self.pos, record_end_exclusive);
            let read_group = state.read_group;
            let agreed_with_reference = self.read_agreed_with_reference(state);
            let existing = rows.iter().position(|row| {
                row.key.bases == bases
                    && row.key.read_witness == read_witness
                    && row.key.read_group == read_group
            });
            let row = match existing {
                Some(index) => &mut rows[index],
                None => {
                    rows.push(ObservationRow {
                        key: ObservationKey {
                            bases: bases.to_vec(),
                            read_witness,
                            read_group,
                        },
                        support: AlleleSupportStats::default(),
                        chain_ids: Vec::new(),
                    });
                    // PANIC-FREE: pushed on the line above.
                    rows.last_mut().expect("just pushed")
                }
            };
            add_contribution(&mut row.support, &state.contribution);
            if !agreed_with_reference {
                row.chain_ids.push(state.chain_id);
            }
        }
        for row in &mut rows {
            row.chain_ids.sort_unstable();
            row.chain_ids.dedup();
        }
        rows
    }

    /// Convert into the finished locus ng emits: the region, the reference bytes under it,
    /// one [`SequenceObservation`] per row, and the two per-record counters.
    ///
    /// **Coverage is resolved here, and here is the only place it can be.** A read's
    /// `witnessed` extent is absolute; what it *means* — complete witness, or one run short
    /// of it — is relative to a footprint that grows until the record closes. Resolving at
    /// fold time would answer against a footprint the record may still outgrow, and no
    /// re-fold would come back to correct it, because the read may have expired in between
    /// (spec §4).
    ///
    /// Chain ids are emitted per row and **absent for a read that agreed with the reference
    /// across everything it witnessed** — see
    /// [`read_agreed_with_reference`](Self::read_agreed_with_reference) for why production's
    /// positional rule does not survive rows that split.
    ///
    /// `placed_start` is **gone**, as A1 said it would be: ng's stats never carried it, the
    /// `PileupRecord` boundary did, and this step removed that boundary. Nothing consumes
    /// it, and it is a pure function of the read's start against the anchor, so a later
    /// consumer re-derives it without changing the fold (spec §6).
    pub fn finalise(self) -> (SampleLocusObservations, RecordWitness) {
        let record_pos = self.pos;
        let record_end_exclusive = self.footprint_end_exclusive();
        let chrom_id = self.chrom_id;
        // Every field named, no `..Default::default()`: a field added to `RecordWitness`
        // would otherwise compile here and arrive as a silent `0`.
        let mut witness = RecordWitness {
            reads_complete: 0,
            reads_partially_observed: 0,
            reads_without_observation: self.reads_without_observation.len() as u32,
            // **Resolved here, not counted in the walk.** A read truncated at one position
            // of this footprint may have folded at another, and a read that folds does so
            // with its whole window — so only the ones still absent from `folded_reads` at
            // the end were actually discarded (spec §6).
            //
            // **And absent for the cap's reason, not for A5's.** "Not in `folded_reads`"
            // has two causes: the cap kept the read out, or it folded and then lost its row
            // when its witness turned out non-contiguous. Counting the second here reports
            // one read in *both* `reads_without_observation` and `reads_discarded_by_cap`,
            // which double-counts it and tells a model the support is a subsample when the
            // truth is that a read said nothing usable. Measured at 240 records in ~506,000
            // before this exclusion.
            reads_discarded_by_cap: self
                .reads_discarded_by_cap
                .iter()
                .filter(|read_id| {
                    !self.folded_reads.contains_key(read_id)
                        && !self.reads_without_observation.contains(read_id)
                })
                .count() as u32,
        };
        // **A3's eviction, checked where the buckets still exist — and D1 is why it is
        // here.** The property is "no bucket survives that no read is folded into", and it
        // is *invisible in the emitted locus*: `observation_rows` derives rows from
        // `folded_reads`, so a stranded bucket produces no row and leaves no trace. A parity
        // test asserted it on the emitted records, which was still meaningful while the walk
        // emitted `PileupRecord`s and stopped being so at B2; D1 mutated the code the test
        // named — moving `evict_unsupported_alleles` above the contributor fold loop, which
        // strands every bucket that loop empties — and **the whole 198-test module stayed
        // green.** That is the eleventh test on this branch that could not fail.
        //
        // A `debug_assert` rather than a test, because the invariant lives on a structure no
        // test outside this file can reach, and this way every walk in the suite checks it —
        // the census alone runs it over ~257,000 loci, and the `soak` profile keeps it armed
        // at release speed. The mutation above now fails there instead of nowhere.
        debug_assert!(
            self.alleles
                .iter()
                .skip(1)
                .all(|allele| allele.support.num_obs > 0),
            "a bucket no read is folded into survived to finalise: the eviction did not run \
             after the fold loop that emptied it. {:?}",
            self.alleles,
        );
        // **The rows come from the reads, not from the bucket totals** (B1). See
        // `observation_rows` for why, and for what makes the two agree at one read group.
        let mut rows = self.observation_rows(record_end_exclusive);
        for state in self.folded_reads.values() {
            match coverage_of(state.witnessed, record_pos, record_end_exclusive) {
                ReadWitness::Complete => witness.reads_complete += 1,
                ReadWitness::Observed { .. } => witness.reads_partially_observed += 1,
            }
        }
        // Every folded read is resolved exactly once and lands in exactly one of the two
        // classes; the no-observation set is disjoint from `folded_reads` by construction,
        // its members having been removed when they produced nothing.
        debug_assert_eq!(
            witness.reads_complete + witness.reads_partially_observed,
            self.folded_reads.len() as u32,
            "every folded read resolves to exactly one coverage class",
        );

        // **Sorted into a canonical order — and this is *not* what makes the output
        // deterministic, which is worth stating because the first version of this comment
        // said it was.** Mutation-tested both ways: deleting this sort leaves
        // `ng_emits_the_same_bytes_in_a_second_process` green, while deleting
        // `observation_rows`' `ids.sort_unstable()` fails it. Determinism is already won
        // upstream, by taking the reads in `read_id` order; first-seen row order inherits it.
        //
        // What the sort buys is that row order is a function of the row's **own identity**
        // rather than of which read happened to arrive first — so the emitted order does not
        // change when an unrelated read is added, and two loci with the same evidence
        // present it the same way. That is what a consumer diffing output wants, and it is
        // what the STR generator sorts for. Belt and braces on determinism; the belt is
        // upstream.
        rows.sort_by(|a, b| {
            a.key
                .bases
                .cmp(&b.key.bases)
                .then_with(|| {
                    coverage_order(a.key.read_witness).cmp(&coverage_order(b.key.read_witness))
                })
                .then_with(|| a.key.read_group.0.cmp(&b.key.read_group.0))
        });

        let observations = rows
            .into_iter()
            .map(|row| {
                // Exhaustively destructured on the way out, in the direction that can lose
                // information: a field added to ng's stats stops this compiling instead of
                // being silently dropped at the boundary.
                let AlleleSupportStats {
                    num_obs,
                    q_sum,
                    fwd,
                    placed_left,
                    mapq_sum,
                    mapq_sum_sq,
                } = row.support;
                SequenceObservation {
                    bases: row.key.bases.into_boxed_slice(),
                    read_witness: row.key.read_witness,
                    read_group: row.key.read_group,
                    num_obs,
                    num_fwd: fwd,
                    q_sum,
                    mapq_sum,
                    mapq_sum_sq,
                    placed_left,
                    chain_ids: row.chain_ids,
                }
            })
            .collect();

        // `alleles[0]` is the REF bucket and its bytes are the record's reference sequence
        // — there is no separate copy, which is why this moves rather than clones.
        let reference_bases = self
            .alleles
            .into_iter()
            .next()
            .expect("alleles[0] is created with the record and never evicted")
            .seq
            .into_boxed_slice();

        let locus = SampleLocusObservations {
            // `GenomeRegion` is 1-based **inclusive**, so the end is the last covered
            // position rather than one past it (spec §6).
            region: GenomeRegion {
                contig: ContigId(chrom_id),
                start: Position(u64::from(record_pos)),
                end: Position(u64::from(record_end_exclusive.saturating_sub(1))),
            },
            reference_bases,
            observations,
            reads_without_observation: witness.reads_without_observation,
            reads_discarded_by_cap: witness.reads_discarded_by_cap,
            kind: LocusKind::Generic,
        };
        (locus, witness)
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
    /// haplotype string built by [`apply_events_into`]. Hoisted
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
    /// Reusable scratch for the read ids `widen` re-folds. Held rather than allocated
    /// per widen, and **sorted** before use: `folded_reads` is an `AHashMap` with a
    /// per-process seed, so folding in its iteration order would make bucket *creation*
    /// order — and therefore the emitted allele order — differ run to run.
    refold_ids_buf: Vec<u32>,
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
            refold_ids_buf: Vec::new(),
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
        self.refold_ids_buf.clear();
    }

    /// The lowest anchor position still open, or `None` when nothing is —
    /// **ng's addition, C2 (plan 3)**, and the whole open-record half of the
    /// region walk's stop rule.
    ///
    /// The table is keyed by anchor, so "is any record anchored at or before
    /// the region's end still open?" is one `BTreeMap` first-key lookup rather
    /// than a scan. It has to be the *anchor* and not the footprint: a record
    /// anchored inside the region is the region's to emit however far past the
    /// boundary its footprint runs (spec §2).
    pub(super) fn first_open_anchor(&self) -> Option<u32> {
        self.records.keys().next().copied()
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

    /// Widen the record at `key` so its REF span covers up to `new_end_exclusive`, fetching
    /// the additional reference bases.
    ///
    /// **Only `alleles[0]` grows** — the record's own reference sequence. Every other bucket
    /// holds what some read witnessed, and a read's witness does not change because the
    /// window around it did (A3, spec §4). Production appends the new bases to *every*
    /// bucket, and this doc comment said so until A3 inverted the behaviour and left the
    /// sentence standing; the two paragraphs below the fetch carry the full argument.
    ///
    /// **Every live folded read is then re-placed against the wider window**
    /// (`refold_live_reads`), which is why this takes `active_reads` and `contributors`.
    /// Production re-folds only the contributors at the current walker position, so a read
    /// sitting inside its own deletion — live, in this record, no event anchored here —
    /// would otherwise keep an extent measured against a footprint the record has outgrown,
    /// and A4 resolves coverage from exactly that.
    ///
    /// Returns `true` when the record actually widened; `false` when `new_end_exclusive` was
    /// already covered (no-op). Callers use the bool to count real widen events without
    /// conflating them with fresh `open_new` calls.
    fn widen(
        &mut self,
        key: u32,
        new_end_exclusive: u32,
        reference: &dyn RefSeq,
        active_reads: &ActiveReads,
        contributors: &[ReadContribution],
    ) -> Result<bool, WalkerError> {
        // Split-borrow so the fetch can write into `widen_bases_buf` while
        // `rec` holds a mutable borrow of `records`. `max_record_span` is
        // read through the same destructure rather than through `self`,
        // which the outstanding `records` borrow would otherwise block.
        let Self {
            records,
            allele_seq_buf,
            widen_bases_buf,
            refold_ids_buf,
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

        // **Only the REF bucket grows.** `alleles[0]` is the record's own reference
        // sequence and genuinely does get longer; every other bucket holds what some
        // read witnessed, and a read's witness does not change because the window
        // around it did.
        //
        // Production appends `extra_bases` to *every* bucket, with a 25-line comment
        // above the loop proving that this reproduces what a re-fold would emit
        // ([open_record.rs:390-415](../../../../src/pileup/walker/open_record.rs#L390)).
        // That is exactly the mechanism the no-fabrication rule has to remove: it is
        // what lets an **expired** read's bucket keep growing with reference bases the
        // read never saw, retroactively, after it has left the active set. A live read
        // re-folds against the wider window at the next position it has an event
        // inside the footprint and lands wherever its bases put it; an expired one
        // keeps a bucket whose bases already say exactly what it saw (spec §4).
        //
        // **This is what makes the rule implementable at all** — production cannot
        // express it, because the bases live on the shared bucket and its
        // `FoldedReadState` holds none of its own.
        rec.alleles[0].seq.extend_from_slice(extra_bases);

        // **Every live read that had folded into this record folds again, against the
        // wider window** (owner, 2026-07-29). Production's `process_position` re-folds
        // only the *contributors* at the current walker position, and a read sitting
        // inside its own deletion has no event anchored there — so it is live, it is in
        // this record, and it does not re-fold. Production hides that by appending the
        // reference bases to every bucket; with REF-only widening it would leave the
        // read's `witnessed` extent pinned to the pre-widen footprint, and `finalise`
        // resolves coverage from exactly that. The result would be a **wrong depth**, at
        // the long-deletion loci this port exists to fix, with no error.
        //
        // Spec §4 already asserts "a live read re-folds against the wider window"; this
        // is what makes the assertion true rather than nearly true.
        //
        // **An expired read is deliberately not re-folded** — it cannot be, since its
        // cursor is gone with it, and it should not be: its bucket already says exactly
        // what it saw, and extending it is the retroactive fabrication this rule removes.
        refold_live_reads(
            rec,
            allele_seq_buf,
            refold_ids_buf,
            active_reads,
            contributors,
        );

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
/// 4. Each event's **reach** — `anchor_pos() + footprint_span()`,
///    the first position past what it witnesses — is
///    non-decreasing. Preconditions 2 and 3 do *not* imply this:
///    two `Deletion`s anchored at *a* with `deleted_len` 5 then 2
///    satisfy both and reach *a+6* then *a+3*. The cursor's CIGAR
///    walk cannot emit that — a read deletes a run once, and the
///    op after a D op starts at the far end of it — so this is a
///    demand on new callers rather than a property of the current
///    ones, and it is what makes the run below grow monotonically.
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

    // Precondition 4, asserted separately from the sort because sorting does not
    // imply it: a longer deletion behind a shorter one at the same anchor is
    // sorted and reaches backwards. This is the assertion the `run_end.max(…)`
    // below defends against in release builds.
    debug_assert!(
        events.windows(2).all(|w| {
            let a = w[0].anchor_pos().saturating_add(w[0].footprint_span());
            let b = w[1].anchor_pos().saturating_add(w[1].footprint_span());
            a <= b
        }),
        "apply_events_into: each event's reach (anchor + footprint_span) must be non-decreasing; got {events:?}",
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
                    // `.max()`, not a plain `event_end`: under precondition 4 the
                    // two are equal — which is why no test can tell them apart,
                    // and why deleting the `.max()` fails nothing. It is the
                    // release build's defence for the input that precondition
                    // rules out (two deletions at one anchor, the longer first),
                    // where a plain assignment would shrink the run below what
                    // the read witnessed. It defends the **stated** precondition,
                    // not the cursor's output: the cursor cannot produce that
                    // input, so the debug assertion above is where a caller that
                    // could learns so.
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
/// [`apply_events_into`]: at the same anchor, Match must come
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

/// The three record fields a fold mutates, split-borrowed out of an
/// [`OpenPileupRecord`] so the caller can hold `allele_seq_buf` — which lives on the
/// *table*, not the record — mutably at the same time.
///
/// Grouped rather than passed as three parameters because they always travel together
/// and always come from the same destructure: a fold that could touch the buckets
/// without being able to record a read that produced no observation would be the A5
/// defect, spelled in the type.
struct RecordFoldState<'a> {
    alleles: &'a mut Vec<OpenAllele>,
    folded_reads: &'a mut AHashMap<u32, FoldedReadState>,
    reads_without_observation: &'a mut Vec<u32>,
    reads_discarded_by_cap: &'a mut Vec<u32>,
}

/// Fold one read into one open record, exactly once: build the haplotype its events
/// present, place it in the matching bucket, and record what it contributed.
///
/// **Idempotent by construction** — the read's prior contribution is subtracted from
/// whichever bucket it was in before the new one is added, which is what makes "each
/// (record, read) pair folds exactly once over the record's lifetime" hold across
/// re-folds. That invariant is the single thing most likely to be lost in a port: a
/// six-base footprint would otherwise count a spanning read six times.
///
/// Its sibling [`refold_live_reads`] deliberately does **not** call this: a widen re-places
/// a read without re-deciding its quality, which this function would recompute.
///
/// `bq_fallback` is `ln_bq_for_read`'s answer for an **empty** window, which no caller can
/// reach: a read with no events in the window yields no observation a line below. The
/// contributor path passes the contributor's walker-position BQ anyway, because that is
/// what production passes and the parity claim is about the walk.
fn fold_read_into_record(
    record: &mut RecordFoldState<'_>,
    allele_seq_buf: &mut Vec<u8>,
    rec_pos: u32,
    active: &super::active_read_set::ActiveRead,
    window_events: &[ReadEvent],
    bq_fallback: u8,
) {
    let RecordFoldState {
        alleles,
        folded_reads,
        reads_without_observation,
        // The fold does not decide what the cap removed — that is registered by
        // `process_position` after this loop, from ids the walk hands down.
        reads_discarded_by_cap: _,
    } = record;
    let Some(witnessed) =
        apply_events_into(allele_seq_buf, rec_pos, &alleles[0].seq, window_events)
    else {
        // A5: **which** reads these were, as a per-record set. The path is reached at
        // every position the record is affected at, so a counter here would multiply
        // by the footprint length (spec §4).
        note_no_observation(reads_without_observation, active.read_id);
        // A non-contiguous witness — an interior `N`, or a ref-skip — yields no
        // observation (spec §6). **The prior contribution has to come off first**, and
        // this is not a corner: a read that folded contiguously *becomes* non-contiguous
        // when the window widens right across an interior gap, so a bare early return
        // would strand a live contribution in a bucket for a read that now has no row,
        // breaking `chain_ids.len() <= num_obs` silently and only on multi-base records
        // (spec §4).
        if let Some(prev) = folded_reads.remove(&active.read_id) {
            subtract_contribution(&mut alleles[prev.allele_index].support, &prev.contribution);
        }
        return;
    };

    // PANIC-FREE, and the assert states the invariant rather than defending against it:
    // the record's left edge is its anchor and never moves, and a widen only extends the
    // right, so a hole inside the footprint stays a hole. A read that once yielded no
    // observation cannot fold successfully later — which is why nothing removes it from
    // the set, and why a removal path would be code no input can reach.
    debug_assert!(
        !reads_without_observation.contains(&active.read_id),
        "read {} folded after having been recorded as witnessing nothing contiguous",
        active.read_id,
    );

    let ln_q = ln_bq_for_read(window_events, bq_fallback).max(active.read.mq_log_err);
    let mapq = u32::from(active.read.mapq);
    let new_contribution = AlleleSupportStats {
        num_obs: 1,
        q_sum: ln_q,
        fwd: u32::from(!active.read.is_reverse_strand),
        placed_left: u32::from(active.read.alignment_start < rec_pos),
        mapq_sum: mapq,
        mapq_sum_sq: (mapq as u64) * (mapq as u64),
    };

    if let Some(prev) = folded_reads.remove(&active.read_id) {
        subtract_contribution(&mut alleles[prev.allele_index].support, &prev.contribution);
    }
    // Borrowed lookup; only `clone()` the bytes when adding a genuinely new allele. In
    // SNP/REF steady state every contributor lands in the same existing bucket, so the
    // clone path never fires.
    let new_index = match find_allele_index(alleles, allele_seq_buf) {
        Some(index) => index,
        None => {
            alleles.push(OpenAllele::new(allele_seq_buf.clone()));
            alleles.len() - 1
        }
    };
    add_contribution(&mut alleles[new_index].support, &new_contribution);
    folded_reads.insert(
        active.read_id,
        FoldedReadState {
            allele_index: new_index,
            contribution: new_contribution,
            chain_id: active.chain_id,
            read_group: active.read.read_group,
            witnessed,
        },
    );
}

/// Re-place every **live** read already folded into this record against the window it now
/// has. Called from `widen`, immediately after `alleles[0]` grows — see the comment there
/// for why the contributor-only re-fold production does is not enough.
///
/// # It moves the read; it does not re-decide anything
///
/// The read's **`contribution` is carried across untouched** and only its bases, its
/// bucket and its `witnessed` extent are recomputed. That is not an optimisation, it is
/// the correctness argument: a widen changes *which positions the record covers*, not
/// *what quality evidence a read carried*. `q_sum` in particular encodes decisions the
/// **walk** took at earlier positions and replayed into the fold — a mate-overlap loser's
/// zeroed Match, an agree-case keeper's summed BQ — and those live on the `ReadContribution`
/// of the position they were taken at, which is long gone by the time a later widen fires.
/// Recomputing quality here would silently undo every reconciliation made earlier in the
/// record's footprint, which is exactly the failure spec §5 names: *"forget the replay and
/// reconciliation silently applies at one position out of a record's whole footprint"*.
///
/// A **contributor at this position** is skipped outright: the fold loop is about to
/// re-fold it *with* its mate-overlap replay, so re-placing it here is work whose result is
/// overwritten a moment later.
///
/// **Unpinned, and deliberately so — the honest state of this line.** Mutating the skip away
/// leaves the whole suite green, and the reason is the carry above: a contributor re-placed
/// here keeps its old contribution, and the fold loop then subtracts that and adds the
/// recomputed one, landing on the same final state. The residue is bucket *creation* order,
/// which decides the order alleles are emitted in, and buckets left transiently empty, which
/// `evict_unsupported_alleles` clears at the end of the same call. So the skip earns its
/// place on cost, not on correctness — it was load-bearing before the contribution was
/// carried, and it is kept because re-placing every contributor twice per widen is pure
/// waste. Do not read the absence of a failing test here as the absence of a reason.
///
/// An **expired** read is skipped because it cannot be reached — its cursor went with it —
/// and should not be: its bucket already says exactly what it saw, and extending it is the
/// retroactive fabrication this whole rule removes.
///
/// Reads are taken in **`read_id` order**, not `folded_reads` iteration order: that map is
/// an `AHashMap` with a per-process seed, and bucket *creation* order decides the emitted
/// allele order, so working in hash order would make the output differ run to run.
fn refold_live_reads(
    rec: &mut OpenPileupRecord,
    allele_seq_buf: &mut Vec<u8>,
    ids: &mut Vec<u32>,
    active_reads: &ActiveReads,
    contributors: &[ReadContribution],
) {
    ids.clear();
    ids.extend(rec.folded_reads.keys().copied());
    ids.sort_unstable();

    let rec_pos = rec.pos;
    let rec_end = rec.footprint_end_exclusive();
    let OpenPileupRecord {
        alleles,
        folded_reads,
        reads_without_observation,
        ..
    } = rec;
    for &read_id in ids.iter() {
        if contributors.iter().any(|c| c.read_id == read_id) {
            continue;
        }
        let Some(active) = active_reads.get_by_read_id(read_id) else {
            continue;
        };
        // PANIC-FREE: `read_id` was just taken from this map's keys, and nothing between
        // then and here removes an entry other than this loop — which `continue`s past
        // any id it has already handled, each id appearing once.
        let previous = *folded_reads
            .get(&read_id)
            .expect("the id came from this record's own fold state");
        let window = active
            .cursor
            .events_overlapping(rec_pos, rec_end, &active.read);

        let Some(witnessed) = apply_events_into(allele_seq_buf, rec_pos, &alleles[0].seq, &window)
        else {
            // The wider window opened a hole in what was one run — see
            // `fold_read_into_record` for why the contribution has to come off, and
            // why the read is recorded rather than merely dropped. **This is the
            // arrival that makes the case ordinary rather than a corner:** a read
            // folds contiguously, the record widens right across an interior gap, and
            // the read that had a row a moment ago now has none.
            note_no_observation(reads_without_observation, read_id);
            folded_reads.remove(&read_id);
            subtract_contribution(
                &mut alleles[previous.allele_index].support,
                &previous.contribution,
            );
            continue;
        };

        let new_index = match find_allele_index(alleles, allele_seq_buf) {
            Some(index) => index,
            None => {
                alleles.push(OpenAllele::new(allele_seq_buf.clone()));
                alleles.len() - 1
            }
        };
        // Guarded, and not only to save work: subtracting and re-adding the same `f64`
        // into the same bucket is not the identity, and would put rounding noise into
        // `q_sum` at every widen for every read that did not move.
        if new_index != previous.allele_index {
            subtract_contribution(
                &mut alleles[previous.allele_index].support,
                &previous.contribution,
            );
            add_contribution(&mut alleles[new_index].support, &previous.contribution);
        }
        // **Rebuilt from an exhaustive destructure, not patched by assigning two fields.**
        // A re-place has to decide *every* field of the state it leaves behind. Assigning
        // `allele_index` and `witnessed` and letting the rest ride means a field added to
        // `FoldedReadState` errors at the fold literal in `fold_read_into_record` and is
        // silently carried stale here — which is precisely how `witnessed` itself went
        // wrong before A3, and the field's own doc records the consequence as "a wrong
        // depth with no error". The destructure makes the compiler ask the question at both
        // sites.
        let FoldedReadState {
            allele_index: _,
            contribution,
            chain_id,
            read_group,
            witnessed: _,
        } = previous;
        let state = folded_reads
            .get_mut(&read_id)
            .expect("the id came from this record's own fold state");
        *state = FoldedReadState {
            // The two the re-place decides.
            allele_index: new_index,
            witnessed,
            // The rest are facts about the read, which a widen does not change — see the
            // "it moves the read; it does not re-decide anything" note above.
            contribution,
            chain_id,
            read_group,
        };
    }
}

/// Drop every non-REF bucket no read is folded into, remapping the fold state.
///
/// **Required by A3's REF-only widening, not merely tidy.** Production's
/// append-to-every-bucket kept a re-folding read landing back in its *existing* bucket;
/// without it, a read that re-folds after a widen lands somewhere new and leaves the old
/// bucket behind at `num_obs == 0`. `find_allele_index` is a **linear scan with a full
/// byte compare**, run once per (record, contributor) at every position of the footprint,
/// and the comment above it — "records typically carry ≤ a few alleles" — stops being true
/// exactly at the long-deletion loci this port exists to fix (spec §7).
///
/// Called at the end of the fold rather than inside `widen`, because the empties are
/// created *by* the re-fold that a widen triggers, not by the widen itself.
///
/// `alleles[0]` is never evicted: it is the REF sequence, it is what `ref_span()`
/// measures, and production creates it with zero observations by design. The
/// counter-pressure spec §7 names is paid right here — a positional `allele_index` on
/// `FoldedReadState` has to be remapped, which is why a mapping is built rather than the
/// buckets simply retained.
fn evict_unsupported_alleles(
    alleles: &mut Vec<OpenAllele>,
    folded_reads: &mut AHashMap<u32, FoldedReadState>,
) {
    if alleles.len() < 2 || alleles[1..].iter().all(|a| a.support.num_obs > 0) {
        return;
    }
    // `None` is "evicted", rather than a `usize::MAX` sentinel. The sentinel makes the
    // evicted state *representable* as an index, so a missed remap is a plausible-looking
    // subscript rather than a type error — and this project has hit that trap in release
    // twice already ([locus_generation/mod.rs:79](super::super)). With `Option` the only way
    // past it is an `expect` that says what went wrong.
    let mut mapping: Vec<Option<usize>> = Vec::with_capacity(alleles.len());
    let mut index = 0usize;
    alleles.retain(|allele| {
        let keep = index == 0 || allele.support.num_obs > 0;
        index += 1;
        // `kept` was a third counter tracking what `mapping` already knows.
        mapping.push(keep.then(|| mapping.iter().filter(|slot| slot.is_some()).count()));
        keep
    });
    for state in folded_reads.values_mut() {
        // PANIC-FREE, and the message says why: every bucket a folded read points at holds
        // that read's own `num_obs`, so it is never one of the evicted ones. Reaching the
        // `expect` would mean the fold state and the bucket totals had drifted apart.
        state.allele_index = mapping[state.allele_index].expect(
            "a folded read pointed at a bucket with no support — the fold state and the \
             bucket totals have drifted",
        );
    }
}

/// Find or create the allele bucket inside a record matching `seq`.
/// Returns the bucket index. Linear scan is fine — records
/// typically carry ≤ a few alleles.
///
/// Takes `&mut Vec<OpenAllele>` rather than `&mut OpenPileupRecord`
/// so the caller can hold an independent borrow on the record's
/// `ref_seq` (e.g. for [`apply_events_into`]) at the same time.
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
    truncated_by_cap: &[u32],
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
                if event_end > cur_end
                    && open.widen(k, event_end, reference, active_reads, contributors)?
                {
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
        refold_ids_buf: _,
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
        // ends inside each `apply_events_into` call (which writes into
        // the table's scratch buffer), leaving the rest of the inner loop
        // free to mutate `alleles[other_idx]`. This is the same
        // disjoint-borrow shape Mi6 (`pileup_2026-05-09.md`)
        // introduced, with `ref_seq` collapsed onto `alleles[0]` —
        // the bytes are no longer duplicated on the record.
        let OpenPileupRecord {
            alleles,
            folded_reads,
            reads_without_observation,
            reads_discarded_by_cap,
            ..
        } = rec;
        // `saturating_add`, matching `footprint_end_exclusive()` — every other reader of
        // this quantity goes through that method, and Mi8 added the saturation there
        // because a wrapped end on a multi-Gbp chromosome yields an empty event window
        // and therefore a silently unfolded read. Open-coding a plain `+` here opted out
        // of the same defence for no reason.
        let rec_end = rec_pos.saturating_add(alleles[0].seq.len() as u32);
        let mut fold_state = RecordFoldState {
            alleles,
            folded_reads,
            reads_without_observation,
            reads_discarded_by_cap,
        };

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

            fold_read_into_record(
                &mut fold_state,
                allele_seq_buf,
                rec_pos,
                active_read,
                &window_events,
                contrib.bq_baq_at_walker_pos,
            );
        }

        // **The reads the cap removed at this position, registered against every record
        // affected here** (B3). Candidates, not the count: `finalise` keeps only those
        // still absent from `folded_reads`, because a read truncated at one position of a
        // footprint may fold at another — and if it folds at all, it folds with its whole
        // window. Registered *after* the fold loop, so a read that both contributed and was
        // truncated at this position could not be listed against a record it just folded
        // into. That cannot happen today (truncation removes it from `contributors` before
        // the fold ever sees it), and the ordering keeps it true if the cap ever moves.
        for &read_id in truncated_by_cap {
            note_no_observation(fold_state.reads_discarded_by_cap, read_id);
        }

        // After every contributor has folded, not before: the empty buckets are created
        // *by* the re-folds a widen triggers, and by the no-observation path above.
        // Moving this above the loop strands every bucket the fold loop empties, and no
        // fixture in this module catches that — the `debug_assert!` at the top of
        // `finalise` is what does, every walk in the suite. (It used to say a parity test
        // on the emitted records caught it; D1 mutated the code and found that it did not.)
        evict_unsupported_alleles(fold_state.alleles, fold_state.folded_reads);
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
    // **No copies of the read's own scalars.** Production carried `mq_log_err`, `mapq`,
    // `is_reverse_strand` and the read group here, duplicating fields the fold could read
    // off the active read it already looks up. A3's `widen` re-fold has no contributor to
    // read them from, so the fold takes them from the read directly — which removes the
    // duplicate rather than adding a second way to reach it.
    /// Active-set local id of the contributing read. Keys the
    /// per-record `folded_reads` map (so re-folds subtract the
    /// prior contribution) and is also how the fold looks the
    /// read up against the `ActiveReads` to query window events.
    pub read_id: u32,
    pub chain_id: ChainId,
    /// Events whose anchor *is* this walker_pos (used by step 3
    /// to identify candidate records). At most 2 events anchor at
    /// any walker_pos (one Match plus at most one indel), so the
    /// SmallVec keeps this list off the heap on every step.
    pub events_at_pos: super::cigar_cursor::EventsAt,
    /// BAQ-capped BQ at this walker position (the Match-event's
    /// quality, used as the fallback when the cursor's window
    /// returns no events — a clean REF read).
    pub bq_baq_at_walker_pos: u8,
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

    /// **`widen` grows the REF bucket and nothing else** (A3). Production appends the
    /// new reference bases to *every* bucket; ng appends them to `alleles[0]` alone,
    /// because that is the record's own reference sequence and the others hold what some
    /// read witnessed. The empty active set is honest here: no read has folded into this
    /// record, so there is nothing to re-fold.
    #[test]
    fn widen_extends_the_ref_bucket_and_leaves_the_others_alone() {
        let mut t = OpenPileupRecordTable::new();
        let f = fa("ACGTAC");
        let active = ActiveReads::new();
        // Open at pos 1 with span 1 ("A"), and give it a second bucket by hand.
        t.open_new(0, 1, 1, &f).unwrap();
        let rec = t.records.get_mut(&1).unwrap();
        rec.alleles.push(OpenAllele::new(b"T".to_vec()));
        rec.alleles[1].support.num_obs = 1;
        // Now widen to span 3 ("ACG").
        t.widen(1, 4, &f, &active, &[]).unwrap();
        let rec = t.records.get(&1).unwrap();
        assert_eq!(rec.alleles[0].seq, b"ACG", "the REF bucket grows");
        assert_eq!(
            rec.alleles[1].seq, b"T",
            "an allele bucket holds what a read witnessed, and the window growing around \
             it does not change that — production appends `CG` here"
        );
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

    // --- Precondition 4: the events' reach is non-decreasing -----
    //
    // Two `Deletion`s anchored at the same position, the longer
    // first, satisfy the sort precondition and reach backwards:
    // `deleted_len` 5 then 2 reach 106 then 103. The cursor cannot
    // emit that — a read deletes a run once — which is why the
    // `run_end.max(event_end)` that survives it is not reachable
    // from any real input, and why nothing failed when a reviewer
    // deleted the `.max()`. Both ends of the contract are pinned,
    // as `subtract_contribution`'s pair below does: in debug the
    // input trips the assertion, in release the `.max()` keeps the
    // run at what the read witnessed.
    #[cfg(debug_assertions)]
    #[test]
    #[should_panic(expected = "reach (anchor + footprint_span) must be non-decreasing")]
    fn apply_events_two_deletions_at_one_anchor_trip_the_reach_assertion_in_debug() {
        let ref_seq = b"ACGTACGTAC"; // positions 100..=109
        let events = vec![
            ReadEvent::Deletion {
                anchor_ref_pos: 100,
                deleted_len: 5,
                bq_proxy: 30,
            },
            ReadEvent::Deletion {
                anchor_ref_pos: 100,
                deleted_len: 2,
                bq_proxy: 30,
            },
        ];
        let _ = apply_events(100, ref_seq, &events);
    }

    #[cfg(not(debug_assertions))]
    #[test]
    fn apply_events_a_shorter_deletion_behind_a_longer_one_cannot_shrink_the_run_in_release() {
        let ref_seq = b"ACGTACGTAC"; // positions 100..=109
        let events = vec![
            ReadEvent::Deletion {
                anchor_ref_pos: 100,
                deleted_len: 5,
                bq_proxy: 30,
            },
            ReadEvent::Deletion {
                anchor_ref_pos: 100,
                deleted_len: 2,
                bq_proxy: 30,
            },
        ];
        let (bases, witnessed) = apply_events(100, ref_seq, &events).expect("one run");
        assert_eq!(bases, b"A", "the anchor base once, from the first deletion");
        assert_eq!(
            witnessed,
            RefSpan {
                start: 100,
                end: 105
            },
            "the run keeps the longer deletion's reach; a plain assignment would report 102",
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

        // On ng's own locus type — the back-projection this used to go through is gone.
        let records: Vec<_> = run(reads, MockFasta::new(contig), &WalkerConfig::default())
            .map(|record| record.expect("the walk succeeds"))
            .collect();
        let widened = records
            .iter()
            .find(|locus| locus.region.start.get() == 5)
            .expect("the opener's deletion opens a record at 5");
        assert_eq!(
            widened.reference_bases.len(),
            16,
            "the widener's deletion must grow the record to 5..=20, or this fixture does \
             not reach the path it exists for"
        );

        let folded: u32 = widened.observations.iter().map(|row| row.num_obs).sum();
        assert_eq!(
            folded, 2,
            "only `opener` and `widener` witnessed this record as one run; `wide` folded \
             into it before the widen and must have been subtracted back out when its \
             witness split in two. Record: {widened:?}"
        );
    }

    // -----------------------------------------------------------------
    // A3 — what `widen` does to the reads already in the record
    //
    // These drive `process_position` against the table directly rather than through
    // `run`, because the three properties A3 introduces are **not visible on a
    // `PileupRecord`**: `witnessed` is consumed at A4, a carried-over `contribution` is
    // indistinguishable from a recomputed one unless the two differ, and an evicted
    // bucket is one a comparison would drop anyway. Mutating each of the three left the
    // whole suite green, which is how they came to be written.
    // -----------------------------------------------------------------

    use super::super::chain_id_allocator::ChainIdAllocator;
    use super::super::{CigarOp, MateRole, PreparedRead};

    /// Admit reads into an active set and return it with its allocator.
    fn admitted(reads: Vec<PreparedRead>) -> (ActiveReads, ChainIdAllocator) {
        let mut active = ActiveReads::new();
        let mut chain_ids = ChainIdAllocator::new();
        for read in reads {
            active
                .admit(read, &mut chain_ids)
                .expect("the fixture reads are well formed");
        }
        (active, chain_ids)
    }

    /// The contributors at `pos`: every active read with an event anchored there, in
    /// admission order — what `genome_walk` builds, minus the mate-overlap resolution
    /// these fixtures do not need.
    fn contributors_at(active: &ActiveReads, pos: u32) -> Vec<ReadContribution> {
        active
            .iter()
            .filter_map(|entry| {
                let events_at_pos = entry.cursor.events_at(pos, &entry.read);
                if events_at_pos.is_empty() {
                    return None;
                }
                let bq = events_at_pos
                    .iter()
                    .find_map(|event| match event {
                        ReadEvent::Match { bq_baq, .. } => Some(*bq_baq),
                        _ => None,
                    })
                    .unwrap_or(0);
                Some(ReadContribution {
                    read_id: entry.read_id,
                    chain_id: entry.chain_id,
                    events_at_pos,
                    bq_baq_at_walker_pos: bq,
                    alignment_start: entry.read.alignment_start,
                    mate_role: entry.read.mate_role,
                    bq_zero_in_window: false,
                    bq_override_at_walker_pos: None,
                })
            })
            .collect()
    }

    fn plain_read(
        qname: &str,
        start: u32,
        end: u32,
        cigar: Vec<CigarOp>,
        seq: Vec<u8>,
    ) -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start: start,
            alignment_end: end,
            cigar,
            bq_baq: vec![30; seq.len()],
            seq,
            mq_log_err: -3.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: std::sync::Arc::from(qname),
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
            read_group: crate::ng::types::ReadGroupId(0),
        }
    }

    /// The contig every widen fixture below walks, 25 bases.
    const WIDEN_CONTIG: &str = "ACGTACGTACGTACGTACGTACGTA";

    /// The three reads the widen fixtures share.
    ///
    /// `spanner` matches every position; `opener`'s two-base deletion opens a record at 5
    /// spanning 5..=7; `widener`'s thirteen-base deletion anchored at 7 grows that record
    /// to 5..=20. At position 7 `spanner` and `widener` are contributors and `opener` is
    /// **not** — it is inside its own deletion — so the fixture holds both a read the fold
    /// loop re-folds and a read only `refold_live_reads` reaches.
    ///
    /// **The widener's matched base is `T` where the opener's is `A`, and that is
    /// load-bearing.** Both bases sit at a position where the reference reads `G`, so
    /// either is a mismatch and either opens a bucket. But with both `A` the widener folds
    /// into the very bucket the opener's re-fold has just left, refilling it in the same
    /// call — no bucket ever reaches `num_obs == 0`, and the eviction test below passes
    /// against an implementation that never evicts anything.
    fn widen_fixture_reads() -> Vec<PreparedRead> {
        vec![
            plain_read(
                "spanner",
                1,
                25,
                vec![CigarOp::Match(25)],
                WIDEN_CONTIG.as_bytes().to_vec(),
            ),
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(18)],
                vec![b'A'; 19],
            ),
            plain_read(
                "widener",
                7,
                24,
                vec![CigarOp::Match(1), CigarOp::Deletion(13), CigarOp::Match(4)],
                b"TAAAA".to_vec(),
            ),
        ]
    }

    /// The read id of the fixture read named `qname`.
    fn read_id_of(active: &ActiveReads, qname: &str) -> u32 {
        active
            .iter()
            .find(|entry| &*entry.read.qname == qname)
            .unwrap_or_else(|| panic!("{qname} is still active"))
            .read_id
    }

    /// [`widen_fixture_reads`] walked at 5 and then at 7, leaving one record at 5 that has
    /// widened once.
    fn widened_record() -> (OpenPileupRecordTable, ActiveReads) {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(widen_fixture_reads());
        let mut open = OpenPileupRecordTable::new();
        for pos in [5u32, 7] {
            let contributors = contributors_at(&active, pos);
            process_position(&mut open, pos, 0, &contributors, &[], &active, &reference)
                .expect("the fixture walks cleanly");
        }
        (open, active)
    }

    /// **`witnessed` follows the record when it widens** — the property option (b) exists
    /// for, and the only observable effect it has.
    ///
    /// `opener` is inside its own deletion at position 7, so it is not a contributor
    /// there and the fold loop never revisits it. Without the live re-fold its extent
    /// would stay `5..=7`, measured against a footprint the record has outgrown — and
    /// A4 resolves `ReadWitness` from exactly that, so the read would be reported as
    /// having seen three positions of sixteen when its deletion witnessed all of them.
    /// A wrong depth, with no error.
    #[test]
    fn widening_updates_the_witnessed_extent_of_a_read_that_is_not_a_contributor() {
        let (open, active) = widened_record();
        let record = open.records.get(&5).expect("the record at 5");
        assert_eq!(
            record.alleles[0].seq.len(),
            16,
            "the widener's deletion must grow the record to 5..=20, or this fixture does \
             not reach the path it exists for"
        );
        let opener = record
            .folded_reads
            .get(&read_id_of(&active, "opener"))
            .expect("the opener folded into this record at position 5");
        assert_eq!(
            opener.witnessed,
            RefSpan { start: 5, end: 20 },
            "the opener's deletion witnesses 6 and 7 and its matches witness 8..=20, so \
             its extent covers the widened footprint whole"
        );
    }

    /// **The re-placement carries the read's contribution; it does not recompute it.**
    ///
    /// A widen changes which positions the record covers, not what quality evidence a
    /// read carried. `q_sum` in particular encodes decisions the walk took at *earlier*
    /// positions and replayed into the fold — a mate-overlap loser's zeroed Match, an
    /// agree-case keeper's summed BQ — and those live on the `ReadContribution` of the
    /// position they were taken at, which is gone by the time a later widen fires.
    ///
    /// Here the opener's contribution is stamped with a quality no recomputation could
    /// produce from its events, and it has to survive the widen unchanged.
    #[test]
    fn widening_carries_a_reads_contribution_rather_than_recomputing_it() {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(widen_fixture_reads());
        let mut open = OpenPileupRecordTable::new();
        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference).expect("opens");

        let opener_id = read_id_of(&active, "opener");
        // A quality the walk decided and the events cannot reproduce.
        const RECONCILED: f64 = -0.5;
        {
            let record = open.records.get_mut(&5).expect("the record at 5");
            let state = record
                .folded_reads
                .get_mut(&opener_id)
                .expect("the opener folded");
            let was = state.contribution.q_sum;
            record.alleles[state.allele_index].support.q_sum += RECONCILED - was;
            state.contribution.q_sum = RECONCILED;
        }

        let contributors = contributors_at(&active, 7);
        process_position(&mut open, 7, 0, &contributors, &[], &active, &reference).expect("widens");

        let record = open.records.get(&5).expect("the record at 5");
        assert_eq!(
            record.alleles[0].seq.len(),
            16,
            "the record must have widened"
        );
        let opener = record
            .folded_reads
            .get(&opener_id)
            .expect("the opener is still folded");
        assert_eq!(
            opener.contribution.q_sum, RECONCILED,
            "a widen re-places the read; recomputing its quality here would silently undo \
             every reconciliation the walk made earlier in this record's footprint"
        );
    }

    /// **A bucket no read is folded into is evicted, and the fold state is remapped.**
    ///
    /// With REF-only widening a re-folding read lands in a new bucket and leaves its old
    /// one behind at `num_obs == 0`. `find_allele_index` is a linear scan with a full
    /// byte compare, run once per (record, contributor) at every position of the
    /// footprint, so those accumulate against exactly the long-deletion loci this port
    /// exists to fix.
    ///
    /// The assertion is **which bucket each read sits in**, by bytes, because that is the
    /// only form that survives both ways of getting this wrong. Skipping the eviction
    /// leaves the opener's pre-widen `b"A"` bucket behind at `num_obs == 0`; skipping the
    /// *remap* leaves the opener's `allele_index` one bucket too high — pointing at the
    /// widener's, which does have support, so "every folded read points at a bucket with
    /// observations" would wave it straight through.
    #[test]
    fn widening_evicts_the_buckets_its_re_folds_emptied() {
        let (open, active) = widened_record();
        let record = open.records.get(&5).expect("the record at 5");
        assert_eq!(
            record.alleles[0].seq.len(),
            16,
            "the record must have widened"
        );

        let buckets: Vec<_> = record
            .alleles
            .iter()
            .map(|a| {
                (
                    String::from_utf8_lossy(&a.seq).to_string(),
                    a.support.num_obs,
                )
            })
            .collect();
        let bucket_of = |qname: &str| -> &OpenAllele {
            let state = record
                .folded_reads
                .get(&read_id_of(&active, qname))
                .unwrap_or_else(|| panic!("{qname} folded into this record"));
            &record.alleles[state.allele_index]
        };

        assert_eq!(
            bucket_of("spanner").seq.as_slice(),
            &WIDEN_CONTIG.as_bytes()[4..20],
            "the spanner matched the widened footprint whole, so it is in REF: {buckets:?}"
        );
        assert_eq!(
            bucket_of("opener").seq.as_slice(),
            [b'A'; 14].as_slice(),
            "the opener's anchor base plus its thirteen matches beyond the deletion: \
             {buckets:?}"
        );
        assert_eq!(
            bucket_of("widener").seq.as_slice(),
            b"T".as_slice(),
            "the widener witnessed one base and deleted the rest: {buckets:?}"
        );
        assert_eq!(
            record.alleles.len(),
            3,
            "REF, the opener's and the widener's — the opener's pre-widen b\"A\" bucket \
             was emptied by its re-fold and must be gone: {buckets:?}"
        );
        assert!(
            record.alleles[1..].iter().all(|a| a.support.num_obs > 0),
            "every non-REF bucket must support a read: {buckets:?}"
        );
    }

    // -----------------------------------------------------------------
    // A4 — coverage resolved at `finalise`, against the final footprint.
    // -----------------------------------------------------------------

    /// A witness that tiles the footprint is a complete one.
    #[test]
    fn coverage_of_a_witness_covering_the_whole_footprint_is_complete() {
        assert_eq!(
            coverage_of(RefSpan { start: 5, end: 20 }, 5, 21),
            ReadWitness::Complete
        );
    }

    /// Flush with the left border and short of the right — freebayes' prefix constraint,
    /// which `is_flush_left` has to keep being able to read off the run.
    #[test]
    fn coverage_of_a_witness_flush_left_reports_a_zero_offset() {
        let coverage = coverage_of(RefSpan { start: 5, end: 7 }, 5, 21);
        assert_eq!(
            coverage,
            ReadWitness::Observed {
                offset_in_locus: 0,
                positions_covered: 3,
            }
        );
        assert!(coverage.is_flush_left());
        assert!(!coverage.is_flush_right(LocusLen::from_positions(16)));
    }

    /// Flush with the right border — the suffix constraint, and the offset is derived
    /// rather than assumed.
    #[test]
    fn coverage_of_a_witness_flush_right_reports_the_offset_it_starts_at() {
        let coverage = coverage_of(RefSpan { start: 18, end: 20 }, 5, 21);
        assert_eq!(
            coverage,
            ReadWitness::Observed {
                offset_in_locus: 13,
                positions_covered: 3,
            }
        );
        assert!(!coverage.is_flush_left());
        assert!(coverage.is_flush_right(LocusLen::from_positions(16)));
    }

    /// **A run flush with neither border**, which is what the generic path mints and
    /// neither `from_left` nor `from_right` can express — the case
    /// [`ReadWitness`](super::super::ReadWitness)'s own note said would only be
    /// knowable when this generator produced its first run.
    #[test]
    fn coverage_of_an_interior_witness_is_flush_with_neither_border() {
        let coverage = coverage_of(RefSpan { start: 9, end: 12 }, 5, 21);
        assert_eq!(
            coverage,
            ReadWitness::Observed {
                offset_in_locus: 4,
                positions_covered: 4,
            }
        );
        assert!(!coverage.is_flush_left());
        assert!(!coverage.is_flush_right(LocusLen::from_positions(16)));
    }

    /// **An extent reaching past either border is clamped into the footprint, not
    /// believed.** `events_overlapping` does not clip a deletion — one anchored before
    /// the record comes back whole, so its run can start below `record_pos` and end past
    /// `record_end` (spec §8). Unclamped, a deletion spanning the record would report an
    /// enormous `positions_covered`, or an `offset_in_locus` that underflowed.
    #[test]
    fn coverage_of_clamps_an_extent_that_overruns_the_footprint() {
        assert_eq!(
            coverage_of(RefSpan { start: 1, end: 40 }, 5, 21),
            ReadWitness::Complete,
            "a witness swallowing the footprint witnessed all of it, and nothing more"
        );
        assert_eq!(
            coverage_of(RefSpan { start: 1, end: 7 }, 5, 21),
            ReadWitness::Observed {
                offset_in_locus: 0,
                positions_covered: 3,
            },
            "the run starts at the record's own anchor, never before it"
        );
        assert_eq!(
            coverage_of(RefSpan { start: 18, end: 40 }, 5, 21),
            ReadWitness::Observed {
                offset_in_locus: 13,
                positions_covered: 3,
            },
            "and it ends at the record's own end"
        );
    }

    /// **A read that was a complete witness becomes `Observed` when the record widens
    /// under it, with nothing about the read having changed** — the whole reason
    /// coverage is resolved at `finalise` and not at the fold (spec §4, plan A4).
    ///
    /// `shortie` matches positions 5..=7 and stops. When the record spans 5..=7 that is
    /// every position it has. `widener`'s deletion then grows the record to 5..=20, and
    /// the same three positions are now three of sixteen. The read is not consulted; its
    /// `witnessed` extent is byte for byte what it was.
    #[test]
    fn a_complete_witness_becomes_observed_when_the_record_widens_under_it() {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(vec![
            plain_read(
                "shortie",
                5,
                7,
                vec![CigarOp::Match(3)],
                WIDEN_CONTIG.as_bytes()[4..7].to_vec(),
            ),
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(18)],
                vec![b'A'; 19],
            ),
            plain_read(
                "widener",
                7,
                24,
                vec![CigarOp::Match(1), CigarOp::Deletion(13), CigarOp::Match(4)],
                b"TAAAA".to_vec(),
            ),
        ]);
        let mut open = OpenPileupRecordTable::new();

        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference).expect("opens");
        let shortie_id = read_id_of(&active, "shortie");
        let narrow = {
            let record = open.records.get(&5).expect("the record at 5");
            assert_eq!(
                record.ref_span(),
                3,
                "the opener's deletion must open the record at 5..=7"
            );
            let state = record
                .folded_reads
                .get(&shortie_id)
                .expect("the shortie folded at 5");
            assert_eq!(
                coverage_of(
                    state.witnessed,
                    record.pos,
                    record.footprint_end_exclusive()
                ),
                ReadWitness::Complete,
                "against a 5..=7 footprint the shortie witnessed everything"
            );
            state.witnessed
        };

        let contributors = contributors_at(&active, 7);
        process_position(&mut open, 7, 0, &contributors, &[], &active, &reference).expect("widens");

        let record = open.records.get(&5).expect("the record at 5");
        assert_eq!(
            record.ref_span(),
            16,
            "the record must have widened to 5..=20"
        );
        let state = record
            .folded_reads
            .get(&shortie_id)
            .expect("the shortie is still folded");
        assert_eq!(
            state.witnessed, narrow,
            "the read saw exactly what it saw; a widen is not news about the read"
        );
        assert_eq!(
            coverage_of(
                state.witnessed,
                record.pos,
                record.footprint_end_exclusive()
            ),
            ReadWitness::Observed {
                offset_in_locus: 0,
                positions_covered: 3,
            },
            "the same extent, resolved against the footprint the record ended with"
        );
    }

    /// **`finalise` resolves every folded read against the record's final footprint**,
    /// and reports what it found. Same walk as the test above: the shortie stopped at 7,
    /// the widener started there, and only the opener's deletion carried it across the
    /// whole widened footprint.
    ///
    /// A `finalise` that answered from the *fold-time* footprint would call the shortie
    /// complete, and the number it feeds — a depth — would be wrong with no error.
    #[test]
    fn finalise_counts_the_witnesses_against_the_footprint_the_record_ended_with() {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(vec![
            plain_read(
                "shortie",
                5,
                7,
                vec![CigarOp::Match(3)],
                WIDEN_CONTIG.as_bytes()[4..7].to_vec(),
            ),
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(18)],
                vec![b'A'; 19],
            ),
            plain_read(
                "widener",
                7,
                24,
                vec![CigarOp::Match(1), CigarOp::Deletion(13), CigarOp::Match(4)],
                b"TAAAA".to_vec(),
            ),
        ]);
        let mut open = OpenPileupRecordTable::new();
        for pos in [5u32, 7] {
            let contributors = contributors_at(&active, pos);
            process_position(&mut open, pos, 0, &contributors, &[], &active, &reference)
                .expect("the fixture walks cleanly");
        }

        let record = open.records.remove(&5).expect("the record at 5");
        let (_record, witness) = record.finalise();
        assert_eq!(
            witness,
            RecordWitness {
                reads_complete: 1,
                reads_partially_observed: 2,
                reads_without_observation: 0,
                reads_discarded_by_cap: 0,
            },
            "the opener's deletion witnessed 5..=20 whole; the shortie stopped at 7 and \
             the widener started there"
        );
    }

    // -----------------------------------------------------------------
    // B1 — the row identity: bases, coverage, read group.
    // -----------------------------------------------------------------

    /// A record at 7 spanning 7..=20, with `shortie` and `deleter` both showing the single
    /// base `G` — and showing it from **different amounts of the locus**. `shortie`'s
    /// alignment ends at 7, so it witnessed one position of fourteen; `deleter`'s deletion
    /// runs to 20, so it witnessed every one of them and still emits one base.
    ///
    /// The two share a *bucket*, because a bucket is keyed on bases alone. Whether they
    /// share a **row** is what B1 decides.
    fn same_bases_different_coverage() -> (OpenPileupRecordTable, ActiveReads) {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(vec![
            plain_read(
                "shortie",
                5,
                7,
                vec![CigarOp::Match(3)],
                WIDEN_CONTIG.as_bytes()[4..7].to_vec(),
            ),
            plain_read(
                "deleter",
                5,
                24,
                vec![CigarOp::Match(3), CigarOp::Deletion(13), CigarOp::Match(4)],
                WIDEN_CONTIG.as_bytes()[4..7]
                    .iter()
                    .copied()
                    .chain(std::iter::repeat_n(b'A', 4))
                    .collect(),
            ),
        ]);
        let mut open = OpenPileupRecordTable::new();
        let contributors = contributors_at(&active, 7);
        process_position(&mut open, 7, 0, &contributors, &[], &active, &reference)
            .expect("the fixture walks cleanly");
        (open, active)
    }

    /// **A complete witness and a partial one of the same bases are different evidence, and
    /// stay different rows** (spec §3, §6) — even though they share one allele bucket,
    /// because a bucket is keyed on bases alone and cannot tell them apart.
    ///
    /// This is what "rows are re-derived per read, not read off the bucket totals" buys.
    /// Read the bucket and there is one observation of `G` with `num_obs == 2`; read the
    /// reads and there are two, one of which saw one position of fourteen and is a lower
    /// bound on nothing more than that.
    #[test]
    fn rows_split_a_complete_witness_from_a_partial_one_of_the_same_bases() {
        let (open, _active) = same_bases_different_coverage();
        let record = open.records.get(&7).expect("the record at 7");
        assert_eq!(
            record.ref_span(),
            14,
            "the deleter's deletion must widen the record to 7..=20"
        );

        let bucket = record
            .alleles
            .iter()
            .find(|allele| allele.seq == b"G")
            .expect("both reads showed a lone `G`");
        assert_eq!(
            bucket.support.num_obs, 2,
            "and they share one bucket, which is the premise of this test"
        );

        let mut rows = record.observation_rows(record.footprint_end_exclusive());
        rows.sort_by_key(|row| match row.key.read_witness {
            ReadWitness::Complete => 0u16,
            ReadWitness::Observed {
                positions_covered, ..
            } => positions_covered,
        });
        let split: Vec<_> = rows
            .iter()
            .filter(|row| row.key.bases == b"G")
            .map(|row| (row.key.read_witness, row.support.num_obs))
            .collect();
        assert_eq!(
            split,
            vec![
                (ReadWitness::Complete, 1),
                (
                    ReadWitness::Observed {
                        offset_in_locus: 0,
                        positions_covered: 1,
                    },
                    1
                ),
            ],
            "one bucket, two rows: the deleter witnessed all fourteen positions and the \
             shortie one. Rows: {rows:?}"
        );
    }

    /// **Two read groups supporting one sequence are two rows, and they sum to what one
    /// group would have shown** (spec §13.1, fixture 5).
    ///
    /// The sum is the half that matters: splitting is only free if nothing is lost, and a
    /// per-group model reading these two rows must be able to recover the merged answer
    /// exactly — that is what makes the grain the *consumer's* choice rather than this
    /// step's guess.
    #[test]
    fn rows_split_when_two_read_groups_support_one_sequence() {
        let reference = fa(WIDEN_CONTIG);
        // Two reads, identical in every way a row's identity can see except the lane.
        let mut first = plain_read("lane0", 5, 5, vec![CigarOp::Match(1)], b"T".to_vec());
        let mut second = plain_read("lane1", 5, 5, vec![CigarOp::Match(1)], b"T".to_vec());
        first.read_group = crate::ng::types::ReadGroupId(0);
        second.read_group = crate::ng::types::ReadGroupId(1);
        let (active, _ids) = admitted(vec![first, second]);
        let mut open = OpenPileupRecordTable::new();
        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference).expect("walks");

        let record = open.records.get(&5).expect("the record at 5");
        let rows = record.observation_rows(record.footprint_end_exclusive());
        let mut alt: Vec<_> = rows.iter().filter(|row| row.key.bases == b"T").collect();
        alt.sort_by_key(|row| row.key.read_group.0);
        assert_eq!(
            alt.len(),
            2,
            "one mismatched base from two lanes is two rows: {rows:?}"
        );
        assert_eq!(
            (alt[0].key.read_group.0, alt[1].key.read_group.0),
            (0, 1),
            "and each names its own lane"
        );
        assert_eq!(
            (alt[0].support.num_obs, alt[1].support.num_obs),
            (1, 1),
            "one read each"
        );
        let bucket = record
            .alleles
            .iter()
            .find(|allele| allele.seq == b"T")
            .expect("both reads showed `T`");
        assert_eq!(
            alt[0].support.num_obs + alt[1].support.num_obs,
            bucket.support.num_obs,
            "the two rows must sum to the single-group total, or the split loses evidence"
        );
        assert_eq!(
            alt[0].support.mapq_sum + alt[1].support.mapq_sum,
            bucket.support.mapq_sum,
            "and the quality moments too — a count beside a merged observation is exactly \
             what carrying the group is meant to replace"
        );
    }

    /// **At one read group the rows are the buckets** — "free at one read group" stated as
    /// a test rather than as a hope (plan B1).
    ///
    /// The general fixture's reads all carry `ReadGroupId(0)`, so every split must come
    /// from coverage, never from the group; and where coverage does not split either, the
    /// row count equals the bucket count with equal support. This is the property that
    /// keeps the stage-1 differential green across B1.
    #[test]
    fn rows_are_the_buckets_when_one_read_group_witnesses_completely() {
        let (open, _active) = widened_record();
        let record = open.records.get(&5).expect("the record at 5");
        let rows = record.observation_rows(record.footprint_end_exclusive());

        for row in &rows {
            assert_eq!(
                row.key.read_group,
                crate::ng::types::ReadGroupId(0),
                "this fixture has one lane, so no row may name another"
            );
        }
        // Every bucket that supports a read has at least one row, and the rows over a
        // bucket sum to it. Matched on **bases**, which is a bijection with the bucket —
        // `find_allele_index` never creates two buckets with the same bytes.
        for allele in record.alleles.iter() {
            if allele.support.num_obs == 0 {
                continue;
            }
            let over_bucket: u32 = rows
                .iter()
                .filter(|row| row.key.bases == allele.seq)
                .map(|row| row.support.num_obs)
                .sum();
            assert_eq!(
                over_bucket,
                allele.support.num_obs,
                "bucket {:?} has {} observations but its rows carry {over_bucket}",
                String::from_utf8_lossy(&allele.seq),
                allele.support.num_obs,
            );
        }
    }

    /// **Reads that share a row identity are merged into one row** — the half of B1 the
    /// three fixtures above do not test, and the half the differential is structurally
    /// blind to.
    ///
    /// `to_pileup_record` — the back-projection the suite used until Milestone D, now deleted —
    /// merged rows back together by bases before the two walkers were compared, so an
    /// `observation_rows` that emitted one row per read — every `num_obs == 1` — projected to
    /// exactly the same `PileupRecord` and left the whole suite green, at 20,000 soak cases as
    /// well. The projection undoes precisely the
    /// defect. So the merge has to be asserted here, on ng's own type, or not at all.
    #[test]
    fn rows_merge_the_reads_that_share_an_identity() {
        let reference = fa(WIDEN_CONTIG);
        // Three reads, same lane, same single mismatched base, all complete witnesses of a
        // one-base record: one identity, three reads.
        let (active, _ids) = admitted(
            (0..3)
                .map(|index| {
                    plain_read(
                        &format!("same{index}"),
                        5,
                        5,
                        vec![CigarOp::Match(1)],
                        b"T".to_vec(),
                    )
                })
                .collect(),
        );
        let mut open = OpenPileupRecordTable::new();
        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference)
            .expect("the fixture walks cleanly");

        let record = open.records.remove(&5).expect("the record at 5");
        let rows = record.observation_rows(record.footprint_end_exclusive());
        let alt: Vec<_> = rows.iter().filter(|row| row.key.bases == b"T").collect();
        assert_eq!(
            alt.len(),
            1,
            "three reads with one identity are one row, not three: {rows:?}"
        );
        assert_eq!(
            alt[0].support.num_obs, 3,
            "and the row carries all three, which is what a row *is*"
        );
        assert_eq!(
            alt[0].chain_ids.len(),
            3,
            "each read keeps its own identity inside the row it shares"
        );
    }

    /// **The emitted region covers exactly the record's footprint, inclusive of both ends.**
    ///
    /// Unpinned by anything before this: `to_pileup_record` — the back-projection the suite used
    /// until Milestone D, now deleted — carried the contig and the start and **discarded the
    /// end**, so the 44 inherited tests and the whole differential were blind to it — dropping the `saturating_sub(1)`, or replacing the end with
    /// `Position(0)`, left the entire library green. `region.len()` is what sizes
    /// `num_obs_along_locus`'s depth vector and what defines flush-right, so a wrong end is
    /// a wrong depth profile everywhere downstream.
    #[test]
    fn the_emitted_region_covers_the_footprint_inclusively() {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(vec![plain_read(
            "solo",
            5,
            5,
            vec![CigarOp::Match(1)],
            b"A".to_vec(),
        )]);
        let mut open = OpenPileupRecordTable::new();
        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference).expect("walks");
        let (one_base, _) = open.records.remove(&5).expect("the record at 5").finalise();
        assert_eq!(one_base.region.start.get(), 5);
        assert_eq!(
            one_base.region.end.get(),
            5,
            "a one-base record starts and ends on the same position — `GenomeRegion` is \
             inclusive of both ends (spec §6)"
        );
        assert_eq!(one_base.region.len(), 1);
        assert_eq!(
            one_base.region.len() as usize,
            one_base.reference_bases.len(),
            "the region and the reference bytes describe the same stretch"
        );

        let (widened, _) = {
            let (mut open, _active) = widened_record();
            open.records.remove(&5).expect("the record at 5").finalise()
        };
        assert_eq!(
            (widened.region.start.get(), widened.region.end.get()),
            (5, 20),
            "the widened record runs 5..=20"
        );
        assert_eq!(widened.region.len(), 16);
        assert_eq!(
            widened.region.len() as usize,
            widened.reference_bases.len(),
            "and its reference bytes are the same sixteen"
        );
    }

    /// **A read that departed from the reference keeps its chain id even on a partial row,
    /// and one that agreed carries none even when its row is not the REF row** — the
    /// per-read rule B2 replaced production's positional one with.
    ///
    /// Asserted nowhere before this, and the reason is structural rather than a missing
    /// fixture: replacing the whole rule with production's `allele_index == 0` left 158/158
    /// green at 10,000 cases, and so did making every partial row lose its ids. The
    /// differential compares chain ids by equality only on the complete-reads fixture, where
    /// the two rules coincide by construction, and the census that does see partial rows
    /// never looks at ids at all.
    ///
    /// `shortie` matched positions 5..=7 and stopped, so its row is a *partial* one whose
    /// bases equal the reference over what it saw — production would give it an id, and ng
    /// must not. `deleter` witnessed the whole footprint through a deletion, so its bases
    /// are not the reference and it keeps its id.
    #[test]
    fn only_the_reads_that_departed_from_the_reference_carry_a_chain_id() {
        let (open, active) = same_bases_different_coverage();
        let record = open.records.get(&7).expect("the record at 7");
        let rows = record.observation_rows(record.footprint_end_exclusive());

        // The shortie agreed with the reference across its one witnessed position; the
        // deleter's fourteen-position witness emits a single base that is not those
        // fourteen reference bytes.
        let shortie_row = rows
            .iter()
            .find(|row| matches!(row.key.read_witness, ReadWitness::Observed { .. }))
            .expect("the shortie's partial row");
        let deleter_row = rows
            .iter()
            .find(|row| row.key.read_witness == ReadWitness::Complete && row.key.bases == b"G")
            .expect("the deleter's complete row");
        assert_eq!(
            shortie_row.chain_ids,
            Vec::new(),
            "the shortie agreed with the reference across everything it witnessed, so it \
             carries no id — production's positional rule would have given it one, because \
             its partial row is not `alleles[0]`"
        );
        assert_eq!(
            deleter_row.chain_ids.len(),
            1,
            "the deleter deleted thirteen reference bases, so it departed and keeps its id"
        );
        let _ = read_id_of(&active, "deleter");
    }

    // -----------------------------------------------------------------
    // B3 — reads a depth cap discarded, per record.
    // -----------------------------------------------------------------

    /// The two halves the walk hands `process_position`: the contributors that survived the
    /// column-depth cap, and the read ids it removed. Splitting them here rather than
    /// truncating a list keeps the fixtures explicit about *which* read the cap took.
    fn contributors_at_with_cap(
        active: &ActiveReads,
        pos: u32,
        capped: &[u32],
    ) -> (Vec<ReadContribution>, Vec<u32>) {
        let (kept, removed): (Vec<_>, Vec<_>) = contributors_at(active, pos)
            .into_iter()
            .partition(|contrib| !capped.contains(&contrib.read_id));
        (kept, removed.into_iter().map(|c| c.read_id).collect())
    }

    /// A record at 5 spanning 5..=9, opened by `opener`'s four-base deletion, with `capped`
    /// matching every one of those positions. Walked at 5, 6, 8 and 9, with `capped`
    /// removed by the cap at whichever of those `capped_at` names.
    fn record_with_a_capped_read(capped_at: &[u32]) -> (OpenPileupRecordTable, ActiveReads) {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(vec![
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(16)],
                vec![b'A'; 17],
            ),
            plain_read(
                "capped",
                5,
                9,
                vec![CigarOp::Match(5)],
                WIDEN_CONTIG.as_bytes()[4..9].to_vec(),
            ),
        ]);
        let capped_id = read_id_of(&active, "capped");
        let mut open = OpenPileupRecordTable::new();
        for pos in [5u32, 6, 8, 9] {
            let capped_here: &[u32] = if capped_at.contains(&pos) {
                std::slice::from_ref(&capped_id)
            } else {
                &[]
            };
            let (contributors, truncated) = contributors_at_with_cap(&active, pos, capped_here);
            process_position(
                &mut open,
                pos,
                0,
                &contributors,
                &truncated,
                &active,
                &reference,
            )
            .expect("the fixture walks cleanly");
        }
        (open, active)
    }

    /// **A read the cap removed at every position it had events is reported as discarded** —
    /// which is what "the support counts are a subsample, not the depth" means (spec §6).
    #[test]
    fn a_read_the_cap_removed_everywhere_is_reported_as_discarded() {
        let (mut open, _active) = record_with_a_capped_read(&[5, 6, 8, 9]);
        let record = open.records.remove(&5).expect("the record at 5");
        assert_eq!(
            record.ref_span(),
            5,
            "the opener's deletion must open the record at 5..=9"
        );
        let (locus, witness) = record.finalise();
        assert_eq!(
            witness.reads_discarded_by_cap, 1,
            "`capped` had events at every position of this footprint and the cap took it at \
             all of them, so it folded nowhere"
        );
        assert_eq!(
            locus.reads_discarded_by_cap, 1,
            "and the locus reports it, which is what a model reads to know the counts are a \
             subsample"
        );
    }

    /// **A read that folded and then lost its row to a hole is not counted twice.**
    ///
    /// "Absent from `folded_reads`" has two causes — the cap kept the read out, or its
    /// witness turned out non-contiguous and A5's path removed it — and counting the second
    /// here reported one read in *both* `reads_without_observation` and
    /// `reads_discarded_by_cap`. Measured at 240 records in ~506,000 before the exclusion.
    /// The two counters mean different things to a model: one says the support is a
    /// subsample of the depth, the other says a read covered the locus and said nothing
    /// usable.
    #[test]
    fn a_read_that_lost_its_row_to_a_hole_is_not_also_counted_as_capped() {
        let reference = fa(WIDEN_CONTIG);
        // `holey` matches 5..=9 with an `N` at 7, so it folds at 5 and then yields no
        // observation — A5's path. It is *also* named to the cap, which the walk would do
        // if the column were over depth.
        let (active, _ids) = admitted(vec![
            plain_read("holey", 5, 9, vec![CigarOp::Match(5)], b"ACNTA".to_vec()),
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(16)],
                vec![b'A'; 17],
            ),
        ]);
        let holey_id = read_id_of(&active, "holey");
        let mut open = OpenPileupRecordTable::new();
        for pos in [5u32, 6, 8, 9] {
            let contributors = contributors_at(&active, pos);
            process_position(
                &mut open,
                pos,
                0,
                &contributors,
                std::slice::from_ref(&holey_id),
                &active,
                &reference,
            )
            .expect("the fixture walks cleanly");
        }

        let record = open.records.remove(&5).expect("the record at 5");
        let (locus, witness) = record.finalise();
        assert_eq!(
            witness.reads_without_observation, 1,
            "the hole at 7 is why `holey` has no row"
        );
        assert_eq!(
            witness.reads_discarded_by_cap, 0,
            "and it must not *also* be reported as discarded by the cap — one read, one \
             reason, or a model reads the support as a subsample when it is not"
        );
        assert_eq!(locus.reads_without_observation, 1);
        assert_eq!(locus.reads_discarded_by_cap, 0);
    }

    /// **A read the cap removed at one position but not another is *not* discarded** — the
    /// correction that makes this quantity worth having.
    ///
    /// Production counts truncated *positions*, run-wide. Counting them per record would
    /// flag this record, whose support is complete: `capped` survived the cap at 6 and, like
    /// every read that folds at all, folded with its **whole window**. Its evidence is not a
    /// subsample of anything. That is why the fold keeps a membership list and `finalise`
    /// resolves it against `folded_reads`, rather than incrementing a counter at the cap.
    #[test]
    fn a_read_the_cap_removed_at_only_some_positions_is_not_discarded() {
        let (mut open, active) = record_with_a_capped_read(&[5, 8, 9]);
        let record = open.records.remove(&5).expect("the record at 5");
        assert!(
            record
                .folded_reads
                .contains_key(&read_id_of(&active, "capped")),
            "`capped` must have folded at 6, or this fixture is the other test"
        );
        let (locus, witness) = record.finalise();
        assert_eq!(
            witness.reads_discarded_by_cap, 0,
            "it was truncated at three positions of four and folded anyway; a per-position \
             counter would have reported three"
        );
        assert_eq!(locus.reads_discarded_by_cap, 0);
    }

    // -----------------------------------------------------------------
    // A5 — the no-observation path: which reads, as a set.
    // -----------------------------------------------------------------

    /// **A read that witnesses nothing contiguous is recorded once, not once per
    /// position** — the reason the field is a set of read ids and not a counter
    /// (spec §4).
    ///
    /// The path is reached at *every* position the record is affected at. `holey`
    /// matches 5..=9 with an `N` at 7, so it is a contributor at four of the record's
    /// five positions and yields no observation at every one of them. A counter
    /// incremented there would report **four** reads without an observation where
    /// there is one — a number multiplied by the footprint's length, on the one path
    /// with no inherited test to catch it.
    #[test]
    fn a_read_witnessing_nothing_contiguous_is_recorded_once_not_once_per_position() {
        let reference = fa(WIDEN_CONTIG);
        // `ACNTA` over 5..=9: the reference reads `ACGTA` there, and the `N` at 7
        // makes the cursor emit nothing — a hole the read's *alignment span* is blind
        // to, which is why coverage comes from the events.
        let (active, _ids) = admitted(vec![
            plain_read("holey", 5, 9, vec![CigarOp::Match(5)], b"ACNTA".to_vec()),
            plain_read(
                "opener",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(16)],
                vec![b'A'; 17],
            ),
        ]);
        let mut open = OpenPileupRecordTable::new();
        // 7 is absent on purpose: the `N` means `holey` has no event there, so the
        // walker would not stop for it either.
        for pos in [5u32, 6, 8, 9] {
            let contributors = contributors_at(&active, pos);
            process_position(&mut open, pos, 0, &contributors, &[], &active, &reference)
                .expect("the fixture walks cleanly");
        }

        let record = open.records.remove(&5).expect("the record at 5");
        assert_eq!(
            record.ref_span(),
            5,
            "the opener's deletion must open the record at 5..=9, or the hole at 7 is \
             not inside the footprint and the fixture reaches nothing"
        );
        assert_eq!(
            record.reads_without_observation.len(),
            1,
            "the same read at four positions is one read: {:?}",
            record.reads_without_observation
        );
        let (_record, witness) = record.finalise();
        assert_eq!(
            witness.reads_without_observation, 1,
            "a counter incremented on the path would report one per affected position"
        );
        assert_eq!(
            witness.reads_complete, 1,
            "the opener's deletion witnessed the whole footprint"
        );
    }

    /// **A read whose witness splits when the record widens is recorded, not merely
    /// dropped.**
    ///
    /// `holed`'s deletion covers 6 and 7, so it witnesses 5..=7 as one run and — being
    /// inside that deletion at the widening position — is not a contributor there.
    /// Only [`refold_live_reads`] reaches it. `widener` then grows the record to
    /// 5..=20, which brings `holed`'s `N` at 11 inside the footprint and splits its
    /// witness in two. It leaves its bucket, and it has to leave a record of itself
    /// behind: a read that had a row a moment ago now has none, and nothing else in
    /// the output says so.
    #[test]
    fn a_read_whose_witness_splits_at_a_widen_is_recorded_not_merely_dropped() {
        let reference = fa(WIDEN_CONTIG);
        // Match(1) at 5, deletion over 6..=7, then matches from 8. Position 11 is the
        // fifth of those matched bases and is `N`.
        let mut holed_seq = vec![b'A'; 19];
        holed_seq[4] = b'N';
        let (active, _ids) = admitted(vec![
            plain_read(
                "holed",
                5,
                25,
                vec![CigarOp::Match(1), CigarOp::Deletion(2), CigarOp::Match(18)],
                holed_seq,
            ),
            plain_read(
                "widener",
                7,
                24,
                vec![CigarOp::Match(1), CigarOp::Deletion(13), CigarOp::Match(4)],
                b"TAAAA".to_vec(),
            ),
        ]);
        let mut open = OpenPileupRecordTable::new();

        let contributors = contributors_at(&active, 5);
        process_position(&mut open, 5, 0, &contributors, &[], &active, &reference).expect("opens");
        let holed_id = read_id_of(&active, "holed");
        assert!(
            open.records
                .get(&5)
                .expect("the record at 5")
                .folded_reads
                .contains_key(&holed_id),
            "before the widen `holed` witnessed 5..=7 as one run and has a row"
        );

        let contributors = contributors_at(&active, 7);
        assert!(
            !contributors.iter().any(|c| c.read_id == holed_id),
            "`holed` must be inside its own deletion at 7, or the fold loop reaches it \
             and `refold_live_reads` is not the path under test"
        );
        process_position(&mut open, 7, 0, &contributors, &[], &active, &reference).expect("widens");

        let record = open.records.remove(&5).expect("the record at 5");
        assert_eq!(
            record.ref_span(),
            16,
            "the record must have widened to 5..=20"
        );
        assert_eq!(
            record.reads_without_observation,
            vec![holed_id],
            "the widen brought the hole at 11 inside the footprint"
        );
        assert!(
            !record.folded_reads.contains_key(&holed_id),
            "and the read has no row left"
        );
        let (emitted, witness) = record.finalise();
        assert_eq!(witness.reads_without_observation, 1);
        assert_eq!(
            emitted
                .observations
                .iter()
                .map(|observation| observation.num_obs)
                .sum::<u32>(),
            1,
            "only `widener` still supports this record — `holed`'s contribution came \
             off the bucket it was in: {emitted:?}"
        );
    }

    /// **The REF bucket is never evicted, whatever its support.**
    ///
    /// `alleles[0]` is the record's own reference sequence: it is what `ref_span()`
    /// measures, what `widen` extends, and what `finalise` emits as the record's reference
    /// bytes. Production creates it with zero observations by design, and a record every
    /// read disagrees with keeps it at zero. Drop the `index == 0` guard and such a record
    /// loses its reference bytes — and with them its span — the moment any other bucket
    /// empties.
    ///
    /// Same walk as [`widened_record`] minus the spanner, which is the only read that
    /// matches the reference. The opener's re-fold still empties its pre-widen bucket, so
    /// the eviction genuinely runs rather than returning early.
    #[test]
    fn eviction_keeps_the_ref_bucket_even_with_no_observations() {
        let reference = fa(WIDEN_CONTIG);
        let (active, _ids) = admitted(
            widen_fixture_reads()
                .into_iter()
                .filter(|read| &*read.qname != "spanner")
                .collect(),
        );
        let mut open = OpenPileupRecordTable::new();
        for pos in [5u32, 7] {
            let contributors = contributors_at(&active, pos);
            process_position(&mut open, pos, 0, &contributors, &[], &active, &reference)
                .expect("the fixture walks cleanly");
        }

        let record = open.records.get(&5).expect("the record at 5");
        let buckets: Vec<_> = record
            .alleles
            .iter()
            .map(|a| {
                (
                    String::from_utf8_lossy(&a.seq).to_string(),
                    a.support.num_obs,
                )
            })
            .collect();
        assert_eq!(
            record.alleles.len(),
            3,
            "REF, the opener's and the widener's — one bucket must have been evicted, or \
             this fixture returns early and never reaches the guard it exists for: \
             {buckets:?}"
        );
        assert_eq!(
            record.alleles[0].support.num_obs, 0,
            "no read here matches the reference, or the guard is not under test: {buckets:?}"
        );
        assert_eq!(
            record.alleles[0].seq.as_slice(),
            &WIDEN_CONTIG.as_bytes()[4..20],
            "the REF bucket must still hold the record's reference bytes: {buckets:?}"
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
