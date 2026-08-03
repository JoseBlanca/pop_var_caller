//! ng step 1 — read filtering: the whole-read keep/drop prelude. It takes the
//! reads of one sample's alignment file (BAM/CRAM) and yields the subset worth
//! carrying forward, plus a running tally of what was dropped and why. Every
//! decision is per-read, locus-independent, and content-preserving — filtering
//! *selects* reads, it never rewrites them.
//!
//! Design: `doc/devel/ng/spec/read_filtering.md` (the "why"),
//! `doc/devel/ng/arch/read_filtering.md` (types & interfaces).
//!
//! Read filtering is a **port** of the production filter stack in
//! [`crate::bam::alignment_input`]: it reuses that module's pure predicates
//! (`read_exceeds_mismatch_fraction`, `cigar_is_bad`) and its `FLAG_*` /
//! `DEFAULT_*` constants, and supplies only its own driver and config. The
//! `RecordBuf` → read decode is **no longer** production's: ng needs the read to
//! carry its read group, so it owns [`AlignedRead`] and the assembly that builds
//! it ([`crate::ng::read::aligned_read`]), while still calling production's
//! `compute_adaptor_boundary` and `cigar_to_ops` rather than reproducing them.
//!
//! The whole step is here: the step-1-local types, the two-phase `verdict_*`
//! cascade, the `RecordSource` input seam, and the `ReadFilter` iterator that
//! drives them into a stream of kept reads with a running drop tally. The read
//! it filters — [`RawAlignedRead`] undecoded, [`AlignedRead`] decoded, and the
//! conversion between them — lives in [`crate::ng::read::aligned_read`].

use crate::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR, FLAG_DUPLICATE, FLAG_QC_FAIL, FLAG_SECONDARY, FLAG_SUPPLEMENTARY,
    FLAG_UNMAPPED, cigar_is_bad, cigar_ref_span, read_exceeds_mismatch_fraction,
};
use crate::ng::read::aligned_read::{AlignedRead, RawAlignedRead};
use crate::ng::read::input::read_groups::{ReadGroupResolution, RecordOwner};
use crate::ng::ref_seq::{RawRefSeq, RefSeqError};
use crate::ng::types::{BaseQual, Bp, ContigId, MapQual, MismatchFraction, ReadGroupId};
// **No `noodles_bam`, `noodles_cram` or `noodles_fasta` here, and their absence is the
// point.** This module knows what the filtering policy is; it does not know what a BAM or a
// CRAM is. It did until 2026-08-03, when the two whole-file `RecordSource` implementations
// that lived below were deleted — see the note where they were. What a *record* is left the
// same day, to `aligned_read.rs`, with the two raw types.
// **No `noodles_sam as sam` here any more, and its absence is the point.** The last thing in
// this module that knew what a SAM *header* was, was `RecordSource::header` — the contig
// probe's input — and B1 took both. What is left needs `RecordBuf` alone, to read a flag and a
// mapping quality off — and the tests do not build a header either, for the same reason.
use noodles_sam::alignment::RecordBuf;
use std::io;

/// The filtering policy: which filters are active and their thresholds. Minimal
/// by design — one field per active filter, no dormant levers (downsampling,
/// read pooling enter only when they enter the pipeline). [`Default`] is the
/// production policy the lab runs with, its thresholds the reused `DEFAULT_*`
/// constants from [`crate::bam::alignment_input`]. Mirrors the filtering subset
/// of the existing `AlignmentMergedReaderConfig`.
///
/// `Option<T>` is the "no threshold" state, never a sentinel: `None` means *no
/// minimum*, `Some(q)` means *drop below `q`* — so `Some(0)` (drop nothing, but
/// a threshold is set) stays structurally distinct from `None` (no threshold).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReadFilterConfig {
    /// `None` = no minimum; `Some(q)` = drop reads with MAPQ `< q` (filter #2).
    pub min_mapq: Option<MapQual>,
    /// `None` = no minimum; `Some(n)` = drop reads shorter than `n` bp (#7).
    pub min_read_length: Option<Bp>,
    /// Drop reads flagged QC-fail (`FLAG_QC_FAIL`) — filter #6.
    pub drop_qc_fail: bool,
    /// Drop reads flagged as PCR/optical duplicates (`FLAG_DUPLICATE`) — #1.
    pub drop_duplicate: bool,
    /// `None` = filter #8 off (no reference access at all); `Some(x)` = drop a
    /// read whose quality-clearing `M`-op mismatch fraction exceeds `x`.
    pub max_read_mismatch_fraction: Option<MismatchFraction>,
    /// BQ floor below which a mismatch does not count toward filter #8. Only
    /// meaningful when `max_read_mismatch_fraction` is `Some`.
    pub mismatch_bq_floor: BaseQual,
}

impl Default for ReadFilterConfig {
    fn default() -> Self {
        Self {
            min_mapq: Some(MapQual(DEFAULT_MIN_MAPQ)),
            // `DEFAULT_MIN_READ_LENGTH` is production's (`src/bam/`, frozen) and is
            // `u32`; ng's `Bp` is `u64` since B2 (spec §4). Widening at ng's own
            // boundary is lossless, and production does not move.
            min_read_length: Some(Bp(u64::from(DEFAULT_MIN_READ_LENGTH))),
            drop_qc_fail: true,
            drop_duplicate: true,
            // PANIC-FREE: the default fraction is a known-good in-range constant
            // (0.10, `DEFAULT_MAX_READ_MISMATCH_FRACTION`), so the checked
            // constructor cannot fail here. Guarded by the
            // `default_config_reproduces_the_production_filter_policy` test, which
            // exercises this exact path.
            max_read_mismatch_fraction: Some(
                MismatchFraction::try_new(DEFAULT_MAX_READ_MISMATCH_FRACTION)
                    .expect("DEFAULT_MAX_READ_MISMATCH_FRACTION is in [0, 1]"),
            ),
            mismatch_bq_floor: BaseQual(DEFAULT_MISMATCH_BQ_FLOOR),
        }
    }
}

/// The verdict for one read. `Keep` carries the read on; `Drop` records which
/// filter fired — the first one, per the hit-rate order — for the tally.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterVerdict {
    Keep,
    Drop(DropReason),
}

/// Which filter dropped a read. Variant names line up 1:1 with the
/// [`ReadFilterCounts`] fields, so a drop maps to exactly one counter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    Duplicate,
    LowMapq,
    Supplementary,
    Secondary,
    Unmapped,
    QcFail,
    TooShort,
    HighMismatchFraction,
    BadCigar,
}

/// A per-sample tally of the filtering pass — one counter per drop reason, plus
/// the kept count. The ng port of `FilterCounts`. Surfacing every drop is the
/// "no silent caps" discipline: a read that vanished must be accounted for. It
/// is a **running** tally — readable at any point, final once the input is
/// exhausted.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReadFilterCounts {
    pub kept: u64,
    pub duplicate: u64,
    pub low_mapq: u64,
    pub supplementary: u64,
    pub secondary: u64,
    pub unmapped: u64,
    pub qc_fail: u64,
    pub too_short: u64,
    pub high_mismatch_fraction: u64,
    pub bad_cigar: u64,
    /// Records skipped because they belong to **another sample** — a file's read
    /// groups need not all be this open's.
    ///
    /// **Not a drop, and deliberately outside every other counter here.** The
    /// rest of this struct answers "how did this read group behave?", and a read
    /// belonging to someone else says nothing about that: counting it as a drop
    /// would make a shared file look like a low-quality one. It is kept so the
    /// records a source saw still add up.
    pub other_sample: u64,
}

impl ReadFilterCounts {
    /// Add another tally into this one, counter by counter.
    ///
    /// Destructured exhaustively rather than field-by-field, so a counter added
    /// to this struct later must be routed here explicitly or this stops
    /// compiling — the same guard `record_drop`'s exhaustive `match` gives the
    /// `DropReason` mapping. A tally that silently stopped summing one reason
    /// would under-report drops without failing anything.
    pub(crate) fn add(&mut self, other: &Self) {
        let Self {
            kept,
            duplicate,
            low_mapq,
            supplementary,
            secondary,
            unmapped,
            qc_fail,
            too_short,
            high_mismatch_fraction,
            bad_cigar,
            other_sample,
        } = other;

        self.kept += kept;
        self.duplicate += duplicate;
        self.low_mapq += low_mapq;
        self.supplementary += supplementary;
        self.secondary += secondary;
        self.unmapped += unmapped;
        self.qc_fail += qc_fail;
        self.too_short += too_short;
        self.high_mismatch_fraction += high_mismatch_fraction;
        self.bad_cigar += bad_cigar;
        self.other_sample += other_sample;
    }

    /// Tally one drop against its counter. The exhaustive `match` is the guard for
    /// the documented `DropReason` ↔ `ReadFilterCounts` 1:1 mapping: adding a
    /// `DropReason` variant without a counter here is a compile error, so the two
    /// cannot silently desync (mirrors production `FilterCounts::record_drop`).
    fn record_drop(&mut self, reason: DropReason) {
        match reason {
            DropReason::Duplicate => self.duplicate += 1,
            DropReason::LowMapq => self.low_mapq += 1,
            DropReason::Supplementary => self.supplementary += 1,
            DropReason::Secondary => self.secondary += 1,
            DropReason::Unmapped => self.unmapped += 1,
            DropReason::QcFail => self.qc_fail += 1,
            DropReason::TooShort => self.too_short += 1,
            DropReason::HighMismatchFraction => self.high_mismatch_fraction += 1,
            DropReason::BadCigar => self.bad_cigar += 1,
        }
    }
}

/// Phase one of the cascade — the flag/MAPQ filters (#1–#6), decided on an
/// undecoded record's `flag` and `mapq` alone. Reference-free and decode-free:
/// `Keep` means "decode and continue to phase two", `Drop` is charged to the
/// first filter that fires. Order is identical to production's
/// `classify_pre_decode` (hit-rate-ordered: duplicate, low-MAPQ, supplementary,
/// secondary, unmapped, QC-fail).
///
/// `mapq` is already resolved: SAM's "unavailable" (`0xFF`) is mapped to
/// `MapQual(0)` by the record source (Milestone C), so a non-zero `min_mapq`
/// drops it — matching production. `flag` is the raw SAM bitfield
/// (`AlignedRead.flag`), tested against the reused `FLAG_*` constants.
fn verdict_pre_decode(flag: u16, mapq: MapQual, config: &ReadFilterConfig) -> FilterVerdict {
    // 1. Duplicate — a PCR/optical copy of another molecule (toggle).
    if config.drop_duplicate && (flag & FLAG_DUPLICATE) != 0 {
        return FilterVerdict::Drop(DropReason::Duplicate);
    }
    // 2. Low MAPQ — the aligner is unsure of the placement (toggle via threshold).
    if let Some(min) = config.min_mapq
        && mapq < min
    {
        return FilterVerdict::Drop(DropReason::LowMapq);
    }
    // 3. Supplementary — a chunk of a chimeric read (unconditional).
    if (flag & FLAG_SUPPLEMENTARY) != 0 {
        return FilterVerdict::Drop(DropReason::Supplementary);
    }
    // 4. Secondary — a duplicate projection of a primary alignment (unconditional).
    if (flag & FLAG_SECONDARY) != 0 {
        return FilterVerdict::Drop(DropReason::Secondary);
    }
    // 5. Unmapped — no placement, so no allele evidence (unconditional).
    if (flag & FLAG_UNMAPPED) != 0 {
        return FilterVerdict::Drop(DropReason::Unmapped);
    }
    // 6. QC fail — the sequencer/pipeline flagged the read (toggle).
    if config.drop_qc_fail && (flag & FLAG_QC_FAIL) != 0 {
        return FilterVerdict::Drop(DropReason::QcFail);
    }
    FilterVerdict::Keep
}

/// Phase two of the cascade — the decode-dependent filters, run on a decoded
/// [`AlignedRead`] only after it clears phase one. Evaluated **cheapest-first**:
///
/// 1. **#7 too-short** — one length compare, no reference.
/// 2. **#9 bad-CIGAR** — a pure CIGAR scan, no reference.
/// 3. **#8 high-mismatch** — a reference fetch plus a per-base walk, and only
///    when `max_read_mismatch_fraction` is `Some`.
///
/// This is a deliberate reordering relative to the spec's #7/#8/#9 *table*: it
/// puts the one reference-touching filter (#8) last, so a read dropped for
/// being too short or for a malformed CIGAR never pays the reference fetch and
/// base walk. It honours the spec's stated "cheapest, most-often-firing first"
/// principle (spec §3), charges a both-failing read to the root cause
/// (`BadCigar`) rather than the symptom (`HighMismatchFraction`), and leaves the
/// keep/drop *set* unchanged (the filters are independent — a read failing any
/// is dropped regardless of order). The mismatch fraction is measured on the
/// read's own (un-left-aligned) CIGAR: left-alignment only shifts indels across
/// equal bases, so it does not change the match/mismatch tally, and it is
/// deferred to `pileup/` anyway (spec §6).
///
/// `ref_buf` is a caller-owned scratch buffer the reference bytes are read into
/// (reused across reads, so #8 costs no per-read allocation). It is written only
/// when #8 actually runs. A [`RefSeqError`] from the fetch is **fatal to the
/// run**, propagated rather than swallowed into a drop or a keep. This includes
/// an `OutOfBounds` window running past the contig end: a validly-aligned read
/// never covers reference positions the contig does not have, so an
/// out-of-bounds fetch signals a malformed record — and the fatal model treats
/// it, like a truncated file, as corrupt input to fail loudly on rather than
/// filter around (spec §7).
fn verdict_post_decode(
    read: &AlignedRead,
    reference: &impl RawRefSeq,
    config: &ReadFilterConfig,
    ref_buf: &mut Vec<u8>,
) -> Result<FilterVerdict, RefSeqError> {
    // #7 — too short (decoded SEQ length). Cheapest: no CIGAR walk, no reference.
    if let Some(min) = config.min_read_length
        && (read.seq.len() as u64) < min.get()
    {
        return Ok(FilterVerdict::Drop(DropReason::TooShort));
    }

    // #9 — bad CIGAR (adjacent I/D, or a boundary deletion). Cheap: a pure scan
    // of the aligner's CIGAR, no reference. Run before #8 so a malformed read
    // never pays the reference fetch below.
    if cigar_is_bad(&read.cigar) {
        return Ok(FilterVerdict::Drop(DropReason::BadCigar));
    }

    // #8 — high mismatch fraction. The only reference-dependent filter, and the
    // only post-decode one that allocates work; skipped entirely when disabled.
    if let Some(max) = config.max_read_mismatch_fraction {
        let ref_span = cigar_ref_span(&read.cigar);
        // PANIC-FREE: **ids stay `u32`, coordinates are now `u64`** (spec §4 / B2).
        // `ref_id` indexes the `u32` contig table, so that conversion is the only
        // one that can fail, and a value that did not fit would be a corrupt
        // record — failing loudly is the intended response under the fatal error
        // model, not a silent `as` truncation that would fetch the wrong window and
        // mis-verdict the read. The position and span now *widen* into the RefSeq
        // boundary, which cannot fail at all: B2 moved that surface to `u64`, so the
        // `u32::try_from(read.pos)` this replaced — which could reject a legal
        // position on a >4 Gb contig — is simply gone.
        let contig = ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32"));
        // Reads raw (un-canonicalised) bytes, matching production's
        // `RawContigRefCache` path so the ported filter behaves identically.
        // A zero span yields an empty slice → `read_exceeds_mismatch_fraction`
        // has no comparable bases and keeps the read (same as production).
        reference.fetch_raw_into(contig, read.pos, u64::from(ref_span), ref_buf)?;
        if read_exceeds_mismatch_fraction(
            &read.cigar,
            &read.seq,
            &read.qual,
            ref_buf,
            config.mismatch_bq_floor.get(),
            max.get(),
        ) {
            return Ok(FilterVerdict::Drop(DropReason::HighMismatchFraction));
        }
    }

    Ok(FilterVerdict::Keep)
}

// ---------------------------------------------------------------------
// The record-source seam (input edge)
// ---------------------------------------------------------------------

/// The filter's input: fills a caller-owned buffer with the next record, reusing
/// its allocations. Modelled as a source that *fills* a buffer rather than an
/// `Iterator<Item = RawAlignedRead>` precisely so one buffer survives the whole pass
/// (a std iterator's owned `Item` cannot borrow a reused buffer — the
/// lending-iterator problem, spec §5). `Ok(true)` = filled, `Ok(false)` = end of
/// input; an `Err` is **fatal to the run**.
pub(crate) trait RecordSource {
    /// The reused buffer type — a [`RawAlignedRead`] the source refills in place. The
    /// `Default` bound is load-bearing: the [`ReadFilter`] iterator seeds *one*
    /// buffer with `Default::default()` and hands the same `&mut` to
    /// [`Self::read_next`] on every read, so the whole pass allocates one record.
    type Record: RawAlignedRead + Default;
    // **`header()` was here, and it went with the contig probe (2026-08-03, B1).** Its only
    // caller was `ReadFilter::new`, which read the `@SQ` list to know how many contigs to
    // fetch a window for. The check that replaced it compares contig tables one layer up, at
    // `AlignmentFile::cursor`, where the file's list is already in hand — so no source is
    // asked for a header any more. Note for C3, which turns this trait's methods into inherent
    // ones on `RegionRawAlignedReads`: arch §3.3 lists `header()` among them, and it should be
    // re-added only if a caller appears, because none exists today.
    /// Fill `buf` with the next record, reusing its allocations. `Ok(true)` =
    /// filled; `Ok(false)` = end of input, after which `buf` holds an unspecified
    /// (stale) record the caller must not read; `Err` is fatal to the run.
    fn read_next(&mut self, buf: &mut Self::Record) -> io::Result<bool>;

    /// How many records this source skipped because they belong to another
    /// sample.
    ///
    /// Read when the stream is taken apart and folded into the tally beside the
    /// drops — **as its own counter, never as one of them** (see
    /// [`ReadFilterCounts::other_sample`]).
    ///
    /// **Required, with no default, deliberately.** A default returning `0`
    /// would be the obvious convenience and is a trap: a source that skips
    /// records and forgets to override it reports a *wrong* number rather than
    /// an absent one, and the records vanish from the accounting with nothing
    /// failing. That is not hypothetical — it happened twice while this was
    /// being written, once to a delegating wrapper that inherited the default
    /// and once to a source that incremented a counter nothing could read. A
    /// source with nothing to skip returns `0` and says why; every other one is
    /// a compile error until it answers.
    fn other_sample_records(&self) -> u64;
}

/// The read group one filled record belongs to, or a fatal error naming it.
///
/// **The `Sole` arm returns before touching the record**, which is not merely an
/// optimisation. A file declaring one read group is read without its records'
/// `RG` tags being consulted at all — that is what lets a file re-headered
/// without rewriting its records be read, and it means a malformed tag in such a
/// file cannot end a run over a value nothing needs. It is also the universal
/// path: every record would otherwise pay a linear scan of its auxiliary fields,
/// since noodles' `Data::get` is a find, not a map lookup.
///
/// It runs where the record is filled rather than at decode because the answer is
/// needed before anything can be attributed to a library.
///
/// Errors name the read. Nothing above this point can: the failure travels up as
/// a source failure and then a filter failure, and neither knows which record was
/// being read.
pub(crate) fn resolve_read_group(
    record: &RecordBuf,
    resolution: &ReadGroupResolution,
) -> io::Result<RecordOwner> {
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;

    // Before the tag is read, not after: see the note above.
    if let ReadGroupResolution::Sole(id) = resolution {
        return Ok(RecordOwner::Mine(*id));
    }

    let tag = match record.data().get(&Tag::READ_GROUP) {
        Some(Value::String(name)) => Some(name.as_ref()),
        // A non-string `RG` is a malformed record, not an absent tag; reading it
        // as absent would report a broken file as a missing tag.
        Some(_) => {
            return Err(unreadable_record(
                record,
                "has an RG tag that is not a string",
            ));
        }
        None => None,
    };

    resolution
        .owner_of(tag)
        .map_err(|unresolved| unreadable_record(record, &unresolved.to_string()))
}

/// A fatal per-record error that names the record it is about.
///
/// Records need not be named, and an unnamed one is exactly the sort that turns
/// up in a malformed file — so it says so rather than rendering an empty pair of
/// quotes and leaving the reader to guess whether the name was blank or the
/// message was broken.
fn unreadable_record(record: &RecordBuf, problem: &str) -> io::Error {
    let name = match record.name() {
        Some(name) => format!("read '{}'", String::from_utf8_lossy(name.as_ref())),
        None => "an unnamed read".to_string(),
    };
    io::Error::new(io::ErrorKind::InvalidData, format!("{name} {problem}"))
}

// **`NoodlesRawAlignedRead` was here, and it moved to `read/aligned_read.rs` with the
// `RawAlignedRead` trait it implements (2026-08-03, `read_filtering_stages.md` §6).** A raw
// aligned read and a decoded one are one thing in two states, and the conversion between them
// already lived there. What is left here is the keep-or-drop rules.

// **`BamRecordSource` and `CramRecordSource` were here, and they are gone (2026-08-03).**
//
// They were whole-file `RecordSource` implementations: open a BAM or a CRAM, hand over every
// record in file order, let the filter decide. They were built with this seam
// (`read_filtering.md`), before there was any other way to read a file.
//
// **Two things retired them.** They had no caller — only this module's own tests drove them —
// and, more to the point, *a filter module has no business opening files*. Finding and
// unpacking records is `read/input/aligned_reads_reader/`'s job, and since the alignment cursor
// (`spec/alignment_cursor.md`) there is exactly one shape for it: an `AlignedReadsReader` positions,
// `RegionRawAlignedReads` narrows, and this filter consumes what they hand over. Keeping a second,
// unused file reader here left ng with two ways to read a BAM and only one of them reachable.
//
// **What became of what they proved.** The four tests whose subject was the sources themselves
// went with them — raw records out of a real file, buffer reuse without leaking the previous
// record, and CRAM agreement across container boundaries are all asserted on the cursor path
// (`aligned_reads_reader/bam.rs`, `aligned_reads_reader/cram.rs`, and the run-of-regions oracles in
// `open_bam.rs`). The two that tested *this filter* against a hand-counted drop tally over a
// real BAM and a real CRAM were kept and re-pointed at an `AlignedReadsReader`, which is where they
// belong: they were never about the source.
//
// **A whole-file pass is therefore not available in ng, and that is now explicit.** Point a
// cursor at a whole chromosome — which is what the linear-scan oracles do. If a genuine
// whole-file need appears (a coverage histogram, an unindexed input), it arrives as a
// `AlignedReadsReader` arm, beside the others, not as a type in the filter.

// ---------------------------------------------------------------------
// The ReadFilter iterator (the driver)
// ---------------------------------------------------------------------

/// A fatal, run-level error from a [`ReadFilter`] pass. It is yielded **in the
/// iterator's item stream** — a fatal condition (a failed record read, a failed
/// decode, or a reference-fetch failure) makes `next()` return `Some(Err(..))`
/// once and then `None` — so the caller cannot mistake it for a clean end of
/// input: `let read = read?;` propagates it. It is never folded into a per-read
/// drop or a silent EOF.
#[derive(Debug, thiserror::Error)]
pub(crate) enum ReadFilterError {
    /// The record source failed to read the next record (e.g. a truncated file).
    #[error("reading the next alignment record failed")]
    Source(#[source] io::Error),
    /// A record that cleared the pre-decode gate failed to decode — a corrupt
    /// record (the unmapped flag clear yet no position).
    #[error("decoding an alignment record failed")]
    Decode(#[source] io::Error),
    /// Filter #8's reference fetch failed — corrupt input, or a read reaching past a contig's
    /// end.
    ///
    /// It can also mean the caller's contig guarantee did not hold, since
    /// `with_validated_contigs` takes that on trust. Since B1 that is much narrower than it
    /// was: `AlignmentFile::cursor` compares the two contig tables *and* proves the accessor
    /// can serve this cursor's own contig before building anything, so what is left to fail
    /// here is a contig the cursor is **not** for — reachable only through a caller that builds
    /// a filter itself.
    #[error("reference access failed during filtering")]
    Reference(#[source] RefSeqError),
}

/// Filters one sample's reads, lazily, as an
/// `Iterator<Item = Result<AlignedRead, ReadFilterError>>`. Each `next()` reads
/// the next record into the single reused buffer, runs the pre-decode cascade
/// (#1–#6) on its flag/MAPQ, decodes only on survival, runs the decode-dependent
/// cascade (#7, #9, #8), tallies every drop, and returns the first read that
/// passes as `Ok(read)`. A read dropped pre-decode builds no `AlignedRead` at all.
///
/// **The item is a `Result`** (unlike the spec's original `Item = AlignedRead`,
/// revised so a fatal error cannot be silently lost): a fatal condition yields
/// `Some(Err(_))` once and then `None`, so `for read in &mut filter { let read =
/// read?; … }` surfaces it with no separate bookkeeping. This matches the noodles
/// readers this sits on (`records()` is `Item = io::Result<_>`). The iterator is
/// fused — after a clean `None` or a fatal `Err`, further `next()`s return `None`.
///
/// `counts()` is a **running** tally, readable at any point and final once the
/// iterator is exhausted (iterate by `&mut` so the filter — and its counts —
/// outlives the loop).
pub(crate) struct ReadFilter<S: RecordSource, R> {
    source: S,
    /// The single record buffer reused across every read.
    record_buf: S::Record,
    /// Raw reference bytes for filter #8; see `ref_seq.md`.
    reference: R,
    config: ReadFilterConfig,
    /// One tally per read group met, in first-seen order.
    ///
    /// A `Vec` scanned linearly rather than a map: a file declares a handful of
    /// read groups, and the scan is a few integer compares on a path that has
    /// just decoded a record. Per read group rather than per *file* because a
    /// drop rate is a read group's property — one bad run shows up as one read
    /// group with an anomalous MAPQ or mismatch rate, and a per-file tally over
    /// a file that holds several would average that away.
    counts: Vec<ReadGroupCounts>,
    /// Reused scratch for #8's reference fetch (touched only when #8 runs).
    ref_buf: Vec<u8>,
    /// Whether this filter is still yielding, and if not, **why it stopped**.
    /// Makes the iterator fused: anything but `Running` returns `None`.
    state: FilterState,
}

/// Whether a [`ReadFilter`] is still yielding reads, and if not, why it stopped.
///
/// **The two ways of stopping have to be told apart**, which a single `done` flag could not.
/// A long-lived filter — one the alignment cursor keeps for a whole chromosome and repoints
/// at each region — reaches `EndOfInput` at the end of *every* region, because the region
/// narrowing below it reports a region boundary the only way a `RecordSource` can, as
/// `Ok(false)`. Repositioning must be able to undo that. It must **not** undo `Failed`: a
/// filter that met a corrupt record or an I/O fault yields the error once and then stays
/// quiet, and resurrecting it would start the next region reading from a file already known
/// to be broken (owner, 2026-08-02).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FilterState {
    /// Reads are still coming.
    Running,
    /// The source reported a clean end of input. Undone by
    /// [`ReadFilter::restart_after_end_of_input`].
    EndOfInput,
    /// A fatal error was yielded. **Permanent** — the only way on is a new filter.
    Failed,
}

/// One read group's step-1 tally, keyed by the read group it belongs to.
///
/// `None` keys the records whose source never stamped a read group — which
/// [`RawAlignedRead::decode`] refuses, so they are counted apart rather than charged
/// to an arbitrary group.
pub type ReadGroupCounts = (Option<ReadGroupId>, ReadFilterCounts);

/// The two buffers a [`ReadFilter`] pass reuses: the record buffer it refills
/// per read, and the scratch filter #8 fetches reference bases into.
///
/// **A vestige of the per-region query, and it should collapse.** They were grouped into one
/// value so that path could keep them beside each pooled reader and lend them to a fresh
/// per-region filter, which is what kept ~10⁶ region queries from allocating two buffers
/// each. Milestone F deleted that path along with the method that handed them back, so the
/// only value anything passes here now is `Default::default()`, from the two remaining call
/// sites: a
/// filter lives as long as the cursor it belongs to and reuses its own buffers for every
/// region. Folding the parameter away is a `filtering.rs` cleanup, recorded rather than done
/// inside a step whose subject is the read path.
#[derive(Debug, Default)]
pub(crate) struct ReadFilterBuffers<Rec> {
    /// The single record buffer refilled by every [`RecordSource::read_next`].
    pub(crate) record_buf: Rec,
    /// Reused scratch for filter #8's reference fetch.
    pub(crate) ref_buf: Vec<u8>,
}

impl<S: RecordSource, R: RawRefSeq> ReadFilter<S, R> {
    /// Seed a filter over a source whose contigs are **already known to match the
    /// reference**, taking its two reused buffers from the caller.
    ///
    /// # There used to be a second constructor, and its body was a loop
    ///
    /// `ReadFilter::new` fetched a zero-length window on every contig of the source's `@SQ`
    /// list to prove each one resolved in `reference` — ~2,580 opens on GRCh38, once per
    /// cursor. It was deleted on 2026-08-03, when `AlignmentFile::cursor` started **comparing
    /// the two contig tables** instead (`read_filtering_stages.md` §9 Q2). The comparison is
    /// microseconds rather than thousands of file opens, and it proves strictly more: the loop
    /// only showed that each contig *resolves*, so an accessor whose contigs had the right
    /// names and the wrong lengths sailed through it.
    ///
    /// # What the caller must have proved
    ///
    /// That the source's `@SQ` list is *equal to* the reference's contig table — name, length
    /// and **order** included — which is what both `AlignmentFile::open`'s gate and
    /// `AlignmentFile::cursor`'s check establish.
    ///
    /// Handing over a source whose contigs were never checked has two failure modes, and **the
    /// quiet one is worse.** A contig missing from the reference surfaces as a mid-stream
    /// [`ReadFilterError::Reference`] — loud, and recoverable. But a *permuted* `@SQ` list
    /// resolves on every fetch, so filter #8 compares each read against the wrong contig's
    /// bases and drops and keeps the wrong reads with no error at all. That is why the
    /// precondition is equality and not resolvability, and why the deleted loop was never
    /// enough on its own.
    ///
    /// The guarantee is documented rather than encoded in a typestate: the caller that needs
    /// this is the one that just ran the check, and a marker type would buy proof only for
    /// callers who already have it.
    ///
    /// The `buffers` parameter is a vestige of the per-region query that once lent them
    /// (see [`ReadFilterBuffers`]); every caller now passes `Default::default()`.
    pub(crate) fn with_validated_contigs(
        source: S,
        reference: R,
        config: ReadFilterConfig,
        buffers: ReadFilterBuffers<S::Record>,
    ) -> Self {
        Self {
            source,
            record_buf: buffers.record_buf,
            reference,
            config,
            counts: Vec::new(),
            ref_buf: buffers.ref_buf,
            state: FilterState::Running,
        }
    }

    /// Undo a **clean** end of input, so a repositioned source can be read again.
    ///
    /// The other half of [`source_mut`](Self::source_mut), and the half that is easy to miss:
    /// pointing the source at a new region is not enough, because reaching the end of the
    /// *previous* region already fused this filter. A region boundary arrives here as an
    /// ordinary end of input — `Ok(false)` is the only way a [`RecordSource`] can report one
    /// — so without this the first region a cursor drains silences it for the whole
    /// chromosome. Measured before it existed: two regions through one filter gave 2 reads,
    /// then 0.
    ///
    /// **A failed filter is not restarted**, and that is the point of separating the two
    /// states. A corrupt record or an I/O fault is yielded once and then the filter stays
    /// quiet; carrying on would mean reading the next region out of a file already known to
    /// be broken. Answers whether it restarted anything, so a caller that cares can tell
    /// "there was nothing to undo" from "this filter is finished for good".
    pub(crate) fn restart_after_end_of_input(&mut self) -> bool {
        if self.state == FilterState::EndOfInput {
            self.state = FilterState::Running;
            return true;
        }
        false
    }

    /// Whether this filter stopped because of a **fatal error**, as opposed to running out of
    /// input or still running.
    ///
    /// A caller that keeps one filter across many regions has to be able to tell the two
    /// stops apart: an end of input is this region's business, while a failure is the whole
    /// file's and no later region can be served from it. Without this the distinction exists
    /// inside the filter and nowhere a caller can see it, and the caller reports later
    /// regions as **empty** rather than as impossible.
    pub(crate) fn has_failed(&self) -> bool {
        self.state == FilterState::Failed
    }

    /// The source underneath, to reposition it.
    ///
    /// **The one thing a long-lived filter needs that a per-query one did not.** Until the
    /// alignment cursor, a filter was built per region and thrown away, so the only way to
    /// reach its source was `into_parts` — which consumed the filter, and therefore *forced*
    /// a new one per region. (That method went with the per-region path at Milestone F, which
    /// is why this names it rather than linking it.) A cursor keeps one filter for as long as it
    /// reads a chromosome and points the layer below at each new region through it
    /// (`arch/alignment_cursor.md` §2.3).
    ///
    /// **This accessor alone is not enough**, and the clean case is the one that bites: a
    /// region's end reaches the filter as an ordinary end of input, so the first region a
    /// cursor drains would fuse it for the rest of the chromosome — measured at 2 reads, then
    /// 0. Pair it with [`restart_after_end_of_input`](Self::restart_after_end_of_input),
    /// which is what undoes that and refuses to undo a failure.
    ///
    /// The tally is deliberately **not** reset either: `counts()` is then a whole-cursor
    /// total rather than whichever region happened to be last.
    pub(crate) fn source_mut(&mut self) -> &mut S {
        &mut self.source
    }

    /// The reference this filter reads for the mismatch-fraction check.
    ///
    /// **Shared, not mutable, and that is the point.** The one thing a caller needs it for is
    /// *releasing the bases the walk has gone past*, which is a `&self` operation
    /// ([`EvictableRefSeq`](crate::ng::ref_seq::EvictableRefSeq)). A filter that lives as long
    /// as a cursor reads the reference once per surviving read, forward, and never releases
    /// any of it on its own — on a densely covered chromosome that is one byte held for every
    /// base walked, about 250 MB on human chromosome 1.
    pub(crate) fn reference(&self) -> &R {
        &self.reference
    }

    // `into_parts` is gone. It handed the source, the two buffers and the final
    // tally back to whoever lent them, which was the per-region query: a pooled
    // reader reclaimed them as each region stream ended. Milestone F deleted that
    // path, and a filter now lives as long as the cursor it belongs to — so its
    // buffers are reused by construction and its tally is read in place, through
    // [`counts`](Self::counts), rather than folded somewhere before the filter dies.

    /// The tally of the read group the buffered record belongs to, created on
    /// first sight.
    ///
    /// Keyed on the buffer rather than passed in, because the pre-decode filters
    /// drop a record before anything owns it — the buffer is the only thing that
    /// knows, and it knows because its source resolved the read group before
    /// handing it over.
    fn tally_for_current_record(&mut self) -> &mut ReadFilterCounts {
        let read_group = self.record_buf.read_group();
        let index = match self.counts.iter().position(|(id, _)| *id == read_group) {
            Some(index) => index,
            None => {
                self.counts.push((read_group, ReadFilterCounts::default()));
                self.counts.len() - 1
            }
        };
        &mut self.counts[index].1
    }

    /// One tally per read group this pass met, in first-seen order, plus what
    /// the source skipped as another sample's.
    ///
    /// **Never summed here.** A drop rate is a read group's property: one bad
    /// run is one read group with an anomalous rate, and adding them up erases
    /// exactly that. A caller wanting a total adds them itself.
    ///
    /// The skipped-record count belongs to no single read group — a foreign
    /// record is never yielded to this filter, and its own group is not one of
    /// these — so it is reported once, against the first entry.
    pub fn counts(&self) -> Vec<ReadGroupCounts> {
        let mut counts = self.counts.clone();
        match counts.first_mut() {
            Some((_, first)) => first.other_sample = self.source.other_sample_records(),
            // Every record was another sample's, so the filter met none at all.
            None => counts.push((
                None,
                ReadFilterCounts {
                    other_sample: self.source.other_sample_records(),
                    ..ReadFilterCounts::default()
                },
            )),
        }
        counts
    }

    /// Mark the iterator finished and yield a fatal error. Shared by the three
    /// fatal arms of `next()`.
    fn fail(&mut self, error: ReadFilterError) -> Option<Result<AlignedRead, ReadFilterError>> {
        self.state = FilterState::Failed;
        Some(Err(error))
    }
}

impl<S: RecordSource, R: RawRefSeq> Iterator for ReadFilter<S, R> {
    type Item = Result<AlignedRead, ReadFilterError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Fused: once a clean EOF or a fatal error has been reported, stay stopped.
        if self.state != FilterState::Running {
            return None;
        }
        loop {
            match self.source.read_next(&mut self.record_buf) {
                Ok(true) => {}
                Ok(false) => {
                    self.state = FilterState::EndOfInput;
                    return None; // clean end of input
                }
                Err(error) => return self.fail(ReadFilterError::Source(error)),
            }

            // Phase one — flag/MAPQ, before decode. Exhaustive so a future
            // `FilterVerdict` variant cannot silently fall through to decode.
            match verdict_pre_decode(self.record_buf.flag(), self.record_buf.mapq(), &self.config) {
                FilterVerdict::Keep => {}
                FilterVerdict::Drop(reason) => {
                    self.tally_for_current_record().record_drop(reason);
                    continue;
                }
            }

            // Decode only the pre-decode survivors.
            let read = match self.record_buf.decode() {
                Ok(read) => read,
                Err(error) => return self.fail(ReadFilterError::Decode(error)),
            };

            // Phase two — length / CIGAR / mismatch, on the decoded read.
            match verdict_post_decode(&read, &self.reference, &self.config, &mut self.ref_buf) {
                Ok(FilterVerdict::Keep) => {
                    self.tally_for_current_record().kept += 1;
                    return Some(Ok(read));
                }
                Ok(FilterVerdict::Drop(reason)) => {
                    self.tally_for_current_record().record_drop(reason);
                    continue;
                }
                Err(error) => return self.fail(ReadFilterError::Reference(error)),
            }
        }
    }
}

impl<S: RecordSource, R: RawRefSeq> std::iter::FusedIterator for ReadFilter<S, R> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FLAG_PAIRED;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::pileup::walker::CigarOp;

    #[test]
    fn default_config_reproduces_the_production_filter_policy() {
        let config = ReadFilterConfig::default();
        assert_eq!(config.min_mapq, Some(MapQual(DEFAULT_MIN_MAPQ)));
        assert_eq!(
            config.min_read_length,
            Some(Bp(u64::from(DEFAULT_MIN_READ_LENGTH)))
        );
        assert!(config.drop_qc_fail);
        assert!(config.drop_duplicate);
        assert_eq!(
            config.max_read_mismatch_fraction.map(MismatchFraction::get),
            Some(DEFAULT_MAX_READ_MISMATCH_FRACTION)
        );
        assert_eq!(
            config.mismatch_bq_floor,
            BaseQual(DEFAULT_MISMATCH_BQ_FLOOR)
        );
    }

    #[test]
    fn counts_default_is_all_zero() {
        // Explicit all-zero literal (no `..`): pins every counter to 0 and forces
        // this test to be revisited if a counter field is ever added.
        assert_eq!(
            ReadFilterCounts::default(),
            ReadFilterCounts {
                kept: 0,
                duplicate: 0,
                low_mapq: 0,
                supplementary: 0,
                secondary: 0,
                unmapped: 0,
                qc_fail: 0,
                too_short: 0,
                high_mismatch_fraction: 0,
                bad_cigar: 0,
                other_sample: 0,
            }
        );
    }

    // ----- verdict_pre_decode (#1–#6) -------------------------------------

    /// A record's flag/MAPQ pair; the pre-decode cascade needs nothing else.
    fn pre(flag: u16, mapq: u8, config: &ReadFilterConfig) -> FilterVerdict {
        verdict_pre_decode(flag, MapQual(mapq), config)
    }

    #[test]
    fn pre_decode_keeps_a_clean_primary_read() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(pre(0, 60, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn low_mapq_boundary_keeps_at_threshold_drops_one_below() {
        let cfg = ReadFilterConfig::default(); // min_mapq = Some(20)
        assert_eq!(pre(0, 20, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(0, 19, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
        // Unavailable MAPQ arrives as 0 (resolved by the record source) → dropped.
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Drop(DropReason::LowMapq));
    }

    #[test]
    fn no_mapq_minimum_keeps_any_quality() {
        let cfg = ReadFilterConfig {
            min_mapq: None,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(0, 0, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn each_flag_bit_drops_to_its_own_bucket() {
        let cfg = ReadFilterConfig::default();
        assert_eq!(
            pre(FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    #[test]
    fn duplicate_and_qc_fail_toggles_off_keep_those_reads() {
        let cfg = ReadFilterConfig {
            drop_duplicate: false,
            drop_qc_fail: false,
            ..ReadFilterConfig::default()
        };
        assert_eq!(pre(FLAG_DUPLICATE, 60, &cfg), FilterVerdict::Keep);
        assert_eq!(pre(FLAG_QC_FAIL, 60, &cfg), FilterVerdict::Keep);
        // Supplementary/secondary/unmapped have no toggle — still dropped.
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
    }

    #[test]
    fn pre_decode_attribution_charges_the_first_firing_filter() {
        let cfg = ReadFilterConfig::default();
        // Duplicate (1) + unmapped (5) both set → charged to duplicate (earlier).
        assert_eq!(
            pre(FLAG_DUPLICATE | FLAG_UNMAPPED, 60, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        // Supplementary (3) + secondary (4) → supplementary (earlier).
        assert_eq!(
            pre(FLAG_SUPPLEMENTARY | FLAG_SECONDARY, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        // Low MAPQ (2) beats unmapped (5) when both apply.
        assert_eq!(
            pre(FLAG_UNMAPPED, 5, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
    }

    #[test]
    fn pre_decode_charges_filters_in_full_cascade_order() {
        // Every filter would fire; peel them off one at a time and confirm the
        // read is charged to each in the exact cascade order
        // (duplicate → low-MAPQ → supplementary → secondary → unmapped → QC-fail).
        let cfg = ReadFilterConfig::default();
        let all_flags =
            FLAG_DUPLICATE | FLAG_SUPPLEMENTARY | FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL;
        // mapq 0 also fails #2 throughout, so every stage after #1 has low-MAPQ
        // waiting behind it — proving each earlier filter really wins.
        assert_eq!(
            pre(all_flags, 0, &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 0, &cfg),
            FilterVerdict::Drop(DropReason::LowMapq)
        );
        // From here MAPQ is fine (60), so the next flag in order wins.
        assert_eq!(
            pre(all_flags & !FLAG_DUPLICATE, 60, &cfg),
            FilterVerdict::Drop(DropReason::Supplementary)
        );
        assert_eq!(
            pre(FLAG_SECONDARY | FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Secondary)
        );
        assert_eq!(
            pre(FLAG_UNMAPPED | FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::Unmapped)
        );
        assert_eq!(
            pre(FLAG_QC_FAIL, 60, &cfg),
            FilterVerdict::Drop(DropReason::QcFail)
        );
    }

    // ----- verdict_post_decode (#7, #9, #8) -------------------------------

    /// Build a mapped read at contig 0, pos 1, with the given decoded sequence,
    /// per-base qualities, and CIGAR. Flag/MAPQ are already-passed values.
    fn mapped(seq: &[u8], qual: &[u8], cigar: Vec<CigarOp>) -> AlignedRead {
        AlignedRead {
            qname: b"read".to_vec(),
            flag: FLAG_PAIRED,
            ref_id: 0,
            pos: 1,
            mapq: 60,
            cigar,
            seq: seq.to_vec(),
            qual: qual.to_vec(),
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    /// A single contig of `n` adenines — the reference for the post-decode
    /// tests, where a `Match` read of `T`s mismatches every aligned base.
    fn poly_a_ref(n: usize) -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![vec![b'A'; n]])
    }

    /// Post-decode config isolating one filter at a time: #7 off, toggles off,
    /// #8 set to `max` (`None` disables it).
    fn post_config(max: Option<f32>) -> ReadFilterConfig {
        ReadFilterConfig {
            min_mapq: None,
            min_read_length: None,
            drop_qc_fail: false,
            drop_duplicate: false,
            max_read_mismatch_fraction: max.map(|x| MismatchFraction::try_new(x).unwrap()),
            mismatch_bq_floor: BaseQual(0),
        }
    }

    fn post(
        read: &AlignedRead,
        reference: &impl RawRefSeq,
        config: &ReadFilterConfig,
    ) -> FilterVerdict {
        let mut buf = Vec::new();
        verdict_post_decode(read, reference, config, &mut buf).unwrap()
    }

    #[test]
    fn too_short_boundary_keeps_at_threshold_drops_one_below() {
        let reference = poly_a_ref(40);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None) // #8 already disabled by post_config(None)
        };
        let at = mapped(&[b'A'; 30], &[40; 30], vec![CigarOp::Match(30)]);
        let below = mapped(&[b'A'; 29], &[40; 29], vec![CigarOp::Match(29)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&below, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn bad_cigar_drops_the_two_ill_formed_shapes() {
        let reference = poly_a_ref(40);
        let cfg = post_config(None);
        // Adjacent insertion/deletion pair.
        let adjacent_indel = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![
                CigarOp::Match(4),
                CigarOp::Insertion(1),
                CigarOp::Deletion(1),
                CigarOp::Match(4),
            ],
        );
        // Leading deletion (boundary deletion).
        let boundary_deletion = mapped(
            b"AAAAAAAA",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        assert_eq!(
            post(&adjacent_indel, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        assert_eq!(
            post(&boundary_deletion, &reference, &cfg),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
    }

    #[test]
    fn high_mismatch_boundary_keeps_at_threshold_drops_above() {
        // ref = 10×A; a Match(10) read with k mismatches has fraction k/10.
        let reference = poly_a_ref(10);
        let cfg = post_config(Some(0.10));
        // 1/10 = 0.10, not > 0.10 → kept (boundary is exclusive).
        let at = mapped(b"TAAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        // 2/10 = 0.20 > 0.10 → dropped.
        let above = mapped(b"TTAAAAAAAA", &[40; 10], vec![CigarOp::Match(10)]);
        assert_eq!(post(&at, &reference, &cfg), FilterVerdict::Keep);
        assert_eq!(
            post(&above, &reference, &cfg),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn low_quality_mismatches_do_not_count_toward_the_fraction() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            mismatch_bq_floor: BaseQual(10),
            ..post_config(Some(0.0))
        };
        // Two mismatches, both below the BQ floor → neither counts → kept even
        // at a zero threshold.
        let read = mapped(
            b"TTAAAAAAAA",
            &[5, 5, 40, 40, 40, 40, 40, 40, 40, 40],
            vec![CigarOp::Match(10)],
        );
        assert_eq!(post(&read, &reference, &cfg), FilterVerdict::Keep);
    }

    #[test]
    fn mismatch_filter_disabled_makes_no_reference_access() {
        // Empty reference: any fetch errors with UnknownContig. With #8 disabled
        // the fetch never happens, so a high-mismatch read is kept.
        let empty = InMemoryRefSeq::from_contigs(Vec::new());
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);

        let disabled = post_config(None);
        assert_eq!(post(&read, &empty, &disabled), FilterVerdict::Keep);

        // Enabling #8 on the same empty reference proves a fetch is attempted.
        let enabled = post_config(Some(0.10));
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &empty, &enabled, &mut buf),
            Err(RefSeqError::UnknownContig(_))
        ));
    }

    #[test]
    fn high_mismatch_fetch_past_contig_end_is_fatal() {
        // A read whose reference span runs past the contig end (Match(10) at
        // pos 1 on a 5-base contig). Under the fatal error model the OutOfBounds
        // fetch propagates as `Err`, not a per-read drop or keep — a
        // validly-aligned read cannot cover positions the contig lacks, so this
        // signals a malformed record.
        let reference = poly_a_ref(5);
        let read = mapped(b"TTTTTTTTTT", &[40; 10], vec![CigarOp::Match(10)]);
        let mut buf = Vec::new();
        assert!(matches!(
            verdict_post_decode(&read, &reference, &post_config(Some(0.10)), &mut buf),
            Err(RefSeqError::OutOfBounds { .. })
        ));
    }

    #[test]
    fn bad_cigar_is_charged_before_high_mismatch() {
        // A read that is BOTH a boundary deletion (#9) and, on its M positions,
        // fully mismatched against the reference (#8 at a zero threshold).
        let reference = poly_a_ref(10);
        let both = mapped(
            b"TTTTTTTT",
            &[40; 8],
            vec![CigarOp::Deletion(2), CigarOp::Match(8)],
        );
        // #9 fires first → charged BadCigar, and #8's reference walk is skipped.
        assert_eq!(
            post(&both, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::BadCigar)
        );
        // The same sequence with a well-formed CIGAR does reach #8 and drops there,
        // proving the read would have failed #8 too — attribution is the only
        // difference the ordering makes.
        let good_cigar = mapped(b"TTTTTTTT", &[40; 8], vec![CigarOp::Match(8)]);
        assert_eq!(
            post(&good_cigar, &reference, &post_config(Some(0.0))),
            FilterVerdict::Drop(DropReason::HighMismatchFraction)
        );
    }

    #[test]
    fn too_short_is_charged_before_bad_cigar() {
        let reference = poly_a_ref(10);
        let cfg = ReadFilterConfig {
            min_read_length: Some(Bp(30)),
            ..post_config(None)
        };
        // Short AND a boundary deletion → charged TooShort (#7 before #9).
        let read = mapped(
            b"AAAAA",
            &[40; 5],
            vec![CigarOp::Deletion(2), CigarOp::Match(5)],
        );
        assert_eq!(
            post(&read, &reference, &cfg),
            FilterVerdict::Drop(DropReason::TooShort)
        );
    }

    #[test]
    fn zero_reference_span_read_is_kept_by_the_mismatch_filter() {
        // An all-soft-clip read has ref_span 0 → empty slice → no comparable
        // bases → kept (matches production's skip-on-no-span behaviour).
        let reference = poly_a_ref(10);
        let read = mapped(b"TTTT", &[40; 4], vec![CigarOp::SoftClip(4)]);
        assert_eq!(
            post(&read, &reference, &post_config(Some(0.0))),
            FilterVerdict::Keep
        );
    }

    // ----- the record-source seam (C) -------------------------------------

    /// A trivial in-memory `RawAlignedRead`: it carries a whole `AlignedRead` and reads
    /// its `flag`/`mapq` back off it, so the cascade can be driven with no BAM.
    /// `decode_fails` makes `decode` return an error (to exercise the fatal
    /// decode-error path).
    #[derive(Clone)]
    struct FakeRecord {
        read: AlignedRead,
        decode_fails: bool,
    }

    impl Default for FakeRecord {
        fn default() -> Self {
            Self {
                read: mapped(b"", &[], Vec::new()),
                decode_fails: false,
            }
        }
    }

    /// The single read group's tally, for a fixture that has exactly one.
    ///
    /// Every fixture here declares one read group (or none at all, for the fake
    /// sources), so a pass produces one entry. Asserting that rather than taking
    /// the first is the point: a fixture that quietly grew a second read group
    /// would otherwise have half its drops silently ignored by these tests.
    fn only_tally<S: RecordSource, R: RawRefSeq>(filter: &ReadFilter<S, R>) -> ReadFilterCounts {
        let counts = filter.counts();
        assert_eq!(
            counts.len(),
            1,
            "this fixture is expected to meet exactly one read group"
        );
        counts[0].1.clone()
    }

    impl RawAlignedRead for FakeRecord {
        /// Unstamped: this fake has no source that resolves read groups, so its
        /// drops are charged to no read group — which is what `None` keys.
        fn read_group(&self) -> Option<ReadGroupId> {
            None
        }

        fn flag(&self) -> u16 {
            self.read.flag
        }
        fn mapq(&self) -> MapQual {
            MapQual(self.read.mapq)
        }
        fn decode(&self) -> io::Result<AlignedRead> {
            if self.decode_fails {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "fake decode failure",
                ));
            }
            Ok(self.read.clone())
        }
    }

    /// A `RecordSource` yielding a fixed slice of `FakeRecord`s, refilling the
    /// caller's buffer each call.
    struct FakeSource {
        records: Vec<FakeRecord>,
        next_index: usize,
    }

    impl FakeSource {
        /// **No header.** These doubles carried one only to answer
        /// `RecordSource::header`, which went with the contig probe at B1 — nothing asks a
        /// source for a header any more.
        fn new(records: Vec<FakeRecord>) -> Self {
            Self {
                records,
                next_index: 0,
            }
        }

        /// Send the source back to its first record — what a real source does when it is
        /// pointed at a new region, in the smallest form that is observable from outside.
        fn rewind(&mut self) {
            self.next_index = 0;
        }

        /// How many of the scripted records this source has handed over — **not** a
        /// reference coordinate, which is what `pos` means everywhere else in this file.
        /// Lets a test say *where* the source was repositioned to rather than only that
        /// something changed.
        fn records_consumed(&self) -> usize {
            self.next_index
        }
    }

    impl RecordSource for FakeSource {
        type Record = FakeRecord;
        /// Nothing to skip: this fake carries no read-group resolution, so it
        /// can never meet another sample's record.
        fn other_sample_records(&self) -> u64 {
            0
        }

        fn read_next(&mut self, buf: &mut FakeRecord) -> io::Result<bool> {
            match self.records.get(self.next_index) {
                Some(record) => {
                    *buf = record.clone();
                    self.next_index += 1;
                    Ok(true)
                }
                None => Ok(false),
            }
        }
    }

    /// A [`fake`] record under a name, so two records in one script can be told apart.
    fn named_fake(qname: &[u8], flag: u16, mapq: u8) -> FakeRecord {
        let mut record = fake(flag, mapq);
        record.read.qname = qname.to_vec();
        record
    }

    /// A 30 bp all-`A` mapped read carrying the given flag/MAPQ — long enough to
    /// clear the default `min_read_length`, and matching a poly-A reference.
    fn fake(flag: u16, mapq: u8) -> FakeRecord {
        let mut read = mapped(&[b'A'; 30], &[30u8; 30], vec![CigarOp::Match(30)]);
        read.flag = flag;
        read.mapq = mapq;
        FakeRecord {
            read,
            decode_fails: false,
        }
    }

    #[test]
    fn fake_source_drives_the_seam() {
        // Mimic what the Milestone-D iterator will do: read_next → pre-decode →
        // decode survivors → post-decode, over the fake source.
        let reference = poly_a_ref(30);
        let cfg = ReadFilterConfig::default();
        let mut source = FakeSource::new(vec![fake(FLAG_DUPLICATE, 60), fake(FLAG_PAIRED, 60)]);
        let mut buf = FakeRecord::default();
        let mut ref_buf = Vec::new();

        // Record 1 — duplicate: dropped pre-decode, never decoded.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(buf.flag(), FLAG_DUPLICATE);
        assert_eq!(buf.mapq(), MapQual(60));
        assert_eq!(
            verdict_pre_decode(buf.flag(), buf.mapq(), &cfg),
            FilterVerdict::Drop(DropReason::Duplicate)
        );

        // Record 2 — clean: pre-decode keep → decode → post-decode keep.
        assert!(source.read_next(&mut buf).unwrap());
        assert_eq!(
            verdict_pre_decode(buf.flag(), buf.mapq(), &cfg),
            FilterVerdict::Keep
        );
        let read = buf.decode().unwrap();
        assert_eq!(
            verdict_post_decode(&read, &reference, &cfg, &mut ref_buf).unwrap(),
            FilterVerdict::Keep
        );

        // End of input.
        assert!(!source.read_next(&mut buf).unwrap());
    }

    // -----------------------------------------------------------------
    // Reaching the source through the filter (alignment cursor, A4)
    // -----------------------------------------------------------------

    /// **The accessor the cursor is built on.** A filter used to be made per region and
    /// consumed to get its source back, which is what forced a new filter per region in the
    /// first place. A cursor keeps one filter for a whole chromosome and points the layer
    /// below at each new region *through* it, so the reach has to be by `&mut`.
    #[test]
    fn repositioning_through_the_filter_reaches_the_source() {
        let mut filter = filter_over(
            FakeSource::new(vec![
                named_fake(b"first", FLAG_PAIRED, 60),
                named_fake(b"second", FLAG_PAIRED, 60),
            ]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert_eq!(filter.source_mut().records_consumed(), 0);
        let first = filter.next().expect("a read").expect("it passes");
        assert_eq!(
            filter.source_mut().records_consumed(),
            1,
            "the filter consumed one record"
        );

        // The reposition, through the filter, without consuming it.
        filter.source_mut().rewind();
        assert_eq!(filter.source_mut().records_consumed(), 0);

        // Named records, because `fake` builds every one at the same position with the same
        // flags: two byte-identical reads made this assertion unable to fail.
        let again = filter.next().expect("a read").expect("it passes");
        assert_eq!(first.qname, b"first");
        assert_eq!(
            again.qname, b"first",
            "after rewinding the source the filter re-read the record after the first, not \
             the first",
        );
    }

    /// **The half of the contract the cursor will have to work around, pinned here so it is
    /// discovered now rather than at the first cursor.**
    ///
    /// `ReadFilter` is a *fused* iterator: a clean end of input sets `done`, and `next()`
    /// returns `None` from then on. `source_mut` reaches the source but does not reset that
    /// flag — so a caller that drains a region and then repositions gets nothing, however
    /// well the source was repositioned. Arch §2.3 describes `move_to_region` as
    /// `self.filter.source_mut().move_to(region)` and does not mention it.
    ///
    /// Whether the cursor drains its regions, or the flag gains a way to be cleared, is a
    /// design question for Milestone B — not something to settle in an accessor.
    #[test]
    fn a_filter_that_has_reached_the_end_stays_finished_after_a_reposition() {
        let mut filter = filter_over(
            FakeSource::new(vec![fake(FLAG_PAIRED, 60)]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(filter.next().is_some(), "the one record passes the filter");
        assert!(filter.next().is_none(), "and then the source is exhausted");

        filter.source_mut().rewind();
        assert_eq!(
            filter.source_mut().records_consumed(),
            0,
            "the source really did rewind"
        );
        assert!(
            filter.next().is_none(),
            "the filter is fused: repositioning its source does not restart it",
        );
    }

    /// The tally is **not** reset by a reposition, which is the behaviour a cursor wants:
    /// `counts()` then totals a whole chromosome rather than whichever region happened to be
    /// last. Stated as a test because the opposite is just as plausible a design.
    #[test]
    fn repositioning_the_source_does_not_reset_the_running_tally() {
        let mut filter = filter_over(
            FakeSource::new(vec![fake(FLAG_DUPLICATE, 60), fake(FLAG_PAIRED, 60)]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(
            filter.next().is_some(),
            "the second record passes; the first is a duplicate"
        );
        let after_one_pass = filter.counts();
        assert!(
            !after_one_pass.is_empty(),
            "the duplicate drop must have been tallied, or this test pins nothing",
        );

        filter.source_mut().rewind();
        assert_eq!(
            filter.counts(),
            after_one_pass,
            "repositioning the source discarded the tally",
        );

        // And it keeps **adding**, which is the half the doc comment's claim rests on:
        // reading the script a second time tallies it a second time, so `counts()` totals a
        // whole cursor rather than whichever region happened to be last. Non-erasure alone
        // would also be satisfied by a tally that stopped counting.
        assert!(filter.next().is_some(), "the rewound script replays");
        let after_two_passes = only_tally(&filter);
        assert_eq!(
            after_two_passes.duplicate, 2,
            "both passes' duplicate drops must be counted",
        );
        assert_eq!(
            after_two_passes.kept, 2,
            "both passes' keeps must be counted"
        );
    }

    /// **The guard that is the whole reason there are three states and not two.** A filter
    /// stopped by a fatal error must not be restarted: carrying on would read the next region
    /// out of a file already known to be corrupt or unreadable, and the error has been
    /// reported exactly once. A reviewer deleted this guard — making the restart
    /// unconditional — and the entire 1,503-test suite stayed green.
    #[test]
    fn a_failed_filter_is_not_restarted_and_says_so() {
        let mut undecodable = fake(FLAG_PAIRED, 60);
        undecodable.decode_fails = true;
        let mut filter = filter_over(
            FakeSource::new(vec![undecodable]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(matches!(
            filter.next(),
            Some(Err(ReadFilterError::Decode(_)))
        ));
        assert!(
            !filter.restart_after_end_of_input(),
            "a failed filter reports that there was nothing to restart",
        );

        filter.source_mut().rewind();
        assert!(
            filter.next().is_none(),
            "and it stays finished — the only way on is a new filter",
        );
    }

    /// The other side: a filter stopped by a clean end of input **is** restarted, and says so.
    /// Without this the two states are indistinguishable from the outside.
    #[test]
    fn a_filter_stopped_cleanly_is_restarted_and_says_so() {
        let mut filter = filter_over(
            FakeSource::new(vec![fake(FLAG_PAIRED, 60)]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(filter.next().is_some());
        assert!(filter.next().is_none(), "the source is exhausted");

        assert!(
            filter.restart_after_end_of_input(),
            "a cleanly finished filter reports that it restarted",
        );
        filter.source_mut().rewind();
        assert!(
            filter.next().is_some(),
            "and it yields the rewound script again",
        );
    }

    /// A running filter has nothing to restart, and must not claim it did.
    #[test]
    fn a_running_filter_reports_that_there_was_nothing_to_restart() {
        let mut filter = filter_over(
            FakeSource::new(vec![fake(FLAG_PAIRED, 60)]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(!filter.restart_after_end_of_input());
        assert!(filter.next().is_some(), "and it is untouched by the asking");
    }

    /// A filter finished by a **fatal error** stays finished too, and for the same reason —
    /// the flag is one flag. Documented as an invariant by `source_mut`; untested until now,
    /// which is how a "needs a fresh filter" instruction quietly becomes optional.
    #[test]
    fn a_filter_stopped_by_a_fatal_error_stays_finished_after_a_reposition() {
        let mut undecodable = fake(FLAG_PAIRED, 60);
        undecodable.decode_fails = true;
        let mut filter = filter_over(
            FakeSource::new(vec![undecodable]),
            poly_a_ref(30),
            ReadFilterConfig::default(),
        );

        assert!(
            matches!(filter.next(), Some(Err(ReadFilterError::Decode(_)))),
            "the record is built to fail its decode",
        );
        assert!(filter.next().is_none(), "and the filter fuses on it");

        filter.source_mut().rewind();
        assert_eq!(filter.source_mut().records_consumed(), 0);
        assert!(
            filter.next().is_none(),
            "a filter stopped by a fatal error is not restarted by repositioning its source",
        );
    }

    // **`contig_header` and `one_contig_header` went too (2026-08-03, B1), and their going is
    // worth a line.** Every test in this module built a SAM header, for one reason: the
    // `RecordSource` it constructed had to answer `header()` so the contig probe could count
    // `@SQ` lines. With the probe replaced by a contig-table comparison one layer up, no test
    // here needs a header at all — which is the same subtraction the module doc claims, seen
    // from the tests: **this module no longer knows what a header is.**

    // **The three `NoodlesRawAlignedRead` tests moved to `read/aligned_read.rs` with the type
    // (2026-08-03).** Their subject is the adapter's flag/MAPQ reads and its refusal to decode
    // an unstamped buffer, not any keep-or-drop rule.

    // **The four tests of `BamRecordSource`/`CramRecordSource` went with them (2026-08-03).**
    // Their subject was the sources: raw flag/MAPQ off a real file, the reused buffer not
    // leaking the previous record, and CRAM agreeing with noodles across container boundaries.
    // All three properties are asserted on the read path that survived —
    // `aligned_reads_reader::bam`'s `a_record_comes_out_of_the_bam_arm_raw` and its sibling walk
    // tests, `aligned_reads_reader::in_memory`'s buffer tests, and the run-of-regions oracles in
    // `open_bam` for the CRAM container walk.

    // **The drop-tally fixture moved to `read/input/cursor.rs` (2026-08-03).**
    //
    // Two tests here drove a fixture with one read per drop reason through a real BAM and a
    // real CRAM and asserted the tally by hand — including the one assertion that pins ng's
    // **#9-before-#8** order, since its last record fails both. They were the only tests worth
    // keeping when `BamRecordSource`/`CramRecordSource` were deleted, and they were never
    // about the source: what they pin is *this* filter's accounting.
    //
    // They live where the chain they need now composes — `AlignedReadsReader` →
    // `RegionRawAlignedReads` → `ReadFilter` → `AlignmentCursor` — as
    // `a_walk_charges_every_drop_reason_by_hand_count`.
    // A file is not needed to state the property and, since this module no longer knows what a
    // BAM is, could not be opened from here anyway.

    // **Three tests of the contig probe went with it (2026-08-03, B1), and one of them owed a
    // successor.**
    //
    // `read_filter_new_rejects_a_contig_missing_from_the_reference` pinned a real property —
    // a reference that does not match the file's contigs is refused **up front**, not
    // mid-stream. That property did not go away; it **moved up**, to
    // `AlignmentFile::cursor`, which now compares the two contig tables instead of fetching a
    // window per contig. Its successor is `a_cursor_refuses_an_accessor_over_a_different_
    // contig_table` in `read/input/open_bam.rs`, and it is stronger: the probe only proved
    // each contig *resolved*, so it accepted an accessor whose contigs had the right names
    // and the wrong lengths, which the successor rejects.
    //
    // The other two — `probe_free_constructor_filters_identically_to_new` and
    // `probe_free_constructor_skips_the_contig_probe_new_would_fail` — have **no successor and
    // should not have one**. Both compared the two constructors, and there is one constructor
    // now; "the probe is skipped" is not a property that can be stated about a probe that no
    // longer exists.

    /// Build a filter the way the cursor does — the stand-in for `ReadFilter::new`, which was
    /// deleted with the contig probe at B1.
    ///
    /// Infallible, which is the point: construction had exactly one failure mode and it was
    /// the probe. Every test below is about what the filter *does*, not about building one, so
    /// they take it through here rather than restating the buffers each time.
    fn filter_over<S: RecordSource, R: RawRefSeq>(
        source: S,
        reference: R,
        config: ReadFilterConfig,
    ) -> ReadFilter<S, R> {
        ReadFilter::with_validated_contigs(source, reference, config, ReadFilterBuffers::default())
    }

    /// A `RecordSource` whose `read_next` always fails — to drive the fatal path.
    struct ErroringSource;

    impl RecordSource for ErroringSource {
        type Record = FakeRecord;
        /// Nothing to skip: this fake carries no read-group resolution, so it
        /// can never meet another sample's record.
        fn other_sample_records(&self) -> u64 {
            0
        }

        fn read_next(&mut self, _buf: &mut FakeRecord) -> io::Result<bool> {
            Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated"))
        }
    }

    #[test]
    fn read_filter_source_read_error_is_fatal() {
        let source = ErroringSource;
        let mut filter = filter_over(source, poly_a_ref(100), ReadFilterConfig::default());
        // The read error surfaces as the yielded item …
        assert!(matches!(
            filter.next(),
            Some(Err(ReadFilterError::Source(_)))
        ));
        // … and the iterator is fused (stays stopped) afterward.
        assert!(filter.next().is_none());
    }

    #[test]
    fn read_filter_decode_error_is_fatal() {
        // A record that clears the pre-decode gate but fails to decode.
        let mut record = fake(FLAG_PAIRED, 60);
        record.decode_fails = true;
        let source = FakeSource::new(vec![record]);
        let mut filter = filter_over(source, poly_a_ref(100), ReadFilterConfig::default());
        assert!(matches!(
            filter.next(),
            Some(Err(ReadFilterError::Decode(_)))
        ));
        assert!(filter.next().is_none());
    }

    #[test]
    fn read_filter_reference_error_mid_stream_is_fatal() {
        // A read whose #8 reference window runs past the (short) contig end. The
        // contig resolves in `new` (a zero-length probe), so this fails only
        // in-loop → fatal Reference error.
        let read = mapped(&[b'A'; 50], &[30u8; 50], vec![CigarOp::Match(50)]);
        let source = FakeSource::new(vec![FakeRecord {
            read,
            decode_fails: false,
        }]);
        // contig length 10 < the read's 50-base reference span.
        let mut filter = filter_over(source, poly_a_ref(10), ReadFilterConfig::default());
        assert!(matches!(
            filter.next(),
            Some(Err(ReadFilterError::Reference(
                RefSeqError::OutOfBounds { .. }
            )))
        ));
        assert!(filter.next().is_none());
    }

    #[test]
    fn read_filter_charges_an_unmapped_read_end_to_end() {
        // The #5 counter the BAM/CRAM fixture omits: an unmapped read that clears
        // #2 (MAPQ 60) is charged to `unmapped` through the full iterator.
        let unmapped = fake(FLAG_PAIRED | FLAG_UNMAPPED, 60);
        let clean = fake(FLAG_PAIRED, 60);
        let source = FakeSource::new(vec![unmapped, clean]);
        let mut filter = filter_over(source, poly_a_ref(30), ReadFilterConfig::default());
        let kept: Vec<AlignedRead> = (&mut filter).collect::<Result<_, _>>().unwrap();
        assert_eq!(kept.len(), 1, "only the clean read survives");
        assert_eq!(only_tally(&filter).unmapped, 1);
        assert_eq!(only_tally(&filter).kept, 1);
    }

    #[test]
    fn read_filter_over_an_empty_source_yields_nothing_and_zero_counts() {
        let source = FakeSource::new(Vec::new());
        let mut filter = filter_over(source, poly_a_ref(30), ReadFilterConfig::default());
        assert!(filter.next().is_none());
        assert_eq!(only_tally(&filter), ReadFilterCounts::default());
    }

    #[test]
    fn read_filter_counts_is_a_running_tally_before_exhaustion() {
        // A duplicate (dropped) then a clean read (kept): the first `next()`
        // returns the kept read, and `counts()` already reflects both.
        let source = FakeSource::new(vec![
            fake(FLAG_PAIRED | FLAG_DUPLICATE, 60),
            fake(FLAG_PAIRED, 60),
        ]);
        let mut filter = filter_over(source, poly_a_ref(30), ReadFilterConfig::default());
        assert!(matches!(filter.next(), Some(Ok(_))));
        assert_eq!(only_tally(&filter).duplicate, 1);
        assert_eq!(only_tally(&filter).kept, 1);
        assert!(filter.next().is_none());
        assert_eq!(only_tally(&filter).kept, 1);
    }

    // -----------------------------------------------------------------
    // The probe-free constructor and the buffer hand-off (read_input A3)
    // -----------------------------------------------------------------

    // `mixed_batch` was here, and it went with `probe_free_constructor_filters_identically_to_
    // _new` (2026-08-03, B1) — its only caller. It existed to feed the *same* records through
    // both constructors and prove they agreed, and there is one constructor now.

    // Two tests lived here and both went with `into_parts` at Milestone F:
    // `into_parts_returns_the_buffers_with_their_allocations_and_the_tally` and
    // `with_validated_contigs_adopts_the_lent_buffers`. Each drove the lend-and-reclaim
    // protocol and nothing else, and nothing lends a filter its buffers any more — one filter
    // serves a whole chromosome, so it allocates its own once and reuses them for every
    // region.
}
