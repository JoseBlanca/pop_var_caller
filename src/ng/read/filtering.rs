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
//! **What is here is the keep-or-drop rules and the thresholds they use.** The two verdicts —
//! [`verdict_on_raw_read`] on the flag and the mapping quality, [`verdict_on_aligned_read`] on
//! the length, the CIGAR and the mismatch fraction — plus the config they read, the drop
//! reasons, and the tally's counters.
//!
//! **The loop that calls them is not here, and its absence is the design** (spec §5). It moved
//! to [`AlignmentCursor`](crate::ng::read::input::cursor::AlignmentCursor) on 2026-08-03, with the
//! reused buffer, the reference and the running tally; before that, a `ReadFilter` iterator drove
//! the two verdicts and the conversion between them, and could not tell why its source had
//! stopped. The cursor causes region ends, so it never has to ask.
//!
//! The read these rules judge — [`RawAlignedRead`](crate::ng::read::aligned_read::RawAlignedRead)
//! undecoded, [`AlignedRead`] decoded, and the
//! conversion between them — lives in [`crate::ng::read::aligned_read`].

use crate::bam::alignment_input::{
    DEFAULT_MAX_READ_MISMATCH_FRACTION, DEFAULT_MIN_MAPQ, DEFAULT_MIN_READ_LENGTH,
    DEFAULT_MISMATCH_BQ_FLOOR, FLAG_DUPLICATE, FLAG_QC_FAIL, FLAG_SECONDARY, FLAG_SUPPLEMENTARY,
    FLAG_UNMAPPED, cigar_is_bad, cigar_ref_span, read_exceeds_mismatch_fraction,
};
use crate::ng::read::aligned_read::AlignedRead;
use crate::ng::ref_seq::{RawRefSeq, RefSeqError};
use crate::ng::types::{BaseQual, Bp, ContigId, MapQual, MismatchFraction};
// **No `noodles_bam`, `noodles_cram` or `noodles_fasta` here, and their absence is the
// point.** This module knows what the filtering policy is; it does not know what a BAM or a
// CRAM is. It did until 2026-08-03, when the two whole-file `RecordSource` implementations
// that lived below were deleted — see the note where they were. What a *record* is left the
// same day, to `aligned_read.rs`, with the two raw types.
// **No `noodles_sam as sam` here any more, and its absence is the point.** The last thing in
// this module that knew what a SAM *header* was, was the deleted source trait's `header` — the
// contig probe's input — and B1 took both. Since C3 the module names no noodles type at all: the
// two verdicts read a flag, a mapping quality and an `AlignedRead`, and nothing else.

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

/// Every provenance key [`ReadFilterConfig::provenance_parameters`] can write, so a
/// caller replacing an earlier recording can clear the whole family first — the key set
/// shrinks when the mismatch filter turns off, and a stale floor beside an `"off"` is
/// exactly the misreading the conditional key exists to prevent.
pub const READ_FILTER_PROVENANCE_KEYS: [&str; 6] = [
    "read-filter-min-mapq",
    "read-filter-min-read-length-bp",
    "read-filter-drop-qc-fail",
    "read-filter-drop-duplicates",
    "read-filter-max-read-mismatch-fraction",
    "read-filter-mismatch-bq-floor",
];

/// How an off filter is spelled in the header: a value, not an absent key — an absent
/// key would read as "unrecorded", which is a different fact.
const OFF: &str = "off";

impl ReadFilterConfig {
    /// This policy as psp-header provenance parameters — recorded, never compared
    /// (spec `run_streaming.md` §6.1): the census names the psp it was built from by
    /// digesting the header, filters included, so a header that omitted them could not
    /// tell a census built from filtered reads apart from one built from all of them.
    ///
    /// One key per **configurable** filter, each readable by eye; an off filter is the
    /// string `"off"`. The one absent-when-meaningless key is the mismatch base-quality
    /// floor, which exists only while the mismatch-fraction filter is on. The
    /// unconditional drops — secondary, supplementary and unmapped records, and a read
    /// whose CIGAR will not decode — are not settings, so no key records them.
    pub fn provenance_parameters(&self) -> Vec<(String, crate::ng::psp::ParameterValue)> {
        use crate::ng::psp::ParameterValue;

        // Exhaustive, so a filter added to this policy fails to compile here rather
        // than going unrecorded — which would make two differently-filtered walks
        // digest identically in the census, silently.
        let Self {
            min_mapq,
            min_read_length,
            drop_qc_fail,
            drop_duplicate,
            max_read_mismatch_fraction,
            mismatch_bq_floor,
        } = *self;

        let mut entries = vec![
            (
                "read-filter-min-mapq".to_string(),
                match min_mapq {
                    Some(quality) => ParameterValue::Integer(i64::from(quality.get())),
                    None => ParameterValue::String(OFF.to_string()),
                },
            ),
            (
                "read-filter-min-read-length-bp".to_string(),
                match min_read_length {
                    // A length is written as a number while one fits TOML's signed
                    // integer, and as its own digits in a string past that —
                    // unreachable from any real configuration, but a writer must never
                    // build a header its own reader refuses.
                    Some(length) => match i64::try_from(length.get()) {
                        Ok(length) => ParameterValue::Integer(length),
                        Err(_) => ParameterValue::String(length.get().to_string()),
                    },
                    None => ParameterValue::String(OFF.to_string()),
                },
            ),
            (
                "read-filter-drop-qc-fail".to_string(),
                ParameterValue::Boolean(drop_qc_fail),
            ),
            (
                "read-filter-drop-duplicates".to_string(),
                ParameterValue::Boolean(drop_duplicate),
            ),
        ];
        match max_read_mismatch_fraction {
            Some(fraction) => {
                entries.push((
                    "read-filter-max-read-mismatch-fraction".to_string(),
                    ParameterValue::Float(crate::ng::psp::header::wire_float_of(fraction.get())),
                ));
                entries.push((
                    "read-filter-mismatch-bq-floor".to_string(),
                    ParameterValue::Integer(i64::from(mismatch_bq_floor.get())),
                ));
            }
            None => entries.push((
                "read-filter-max-read-mismatch-fraction".to_string(),
                ParameterValue::String(OFF.to_string()),
            )),
        }
        entries
    }
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
    pub(in crate::ng::read) fn record_drop(&mut self, reason: DropReason) {
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
pub(in crate::ng::read) fn verdict_on_raw_read(
    flag: u16,
    mapq: MapQual,
    config: &ReadFilterConfig,
) -> FilterVerdict {
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
pub(in crate::ng::read) fn verdict_on_aligned_read(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::FLAG_PAIRED;
    use crate::ng::psp::ParameterValue;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::ReadGroupId;
    use crate::pileup::walker::CigarOp;

    /// Every configurable filter lands under its own key with its **exact** value, from
    /// a policy whose two booleans differ — so transposing the qc-fail and duplicate
    /// sources, or nudging any number, fails here rather than passing everywhere (a
    /// mutation pass found exactly that: four value mutations survived a test that
    /// pinned only two of the six).
    #[test]
    fn provenance_parameters_pin_every_filter_value() {
        let policy = ReadFilterConfig {
            min_mapq: Some(MapQual(17)),
            min_read_length: Some(Bp(31)),
            drop_qc_fail: true,
            drop_duplicate: false,
            max_read_mismatch_fraction: Some(
                MismatchFraction::try_new(0.25).expect("a quarter is in [0, 1]"),
            ),
            mismatch_bq_floor: BaseQual(13),
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            policy.provenance_parameters().into_iter().collect();

        assert_eq!(
            entries.get("read-filter-min-mapq"),
            Some(&ParameterValue::Integer(17))
        );
        assert_eq!(
            entries.get("read-filter-min-read-length-bp"),
            Some(&ParameterValue::Integer(31))
        );
        assert_eq!(
            entries.get("read-filter-drop-qc-fail"),
            Some(&ParameterValue::Boolean(true))
        );
        assert_eq!(
            entries.get("read-filter-drop-duplicates"),
            Some(&ParameterValue::Boolean(false))
        );
        assert_eq!(
            entries.get("read-filter-max-read-mismatch-fraction"),
            Some(&ParameterValue::Float(0.25))
        );
        assert_eq!(
            entries.get("read-filter-mismatch-bq-floor"),
            Some(&ParameterValue::Integer(13))
        );
        assert_eq!(entries.len(), 6, "six keys, nothing extra");
        for key in entries.keys() {
            assert!(
                READ_FILTER_PROVENANCE_KEYS.contains(&key.as_str()),
                "{key} is missing from the published key family"
            );
        }
    }

    /// An off filter is spelled `"off"`, never omitted — an absent key would read as
    /// "unrecorded", a different fact — and the one deliberate omission is the
    /// base-quality floor beside an off mismatch filter, which would invite misreading.
    #[test]
    fn provenance_parameters_spell_an_off_filter_rather_than_omitting_it() {
        let all_off = ReadFilterConfig {
            min_mapq: None,
            min_read_length: None,
            max_read_mismatch_fraction: None,
            ..ReadFilterConfig::default()
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            all_off.provenance_parameters().into_iter().collect();
        for key in [
            "read-filter-min-mapq",
            "read-filter-min-read-length-bp",
            "read-filter-max-read-mismatch-fraction",
        ] {
            assert_eq!(
                entries.get(key),
                Some(&ParameterValue::String("off".to_string())),
                "{key} must say off, not vanish"
            );
        }
        assert!(
            !entries.contains_key("read-filter-mismatch-bq-floor"),
            "a floor beside an off mismatch filter would invite misreading"
        );
        assert_eq!(entries.len(), 5, "five keys while the floor is meaningless");
    }

    /// A length past TOML's signed integer is recorded as its own digits in a string —
    /// unreachable from any real policy, but a wrapping-cast regression here would
    /// record a negative length silently.
    #[test]
    fn a_read_length_past_toml_integers_is_recorded_as_its_own_digits() {
        let policy = ReadFilterConfig {
            min_read_length: Some(Bp(u64::MAX)),
            ..ReadFilterConfig::default()
        };
        let entries: std::collections::BTreeMap<String, ParameterValue> =
            policy.provenance_parameters().into_iter().collect();
        assert_eq!(
            entries.get("read-filter-min-read-length-bp"),
            Some(&ParameterValue::String("18446744073709551615".to_string()))
        );
    }

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

    // ----- verdict_on_raw_read (#1–#6) -------------------------------------

    /// A record's flag/MAPQ pair; the pre-decode cascade needs nothing else.
    fn pre(flag: u16, mapq: u8, config: &ReadFilterConfig) -> FilterVerdict {
        verdict_on_raw_read(flag, MapQual(mapq), config)
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

    // ----- verdict_on_aligned_read (#7, #9, #8) -------------------------------

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
        verdict_on_aligned_read(read, reference, config, &mut buf).unwrap()
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
            verdict_on_aligned_read(&read, &empty, &enabled, &mut buf),
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
            verdict_on_aligned_read(&read, &reference, &post_config(Some(0.10)), &mut buf),
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
