//! **Filtering on flag and mapping quality, with no reference anywhere in the scope** — spec §5's
//! capability, and the first of Milestone D's two tests output identity cannot see (plan D1).
//!
//! Nothing in ng needs it today. It is a capability the shape used to foreclose: filtering meant
//! building a `ReadFilter`, which until B1 cost a zero-length probe fetch for every contig in the
//! header, and until C2 still required a `RawRefSeq` to construct. Spec §5 names three callers
//! that would want it — a coverage histogram, an insert-size pass, a read-group pre-pass — and
//! what they would write is this file: a reader, the narrowing, the first filter, a tally.
//!
//! # What is actually pinned here, measured rather than argued
//!
//! **The property is over-determined, and D1's first draft added nothing.** That draft called
//! `verdict_on_raw_read` twice inside `filtering.rs` and asserted the verdicts. Review measured its
//! unique detection power at zero, and re-measuring confirms it: the checkpoint's mutation — a
//! reference parameter on the first filter — already fails at `cursor.rs:706`, which is *not a
//! test*, so `cargo build` breaks with none of this present. Narrowing the filter's visibility back
//! to its own module fails there too. Every layer this file drives already has a reference-free
//! test somewhere, written for other reasons: `region_raw_aligned_reads.rs`'s tests build the
//! reader and the narrowing without a reference, and `filtering.rs`'s own tests call the filter
//! without one.
//!
//! **So there is exactly one thing here that nothing else in the tree catches, and it is the
//! config.** Review broke spec §5's property outright with an earlier draft of this file green:
//! put a working reference behind `ReadFilterConfig`'s `Default`, have the first filter read it,
//! and every caller that says `ReadFilterConfig::default()` — or `..Default::default()` — notices
//! nothing. Measured here, with that field added, `Default` supplying it, and the tree's one other
//! exhaustive literal (`filtering.rs`'s `post_config`) repaired the way its author would repair it,
//! since a reference is already in scope in that module:
//!
//! ```text
//! error[E0063]: missing field `reference` in initializer of `filtering::ReadFilterConfig`
//!   --> src/ng/read/reference_free_first_filter.rs:100:5
//! error: could not compile `pop_var_caller` (lib test) due to 1 previous error
//! ```
//!
//! One error, in this file, and nowhere else. **That is why the config below is written out field
//! by field**, and why this module has to be somewhere a reference cannot be produced to repair it.
//!
//! **The pass is what makes that scope real.** A helper building a config proves nothing on its
//! own; it has to have a caller, and the honest caller is the whole pass spec §5 describes —
//! reader, narrowing, first filter, tally — with no `RawRefSeq` named, imported or constructible.
//! The pass also carries the ordinary behavioural weight: a flag routed to the wrong counter fails
//! it (verified by setting the secondary record's flag to `SUPPLEMENTARY`), though `filtering.rs`'s
//! own tests catch that too.
//!
//! **This module is a sibling of `filtering`, not a child.** `verdict_on_raw_read` and
//! `ReadFilterCounts::record_drop` are `pub(in crate::ng::read)`; a test inside `filtering.rs` sees
//! them however private they become. From here, narrowing either back to its module breaks this
//! file — co-detected by `cursor.rs`, so not a unique alarm, but it is the difference between
//! testing a filter and testing a *caller* of one.
//!
//! **What it still cannot do**, stated so nobody reads more into it: nothing stops a future author
//! reaching `crate::ng::ref_seq` by path from inside this file. The guard is that doing so is a
//! visible edit *to this file*, beside the paragraph saying why it must not happen — where the
//! config mechanism above catches a change made from anywhere else.
//!
//! # One thing this does not prove, and the design should say so
//!
//! Spec §5 says "the reference stops being a precondition for filtering at all". That is true of
//! the **filter**, and of the reader and the narrowing below it — which is what this file drives.
//! It is not true of `AlignmentCursor<R: RawRefSeq>` or `AlignmentFile::cursor`, whose bound is
//! unconditional, so §5's three callers still cannot reach a *file's* reads without producing a
//! reference. Raised at Checkpoint D; this file tests what is available and claims nothing more.

use crate::ng::read::aligned_read::{NoodlesRawAlignedRead, RawAlignedRead};
use crate::ng::read::filtering::{
    FilterVerdict, ReadFilterConfig, ReadFilterCounts, verdict_on_raw_read,
};
use crate::ng::read::input::aligned_reads_reader::{
    AlignedReadsReader, InMemoryAlignedReadsReader,
};
use crate::ng::read::input::region_raw_aligned_reads::RegionRawAlignedReads;
use crate::ng::read::input::test_fixtures::{
    bam_header, fixture_read_group, matching_contigs, read_named_with_length,
};
use crate::ng::types::{BaseQual, ContigId, GenomeRegion, MapQual, Position};
use noodles_sam::alignment::RecordBuf;
use noodles_sam::alignment::record::Flags;

/// A flag, a mapping quality and the thresholds — and nothing else.
///
/// **The signature is the property**, so a coercion is where it is written down: an added
/// parameter fails this line with `E0308`, and no import repairs it. The same mutation fails the
/// tree's ordinary call sites too, with `E0061` — this line is not the alarm, it is the statement,
/// and it is the only one that says the *whole* signature rather than its arity.
///
/// `#[cfg(test)]`-only, like this whole file, because Milestone D adds no shipping code. Beside
/// `verdict_on_raw_read` it would make `cargo build` carry the check as well; that is a change to a
/// non-test file and belongs to whoever wants it, not to a step scoped as tests.
const _FIRST_FILTER_TAKES_NO_REFERENCE: fn(u16, MapQual, &ReadFilterConfig) -> FilterVerdict =
    verdict_on_raw_read;

/// What a caller that never touches a reference would set: the flag and mapping-quality
/// thresholds, and **no** threshold for either filter that reads a decoded read.
///
/// **Written out field by field, and that is the one mechanism in this file with detection power
/// nothing else has** — the module doc has the measurement. `..Default::default()` here would
/// re-open the hole; `ReadFilterConfig::default()` would hide it completely.
fn config_of_a_caller_with_no_reference() -> ReadFilterConfig {
    ReadFilterConfig {
        min_mapq: Some(MapQual(20)),
        // `None`, not a value: this caller runs the first filter only, and #7 reads a length the
        // conversion produces.
        min_read_length: None,
        drop_qc_fail: true,
        drop_duplicate: true,
        // The threshold the one reference-touching filter uses. `None` is what makes this caller
        // reference-free in fact and not merely in scope — and the field still has to be named,
        // because the config is one struct for both filters.
        max_read_mismatch_fraction: None,
        mismatch_bq_floor: BaseQual(13),
    }
}

/// A 40-base read at `start` on contig 0, carrying `flags`.
///
/// 40 rather than the fixture's 10 so these are records a *whole* pass would also keep; nothing
/// here runs filter #7, which is the one that would notice.
fn read_flagged(qname: &str, start: usize, flags: Flags) -> RecordBuf {
    let mut record = read_named_with_length(qname, 0, start, 40);
    *record.flags_mut() = flags;
    record
}

/// **A whole pass over raw aligned reads with no reference in existence** — the shape spec §5's
/// three callers would write, and the caller that makes this file's scope a real one rather than a
/// claim about imports.
///
/// The six records cover one keep and five of the six drop reasons. `LowMapq` is left out because
/// the fixture builds every record at MAPQ 60; the mapping-quality threshold is exercised by
/// `filtering.rs`'s own boundary test, which is where a threshold belongs.
#[test]
fn a_pass_over_raw_reads_filters_and_charges_every_drop_with_no_reference() {
    let config = config_of_a_caller_with_no_reference();

    // Ascending positions: the narrowing stops at the first record beginning past the region, so
    // an out-of-order script would silently shorten the pass.
    let script = vec![
        read_flagged("kept", 10, Flags::empty()),
        read_flagged("duplicate", 11, Flags::DUPLICATE),
        read_flagged("supplementary", 12, Flags::SUPPLEMENTARY),
        read_flagged("secondary", 13, Flags::SECONDARY),
        // Unmapped **and placed**, which is the only way an unmapped record reaches a filter at
        // all: the narrowing compares a footprint against the region, so a record without an
        // alignment start is discarded a layer below and charged to nothing.
        read_flagged("unmapped", 14, Flags::UNMAPPED),
        read_flagged("qc-fail", 15, Flags::QC_FAIL),
    ];

    let mut reads = RegionRawAlignedReads::new(
        AlignedReadsReader::InMemory(InMemoryAlignedReadsReader::new(
            bam_header(&matching_contigs()),
            script,
        )),
        ContigId(0),
        fixture_read_group(),
    );
    reads
        .jump_to(GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(100),
        })
        .expect("this script has no scripted seek failure");

    let mut buffer = NoodlesRawAlignedRead::default();
    let mut counts = ReadFilterCounts::default();
    while reads
        .read_next(&mut buffer)
        .expect("this script has no scripted fault")
    {
        // The drop is charged to a read group the narrowing stamped, before anything decoded —
        // which is why the tally can be kept on this side of the conversion at all (spec §4).
        assert!(
            buffer.read_group().is_some(),
            "a raw aligned read reaches the first filter already attributable"
        );
        match verdict_on_raw_read(buffer.flag(), buffer.mapq(), &config) {
            FilterVerdict::Keep => counts.kept += 1,
            FilterVerdict::Drop(reason) => counts.record_drop(reason),
        }
    }

    assert_eq!(
        counts,
        ReadFilterCounts {
            kept: 1,
            duplicate: 1,
            supplementary: 1,
            secondary: 1,
            unmapped: 1,
            qc_fail: 1,
            ..ReadFilterCounts::default()
        },
        "a reference-free pass must keep one read and charge each of the five flagged drops"
    );
}
