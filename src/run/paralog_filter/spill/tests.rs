//! What the codec has to survive: an absent sample, a cohort of none, a cohort of a thousand,
//! a record that already carries a filter, a record with no annotations at all, and every way
//! a file can be cut or corrupted.
//!
//! **The comparison is by bit pattern throughout.** An absent sample is the `NaN` pair, and
//! `NaN == NaN` is false, so an equality-based round-trip test would fail on the one input the
//! file exists to carry — while a test that compared the two floats numerically at all would
//! pass on a codec that turned every absence into a zero (spec §6 trap 4).
//!
//! **[`holds_the_same_bits`] is what nine round trips assert, so it has its own negative
//! test.** A comparator that stopped discriminating would disarm all of them at once, and four
//! of the nine would then assert nothing but *the decode did not panic*.
//!
//! The tests that corrupt bytes by index all work on [`a_tiny_record`], whose body and fragments
//! are written out in [`the_encoding_matches_the_layout_the_spec_fixes`]; the offsets they use are
//! read off that list rather than counted from a fixture that could change under them. **They
//! corrupt the uncompressed fields and then frame them** ([`framed`]), because a byte changed
//! inside a compressed frame reaches no field — it is the frame's checksum that refuses it, which
//! [`a_damaged_frame_is_refused_and_not_read_as_a_record`] tests on its own.

use std::io::Cursor;

use proptest::prelude::*;

use super::{
    EncodedSpillEntry, MAX_BODY_BYTES, MAX_FRAGMENTS_BYTES, SpillEntry, SpillError, SpillReader,
    SpillReaderWithFragments, SpillWriter, SpilledSample, encode_body, encode_fragments,
    max_frame_bytes,
};
use crate::psp::varint::encode_u64_leb128;
use crate::read_phase::capture::LARGEST_ERROR_PHRED;
use crate::read_phase::fragments::{FragmentAllele, RecordFragments};
use crate::run::paralog_filter::WindowCoverage;
use crate::types::{ContigId, Position};

/// A window pair.
fn a_window(gc_fraction: f32, mean_depth: f32) -> WindowCoverage {
    WindowCoverage {
        gc_fraction,
        mean_depth,
    }
}

/// The `NaN` pair a sample with no usable window at this locus carries.
fn no_window() -> WindowCoverage {
    a_window(f32::NAN, f32::NAN)
}

/// A one-position sample with a window and both allele counts.
fn a_sample_with_a_window(
    gc_fraction: f32,
    mean_depth: f32,
    ref_reads: u32,
    alt_reads: u32,
) -> SpilledSample {
    SpilledSample {
        window: a_window(gc_fraction, mean_depth),
        ref_reads,
        alt_reads,
    }
}

/// A one-position sample with no usable window — the `NaN` pair — but reads all the same.
fn a_sample_without_a_window(ref_reads: u32, alt_reads: u32) -> SpilledSample {
    SpilledSample {
        window: no_window(),
        ref_reads,
        alt_reads,
    }
}

/// A biallelic SNP two samples covered, one of them with no window.
fn a_snp_two_samples_covered() -> SpillEntry {
    SpillEntry {
        contig: ContigId(3),
        position: Position(1_000_000),
        is_repeat_tract: false,
        line: b"SL4.0ch04\t1000000\t.\tA\tG\t42.5\tPASS\tAF=0.5;AC=1\tGT:GQ:AD\t0/1:30:5,5\t0/0:99:8,0"
            .to_vec(),
        samples: vec![
            a_sample_with_a_window(0.41, 6.25, 5, 5),
            a_sample_without_a_window(8, 0),
        ],
    }
}

/// A repeat tract two samples covered, one of them with no window. Its samples carry read counts
/// as every record's do (spec §3.2).
fn a_tract_two_samples_covered() -> SpillEntry {
    SpillEntry {
        contig: ContigId(3),
        position: Position(2_000_000),
        is_repeat_tract: true,
        line: b"SL4.0ch04\t2000000\t.\tATAT\tAT\t31.0\tPASS\tAN=4;DP=12;STR;RU=AT;PERIOD=2\tGT:GQ:DP:AD:REPCN\t1/1:30:6:0,6:2\t0/0:99:6:6,0:4"
            .to_vec(),
        samples: vec![
            a_sample_with_a_window(0.38, 5.75, 0, 6),
            a_sample_without_a_window(6, 0),
        ],
    }
}

/// The smallest entry with every field non-trivial, whose encoding is short enough to write
/// out byte by byte. Every test that corrupts a byte by index uses this one.
fn a_tiny_record() -> SpillEntry {
    SpillEntry {
        contig: ContigId(1),
        position: Position(300),
        is_repeat_tract: false,
        line: b"AB".to_vec(),
        samples: vec![a_sample_with_a_window(0.5, 2.0, 3, 130)],
    }
}

/// The fragments [`a_tiny_record`] is written with: one, of its one sample.
fn a_tiny_records_fragments() -> RecordFragments {
    RecordFragments::from_entries(vec![FragmentAllele {
        sample: 0,
        fragment: 5,
        allele: 1,
        error_phred: 30,
    }])
}

/// [`a_tiny_record`]'s entry as the writer writes it, with its fragments.
fn a_tiny_records_bytes() -> Vec<u8> {
    encoded_bytes_with_fragments(&[(a_tiny_record(), a_tiny_records_fragments())])
}

/// `entry`'s body, uncompressed.
fn raw_body(entry: &SpillEntry) -> Vec<u8> {
    let mut out = Vec::new();
    encode_body(entry, &mut out);
    out
}

/// `fragments`, uncompressed, as the fragments' frame holds them.
fn raw_fragments(fragments: &RecordFragments) -> Vec<u8> {
    let mut out = Vec::new();
    encode_fragments(fragments, &mut out);
    out
}

/// **An entry made of these uncompressed bytes**: `body` as one frame behind its length, then
/// `fragments` likewise, or a length of 0 where there are none. What the corruption tests read, so
/// that the bytes they change reach the field decoders.
fn framed(body: &[u8], fragments: Option<&[u8]>) -> Vec<u8> {
    let mut out = Vec::new();
    push_frame(body, &mut out);
    match fragments {
        Some(fragments) => push_frame(fragments, &mut out),
        None => encode_u64_leb128(0, &mut out),
    }
    out
}

/// Compress `raw` as one frame, as the writer does, and append it behind its length.
fn push_frame(raw: &[u8], out: &mut Vec<u8>) {
    let frame = zstd::bulk::compress(raw, 1).expect("zstd compresses any bytes");
    encode_u64_leb128(frame.len() as u64, out);
    out.extend_from_slice(&frame);
}

/// One entry's frames, inflated: its body, and its fragments or `None` for a length of 0.
fn unframed(entry: &[u8]) -> (Vec<u8>, Option<Vec<u8>>) {
    let mut at = entry;
    let body = next_frame(&mut at).expect("an entry starts with its body");
    let fragments = next_frame(&mut at);
    assert!(at.is_empty(), "one entry, and nothing after it");
    (body, fragments)
}

/// The frame at the front of `at`, inflated, or `None` for a length of 0; moves past it.
fn next_frame(at: &mut &[u8]) -> Option<Vec<u8>> {
    let (length, read) =
        crate::psp::varint::decode_u64_leb128(at).expect("a frame starts with its length");
    *at = &at[read..];
    if length == 0 {
        return None;
    }
    let (frame, rest) = at.split_at(length as usize);
    *at = rest;
    Some(zstd::decode_all(frame).expect("the writer's frames inflate"))
}

/// Offsets into [`a_tiny_record`]'s body, read off the byte list in
/// [`the_encoding_matches_the_layout_the_spec_fixes`].
const TRACT_FLAG_OFFSET_IN_TINY_RECORD: usize = 3;

/// Whether two entries hold the same bytes and the same bit patterns — the only comparison
/// that is meaningful over fields that can be `NaN`.
///
/// **Destructured, not field-accessed**: this is what "came back unchanged" means for every
/// round-trip test in this module, so a field either struct gains has to be answered for here
/// or it drops out of all of them at once. The precedents are the two this module's own header
/// cites — `var_calling::types`'s `LocusWindowCoverage` and `cohort_merge`'s `render`.
fn holds_the_same_bits(left: &SpillEntry, right: &SpillEntry) -> bool {
    let SpillEntry {
        contig,
        position,
        is_repeat_tract,
        line,
        samples,
    } = left;
    let SpillEntry {
        contig: other_contig,
        position: other_position,
        is_repeat_tract: other_is_repeat_tract,
        line: other_line,
        samples: other_samples,
    } = right;

    contig == other_contig
        && position == other_position
        && is_repeat_tract == other_is_repeat_tract
        && line == other_line
        && samples_hold_the_same_bits(samples, other_samples)
}

/// The samples, compared by bit pattern.
fn samples_hold_the_same_bits(left: &[SpilledSample], right: &[SpilledSample]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .zip(right)
            .all(|(one, other)| one_position_sample_holds_the_same_bits(one, other))
}

/// One one-position sample's four numbers, by bit pattern. Destructured for the same reason.
fn one_position_sample_holds_the_same_bits(one: &SpilledSample, other: &SpilledSample) -> bool {
    let SpilledSample {
        window,
        ref_reads,
        alt_reads,
    } = one;
    let SpilledSample {
        window: other_window,
        ref_reads: other_ref_reads,
        alt_reads: other_alt_reads,
    } = other;

    window_holds_the_same_bits(window, other_window)
        && ref_reads == other_ref_reads
        && alt_reads == other_alt_reads
}

/// The window pair, by bit pattern — the comparison the `NaN` absence needs.
fn window_holds_the_same_bits(window: &WindowCoverage, other: &WindowCoverage) -> bool {
    let WindowCoverage {
        gc_fraction,
        mean_depth,
    } = window;
    let WindowCoverage {
        gc_fraction: other_gc_fraction,
        mean_depth: other_mean_depth,
    } = other;

    gc_fraction.to_bits() == other_gc_fraction.to_bits()
        && mean_depth.to_bits() == other_mean_depth.to_bits()
}

/// The bytes the writer produces for these entries, each written without fragments.
fn encoded_bytes(entries: &[SpillEntry]) -> Vec<u8> {
    let entries: Vec<(SpillEntry, RecordFragments)> = entries
        .iter()
        .map(|entry| (entry.clone(), RecordFragments::default()))
        .collect();
    encoded_bytes_with_fragments(&entries)
}

/// The bytes the writer produces for these entries, each with its fragments.
fn encoded_bytes_with_fragments(entries: &[(SpillEntry, RecordFragments)]) -> Vec<u8> {
    let mut writer = SpillWriter::new(Vec::new());
    for (entry, fragments) in entries {
        writer
            .append(entry, fragments)
            .expect("a Vec sink never refuses bytes");
    }
    assert_eq!(writer.entries_written(), entries.len() as u64);
    writer.finish().expect("a Vec sink never refuses a flush")
}

/// Encode the entries and read them back.
fn round_trip(entries: &[SpillEntry]) -> Vec<SpillEntry> {
    SpillReader::new(Cursor::new(encoded_bytes(entries)), entries.len() as u64)
        .map(|entry| entry.expect("every entry the writer wrote reads back"))
        .collect()
}

/// Encode the entries with their fragments and read both back with pass three's reader.
fn round_trip_with_fragments(
    entries: &[(SpillEntry, RecordFragments)],
) -> Vec<(SpillEntry, RecordFragments)> {
    SpillReaderWithFragments::new(
        Cursor::new(encoded_bytes_with_fragments(entries)),
        entries.len() as u64,
    )
    .map(|entry| entry.expect("every entry the writer wrote reads back"))
    .collect()
}

/// Assert one entry survives the file unchanged, and hand back what came out of it.
fn assert_round_trips(entry: &SpillEntry) -> SpillEntry {
    let mut read = round_trip(std::slice::from_ref(entry));
    assert_eq!(read.len(), 1, "one entry in, one entry out");
    let came_back = read.remove(0);
    assert!(
        holds_the_same_bits(entry, &came_back),
        "the entry came back changed:\n  wrote {entry:?}\n  read  {came_back:?}"
    );
    came_back
}

/// The field a truncation blames, reached through the sample wrapper where there is one.
fn truncated_field(error: &SpillError) -> &'static str {
    match error {
        SpillError::Truncated { field } => field,
        SpillError::InSample { source, .. } => truncated_field(source),
        other => panic!("expected a truncation, got {other:?}"),
    }
}

/// A sink that refuses every byte and every flush, so the `Write` and `Flush` variants are
/// reachable from a test. Without it nothing in this module ever meets a stream that fails.
struct RefusesEverything;

impl std::io::Write for RefusesEverything {
    fn write(&mut self, _bytes: &[u8]) -> std::io::Result<usize> {
        Err(std::io::Error::other("the disk is full"))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Err(std::io::Error::other("the disk is full"))
    }
}

// ---------------------------------------------------------------------------
// The comparator every round trip rests on
// ---------------------------------------------------------------------------

#[test]
fn holds_the_same_bits_separates_entries_that_differ_in_any_one_field() {
    let base = a_snp_two_samples_covered();
    assert!(
        holds_the_same_bits(&base, &base.clone()),
        "an entry matches itself"
    );

    let differing: Vec<(&str, SpillEntry)> = vec![
        (
            "contig",
            SpillEntry {
                contig: ContigId(4),
                ..base.clone()
            },
        ),
        (
            "position",
            SpillEntry {
                position: Position(1_000_001),
                ..base.clone()
            },
        ),
        (
            "is_repeat_tract",
            SpillEntry {
                is_repeat_tract: true,
                ..base.clone()
            },
        ),
        (
            "line",
            SpillEntry {
                line: b"X".to_vec(),
                ..base.clone()
            },
        ),
        (
            "the sample count",
            SpillEntry {
                samples: vec![a_sample_with_a_window(0.41, 6.25, 5, 5)],
                ..base.clone()
            },
        ),
        (
            "gc_fraction",
            SpillEntry {
                samples: vec![
                    a_sample_with_a_window(0.42, 6.25, 5, 5),
                    a_sample_without_a_window(8, 0),
                ],
                ..base.clone()
            },
        ),
        (
            "mean_depth",
            SpillEntry {
                samples: vec![
                    a_sample_with_a_window(0.41, 6.26, 5, 5),
                    a_sample_without_a_window(8, 0),
                ],
                ..base.clone()
            },
        ),
        (
            "an absent window against a zeroed one",
            SpillEntry {
                samples: vec![
                    a_sample_with_a_window(0.41, 6.25, 5, 5),
                    a_sample_with_a_window(0.0, 0.0, 8, 0),
                ],
                ..base.clone()
            },
        ),
        (
            "ref_reads",
            SpillEntry {
                samples: vec![
                    a_sample_with_a_window(0.41, 6.25, 6, 5),
                    a_sample_without_a_window(8, 0),
                ],
                ..base.clone()
            },
        ),
        (
            "alt_reads",
            SpillEntry {
                samples: vec![
                    a_sample_with_a_window(0.41, 6.25, 5, 6),
                    a_sample_without_a_window(8, 0),
                ],
                ..base.clone()
            },
        ),
    ];

    for (what_differs, other) in differing {
        assert!(
            !holds_the_same_bits(&base, &other),
            "the comparator agreed on two entries differing in {what_differs}, so every round \
             trip that rests on it is asserting only that decoding did not panic"
        );
    }
}

// ---------------------------------------------------------------------------
// Round trips
// ---------------------------------------------------------------------------

#[test]
fn a_snp_with_a_sample_that_has_no_window_comes_back_with_the_absence_intact() {
    let entry = a_snp_two_samples_covered();
    let came_back = assert_round_trips(&entry);

    let absent = came_back
        .samples
        .get(1)
        .expect("two samples went in")
        .window;
    assert!(
        absent.gc_fraction.is_nan() && absent.mean_depth.is_nan(),
        "the absent sample came back as {absent:?}, and an absent sample that reads as a \
         number is a sample with no evidence looking like one with average coverage"
    );
}

#[test]
fn a_nan_comes_back_as_the_same_nan_and_not_merely_as_a_nan() {
    // A NaN with a payload of its own: a codec that reconstructed absence by writing
    // `f32::NAN` on the way out would pass an `is_nan()` test and fail this one.
    let odd_nan = f32::from_bits(0x7F80_0001);
    let entry = SpillEntry {
        samples: vec![SpilledSample {
            window: WindowCoverage {
                gc_fraction: odd_nan,
                mean_depth: f32::NAN,
            },
            ref_reads: 0,
            alt_reads: 0,
        }],
        ..a_snp_two_samples_covered()
    };

    let came_back = assert_round_trips(&entry);

    assert_eq!(
        came_back
            .samples
            .first()
            .expect("one sample went in")
            .window
            .gc_fraction
            .to_bits(),
        0x7F80_0001,
        "the float came back through arithmetic rather than as its bits"
    );
}

#[test]
fn a_record_with_no_samples_at_all_round_trips() {
    let entry = SpillEntry {
        samples: Vec::new(),
        ..a_snp_two_samples_covered()
    };
    assert!(assert_round_trips(&entry).samples.is_empty());
}

#[test]
fn the_columns_the_verdict_will_touch_are_carried_opaquely() {
    // The codec writes the line with `extend_from_slice` and reads it back by length, so
    // nothing branches on its contents: these two shapes cannot take a path the other round
    // trips do not. They are here because the plan names them, and what they prove is that a
    // record arriving with a `FILTER` already set, or with no annotations at all, is handed to
    // pass three exactly as it was written — the patch that rewrites those two columns is a
    // later step's.
    let lines: [&[u8]; 2] = [
        b"SL4.0ch04\t1000000\t.\tA\tG\t12.0\tEMNoConv\tAF=0.5;AC=1\tGT:GQ:AD\t0/1:12:3,3",
        b"SL4.0ch04\t1000000\t.\tA\t.\t0\tPASS\t.\tGT\t./.\t./.",
    ];

    for line in lines {
        let entry = SpillEntry {
            line: line.to_vec(),
            ..a_snp_two_samples_covered()
        };
        assert_eq!(assert_round_trips(&entry).line, line);
    }
}

#[test]
fn a_repeat_tract_round_trips_with_its_tract_flag_and_its_read_counts() {
    // Both facts must come back: the tract flag, which the writer's ordering rule reads, and the
    // counts, which the scorer reads (spec §3.2).
    let entry = SpillEntry {
        contig: ContigId(0),
        position: Position(1),
        is_repeat_tract: true,
        line: b"SL4.0ch00\t1\t.\tATATAT\tATATATAT\t61.0\tPASS\tAF=0.25\tGT:GQ:AD\t0/1:20:4,4"
            .to_vec(),
        samples: vec![
            a_sample_with_a_window(0.19, 3.5, 4, 4),
            a_sample_without_a_window(0, 3),
        ],
    };
    let came_back = assert_round_trips(&entry);
    assert!(came_back.is_repeat_tract);
    assert_eq!(
        (
            came_back.samples[0].ref_reads,
            came_back.samples[0].alt_reads
        ),
        (4, 4)
    );
}

#[test]
fn an_empty_line_round_trips() {
    // Not something the encoder produces, but the length prefix has to mean zero when it
    // says zero rather than run to the end of the entry.
    let entry = SpillEntry {
        line: Vec::new(),
        ..a_snp_two_samples_covered()
    };
    assert_round_trips(&entry);
}

#[test]
fn values_at_the_edges_of_their_fields_round_trip() {
    // `-0.0` is the one a bit-pattern codec can lose in a way `==` hides, which is the same
    // argument the module makes for `NaN`: `-0.0 == 0.0` is true.
    let entry = SpillEntry {
        contig: ContigId(u32::MAX),
        position: Position(u64::MAX),
        samples: vec![
            a_sample_with_a_window(f32::MAX, f32::MIN_POSITIVE, u32::MAX, u32::MAX),
            a_sample_with_a_window(-0.0, f32::INFINITY, 0, 0),
            a_sample_with_a_window(f32::NEG_INFINITY, 0.0, 1, 1),
        ],
        ..a_snp_two_samples_covered()
    };
    let came_back = assert_round_trips(&entry);
    assert_eq!(
        came_back
            .samples
            .get(1)
            .expect("three samples went in")
            .window
            .gc_fraction
            .to_bits(),
        (-0.0f32).to_bits(),
        "a negative zero came back as a positive one"
    );
}

#[test]
fn a_line_longer_than_one_varint_byte_round_trips() {
    // Every other fixture's line is under 128 bytes, so its length prefix is one byte. A real
    // record's is not: a 63-sample tomato line runs to several hundred bytes, and at a
    // thousand samples it is kilobytes.
    let long_line = vec![b'X'; 300];
    let entry = SpillEntry {
        line: long_line.clone(),
        ..a_tiny_record()
    };
    assert_eq!(assert_round_trips(&entry).line, long_line);
}

#[test]
fn a_cohort_of_a_thousand_samples_round_trips() {
    // The sample count's prefix is one byte below 128 and two above it, and every other
    // fixture has 0, 1 or 2 samples.
    let entry = SpillEntry {
        samples: (0..1000u32)
            .map(|index| {
                if index % 3 == 0 {
                    a_sample_without_a_window(index, 0)
                } else {
                    a_sample_with_a_window(0.4, index as f32, index, index + 1)
                }
            })
            .collect(),
        ..a_tiny_record()
    };
    let came_back = assert_round_trips(&entry);
    assert_eq!(came_back.samples.len(), 1000);
    let samples = &came_back.samples;
    assert_eq!(samples[999].ref_reads, 999);
}

#[test]
fn a_stream_of_records_comes_back_in_the_order_it_was_written() {
    let first = a_snp_two_samples_covered();
    let second = SpillEntry {
        position: Position(1_000_001),
        is_repeat_tract: true,
        line: b"SL4.0ch04\t1000001\t.\tAT\tATAT\t9.0\tPASS\t.\tGT\t0/1\t0/0".to_vec(),
        ..a_snp_two_samples_covered()
    };
    let third = SpillEntry {
        contig: ContigId(4),
        position: Position(7),
        samples: Vec::new(),
        ..a_snp_two_samples_covered()
    };

    let read = round_trip(&[first.clone(), second.clone(), third.clone()]);

    assert_eq!(read.len(), 3);
    assert!(holds_the_same_bits(&first, &read[0]));
    assert!(holds_the_same_bits(&second, &read[1]));
    assert!(holds_the_same_bits(&third, &read[2]));
}

/// **The fragments ride the spill, and only pass three's reader reads them** (read phase spec
/// §4.1): the reader pass three opens gets them back exactly, beside their entry; pass two's
/// steps over them by their length and still meets the next entry where it starts.
#[test]
fn fragments_are_read_back_by_pass_three_and_stepped_over_by_pass_two() {
    let deep = a_snp_two_samples_covered();
    let deep_fragments = RecordFragments::from_entries(
        (0..40)
            .map(|at| FragmentAllele {
                sample: at % 2,
                fragment: 1_000 + at * 7,
                allele: (at % 2) as u16,
                error_phred: 20 + (at % 30) as u8,
            })
            .collect(),
    );
    let after = a_tract_two_samples_covered();
    let written = [
        (deep.clone(), deep_fragments.clone()),
        (after.clone(), RecordFragments::default()),
    ];

    let read = round_trip_with_fragments(&written);
    assert!(holds_the_same_bits(&deep, &read[0].0));
    assert_eq!(read[0].1, deep_fragments);
    assert!(holds_the_same_bits(&after, &read[1].0));
    assert!(read[1].1.is_empty());

    let skipped: Vec<SpillEntry> =
        SpillReader::new(Cursor::new(encoded_bytes_with_fragments(&written)), 2)
            .collect::<Result<_, _>>()
            .expect("pass two reads both");
    assert!(holds_the_same_bits(&deep, &skipped[0]));
    assert!(holds_the_same_bits(&after, &skipped[1]));
}

#[test]
fn an_empty_file_yields_no_entries() {
    let mut reader = SpillReader::new(Cursor::new(Vec::new()), 0);
    assert!(reader.next_entry().is_none());
}

// ---------------------------------------------------------------------------
// The bytes themselves
// ---------------------------------------------------------------------------

#[test]
fn the_encoding_matches_the_layout_the_spec_fixes() {
    // Pins the field order and each field's width inside each frame, and the frames' order. A
    // codec whose encoder and decoder agree with each other but not with §3.4 passes every round
    // trip above and fails this.
    assert_eq!(
        raw_body(&a_tiny_record()),
        vec![
            0x01, // contig 1, one varint byte                          — offset 0
            0xAC, 0x02, // position 300, two varint bytes               — offsets 1, 2
            0x00, // is_repeat_tract = false                           — offset 3
            0x02, // the line is two bytes long                         — offset 4
            b'A', b'B', // the line                                     — offsets 5, 6
            0x01, // one sample                                         — offset 7
            0x00, 0x00, 0x00, 0x3F, // gc_fraction 0.5, an f32's bits little-endian
            0x00, 0x00, 0x00, 0x40, // mean_depth 2.0
            0x03, // ref_reads 3
            0x82, 0x01, // alt_reads 130, two varint bytes
        ]
    );
    assert_eq!(
        raw_fragments(&a_tiny_records_fragments()),
        vec![
            0x01, // one fragment allele                                — offset 0
            0x00, // sample 0, a step of 0 from the start               — offset 1
            0x05, // fragment 5, whole at its sample's first entry      — offset 2
            0x01, // allele 1                                           — offset 3
            0x1E, // error Phred 30                                     — offset 4
        ]
    );
    // And the entry the writer writes is the body's frame, then the fragments' frame.
    assert_eq!(
        unframed(&a_tiny_records_bytes()),
        (
            raw_body(&a_tiny_record()),
            Some(raw_fragments(&a_tiny_records_fragments()))
        )
    );
}

/// **The fragments are written as steps**: a sample as the step from the entry before, and a
/// fragment as the step from the one before within its sample — whole again at each sample's
/// first entry, which is what makes the steps short.
#[test]
fn fragments_are_written_as_steps_that_restart_at_each_sample() {
    let fragments = RecordFragments::from_entries(vec![
        FragmentAllele {
            sample: 0,
            fragment: 1_000,
            allele: 0,
            error_phred: 30,
        },
        FragmentAllele {
            sample: 0,
            fragment: 1_003,
            allele: 1,
            error_phred: 30,
        },
        FragmentAllele {
            sample: 2,
            fragment: 900,
            allele: 0,
            error_phred: 20,
        },
    ]);
    assert_eq!(
        raw_fragments(&fragments),
        vec![
            0x03, // three fragment alleles
            // sample 0; fragment 1,000 whole; allele 0; Phred 30
            0x00, 0xE8, 0x07, 0x00, 0x1E,
            // sample 0 again; fragment 1,003 as a step of 3
            0x00, 0x03, 0x01, 0x1E,
            // sample 2, a step of 2; fragment 900 whole again; Phred 20
            0x02, 0x84, 0x07, 0x00, 0x14,
        ]
    );
}

/// **A record with no fragments writes no frame for them**, only a length of 0 — so a run whose
/// records carry none pays one byte an entry for them, not an empty frame's dozen.
#[test]
fn a_record_with_no_fragments_writes_a_length_of_zero_and_no_frame() {
    let written = encoded_bytes(&[a_tiny_record()]);
    assert_eq!(written.last(), Some(&0x00));
    assert_eq!(unframed(&written), (raw_body(&a_tiny_record()), None));
    let (entry, fragments) =
        read_with_fragments(written).expect("pass three reads an entry without fragments");
    assert!(holds_the_same_bits(&entry, &a_tiny_record()));
    assert!(fragments.is_empty());
}

/// **Real entries shrink**: a cohort's line repeats the same few genotype strings, and the frame
/// is where that is taken out. Pinned loosely — a third — so that a frame written uncompressed by
/// mistake fails, and a zstd version that compresses a little differently does not.
#[test]
fn a_cohort_whose_genotypes_repeat_is_written_in_a_fraction_of_its_bytes() {
    let samples = 1_000;
    let mut line = b"SL4.0ch01\t1000\t.\tA\tG\t50\tPASS\t.\tGT:AD:DP".to_vec();
    for at in 0..samples {
        line.extend_from_slice(if at % 3 == 0 {
            b"\t0/0:5,0:5"
        } else {
            b"\t1/1:0,4:4"
        });
    }
    let entry = SpillEntry {
        line,
        samples: vec![a_sample_with_a_window(0.4, 5.0, 5, 0); samples],
        ..a_tiny_record()
    };
    let written = encoded_bytes(std::slice::from_ref(&entry)).len();
    let raw = raw_body(&entry).len();
    assert!(
        written * 3 < raw,
        "{raw} bytes of body were written as {written}"
    );
}

// ---------------------------------------------------------------------------
// Files that are cut, and files that are corrupt
// ---------------------------------------------------------------------------

#[test]
fn a_body_cut_anywhere_inside_names_the_field_the_bytes_ran_out_in() {
    // A body that fails is not enough: the field name is the whole diagnostic a spill failure
    // carries, so every cut has to blame the right one. The offsets are a_tiny_record's body,
    // from the byte list above; each cut body is framed whole, so it is the field decoders that
    // meet the end.
    let expected: &[(usize, &str)] = &[
        (0, "contig"),
        (1, "position"),
        (2, "position"),
        (3, "is_repeat_tract"),
        (4, "line_length"),
        (5, "line"),
        (6, "line"),
        (7, "sample_count"),
        (8, "gc_fraction"),
        (9, "gc_fraction"),
        (10, "gc_fraction"),
        (11, "gc_fraction"),
        (12, "mean_depth"),
        (13, "mean_depth"),
        (14, "mean_depth"),
        (15, "mean_depth"),
        (16, "ref_reads"),
        (17, "alt_reads"),
        (18, "alt_reads"),
    ];

    let whole = raw_body(&a_tiny_record());
    assert_eq!(
        whole.len(),
        19,
        "the byte list moved and the offsets below move with it"
    );

    for &(cut, field) in expected {
        let mut reader = SpillReader::new(Cursor::new(framed(&whole[..cut], None)), 1);
        match reader.next_entry() {
            Some(Err(error)) => assert_eq!(
                truncated_field(&error),
                field,
                "a body cut after {cut} bytes blamed {} rather than {field}",
                truncated_field(&error)
            ),
            other => panic!("a body cut after {cut} bytes gave {other:?} rather than failing"),
        }
    }
}

/// **A file cut anywhere inside an entry is refused, naming the part it ran out in**: a frame's
/// length, or the frame itself. A spill that lost its tail mid-entry is a partial final write.
#[test]
fn a_file_cut_anywhere_inside_an_entry_names_the_frame_the_bytes_ran_out_in() {
    let whole = a_tiny_records_bytes();
    // The body's frame is under 128 bytes, so its length is one byte.
    let body_frame = usize::from(whole[0]);
    let fragments_length_at = 1 + body_frame;
    for cut in 1..whole.len() {
        let field = if cut < fragments_length_at {
            "body"
        } else if cut == fragments_length_at {
            "fragments_bytes"
        } else {
            "fragments"
        };
        match read_with_fragments(whole[..cut].to_vec()) {
            Err(error) => assert_eq!(truncated_field(&error), field, "a cut after {cut} bytes"),
            other => panic!("a file cut after {cut} bytes gave {other:?} rather than failing"),
        }
    }
}

#[test]
fn a_failure_inside_the_samples_names_which_sample() {
    // At a cohort of three thousand a message that names only the field points at three
    // thousand places at once.
    let entry = SpillEntry {
        samples: vec![a_sample_with_a_window(0.5, 2.0, 1, 1); 3],
        ..a_tiny_record()
    };
    let whole = raw_body(&entry);
    let head = raw_body(&SpillEntry {
        samples: Vec::new(),
        ..entry
    })
    .len();
    let bytes_per_sample = (whole.len() - head) / 3;
    let two_bytes_into_the_second_sample = head + bytes_per_sample + 2;

    let mut reader = SpillReader::new(
        Cursor::new(framed(&whole[..two_bytes_into_the_second_sample], None)),
        1,
    );

    match reader.next_entry() {
        Some(Err(SpillError::InSample { index, source })) => {
            assert_eq!(index, 1);
            assert_eq!(truncated_field(&source), "gc_fraction");
        }
        other => panic!("expected a failure naming its sample, got {other:?}"),
    }
}

#[test]
fn a_tract_flag_byte_that_is_neither_zero_nor_one_is_refused() {
    let mut body = raw_body(&a_tiny_record());
    body[TRACT_FLAG_OFFSET_IN_TINY_RECORD] = 2;
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::NotABoolean { field, byte })) => {
            assert_eq!(field, "is_repeat_tract");
            assert_eq!(byte, 2);
        }
        other => panic!("expected a refused flag byte, got {other:?}"),
    }
}

/// **A body holding bytes past its last sample is refused**: the frame's size and the body's own
/// counts disagree, so one of them is damaged.
#[test]
fn a_body_longer_than_its_fields_is_refused() {
    let mut body = raw_body(&a_tiny_record());
    body.extend_from_slice(&[0x00, 0x00, 0x00]);
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::SectionLongerThanItsEntries { field, unread })) => {
            assert_eq!((field, unread), ("body", 3));
        }
        other => panic!("a body with bytes past its samples gave {other:?}"),
    }
}

#[test]
fn a_contig_too_large_for_the_field_is_refused_rather_than_wrapped() {
    // A varint holding 2^32 — one past the largest contig index there can be. Written by
    // hand because the encoder cannot produce it.
    let mut body = Vec::new();
    encode_u64_leb128(1 << 32, &mut body);
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!(field, "contig");
            assert_eq!(value, 1 << 32);
        }
        other => panic!("expected a refused contig, got {other:?}"),
    }
}

#[test]
fn a_read_count_too_large_for_its_field_is_refused() {
    // `ref_reads` and `alt_reads` are `u32` behind a `u64` varint; a decoder that read them
    // with `read_varint` instead of `read_u32` would wrap silently.
    let mut body = Vec::new();
    encode_u64_leb128(0, &mut body); // contig
    encode_u64_leb128(1, &mut body); // position
    body.push(0); // not a tract, so the samples carry read counts
    encode_u64_leb128(0, &mut body); // an empty line
    encode_u64_leb128(1, &mut body); // one sample
    body.extend_from_slice(&0.5f32.to_bits().to_le_bytes());
    body.extend_from_slice(&1.0f32.to_bits().to_le_bytes());
    encode_u64_leb128(1 << 32, &mut body); // ref_reads, one past the field
    encode_u64_leb128(0, &mut body);

    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::InSample { index, source })) => {
            assert_eq!(index, 0);
            match *source {
                SpillError::OutOfRange { field, value } => {
                    assert_eq!(field, "ref_reads");
                    assert_eq!(value, 1 << 32);
                }
                other => panic!("expected a refused read count, got {other:?}"),
            }
        }
        other => panic!("expected a refused read count, got {other:?}"),
    }
}

#[test]
fn an_over_long_varint_is_refused() {
    // Eleven continuation bytes: longer than any `u64` needs, so corruption rather than a
    // value — inside a frame, and as a frame's own length.
    let mut reader = SpillReader::new(Cursor::new(framed(&[0xFF; 11], None)), 1);
    match reader.next_entry() {
        Some(Err(SpillError::OverlongVarint { field })) => assert_eq!(field, "contig"),
        other => panic!("expected an over-long varint, got {other:?}"),
    }

    let mut reader = SpillReader::new(Cursor::new(vec![0xFF; 11]), 1);
    match reader.next_entry() {
        Some(Err(SpillError::OverlongVarint { field })) => assert_eq!(field, "body_bytes"),
        other => panic!("expected an over-long frame length, got {other:?}"),
    }
}

#[test]
fn a_ten_byte_varint_holding_more_than_a_u64_is_refused_rather_than_truncated() {
    // Ten bytes is a legal LEB128 width, but the tenth byte contributes `data << 63`, so all
    // but its lowest bit falls off the top of the `u64`. The psp primitive returns `Ok` with
    // the high bits dropped; this codec refuses corruption rather than absorbing it.
    let mut body = vec![0x80; 9];
    body.push(0x7F);
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OverlongVarint { field })) => assert_eq!(field, "contig"),
        other => panic!("expected a refused varint, got {other:?}"),
    }
}

#[test]
fn a_line_longer_than_the_ceiling_is_refused_before_it_is_read() {
    // Without the ceiling the length is believed and every byte the body still holds is
    // appended before the short-read check fires.
    let mut body = Vec::new();
    encode_u64_leb128(0, &mut body); // contig
    encode_u64_leb128(1, &mut body); // position
    body.push(0); // not a tract, so generic-locus rows
    encode_u64_leb128(u64::from(u32::MAX), &mut body); // a line of four billion bytes
    body.extend_from_slice(&vec![b'X'; 4096]);

    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!(field, "line_length");
            assert_eq!(value, u64::from(u32::MAX));
        }
        other => panic!("expected a refused line length, got {other:?}"),
    }
}

/// **A frame's length past what any frame can be is refused before a byte of the frame is
/// read**, so one damaged length cannot make the reader hold the rest of a multi-terabyte file.
#[test]
fn a_frame_length_past_the_ceiling_is_refused_before_the_frame_is_read() {
    let mut bytes = Vec::new();
    let too_long = max_frame_bytes(MAX_BODY_BYTES) as u64 + 1;
    encode_u64_leb128(too_long, &mut bytes);
    bytes.extend_from_slice(&[0u8; 4096]);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!((field, value), ("body_bytes", too_long));
        }
        other => panic!("expected a refused frame length, got {other:?}"),
    }
}

/// **A frame that claims to inflate past the ceiling is refused before anything is allocated for
/// it**: the frame records its size, and that is checked first. Zeros compress to almost
/// nothing, so a small frame can claim a large size.
#[test]
fn a_frame_claiming_more_than_the_ceiling_is_refused_before_it_is_inflated() {
    let too_large = vec![0u8; MAX_BODY_BYTES + 1];
    let mut reader = SpillReader::new(Cursor::new(framed(&too_large, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!((field, value), ("body", MAX_BODY_BYTES as u64 + 1));
        }
        other => panic!("expected a refused frame size, got {other:?}"),
    }
}

/// **A frame damaged inside is refused, not read as a record**: the byte changed is in the
/// compressed data, where no field decoder sees it, and it is the frame's checksum or zstd's own
/// checks that catch it.
#[test]
fn a_damaged_frame_is_refused_and_not_read_as_a_record() {
    let whole = a_tiny_records_bytes();
    let body_frame = usize::from(whole[0]);
    // Every byte of the body's frame after its four-byte magic number, flipped one at a time.
    for at in 1 + 4..1 + body_frame {
        let mut bytes = whole.clone();
        bytes[at] ^= 0x55;
        match read_with_fragments(bytes) {
            Err(SpillError::NotAFrame { field, .. }) => assert_eq!(field, "body", "byte {at}"),
            // A damaged frame header can also make the frame claim a size past the ceiling.
            Err(SpillError::OutOfRange { field, .. }) => assert_eq!(field, "body", "byte {at}"),
            other => panic!("a frame damaged at byte {at} gave {other:?}"),
        }
    }
}

/// **A frame that does not record its size is refused**: the writer always records it, and the
/// reader needs it to bound what it allocates.
#[test]
fn a_frame_that_does_not_record_its_size_is_refused() {
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), 1).expect("a zstd encoder over a Vec");
    encoder
        .include_contentsize(false)
        .expect("the flag is accepted");
    std::io::Write::write_all(&mut encoder, &raw_body(&a_tiny_record())).expect("a Vec");
    let frame = encoder.finish().expect("the frame ends");
    let mut bytes = Vec::new();
    encode_u64_leb128(frame.len() as u64, &mut bytes);
    bytes.extend_from_slice(&frame);
    encode_u64_leb128(0, &mut bytes);

    match SpillReader::new(Cursor::new(bytes), 1).next_entry() {
        Some(Err(SpillError::NotAFrame { field, reason })) => {
            assert_eq!(field, "body");
            assert_eq!(reason, "the frame does not record its size");
        }
        other => panic!("a frame without its size gave {other:?}"),
    }
}

#[test]
fn a_sample_count_larger_than_the_body_is_an_error_and_not_an_allocation() {
    // Nine hundred thousand samples — under the ceiling, so the decoder believes the count —
    // and a body that ends immediately. It must fail on the bytes it does not have rather
    // than reserve for the count it was told: without the reservation cap this allocates
    // 900,000 x 16 bytes before reading one.
    let mut body = vec![
        0x00, // contig 0
        0x01, // position 1
        0x00, // not a tract, so generic-locus rows
        0x00, // an empty line
    ];
    encode_u64_leb128(900_000, &mut body);
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::InSample { index, source })) => {
            assert_eq!(index, 0);
            assert_eq!(truncated_field(&source), "gc_fraction");
        }
        other => panic!("expected a truncated sample, got {other:?}"),
    }
}

#[test]
fn a_sample_count_past_the_ceiling_is_refused_before_the_samples_are_built() {
    // The reservation cap bounds what is *reserved*; without a cap on the count itself the
    // loop keeps decoding until the body runs out. The bytes below are a head followed by a
    // count of ten million and enough sample bytes that a decoder without the ceiling would
    // build them all.
    let mut body = vec![
        0x00, // contig 0
        0x01, // position 1
        0x00, // not a tract, so generic-locus rows
        0x00, // an empty line
    ];
    encode_u64_leb128(10_000_000, &mut body);
    body.extend_from_slice(&vec![0u8; 4096]);
    let mut reader = SpillReader::new(Cursor::new(framed(&body, None)), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!(field, "sample_count");
            assert_eq!(value, 10_000_000);
        }
        other => panic!("expected a refused sample count, got {other:?}"),
    }
}

#[test]
fn a_file_holding_more_records_than_were_written_is_refused() {
    // The completeness check runs in both directions: a file that grew is as wrong as one that
    // lost its tail, and a comparison written as "fewer than expected" would miss this.
    let bytes = encoded_bytes(&[a_tiny_record(), a_tiny_record()]);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    assert!(reader.next_entry().expect("the first entry").is_ok());
    match reader.next_entry() {
        Some(Ok(_)) => match reader.next_entry() {
            Some(Err(SpillError::TheWrongNumberOfRecords { expected, read })) => {
                assert_eq!(expected, 1);
                assert_eq!(read, 2);
            }
            other => panic!("expected the extra record to be refused, got {other:?}"),
        },
        other => panic!("expected a second entry, got {other:?}"),
    }
}

#[test]
fn a_source_that_fails_mid_frame_says_the_read_failed_and_not_that_the_file_ended() {
    // A device failure and a file that ends are different things.
    struct FailsAfterTheStart {
        remaining: Vec<u8>,
    }

    impl std::io::Read for FailsAfterTheStart {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            if self.remaining.is_empty() {
                return Err(std::io::Error::other("the device went away"));
            }
            let taken = out.len().min(self.remaining.len());
            out[..taken].copy_from_slice(&self.remaining[..taken]);
            self.remaining.drain(..taken);
            Ok(taken)
        }
    }

    // The body's length and the first half of its frame, then the device fails.
    let whole = a_tiny_records_bytes();
    let head = whole[..1 + usize::from(whole[0]) / 2].to_vec();
    let source = std::io::BufReader::with_capacity(1, FailsAfterTheStart { remaining: head });
    let mut reader = SpillReader::new(source, 1);

    match reader.next_entry() {
        Some(Err(SpillError::Read { field, .. })) => assert_eq!(field, "body"),
        other => panic!("expected a read failure inside the frame, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// The reader stops, and the writer reports its sink
// ---------------------------------------------------------------------------

#[test]
fn the_reader_stops_after_a_decode_error() {
    // A stream that failed once is not at an entry boundary any more, so bytes from the middle
    // of an entry would decode — plausibly — into a record that was never written.
    let mut bytes = framed(&[0x01, 0x01, 0x02], None); // contig, position, then a flag byte of 2
    bytes.extend_from_slice(&a_tiny_records_bytes());
    let mut reader = SpillReader::new(Cursor::new(bytes), 2);

    assert!(matches!(
        reader.next_entry(),
        Some(Err(SpillError::NotABoolean { .. }))
    ));
    let after = reader.next_entry();
    assert!(
        after.is_none(),
        "the reader kept going after an error and produced {after:?}"
    );
}

#[test]
fn append_reports_a_sink_that_refused_the_bytes_and_does_not_count_the_entry() {
    let mut writer = SpillWriter::new(RefusesEverything);
    assert!(matches!(
        writer.append(&a_tiny_record(), &a_tiny_records_fragments()),
        Err(SpillError::Write { .. })
    ));
    assert_eq!(writer.entries_written(), 0);
}

#[test]
fn finish_reports_a_flush_that_failed() {
    // A `finish` that swallowed this would report a complete spill for a truncated one, and
    // pass two would score a cohort missing its last records.
    assert!(matches!(
        SpillWriter::new(RefusesEverything).finish(),
        Err(SpillError::Flush { .. })
    ));
}

// ---------------------------------------------------------------------------
// The round-trip law, over generated entries
// ---------------------------------------------------------------------------

/// Floats the fixtures reach for by hand, plus any bit pattern at all.
fn any_float() -> impl Strategy<Value = f32> {
    prop_oneof![
        Just(f32::NAN),
        Just(f32::from_bits(0x7F80_0001)),
        Just(0.0f32),
        Just(-0.0f32),
        Just(f32::INFINITY),
        Just(f32::NEG_INFINITY),
        any::<u32>().prop_map(f32::from_bits),
    ]
}

/// An entry of any shape the writer will accept — lines and cohorts on both sides of the
/// 128-byte varint boundary, and both row shapes.
///
/// **The shape follows the tract flag** rather than being drawn beside it, which is the writer's
/// rule and now also the file's: there is one flag, and it decides.
fn any_entry() -> impl Strategy<Value = (SpillEntry, RecordFragments)> {
    (
        any::<u32>(),
        any::<u64>(),
        any::<bool>(),
        prop::collection::vec(any::<u8>(), 0..400),
        prop::collection::vec(
            (any_float(), any_float(), any::<u32>(), any::<u32>()),
            0..300,
        ),
        // Fragments: sample indices are reduced into the entry's sample count below; fragment
        // numbers and alleles span every varint length their types allow.
        prop::collection::vec(
            (
                any::<u32>(),
                any::<u32>(),
                any::<u16>(),
                0..=LARGEST_ERROR_PHRED,
            ),
            0..40,
        ),
    )
        .prop_map(
            |(contig, position, is_repeat_tract, line, samples, drawn)| {
                let samples: Vec<SpilledSample> = samples
                    .into_iter()
                    .map(|(gc, depth, ref_reads, alt_reads)| {
                        a_sample_with_a_window(gc, depth, ref_reads, alt_reads)
                    })
                    .collect();
                let sample_count = samples.len() as u32;
                let fragments = RecordFragments::from_entries(if sample_count == 0 {
                    Vec::new()
                } else {
                    drawn
                        .into_iter()
                        .map(|(sample, fragment, allele, error_phred)| FragmentAllele {
                            sample: sample % sample_count,
                            fragment,
                            allele,
                            error_phred,
                        })
                        .collect()
                });
                (
                    SpillEntry {
                        contig: ContigId(contig),
                        position: Position(position),
                        is_repeat_tract,
                        line,
                        samples,
                    },
                    fragments,
                )
            },
        )
}

proptest! {
    #[test]
    fn any_stream_of_entries_comes_back_bit_for_bit(
        entries in prop::collection::vec(any_entry(), 1..8)
    ) {
        let read = round_trip_with_fragments(&entries);
        prop_assert_eq!(read.len(), entries.len());
        for ((wrote, fragments), (came_back, fragments_back)) in entries.iter().zip(&read) {
            prop_assert!(holds_the_same_bits(wrote, came_back));
            prop_assert_eq!(fragments, fragments_back);
        }
        // Pass two's reader steps over every section and meets every entry where it starts.
        let bare: Vec<SpillEntry> = entries.iter().map(|(entry, _)| entry.clone()).collect();
        let skipped: Vec<SpillEntry> = SpillReader::new(
            Cursor::new(encoded_bytes_with_fragments(&entries)),
            entries.len() as u64,
        )
        .collect::<Result<_, _>>()
        .expect("pass two reads every entry");
        prop_assert_eq!(skipped.len(), bare.len());
        for (wrote, came_back) in bare.iter().zip(&skipped) {
            prop_assert!(holds_the_same_bits(wrote, came_back));
        }
    }
}

/// **An entry encoded ahead of time writes the bytes `append` writes**, for both row shapes and a
/// record with no samples — so moving the encoding to the thread that called the locus changes
/// nothing in the file.
#[test]
fn an_entry_encoded_ahead_of_time_writes_what_append_writes() {
    let entries = [
        (a_snp_two_samples_covered(), RecordFragments::default()),
        (a_tract_two_samples_covered(), a_tiny_records_fragments()),
        (a_tiny_record(), a_tiny_records_fragments()),
    ];
    let mut writer = SpillWriter::new(Vec::new());
    for (entry, fragments) in &entries {
        let encoded = EncodedSpillEntry::of(entry, fragments).expect("a well-formed entry encodes");
        writer
            .append_encoded(&encoded)
            .expect("a Vec takes every byte");
    }
    assert_eq!(writer.entries_written(), entries.len() as u64);
    let written = writer.finish().expect("a Vec flushes");
    assert_eq!(written, encoded_bytes_with_fragments(&entries));
}

// ---------------------------------------------------------------------------
// A damaged fragments frame, read by pass three
// ---------------------------------------------------------------------------

/// What pass three's reader makes of these bytes, holding one entry.
fn read_with_fragments(bytes: Vec<u8>) -> Result<(SpillEntry, RecordFragments), SpillError> {
    SpillReaderWithFragments::new(Cursor::new(bytes), 1)
        .next()
        .expect("the reader yields something for one entry's bytes")
}

/// **Each kind of damage to the fragments is refused, naming its field** (review Mi2): a count
/// the frame cannot hold, bytes past the last entry, an allele past `u16`, a frame larger than any
/// real one, a Phred above the largest stored, a sample the entry does not have, and steps that
/// run past the end of a `u32`. Read as data, the sample and the Phred would change or drop phase
/// silently. The offsets are [`a_tiny_records_fragments`]'s, from the byte list in
/// [`the_encoding_matches_the_layout_the_spec_fixes`].
#[test]
fn a_damaged_fragments_frame_is_refused() {
    let body = raw_body(&a_tiny_record());
    let section = raw_fragments(&a_tiny_records_fragments());
    let with = |edit: &dyn Fn(&mut Vec<u8>)| {
        let mut bytes = section.clone();
        edit(&mut bytes);
        read_with_fragments(framed(&body, Some(&bytes)))
    };

    // A count of 2 in a five-byte frame: two entries need at least eight bytes.
    match with(&|bytes| bytes[0] = 2) {
        Err(SpillError::OutOfRange { field, value: 2 }) => assert_eq!(field, "fragment_count"),
        other => panic!("a count the frame cannot hold gave {other:?}"),
    }
    // A count of 0, leaving the four bytes of the one entry unread.
    match with(&|bytes| bytes[0] = 0) {
        Err(SpillError::SectionLongerThanItsEntries { field, unread: 4 }) => {
            assert_eq!(field, "fragments");
        }
        other => panic!("bytes past the last entry gave {other:?}"),
    }
    // An allele of 65,536: three varint bytes in place of one.
    match with(&|bytes| {
        bytes.splice(3..4, [0x80, 0x80, 0x04]);
    }) {
        Err(SpillError::OutOfRange { field, value }) => {
            assert_eq!((field, value), ("fragment_allele", 65_536));
        }
        other => panic!("an allele past u16 gave {other:?}"),
    }
    // A Phred of 61, one above the largest a read allele stores.
    match with(&|bytes| bytes[4] = LARGEST_ERROR_PHRED + 1) {
        Err(SpillError::OutOfRange { field, value: 61 }) => {
            assert_eq!(field, "fragment_error_phred");
        }
        other => panic!("a Phred above the largest stored gave {other:?}"),
    }
    // Sample 1, a step of 1, in an entry with one sample.
    match with(&|bytes| bytes[1] = 1) {
        Err(SpillError::OutOfRange { field, value: 1 }) => assert_eq!(field, "fragment_sample"),
        other => panic!("a sample the entry does not have gave {other:?}"),
    }
    // A sample step of 2^32, past any sample a `u32` holds.
    match with(&|bytes| {
        let mut step = Vec::new();
        encode_u64_leb128(1 << 32, &mut step);
        bytes.splice(1..2, step);
    }) {
        Err(SpillError::OutOfRange { field, value }) => {
            assert_eq!((field, value), ("fragment_sample", 1 << 32));
        }
        other => panic!("a sample step past u32 gave {other:?}"),
    }
    // A fragment of 2^32, past any fragment a `u32` holds.
    match with(&|bytes| {
        let mut step = Vec::new();
        encode_u64_leb128(1 << 32, &mut step);
        bytes.splice(2..3, step);
    }) {
        Err(SpillError::OutOfRange { field, value }) => {
            assert_eq!((field, value), ("fragment", 1 << 32));
        }
        other => panic!("a fragment past u32 gave {other:?}"),
    }
    // And the undamaged bytes read back whole.
    let (entry, fragments) =
        read_with_fragments(framed(&body, Some(&section))).expect("the undamaged entry");
    assert!(holds_the_same_bits(&entry, &a_tiny_record()));
    assert_eq!(fragments, a_tiny_records_fragments());
}

/// **A fragments frame larger than any real one is refused**, by its compressed length before it
/// is read, and by the size it records before it is inflated.
#[test]
fn a_fragments_frame_past_the_ceilings_is_refused() {
    let body = raw_body(&a_tiny_record());

    let mut bytes = framed(&body, None);
    bytes.pop(); // the fragments' length of 0
    let too_long = max_frame_bytes(MAX_FRAGMENTS_BYTES) as u64 + 1;
    encode_u64_leb128(too_long, &mut bytes);
    match read_with_fragments(bytes) {
        Err(SpillError::OutOfRange { field, value }) => {
            assert_eq!((field, value), ("fragments_bytes", too_long));
        }
        other => panic!("an over-long fragments frame gave {other:?}"),
    }

    let too_large = vec![0u8; MAX_FRAGMENTS_BYTES + 1];
    match read_with_fragments(framed(&body, Some(&too_large))) {
        Err(SpillError::OutOfRange { field, value }) => {
            assert_eq!(
                (field, value),
                ("fragments", MAX_FRAGMENTS_BYTES as u64 + 1)
            );
        }
        other => panic!("a fragments frame claiming too much gave {other:?}"),
    }
}

/// **Pass two steps over the fragments without inflating them**: a fragments frame that is
/// damaged inside, which pass three refuses, is stepped over by pass two all the same, because
/// pass two reads only its length.
#[test]
fn pass_two_steps_over_the_fragments_without_inflating_them() {
    let body = raw_body(&a_tiny_record());
    let mut bytes = framed(&body, None);
    bytes.pop();
    let not_a_frame = [0xAAu8; 12];
    encode_u64_leb128(not_a_frame.len() as u64, &mut bytes);
    bytes.extend_from_slice(&not_a_frame);

    let entry = SpillReader::new(Cursor::new(bytes.clone()), 1)
        .next_entry()
        .expect("one entry")
        .expect("pass two steps over the frame unread");
    assert!(holds_the_same_bits(&entry, &a_tiny_record()));
    assert!(matches!(
        read_with_fragments(bytes),
        Err(SpillError::NotAFrame {
            field: "fragments",
            ..
        })
    ));
}

/// **A file cut inside the fragments frame is refused by both readers**, each naming the frame:
/// pass two's while stepping over it, pass three's while reading it.
#[test]
fn a_file_cut_inside_the_fragments_frame_is_refused_by_both_readers() {
    let whole = a_tiny_records_bytes();
    let fragments_length_at = 1 + usize::from(whole[0]);
    for cut in fragments_length_at + 1..whole.len() {
        let bytes = whole[..cut].to_vec();
        match SpillReader::new(Cursor::new(bytes.clone()), 1).next_entry() {
            Some(Err(error)) => assert_eq!(truncated_field(&error), "fragments", "{cut}"),
            other => panic!("pass two read a file cut after {cut} bytes as {other:?}"),
        }
        match read_with_fragments(bytes) {
            Err(error) => assert_eq!(truncated_field(&error), "fragments", "{cut}"),
            other => panic!("pass three read a file cut after {cut} bytes as {other:?}"),
        }
    }
}

/// **The ceiling on a frame's length is above zstd's own bound** for the largest body and the
/// largest fragments, so no frame the writer can make is refused by it.
#[test]
fn the_frame_ceilings_are_above_what_zstd_can_write() {
    for inflated in [0, 1, 1_000, MAX_BODY_BYTES, MAX_FRAGMENTS_BYTES] {
        assert!(
            max_frame_bytes(inflated) >= zstd::zstd_safe::compress_bound(inflated),
            "{inflated}"
        );
    }
}
