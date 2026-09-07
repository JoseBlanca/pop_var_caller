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
//! The tests that corrupt bytes by index all work on [`a_tiny_record`], whose whole encoding is
//! written out in [`the_encoding_matches_the_layout_the_spec_fixes`]; the offsets they use are
//! read off that list rather than counted from a fixture that could change under them.

use std::io::Cursor;

use proptest::prelude::*;

use super::{
    GenericLocusSample, RepeatTractSample, SpillEntry, SpillError, SpillReader, SpillWriter,
    SpilledSamples,
};
use crate::ng::run::paralog_filter::WindowCoverage;
use crate::ng::types::{ContigId, Position};
use crate::psp::varint::encode_u64_leb128;

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
) -> GenericLocusSample {
    GenericLocusSample {
        window: a_window(gc_fraction, mean_depth),
        ref_reads,
        alt_reads,
    }
}

/// A one-position sample with no usable window — the `NaN` pair — but reads all the same.
fn a_sample_without_a_window(ref_reads: u32, alt_reads: u32) -> GenericLocusSample {
    GenericLocusSample {
        window: no_window(),
        ref_reads,
        alt_reads,
    }
}

/// A wide-locus sample, which carries a window and nothing else.
fn a_tract_sample(gc_fraction: f32, mean_depth: f32) -> RepeatTractSample {
    RepeatTractSample {
        window: a_window(gc_fraction, mean_depth),
    }
}

/// A wide-locus sample with no usable window — the case that must survive by bit pattern.
fn a_tract_sample_without_a_window() -> RepeatTractSample {
    RepeatTractSample {
        window: no_window(),
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
        samples: SpilledSamples::GenericLocus(vec![
            a_sample_with_a_window(0.41, 6.25, 5, 5),
            a_sample_without_a_window(8, 0),
        ]),
    }
}

/// A repeat tract two samples covered, one of them with no window. **Its samples carry no read
/// counts at all** — the shape spec §3.2 gives a locus spanning more than one base.
fn a_tract_two_samples_covered() -> SpillEntry {
    SpillEntry {
        contig: ContigId(3),
        position: Position(2_000_000),
        is_repeat_tract: true,
        line: b"SL4.0ch04\t2000000\t.\tATAT\tAT\t31.0\tPASS\tAN=4;DP=12;STR;RU=AT;PERIOD=2\tGT:GQ:DP:AD:REPCN\t1/1:30:6:0,6:2\t0/0:99:6:6,0:4"
            .to_vec(),
        samples: SpilledSamples::RepeatTract(vec![
            a_tract_sample(0.38, 5.75),
            a_tract_sample_without_a_window(),
        ]),
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
        samples: SpilledSamples::GenericLocus(vec![a_sample_with_a_window(0.5, 2.0, 3, 130)]),
    }
}

/// Offsets into [`a_tiny_record`]'s encoding, read off the byte list in
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

/// The samples, compared by bit pattern — **and a record that changed shape is not the same
/// record**, so the two variants never compare equal to each other.
fn samples_hold_the_same_bits(left: &SpilledSamples, right: &SpilledSamples) -> bool {
    match (left, right) {
        (SpilledSamples::GenericLocus(one), SpilledSamples::GenericLocus(other)) => {
            one.len() == other.len()
                && one
                    .iter()
                    .zip(other)
                    .all(|(one, other)| one_position_sample_holds_the_same_bits(one, other))
        }
        (SpilledSamples::RepeatTract(one), SpilledSamples::RepeatTract(other)) => {
            one.len() == other.len()
                && one
                    .iter()
                    .zip(other)
                    .all(|(one, other)| wide_sample_holds_the_same_bits(one, other))
        }
        (SpilledSamples::GenericLocus(_) | SpilledSamples::RepeatTract(_), _) => false,
    }
}

/// One one-position sample's four numbers, by bit pattern. Destructured for the same reason.
fn one_position_sample_holds_the_same_bits(
    one: &GenericLocusSample,
    other: &GenericLocusSample,
) -> bool {
    let GenericLocusSample {
        window,
        ref_reads,
        alt_reads,
    } = one;
    let GenericLocusSample {
        window: other_window,
        ref_reads: other_ref_reads,
        alt_reads: other_alt_reads,
    } = other;

    window_holds_the_same_bits(window, other_window)
        && ref_reads == other_ref_reads
        && alt_reads == other_alt_reads
}

/// One wide sample's window, by bit pattern. Destructured for the same reason.
fn wide_sample_holds_the_same_bits(one: &RepeatTractSample, other: &RepeatTractSample) -> bool {
    let RepeatTractSample { window } = one;
    let RepeatTractSample {
        window: other_window,
    } = other;

    window_holds_the_same_bits(window, other_window)
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

/// The bytes the writer produces for these entries.
fn encoded_bytes(entries: &[SpillEntry]) -> Vec<u8> {
    let mut writer = SpillWriter::new(Vec::new());
    for entry in entries {
        writer
            .append(entry)
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
                samples: SpilledSamples::RepeatTract(vec![
                    a_tract_sample(0.41, 6.25),
                    a_tract_sample_without_a_window(),
                ]),
                ..base.clone()
            },
        ),
        (
            "the samples' shape, with the same windows in them",
            SpillEntry {
                samples: SpilledSamples::RepeatTract(vec![
                    a_tract_sample(0.41, 6.25),
                    a_tract_sample_without_a_window(),
                ]),
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
                samples: SpilledSamples::GenericLocus(vec![a_sample_with_a_window(
                    0.41, 6.25, 5, 5,
                )]),
                ..base.clone()
            },
        ),
        (
            "gc_fraction",
            SpillEntry {
                samples: SpilledSamples::GenericLocus(vec![
                    a_sample_with_a_window(0.42, 6.25, 5, 5),
                    a_sample_without_a_window(8, 0),
                ]),
                ..base.clone()
            },
        ),
        (
            "mean_depth",
            SpillEntry {
                samples: SpilledSamples::GenericLocus(vec![
                    a_sample_with_a_window(0.41, 6.26, 5, 5),
                    a_sample_without_a_window(8, 0),
                ]),
                ..base.clone()
            },
        ),
        (
            "an absent window against a zeroed one",
            SpillEntry {
                samples: SpilledSamples::GenericLocus(vec![
                    a_sample_with_a_window(0.41, 6.25, 5, 5),
                    a_sample_with_a_window(0.0, 0.0, 8, 0),
                ]),
                ..base.clone()
            },
        ),
        (
            "ref_reads",
            SpillEntry {
                samples: SpilledSamples::GenericLocus(vec![
                    a_sample_with_a_window(0.41, 6.25, 6, 5),
                    a_sample_without_a_window(8, 0),
                ]),
                ..base.clone()
            },
        ),
        (
            "alt_reads",
            SpillEntry {
                samples: SpilledSamples::GenericLocus(vec![
                    a_sample_with_a_window(0.41, 6.25, 5, 6),
                    a_sample_without_a_window(8, 0),
                ]),
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

    let absent = came_back.samples.window(1).expect("two samples went in");
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
        samples: SpilledSamples::GenericLocus(vec![GenericLocusSample {
            window: WindowCoverage {
                gc_fraction: odd_nan,
                mean_depth: f32::NAN,
            },
            ref_reads: 0,
            alt_reads: 0,
        }]),
        ..a_snp_two_samples_covered()
    };

    let came_back = assert_round_trips(&entry);

    assert_eq!(
        came_back
            .samples
            .window(0)
            .expect("one sample went in")
            .gc_fraction
            .to_bits(),
        0x7F80_0001,
        "the float came back through arithmetic rather than as its bits"
    );
}

#[test]
fn a_wide_locus_carries_an_absent_window_through_by_its_bits_too() {
    // The wide row has one field and it is the one that can be absent, so nothing else in the
    // codec stands between it and the file.
    let odd_nan = f32::from_bits(0x7F80_0001);
    let entry = SpillEntry {
        samples: SpilledSamples::RepeatTract(vec![RepeatTractSample {
            window: WindowCoverage {
                gc_fraction: odd_nan,
                mean_depth: f32::NAN,
            },
        }]),
        ..a_tract_two_samples_covered()
    };

    let came_back = assert_round_trips(&entry);

    assert_eq!(
        came_back
            .samples
            .window(0)
            .expect("one sample went in")
            .gc_fraction
            .to_bits(),
        0x7F80_0001
    );
}

#[test]
fn a_record_with_no_samples_at_all_round_trips() {
    let entry = SpillEntry {
        samples: SpilledSamples::GenericLocus(Vec::new()),
        ..a_snp_two_samples_covered()
    };
    assert!(assert_round_trips(&entry).samples.is_empty());

    let wide = SpillEntry {
        samples: SpilledSamples::RepeatTract(Vec::new()),
        ..a_tract_two_samples_covered()
    };
    assert!(assert_round_trips(&wide).samples.is_empty());
}

#[test]
fn a_tract_comes_back_carrying_its_windows_and_no_read_counts() {
    // The shape itself round-trips: a wide entry must not come back as a one-position entry
    // with zeros in it, which is the confusion the two variants exist to prevent.
    let entry = a_tract_two_samples_covered();

    let came_back = assert_round_trips(&entry);

    assert!(
        matches!(came_back.samples, SpilledSamples::RepeatTract(_)),
        "a wide locus came back as {:?}",
        came_back.samples
    );
    assert_eq!(came_back.samples.len(), 2);
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
fn a_repeat_tract_round_trips_with_its_tract_flag_and_no_read_counts_at_all() {
    let entry = SpillEntry {
        contig: ContigId(0),
        position: Position(1),
        is_repeat_tract: true,
        line: b"SL4.0ch00\t1\t.\tATATAT\tATATATAT\t61.0\tPASS\tAF=0.25\tGT:GQ:AD\t0/1:20:4,4"
            .to_vec(),
        samples: SpilledSamples::RepeatTract(vec![
            a_tract_sample(0.19, 3.5),
            a_tract_sample_without_a_window(),
        ]),
    };
    let came_back = assert_round_trips(&entry);
    assert!(came_back.is_repeat_tract);
    assert!(matches!(came_back.samples, SpilledSamples::RepeatTract(_)));
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
        samples: SpilledSamples::GenericLocus(vec![
            a_sample_with_a_window(f32::MAX, f32::MIN_POSITIVE, u32::MAX, u32::MAX),
            a_sample_with_a_window(-0.0, f32::INFINITY, 0, 0),
            a_sample_with_a_window(f32::NEG_INFINITY, 0.0, 1, 1),
        ]),
        ..a_snp_two_samples_covered()
    };
    let came_back = assert_round_trips(&entry);
    assert_eq!(
        came_back
            .samples
            .window(1)
            .expect("three samples went in")
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
        samples: SpilledSamples::GenericLocus(
            (0..1000u32)
                .map(|index| {
                    if index % 3 == 0 {
                        a_sample_without_a_window(index, 0)
                    } else {
                        a_sample_with_a_window(0.4, index as f32, index, index + 1)
                    }
                })
                .collect(),
        ),
        ..a_tiny_record()
    };
    let came_back = assert_round_trips(&entry);
    assert_eq!(came_back.samples.len(), 1000);
    let SpilledSamples::GenericLocus(samples) = &came_back.samples else {
        panic!("a one-position entry came back wide");
    };
    assert_eq!(samples[999].ref_reads, 999);
}

#[test]
fn a_cohort_of_a_thousand_at_a_wide_locus_round_trips() {
    // Same crossing of the varint boundary, on the shape that carries eight bytes a sample
    // instead of ten or more.
    let entry = SpillEntry {
        samples: SpilledSamples::RepeatTract(
            (0..1000u32)
                .map(|index| {
                    if index % 3 == 0 {
                        a_tract_sample_without_a_window()
                    } else {
                        a_tract_sample(0.4, index as f32)
                    }
                })
                .collect(),
        ),
        ..a_tract_two_samples_covered()
    };
    let came_back = assert_round_trips(&entry);
    assert_eq!(came_back.samples.len(), 1000);
}

#[test]
fn a_stream_of_records_comes_back_in_the_order_it_was_written() {
    let first = a_snp_two_samples_covered();
    let second = SpillEntry {
        position: Position(1_000_001),
        is_repeat_tract: true,
        line: b"SL4.0ch04\t1000001\t.\tAT\tATAT\t9.0\tPASS\t.\tGT\t0/1\t0/0".to_vec(),
        samples: SpilledSamples::RepeatTract(vec![
            a_tract_sample(0.41, 6.25),
            a_tract_sample_without_a_window(),
        ]),
        ..a_snp_two_samples_covered()
    };
    let third = SpillEntry {
        contig: ContigId(4),
        position: Position(7),
        samples: SpilledSamples::GenericLocus(Vec::new()),
        ..a_snp_two_samples_covered()
    };

    let read = round_trip(&[first.clone(), second.clone(), third.clone()]);

    assert_eq!(read.len(), 3);
    assert!(holds_the_same_bits(&first, &read[0]));
    assert!(holds_the_same_bits(&second, &read[1]));
    assert!(holds_the_same_bits(&third, &read[2]));
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
    // Pins the field order and each field's width. A codec whose encoder and decoder agree
    // with each other but not with §3.4 passes every round trip above and fails this.
    assert_eq!(
        encoded_bytes(&[a_tiny_record()]),
        vec![
            0x01, // contig 1, one varint byte                          — offset 0
            0xAC, 0x02, // position 300, two varint bytes               — offsets 1, 2
            0x00, // is_repeat_tract = false, so generic-locus rows      — offset 3
            0x02, // the line is two bytes long                         — offset 4
            b'A', b'B', // the line                                     — offsets 5, 6
            0x01, // one sample                                         — offset 7
            0x00, 0x00, 0x00, 0x3F, // gc_fraction 0.5, an f32's bits little-endian
            0x00, 0x00, 0x00, 0x40, // mean_depth 2.0
            0x03, // ref_reads 3
            0x82, 0x01, // alt_reads 130, two varint bytes
        ]
    );
}

// ---------------------------------------------------------------------------
// Files that are cut, and files that are corrupt
// ---------------------------------------------------------------------------

#[test]
fn a_file_cut_anywhere_inside_a_record_names_the_field_the_bytes_ran_out_in() {
    // A file that fails is not enough: the field name is the whole diagnostic a spill failure
    // carries, so every cut has to blame the right one. The offsets are a_tiny_record's, from
    // the byte list above.
    let expected: &[(usize, &str)] = &[
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

    let whole = encoded_bytes(&[a_tiny_record()]);
    assert_eq!(
        whole.len(),
        19,
        "the byte list moved and the offsets below move with it"
    );

    for &(cut, field) in expected {
        let mut reader = SpillReader::new(Cursor::new(whole[..cut].to_vec()), 1);
        match reader.next_entry() {
            Some(Err(error)) => assert_eq!(
                truncated_field(&error),
                field,
                "a file cut after {cut} bytes blamed {} rather than {field}",
                truncated_field(&error)
            ),
            other => panic!("a file cut after {cut} bytes gave {other:?} rather than failing"),
        }
    }
}

#[test]
fn a_failure_inside_the_samples_names_which_sample() {
    // At a cohort of three thousand a message that names only the field points at three
    // thousand places at once.
    let entry = SpillEntry {
        samples: SpilledSamples::GenericLocus(vec![a_sample_with_a_window(0.5, 2.0, 1, 1); 3]),
        ..a_tiny_record()
    };
    let whole = encoded_bytes(std::slice::from_ref(&entry));
    let head = encoded_bytes(&[SpillEntry {
        samples: SpilledSamples::GenericLocus(Vec::new()),
        ..entry
    }])
    .len();
    let bytes_per_sample = (whole.len() - head) / 3;
    let two_bytes_into_the_second_sample = head + bytes_per_sample + 2;

    let mut reader = SpillReader::new(
        Cursor::new(whole[..two_bytes_into_the_second_sample].to_vec()),
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
    let mut bytes = encoded_bytes(&[a_tiny_record()]);
    bytes[TRACT_FLAG_OFFSET_IN_TINY_RECORD] = 2;
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    match reader.next_entry() {
        Some(Err(SpillError::NotABoolean { field, byte })) => {
            assert_eq!(field, "is_repeat_tract");
            assert_eq!(byte, 2);
        }
        other => panic!("expected a refused flag byte, got {other:?}"),
    }
}

#[test]
fn an_entry_whose_sample_shape_disagrees_with_its_tract_flag_is_refused_by_the_writer() {
    // **There is no reader-side counterpart, and that is by design.** The row shape is read off
    // the tract flag, so the file cannot hold the two disagreeing — which is why the check lives
    // where the caller that built the entry can still be named. Both directions, because either
    // one silently changes what the scorer is handed.
    let a_tract_carrying_read_counts = SpillEntry {
        is_repeat_tract: true,
        ..a_tiny_record()
    };
    let a_generic_locus_carrying_none = SpillEntry {
        is_repeat_tract: false,
        ..a_tract_two_samples_covered()
    };

    for (entry, flag) in [
        (a_tract_carrying_read_counts, true),
        (a_generic_locus_carrying_none, false),
    ] {
        let mut writer = SpillWriter::new(Vec::new());

        match writer.append(&entry) {
            Err(SpillError::SampleShapeDisagreesWithTheTractFlag {
                is_repeat_tract, ..
            }) => assert_eq!(is_repeat_tract, flag),
            other => panic!("expected the writer to refuse the record, got {other:?}"),
        }
        assert_eq!(
            writer.entries_written(),
            0,
            "a record the writer refused is not counted"
        );
    }
}

#[test]
fn a_contig_too_large_for_the_field_is_refused_rather_than_wrapped() {
    // A varint holding 2^32 — one past the largest contig index there can be. Written by
    // hand because the encoder cannot produce it.
    let mut bytes = Vec::new();
    encode_u64_leb128(1 << 32, &mut bytes);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

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
    let mut bytes = Vec::new();
    encode_u64_leb128(0, &mut bytes); // contig
    encode_u64_leb128(1, &mut bytes); // position
    bytes.push(0); // not a tract, so the samples carry read counts
    encode_u64_leb128(0, &mut bytes); // an empty line
    encode_u64_leb128(1, &mut bytes); // one sample
    bytes.extend_from_slice(&0.5f32.to_bits().to_le_bytes());
    bytes.extend_from_slice(&1.0f32.to_bits().to_le_bytes());
    encode_u64_leb128(1 << 32, &mut bytes); // ref_reads, one past the field
    encode_u64_leb128(0, &mut bytes);

    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

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
    // value.
    let bytes = vec![0xFF; 11];
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OverlongVarint { field })) => assert_eq!(field, "contig"),
        other => panic!("expected an over-long varint, got {other:?}"),
    }
}

#[test]
fn a_ten_byte_varint_holding_more_than_a_u64_is_refused_rather_than_truncated() {
    // Ten bytes is a legal LEB128 width, but the tenth byte contributes `data << 63`, so all
    // but its lowest bit falls off the top of the `u64`. The psp primitive returns `Ok` with
    // the high bits dropped; this codec refuses corruption rather than absorbing it.
    let mut bytes = vec![0x80; 9];
    bytes.push(0x7F);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OverlongVarint { field })) => assert_eq!(field, "contig"),
        other => panic!("expected a refused varint, got {other:?}"),
    }
}

#[test]
fn a_line_longer_than_the_ceiling_is_refused_before_it_is_read() {
    // Without the ceiling the length is believed and every byte the file still holds is
    // appended before the short-read check fires — on a spill the size of a genome's VCF,
    // that is the rest of the file in one `Vec`.
    let mut bytes = Vec::new();
    encode_u64_leb128(0, &mut bytes); // contig
    encode_u64_leb128(1, &mut bytes); // position
    bytes.push(0); // not a tract, so generic-locus rows
    encode_u64_leb128(u64::from(u32::MAX), &mut bytes); // a line of four billion bytes
    bytes.extend_from_slice(&vec![b'X'; 4096]);

    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

    match reader.next_entry() {
        Some(Err(SpillError::OutOfRange { field, value })) => {
            assert_eq!(field, "line_length");
            assert_eq!(value, u64::from(u32::MAX));
        }
        other => panic!("expected a refused line length, got {other:?}"),
    }
}

#[test]
fn a_sample_count_larger_than_the_file_is_an_error_and_not_an_allocation() {
    // Nine hundred thousand samples — under the ceiling, so the decoder believes the count —
    // and a file that ends immediately. It must fail on the bytes it does not have rather
    // than reserve for the count it was told: without the reservation cap this allocates
    // 900,000 x 16 bytes before reading one.
    let mut bytes = vec![
        0x00, // contig 0
        0x01, // position 1
        0x00, // not a tract, so generic-locus rows
        0x00, // an empty line
    ];
    encode_u64_leb128(900_000, &mut bytes);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

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
    // loop keeps decoding until the file runs out, so the memory held is the rest of the file
    // rather than one entry. The bytes below are a header followed by a count of ten million
    // and enough sample bytes that a decoder without the ceiling would build them all.
    let mut bytes = vec![
        0x00, // contig 0
        0x01, // position 1
        0x00, // not a tract, so generic-locus rows
        0x00, // an empty line
    ];
    encode_u64_leb128(10_000_000, &mut bytes);
    bytes.extend_from_slice(&vec![0u8; 4096]);
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

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
fn a_source_that_fails_mid_float_says_the_read_failed_and_not_that_the_file_ended() {
    // A device failure and a file that ends are different things, and the field's four bytes
    // are the one place the codec has to tell them apart by hand.
    struct FailsAfterTheHeader {
        remaining: Vec<u8>,
    }

    impl std::io::Read for FailsAfterTheHeader {
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

    let head = vec![
        0x00, // contig 0
        0x01, // position 1
        0x00, // not a tract, so generic-locus rows
        0x00, // an empty line
        0x01, // one sample
    ];
    let source = std::io::BufReader::with_capacity(1, FailsAfterTheHeader { remaining: head });
    let mut reader = SpillReader::new(source, 1);

    match reader.next_entry() {
        Some(Err(SpillError::InSample { index, source })) => {
            assert_eq!(index, 0);
            match *source {
                SpillError::Read { field, .. } => assert_eq!(field, "gc_fraction"),
                other => panic!("expected a read failure, got {other:?}"),
            }
        }
        other => panic!("expected a read failure inside the sample, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// The reader stops, and the writer reports its sink
// ---------------------------------------------------------------------------

#[test]
fn the_reader_stops_after_a_decode_error() {
    // A stream that failed once is not at an entry boundary any more, so bytes from the middle
    // of an entry would decode — plausibly — into a record that was never written.
    let mut bytes = vec![0x01, 0x01, 0x02]; // contig, position, then a flag byte of 2
    bytes.extend_from_slice(&encoded_bytes(&[a_tiny_record()]));
    let mut reader = SpillReader::new(Cursor::new(bytes), 1);

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
        writer.append(&a_tiny_record()),
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
fn any_entry() -> impl Strategy<Value = SpillEntry> {
    (
        any::<u32>(),
        any::<u64>(),
        any::<bool>(),
        prop::collection::vec(any::<u8>(), 0..400),
        prop::collection::vec(
            (any_float(), any_float(), any::<u32>(), any::<u32>()),
            0..300,
        ),
    )
        .prop_map(|(contig, position, is_repeat_tract, line, samples)| {
            let samples = if is_repeat_tract {
                SpilledSamples::RepeatTract(
                    samples
                        .into_iter()
                        .map(|(gc, depth, _, _)| a_tract_sample(gc, depth))
                        .collect(),
                )
            } else {
                SpilledSamples::GenericLocus(
                    samples
                        .into_iter()
                        .map(|(gc, depth, ref_reads, alt_reads)| {
                            a_sample_with_a_window(gc, depth, ref_reads, alt_reads)
                        })
                        .collect(),
                )
            };
            SpillEntry {
                contig: ContigId(contig),
                position: Position(position),
                is_repeat_tract,
                line,
                samples,
            }
        })
}

proptest! {
    #[test]
    fn any_stream_of_entries_comes_back_bit_for_bit(
        entries in prop::collection::vec(any_entry(), 1..8)
    ) {
        let read = round_trip(&entries);
        prop_assert_eq!(read.len(), entries.len());
        for (wrote, came_back) in entries.iter().zip(&read) {
            prop_assert!(holds_the_same_bits(wrote, came_back));
        }
    }
}
