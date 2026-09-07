//! The one thing this module owns: **the file cannot outlive its run**, and nothing else may
//! write over it while it lives.
//!
//! Three exits are tested separately because they are three different mechanisms — a value
//! going out of scope at the end of a block, a value going out of scope because `?` returned
//! early, and a value going out of scope because a panic is unwinding past it. A guard written
//! as a call rather than as a `Drop` would pass the first and fail the other two, which are
//! exactly the runs a leftover file is most likely to come from.

use std::fs;
use std::panic::{self, AssertUnwindSafe};
use std::path::PathBuf;

use super::{SpillFile, SpillFileError};
use crate::ng::run::RunError;
use crate::ng::run::paralog_filter::{
    OnePositionSample, SpillEntry, SpillError, SpillWriter, SpilledSamples, WindowCoverage,
};
use crate::ng::types::{ContigId, GenomeRegion, Position};

/// A scratch directory inside the project, per `CLAUDE.md` — never the system temp.
fn scratch(name: &str) -> PathBuf {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_spill_file_tests")
        .join(name);
    fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

/// The output path a run would write, inside its own scratch directory.
///
/// **Any spill left by an earlier run of this test is removed first.** A test that asserts a
/// file is absent must not be able to pass, or fail, on what a previous run left behind — and
/// the failing runs are exactly the ones that leave something.
fn an_output_path(test_name: &str) -> PathBuf {
    let output = scratch(test_name).join("cohort.vcf");
    let _ = fs::remove_file(SpillFile::beside(&output).path());
    output
}

/// One record, distinguishable from the next by its position.
fn a_record_at(position: u64) -> SpillEntry {
    SpillEntry {
        contig: ContigId(2),
        position: Position(position),
        is_repeat_tract: false,
        line: format!("SL4.0ch03\t{position}\t.\tA\tG\t42.5\tPASS\tAF=0.5\tGT:AD\t0/1:5,5")
            .into_bytes(),
        samples: SpilledSamples::OnePosition(vec![
            OnePositionSample {
                window: WindowCoverage {
                    gc_fraction: 0.41,
                    mean_depth: 6.25,
                },
                ref_reads: 5,
                alt_reads: 5,
            },
            OnePositionSample {
                window: WindowCoverage {
                    gc_fraction: f32::NAN,
                    mean_depth: f32::NAN,
                },
                ref_reads: 0,
                alt_reads: 0,
            },
        ]),
    }
}

/// How long these records are once encoded, computed in memory.
///
/// **Not by writing a second `SpillFile` at the same path**, which is what an earlier version of
/// the lost-tail test did — and since two spills for one output now refuse each other, that
/// route measured nothing and truncated the file it was measuring.
fn encoded_length(records: &[SpillEntry]) -> u64 {
    let mut writer = SpillWriter::new(Vec::new());
    for record in records {
        writer.append(record).expect("a Vec sink never refuses");
    }
    writer.finish().expect("a Vec sink never refuses").len() as u64
}

/// Any `RunError`, for the early-return path. Which one it is does not matter; that a `?` on it
/// takes the spill out of scope does.
fn a_run_error() -> RunError {
    RunError::RecordNotWritten {
        locus: GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(2),
        },
        source: Box::new(std::io::Error::other("something the run could not do")),
    }
}

// ---------------------------------------------------------------------------
// Where the file is, and when it appears
// ---------------------------------------------------------------------------

#[test]
fn the_path_is_the_whole_output_path_with_a_suffix_added() {
    // Appended rather than substituted, so an output with no extension, one extension or
    // several all give one predictable answer and none of them can collide with the output.
    let directory = scratch("the_path_is_the_whole_output_path_with_a_suffix_added");
    let cases = ["cohort.vcf", "cohort", "cohort.vcf.gz"];

    for name in cases {
        let output = directory.join(name);
        assert_eq!(
            SpillFile::beside(&output).path(),
            directory.join(format!("{name}.paralog-spill.tmp")),
        );
    }
}

#[test]
fn the_file_does_not_exist_until_the_first_record_is_written() {
    let output = an_output_path("the_file_does_not_exist_until_the_first_record_is_written");
    let mut spill = SpillFile::beside(&output);

    assert!(
        !spill.path().exists(),
        "naming the spill must not create it: a run with the filter off names nothing and \
         leaves nothing"
    );

    spill.append(&a_record_at(100)).expect("the first record");

    assert!(spill.path().exists(), "the first record creates the file");
}

// ---------------------------------------------------------------------------
// The exits
// ---------------------------------------------------------------------------

#[test]
fn the_file_is_gone_when_the_run_ends_normally() {
    let output = an_output_path("the_file_is_gone_when_the_run_ends_normally");
    let path = SpillFile::beside(&output).path().to_path_buf();

    {
        let mut spill = SpillFile::beside(&output);
        spill.append(&a_record_at(100)).expect("a record");
        spill.finish_writing().expect("the flush");
        assert!(path.exists(), "the file is there while the run holds it");
    }

    assert!(!path.exists(), "the file is gone when the run ends");
}

#[test]
fn the_file_is_gone_when_the_run_returns_an_error() {
    let output = an_output_path("the_file_is_gone_when_the_run_returns_an_error");
    let path = SpillFile::beside(&output).path().to_path_buf();

    // A run that opens the spill and then fails on its way to the end, exactly as an `?` on a
    // `RunError` would leave it.
    fn a_run_that_fails(output: &std::path::Path) -> Result<(), RunError> {
        let mut spill = SpillFile::beside(output);
        spill.append(&a_record_at(100)).expect("a record");
        Err(a_run_error())?;
        unreachable!("the run returned above")
    }

    let outcome = a_run_that_fails(&output);

    assert!(outcome.is_err());
    assert!(
        !path.exists(),
        "the file survived a run that returned an error, which is one of the two exits a \
         guard written as a call would miss"
    );
}

#[test]
fn the_file_is_gone_when_a_panic_unwinds_through_the_run() {
    let output = an_output_path("the_file_is_gone_when_a_panic_unwinds_through_the_run");
    let path = SpillFile::beside(&output).path().to_path_buf();

    // **The panic hook is left alone.** Replacing it to keep this deliberate panic out of the
    // log would replace it for the whole process, and the lib binary runs several thousand
    // tests on many threads: any other test failing inside this window would report its name
    // and no message at all. One backtrace in the log is the cheaper price.
    let outcome = panic::catch_unwind(AssertUnwindSafe(|| {
        let mut spill = SpillFile::beside(&output);
        spill.append(&a_record_at(100)).expect("a record");
        panic!("a deliberate panic: the run fell over with its spill open");
    }));

    assert!(outcome.is_err(), "the panic was caught");
    assert!(
        !path.exists(),
        "the file survived a panic unwinding through the run, which is the other exit a guard \
         written as a call would miss"
    );
}

#[test]
fn the_file_is_gone_even_when_the_writer_was_never_finished() {
    let output = an_output_path("the_file_is_gone_even_when_the_writer_was_never_finished");
    let path = SpillFile::beside(&output).path().to_path_buf();

    {
        let mut spill = SpillFile::beside(&output);
        spill.append(&a_record_at(100)).expect("a record");
        // No `finish_writing`: the run ends with the writer still open.
    }

    assert!(!path.exists());
}

#[test]
fn a_spill_nothing_was_written_to_creates_nothing_and_removes_nothing() {
    let output = an_output_path("a_spill_nothing_was_written_to_creates_nothing_and_removes");
    let path = SpillFile::beside(&output).path().to_path_buf();

    // A file that is not this run's, sitting where the spill would go. A `SpillFile` that was
    // never written to must leave it exactly as it found it — its `Drop` has nothing of its own
    // to remove, and removing somebody else's file would destroy a leftover that is evidence.
    fs::write(&path, b"not this run's bytes").expect("the planted file");

    {
        let spill = SpillFile::beside(&output);
        assert_eq!(spill.entries_written(), 0);
        assert!(
            matches!(spill.read(), Err(SpillFileError::PassOneHasNotEnded { .. })),
            "a spill nothing was written to has not finished pass one"
        );
    }

    assert_eq!(
        fs::read(&path).expect("the planted file survives"),
        b"not this run's bytes",
        "a spill that was never written to removed a file it did not create"
    );
    fs::remove_file(&path).expect("the planted file is the test's to clean up");
}

// ---------------------------------------------------------------------------
// Nothing else may write over it
// ---------------------------------------------------------------------------

#[test]
fn something_already_at_the_path_is_refused_rather_than_overwritten() {
    // Either it is a leftover from a run that was killed — evidence — or it is another run
    // writing the same output right now. Truncating loses one or the other, and returns `Ok`.
    let output = an_output_path("something_already_at_the_path_is_refused");
    let path = SpillFile::beside(&output).path().to_path_buf();
    fs::write(&path, b"an earlier run's spill").expect("the leftover");

    let mut spill = SpillFile::beside(&output);

    match spill.append(&a_record_at(100)) {
        Err(SpillFileError::SomethingIsAlreadyThere { path: named }) => {
            assert_eq!(named, path);
        }
        other => panic!("expected the leftover to be refused, got {:?}", other.err()),
    }
    assert_eq!(
        fs::read(&path).expect("the leftover survives"),
        b"an earlier run's spill",
        "the leftover was destroyed by the run that refused it"
    );
    drop(spill);
    fs::remove_file(&path).expect("the leftover is the test's to clean up");
}

#[test]
fn a_second_spill_for_one_output_is_refused_rather_than_clobbering_the_first() {
    let output = an_output_path("a_second_spill_for_one_output_is_refused");
    let mut first = SpillFile::beside(&output);
    first.append(&a_record_at(100)).expect("a record");
    let bytes_before = fs::metadata(first.path()).expect("the file").len();

    let mut second = SpillFile::beside(&output);
    let outcome = second.append(&a_record_at(200));

    assert!(matches!(
        outcome,
        Err(SpillFileError::SomethingIsAlreadyThere { .. })
    ));
    drop(second);
    assert_eq!(
        fs::metadata(first.path()).expect("the file").len(),
        bytes_before,
        "the second spill truncated the first's file, or its Drop unlinked it"
    );
}

#[test]
fn appending_after_pass_one_has_ended_is_refused() {
    // Silently reopening would empty the file — the whole spill lost, reported later as a lost
    // *tail* — and the record count is already settled.
    let output = an_output_path("appending_after_pass_one_has_ended_is_refused");
    let mut spill = SpillFile::beside(&output);
    spill.append(&a_record_at(100)).expect("a record");
    spill.finish_writing().expect("the flush");
    let bytes_after_pass_one = fs::metadata(spill.path()).expect("the file").len();

    match spill.append(&a_record_at(200)) {
        Err(SpillFileError::AppendedAfterPassOneEnded { path }) => {
            assert_eq!(path, spill.path());
        }
        other => panic!("expected the append to be refused, got {:?}", other.err()),
    }

    assert_eq!(
        spill.entries_written(),
        1,
        "the refused record is not counted"
    );
    assert_eq!(
        fs::metadata(spill.path()).expect("the file").len(),
        bytes_after_pass_one,
        "the refused append emptied the file"
    );
}

#[test]
fn a_record_the_codec_refused_is_not_counted() {
    // A repeat tract marked as a biallelic SNP: spec §3.2 excludes it, and the codec says so.
    let output = an_output_path("a_record_the_codec_refused_is_not_counted");
    let mut spill = SpillFile::beside(&output);
    let impossible = SpillEntry {
        is_repeat_tract: true,
        ..a_record_at(100)
    };

    assert!(matches!(
        spill.append(&impossible),
        Err(SpillFileError::Write { .. })
    ));
    assert_eq!(spill.entries_written(), 0);

    spill.finish_writing().expect("the flush");
    assert_eq!(
        spill.read().expect("a cursor").count(),
        0,
        "a spill that counted a record it never wrote reports the wrong number"
    );
}

// ---------------------------------------------------------------------------
// Reading it back
// ---------------------------------------------------------------------------

#[test]
fn the_records_come_back_from_the_file_bit_for_bit() {
    let output = an_output_path("the_records_come_back_from_the_file_bit_for_bit");
    let mut spill = SpillFile::beside(&output);
    let written: Vec<SpillEntry> = (0..3).map(|index| a_record_at(100 + index)).collect();

    for record in &written {
        spill.append(record).expect("a record");
    }
    spill.finish_writing().expect("the flush");

    let read: Vec<SpillEntry> = spill
        .read()
        .expect("a cursor")
        .map(|entry| entry.expect("an entry"))
        .collect();

    assert_eq!(spill.entries_written(), 3);
    assert_eq!(read.len(), 3);
    for (wrote, came_back) in written.iter().zip(&read) {
        assert_eq!(wrote.position, came_back.position);
        assert_eq!(wrote.line, came_back.line);
        assert!(
            came_back
                .samples
                .window(1)
                .expect("two samples went in")
                .gc_fraction
                .is_nan(),
            "the absent sample came back as a number after a round trip through the filesystem"
        );
    }
}

#[test]
fn the_file_can_be_read_twice() {
    // Pass two scores the entries and pass three writes them; each needs its own cursor over
    // the same file, from the beginning.
    let output = an_output_path("the_file_can_be_read_twice");
    let mut spill = SpillFile::beside(&output);
    spill.append(&a_record_at(100)).expect("a record");
    spill.append(&a_record_at(101)).expect("a record");
    spill.finish_writing().expect("the flush");

    let first_pass: Vec<u64> = spill
        .read()
        .expect("a cursor")
        .map(|entry| entry.expect("an entry").position.get())
        .collect();
    let second_pass: Vec<u64> = spill
        .read()
        .expect("a second cursor")
        .map(|entry| entry.expect("an entry").position.get())
        .collect();

    assert_eq!(first_pass, vec![100, 101]);
    assert_eq!(second_pass, first_pass, "the second cursor starts over");
}

#[test]
fn a_run_that_called_nothing_reads_an_empty_spill_rather_than_a_missing_file() {
    let output = an_output_path("a_run_that_called_nothing_reads_an_empty_spill");
    let mut spill = SpillFile::beside(&output);

    spill.finish_writing().expect("the flush ends pass one");

    let read: Vec<SpillEntry> = spill
        .read()
        .expect("a cursor over an empty spill")
        .map(|entry| entry.expect("no entry can fail"))
        .collect();

    assert!(read.is_empty());
    assert_eq!(spill.entries_written(), 0);
}

#[test]
fn finishing_pass_one_twice_leaves_the_records_where_they_are() {
    // The second call must not reopen the file: reopening would empty it, and the whole spill
    // would be lost as what looks later like a missing tail.
    let output = an_output_path("finishing_pass_one_twice_leaves_the_records_where_they_are");
    let mut spill = SpillFile::beside(&output);
    spill.append(&a_record_at(100)).expect("a record");
    spill.finish_writing().expect("the first flush");
    let bytes = fs::metadata(spill.path()).expect("the file").len();

    spill.finish_writing().expect("the second flush");

    assert_eq!(fs::metadata(spill.path()).expect("the file").len(), bytes);
    assert_eq!(spill.read().expect("a cursor").count(), 1);
}

#[test]
fn reading_before_pass_one_has_ended_is_refused() {
    let output = an_output_path("reading_before_pass_one_has_ended_is_refused");
    let mut spill = SpillFile::beside(&output);

    // Never opened: the answer is "pass one has not finished", not "the file is missing".
    match spill.read() {
        Err(SpillFileError::PassOneHasNotEnded { path }) => assert_eq!(path, spill.path()),
        other => panic!("expected the read to be refused, got {:?}", other.err()),
    }

    // Open but not flushed: at 64 KiB of buffering, reading here usually sees nothing at all.
    spill.append(&a_record_at(100)).expect("a record");
    match spill.read() {
        Err(SpillFileError::PassOneHasNotEnded { path }) => assert_eq!(path, spill.path()),
        other => panic!("expected the read to be refused, got {:?}", other.err()),
    }
}

#[test]
fn a_spill_whose_tail_was_lost_is_refused_rather_than_read_short() {
    // The failure this catches: a file cut on an entry boundary is byte for byte the prefix of
    // a complete one, so without the writer's count the reader cannot tell them apart — and
    // pass three would write a VCF short by its last records, with no panic and no message.
    let output = an_output_path("a_spill_whose_tail_was_lost_is_refused_rather_than_read_short");
    let three: Vec<SpillEntry> = (0..3).map(|index| a_record_at(100 + index)).collect();
    let mut spill = SpillFile::beside(&output);
    for record in &three {
        spill.append(record).expect("a record");
    }
    spill.finish_writing().expect("the flush");

    // The boundary is computed in memory. Measuring it by writing a second `SpillFile` at the
    // same path would truncate the file being measured — which is what an earlier version of
    // this test did, leaving its `set_len` shortening the file by nothing at all.
    let boundary = encoded_length(&three[..2]);
    let whole = fs::metadata(spill.path()).expect("the file").len();
    assert!(
        boundary < whole,
        "two records must encode shorter than three: {boundary} against {whole}"
    );

    fs::OpenOptions::new()
        .write(true)
        .open(spill.path())
        .expect("the file")
        .set_len(boundary)
        .expect("the truncation");

    let outcome: Vec<Result<SpillEntry, SpillError>> = spill.read().expect("a cursor").collect();

    assert_eq!(outcome.len(), 3, "two entries, then the refusal");
    assert!(outcome[0].is_ok() && outcome[1].is_ok());
    match &outcome[2] {
        Err(SpillError::TheWrongNumberOfRecords { expected, read }) => {
            assert_eq!(*expected, 3);
            assert_eq!(*read, 2);
        }
        other => panic!("expected the short file to be refused, got {other:?}"),
    }
}

#[test]
fn a_failure_off_the_reader_can_be_given_the_file_it_happened_on() {
    // The reader is generic over its stream and carries no path; spec §5 asks that a read
    // failure name the file, and this is the join.
    let output = an_output_path("a_failure_off_the_reader_can_be_given_the_file");
    let spill = SpillFile::beside(&output);

    let named = spill.naming_this_file(SpillError::Truncated { field: "line" });

    match named {
        SpillFileError::Read { path, source } => {
            assert_eq!(path, spill.path());
            assert!(matches!(source, SpillError::Truncated { field: "line" }));
        }
        other => panic!("expected a read failure naming the file, got {other:?}"),
    }
}

#[test]
fn a_spill_that_cannot_be_created_names_the_path_and_says_what_the_filesystem_said() {
    let output = scratch("a_spill_that_cannot_be_created")
        .join("no-such-directory")
        .join("cohort.vcf");
    let mut spill = SpillFile::beside(&output);

    match spill.append(&a_record_at(100)) {
        Err(SpillFileError::Create { path, source }) => {
            assert_eq!(path, spill.path());
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("expected the create to fail, got {:?}", other.err()),
    }
}
