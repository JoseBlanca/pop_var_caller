//! The windows a calling run read at each locus, written to a file so that something outside the
//! run can check them.
//!
//! **Why a run has to be asked rather than read.** The window-coverage measurement changes no VCF
//! byte — that is the standing oracle for the whole of this work — so a run whose windows were all
//! absent, or all one cover behind, produces exactly the output a correct one produces. The failure
//! plan step C3 exists to prevent is precisely of that shape: a cover that stops at its region's
//! last base leaves that region's last centres unfinalised, and every sample reads as absent at the
//! loci nearest each region boundary, silently. Nothing in the run's own output can show it.
//!
//! So a run set up for the comparison writes down what it read, and
//! [`examples/ng_window_coverage_probe.rs`](../../../../../examples/ng_window_coverage_probe.rs)
//! walks the same stored files end to end — one pass, no covers, no eviction — and checks the two
//! agree bit for bit (spec
//! [`window_coverage.md`](../../../../../doc/devel/ng/spec/window_coverage.md) §10).
//!
//! **Both halves of the row format live here**, [`write_the_row`] beside [`read_a_row`], because
//! the reader is in another crate target: a column reordered on one side alone would still parse
//! on the other, and the comparison would then report every locus as one the run had no window
//! for — which is the shape of a *pass* in this measurement, not of a failure.
//!
//! **One row a sample per record the run writes, not per locus the merge builds.** A locus that
//! establishes no variant becomes no record and leaves no row, so it is outside what the
//! comparison checks. That is the right set — the record is where the filter will read the pair,
//! and spec §10 asks for "every written record" — but a reader chasing a boundary failure should
//! know the built-but-unwritten loci are not in the file.
//!
//! **It is off unless [`PATH_VARIABLE`] names a file**, and off it costs one already-resolved
//! `OnceLock` read per written record. It is not a debugging aid left lying about: it is the only way
//! the look-ahead's failure is observable, and deleting it would leave that step proved by nothing.
//!
//! **The file is a process's, resolved once.** Two calling runs inside one process would interleave
//! into it, and a row carries a sample index but no run identity, so the comparison would silently
//! read one run's windows as the other's. Nothing does that today — a run is a command — and a
//! caller that wants to must give each run its own file.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

use crate::ng::types::{ContigId, GenomePosition, Position};
use crate::ng::window_coverage::WindowCoverage;

/// The environment variable naming the file to write. Anything else about the run is unchanged.
pub const PATH_VARIABLE: &str = "NG_WINDOW_COVERAGE_FILE";

/// One row: which sample, which locus, and what that sample's window there was — or that it had
/// none.
///
/// **`None` and an absent window are two different rows, and no run writes the first.** `None` is
/// a run holding nothing at that position at all; an absent window is one that was finalised and
/// reports nothing. A cohort locus collapses the two before a record is written, so what reaches
/// this type from a run is always `Some` — measured on the tomato slice, 0 of 13,866 rows carry
/// the `-` form. The distinction is the format's, kept because reading it back is what a
/// comparison does.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecordedWindow {
    /// The sample's index in the run's sample order.
    pub sample: usize,
    /// The locus's first base.
    pub at: GenomePosition,
    /// What the run's window was there, or `None` where it held none.
    pub window: Option<WindowCoverage>,
}

/// The open file, or `None` where the variable was not set. Resolved once per process.
///
/// **A `Mutex` around it rather than a bare handle**, because a `static` has to be `Sync` and the
/// writer is reached through a shared reference. The record loop is one thread, so it is never
/// contended; it is taken once per written record, beside the VCF write that record already
/// causes.
fn sink() -> Option<&'static Mutex<BufWriter<File>>> {
    static SINK: OnceLock<Option<Mutex<BufWriter<File>>>> = OnceLock::new();
    SINK.get_or_init(|| {
        let path = std::env::var_os(PATH_VARIABLE)?;
        let file = File::create(&path).unwrap_or_else(|why| {
            panic!(
                "{PATH_VARIABLE} names {} and it could not be created: {why}",
                std::path::Path::new(&path).display()
            )
        });
        Some(Mutex::new(BufWriter::new(file)))
    })
    .as_ref()
}

/// One row as it is written: sample, contig, position, then the two fields' **bit patterns**, or
/// `-` in both where the run had no window.
///
/// Bits rather than digits because an absent window is a pair of `NaN`s and a printed `NaN` cannot
/// be told from another one (spec §6 trap 4), and because a decimal rendering would compare equal
/// where two values differ in the last bit.
pub fn write_the_row(row: RecordedWindow, into: &mut String) {
    use core::fmt::Write as _;

    let (gc, mean_depth) = match row.window {
        None => (
            ABSENT_FROM_THE_RUN.to_owned(),
            ABSENT_FROM_THE_RUN.to_owned(),
        ),
        Some(window) => (
            window.gc_fraction.to_bits().to_string(),
            window.mean_depth.to_bits().to_string(),
        ),
    };
    writeln!(
        into,
        "{}\t{}\t{}\t{gc}\t{mean_depth}",
        row.sample,
        row.at.contig.get(),
        row.at.position.get(),
    )
    .expect("writing into a String cannot fail");
}

/// [`write_the_row`]'s inverse, for the comparison in the other crate target.
///
/// # Errors
///
/// A row that is not five tab-separated fields, or whose numbers will not parse — naming the
/// column, because a reader looking at a malformed file needs to know which one.
pub fn read_a_row(line: &str) -> Result<RecordedWindow, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    let [sample, contig, position, gc, mean_depth] = fields[..] else {
        return Err(format!(
            "expected five tab-separated fields, got {}: {line:?}",
            fields.len()
        ));
    };
    let number = |column: &str, field: &str| -> Result<u64, String> {
        field
            .parse()
            .map_err(|_| format!("{column} is not a number: {field:?}"))
    };
    let window = match (gc, mean_depth) {
        (ABSENT_FROM_THE_RUN, ABSENT_FROM_THE_RUN) => None,
        _ => Some(WindowCoverage {
            // `as u32` cannot lose anything a `write_the_row` produced: both fields are `f32`
            // bit patterns, so they are below `u32::MAX` by construction. A hand-edited file with
            // a larger number is a malformed file, and it becomes some other float rather than an
            // error — which the round-trip test below is what stands against.
            gc_fraction: f32::from_bits(number("the GC bits", gc)? as u32),
            mean_depth: f32::from_bits(number("the mean-depth bits", mean_depth)? as u32),
        }),
    };
    Ok(RecordedWindow {
        sample: number("the sample", sample)? as usize,
        at: GenomePosition {
            contig: ContigId(number("the contig", contig)? as u32),
            position: Position(number("the position", position)?),
        },
        window,
    })
}

/// What both halves write and read where the run held no window at all.
const ABSENT_FROM_THE_RUN: &str = "-";

/// One locus's rows: one per entry of `windows_of_every_sample`, which is dense over the run's
/// samples and in the run's own sample order.
///
/// **A sample with no window is written down rather than omitted**, as an *absent* pair — two
/// `NaN`s. That is what "no window" looks like by the time it reaches here: the merge's own
/// three-way distinction collapses on the way to the record
/// ([`CohortObservation`](super::build::CohortObservation)'s own `window_coverage`), so **no run
/// writes the row's `-` form**, which survives on the reading side alone. What the comparison
/// reads is a row present against a row missing, and a row is missing only where a locus became
/// no record.
///
/// **Split from the write below because the write is unreachable under test**: the file is named
/// by an environment variable and no fixture sets one, so without this the sample index, the
/// position and the absent rows would all be untested — and each of the three is silently wrong
/// in a way the comparison would blame on the merge rather than on this.
fn rows_for(at: GenomePosition, windows_of_every_sample: &[WindowCoverage]) -> String {
    let mut rows = String::new();
    for (sample, window) in windows_of_every_sample.iter().enumerate() {
        write_the_row(
            RecordedWindow {
                sample,
                at,
                window: Some(*window),
            },
            &mut rows,
        );
    }
    rows
}

/// Write every sample's window at `at` — the first base of a locus the run has just written a
/// record for.
///
/// # Panics
///
/// A write or a flush that fails ends the process rather than the run. That is the cost of a
/// facility with no error channel into the record loop, and it is bounded by being opt-in: a run
/// that names no file cannot reach it.
pub fn record_the_windows_at(at: GenomePosition, windows_of_every_sample: &[WindowCoverage]) {
    let Some(sink) = sink() else {
        return;
    };
    let rows = rows_for(at, windows_of_every_sample);
    let mut file = sink.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    file.write_all(rows.as_bytes())
        .unwrap_or_else(|why| panic!("writing the window-coverage record: {why}"));
    // **Flushed per locus.** The run holds the writer for its whole life and never drops it in a
    // place that would flush; a buffered tail lost at exit would silently shorten the comparison.
    file.flush()
        .unwrap_or_else(|why| panic!("flushing the window-coverage record: {why}"));
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The two halves of the format, held together.** They are compiled into two different
    /// crate targets and a column reordered on one side alone still parses on the other, so this
    /// is the only thing that says they describe one file.
    #[test]
    fn a_row_survives_being_written_and_read_back() {
        let cases = [
            RecordedWindow {
                sample: 3,
                at: GenomePosition {
                    contig: ContigId(11),
                    position: Position(13_906_474),
                },
                window: Some(WindowCoverage {
                    gc_fraction: 0.371,
                    mean_depth: 14.375,
                }),
            },
            RecordedWindow {
                sample: 0,
                at: GenomePosition {
                    contig: ContigId(0),
                    position: Position(1),
                },
                // The floor silenced it: finalised, and reporting nothing.
                window: Some(WindowCoverage::absent()),
            },
            RecordedWindow {
                sample: 5,
                at: GenomePosition {
                    contig: ContigId(2),
                    position: Position(500),
                },
                // The run held none at all — a different fact from the one above.
                window: None,
            },
        ];
        for row in cases {
            let mut written = String::new();
            write_the_row(row, &mut written);
            let read = read_a_row(written.trim_end_matches('\n')).expect("its own row parses");
            assert_eq!(read.sample, row.sample);
            assert_eq!(read.at, row.at);
            // `WindowCoverage`'s own `PartialEq` is **bitwise** (`window_coverage/mod.rs`), which
            // is what the second case needs: an absent window is a pair of `NaN`s, and field-wise
            // `==` on floats would pass it however the bits came back.
            assert_eq!(
                read.window, row.window,
                "the window did not survive the round trip: {written:?}",
            );
        }
    }

    /// **One row a sample, naming that sample and that locus** — the three fields the comparison
    /// keys on, each silently wrong in its own way: an index off by one compares every sample
    /// against its neighbour's windows, a position off by one finds no centre at all, and a
    /// silenced sample omitted turns "the run had nothing here" into "the run never reached this
    /// locus", which is one of the two halves of the failure this whole facility exists to see.
    #[test]
    fn a_locus_writes_one_row_a_sample_naming_the_sample_and_the_locus() {
        let at = GenomePosition {
            contig: ContigId(1),
            position: Position(400),
        };
        let measured = WindowCoverage {
            gc_fraction: 0.5,
            mean_depth: 3.25,
        };
        // Sample 0 was measured here; sample 1 had no usable window and carries the absent pair.
        let windows_of_every_sample = [measured, WindowCoverage::absent()];

        let rows: Vec<RecordedWindow> = rows_for(at, &windows_of_every_sample)
            .lines()
            .map(|line| read_a_row(line).expect("its own rows parse"))
            .collect();

        assert_eq!(rows.len(), 2, "one row a sample, measured or not");
        assert_eq!(rows[0].sample, 0);
        assert_eq!(rows[0].at, at);
        assert_eq!(rows[0].window, Some(measured));
        assert_eq!(rows[1].sample, 1);
        assert_eq!(rows[1].at, at);
        assert_eq!(
            rows[1].window,
            Some(WindowCoverage::absent()),
            "a sample with no usable window is written down, not left out",
        );
    }

    /// A malformed row names the column, because that is what a reader looking at the file needs.
    #[test]
    fn a_row_that_will_not_parse_says_which_column_failed() {
        assert!(
            read_a_row("0\t1\t2\t3")
                .expect_err("four fields")
                .contains("five"),
            "a short row must say how many fields were expected",
        );
        assert!(
            read_a_row("0\tchr1\t2\t3\t4")
                .expect_err("a contig name where an id belongs")
                .contains("the contig"),
            "the failing column must be named",
        );
    }

    /// **Off, the recorder must not touch the windows at all**, which is what makes it free in an
    /// ordinary run and what lets it be called from inside the record loop.
    #[test]
    fn with_no_file_named_nothing_is_written_and_nothing_is_read() {
        // The variable is process-wide and this suite runs threaded, so the test asserts the
        // resolved sink rather than setting the variable: no fixture in this repository sets it,
        // so `sink()` is `None` for the whole run of the suite.
        assert!(
            std::env::var_os(PATH_VARIABLE).is_none(),
            "the test suite must not run with {PATH_VARIABLE} set — the recorder would write a \
             file per locus of every fixture",
        );
        assert!(sink().is_none());
        record_the_windows_at(
            GenomePosition {
                contig: ContigId(0),
                position: Position(1),
            },
            &[WindowCoverage::absent()],
        );
    }
}
