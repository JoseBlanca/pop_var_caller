//! **One sample's [`SampleHistogram`] as one line of text, and back** — the form the
//! hidden-duplication filter's resume file keeps (`coverage-histograms.tsv`,
//! [`crate::run::paralog_filter::resume`]) and the form a calling run writes beside its window rows
//! when asked to (`crate::run::cohort_merge::recorded_windows`).
//!
//! **Both directions in one file, beside the type they encode**, because the resume path reads
//! back what an earlier run wrote: a field reordered on one side alone would still parse on the
//! other and come back as a different histogram. The line format is fixed by the files already
//! on disk, so a change here is a change to what every existing resume point means.

use super::{CoverageByGcHistogram, SampleHistogram};

/// One sample's line, whichever of the four answers it is.
///
/// **Exhaustively matched, so a fifth silence is a compile error here** rather than a line that
/// says nothing about it.
///
/// `pub` because two modules outside this one write histograms through it: the resume file and
/// the run's recorded windows, which must agree on what a histogram's line is.
pub fn write_the_histogram(sample: usize, histogram: &SampleHistogram, into: &mut String) {
    use core::fmt::Write as _;

    let write = |into: &mut String, what: &str| {
        writeln!(into, "{sample}\t{what}").expect("writing into a String cannot fail");
    };
    match histogram {
        SampleHistogram::NoWindowFinalised => write(into, "no-window-finalised"),
        SampleHistogram::EveryWindowUnderTheFloor => write(into, "every-window-under-the-floor"),
        SampleHistogram::MedianDepthNotPositive => write(into, "median-depth-not-positive"),
        SampleHistogram::Fitted(fitted) => {
            // Destructured, so a field the histogram gains is a compile error here instead of a
            // comparison that quietly stops covering it.
            let CoverageByGcHistogram {
                window_bp,
                gc_bins,
                depth_bin_width,
                depth_bins,
                windows_folded,
                windows_under_the_floor,
                depth_summed,
                depth_rises,
                counts,
            } = fitted;
            // **The cells must be the scheme's own count**, because the line writes the scheme
            // and then whatever the vector holds: a `counts` shorter or longer than
            // `gc_bins × (depth_bins + 1)` produces a line no reader can index by bin, and two
            // such lines compare *equal* when both are wrong the same way — which is what a
            // comparison of two such files would then report as a pass. The accumulator
            // allocates exactly this and never resizes, but every field here is `pub`.
            let cells = u64::from(*gc_bins) * (u64::from(*depth_bins) + 1);
            assert_eq!(
                counts.len() as u64,
                cells,
                "sample {sample}'s histogram holds {} cells and its scheme describes {cells} \
                 ({gc_bins} GC bins over {depth_bins} depth bins and one overflow)",
                counts.len(),
            );
            // **The width by bit pattern**, because it is fitted from a sample's own median and
            // two runs that differ in its last bit cut every cell at a different depth.
            write!(
                into,
                "{sample}\tfitted\t{window_bp}\t{gc_bins}\t{}\t{depth_bins}\t\
                 {windows_folded}\t{windows_under_the_floor}\t{depth_summed}\t{depth_rises}",
                depth_bin_width.to_bits(),
            )
            .expect("writing into a String cannot fail");
            for count in counts {
                write!(into, "\t{count}").expect("writing into a String cannot fail");
            }
            into.push('\n');
        }
    }
}

/// **One line [`write_the_histogram`] wrote, read back**: the sample's index and its histogram,
/// the same value bit for bit — the bin width comes back from its bit pattern, and every cell
/// from its own column.
///
/// **What a resumed run fits the hidden-duplication filter from**
/// ([`crate::run::paralog_filter::resume`]): the histograms are finished when calling ends and
/// nothing else can rebuild them without reading every sample's reads again, so a run that stops
/// later keeps them in this form beside its spill.
///
/// # Errors
///
/// A line whose kind is not one of the four, whose fields are missing or will not parse, or whose
/// cells are not the count its scheme describes — naming what was wrong.
pub fn read_the_histogram(line: &str) -> Result<(usize, SampleHistogram), String> {
    let fields: Vec<&str> = line.split('\t').collect();
    let [sample, kind, rest @ ..] = &fields[..] else {
        return Err(format!("expected a sample and a kind, got {line:?}"));
    };
    let number = |column: &str, field: &str| -> Result<u64, String> {
        field
            .parse()
            .map_err(|_| format!("{column} is not a number: {field:?}"))
    };
    let small = |column: &str, field: &str| -> Result<u32, String> {
        u32::try_from(number(column, field)?)
            .map_err(|_| format!("{column} is too large: {field:?}"))
    };
    let sample = usize::try_from(number("the sample", sample)?)
        .map_err(|_| format!("the sample is too large: {sample:?}"))?;
    let silent = |histogram: SampleHistogram| -> Result<(usize, SampleHistogram), String> {
        if rest.is_empty() {
            Ok((sample, histogram))
        } else {
            Err(format!(
                "a {kind} line carries nothing after its kind: {line:?}"
            ))
        }
    };
    match *kind {
        "no-window-finalised" => silent(SampleHistogram::NoWindowFinalised),
        "every-window-under-the-floor" => silent(SampleHistogram::EveryWindowUnderTheFloor),
        "median-depth-not-positive" => silent(SampleHistogram::MedianDepthNotPositive),
        "fitted" => {
            let [
                window_bp,
                gc_bins,
                width_bits,
                depth_bins,
                windows_folded,
                windows_under_the_floor,
                depth_summed,
                depth_rises,
                cells @ ..,
            ] = rest
            else {
                return Err(format!("a fitted line needs eight scheme fields: {line:?}"));
            };
            let gc_bins = small("the GC bins", gc_bins)?;
            let depth_bins = small("the depth bins", depth_bins)?;
            let expected = u64::from(gc_bins) * (u64::from(depth_bins) + 1);
            if cells.len() as u64 != expected {
                return Err(format!(
                    "sample {sample}'s line holds {} cells and its scheme describes {expected}",
                    cells.len()
                ));
            }
            let counts = cells
                .iter()
                .map(|cell| small("a cell", cell))
                .collect::<Result<Vec<u32>, String>>()?;
            Ok((
                sample,
                SampleHistogram::Fitted(CoverageByGcHistogram {
                    window_bp: small("the window length", window_bp)?,
                    gc_bins,
                    depth_bin_width: f64::from_bits(number("the bin width's bits", width_bits)?),
                    depth_bins,
                    windows_folded: number("the windows folded", windows_folded)?,
                    windows_under_the_floor: number(
                        "the windows under the floor",
                        windows_under_the_floor,
                    )?,
                    depth_summed: number("the depth summed", depth_summed)?,
                    depth_rises: number("the depth rises", depth_rises)?,
                    counts,
                }),
            ))
        }
        other => Err(format!("unknown histogram kind {other:?}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **A histogram comes back from its line as the value that was written**, for each of the
    /// four kinds — the fitted one with a bin width whose last bits a decimal rendering would
    /// round, and empty cells among the full ones.
    #[test]
    fn a_histogram_survives_being_written_and_read_back() {
        let fitted = SampleHistogram::Fitted(CoverageByGcHistogram {
            window_bp: 1_000,
            gc_bins: 2,
            depth_bin_width: 0.1 + 0.2,
            depth_bins: 2,
            windows_folded: 17,
            windows_under_the_floor: 3,
            depth_summed: 412,
            depth_rises: 9,
            counts: vec![0, 4, 0, 7, 6, 0],
        });
        let cases = [
            fitted,
            SampleHistogram::NoWindowFinalised,
            SampleHistogram::EveryWindowUnderTheFloor,
            SampleHistogram::MedianDepthNotPositive,
        ];
        for (sample, histogram) in cases.into_iter().enumerate() {
            let mut written = String::new();
            write_the_histogram(sample, &histogram, &mut written);
            let read =
                read_the_histogram(written.trim_end_matches('\n')).expect("its own line parses");
            if let (SampleHistogram::Fitted(back), SampleHistogram::Fitted(sent)) =
                (&read.1, &histogram)
            {
                assert_eq!(
                    back.depth_bin_width.to_bits(),
                    sent.depth_bin_width.to_bits()
                );
            }
            assert_eq!(read, (sample, histogram), "from {written:?}");
        }
    }

    /// **A line cut short, or with a cell too many, is refused** rather than read as a histogram
    /// of another shape.
    #[test]
    fn a_histogram_line_whose_cells_do_not_match_its_scheme_is_refused() {
        assert!(read_the_histogram("0\tfitted\t1000\t2\t0\t2\t1\t0\t1\t0\t1\t2\t3").is_err());
        assert!(read_the_histogram("0\tno-window-finalised\t5").is_err());
        assert!(read_the_histogram("0\tsomething-else").is_err());
    }

    /// **A histogram's line carries everything a comparison needs to see a difference** — the
    /// fitted bin width by bit pattern, both window counts, and every cell including the empty
    /// ones. A sparse rendering would let two histograms of different widths agree wherever both
    /// were empty, which is most of the matrix.
    #[test]
    fn a_fitted_histogram_line_carries_the_bin_width_by_bits_and_every_cell() {
        let histogram = SampleHistogram::Fitted(CoverageByGcHistogram {
            window_bp: 500,
            gc_bins: 2,
            depth_bin_width: 0.25,
            depth_bins: 2,
            windows_folded: 7,
            windows_under_the_floor: 3,
            depth_summed: 900,
            depth_rises: 6,
            // Two GC rows of three cells each — two regular bins and the overflow.
            counts: vec![1, 2, 0, 0, 4, 0],
        });
        let mut written = String::new();
        write_the_histogram(4, &histogram, &mut written);

        // Where the cells begin: the sample index, the word `fitted`, and the eight scheme and
        // count fields come first.
        const CELLS_START_AT: usize = 10;

        let fields: Vec<&str> = written.trim_end().split('\t').collect();
        assert_eq!(fields[0], "4", "the sample index keys the line");
        assert_eq!(fields[1], "fitted");
        assert_eq!(
            &fields[2..CELLS_START_AT],
            ["500", "2", "4598175219545276416", "2", "7", "3", "900", "6"]
        );
        assert_eq!(
            fields[CELLS_START_AT..],
            ["1", "2", "0", "0", "4", "0"],
            "every cell, in row-major order, including the empty ones",
        );
        assert_eq!(
            fields[4].parse::<u64>().map(f64::from_bits),
            Ok(0.25),
            "the bin width is written as the bits it is, so two runs a last bit apart differ here",
        );
    }

    /// **Each silence is its own word**, because a consumer that acts on the difference is the
    /// whole reason `SampleHistogram` is an enum rather than an `Option`.
    #[test]
    fn a_sample_with_no_histogram_says_which_silence_it_was() {
        for (histogram, word) in [
            (SampleHistogram::NoWindowFinalised, "no-window-finalised"),
            (
                SampleHistogram::EveryWindowUnderTheFloor,
                "every-window-under-the-floor",
            ),
            (
                SampleHistogram::MedianDepthNotPositive,
                "median-depth-not-positive",
            ),
        ] {
            let mut written = String::new();
            write_the_histogram(0, &histogram, &mut written);
            assert_eq!(written, format!("0\t{word}\n"));
        }
    }

    /// **A histogram whose cells do not fit its own scheme is refused rather than written**, for
    /// the reason the assertion gives: a line no reader can index by bin compares equal to another
    /// wrong the same way, which a comparison of two such files would report as a pass.
    #[test]
    #[should_panic(expected = "its scheme describes")]
    fn a_histogram_whose_cells_do_not_fit_its_scheme_is_refused() {
        let histogram = SampleHistogram::Fitted(CoverageByGcHistogram {
            window_bp: 500,
            gc_bins: 2,
            depth_bin_width: 0.25,
            depth_bins: 2,
            windows_folded: 0,
            windows_under_the_floor: 0,
            depth_summed: 0,
            depth_rises: 0,
            // Two GC rows of three cells each is six; this is five.
            counts: vec![0; 5],
        });
        write_the_histogram(0, &histogram, &mut String::new());
    }
}
