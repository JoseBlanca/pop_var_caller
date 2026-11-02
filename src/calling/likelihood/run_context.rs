//! **A read's misread rate after a long run of one base, in the direction it was sequenced**
//! (`doc/devel/ng/spec/str_context_error.md`). A prototype: the rates are fixed from one
//! measurement, HG002 chromosome 20 from 10 to 20 Mb at 30× (`tmp/homopolymer/` on branch
//! `flank-indel-junction`), not fitted per read group.
//!
//! A read sequenced through a long run of one base misreads the bases after it: most often the
//! next base, as one more copy of the run, and every base for at least 20 after it somewhat more
//! than elsewhere. Each rate here is the **excess** measured: reads sequenced through the run minus
//! reads reaching the same base before the run. It is added to the read group's fitted rate,
//! because the fitted rate is the background both groups of reads share and the excess is what
//! the run adds on top of it.
//!
//! **The context is read off the candidate**, with the locus's flanks around it: the question is
//! how likely this read is to have misread the candidate, and the candidate is what the read was
//! copied from.

use crate::float;
use crate::types::ErrorRate;

/// The direction a read was sequenced in, on the reference: a forward read from left to right,
/// a reverse read from right to left.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadDirection {
    Forward,
    Reverse,
}

/// The shortest run whose bases after it are charged more than the fitted rate. Below it the
/// measured excess is under 1 in 1,000.
pub const SHORTEST_RUN_CHARGED: usize = 8;

/// How many bases after a run its excess reaches. The measurement stopped at 20, flat to the end.
const BASES_CHARGED_AFTER_A_RUN: usize = 20;

/// The excess at the base right after a run, by the run's length: the chance of reading it as one
/// more copy of the run, and the chance of reading it as any other wrong base.
fn excess_at_the_next_base(run_len: usize) -> (f64, f64) {
    match run_len {
        0..=7 => (0.0, 0.0),
        8 => (0.0010, 0.0001),
        9 => (0.0004, 0.0008),
        10 => (0.0068, 0.0005),
        11 => (0.0157, 0.0018),
        12 => (0.0193, 0.0034),
        13 => (0.0316, 0.0049),
        14 => (0.0442, 0.0049),
        15 => (0.0697, 0.0103),
        16 => (0.0902, 0.0091),
        17 => (0.1031, 0.0134),
        18 => (0.1280, 0.0061),
        19 => (0.1438, 0.0086),
        20..=22 => (0.1702, 0.0085),
        _ => (0.2036, 0.0218),
    }
}

/// The excess at each of the 2nd to 20th bases after a run, any wrong base, by the run's length.
fn excess_further_on(run_len: usize) -> f64 {
    match run_len {
        0..=7 => 0.0,
        8..=11 => 0.0020,
        12..=15 => 0.0165,
        _ => 0.0610,
    }
}

/// **The longest run in `bases`** — what decides whether a candidate needs the per-direction
/// comparison at all.
#[must_use]
pub fn longest_run(bases: &[u8]) -> usize {
    let mut longest = 0;
    let mut current = 0;
    let mut previous = None;
    for &base in bases {
        current = if previous == Some(base) {
            current + 1
        } else {
            1
        };
        previous = Some(base);
        longest = longest.max(current);
    }
    longest
}

/// **How probable `observed`'s letters are against `candidate`**, both over the same bases, for a
/// read sequenced in `direction`, with `left_flank` and `right_flank` the bases around them. Bases
/// no long run precedes are scored at `rate` exactly as the flat comparison scores them.
///
/// `None` when the two differ in length, which the flat comparison answers with zero.
#[must_use]
pub fn letters_after_runs(
    observed: &[u8],
    candidate: &[u8],
    (left_flank, right_flank): (&[u8], &[u8]),
    direction: ReadDirection,
    rate: ErrorRate,
) -> Option<f64> {
    if observed.len() != candidate.len() {
        return None;
    }
    let epsilon = rate.get();
    let match_ln = float::ln(1.0 - epsilon);
    let mismatch_ln = float::ln(epsilon / 3.0);
    // The sequence in the order the read met it, and where the candidate starts in it.
    let (met, first): (Vec<u8>, usize) = match direction {
        ReadDirection::Forward => ([left_flank, candidate].concat(), left_flank.len()),
        ReadDirection::Reverse => (
            right_flank
                .iter()
                .rev()
                .chain(candidate.iter().rev())
                .copied()
                .collect(),
            right_flank.len(),
        ),
    };
    // The runs met so far, as (length, base, the step just past its end), the last few kept.
    let mut finished: Vec<(usize, u8, usize)> = Vec::new();
    let mut run_base = None;
    let mut run_len = 0;
    let mut total_ln = 0.0;
    for (step, &base) in met.iter().enumerate() {
        if run_base == Some(base) {
            run_len += 1;
        } else {
            if let Some(previous) = run_base {
                finished.push((run_len, previous, step));
            }
            run_base = Some(base);
            run_len = 1;
        }
        finished.retain(|&(_, _, end)| step - end < BASES_CHARGED_AFTER_A_RUN);
        if step < first {
            continue;
        }
        // The candidate base this step compares, and the read's.
        let at = step - first;
        let (candidate_base, read_base) = match direction {
            ReadDirection::Forward => (candidate[at], observed[at]),
            ReadDirection::Reverse => {
                let back = candidate.len() - 1 - at;
                (candidate[back], observed[back])
            }
        };
        // The longest run ending within the bases charged, the nearest of equal ones.
        let charged = finished
            .iter()
            .filter(|&&(len, _, _)| len >= SHORTEST_RUN_CHARGED)
            .max_by_key(|&&(len, _, end)| (len, end));
        total_ln += match charged {
            None => {
                if read_base == candidate_base {
                    match_ln
                } else {
                    mismatch_ln
                }
            }
            Some(&(len, run, end)) if step == end => {
                let (as_the_run, as_another) = excess_at_the_next_base(len);
                if read_base == candidate_base {
                    float::ln((1.0 - epsilon - as_the_run - as_another).max(epsilon))
                } else if read_base == run {
                    float::ln(epsilon / 3.0 + as_the_run)
                } else {
                    float::ln(epsilon / 3.0 + as_another / 2.0)
                }
            }
            Some(&(len, _, _)) => {
                let excess = excess_further_on(len);
                if read_base == candidate_base {
                    float::ln((1.0 - epsilon - excess).max(epsilon))
                } else {
                    float::ln((epsilon + excess) / 3.0)
                }
            }
        };
    }
    Some(float::exp(total_ln))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rate() -> ErrorRate {
        ErrorRate::try_new(0.002).unwrap()
    }

    /// Where no run of eight or more comes first, the comparison is the flat one, bit for bit.
    #[test]
    fn without_a_long_run_it_is_the_flat_comparison() {
        let candidate = b"ACGTACGTAAAACGT";
        let observed = b"ACGTACCTAAAACGT";
        let flat = (14.0 * float::ln(0.998) + float::ln(0.002 / 3.0)).exp();
        for direction in [ReadDirection::Forward, ReadDirection::Reverse] {
            let scored =
                letters_after_runs(observed, candidate, (b"GGC", b"TCA"), direction, rate())
                    .unwrap();
            assert!(
                (scored - flat).abs() <= 1e-12 * flat,
                "{scored} against {flat}"
            );
        }
    }

    /// A16 G A10, the `G` read as `A`: about a tenth for a forward read, which crossed the sixteen
    /// `A`s first, and much less for a reverse read, which crossed only the ten.
    #[test]
    fn the_base_after_a_long_run_is_cheap_to_misread_only_from_the_run_s_side() {
        let candidate = [&[b'A'; 16][..], b"G", &[b'A'; 10]].concat();
        let observed = [b'A'; 27];
        let forward = letters_after_runs(
            &observed,
            &candidate,
            (b"", b""),
            ReadDirection::Forward,
            rate(),
        )
        .unwrap();
        let reverse = letters_after_runs(
            &observed,
            &candidate,
            (b"", b""),
            ReadDirection::Reverse,
            rate(),
        )
        .unwrap();
        // The misread: ε/3 + 0.0902 after sixteen against ε/3 + 0.0068 after ten, about twelvefold.
        // The bases after the `G` are charged too: the forward read's ten `A`s at 0.061 each
        // (0.937¹⁰ = 0.52), the reverse read's sixteen at 0.002 (0.94) — so about 6.7.
        let ratio = forward / reverse;
        assert!(ratio > 6.0 && ratio < 7.5, "{ratio}");
        let flat = (26.0 * float::ln(0.998) + float::ln(0.002 / 3.0)).exp();
        assert!(
            reverse > flat * 5.0,
            "reverse {reverse} against flat {flat}"
        );
    }

    #[test]
    fn the_longest_run_is_counted() {
        assert_eq!(longest_run(b""), 0);
        assert_eq!(longest_run(b"ACGT"), 1);
        assert_eq!(longest_run(b"AACCCCGTTT"), 4);
    }
}
