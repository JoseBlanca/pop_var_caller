//! **How much of a site's confidence the *shape* of its variant reads does not support.**
//!
//! The site quality beside this ([`score_uncorrected_site_quality`](super::score_uncorrected_site_quality))
//! asks how unlikely it is that nobody in the cohort carries a non-reference allele, and it
//! answers from the *amount* of variant evidence. That is right at a real site and wrong at an
//! artifact — reads mapped in from a paralogous copy, or a recurrent context error, recur at a
//! steady fraction of the depth, so their number grows with coverage and **the caller gets more
//! confident about a false site the deeper it is sequenced**. Measured on GIAB HG002 by the
//! shipping caller in June 2026: the median quality of its false-positive SNPs went 1 → 3 → 150
//! from 5× to 301×, while freebayes held its own near zero at every depth
//! (`doc/devel/ng/spec/calling_quality.md` §6.1).
//!
//! **The fix is to judge the shape rather than the amount**, because shape is the thing that gets
//! *clearer* with depth at an artifact. Two tests, each a Phred subtracted from the site quality,
//! summed as independent evidence that the site is an artifact (§6.2):
//!
//! - **allele balance** — does the variant-read fraction match what the *called genotypes* imply?
//! - **strand and read position** — are the variant reads a fair sample of the site's reads, or do
//!   they pile on one strand or at one end?
//!
//! # What is here, and what is not
//!
//! Everything in this module is a function of the nine pooled counts
//! ([`ArtifactTestCounts`](super::ArtifactTestCounts)) the worker gathered while the evidence and
//! the genotypes were both in hand. **The stage that calls them is not here.** Which stage runs
//! the correction, the in-place overwrite of the called locus's one quality field, and the
//! emission threshold that reads it all belong to the output stream's own document (§3.4), which
//! is unwritten. This module supplies the arithmetic and nothing else.

use crate::genetics::lgamma;
use crate::ng::types::Phred;

/// **How much the site quality each artifact test took away.**
///
/// Recorded beside the corrected quality so the uncorrected one stays recoverable as their sum,
/// **without a second quality field for anything to read by mistake** (§3.5). That distinction
/// is not fussiness: between the correction shipping and its repair, production compared the
/// *uncorrected* number against its emission threshold while writing the *corrected* one into the
/// file, so sites were emitted `PASS` carrying a `QUAL` of 0 — 40 false positives at 30× on GIAB
/// HG002 and 64 at 50×, taken to 14 and 14 by routing both through one function.
///
/// **The sum is not always the baseline**, and the exception is the uninteresting one: where the
/// two penalties exceed the quality they are taken from, the corrected value floors at zero and
/// the arithmetic above the floor is lost. A site whose penalties exceed its baseline is one no
/// threshold would keep, so what is unrecoverable is exactly what nobody needs.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ArtifactPenalties {
    /// What the observed split between reference and alternative reads costs, against the
    /// split the called genotypes lead you to expect.
    pub allele_balance: Phred,
    /// What the alternative reads' strand and within-read placement cost, against the
    /// reference reads' own at the same site. The larger of the two, scaled by
    /// [`BIAS_RAMP_NO_POWER_BELOW`] / [`BIAS_RAMP_FULL_POWER_AT`].
    pub strand_and_read_position: Phred,
}

/// **The alternative-read count at or below which the strand and read-position test is charged
/// nothing**, with [`BIAS_RAMP_FULL_POWER_AT`] the count at which it is charged in full and a
/// linear ramp between.
///
/// **The test has no power at two or three variant reads.** Three reads land on one strand by
/// chance often enough that the test called genuine low-coverage heterozygotes biased and charged
/// them a flat 10–17 Phred — harmless where the baseline is in the hundreds, lethal at 5×. Adding
/// this ramp restored GIAB recall (5×: 0.640 → 0.706; 10×: 0.885 → 0.913) while holding the
/// medium-depth false-positive floor at 14 rather than the 28 and 43 that simply removing the test
/// gives (§6.2).
///
/// **Soft, and read off one sample's distributions.** The endpoints come from the GIAB HG002
/// alternative-read distributions in June 2026 — killed real heterozygotes carried 2–3 variant
/// reads, medium-depth artifacts 5 or more — at **one human sample**. Whether they hold on a
/// cohort at three reads a position is spec §13's Q1, open: a cohort's *pooled* variant-read count
/// crosses seven for reasons that have nothing to do with one sample's power. The ruling is to
/// port them unchanged and measure before moving anything.
///
/// **A typed constant and not production's `PVC_BIAS_RAMP` environment variable.** This repository
/// configures with types, and an environment variable that silently changes a published quality is
/// the shape §3.5 exists to prevent.
pub const BIAS_RAMP_NO_POWER_BELOW: f64 = 3.0;

/// The alternative-read count at which the strand and read-position test is charged in full —
/// see [`BIAS_RAMP_NO_POWER_BELOW`], which carries the provenance of both.
pub const BIAS_RAMP_FULL_POWER_AT: f64 = 7.0;

/// **The allele-balance test is skipped where the called genotypes expect this share of the
/// reads or more to carry the alternative.**
///
/// A homozygous-variant sample's handful of reference reads is sequencing error; a binomial
/// against a probability near one reads that as a deficit and charges for it. Inherited from
/// production, which applies the same guard at the same value (§6.2).
pub const ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE: f64 = 0.9;

/// **The ramp runs upward, and it is the compiler that says so rather than a test.** A
/// transposition would charge the full penalty at three alternative reads and nothing at seven —
/// the failure the ramp was added to prevent, in reverse — and it is a property of two literals,
/// so it can be settled where they are written instead of at run time.
const _: () = assert!(
    BIAS_RAMP_NO_POWER_BELOW >= 0.0 && BIAS_RAMP_NO_POWER_BELOW < BIAS_RAMP_FULL_POWER_AT,
    "the strand-bias ramp charges nothing below BIAS_RAMP_NO_POWER_BELOW alternative reads and \
     in full at BIAS_RAMP_FULL_POWER_AT, so the first must be the smaller and neither can be a \
     negative count of reads"
);

/// **The allele-balance guard is a share of reads**, so it lives strictly inside `(0, 1)`: at
/// zero the test would never run and at one it would never be skipped, and both are values a
/// typo reaches.
const _: () = assert!(
    ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE > 0.0 && ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE < 1.0,
    "ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE is an expected alternative-read share, so it lies \
     strictly between zero and one"
);

// ---------------------------------------------------------------------------------------
// The two-sided binomial tail — what both artifact tests charge with
// ---------------------------------------------------------------------------------------

/// **How surprising is it that only `observed` of `total` reads did the thing, when `expected`
/// of them should have?** — as a Phred, and never below zero.
///
/// This is what both artifact tests are: the allele-balance one asks it of the alternative
/// reads against the split the called genotypes imply, and the strand one asks it of the
/// alternative reads' forward-strand and placed-left counts against the reference reads' own
/// (§6.2). Everything else those two tests do is choosing which numbers to pass in.
///
/// **"Surprising" here is the two-sided Sterne tail**: the total probability of every outcome no
/// more likely than the one seen — both flanks, not just the near one. A one-sided tail would
/// call a *balanced* split at a lopsided expectation unsurprising, which is the opposite of what
/// the test is for.
///
/// Production's [`tail_phred`](../../../../src/vcf/qual_refine.rs).
pub fn two_sided_binomial_tail_phred(observed: f64, total: f64, expected_share: f64) -> f64 {
    (-10.0 * two_sided_binomial_tail(observed, total, expected_share).log10()).max(0.0)
}

/// **The probability of every outcome no more likely than `observed`**, for `observed` out of
/// `total` at a per-read probability of `expected_share`.
///
/// # One implementation, and it is the closed form
///
/// Production carries two — an exact discrete sum over all `total + 1` outcomes, and this
/// closed-form incomplete-beta tail, switching to the second only above 2,000 reads. **That
/// split is a byte-identity boundary and not an accuracy one**: the beta method is exact at
/// every `total`, and the sum is kept below the boundary purely to hold production's own output
/// bit-for-bit unchanged at the single-sample depths it had already validated
/// ([`qual_refine.rs`](../../../../src/vcf/qual_refine.rs), `EXACT_TAIL_MAX_N`).
///
/// **ng has no such obligation** — its oracle is a differential against production, not
/// byte-identity with it — and at cohort scale `total` is the pooled read count across every
/// sample, which is above 2,000 on any run that matters. So the sum would be dead code exactly
/// where the caller is used. Spec §13's Q2 left the choice to whoever wrote this; it is the
/// closed form alone, and the discrete sum lives on in this module's tests as its oracle
/// (`the_closed_form_agrees_with_the_exact_sum_across_the_grid`).
///
/// The cost is `O(log total)` probability evaluations plus two incomplete-beta evaluations,
/// independent of depth. The sum is `O(total)`.
///
/// # How the two flanks are found
///
/// The binomial is unimodal, so *no more likely than `observed`* is the near tail through
/// `observed` plus a far tail beyond the opposite-flank outcome of equal probability. Each flank
/// is monotone, so that opposite cutoff is one binary search. Production's
/// [`binom_two_sided_p_beta`](../../../../src/vcf/qual_refine.rs), ported unchanged.
fn two_sided_binomial_tail(observed: f64, total: f64, expected_share: f64) -> f64 {
    if total < 1.0 {
        return 1.0;
    }
    let total_reads = total.round() as u64;
    let total_f = total_reads as f64;
    let observed_reads = observed.round().clamp(0.0, total_f) as u64;
    let log_observed = log_binomial_probability(observed_reads as f64, total_f, expected_share);
    // The slack the "no more likely than" comparison is made with. Two outcomes of genuinely
    // equal probability can differ in their last bits, and dropping the far one from the tail
    // would make the answer depend on rounding. Production's value.
    let tolerance = 1e-7;
    // `floor((total + 1) · share)` is always an argmax, so this is the peak.
    let mode = (((total_f + 1.0) * expected_share).floor()).clamp(0.0, total_f) as u64;

    // Observed at the peak: every outcome is "no more likely", so nothing is surprising.
    if observed_reads == mode {
        return 1.0;
    }

    let two_sided = if observed_reads < mode {
        let near = binomial_at_most(observed_reads, total_reads, expected_share);
        // The far flank is `{i >= cutoff}`, the smallest `cutoff` right of the mode whose
        // probability is no greater than the observed one — empty if even the last outcome is
        // more likely than what was seen.
        let far = if log_binomial_probability(total_f, total_f, expected_share)
            > log_observed + tolerance
        {
            0.0
        } else {
            let (mut low, mut high) = (mode, total_reads);
            while low < high {
                let middle = low + (high - low) / 2;
                if log_binomial_probability(middle as f64, total_f, expected_share)
                    <= log_observed + tolerance
                {
                    high = middle;
                } else {
                    low = middle + 1;
                }
            }
            binomial_at_least(low, total_reads, expected_share)
        };
        near + far
    } else {
        let near = binomial_at_least(observed_reads, total_reads, expected_share);
        // The far flank is `{i <= cutoff}`, the largest `cutoff` left of the mode whose
        // probability is no greater than the observed one.
        let far =
            if log_binomial_probability(0.0, total_f, expected_share) > log_observed + tolerance {
                0.0
            } else {
                let (mut low, mut high) = (0_u64, mode);
                while low < high {
                    let middle = low + (high - low).div_ceil(2);
                    if log_binomial_probability(middle as f64, total_f, expected_share)
                        <= log_observed + tolerance
                    {
                        low = middle;
                    } else {
                        high = middle - 1;
                    }
                }
                binomial_at_most(low, total_reads, expected_share)
            };
        near + far
    };
    // **Floored well above zero rather than at it**, because the caller takes a base-ten
    // logarithm of this: a tail of exactly zero is a Phred of infinity, which `Phred` refuses.
    // `1e-300` is 3,000 Phred, far past any threshold, so the floor changes no decision.
    two_sided.clamp(1e-300, 1.0)
}

/// `P(X ≤ at_most)` for `X` binomial on `total` reads at `share` — exact, through the
/// regularised incomplete beta: `P(X ≤ k) = I_{1−p}(n − k, k + 1)`.
fn binomial_at_most(at_most: u64, total: u64, share: f64) -> f64 {
    if at_most >= total {
        return 1.0;
    }
    regularised_incomplete_beta((total - at_most) as f64, (at_most + 1) as f64, 1.0 - share)
}

/// `P(X ≥ at_least)` for `X` binomial on `total` reads at `share`.
fn binomial_at_least(at_least: u64, total: u64, share: f64) -> f64 {
    if at_least == 0 {
        return 1.0;
    }
    1.0 - binomial_at_most(at_least - 1, total, share)
}

/// The log of the binomial probability of exactly `count` of `total` at `share`.
///
/// **`count` outside `[0, total]` is caller error and is refused**, because the middle term
/// would ask [`lgamma`] for the logarithm of a gamma at a non-positive argument — where
/// `libm`'s answer is a finite number from the reflection formula rather than a failure, so the
/// mistake would travel as a plausible probability. Every caller in this module clamps first.
fn log_binomial_probability(count: f64, total: f64, share: f64) -> f64 {
    debug_assert!(
        (0.0..=total).contains(&count),
        "a binomial count lies between zero and the read total, got {count} of {total}"
    );
    if share <= 0.0 {
        return if count == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if share >= 1.0 {
        return if count == total {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    let log_ways = lgamma(total + 1.0) - lgamma(count + 1.0) - lgamma(total - count + 1.0);
    log_ways + count * share.ln() + (total - count) * (1.0 - share).ln()
}

/// The regularised incomplete beta `I_x(a, b)` — which is the exact binomial cumulative
/// distribution in closed form.
///
/// Numerical Recipes' `betai`, over ng's own [`lgamma`] and the continued fraction below.
/// **ng's `lgamma` is `libm`'s and production's is a hand-written Lanczos approximation**, so
/// the two are not bit-identical by construction; they agree to about one part in `1e-13`, and
/// what the agreement costs at the end of the whole correction is measured by the differential
/// this module's port is judged against (spec §11).
fn regularised_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let log_front = lgamma(a + b) - lgamma(a) - lgamma(b) + a * x.ln() + b * (1.0 - x).ln();
    let front = log_front.exp();
    // Each branch is the one the continued fraction converges quickly on; the other is reached
    // through the symmetry `I_x(a, b) = 1 − I_{1−x}(b, a)`.
    if x < (a + 1.0) / (a + b + 2.0) {
        front * incomplete_beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * incomplete_beta_continued_fraction(b, a, 1.0 - x) / b
    }
}

/// Lentz's continued fraction for the incomplete beta (Numerical Recipes' `betacf`).
///
/// It converges to about one part in `1e-15` within a few dozen iterations **whatever the
/// magnitude of `a` and `b`**, which is what makes the tail above cost the same at 30 reads and
/// at 30,000.
fn incomplete_beta_continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITERATIONS: usize = 300;
    const CONVERGED: f64 = 1e-15;
    /// Smallest magnitude the recurrence is allowed to carry, so a denominator that lands on
    /// zero is nudged off it rather than producing an infinity.
    const FLOOR: f64 = 1e-300;

    let a_plus_b = a + b;
    let a_plus_one = a + 1.0;
    let a_minus_one = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - a_plus_b * x / a_plus_one;
    if d.abs() < FLOOR {
        d = FLOOR;
    }
    d = 1.0 / d;
    let mut fraction = d;
    for iteration in 1..=MAX_ITERATIONS {
        let step = iteration as f64;
        let two_steps = 2.0 * step;

        // The even half of the recurrence.
        let term = step * (b - step) * x / ((a_minus_one + two_steps) * (a + two_steps));
        d = 1.0 + term * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + term / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        fraction *= d * c;

        // The odd half.
        let term =
            -(a + step) * (a_plus_b + step) * x / ((a + two_steps) * (a_plus_one + two_steps));
        d = 1.0 + term * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + term / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        let delta = d * c;
        fraction *= delta;
        if (delta - 1.0).abs() < CONVERGED {
            break;
        }
    }
    fraction
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The exact discrete sum, which is what the closed form above is checked against.**
    ///
    /// Every outcome from none to all, keeping the ones no more likely than the observed. This
    /// is production's shipped path below 2,000 reads
    /// ([`binom_two_sided_p`](../../../../src/vcf/qual_refine.rs)); in ng it exists only here,
    /// because it costs `O(total)` and answers nothing the closed form does not (spec §13's Q2).
    /// **It is the oracle precisely because it is the naive definition** — no continued
    /// fraction, no binary search, no unimodality argument, nothing that could be wrong in the
    /// same way as the thing it judges.
    fn exact_two_sided_binomial_tail(observed: f64, total: f64, share: f64) -> f64 {
        if total < 1.0 {
            return 1.0;
        }
        let total_reads = total.round() as u64;
        let total_f = total_reads as f64;
        let observed_reads = observed.round().clamp(0.0, total_f);
        let log_observed = log_binomial_probability(observed_reads, total_f, share);
        let tolerance = 1e-7;
        let mut accumulated = 0.0_f64;
        for outcome in 0..=total_reads {
            let log_probability = log_binomial_probability(outcome as f64, total_f, share);
            if log_probability <= log_observed + tolerance {
                accumulated += log_probability.exp();
            }
        }
        accumulated.clamp(1e-300, 1.0)
    }

    /// **The closed form against the exact sum, over every outcome of every fixture.**
    ///
    /// Read totals spanning one read to 999, expected shares from one in a hundred to
    /// ninety-nine in a hundred — both wider than what the two artifact tests will ask for,
    /// which is the point — and **every** outcome from none to all: 8,155 comparisons.
    ///
    /// **Measured: they agree to `7.0e-13` in the tail probability**, at 746 of 999 reads
    /// against an expected three in four. That is three orders of magnitude below what this
    /// asserts, and on the Phred scale the whole disagreement is `3e-12` — nothing a quality
    /// column could show.
    #[test]
    fn the_closed_form_agrees_with_the_exact_sum_across_the_grid() {
        let mut worst = 0.0_f64;
        let mut worst_case = (0_u64, 0_u64, 0.0_f64);
        for &total in &[1_u64, 2, 3, 5, 10, 37, 100, 999] {
            for &share in &[0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99] {
                for observed in 0..=total {
                    let closed = two_sided_binomial_tail(observed as f64, total as f64, share);
                    let exact = exact_two_sided_binomial_tail(observed as f64, total as f64, share);
                    let difference = (closed - exact).abs();
                    if difference > worst {
                        worst = difference;
                        worst_case = (observed, total, share);
                    }
                }
            }
        }
        assert!(
            worst < 1e-9,
            "the closed form and the exact sum disagree by {worst} at {worst_case:?}"
        );
    }

    /// **A tail a reader can check by hand.** None of ten reads carrying an allele that half of
    /// them should is two outcomes — none and all ten — each `1/1024`, so the two-sided tail is
    /// `2/1024`, and its Phred is `−10·log₁₀(0.001953125)` = 27.09.
    #[test]
    fn a_tail_a_reader_can_check_by_hand() {
        let tail = two_sided_binomial_tail(0.0, 10.0, 0.5);
        assert!(
            (tail - 2.0 / 1024.0).abs() < 1e-12,
            "the tail is {tail} against 2/1024"
        );
        let phred = two_sided_binomial_tail_phred(0.0, 10.0, 0.5);
        assert!(
            (phred - 27.092_699_609_758_306).abs() < 1e-9,
            "the Phred is {phred} against -10 log10(2/1024)"
        );
    }

    /// **The likeliest outcome is charged nothing**, because every other outcome is "no more
    /// likely" and the tail is the whole distribution. This is the branch the two artifact tests
    /// spend most of their life in — a well-behaved locus — so it must cost exactly zero rather
    /// than a small positive number that accumulates over a genome.
    #[test]
    fn the_likeliest_outcome_costs_nothing() {
        assert_eq!(two_sided_binomial_tail(50.0, 100.0, 0.5), 1.0);
        assert_eq!(two_sided_binomial_tail_phred(50.0, 100.0, 0.5), 0.0);
        // Ten reads at nine in ten: the mode is `floor(11 · 0.9)` = 9.
        assert_eq!(two_sided_binomial_tail_phred(9.0, 10.0, 0.9), 0.0);
    }

    /// **Nothing to weigh is not surprising.** A locus with no reads at all, and one with a
    /// single read, both come back charging nothing — the guard that keeps a division by zero
    /// out of the callers rather than a special case they each have to remember.
    #[test]
    fn no_reads_and_one_read_are_charged_nothing() {
        assert_eq!(two_sided_binomial_tail_phred(0.0, 0.0, 0.5), 0.0);
        assert_eq!(two_sided_binomial_tail_phred(0.0, 1.0, 0.5), 0.0);
        assert_eq!(two_sided_binomial_tail_phred(1.0, 1.0, 0.5), 0.0);
    }

    /// **A share of exactly zero or exactly one is answerable rather than a division by zero.**
    /// The two artifact tests clamp before they call, so neither can arrive here in a run — but
    /// the clamp is in the caller and this function is what would produce a `NaN` if it were
    /// ever removed.
    #[test]
    fn a_degenerate_share_still_gives_a_number() {
        assert!(two_sided_binomial_tail_phred(0.0, 20.0, 0.0).is_finite());
        assert!(two_sided_binomial_tail_phred(20.0, 20.0, 1.0).is_finite());
        // Every read where none was expected: the most surprising thing this can be told, and
        // it comes back at the floor rather than as an infinity.
        let impossible = two_sided_binomial_tail_phred(20.0, 20.0, 0.0);
        assert!(
            impossible.is_finite() && impossible > 1000.0,
            "an impossible observation is charged the floor's 3,000 Phred, not an infinity: \
             {impossible}"
        );
    }

    /// **A bigger deviation costs more, and more reads make the same deviation cost more.**
    /// Both are the properties the correction exists for: the second is what makes the test
    /// sharpen with depth where the site quality it is subtracted from also grows with depth
    /// (§6.1).
    #[test]
    fn a_wider_deviation_and_a_deeper_site_both_cost_more() {
        let mild = two_sided_binomial_tail_phred(40.0, 100.0, 0.5);
        let severe = two_sided_binomial_tail_phred(20.0, 100.0, 0.5);
        assert!(
            severe > mild,
            "twenty of a hundred is further from half than forty is, so it must cost more: \
             {severe} against {mild}"
        );

        let shallow = two_sided_binomial_tail_phred(4.0, 10.0, 0.5);
        let deep = two_sided_binomial_tail_phred(400.0, 1000.0, 0.5);
        assert!(
            deep > shallow,
            "the same four-in-ten split costs more when it is seen four hundred times in a \
             thousand: {deep} against {shallow}"
        );
    }

    /// **Two penalties, and neither can be negative or infinite**, because both are [`Phred`].
    /// The type is the check; this pins that the fields keep it rather than becoming bare
    /// `f64`s the day someone finds the conversions tiresome.
    #[test]
    fn a_penalty_pair_is_two_phreds() {
        let penalties = ArtifactPenalties {
            allele_balance: Phred::try_new(12.5).expect("a quality"),
            strand_and_read_position: Phred::try_new(0.0).expect("no penalty"),
        };
        assert_eq!(penalties.allele_balance.get(), 12.5);
        assert_eq!(penalties.strand_and_read_position.get(), 0.0);
    }
}
