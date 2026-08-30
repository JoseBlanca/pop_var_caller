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
use crate::ng::calling::quality::ArtifactTestCounts;
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

// ---------------------------------------------------------------------------------------
// The first test: does the split match what the called genotypes imply?
// ---------------------------------------------------------------------------------------

/// **How much of the site quality the alternative reads' *share* does not support.**
///
/// A single heterozygote should show about half its reads carrying the alternative; a cohort
/// where two samples in sixty carry one copy each should show about that fraction overall. The
/// penalty is the two-sided binomial improbability of the observed split against that
/// expectation.
///
/// **The expectation comes from the called genotypes and never from the fitted allele
/// frequency.** The frequency adapts to the artifact and would excuse it; the genotypes do not
/// (§6.2). It arrives already summed, as
/// [`genotype_expected_alternative_reads`](super::ArtifactTestCounts::genotype_expected_alternative_reads).
///
/// # Two guards, both inherited and both load-bearing
///
/// **Only a deficit is charged.** These artifacts present *fewer* alternative reads than a real
/// call at that frequency would. An excess is a different phenomenon and this test says nothing
/// about it.
///
/// **The test is skipped where the genotypes expect nearly all the reads to be alternative**
/// ([`ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE`]). A homozygous-variant sample's handful of reference
/// reads is sequencing error; a binomial against a probability near one reads that as a deficit
/// and charges for it.
///
/// **No ramp, and that is a measurement rather than a symmetry.** Production's per-record
/// decomposition found this test at zero for true heterozygotes at every depth while growing with
/// depth for false positives; it is the strand test next door that had the power problem (§6.2).
///
/// Production's, inline in [`refine_qual`](../../../../src/vcf/qual_refine.rs).
pub fn allele_balance_penalty(counts: &ArtifactTestCounts) -> Phred {
    // A locus with no alternative reads has nothing to weigh. It cannot arise from a run — the
    // worker hands back no summary at all in that case — but this is where a division by zero
    // would be, and two comparisons are cheaper than the argument that it cannot happen.
    if counts.total_reads < 1.0 || counts.alternative_reads < 1.0 {
        return NO_PENALTY;
    }
    let expected_share = (counts.genotype_expected_alternative_reads / counts.total_reads)
        .clamp(EXPECTED_SHARE_FLOOR, EXPECTED_SHARE_CEILING);
    if expected_share >= ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE {
        return NO_PENALTY;
    }
    if counts.alternative_reads >= expected_share * counts.total_reads {
        return NO_PENALTY;
    }
    as_penalty(two_sided_binomial_tail_phred(
        counts.alternative_reads,
        counts.total_reads,
        expected_share,
    ))
}

/// **The expected alternative share is held strictly inside `(0, 1)`**, because the binomial
/// tail takes a logarithm at both ends.
///
/// **Neither end can bind while the two guards above stand**, and saying so is worth more than
/// the clamp: the deficit rule reaches the tail only where `alternative_reads` is below
/// `expected_share × total_reads` with at least one alternative read, which puts the share above
/// `1 / total_reads` — never at the floor; and [`ALLELE_BALANCE_SKIPPED_AT_OR_ABOVE`] returns
/// before any share reaches [`EXPECTED_SHARE_CEILING`]. Both are production's and both are kept,
/// as the net under two constants somebody may move.
const EXPECTED_SHARE_FLOOR: f64 = 1e-6;

/// The other end of [`EXPECTED_SHARE_FLOOR`]'s clamp, which carries the reasoning for both.
const EXPECTED_SHARE_CEILING: f64 = 0.999;

/// Zero Phred — what both tests return where they have nothing to weigh.
const NO_PENALTY: Phred = Phred::ZERO;

/// A finite non-negative Phred from the tail, which is what
/// [`two_sided_binomial_tail_phred`] always returns: it floors its probability at `1e-300`, so
/// the largest penalty expressible is about 3,000, and it takes `max(0.0)` at the end.
fn as_penalty(phred: f64) -> Phred {
    Phred::try_new(phred as f32).expect(
        "a two-sided binomial tail is a finite non-negative Phred: its probability is floored \
         at 1e-300 before the logarithm and the result takes max(0.0)",
    )
}

// ---------------------------------------------------------------------------------------
// The second test: are the alternative reads a fair sample of the site's reads?
// ---------------------------------------------------------------------------------------

/// **How much of the site quality the alternative reads' *provenance* does not support.**
///
/// A real variant's reads come off both strands and sit at varied places within the read; an
/// artifact's often pile on one strand or at one end. This charges the improbability of the
/// alternative reads' forward-strand share, and of their placed-left share — **the larger of the
/// two**, since either alone is evidence the site is an artifact.
///
/// **The expectation is the reference reads' own share at the same site, not one half.** A locus
/// whose coverage is one-sided for innocent reasons — the edge of a capture target, a repeat that
/// only one orientation maps into — has reference reads that are one-sided too, and comparing the
/// alternative reads against them rather than against a fixed half is what keeps the test from
/// charging the whole locus for its neighbourhood (§6.2).
///
/// **And it is ramped in, because at two or three alternative reads it has no power.** See
/// [`BIAS_RAMP_NO_POWER_BELOW`], which carries what removing the ramp cost and what restoring it
/// bought.
///
/// Production's, inline in [`refine_qual`](../../../../src/vcf/qual_refine.rs).
pub fn strand_and_read_position_penalty(counts: &ArtifactTestCounts) -> Phred {
    if counts.alternative_reads < 1.0 {
        return NO_PENALTY;
    }
    let expected_forward_share =
        reference_share(counts.reference_forward_reads, counts.reference_reads);
    let expected_placed_left_share =
        reference_share(counts.reference_placed_left_reads, counts.reference_reads);
    let unramped = two_sided_binomial_tail_phred(
        counts.alternative_forward_reads,
        counts.alternative_reads,
        expected_forward_share,
    )
    .max(two_sided_binomial_tail_phred(
        counts.alternative_placed_left_reads,
        counts.alternative_reads,
        expected_placed_left_share,
    ));
    as_penalty(unramped * bias_test_power(counts.alternative_reads))
}

/// **What share of the reference reads did the thing** — the expectation the alternative reads
/// are weighed against.
///
/// `0.5` where there are no reference reads at all: with nothing to compare against, an even
/// split is the assumption that charges least. Clamped away from both ends because the binomial
/// tail takes a logarithm there, and because a site where *every* reference read is on one strand
/// would otherwise make a single opposite-strand alternative read infinitely surprising.
fn reference_share(reference_reads_that_did_it: f64, reference_reads: f64) -> f64 {
    if reference_reads > 0.0 {
        (reference_reads_that_did_it / reference_reads)
            .clamp(REFERENCE_SHARE_FLOOR, REFERENCE_SHARE_CEILING)
    } else {
        0.5
    }
}

/// **How much of the raw strand and read-position penalty is charged**, from nothing at
/// [`BIAS_RAMP_NO_POWER_BELOW`] alternative reads or fewer to all of it at
/// [`BIAS_RAMP_FULL_POWER_AT`] or more, linear between.
///
/// This is not a correction for multiple testing or a confidence adjustment; it is a statement
/// that **the test cannot tell a real heterozygote's chance pile-up from an artifact's** until
/// there are enough alternative reads for the pile-up to mean something. Production's
/// `bias_power_factor`.
fn bias_test_power(alternative_reads: f64) -> f64 {
    if alternative_reads <= BIAS_RAMP_NO_POWER_BELOW {
        0.0
    } else if alternative_reads >= BIAS_RAMP_FULL_POWER_AT {
        1.0
    } else {
        (alternative_reads - BIAS_RAMP_NO_POWER_BELOW)
            / (BIAS_RAMP_FULL_POWER_AT - BIAS_RAMP_NO_POWER_BELOW)
    }
}

/// **The reference reads' share is held inside `[0.01, 0.99]`**, which is a wider clamp than the
/// allele-balance test's and, unlike that one, **binds on real data**: a locus every one of whose
/// reference reads is on the forward strand is an ordinary thing to meet, and without the clamp a
/// single reverse-strand alternative read there is charged the tail's whole 3,000-Phred floor.
/// Production's, at the same value.
const REFERENCE_SHARE_FLOOR: f64 = 0.01;

/// The other end of [`REFERENCE_SHARE_FLOOR`]'s clamp.
const REFERENCE_SHARE_CEILING: f64 = 0.99;

// ---------------------------------------------------------------------------------------
// The subtraction
// ---------------------------------------------------------------------------------------

/// **The site quality a file carries** — the worker's baseline, less what the shape of the
/// variant reads does not support.
///
/// The two penalties are summed as **independent** evidence that the site is an artifact: one
/// asks whether there are as many alternative reads as the calls imply, the other whether the
/// ones there are came from everywhere they should have, and a site can fail either without
/// telling you anything about the other.
///
/// **Floored at zero**, because a Phred is not negative — and a site whose penalties exceed its
/// baseline is one no threshold would keep, so the arithmetic lost below the floor is exactly the
/// arithmetic nobody needs.
///
/// # The one rule about the caller, and it is why the penalties come back
///
/// **A called locus carries exactly one quality at every moment** (§3.5). The worker writes the
/// baseline into it and the output stage overwrites that with what this returns; the penalties
/// travel *beside* the corrected quality rather than the baseline travelling beside it, so there
/// is never a second quality field for a threshold to read by mistake. Production kept no
/// corrected value at all and recomputed it at write time, and for sixteen days its emission gate
/// compared the baseline while the corrected number went into the file — 40 sites emitted `PASS`
/// with a written `QUAL` of 0 at 30× on GIAB HG002, and 64 at 50×.
///
/// **A locus with no summary never reaches here.** The worker hands back `None` where the
/// candidate table is the reference alone, or where no read reached an alternative — a quarter of
/// built loci on both benchmarks — and such a locus keeps its baseline unchanged. That branch is
/// the output stage's, which is why this takes the summary by value rather than an `Option`.
pub fn correct_site_quality(
    baseline: Phred,
    counts: &ArtifactTestCounts,
) -> (Phred, ArtifactPenalties) {
    let penalties = ArtifactPenalties {
        allele_balance: allele_balance_penalty(counts),
        strand_and_read_position: strand_and_read_position_penalty(counts),
    };
    let charged = f64::from(penalties.allele_balance.get())
        + f64::from(penalties.strand_and_read_position.get());
    let corrected = (f64::from(baseline.get()) - charged).max(0.0);
    (
        Phred::try_new(corrected as f32)
            .expect("a difference of two finite non-negative Phreds, floored at zero"),
        penalties,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::AlleleId;

    /// A locus's pooled counts, with the two the allele-balance test reads spelled out and the
    /// strand ones set to a fair split so they charge nothing.
    fn balanced_counts(
        reference_reads: f64,
        alternative_reads: f64,
        genotype_expected_alternative_reads: f64,
    ) -> ArtifactTestCounts {
        ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads,
            reference_forward_reads: reference_reads / 2.0,
            reference_placed_left_reads: reference_reads / 2.0,
            alternative_reads,
            alternative_forward_reads: alternative_reads / 2.0,
            alternative_placed_left_reads: alternative_reads / 2.0,
            total_reads: reference_reads + alternative_reads,
            genotype_expected_alternative_reads,
        }
    }

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

    // -----------------------------------------------------------------------------------
    // The allele-balance test
    // -----------------------------------------------------------------------------------

    /// **A heterozygote showing half its reads costs nothing.** Twenty reference reads and
    /// twenty alternative, at genotypes expecting twenty — the ordinary case, and the one the
    /// whole genome spends its life in.
    #[test]
    fn a_split_the_genotypes_predict_costs_nothing() {
        let counts = balanced_counts(20.0, 20.0, 20.0);
        assert_eq!(allele_balance_penalty(&counts).get(), 0.0);
    }

    /// **A deficit costs, and costs about ten times more when the depth is ten times greater.**
    /// One read in five carries the alternative where the called genotypes say half should —
    /// the shape of an artifact that recurs at a steady fraction of the depth. At 50 reads that
    /// is charged **46.2** Phred and at 500 it is **430.8**.
    ///
    /// **That ratio is the whole point of the correction** (§6.1). The site quality this is
    /// subtracted from also grows about linearly with the variant-read count, so a penalty that
    /// did *not* grow with depth would be swamped at 500 reads and the caller would go on
    /// getting more confident about a false site the deeper it was sequenced.
    #[test]
    fn a_deficit_costs_and_costs_more_at_depth() {
        let shallow = f64::from(allele_balance_penalty(&balanced_counts(40.0, 10.0, 25.0)).get());
        let deep = f64::from(allele_balance_penalty(&balanced_counts(400.0, 100.0, 250.0)).get());
        assert!(
            (shallow - 46.22).abs() < 0.01,
            "fifty reads at one in five against a half: {shallow}"
        );
        assert!(
            (deep - 430.81).abs() < 0.01,
            "five hundred reads at one in five against a half: {deep}"
        );
        assert!(
            deep / shallow > 8.0,
            "ten times the depth at the same split has to cost close to ten times as much, or \
             the penalty is swamped by a site quality that does grow with depth: {deep} against \
             {shallow}"
        );
    }

    /// **An excess is charged nothing, and that is a rule rather than an oversight.** These
    /// artifacts present *fewer* alternative reads than a real call at that frequency would; an
    /// excess is a different phenomenon this test says nothing about (§6.2). Thirty-five of
    /// fifty reads where the genotypes expect twenty-five is as far from the expectation as the
    /// charged fixture above, in the other direction, and it costs zero.
    #[test]
    fn an_excess_of_alternative_reads_is_charged_nothing() {
        let counts = balanced_counts(15.0, 35.0, 25.0);
        assert_eq!(allele_balance_penalty(&counts).get(), 0.0);
        // The mirror image is charged, so the fixture is not simply too mild to register.
        assert!(allele_balance_penalty(&balanced_counts(35.0, 15.0, 25.0)).get() > 0.0);
    }

    /// **A cohort of homozygous-variant samples is skipped rather than charged.** Its handful of
    /// reference reads is sequencing error, and a binomial against a probability near one reads
    /// that as a deficit. Ninety-six of a hundred reads carry the alternative where the
    /// genotypes expect all hundred: without the guard that is charged 12.7 Phred, with it
    /// nothing.
    #[test]
    fn a_cohort_the_genotypes_call_homozygous_variant_is_skipped() {
        let counts = balanced_counts(4.0, 96.0, 100.0);
        assert_eq!(allele_balance_penalty(&counts).get(), 0.0);
        // What the guard is worth: the same split weighed at the highest share the guard still
        // lets through.
        let just_under_the_guard = two_sided_binomial_tail_phred(96.0, 100.0, 0.89);
        assert!(
            just_under_the_guard > 10.0,
            "the skipped fixture would be charged {just_under_the_guard} Phred just below the \
             guard, so the guard is doing work"
        );
    }

    /// **Reads nobody's genotype expects are charged nothing, and that is the deficit rule
    /// rather than a gap.**
    ///
    /// Ten alternative reads where every called genotype is homozygous reference looks like the
    /// worst artifact there is, and this test says nothing about it — because *more* alternative
    /// reads than expected is an excess, and only a deficit is charged (§6.2). What catches that
    /// site is the strand test beside this one, or nothing.
    ///
    /// **It also means neither end of the expected-share clamp can bind.** The deficit branch
    /// needs `alternative_reads < expected_share × total_reads` with at least one alternative
    /// read, so the share is above `1 / total_reads` wherever the tail is reached — never at the
    /// floor — and the guard above returns before any share reaches the ceiling. Both are kept
    /// as production has them, as the net under two constants somebody may move.
    #[test]
    fn reads_nobodys_genotype_expects_are_an_excess_and_charged_nothing() {
        let counts = balanced_counts(90.0, 10.0, 0.0);
        assert_eq!(allele_balance_penalty(&counts).get(), 0.0);
    }

    /// **A locus with nothing to weigh is charged nothing.** No summary like this reaches the
    /// test from a run — the worker hands back `None` where no read reached an alternative —
    /// but the guard is where a division by zero would be.
    #[test]
    fn a_locus_with_no_reads_is_charged_nothing() {
        assert_eq!(
            allele_balance_penalty(&balanced_counts(0.0, 0.0, 0.0)).get(),
            0.0
        );
        assert_eq!(
            allele_balance_penalty(&balanced_counts(10.0, 0.0, 5.0)).get(),
            0.0
        );
    }

    // -----------------------------------------------------------------------------------
    // The strand and read-position test
    // -----------------------------------------------------------------------------------

    /// A locus whose reference reads split evenly on both axes, with the alternative reads'
    /// two counts given, so a fixture says only what it is about.
    fn strand_counts(
        reference_reads: f64,
        alternative_reads: f64,
        alternative_forward_reads: f64,
        alternative_placed_left_reads: f64,
    ) -> ArtifactTestCounts {
        ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads,
            reference_forward_reads: reference_reads / 2.0,
            reference_placed_left_reads: reference_reads / 2.0,
            alternative_reads,
            alternative_forward_reads,
            alternative_placed_left_reads,
            total_reads: reference_reads + alternative_reads,
            // Whatever the called genotypes expect is the other test's business; a value that
            // matches the split keeps this fixture from being about two things at once.
            genotype_expected_alternative_reads: alternative_reads,
        }
    }

    /// **Alternative reads drawn as evenly as the reference reads cost nothing.** Twenty of
    /// forty on each axis, against reference reads split down the middle — the ordinary case.
    #[test]
    fn alternative_reads_as_evenly_drawn_as_the_reference_cost_nothing() {
        let counts = strand_counts(60.0, 40.0, 20.0, 20.0);
        assert_eq!(strand_and_read_position_penalty(&counts).get(), 0.0);
    }

    /// **Forty alternative reads all on one strand are charged 111.6 Phred**, where the
    /// reference reads at the same site split evenly. This is the artifact shape the test is
    /// named for.
    #[test]
    fn a_one_strand_pile_up_at_depth_is_charged() {
        let counts = strand_counts(60.0, 40.0, 40.0, 20.0);
        let penalty = f64::from(strand_and_read_position_penalty(&counts).get());
        assert!(
            (penalty - 117.40).abs() < 0.01,
            "forty of forty on the forward strand: {penalty}"
        );
    }

    /// **The larger of the two axes is what is charged.** The same forty reads, this time split
    /// evenly by strand but every one placed left, cost the same 117.4 — so a site whose
    /// artifact shows in only one of the two is caught by that one.
    #[test]
    fn either_axis_alone_is_enough_to_charge_a_site() {
        let by_strand = strand_and_read_position_penalty(&strand_counts(60.0, 40.0, 40.0, 20.0));
        let by_position = strand_and_read_position_penalty(&strand_counts(60.0, 40.0, 20.0, 40.0));
        assert_eq!(by_strand, by_position);
        assert!(by_strand.get() > 100.0);
    }

    /// **Three alternative reads on one strand are charged exactly nothing, and five are
    /// charged half.**
    ///
    /// This is the ramp, and it is the reason the test can be kept at all. Without it, three
    /// reads landing on one strand by chance — which happens one time in four at an even split —
    /// took genuine low-coverage heterozygotes for artifacts and charged them 10 to 17 Phred,
    /// which is harmless where the baseline is in the hundreds and fatal at 5×.
    ///
    /// **Measured: the unramped penalty at five reads is 12.04 Phred and what is charged is
    /// 6.02**, exactly half, because five sits halfway between the ramp's three and seven.
    #[test]
    fn the_ramp_charges_nothing_at_three_reads_and_half_at_five() {
        let three = strand_and_read_position_penalty(&strand_counts(60.0, 3.0, 3.0, 1.5));
        assert_eq!(
            three.get(),
            0.0,
            "three alternative reads have no power, so the test says nothing about them"
        );

        let five =
            f64::from(strand_and_read_position_penalty(&strand_counts(60.0, 5.0, 5.0, 2.5)).get());
        let unramped = two_sided_binomial_tail_phred(5.0, 5.0, 0.5);
        assert!(
            (unramped - 12.041).abs() < 0.01,
            "five of five on one strand against an even expectation: {unramped}"
        );
        assert!(
            (five - unramped / 2.0).abs() < 0.01,
            "five reads sit halfway along the ramp from three to seven, so half the penalty is \
             charged: {five} against {unramped}"
        );

        let seven =
            f64::from(strand_and_read_position_penalty(&strand_counts(60.0, 7.0, 7.0, 3.5)).get());
        assert!(
            (seven - two_sided_binomial_tail_phred(7.0, 7.0, 0.5)).abs() < 0.01,
            "seven reads are the top of the ramp and pay in full: {seven}"
        );
    }

    /// **A site whose reference reads are one-sided does not charge its alternative reads for
    /// it.** Every reference read on the forward strand and every alternative read too: the
    /// expectation is read from the reference, so the alternative reads look exactly as expected
    /// and cost nothing on that axis.
    ///
    /// Comparing against a fixed one half instead would charge this site 117.4 Phred — the
    /// number a caller would take off every locus at the edge of a capture target.
    #[test]
    fn a_one_sided_site_charges_nothing_when_both_alleles_lean_the_same_way() {
        let counts = ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads: 60.0,
            reference_forward_reads: 60.0,
            reference_placed_left_reads: 30.0,
            alternative_reads: 40.0,
            alternative_forward_reads: 40.0,
            alternative_placed_left_reads: 20.0,
            total_reads: 100.0,
            genotype_expected_alternative_reads: 40.0,
        };
        assert_eq!(strand_and_read_position_penalty(&counts).get(), 0.0);

        let against_a_fixed_half = two_sided_binomial_tail_phred(40.0, 40.0, 0.5);
        assert!(
            (against_a_fixed_half - 117.40).abs() < 0.01,
            "what a fixed expectation of one half would have charged: {against_a_fixed_half}"
        );
    }

    /// **A locus with no reference reads falls back to an even expectation.** Nothing to compare
    /// against, so the assumption that charges least is the one taken — and the clamp keeps the
    /// fall-back from ever being zero or one.
    #[test]
    fn no_reference_reads_falls_back_to_an_even_expectation() {
        assert_eq!(reference_share(0.0, 0.0), 0.5);
        let counts = strand_counts(0.0, 40.0, 20.0, 20.0);
        assert_eq!(strand_and_read_position_penalty(&counts).get(), 0.0);
    }

    /// **The reference share is clamped, and unlike the other test's clamp this one binds on
    /// ordinary data.** Every reference read on one strand is a real thing to meet; without the
    /// clamp a single opposite-strand alternative read there is weighed against a probability of
    /// zero and charged the tail's floor.
    #[test]
    fn a_wholly_one_sided_reference_is_clamped_off_the_endpoint() {
        assert_eq!(reference_share(60.0, 60.0), REFERENCE_SHARE_CEILING);
        assert_eq!(reference_share(0.0, 60.0), REFERENCE_SHARE_FLOOR);

        // What the clamp is worth: one alternative read of forty on the other strand, at a site
        // whose every reference read leans one way.
        let counts = ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads: 60.0,
            reference_forward_reads: 60.0,
            reference_placed_left_reads: 30.0,
            alternative_reads: 40.0,
            alternative_forward_reads: 39.0,
            alternative_placed_left_reads: 20.0,
            total_reads: 100.0,
            genotype_expected_alternative_reads: 40.0,
        };
        let penalty = f64::from(strand_and_read_position_penalty(&counts).get());
        assert!(
            penalty.is_finite() && penalty < 10.0,
            "one read of forty against the site's own lean is a small charge, not the tail's \
             floor: {penalty}"
        );
    }

    /// **A locus with no alternative reads is charged nothing.** As with the other test, no such
    /// summary reaches here from a run; this is where the division by zero would be.
    #[test]
    fn a_locus_with_no_alternative_reads_is_charged_nothing_by_the_strand_test() {
        assert_eq!(
            strand_and_read_position_penalty(&strand_counts(60.0, 0.0, 0.0, 0.0)).get(),
            0.0
        );
    }

    // -----------------------------------------------------------------------------------
    // The subtraction
    // -----------------------------------------------------------------------------------

    /// **A clean locus keeps its quality exactly.** Forty alternative reads where the genotypes
    /// expect forty, drawn as evenly as the reference reads: both tests charge nothing and the
    /// baseline comes back unchanged. This is the ordinary case, and *unchanged* has to mean the
    /// same bits rather than nearly — a correction that shaved a fraction of a Phred off every
    /// clean site would move a genome's worth of calls across a threshold.
    #[test]
    fn a_clean_locus_keeps_its_quality_exactly() {
        let baseline = Phred::try_new(742.5).expect("a quality");
        let (corrected, penalties) =
            correct_site_quality(baseline, &strand_counts(60.0, 40.0, 20.0, 20.0));
        assert_eq!(corrected, baseline);
        assert_eq!(penalties.allele_balance.get(), 0.0);
        assert_eq!(penalties.strand_and_read_position.get(), 0.0);
    }

    /// **The two penalties are summed, and the baseline is recoverable from what comes back.**
    /// A locus failing both tests — one read in five carrying the alternative where the
    /// genotypes say half should, *and* every one of them on the same strand — is charged both,
    /// and `corrected + balance + strand` returns the baseline.
    ///
    /// That recoverability is the reason nothing needs a second quality field to hold the
    /// baseline in (§3.5).
    #[test]
    fn a_locus_failing_both_tests_is_charged_both_and_the_baseline_is_recoverable() {
        let counts = ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads: 400.0,
            reference_forward_reads: 200.0,
            reference_placed_left_reads: 200.0,
            alternative_reads: 100.0,
            alternative_forward_reads: 100.0,
            alternative_placed_left_reads: 50.0,
            total_reads: 500.0,
            genotype_expected_alternative_reads: 250.0,
        };
        let baseline = Phred::try_new(900.0).expect("a quality");
        let (corrected, penalties) = correct_site_quality(baseline, &counts);
        assert!(penalties.allele_balance.get() > 0.0);
        assert!(penalties.strand_and_read_position.get() > 0.0);

        let recovered = f64::from(corrected.get())
            + f64::from(penalties.allele_balance.get())
            + f64::from(penalties.strand_and_read_position.get());
        assert!(
            (recovered - 900.0).abs() < 0.01,
            "the baseline is the corrected quality plus the two penalties: {recovered}"
        );
    }

    /// **Penalties larger than the baseline floor at zero rather than going negative.** A
    /// low-quality site that also fails both tests is the case, and a [`Phred`] has no negative
    /// value to give it — so the subtraction is floored where it is done rather than where it is
    /// read.
    ///
    /// The baseline is *not* recoverable here, and that is the deliberate exception: a site whose
    /// penalties exceed its baseline is one no threshold would keep.
    #[test]
    fn penalties_larger_than_the_baseline_floor_at_zero() {
        let counts = ArtifactTestCounts {
            primary_alternative: AlleleId(1),
            reference_reads: 400.0,
            reference_forward_reads: 200.0,
            reference_placed_left_reads: 200.0,
            alternative_reads: 100.0,
            alternative_forward_reads: 100.0,
            alternative_placed_left_reads: 50.0,
            total_reads: 500.0,
            genotype_expected_alternative_reads: 250.0,
        };
        let (corrected, penalties) =
            correct_site_quality(Phred::try_new(30.0).expect("a weak quality"), &counts);
        assert_eq!(corrected.get(), 0.0);
        assert!(
            f64::from(penalties.allele_balance.get())
                + f64::from(penalties.strand_and_read_position.get())
                > 30.0,
            "the fixture has to charge more than the baseline for this to be about the floor"
        );
    }

    /// **What comes back is what the two tests give, not a second computation of them.** A
    /// correction that re-derived either penalty could drift from the function that names it;
    /// this pins that the pair is assembled from the two published functions.
    #[test]
    fn the_pair_that_comes_back_is_what_the_two_tests_charge() {
        let counts = strand_counts(60.0, 40.0, 38.0, 20.0);
        let (_, penalties) =
            correct_site_quality(Phred::try_new(500.0).expect("a quality"), &counts);
        assert_eq!(penalties.allele_balance, allele_balance_penalty(&counts));
        assert_eq!(
            penalties.strand_and_read_position,
            strand_and_read_position_penalty(&counts)
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
