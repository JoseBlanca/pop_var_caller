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

#[cfg(test)]
mod tests {
    use super::*;

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
