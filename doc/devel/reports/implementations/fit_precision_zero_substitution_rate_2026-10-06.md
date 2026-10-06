# Fit precision, after checkpoint E — a repeat-tract substitution rate is never zero

**Date:** 2026-10-06. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md), "Decided at checkpoint
E", item 1. **Branch:** `fit-precision`. Review:
[fit_precision_zero_substitution_rate_2026-10-06.md](../reviews/fit_precision_zero_substitution_rate_2026-10-06.md),
fixes: [fixes_applied_fit_precision_zero_substitution_rate_2026-10-06.md](../reviews/fixes_applied_fit_precision_zero_substitution_rate_2026-10-06.md).

## 1. What changed, and why

Calling scores the reads of a repeat tract under its stratum's substitution rate: how often a read's base disagrees
with the tract's sequence. The rate was the stratum's mismatching bases over its bases compared, so a stratum whose
reads never mismatched got a rate of exactly zero — and then a read with one mismatched base had no likelihood under
any tract length, a prior no read could move. **The owner's ruling at checkpoint E: a small value, never zero.**

`StratumSubstitutionCounts::substitution_rate` now gives, over `n` bases compared with `k` mismatching:

| count | rate |
|---|---|
| `k = 0` | `0.5 / (n + 1)` — half a mismatch |
| `k = n` | `(n + 0.5) / (n + 1)` — half a match short of one |
| otherwise | `k / n`, as before |

Half a count keeps what the count says: a zero over 500 bases becomes 0.001, one over 34 becomes 0.014. A count that
saw both outcomes keeps its own ratio. Two alternatives were put to the owner and not taken: the same half count on every
rate, which on the oracle cohort moved the 136 non-zero rates by a median of 25% and up to 50%; and the stated default
of one in a thousand, which drops what the count says. Since the rate is never zero or one, its binomial standard error
is now always written.

## 2. What moved, measured

**The four-accession oracle cohort** (`scripts/promote_ng_oracle.sh`; `tmp/fit_precision/oracle_e3/` against
`oracle_e2b/`, the run at commit `479e5ddb`):

- **the fitted file**: the 52 rates of 188 that were zero, over 34 to 506 bases compared, become 0.0010 to 0.014; no
  other line moves but the notes;
- **the calls made with the fit**: 4 of 6,706 records change, all repeat tracts: QUAL by at most 0.4, AF by at most
  0.00046. No genotype changes, and no record appears or vanishes;
- the defaults run's file moves in 8 lines of notes; every other checksum is unchanged.

**The cross-platform fixture** ([cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs)): its one
stratum compared 460 bases and found no mismatch, so its three rows go from 0.0 to 0.5 / 461 = 0.0010846, with an error
of 0.0015347; nothing else moves but 8 note lines. Its calls change at the repeat tract chrV:201 only: QUAL 14.3 to
14.9, AF 0.249885 to 0.249902, one sample's GQ 33 to 34, genotypes unchanged. Measured by putting the old rule back,
which reproduces the old calls checksum `2c0a8946…` exactly. Both checksums re-recorded, and the oracle baseline's three
moved lines.

## 3. Tests

- `a_count_with_one_outcome_only_takes_half_a_count_of_the_other` (new): 0 of 459 gives 0.5 / 460, 459 of 459 gives
  459.5 / 460, 1 of 459 gives 1 / 459.
- `every_stratum_with_bases_compared_gets_a_rate_and_its_binomial_error` (renamed from
  `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own`): the fixture's own zero count
  gives half a mismatch over its bases with the binomial error; 5 of 5 gives 5.5 / 6 with its error; 3 of 4,000
  unchanged; a stratum with nothing compared still gets no rate.
- Four mutations by the reviewer, all killed: the zero branch's `n + 1` made `n`, the error's denominator made `n + 1`,
  `== 0` made `<= 1`, the full branch made `(n − 0.5) / n`.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,041 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored.

**Open, for the owner:** a version-1 or hand-edited file can still carry a rate of zero, which calling then uses (the
review's M1).
