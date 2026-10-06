# Fit precision, step A3 — standard errors from the blocks

**Date:** 2026-09-27. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A3.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.2, §3.6 item 2. **Branch:** `fit-precision`.

## 1. What was built

The SNP/indel fit can now turn the information a pass sums (step A2) into a **standard error** for
every parameter — how far its estimate would typically move if the same kind of data were drawn
again — or say why it has none. Nothing calls it yet outside its tests; step A5 does, at the end of a
fit. No fitted number moves.

The information is kept in blocks: the cohort's eight parameters with each other, each sample's three
with each other, and each sample's three with the cohort's. Two samples are never paired, so the whole
matrix has an **arrow** shape, and an arrow is inverted exactly, block by block:

- **the cohort's errors** come from its own block less what every sample's parameters could explain
  of it, `C − Σ_s B_sᵀ A_s⁻¹ B_s`, inverted;
- **a sample's errors** come from its own block inverted plus the cohort's uncertainty carried through
  the parameters the two share, `A_s⁻¹ + A_s⁻¹ B_s V B_sᵀ A_s⁻¹`.

**Each parameter either has an error or says why it has none** (`StandardError`):

| reason | when |
|---|---|
| no information | no position's likelihood depends on it: a sample without reads, a read group without mismapped positions, the duplicated class when the run does not fit it |
| held fixed | the fit does not move it: a sample's homozygote excess in a one-sample cohort |
| not identified | other parameters can mimic its effect: the curvature left to it once the parameters before it are accounted for is below 10⁻⁸ of its own. **Only that parameter is dropped**; the others are inverted without it |
| wider than its range | its error came out wider than the whole interval the fit keeps it in, so the data do not place it anywhere in that interval. It stays in the inversion, so its uncertainty still widens the others' |

## 2. Assumptions and deviations

1. **The spec sketched a different formula; the plan asked for both to be compared and the closer to
   ship.** The sketch reads the cohort's errors from its own block alone and each sample's from
   `A_s − B_s C⁻¹ B_sᵀ` inverted. On the drawn cohorts below it is never closer to the full matrix
   than the exact arrow inverse, and further on the cohort's parameters (up to 8.6% at 20 samples,
   against 0.3%). The arrow inverse ships. Spec §3.2 still shows the sketch; amending it is the
   owner's.
2. **Parameters the data cannot tell apart are dropped one at a time — found in review, not
   anticipated by the spec.** With the duplicated class on, as runs are, some cohort parameters are
   mimicked by others at the smallest cohorts. Measured: at one sample, five of the cohort's eight
   (the density's two shapes and the duplicated class's three) are not identified; at two samples, the
   carrier Beta's second shape. The first version refused every error when the cohort's block did not
   invert — the error rates' included — and, where a pivot of rounding passed its threshold of 10⁻¹²,
   gave a share an error of 13.5 and a carrier shape one of 2 × 10⁶ (measured in review, two samples,
   300,000 positions). The threshold is now 10⁻⁸ of the parameter's own curvature before anything is
   taken from it (`1 − R²`, so units do not matter), and a failing parameter is dropped alone. **The
   errors of the others are then those with it held where it was fitted** — which of a mimicking set
   is dropped depends on their order, the cohort's first.
3. **An error wider than its parameter's whole range is not reported as an error.** Measured at two
   samples: once the carrier's second shape is dropped, the invariant share's error came out at 3.6 on
   [0, 1] and the density's shapes' at 720 and 1,148 on [0.02, 50]; the full matrix over those
   parameters does not invert at all. At 4 and 20 samples on the drawn cohorts the carrier Beta's
   shapes come out wider than their range too, and at 4 samples the duplicated share (0.149 on
   [10⁻⁹, 0.05]). These stay in the inversion — dropping them would shrink every error coupled to
   them, the optimistic direction. Not in the spec; recorded for the checkpoint.
4. **The one-sample rule is one function.** `fits_homozygote_excess(samples)` decides both whether the
   maximisation moves the excess and whether it gets an error; `fit_jointly`'s `Defaulted` marking
   uses it too. Same condition as before (two samples or more), so nothing moves.

## 3. Changes

- [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs) (new):
  `StandardError`, `StandardErrors::of(&InformationSums)`, `invert_identified` (Cholesky that drops a
  parameter whose pivot keeps less than `IDENTIFIED_SHARE` = 10⁻⁸ of its reference curvature, and
  refactors without it), `inverse_of_positive_definite`; only `+ − × ÷` and `sqrt`, fixed order.
  `expect(dead_code)` outside tests until A5.
- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): `fits_homozygote_excess`, used by
  `maximisation` and `fit_jointly`; `mod standard_errors;`.
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) (tests only): the
  comparison against the full matrix; the one- and two-sample test; a helper that rebuilds the fitted
  parameters from a `JointFit` (A4 reuses it).

## 4. Tests

Unit tests of the algebra (`standard_errors.rs`, 13), each against a dense inversion of the whole arrow
to 10⁻¹² relative: the arrow of five samples; no information at the end and in the middle of a block
(a planted slip that writes an error to its position among the kept parameters rather than its own
slot fails both middle-of-block tests); a non-finite diagonal; the one-sample excess held fixed; a
cohort parameter and a sample parameter each mimicked exactly by another, dropped alone; the threshold
keeping `1 − R²` = 10⁻⁶ and dropping 10⁻¹⁰; a parameter in units 10⁻³ as large keeping its error ×10³;
an error wider than its range; no cohort information at all; the Cholesky inverse itself.

**The smallest cohorts** (`the_smallest_cohorts_keep_their_error_rates`, 30,000 positions, 8 reads, the
duplicated class on): at one and two samples every error rate keeps its error (clean about 1.6 × 10⁻⁴
to 4.6 × 10⁻⁴, mismapped 0.006 to 0.020), and no share's error reaches 1.

**Against the full matrix** (`the_block_errors_against_the_full_matrix`): drawn cohorts at 8 reads a
position, 3,000 positions, fitted to their maximum, errors three ways — the blocks, the full matrix of
every position's scores multiplied pairwise (including the products between two samples that the
blocks leave out), and the spec's sketch. `blocks / full − 1`, over the parameters with an error:

| kind | 4 samples | 20 samples | the spec's sketch, 20 samples |
|---|---|---|---|
| shares | −0.203 to +0.023 (3) | −0.003 to 0.000 (4) | −0.086 to −0.005 |
| density shapes | +0.011 to +0.021 | −0.001 to +0.003 | −0.030 to −0.028 |
| clean error rates | −0.084 to −0.041 | −0.011 to −0.004 | −0.012 to −0.004 |
| homozygote excess | −0.022 to +0.089 | −0.065 to −0.026 | −0.067 to −0.028 |
| mismapped error rates | −0.244 to −0.043 (median −0.166) | −0.229 to −0.093 (median −0.145) | the same |

**At 20 samples and 8 reads the blocks match the full matrix, and the samples' own gaps are chance.**
Two samples' scores for the same kind are uncorrelated on average (between +0.000 and +0.001 over 190
pairs), single pairs correlate by chance — root-mean-square 0.020 for the clean rates, 0.065 for the
excesses, 0.121 for the mismapped rates, the order of the gaps — and the review measured the mismapped
rates' median gap fall from −0.145 at 3,000 positions to −0.014 at 30,000.

**At 4 samples, and at 3 reads, the gap is not chance** (measured in review). The mismapped rates'
median gap went −0.166, −0.135, −0.109 at 3,000, 30,000 and 100,000 positions; every pair's clean-rate
scores correlated at about −0.012; at 3 reads the clean rates' blocks sat 0.33 below the full matrix at
3,000 positions and 0.19 below at 30,000. Two samples share a position's unknown frequency and class,
and with few samples that shared part is a real share of each one's information, which the blocks
leave out — **the blocks' errors are then too small**. Which of the two matches the spread of the
estimates is step A4's measurement, and it must be made there: at two to four samples, at 3 reads.

Measured with `cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
(`tmp/fit_precision/module_a3i.log`): 51 passed in the fit module, 13 of them the algebra's.

## 5. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a3b.log`): counts
  in the commit message; the three pre-existing failures only; both cross-platform checksum tests pass
  unchanged.

## 6. Follow-ups, for checkpoint A

- The block approximation understates the samples' errors at small cohorts and low depth (§4). If A4
  confirms it, spec question 3 decides between the full matrix for small cohorts and a correction.
- Spec §3.2 shows the sketch, not the arrow inverse, and has neither of the two new reasons for an
  absent error (§2 items 2 and 3).
- Spec §3.6 item 2 asks for 4, 20 and 63 samples; 63 is not run (the plan's A3 names 4 and 20).
