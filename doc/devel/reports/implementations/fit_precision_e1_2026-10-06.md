# Fit precision, step E1 — each fitted number carries its standard error, where one was computed

**Date:** 2026-10-06. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step E1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §5.1. **Branch:** `fit-precision`.
Review: [fit_precision_e1_2026-10-06.md](../reviews/fit_precision_e1_2026-10-06.md), fixes:
[fixes_applied_fit_precision_e1_2026-10-06.md](../reviews/fixes_applied_fit_precision_e1_2026-10-06.md).

## 1. What was built

`Estimate<T>`, a fitted number with its warrant and evidence count, gains `standard_error: Option<f64>`: how far the
value would typically move if the same kind of data were drawn again, on the value's own scale. Its documentation's old
rule — *no uncertainty interval* — is replaced by the spec's reversal and a list of every reason the field is `None`.

Where the fits already compute an error, it is now carried:

| number | its error | `None` where |
|---|---|---|
| a sample's homozygote excess (`JointFit::hom_excess`) | the SNP/indel fit's own, from that sample's slot | the fit gives none — at one sample it holds the excess at its start |
| a sample's inbreeding coefficient (`fitted_inbreeding_of`) | the excess's: the coefficient is the excess | the excess has none |
| a library's sequencing error rate (`JointFit::sequencing_error_rates`, new) | the fit's error on that library's ordinary-position rate, found through the sample and section the library is | the fit gives none |
| a repeat-tract substitution rate (`parameters_from_the_fit`) | the binomial error √(p(1 − p)/n) over the bases compared | no base mismatched, or every one did: the formula's zero would claim the rate known exactly |

Everything else is `None`: a frequency density and a library's two rates together (several numbers, one error cannot
describe them), a sample's genotype rates (computed from other fitted numbers), supplied and defaulted values, and every
value read from a parameters file, which holds no error until step E2.

- **`JointFit::sequencing_error_rates`** replaces the conversion `parameters_from_the_fit` did itself, so the rate and
  its error leave the fit together. The errors sit in a private map built beside `noise` over the same read groups; a
  read group with rates and no error entry panics rather than reading as *fitted, with no error*.
- **No output changes.** Nothing writes or reads the field yet: calling ignores it, and the parameters file writes it
  from step E2. The library error in particular stops at `RunParameters::assemble`, whose calibration keeps only the
  multiplier and its warrant — E2's to carry.

## 2. Decided here, for the owner at checkpoint E

**A substitution rate counted with one outcome only — no mismatch, or nothing but — carries no error**, rather than an
upper bound such as 3/n. Spec §5.1 says `None` where nothing determined the value; a count that found one outcome
determines a bound, not an error of the binomial's kind. The review found the zero error on the cross-platform
fixture's own stratum (460 bases, none mismatching).

## 3. Tests

| test | what it shows |
|---|---|
| `a_samples_two_libraries_come_back_at_their_own_rates` (extended) | six samples of two libraries, drawn at ordinary-position rates 2 to 4 times apart: each sample's excess carries its own slot's error, and each library's rate from `sequencing_error_rates` carries its own error and the fit's value |
| `one_sample_fits_the_density_and_marks_what_it_could_not_fit` (extended) | at one sample the excess has no error |
| `the_fits_warrant_and_evidence_count_travel_with_the_coefficient` (extended) | an excess's error becomes the coefficient's; an excess with none gives none |
| `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own` (extended) | 3 mismatches in 4,000 bases give 0.000433 (√(p(1 − p)/n)); the fixture's own zero rate and an added stratum of 5 bases all mismatching give none |

**Mutations**, 13 across the two reviewers and 2 more after the fixes: every one that changed behaviour was killed —
the error read from the mismapped class, from the first library or the first sample for all, the read groups'
errors reversed, the excess's or coefficient's error dropped, the binomial divided by n + 1, the library error dropped
in `sequencing_error_rates`, and the one-outcome guard relaxed to admit every base mismatching. Three of the reviewers'
changed no behaviour: two because the only fixture reaching the old conversion had no library errors — the move into
the fit is what tests it now — and one because nothing reads the library error until step E2.

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Touched modules (`cargo test --release --lib -- parameter_estimation::joint::fit:: run::census_fit
  parameter_estimation::tests`): 122 passed, 1 ignored.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,025 passed, 3 failed, 9 ignored — the three failures are
  the pre-existing `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1; the cross-platform
  checksums pass unchanged.
