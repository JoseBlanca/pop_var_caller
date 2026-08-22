# ng genotype prior — B1 review fixes

*Fix-application report, 2026-08-22. Branch `ng-calling-prior`. Applies
[the B1 review](../reviews/ng_calling_prior_b1_2026-08-22.md) to step B1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md).*

## Summary

**The Blocker, all four Majors and every Minor applied. Nothing disputed.** The file went from 350
lines to 432 (`git diff --numstat`: +432 / −5). The arithmetic did not change — the bit-parity test
is inside every run below — and everything else did.

Two of the findings reach past this step and are recorded rather than resolved here: which of
production's two spellings ng should be bit-identical to (§4), and that spec §9's reuse map carries
the same error the review found in this file's prose.

## 1. The Blocker: every fixture ran at a reference concentration of 1

All 25 tests drew their concentration from one helper that pinned the reference entry at 1 with at
most six alleles — the leave-one-out concentration of **one sample at a biallelic site**. What the
primitive is handed is the run's seed plus the cohort's expected allele copies, so at a thousand
diploid samples that entry is near 2,000.

The bit-parity grid now runs the reference entry at **1, 201, 2,001 and 6,001**, the alternative
total to 40, ploidy to 8 and the table to twelve alleles. Parity is the right grid to widen because
it has no tolerance to argue about — the oracle's does not survive out there (§3).

**Re-run against the fixed code, the review's three mutants all die:**

| mutant | before | after |
|---|---|---|
| every concentration entry clamped at 1.0 | passed all 25 | **FAILED, 2 tests** |
| the reference allele's entry hard-coded at 1.0 | passed all 25 | **FAILED, 2 tests** |
| the fold truncated at six alleles | passed all 25 | **FAILED, 1 test** |

## 2. The provenance was wrong about production, in a way that matters downstream

The doc said the ported primitive is "already shared between production's two callers" and that
filling a caller slice is what this port changed. **Both halves were false**, and I checked before
applying: `crate::genetics::dirichlet_multinomial_log_priors` has exactly one shipping caller, the
STR cohort's EM. The SNP/indel engine runs its own copy, `fill_log_indep_per_g_from`, which already
takes a caller's `out` and a caller's `lgamma_alpha` — the same shape this port needed, arrived at
independently and for the same stated reason.

**And the two spellings do not agree to the last bit.** The engine's sums the per-allele terms and
adds the multinomial coefficient last, where the port folds from the coefficient. Measured over 492
genotype values, **112 differ, by at most one unit in the last place**. So the parity test pins this
port to the *shared* primitive, and the GIAB 83.6% → 94.6% measurement the whole design rests on was
taken on the other one. The doc now says all of that, so whoever builds the production differential
meets it in the file rather than in a debugger.

`arch/calling_priors.md` §6's reuse row carries the same correction. **Spec §9 has the error too**
and is not edited here — a spec change is the owner's.

## 3. Three claims the measurements refuted

**The oracle's `1e-12` tolerance is a fact about the grid, not about the primitive.** Measured
disagreement between the two routes: 7.2e-14 at a reference entry of 201, 9.7e-13 at 801, 2.1e-12 at
2,001, 7.6e-12 at 6,001 — so it stops holding somewhere between four hundred and a thousand diploid
samples. Ordinary `lgamma` cancellation rather than a defect, and harmless to genotyping. The oracle
test now says why it stays at the small end and points at parity for the large one.

**Two tests did not guard the zero-count skip from opposite sides; one did.** Its two
concentrations perturbed the affected rows by about `1e-22` against values near `0.69`, six orders
of magnitude below an ulp. Replaced with the reviewer's version, which compares hom-reference at a
one-allele locus against hom-reference at a four-allele one — the skip makes them bit-identical and
adding the zero terms back does not. **Verified both ways**: it passes on the fixed code, and under
the branch removal two tests now fail where one did, one of them without reading production.

**The 2:1 test conflated two ratios.** At `α_ref = 10` the *marginalized* ratio is 19.998:1; the
22:1 the spec also records is the plug-in path's own Hardy–Weinberg ratio at the frequency its
regularisation implies. The doc now separates them, and says plainly that this test cannot guard a
seed builder's choice of `α_ref` — its fixture writes that number three functions away — which is
step D3's job. A third assertion at a reference of 1.5 shows the ratio moving to 2.997 when the pair
does, so the first two are a check on the input rather than an identity.

## 4. Everything else

`concentration_of` now shares the alternative total across the alternatives **so that they sum to
it**; before, they summed to `total × allele_count / 2` — at the six-allele shape with 0.5 that was
1.5, half again as much as the reference entry, and it corrupted three assertion messages. The inner
fold is extracted into `one_genotypes_log_prior`, so the order-of-operations promise the whole port
rests on sits beside the fold it constrains rather than forty lines above it. The monotonicity
assertion is gone, because the `1 − 1.5θ ± 3θ²` windows at the three θ do not overlap and the
tolerance check already forced the ordering — the doc says so rather than leaving a check that could
not fail. "Tracks" is out of a test name (a banned placeholder in this project's writing rules);
`oracle`, `row_for`, `checked`, `untouched`, `previous` and `term` are renamed for what they hold;
and "transcription" is gone, since to this reader it means DNA→RNA.

The debug-only value check is now documented for what it lets through in release: a concentration
entry of zero gives `lgamma(0) = +∞` and a row entry of `−∞`, a `NaN` entry gives a `NaN` row, and
neither is detected here.

**Not applied, with the reason.** The hot-path observation — `lgamma(α_a + k_a)` is recomputed per
genotype though only `alleles × ploidy` distinct values exist — was filed at Low confidence with no
code change and stays that way. It saves nothing at diploid biallelic (4 calls, 4 distinct values),
reaches about 5× at tetraploid with four alleles, and adopting it reopens the `PriorRow` scratch
contract settled at A2. **No criterion bench covers this module**, so there is no evidence either
way; that gap is recorded as a follow-up.

## 5. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | `Finished dev profile … in 6.50s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `test result: ok. 25 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `test result: ok. 22 passed; 0 failed` |
| `cargo test --lib` | 0 | `test result: ok. 4032 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 615.81s` |

Twenty-five in debug and twenty-two in release: the three `#[cfg(debug_assertions)]` value-check
tests in the sibling module are compiled out there, which is the point of running both.

## 6. Follow-ups

1. **Which spelling should ng match?** Recorded in the file and in the review's open questions. The
   recommendation is to stay with the shared primitive the plan names, because the alternative is a
   parity test against a private function in a frozen file.
2. **Spec §9's reuse map** repeats the "two callers" error.
3. **No bench covers `ng::calling::genotype_prior`**, so the one hot-path question has no evidence.
