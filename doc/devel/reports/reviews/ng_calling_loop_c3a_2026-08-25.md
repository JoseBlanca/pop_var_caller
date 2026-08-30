# Review — ng calling loop, C3a: the quality module

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, reviewed at `3ac072f3` + the C3a working-tree diff
**Implementation report:** [C3a](../implementations/ng_calling_loop_c3a_2026-08-25.md)

## 1. Scope

The new module [`src/ng/calling/quality/mod.rs`](../../../../src/ng/calling/quality/mod.rs); four
new `CallingScratch` fields, their sizing and a bundle accessor in
[`src/ng/calling/mod.rs`](../../../../src/ng/calling/mod.rs); and a one-line visibility widening
of `log_sum_exp_2`.

**Four agents, each in its own git worktree**: reliability; a pass whose only brief was *is the
arithmetic right, checked against an independent oracle*; errors + naming; and the diff's own
quantitative and mechanistic claims.

## 2. Verdict

**Request-changes**, all applied. **1 Blocker, 5 Majors, 12 Minors.**

## 3. The Blocker — found three times, by three different oracles

`fold_samples_into_allele_counts` zeroed `next[..ploidy + live]` and wrote exactly that window.
The **next** sample, reading that buffer as `current`, tapped `ploidy` entries further up — the
counts that had just become reachable. Those slots were never written, and `prepare_for_locus`
hands every scratch buffer over filled with `f64::NAN`. The `NaN` multiplied in, survived the
rescaling (`f64::max` returns the other operand), and was converted to `−∞` by the `value > 0.0`
test on the way out.

**From count `ploidy + 1` upward the cohort's allele-count distribution was silently declared
impossible, at every cohort of two or more samples.** At 2 diploid samples the fold gave
`[−40.0, −19.31, 0.0, −∞, −∞]` where the truth is `[−40.0, −19.31, 0.0, −19.31, −40.0]`.

| oracle | disagreement |
|---|---|
| brute-force enumeration of every cohort-wide assignment | 157 of 252 cases; the 95 that agree are all one-sample |
| exact log-domain convolution, independently written | 20 of 24 shapes over ploidy 1/2/3/8 × 1–6 samples; **200 samples: 60.96 against 3316.55** |
| production's `compute_qual_via_exact_af`, same table, same constants | 1 sample agrees to `1.3e-6` Phred; **63 samples: 46.3 against 733.7** |

**The fix is one line**: `next.fill(0.0)` beside the existing `current.fill(0.0)`. Both buffers
are then zero everywhere outside each pass's write window, and stay so, because no pass writes
above its own window. With it, brute force agrees 252 of 252, the log fold to a worst `1.15e-4`
Phred, and production to a worst `1.23e-5`.

**The suite was protecting the bug.** Two tests failed the moment the fold was repaired — one
asserting a quality *below* the ceiling on a fixture whose correct answer is the ceiling, and one
whose whole property rested on the truncated axis.

**Latent in production too** (`posterior_engine.rs:3624`): hidden on a freshly grown
`RecordScratch` because `Vec::resize` zeroes, reachable on a reused one, measured at up to **0.59
Phred of drift at 8 samples**. Production's only test of that path is single-sample — the one
cohort size at which the defect cannot occur. Recorded for its owner; this branch does not touch
production.

## 4. The other Majors

### M1 — three of the diff's own measurements were measurements of the defect

Of 33 claims re-derived, 25 were correct and **all eight failures trace to the Blocker**. The
worst was not a number but an invitation: a doc paragraph concluded that the exact zero-term
override never mattered and that *"a later step trimming it should know that the tests will not
object"* — about the line that, on a working fold, is the difference between 4295.97 Phred and
the ceiling at 50 samples.

### M2 — a single `NaN` beside real probabilities was absorbed, and the message denied it

`genotype_quality`'s check was on the *winner*, and no comparison against a `NaN` is true, so a
`NaN` never wins the fold. Measured: `[0.7, NaN, 0.1]` returned genotype 0 at 5.2288 Phred — an
entirely ordinary call against a row one of whose genotypes had no probability at all. The
assertion's own message claimed *"a NaN anywhere in the row reaches here"*, true only of an
all-`NaN` row, which is the one shape the existing test covered.

**Fixed by checking the total instead**, in the same walk: `NaN` survives addition, and a row
that does not sum to one is caught by the same check.

### M3 — the panic that fires diagnoses the wrong cause

A sample whose whole likelihood row is `−∞` drives `−∞ − −∞` to a `NaN`, which surfaces at the
end as a complaint about *the normalisation* — two steps downstream of the mistake. Fixed with an
assertion in the collapse naming the sample and the cause.

### M4 — a doc comment stated a measurement that does not reproduce

The test for the one-unit-in-the-last-place clamp claimed that removing it makes the call panic.
It does not: `-10·log₁₀(0)` is `+∞` and `f32::clamp` maps `+∞` to the ceiling, so the constructor
never sees an infinity. **The clamp is untested**, and the doc now says so rather than crediting
it with a guarantee something else provides.

### M5 — a new doc block was inserted inside an existing one

`site_quality_buffers_mut`'s documentation landed between `cohort_summing_buffers_mut`'s doc
comment and its function, so the M-step bundle's prose documented the site-quality accessor and
`cohort_summing_buffers_mut` had none at all.

## 5. Three mutations that survived, and the fixtures that hid them

| mutation | why nothing saw it | now caught by |
|---|---|---|
| the collapse strides by `ploidy` instead of `allele_count` | every fixture was diploid **biallelic**, where both are 2 | `the_collapse_strides_by_the_allele_count_and_not_by_the_ploidy` — a triallelic locus, 19.472 against the mutant's 19.991 |
| `log_scale += largest` deleted | every fixture's rows peak at exactly 0.0, so the term adds nothing | `shifting_every_likelihood_by_a_constant_leaves_the_quality_alone` |
| the finite-above-ceiling cap | unreachable while the axis was truncated | `a_finite_quality_above_the_ceiling_is_capped` — 400 samples reach about 34,690 Phred |

## 6. Naming

Applied: the borrowed word *kernel* is gone from the module and from `CallingScratch`
(`copy_count_log_likelihoods`, `copy_counts_per_sample`); `count_axis_*` became
`allele_count_distribution*`; and the two entry points moved to verbs, matching the crate's
`score_one_sample` — `score_best_genotype` and `score_uncorrected_site_quality`.

Not applied, with a reason: `ArtifactTestCounts` keeps the name `spec/calling_quality.md` §10
gives it, though the reviewer is right that its first field is an `AlleleId` rather than a count.

## 7. What the reviews confirmed rather than found

- **The genotype quality is correct** — 1,200 random rows against an independent formula, worst
  disagreement `1.8e-6` Phred, which is the `f32` narrowing.
- **The prior reproduces all fifteen cells of `spec/calling_quality.md` §5.4's table** against an
  independently written closed form.
- **The rescaling invariant holds after every sample**, not just at the end: worst `1.78e-15` in
  the log domain.
- **The edges pass** with the Blocker fixed: ploidy 1 and 8, three and six alleles, one sample,
  `−∞` entries, an alternative concentration of exactly zero.

## 8. Raised for the owner, not fixed here

**`spec/calling_quality.md` §5.1's justification does not reproduce.** It argues that the
rejected `Π_s P(hom-ref)` formula grows with cohort size where the marginal stays bounded.
Measured on the corrected arithmetic, both grow in the same proportion, and in the thin-sample
regime the marginal grows faster. The arithmetic is production's and agrees with it; the question
is about the section's argument. The implementation report carries the numbers.

## 9. After the fixes

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4643 passed; 0 failed; 14 ignored` (4,620 before C3a) |
| `cargo test --release --lib ng::calling --all-features` | `597 passed; 0 failed; 3 ignored` (574 before) |

**Every mutation the review filed is now caught by the test written for it**, each re-run singly:
the Blocker (3 tests), the collapse stride (1), the running log scale (2), the exact zero term
(1). **The nine release-held checks all downgraded together under `--release`: `586 passed; 11
failed`, every check reached.**
