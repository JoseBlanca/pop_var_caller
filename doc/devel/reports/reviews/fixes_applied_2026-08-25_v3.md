# Fix Application Report: ng_calling_loop_b1_2026-08-25.md

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`
**Review:** [ng_calling_loop_b1_2026-08-25.md](ng_calling_loop_b1_2026-08-25.md)

## 1. Executive summary

**3 Blockers, 7 Majors, 14 Minors, 11 Nits. All 3 Blockers fixed, 6 of 7 Majors fixed, 1 Major
carried forward as a decision for step C1.** The Minors and Nits were applied except where two
agents disagreed or where the change belongs to another branch.

### What actually changed in the code, as against in the prose

Four behaviour changes, and one of them is a real bug:

1. **The total-weight check is now release-held**, and it is the function's only `NaN`
   detector. The finiteness check cannot be one: `largest_score` is assigned only through
   `score > largest_score`, which every `NaN` loses, so a `NaN` in a genotype that is not the
   most probable one never reached it. Before this fix, under `--release`, likelihoods
   `[NaN, 0, 0]` returned normally with every posterior entry and every expected copy `NaN` —
   which the M-step would then sum into the cohort's copies and carry to every other sample's
   next prior. It costs one comparison per sample per pass, not per genotype.
2. **A release-held check that the sample's own expected copies are finite and non-negative.**
   Without it the scratch's `UNWRITTEN_SCRATCH_VALUE` sentinel did not survive to be seen:
   `2.0 − NaN` is `NaN`, `f64::max` returns the *other* operand on a `NaN`, so an unwritten row
   collapsed to a zero leave-one-out term and the sample was scored against the bare seed with
   the cohort's evidence silently absent.
3. **The copy-table width check is `debug_assert_eq!`, not `assert_eq!`.** It cannot fire while
   every `GenotypeTableView` comes from `GenotypeTable::build`, which asserts the same identity
   as it builds — so holding it in release would be a check the suite cannot reach and
   therefore cannot keep honest.
4. **`score_one_sample`, `SampleScoringBuffers` and the accessor are `pub(crate)`**, with the
   dead-code lint `expect`ed only outside the test build. `pub` was buying nothing but silence.

Everything else is a test, a name, or a corrected claim.

### Outcome totals

| severity | raised | fixed | carried forward | declined |
|---|---|---|---|---|
| Blocker | 3 | 3 | 0 | 0 |
| Major | 7 | 6 | 1 | 0 |
| Minor | 14 | 11 | 2 | 1 |
| Nit | 11 | 9 | 0 | 2 |

### Validation

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4540 passed; 0 failed; 14 ignored` — B1 now adds **12** tests to A2's 4,528 |
| `cargo test --release --lib ng::calling::inference::summarise_condition` | `11 passed; 0 failed` |

### Every release-held check is now reached by a test that fails in release without it

This is the property the module's whole no-`Result` design rests on
(`spec/calling_em_loop.md` §8), and before this review two of the four checks did not have it.
**Measured, in one run:** all six release-held checks in `score_one_sample` downgraded to
`debug_assert` together, under `--release` —
`test result: FAILED. 5 passed; 6 failed`. Six checks, six failing tests, a one-to-one mapping:

| check | the test that catches its downgrade |
|---|---|
| `genotype_likelihoods.len()` | `a_short_likelihood_row_is_refused` |
| `posterior_row.len()` | `a_short_posterior_row_is_refused` *(new)* |
| `seed_concentration.len()` | `a_seed_of_the_wrong_width_is_refused` *(new)* |
| `sample_expected_copies` finiteness | `a_sample_whose_own_copies_were_never_written_is_refused` *(new)* |
| `largest_score.is_finite()` | `a_row_with_no_usable_score_is_refused` *(renamed)* |
| `total_weight >= 1.0` | `a_nan_below_the_largest_score_is_refused` *(new)* |

And the row-indexing Blocker: with `sample_scoring_buffers_mut` mutated to ignore its argument,
`the_scoring_buffers_hand_back_that_samples_own_rows` fails and nothing else does
(`511 passed; 1 failed`), where before the fix the whole of `ng::calling` stayed green.

Both source files were restored from pristine copies afterwards and re-verified **by hash**,
not by `git diff --stat`.

## 2. Per-finding log

### BL1 — the row indexing had no test — **Fixed**

Added `the_scoring_buffers_hand_back_that_samples_own_rows` in `calling/mod.rs`'s tests, beside
the accessors it exercises rather than in the scoring module, because the defect is the
accessor's. Its doc carries the measurement of why the three-sample scoring test cannot stand
in for it.

### BL2 — the posterior-row check had no test — **Fixed**

Added `a_short_posterior_row_is_refused`, reached with a hand-built `SampleScoringBuffers` for
the reason `a_short_likelihood_row_is_refused` gives. Its doc records that **neither** of the
three-sample test's assertions notices the truncation, and why.

### BL3 / M1 — the `NaN` hole — **Fixed**

`total_weight >= 1.0` promoted from `debug_assert!` to `assert!`, with a comment saying it is
the only `NaN` detector and why the finiteness check cannot be one. The `largest_score` check
is kept for the `±∞` cases it does cover, and its message now says
*"an infinity, or every score was NaN, since a NaN never wins the comparison that picks the
largest"*. `a_non_finite_score_is_refused` is renamed `a_row_with_no_usable_score_is_refused`
and its doc corrected — **the old doc said it tested `NaN` detection and it did not**; it
reaches the check as the `−∞` the maximum started from. `a_nan_below_the_largest_score_is_refused`
is the new test for the case the old name claimed.

### M2 — the one-sample concentration test could not fail independently — **Fixed (doc)**

Kept, with its limit written on it: `seed + (x − x) = seed` is a property of
`fill_sample_concentration` one module away, so any wiring with `f(x, x) = 0` passes; the
discriminating check on the leave-one-out wiring is the hand-computed test, and the doc now
says so with the mutation that shows it.

### M3 — no test ran a second pass — **Fixed**

Added `at_one_sample_a_second_pass_reproduces_the_first`, which reaches spec §13 test 1's fixed
point without the M-step, because at one sample the M-step is `cohort := own`. Its doc records
the limit the reviewer measured: it still passes with `sample_expected_copies.fill(0.0)`
deleted, so it is a check on the concentration path and not on step 4.

### M4 — the sentinel did not survive to be seen — **Fixed (release-held check + test)**

See §1 item 2. The reviewer left the release/debug choice open; it is held in release, because
the failure it prevents is silent and the module's rule is that a silent caller bug is
release-held. `a_sample_whose_own_copies_were_never_written_is_refused` pins it.

### M5 — `pub` as a dead-code silencer — **Fixed**

`pub(crate)` on all three new items. The lint is `expect`ed rather than `allow`ed, and only
under `cfg(not(test))` where it genuinely fires — so when D1 adds a real caller the expectation
becomes unfulfilled, the build fails, and whoever writes that caller deletes the line. A plain
`#[expect]` does not work: it is unfulfilled in the test build, where the tests are callers.

### M6 — step C1's flat pass — **Carried forward, not applied**

The reviewer's `PassPrior::{Flat, LeaveOneOut}` repair is good and compiles green in their
worktree. It is not applied because the plan puts the flat pass at C1 and the shape of C1's
entry point is C1's decision. The measurement is carried in the review §7 and in
`PROJECT_STATUS.md` so C1 does not rediscover it.

### M7 — the bundle's public fields — **Fixed (visibility + honest doc), constructor declined**

The type is now `pub(crate)`, so the exposure is crate-internal. A checking `new` was declined:
the shapes it could check are exactly the ones `score_one_sample`,
`fill_sample_concentration` and `PriorRow::new` already check, and the doc now says where the
shapes are enforced instead of claiming an invariant the type does not hold. The doc's
"made only by" claim — which the patch's own test contradicted, found by two agents — is
replaced with what is true: every *run-time* construction goes through the accessor, and the
tests reach past it on purpose.

### Minors and Nits applied

- **The module doc no longer opens with `Arm A`**, a label nothing in `src/` defines; it leads
  with what the arm does and names the arch's label in passing.
- **`expectation-maximization` is spelled out** — the module defines *E-step* and *M-step* but
  never connected the letters to anything.
- **`prior_row_workspace` → `prior_per_allele_workspace`**, field and both accessors: it is an
  allele-length buffer that sat two lines from the genotype-length `prior_row` and carried a
  genotype-length word. At a diploid biallelic locus those are 2 and 3, so each is a legal
  length for the other.
- **`sample_scoring_buffers` → `sample_scoring_buffers_mut`**, with the `#[inline]` and
  `#[must_use]` its eight siblings carry. It is the mutable half of five accessor pairs at
  once.
- **The bundle carries the sample index**, so the three shape panics can name the sample. At a
  thousand samples that is the first thing the reader of a message wants.
- **`probability` → `genotype_probability`, `genotype_copies` → `copies_per_allele`.**
- **The `if copies != 0` guard now carries its measurement** — 90 of 126 entries zero at six
  alleles, 8–29% against a stub prior, 0–4% against the shipped one, with the reason.
- **The `# Panics` section is rewritten so each check carries its own reason**, because the
  shared justification was false for the seed check: with it removed, every mis-shaped seed is
  still refused in release a few lines later, so it buys the message and not the catch.
- **The triallelic test's doc no longer calls sum-to-ploidy a second invariant.** It follows
  from sum-to-one whenever the fold walks whole genotype rows; what it catches is a wrong
  *stride*.
- **`a_short_likelihood_row_is_refused`'s fixture** had cohort copies below the sample's own,
  which trips a debug-only check elsewhere — so which assertion it pinned depended on statement
  order. Now `[2.0, 1.0]` against the sample's `[1.0, 1.0]`.
- **The function's summary line says it replaces the sample's expected copies and is not
  idempotent** — the destructive read-then-overwrite was invisible at a call site.
- **"Steps 3 and 4 are where the cohort size does not appear"** led a paragraph whose argument
  is about step 1; now *"No line of this function branches on the cohort size, and step 1 is
  why."*
- **The GIAB figure carries its range** — "each sample called on its own", per `CLAUDE.md`.
- **The `&*sample_expected_copies` reborrow is dropped**; `&mut [f64]` coerces at the argument
  position.
- **"the normaliser" and `total_weight` were two words for one quantity**; now one.

### Declined, with reasons

- **`SampleScoringBuffers` → `SampleScoringView`.** The crate's `*Scratch` types all *own*
  their memory and this one borrows, so `*Scratch` would mislead; `*View` usually implies
  read-only and five of these eight borrows are mutable. `Buffers` is the honest third word.
- **`score_one_sample` → `fill_one_samples_scores`.** The `fill_*` prefix names functions that
  fill one buffer; this fills four and reads a fifth. The reviewer's own second option — put
  the destructive input in the summary line — is what was done.
- **`#[non_exhaustive]` on the bundle.** `Cargo.toml` records the crate as internal-only, so it
  is a no-op, and it would suppress the exhaustive destructure that forces a future ninth
  buffer to be handled at the one call site.

## 3. Carried forward

1. **C1's flat-pass entry point** (M6) — a design decision, with the measurement attached.
2. **`homozygous_alleles()` against `homozygous_allele_for`** — one concept, two names, both
   predating this patch.
3. **`SampleAlleleCopies::new` refusing a non-finite entry in release** — `genotype_prior`'s
   owner's call; B1 works around it with a check of its own.
4. **The four `--release` failures in `ng::calling::likelihood`** — another branch's, and the
   thing blocking a release step in CI.
