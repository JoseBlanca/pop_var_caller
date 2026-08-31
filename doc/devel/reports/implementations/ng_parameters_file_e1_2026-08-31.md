# ng parameters file — E1: the defaults as named constants with their origin

**Date:** 2026-08-31
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone E, step E1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3.3, §3.8, §5, §8
**Code:** `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` in
[likelihood/mod.rs](../../../../src/ng/calling/likelihood/mod.rs); the new
[parameters_file/defaults.rs](../../../../src/ng/calling/parameters_file/defaults.rs); the file's
own header and one row comment in
[to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs)

---

## 1. What it does

Spec §8 asks that "the default for every parameter is a named `pub const` in the source with its
origin recorded beside it". Two of the three numbers §8 names already were —
`DEFAULT_OUTLIER_WEIGHT` and `STATED_FLAT_CONCENTRATION`. The third, the base-quality calibration's
scale of one, was a bare `1.0` inside `ReadGroupCalibration::defaulted()`. It is now
**`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`**, named where its reader is, with its origin above it,
and `defaulted()` reads it.

`defaults.rs` is the one place a person can find the whole set. It holds no code: the numbers stay
beside the code that reads them, because a constant re-declared in a summary is a second spelling
that can disagree with the first. What it holds is the inventory and the tests that pin each
default to its warrant at the point of use.

## 2. The inventory is seven rows, not §8's four

§8 names three constants plus contamination's absence. Walking the nine fields of `RunParameters`
found three more things a run with no fit needs, and `to_toml.rs`'s own `origins` module — written
at B3 *"so that step E1 has one list to reconcile against"* — already carried all three:

| what a run needs | what it takes with no fit |
|---|---|
| the base-quality multiplier, per read group | `DEFAULT_ERROR_PROBABILITY_MULTIPLIER`, one |
| the repeat-tract outlier weight, one per run | `DEFAULT_OUTLIER_WEIGHT`, 0.01 |
| the tract ladder's fallback concentration | `STATED_FLAT_CONCENTRATION`, one |
| contamination, per read group | absence — no `[contamination]` section |
| the repeat-tract substitution rate | `DEFAULT_SSR_SUBSTITUTION_RATE`, 0.001, taken at the tract |
| the slippage numbers | **no row** — the tract falls back to `StutterModel::hipstr_shipped` |
| the inbreeding coefficient, per sample | **nothing** |

**The last two are not alike and the header says so.** The slippage numbers are *owed a
measurement* (§8's third bullet, §12 question 1). The inbreeding coefficient is *forbidden a
default*: `parameter_estimation::generic::fallback`'s header states the rule — "it is the parameter
that differs most between an outcrosser and a selfing landrace, and a cohort's diversity divides by
`1 − F`, so a wrong constant would be amplified rather than absorbed" — and `origins`'s own text
for it reads *"a run should not be able to write this line"*.

**The prior's seed is deliberately outside the table.** It has a fallback
(`ExpectedHeterozygosity::SPECIES_FALLBACK`) but what records it is `SeedRegime::FallbackDiversity`,
a rung in the file's `[ordinary_site_prior]` section, not a `warrant` — so it is marked, and not
marked the way the others are. §8 does not discuss the difference.

## 3. ⚑ The rung this step added and then removed, and why it must not come back

E1's first draft added a rung to `validate`: *a `defaulted` base-quality multiplier is
`DEFAULT_ERROR_PROBABILITY_MULTIPLIER` and nothing else*, on the model of the two identical rungs
the outlier weight and the fallback concentration already carry. **It refused a file this caller
had just written.**

**The two rules are not the same rule.** On the outlier weight and the fallback concentration,
`defaulted` means *the run took this compiled-in number*, so the warrant determines the value. On
the multiplier it means something one level up: the multiplier is a fitted error **rate** divided
by the geometric mean of that read group's minted error, and `from_fitted_rate` copies **the
rate's** warrant onto the ratio. The pre-pass's error-rate ladder has a `Defaulted` bottom rung of
its own — `DEFAULT_ERROR_RATE` at 0.001, taken by a read group with fewer than `MIN_SITES_TO_FIT`
sites of its own, no sibling above that floor to borrow from and nothing supplied — so a legitimate
run writes a `defaulted` multiplier of `0.001 / that library's mean minted error`, which is one only
by coincidence.

**Measured**, on the fixture the regression test now carries: two libraries reporting a mean error
of 2.5 × 10⁻⁴ and 5 × 10⁻⁴ give multipliers of **4.0 and 2.0**.

**The reach is the corner `CLAUDE.md` names as the hardest.** `MIN_SITES_TO_FIT` is 10,000 *sites*,
so the ordinary way in is a run over one gene, a panel, a small contig, or a shallow single sample
where no read group clears the floor. And the failure is not a wrong number — it is the run's whole
parameters file becoming unreadable, so the run cannot be re-used or fed back in.

`a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts` is what stops the rung
coming back, and `validate`'s comment says at length why it is absent.

**⚑ It leaves a divergence, and the owner ruled it on 2026-08-31: the code stays and spec §5's third
row is the sentence to correct.** §5 says a read group whose rate could not be fitted gets *"scale
1.0, warrant `Defaulted`"* — take its reads at the quality they claim. It should not: **a library's
real error rate is never its reported sequencing quality**, because the quality scores describe base
calling and the reads also carry mismapping, chimeras and damage. So a library nothing could be
fitted for is charged a stated rate rather than taken at its word, and on a real library that is the
conservative direction — at HG002's measured mean minted error of 2.9055 × 10⁻⁴
([`read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2) the multiplier is **3.44**, so
every read is charged 3.4 times worse than it claimed, 5.4 Phred less confident. `DEFAULT_ERROR_RATE`
at 0.001 is a placeholder and is to be fitted from GIAB, like §8's slippage numbers.

**⚑ The direction was reported backwards before the ruling, and the fixture is why.** The first
draft of the regression test used libraries reporting 0.008 and 0.004 — chosen to make the two
multipliers distinct, not to be realistic; 0.008 is Q21, which no library in this project is
anywhere near. Those gave multipliers *below* one, so the report and the chat summary said a
defaulted rate makes reads **more** confident, which is true of that fixture and false of any real
library. The fixture now reports 2.5 × 10⁻⁴ and 5 × 10⁻⁴, either side of HG002's measured value, and
the test asserts both multipliers are above one — so the fixture cannot be read the wrong way round
again.

## 4. What was absorbed beyond the step's own brief

Four message changes, all found by the review agent that read the produced file as a geneticist,
and all in the two functions this step was already editing:

- **The three `defaulted` refusals now close with one clause, word for word.** The fallback
  concentration's was transposed — *"is `defaulted` at 3.5, and the constant a run falls back to is
  1.0; a number **somebody chose** is `supplied`"* — where the other two say *"a number **you
  changed** is one the run was handed, so change the warrant beside it to `supplied`"*. It stated a
  fact about the format where the others give an instruction. Two tests assert `ends_with` on the
  shared clause now, rather than three substrings both shapes satisfied.
- **Four messages spell their floats the way the writer does.** `to_toml` formats every value with
  `{:?}`, so the file contains `0.0` and `1.0`; four refusals printed `is 0,` and `is 1`, which a
  reader cannot search their own file for. The fixture values in the test were changed to whole
  numbers (3.0, not 3.5) because `Display` and `Debug` agree on 3.5 — measured: reverting the
  `{:?}` passed 1,852 tests before the fixture moved.
- **"every read of the library" became "every read of that read group"** in the zero-multiplier
  refusal, whose key is `by_read_group[...]`. The file spends a paragraph on two lanes of one
  library being different rows; two adjacent messages should not use the two words as synonyms.
- **The row comment for a defaulted multiplier reads "not calibrated: … taken at face value"**
  rather than "no calibration: … used exactly as they came". The reader read the old one as *this
  lane is untrusted*, which is the opposite of what it means.

And two changes to the file's own header, which the step made false and then had to make true:

- it still says **"Two keys do not take every warrant"**, because the rung was removed;
- a new paragraph says **what that checking reaches and what it cannot** — that a number you
  changed still labelled `fitted_here` is accepted silently with its old `observations` beside it;
  that `defaulted` is checked against nothing on every other key; and that on the multiplier there
  *is* a built-in number and the key still is not checked, with the reason.

The second is the answer to a review finding worth recording: a reader refused once for leaving a
`defaulted` warrant behind concludes *the caller checks my warrants*, and then makes the same class
of mistake on a `fitted_here` row with nothing to stop them.

## 5. Tests

Seven added: six in `defaults.rs`, one in `from_run_parameters.rs`.

- `a_read_group_with_no_fitted_rate_is_charged_what_its_reads_were_minted_with` — the constant, the
  warrant, and that `charged_error` at a scale of one is the geometric mean of the reads' own error
  probabilities, bit for bit.
- `the_outlier_weight_a_run_inherits_is_the_existing_callers_number`, and
  `a_run_that_fitted_no_stratum_seeds_every_tract_from_the_flat_rung` — the other two constants,
  their warrants, and that the flat rung hands out no shape of its own.
- `a_file_with_nothing_fitted_for_repeat_tracts_is_still_a_legal_file` — no slippage row, no
  substitution-rate row, no contamination section. **This is what a defaults run will write and
  what E2 rests on**; nothing held it before.
- `a_defaulted_value_that_is_not_the_binarys_own_number_is_refused` and
  `each_of_the_three_constants_is_accepted_beside_a_defaulted_warrant` — the two keys the file's
  reader can hold to a constant, and the third that it cannot.
- `a_run_whose_rates_were_defaulted_writes_a_file_its_own_reader_accepts` — §3 above.

**Mutations.** The review ran fifteen against the first draft, thirteen killed; the two survivors
were the unpinned `{:?}` and are pinned now. Since the fixes: making `from_fitted_rate` stamp
`FittedHere` instead of copying the rate's warrant fails the new regression test (and 22 others);
moving `DEFAULT_ERROR_PROBABILITY_MULTIPLIER` to 1.5 fails 75 tests.

## 6. Validation

In the container, on the committed tree:

- `cargo test --lib` — **5,548 passed, 0 failed, 13 ignored** (5,541 before this step).
- `cargo test --lib ng::calling::parameters_file` — **177 passed** (170 before).
- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo doc --no-deps` — 25 unresolved-link **errors** and 23 redundant-target warnings, both
  unchanged pre-existing baselines. Under `--document-private-items` the new module contributes no
  diagnostic; every one of its links is a full path, because a module whose only imports are inside
  its test block has nothing in scope for a short one.
- `cargo test --all-targets` still exits 101 on the pre-existing panic in `benches/psp_writer_perf.rs`.

## 7. Review

Three agents in isolated worktrees: correctness with mutation testing, design fidelity, and the
produced file read as a geneticist. **One Blocker (§3), four Majors, and about twenty Minors, and
every finding but the Blocker was prose that said something untrue.**

**A correction to how a review runs in a worktree, which two sessions have now had wrong.**
`./scripts/dev.sh` from the main repo bind-mounts only its own `PROJECT_DIR`, so a `cd` into the
worktree fails inside the container without stopping the script, and cargo compiles the **main
repo**. A mutation written into the worktree is therefore compiled away and reports as surviving.
The form that works is the **worktree's own** copy of the script, which computes `PROJECT_DIR` from
its own location:

```
<worktree>/scripts/dev.sh bash -c 'cd <worktree> || exit 9; cargo test --lib ng::calling::parameters_file'
```

Verified both ways in this session. Recorded in `PROJECT_STATUS.md`, where the old form stands.
