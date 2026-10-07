# Code Review: fit_precision_a2
**Date:** 2026-09-27
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A2 of `fit_precision.md` — the information, summed in blocks by an expectation pass
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `fd0b364e..582c8b4c` (`582c8b4c`, "WIP: fit precision A2, for
  review", a review-only commit on no branch).
- **In-scope files:** [src/parameter_estimation/joint/fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
  (`PassModel`, `one_position`'s signature, `Statistics.information`, `expectation_pass`),
  [src/parameter_estimation/joint/fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (`InformationSums`, the rules' slopes once a pass, the scorer split by branch, the read-slope
  table, three new tests), the A2 implementation report.
- **Out of scope:** the rest of `fit.rs`, read only as callers and callees.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A2 and spec §3.2, §3.6 item 5; the prose;
  - *reliability* — tests and mutation (24 mutations);
  - *correctness and determinism* — errors, float portability, concurrency (the chunk joins), extras
    (the cost of a collecting pass); a bit-for-bit comparison of A1's and A2's scorers and of whole
    fits before and after.

## 2. Verdict

**Approve-with-changes.** No fitted number moves: A2's scorer gives the same bits as A1's at every
one of 600 positions in 42 configurations, and whole `fit_jointly` dumps are byte-identical before and
after on three drawn cohorts, one of them a single sample. The joins are deterministic and the diff
adds no platform-dependent arithmetic. All three reviewers found the same gap: the way the final pass
joins its chunks — the path A5 turns on — had no test of its values, and the thread test could not
tell an absent result from a present one.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 36 passed at
  `582c8b4c`; both `cli::cross_platform_digests` tests pass.
- Mutations: design 13 run (5 killed, 2 survived with a behaviour change, 6 no behaviour change);
  reliability 24 run (17 killed, 2 survived with a behaviour change, 1 equivalent, 4 no behaviour
  change); correctness 5 run (3 killed, 1 survived with a behaviour change, 1 no behaviour change).
  Every survivor with a behaviour change is on the ordered join (M1).
- Findings labelled *Needs verification*: 0.

## 4. Open questions and assumptions

1. **Spec §3.2 and plan step A2 still describe a per-sample block of (k + 8)² with k = 1 + 2 × read
   groups** and name `SampleInformationBlock`. What is built has k = 3 and keeps the cohort's block
   once. Only the A2 report records the change; step A3 is specified from the spec and the plan. The
   design documents are the owner's to amend, with the read-group ruling at checkpoint A. Affects Mi5.

## 5. Top 3 priorities

1. **M1** — test the values on the final pass's way of joining, and make the thread test fail on an
   absent result.
2. **Mi1** — `Statistics::absorb` drops the information silently when two chunks disagree on
   keeping it.
3. **Mi6** — the debug checks that the scorer's tables match its rules are never exercised.

## 6. Findings

### Major

#### M1: src/parameter_estimation/joint/fit/information.rs (tests) — no test checks the information on the final pass's way of joining chunks
- **Categories:** design, reliability, correctness (convergent)
- **Confidence:** High
- **Problem:** a pass that keeps per-position lists — the final pass — joins its chunks in position
  order; the others use a fixed halving tree. `a_pass_sums_the_scores_multiplied_pairwise` checked
  values only on the tree. `the_information_is_the_same_bits_at_any_pool_width` compared
  `format!("{:?}", information)` at one and four threads, which is `"None" == "None"` when the
  information is lost. Setting it to `None` on the ordered path, or keeping only the first chunk's
  sums, survived all 36 tests.
- **Suggested fix:** run the value test on both joins; `expect` the information in the thread test.

### Minor

#### Mi1: src/parameter_estimation/joint/fit.rs `Statistics::absorb` — a mismatch in carrying the information is dropped silently
- **Categories:** design, reliability, correctness (convergent) · **Confidence:** High
- **Problem:** `if let (Some, Some)` does nothing when one side is `None`, so the other chunk's sums
  vanish, or the total covers part of the census. Unreachable today; the one way to lose the
  information on the ordered path is exactly this.
- **Suggested fix:** a `match` that is `unreachable!` on a mismatch (or a `debug_assert_eq!`); one
  `Option` deciding whether a position is scored.

#### Mi2: src/parameter_estimation/joint/fit.rs `expectation_pass` — two positional bools, and a new `too_many_arguments` allowance
- **Categories:** design · **Confidence:** High
- **Suggested fix:** a small struct naming what the pass keeps.

#### Mi3: src/parameter_estimation/joint/fit.rs (two doc comments) — "only on a run's last pass"
- **Categories:** design, correctness (convergent) · **Confidence:** High
- **Problem:** nothing collects the information at this commit; from step B1 the passes that refresh
  the errors will too, through the tree join.

#### Mi4: src/parameter_estimation/joint/fit/information.rs — the read-slope table's index written twice, its size a third time
- **Categories:** design · **Confidence:** High
- **Suggested fix:** one index function and one length function.

#### Mi5: spec §3.2, plan step A2 — the design text describes a block shape that was not built
- **Categories:** design · **Confidence:** High
- **Problem:** see open question 1.

#### Mi6: src/parameter_estimation/joint/fit/information.rs `score_position` — the debug checks are never exercised
- **Categories:** reliability; correctness (noted) · **Confidence:** High
- **Problem:** no test hands the scorer tables built for another rule; a check that always passes
  survives. The checks compile out under `--release`.
- **Suggested fix:** a `#[should_panic]` test in debug builds.

### Nits

- `RuleSlopes::describe` returns a `bool`; `describes`.
- The per-branch shares are rebuilt at every position, though documented as the same all pass.
- `own` / `with_cohort` are bare adjectives as field names.
- `ClassAt` and its binding `at` sit beside the unrelated `at_node`.
- `&super::Parameters` written inline while every other item is imported.
- `InformationSums` derives `PartialEq`, which nothing uses.
- `InformationSums::absorb` would truncate on a length mismatch; state the precondition.
- "Two equal strings are two equal sets of bits" holds for finite values.
- "Transposed a block fails it" holds only for the sample-with-cohort block; the other two are
  symmetric.
- Spec §3.6 item 5 asks for 1, 4 and 8 threads; the test checked 1 and 4.
- The report cited review-internal finding codes; and its "450 MB at 2,169 samples" depends on
  kimura's 2,651 read groups (7,478 parameters).
- "At genotype 0 no read's slope depends on the candidate" holds in exact arithmetic, not bit for bit.
- The zero-row skip in `add_position` is a time saving, not a change in what is added.
- The "moves no count" check compares four fields, not every count.

## 7. Out of scope observations

- **The underflow early return in `score_position` is reached by no fixture.** A position whose
  likelihood is `−∞` would otherwise put `∞ × 0 = NaN` into every block. `one_position`'s own comment
  says such positions do not arise from well-formed evidence. Follow-up: a test when one can be
  built cheaply.
- **At one sample the homozygote excess has non-zero information** (90.0), though the fit never moves
  it; step A3 must report it absent. **With the duplicated class off** the cohort block's last three
  rows and columns are zero; A3 must restrict to the parameters with information.
- **The two joins agree to 1.66 × 10⁻¹⁵ relative, not bit for bit**; A5 must not compare them
  exactly.
- **Cost** (correctness, one thread, depth 3): a collecting pass takes 1.41× a plain pass at 4
  samples and 1.72–1.75× at 64. The final pass holds every chunk's sums at once, about 73 MB at
  kimura.

## 8. Missing tests to add now

- The value test on both joins, and the thread test requiring a present result (M1).
- `tables_built_for_another_rule_are_refused` — debug builds; tables built at a shape `a` 0.01 higher
  than the pass's (Mi6).
- The value test on a one-sample cohort (correctness, cross-category).

## 9. What's good

- `PassModel::new` is the old setup code moved without an edit, and `one_position` unpacks it into the
  same local names, so the refactor is visibly arithmetic-free — confirmed bit for bit by the review.
- The information follows the pass's existing fixed chunk bounds and joins; nothing is sized from the
  pool.
- The value test measures each entry against the Cauchy–Schwarz bound of its row and column, so
  cancelling off-diagonal sums are held to a meaningful scale.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
- `scripts/dev.sh cargo test --lib tables_built_for_another_rule_are_refused`
- `scripts/dev.sh cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`;
  `cargo test --all-targets --all-features --no-fail-fast`

Per-category files: `tmp/review_2026-09-27_fit_precision_a2/` (`design.md`, `reliability.md`,
`correctness.md`).
