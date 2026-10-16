# Review: plan step D7 — the allele-length prior's curvature in the standard errors (2026-10-16)

**Scope.** The uncommitted change in worktree `pop_var_caller-str-penalised-fit` on top of HEAD
`85a93015`: `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs`, its `run_log.rs`, the
re-recorded fixture `bench_fixtures/d_tract_prior_reference.jsonl.gz`, spec §2.6 / §2.10.4 and the
plan's D7 note. Judged against spec `doc/devel/ng/spec/str_penalised_fit.md` §2.4, §2.6, §2.10.4, §7
(2026-10-16) and the diagnosis report
`doc/devel/reports/implementations/str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`. One
reviewer, every category of `ai/skills/rust-code-review/` covered directly, no sub-agents. Every
probe ran in an exported copy (`tmp/review_d7/src_copy/`), native cargo, at most 8 threads, nothing
at 2,000 samples.

**Verdict: correct. No Blocker or Major finding.** The added block is exactly minus the second
derivative of the prior's term the fit maximises, in the coordinates `H` uses; the fitted numbers
are bit-identical to HEAD's; the refactors change no output; the tests fail under every mutation
tried. Two Minor findings (one wrong line in the spec's decision record, one degenerate
configuration that silently skips the block) and four Nits.

## What was verified

### 1. The block is minus the second derivative of the prior's term

**What the fit maximises.** Each round of counting adds `q = allele_pseudo_draws +
allele_prior_draws_from_loci(…)` draws to every one of a class's lengths (`penalised_fit.rs:1414-1428`)
and fits the share logits to those counts as a multinomial (`CurveData::evaluate`, Multinomial arm,
`penalised_fit.rs:2215-2238`). The prior's part of that objective is `T(η) = Σ_a q·ln s_a` over all
lengths, the reference's included, with `s` the softmax of the logits and the reference's logit
held at 0.

**The algebra.** `∂ ln s_a/∂η_b = δ_ab − s_b`, so `∂T/∂η_b = q − Q s_b` with `Q = q × alleles`, and
`∂²T/∂η_b∂η_c = −Q s_b (δ_bc − s_c)`. Minus that is `Q (s_b δ_bc − s_b s_c)`, which is what
`allele_prior_curvature` (`penalised_fit.rs:321`) fills, over the non-reference lengths in the
order `softmax_with_reference` assigns logits. It is also the prior's share of the information the
M-step itself uses (`total × (diag(s) − s sᵀ)` with `total` including the prior's draws), so `H`
now holds the curvature of the objective the counting climbs.

**Independent checks.**

- A Python port of the block against central second differences of `T` (step 10⁻⁴), 200 random
  cases at 3, 9 and 19 lengths, the reference at every position, `q` from 0.01 to 5, and one share in
  five pushed to its floor (logit −25): largest error 2.0 × 10⁻⁶ of `max(q, 1)`
  (`tmp/review_d7/scripts/fd_check.py`).
- A reviewer test in the copy (`review_d7_probe_no_second_floor_and_every_class`), on the test
  cohort's fitted period, **every** class rather than the one the shipped test checks, with `q`
  recomputed by hand from the class's last counts: largest error 1.7 × 10⁻⁸ of `max(q, 1)` with
  `allele_prior_tracts = 0` (the bench's method `D`, which keeps only the first floor) and
  9.1 × 10⁻⁸ with the production value 1.

**Edge cases.**

| case | what the code does | right? |
|---|---|---|
| class with no reads | `class_at_point` is `None`: no counts, no block, no shares in `H` | yes, as counting |
| `tracts_with_reads() = 0` at a class with counts | `.max(1)` inside the shared helper, as counting did before | yes, same bits |
| `allele_prior_tracts = 0` | helper returns 0; `q` = first floor only | yes (probe above) |
| `allele_prior_tracts < 0` | helper returns 0, as counting did | same as counting (no validation; pre-existing) |
| a share at its floor | entries go to 0 smoothly; the block stays positive semi-definite | yes (Python check) |
| counts reused from a settled class | `last_counts` holds the clone counting used | yes |
| `max_rounds = 0` | `last_counts` is `None`: block skipped | see Minor finding 2 |

**`Q` held fixed.** The block treats `q` as a constant, as one round of counting does. See Nit 1 for
what that leaves out.

### 2. Counting arithmetic unchanged to the bit

The counting change replaces an inline expression with `allele_prior_draws_from_loci`. The
operations and their order are the same (`allele_prior_tracts * total / tracts.max(1) as f64`, then
`/ draws.len().max(1)`), and the fixture confirms it:

- **HEAD's fixture against the new one, every leaf** (`tmp/review_d7/scripts/cmp_fixture.py`): of
  152,442 values outside `errors` and timings, **none moved**. Of 64,953 error values, 59,443 moved
  and 63 went from absent to present; the 63 are the share-error vectors and three stutter errors of
  run 147 (`p2_n60_d10_pcr_free` seed 2), which had every error but F's withheld at HEAD.
- **F's errors**, read off pairs of samples with their own curvature: `between_locus_share` and
  `concentration` errors bit-identical in all 2,723 classes (`cmp_errors.py`). (No fixture run counts
  F inside the counting, so that branch's F error is not compared here; it moves by design.)

### 3. The refactors change no output

`standard_errors` → `joint_curvature` + `errors_from`, `fit_one_period` → `fitted_period`, and the
run-log tail → `outcomes_and_lines_of_a_motif_length` read as moves. To check that rather than
assume it, the copy was built with an environment switch that skips the new block, and the fixture
grid was re-run:

- **Block skipped, refactored code, against HEAD's fixture**: 144 of 144 runs **exact**, the
  withheld run included — every cell below 2,000 samples, 3 seeds each
  (`tmp/review_d7/runs_omit/`).
- **Block in, against the re-recorded fixture**: 144 of 144 runs exact (`tmp/review_d7/runs_new/`).

So the only change in output is the block's, and the re-recorded fixture is what this code
writes. The 9 runs at 2,000 samples were not re-run (the review's limit).

### 4. Tests, and whether they would catch a wrong block

Mutations in the copy, switched by environment variable, each run against the 169 `ssr_fit`
tests (`tmp/review_d7/mutation_results.txt`):

| mutation | failing tests |
|---|---|
| none | 0 of 169 |
| block omitted from `P` | `a_cohort_whose_reads_alone_curve_upward_gets_its_errors` |
| block added twice | `a_cohort_whose_reads_alone_curve_upward_gets_its_errors` |
| errors use the first floor only (second floor dropped from `q`) | both new tests |
| `Q = q × (alleles − 1)` (reference length left out) | `the_allele_length_prior_curvature_is_the_second_derivative_of_its_term` |
| block scaled by 1.1 | the same |
| off-diagonal entries dropped | the same |

Every mutation is caught. The fallback remains reachable and reported: the built-to-fail matrix
(`a_fit_whose_curvature_cannot_be_inverted`, `penalised_fit.rs:4009`) drives
`UndeterminedErrors::EveryNumberButBetweenLocusShare` and `EveryNumber`, and the run-log test
checks the line's wording and that every error is written as not placed with F kept.

### 5. Tooling

- `cargo fmt --check` on a clean copy: no output, exit 0.
- `cargo clippy --release --features bench-fixtures --all-targets`: no hit in a changed line; in
  `penalised_fit.rs` only the two pre-existing ones at 3115-3116 (test code).
- `cargo test --release --features bench-fixtures --lib ssr_fit`: `169 passed; 0 failed; 5 ignored`.
- Float portability: the block is products and sums in a fixed order; no transcendental function,
  no thread-order dependence.

## Findings

### Minor

#### Spec §7 says "the half draw a class"; the floor is one draw a class

- **Location:** `doc/devel/ng/spec/str_penalised_fit.md:846`
- **Problem:** The owner's 2026-10-16 answer describes §2.4's floor as "one locus's worth of draws
  spread over the lengths, plus the half draw a class". The code adds **one** draw a class:
  `ALLELE_PSEUDO_DRAWS_A_CLASS = 1.0` (`penalised_fit.rs:355`), spread as 1/19 a length at ±9
  (`pseudo_draws_a_length`), and §2.4 (line 188) and the config's doc say one draw a class. The
  config's doc even explains why half a draw a *length* was rejected.
- **Why it matters:** A reader checking the block against §7 would look for a half draw that is
  not there.
- **Fix:** Leave the owner's words and add a correction after them, e.g.
  `*(The floor is one draw a class, `ALLELE_PSEUDO_DRAWS_A_CLASS`; corrected at D7.)*`

#### With `max_rounds = 0` the block is skipped and the fallback is reached

- **Location:** `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:1812`, and the claim at
  `:603` ("no case is known that reaches it").
- **Problem:** The block reads `q` from `self.last_counts[point]`, which is `None` until a round of
  counting has run. With `max_rounds = 0` — a setting the code allows and documents (`:461`, "`None`
  only when no round ran") — the `let … else { continue }` skips every class, so `H` again lacks the
  prior. Measured in the copy on the test cohort: `undetermined_errors =
  Some(EveryNumberButBetweenLocusShare)`, every stutter error absent.
- **Why it matters:** Low: no production path sets `max_rounds` to 0, and errors at unfitted
  starting values mean little. But the silent skip contradicts "the prior the errors account for is
  the prior the fit used", and the fallback's doc says no case reaches it.
- **Fix:** Either `debug_assert!(self.last_counts[point].is_some())` beside the `let … else` with a
  comment that only `max_rounds = 0` gets there, or one line in `UndeterminedErrors`' doc: "except a
  fit that ran no round (`max_rounds` of 0), whose `H` has no allele-length prior".

### Nits

1. **`Q` held fixed is a convention, not the exact Jacobian** (`penalised_fit.rs:317-320`). The
   second floor's `q` is proportional to the class's expected allele draws, which come from the
   E-step at the previous round's shares and concentration (`ClassCounts::draws`, "how many
   independent allele draws", which the Pólya clustering makes depend on F). The fixed point's
   exact derivative would add a term from `∂q/∂θ`. It is at most one locus's worth against the
   class's own draws, and counting holds `q` the same way, so no change is needed; the doc's "it
   comes from the class's expected draws, not from the shares being chosen" could say "held at the
   last round's value, though those draws themselves depend on the numbers". Not measured
   (confidence: Low).
2. **`JointCurvature` holds `A`, `B` and `P`** (`penalised_fit.rs:2911`); only two of the three are
   curvatures. `SandwichMatrices` or `ErrorMatrices` would say what it is.
3. **No shipped test of the block at `allele_prior_tracts = 0`** (the bench's method `D`). The
   reviewer probe above passes there (1.7 × 10⁻⁸); worth keeping as a second case of the
   second-derivative test, or by looping that test over both values.
4. **Plan line 388** still names the old test `a_three_base_cohort_whose_curvature_cannot_be_inverted`
   in the step's text; the *Implemented* note says it was renamed, which is enough. (Out of scope,
   pre-existing: spec §2.4 line 189 cites `penalised_fit.rs:57-75` for the floors, which now sit at
   74-97.)

## Docs checked

- `UndeterminedErrors` doc (`:598-620`): the 56 of 384 runs at 500 samples or fewer and the test
  cohort match the diagnosis report (its lines 12, 25, 153). Accurate, apart from Minor 2.
- `standard_errors` doc: the description of `P` matches the code; the removed sentence ("the half
  allele draw a class adds is left out of `H` and `P`") is no longer true and was rightly removed.
- Spec §2.6's new paragraph: formula, coordinates and "89 of D1's 480 runs" (56 of 384 plus 33 of
  96 at 2,000) consistent with the diagnosis and the implementer's report; the 2,000-sample claim
  was not re-run here.
- Run-log line and its doc: wording matches the test; "no cohort is known to fail" — see Minor 2.

## Scratch

All under `/Users/jose/devel/pop_var_caller-str-penalised-fit/tmp/review_d7/`:

- `src_copy/` — the export with the mutation switches (`REVIEW_D7_MUTATION=omit|double`,
  `REVIEW_D7_HELPER=pseudo_only|no_reference|scale|diagonal_only`) and the probe test
  `review_d7_probe_no_second_floor_and_every_class`; `src_clean/` — unmodified export for fmt/clippy.
- `mutation_results.txt`; `runs_new/`, `runs_omit/` (fixture grid re-runs); `bin/bench_d7`,
  `bin/bench_mut`.
- `fixture_head.jsonl.gz`, `fixture_new.jsonl.gz`; `scripts/cmp_fixture.py`, `scripts/cmp_errors.py`,
  `scripts/fd_check.py`.
- Build trees `target/`, `target_clean/` (large; safe to delete).
