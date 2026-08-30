# Review — ng calling loop, C2: convergence, the cap, and the emitted flag

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, reviewed at `52937a76` + the C2 working-tree diff
**Implementation report:** [C2](../implementations/ng_calling_loop_c2_2026-08-25.md)

## 1. Scope

The uncommitted C2 diff to
[`summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs):
`cohort_copies_have_settled`, `FrequencyLoopOutcome`, `run_frequency_loop`, their tests, and
three edited `expect(dead_code)` reason strings.

**Out of scope:** the rest of that file (B1, B2, C1 — reviewed and committed);
`src/ng/calling/allele_candidates/`; the rest of `src/`.

**Five agents, each in its own git worktree**, detached at `52937a76` with the diff applied as a
patch: reliability, errors, naming, idiomatic + smells, and a dedicated pass over the diff's own
quantitative claims (the skill's step 8a).

## 2. Verdict

**Request-changes**, all applied. **1 Blocker, 4 Majors, 8 Minors, 2 wrong sentences.**

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4613 passed; 0 failed` at review time |
| `cargo test --release --lib ng::calling --all-features` | `567 passed; 0 failed` at review time |

## 4. The findings that mattered

### B1 — the per-sample inbreeding pairing was untested (reliability)

Every C2 fixture reached the loop through a helper that built `vec![outbred(); n]` — **one
coefficient, repeated** — so nothing in the file could see an implementation that scored samples
against a coefficient other than their own. Replacing `inbreeding` with
`inbreeding_by_sample[0]` left the suite green.

The mutant is not a no-op. On the three-sample fixture with coefficients `[0.0, 0.5, 0.9]`:
shipped gives `passes: 3` and cohort copies `[2.883, 3.117]`; the mutant gives `passes: 4` and
`[3.264, 2.736]` — **0.38 copies of the alternative allele out of six chromosomes, an
allele-frequency shift of 0.064**. The function's own doc comment called the length check
"load-bearing rather than defensive" because a mis-paired coefficient is finite, plausible,
wrong and silent; the length check was tested and **the pairing it protects was not**.

Fixed by `each_sample_is_scored_against_its_own_inbreeding_coefficient`, which trades two
samples' coefficients *without moving the samples* — the one move a shared-coefficient
implementation cannot tell from the original.

### M1 — `ploidy` was a second source of truth, and a mismatch was silent (errors)

`run_frequency_loop` took a `Ploidy` argument used on exactly one line, while the
`GenotypeTableView` it also takes already carries the answer — the table was *built* for a
`(ploidy, allele count)` shape. Nothing compared them, and the C3/D1 caller will hold two
(`FrozenParameters::ploidy()` and the table's).

Measured on the three-sample fixture with a diploid table: `ploidy = 2` gives
`passes: 4, converged: true`; **`ploidy = 64` gives `passes: 2, converged: true`** with a
different cohort row, identically in debug and release, with nothing asserting. A too-large
ploidy loosens the threshold by the ratio and claims convergence it did not reach; a too-small
one tightens it and sends loci to the cap, where §6 emits them as the weaker claim.

Fixed by **removing the parameter** rather than asserting agreement — the loop reads
`genotypes.ploidy()`.

### M2 — `.abs()` could not fail at any fixture in the file (reliability)

Deleting it left the suite green, and the reason is structural rather than accidental: expected
copies sum to the cohort's chromosome total on every pass, so **at two alleles the two movements
are exactly equal and opposite** and a signed comparison gives the same verdict as an absolute
one. Every C2 fixture was biallelic, and the division test's hand-built row moved only upward.

From three alleles on the mutant is real: the reference allele falling by more than the
threshold while every alternative rises by less than it — scaled `[−0.0015, +0.0008, +0.0007]`
against `1e-3` — is the ordinary shape of a multi-allelic locus that is still moving, and this
caller meets three and four alleles routinely.

Fixed by `a_fall_larger_than_every_rise_has_not_settled`, plus
`the_loop_settles_at_three_alleles_and_the_copies_still_sum_to_the_chromosomes`, because **no
C2 test ran the loop past two alleles at all**.

### M3 — nothing pinned that the settled test wins over the cap test (reliability)

Swapping the loop's two exit tests left the suite green: no fixture put the cap at the pass on
which its locus settles. At a cap of exactly 4 — the pass the three-sample locus settles on —
shipped reports `converged: true` and the mutant `converged: false`. §6 makes the flag a claim
about the *locus*, so reporting the cap there understates every genotype at the site.

Fixed by `run_frequency_loop_reports_converged_when_the_last_allowed_pass_settles`.

### M4 — two wrong mechanism sentences, both mine (claims pass)

**40 claims re-derived, 38 correct.** Every number about a fixture, a mutation or a test count
was right; both failures were prose explaining *why*.

1. *"`u8` and `usize` both widen to `f64` exactly at every value either can hold"* — `usize`
   does not: `9007199254740993usize as f64 == 9007199254740992.0`, measured. The same sentence's
   second half was correct and contradicted the first.
2. *"Both sides are given a threshold of `1e-300`, so neither can stop early"* — the hand-driven
   side is `for _ in 0..PASSES`. It has no threshold, no cap, and never calls the stopping rule
   at all, which is a **better** reason for the test than the one written.

Both rewritten.

## 5. Minors, and what became of them

| finding | category | outcome |
|---|---|---|
| the too-many-coefficients direction of the length check untested | reliability | `more_inbreeding_coefficients_than_samples_are_refused` |
| `is_finite()` half of the value guards unfalsified — `0.0` and `NaN` are refused by `> 0.0` alone | reliability | `an_infinite_threshold_is_refused`; the infinite *chromosome count* is now unrepresentable (below) |
| `< threshold` vs `<= threshold` untested; `2e-3 / 2.0` is exactly `1e-3` in `f64` | reliability | `a_movement_exactly_at_the_threshold_has_not_settled` |
| two guard tests matched two substrings of one assertion message | reliability | the assertion is now two, each with its own message |
| `FrequencyLoopOutcome` not `#[must_use]` — the only carrier of "did not settle" | errors, idiomatic | applied, to the struct and the predicate |
| `cohort_chromosomes: f64` and `threshold: f64` transposable | naming, idiomatic | the function now takes `Ploidy` and a sample count |
| names drop the crate's word *expected* | naming | renamed `cohort_expected_copies_have_settled`, `previous_/current_expected_copies` |
| `never_settles_before` names a property its fixture lacks | naming, smells | renamed `settles_only_at_a_bitwise_fixed_point` |
| `loop_over`'s `seed` reads as an RNG seed | naming | `seed_concentration` |
| the doc claimed the sentinel row is "the row this function is first handed" | smells | it is not — the initialisation's M-step writes finite copies first; reworded to say the guarantee is for C3's and D1's callers |
| the guard test's comment said a short slice "would panic on an index" | errors | it does not; comment corrected and the fixture given real likelihoods |
| `passes: u32` can hold a zero its doc calls impossible | idiomatic | **not applied** — see below |

**`passes: NonZeroU32` was declined, with a reason.** Its consumer is `LocusInference.passes`,
a `u32` whose constructor asserts `passes > 0` for callers other than this loop. Retyping only
this end would leave two adjacent types in the same plan disagreeing and remove no check. The
right moment is C3, where the two meet.

## 6. What the reviews confirmed rather than found

- **Nothing allocates on any path**, measured with a counting global allocator rather than
  argued: 0 allocations for a 2-pass run, 0 for a 20-pass run, 0 for a direct call to the
  stopping rule. The count is independent of the pass count, which is spec §13 test 7's
  property arriving early.
- **The non-convergence path is genuinely non-fatal.** Replacing the cap's `converged: false`
  return with a `panic!` fails three tests, and `LocusInference::new` asserts nothing a
  `converged = false` locus could fail.
- **A non-finite likelihood reaching the loop is not a gap** — a planted `NaN` panics in release
  at the normaliser's total-weight check, naming a `NaN` score.
- **The doc's quotation of `previous_cohort_expected_copies`' own promise is accurate**, checked
  word for word against `src/ng/calling/mod.rs`.

## 7. Out of scope observations

- `&dyn GenotypePriorModel` is re-dispatched once per sample per pass — up to about 150,000
  virtual calls per locus at the committed range. A performance question owned by B1's
  `score_one_sample`, not by C2.
- `RunnableCallingLoopConfig`'s `Deref<Target = CallingLoopConfig>` (`inference/mod.rs`) is the
  Deref-polymorphism shape the checklist warns about; it landed in A2.

## 8. After the fixes

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4620 passed; 0 failed; 14 ignored` (4,603 before C2 — 17 tests added) |
| `cargo test --release --lib ng::calling --all-features` | `574 passed; 0 failed; 3 ignored` (557 before C2) |

**Every mutation the review filed is now caught by the test written for it**, each re-run
singly against the fixed tree: the hoisted coefficient, the dropped `.abs()`, the reordered exit
tests, the weakened length check, the dropped `is_finite()`, the relaxed `<`, and the deleted
sample-count check — one failing test each, except the reordered exit tests which fail one.

**The release-held checks:** C2 now has five, and downgrading all five to `debug_assert`
together under `--release` gives `567 passed; 7 failed`. Each check is reached:

| check | tests that fail without it |
|---|---|
| the previous row names at least one allele | `a_convergence_test_over_no_alleles_is_refused` |
| the two rows are the same length | `cohort_rows_of_different_lengths_are_refused` |
| the cohort holds at least one sample | `a_convergence_test_over_a_cohort_of_no_samples_is_refused` |
| the threshold is a finite positive fraction | `a_threshold_that_is_not_a_fraction_is_refused`, `an_infinite_threshold_is_refused` |
| one inbreeding coefficient per sample | `fewer_…`, `more_inbreeding_coefficients_than_samples_are_refused` |

**One check was deleted rather than tested**, which is the better outcome: the chromosome count
is no longer a `f64` argument, so an infinite or `NaN` one cannot be expressed and the guard
that refused it is gone. What remains is a `Ploidy` — at least one by construction — and a
sample count.
