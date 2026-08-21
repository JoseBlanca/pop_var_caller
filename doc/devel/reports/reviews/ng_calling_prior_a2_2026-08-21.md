# Code Review: ng_calling_prior_a2
**Date:** 2026-08-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A2 of the genotype-prior plan — the local types and the step-8 seam
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step A2 of
  [`calling_prior.md`](../../ng/impl_plan/calling_prior.md), branch `ng-calling-prior`. One file,
  +525/−0.
- **Reviewed against:** base commit `36bcf213` plus `tmp/a2.patch`. Every sub-agent re-pointed its
  own worktree at that commit, applied the patch, and confirmed
  `pub trait GenotypePriorModel` present before reviewing.
- **In-scope files:**
  [src/ng/calling/genotype_prior/mod.rs](../../../../src/ng/calling/genotype_prior/mod.rs)
- **Deliberately out of scope:** the four sibling stub files, filled by later plan steps; frozen
  production; the three aggregate gates already red on `main`.
- **Categories dispatched:** `reliability`, `errors` and `naming` (always-on); `idiomatic`, because
  the diff mints a borrowing newtype, a trait and a lifetime-carrying API; `smells`, because the
  diff carries a lint suppression and a test module that is half its length. `defaults` and
  `module_structure` were not dispatched — A2 adds no default-acting value and no module.

### 2. Verdict

**Request-changes.** No Blocker, six Major, thirteen Minor. The severity is carried by convergence
rather than by any single finding: **three of the five agents independently identified the same
defect and independently compiled the same class of fix**, which is the strongest signal this
process has produced.

### 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.89s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `9 passed; 0 failed` |
| `cargo test --lib` | 0 | `4016 passed; 0 failed; 11 ignored` |

**And one command that should have been run and was not** — see M1. `cargo clippy --lib` does not
type-check `#[cfg(test)] mod tests`, so the verification quoted with the diff could not see five
denied lints in the diff's own test module.

Findings labelled "Needs verification": **0.**

**Mutation totals across the five agents: 25 run, 12 survived, 2 changed no behaviour** — reliability 4/4/0, errors 7/5/1, idiomatic 10/2/0, smells 4/1/1; the naming agent ran renames and reproductions rather than mutations.

### 4. Open questions and assumptions

1. **How far may the implementation reshape a seam the architecture sketches?** Arch §3.2 gives the
   row function six parameters and §7 records "the prior takes flat slices, not the loop's types —
   decided". Three agents argue, and one states it explicitly, that §7's decision is about the
   *provenance* of the type and not a bar on the module defining an aggregate of its own — §7 gives
   two reasons, nothing allocates and no back-reference into the caller, and both hold for a bundle
   owned here. Affects **M2**, **M3**, **M4**.
2. **What size is the per-allele scratch?** Arch §8 lists the scratch slot sizes as an open item
   while the trait fixes one `f64` per allele. Settled at B1 by counting what the ported primitive
   holds. Affects **Mi-scratch**.

### 5. Top 3 priorities

1. **M2/M3/M4 (one defect, three agents)** — the eight-argument seam and the shape checks no
   implementation can be forced to call are the same fact seen from two sides, and one type fixes
   both.
2. **M1** — five clippy errors in the diff's own test module, invisible to the `--lib` command used
   to verify it.
3. **M5** — the module's headline invariant, "held in release, not debug", is exercised by no
   command anyone runs.

### 6. Findings

#### Major

**M1: mod.rs:473–514 — five denied clippy lints in the test module, and the verification command
cannot see them.** *Category: idiomatic.* `cargo clippy --lib` builds the library target only, so
`#[cfg(test)] mod tests` is never checked; CI and `scripts/precommit-check.sh` both run
`--all-targets`. Adding the test target to the reviewed patch gives `clippy::type_complexity` at
473 and `clippy::redundant_locals` four times at 475/488/501/514 — the `let view = view;` lines,
redundant because `GenotypeTableView` is `Copy`. Stated without overclaiming: `--all-targets` is
already red at the base commit on an unrelated example, so this adds five failures to a failing
gate rather than turning CI red.

**M2: mod.rs:265 — the eight-argument seam and the unenforceable shape check are one defect.**
*Category: smells.* `assert_row_shapes` takes six values; `genotype_log_priors` takes the same six
plus `&self` and `inbreeding`. They are one thing — the row a model fills at one locus for one
sample — spread over eight parameters, with the list written twice and prose keeping the copies in
step. The doc's rebuttal ("bundling would move the count without removing anything") holds for a
*passive* bundle and not for a **checked** one: a bundle whose constructor is the check removes the
`#[allow]`, the `pub` helper and the obligation prose at once. Compiled: a borrow-only bundle gives
a three-argument trait method, clippy-clean with no suppression, dyn-compatible, nothing allocated.

**M3: mod.rs:286 — the same conclusion, reached independently.** *Category: reliability.* Built as
a witness type whose only constructor runs the checks; verified dyn-compatible and clippy-clean
with `new` at six parameters, one under the ceiling. "So the survivor the trait doc records as
unclosable is closable; the doc's claim is true of a method body but not of the goal."

**M4: mod.rs:265–275 — the same conclusion, third derivation, more conservative.** *Category:
idiomatic.* Bundle only the three flat table views, whose lengths are fixed together by the table
that produced them; the count drops to six, the `#[allow]` retires, and two of the four checks move
from once-per-sample-per-pass to once-per-locus. Compiled clean with zero suppressions.

**M5: mod.rs:110, 297–312 — "held in release" is the module's headline ruling and nothing anyone
runs can fail on it.** *Categories: errors, reliability (convergent).* Demoting the
empty-concentration `assert!` to `debug_assert!`, or all four `assert_eq!` to `debug_assert_eq!`,
each leaves the suite green; each dies only under `cargo test --release`. CI runs one test command,
in debug (`.github/workflows/ci.yml:47`), and `[profile.release]` does not arm debug assertions.
The tests are already correct release oracles — 7 of the 9 run there — so only the run is missing.

**M6: mod.rs:366 — neither seam test can fail on a wrong genotype order, the one runtime property
they claim.** *Category: reliability.* The biallelic expected row `[0, ln 2, 0]` is a palindrome, so
an implementation walking the coefficients backwards passes; the trait-object test asserts only a
*count* of zeros. Measured: reversing the stand-in's write order leaves both green while the
triallelic row's `1/1` and `0/2` swap places. Nothing anywhere observes the contents of
`genotype_allele_counts` or `homozygous_alleles`.

**M7: mod.rs:295–317 — an empty `out` passes every check, while an empty concentration is
refused.** *Category: errors.* Measured: with `genotype_count == 0` three checks degenerate to
`0 == 0`. A genotype table has at least the all-reference genotype at any ploidy ≥ 1, so a
zero-length row cannot be a thin locus either — the same wiring-bug class the module refuses one
type earlier and argues about at length.

**M8: mod.rs:296 — a mis-sized `out` makes every message blame a correct array.** *Category:
errors.* `genotype_count` is `out.len()`, so `out` is the yardstick and never the subject.
Measured on a well-formed 3-genotype table with `out` sized 6: the panic reads "one log multinomial
coefficient per genotype / left: 3 / right: 6" — the coefficients are correct and `out` is never
named. The module mandates reused caller-owned buffers, so an untrimmed row buffer is the likeliest
mis-shape at this seam and the one the message misdiagnoses.

#### Minor

**Mi1: the parameter `homozygous_alleles` names the allele axis for a slice indexed by genotype**,
and both the architecture and the plan name it `homozygous_allele_for`, which is also the table's
own field name. *naming.*

**Mi2: `Concentration::new`'s doc claims a match with production that does not hold.** Production
checks `α > 0.0`; this checks `α >= MIN_ALT_CONCENTRATION` (`1e-12`), strictly tighter, so `1e-13`
passes production and panics here. The file's own test states the real rule; the type doc and the
test disagree. *naming.*

**Mi3: `[CandidateAlleles]` is a dead intra-doc link** — `cargo doc` reports "no item named
`CandidateAlleles` in scope", in the one sentence that tells the reader what entry 0 means. The
same file qualifies its two other cross-module links correctly. *naming.*

**Mi4: `assert_row_shapes` names the row, and the row is the one thing it does not check** — `out`
and the concentration are the yardsticks and neither is itself checked. *naming.*

**Mi5: `genotype_log_priors` is a noun for a method that writes into a caller's buffer**, where its
own neighbour `assert_row_shapes` is a verb phrase. *naming.*

**Mi6: `SeedRegime`'s three variants are named on three different axes**, and two of them are the
same `(1, θ)` shape differing only in where θ came from — which the names hide and
`FallbackDiversity`'s doc never states. *naming.*

**Mi7: `data_dominated` is a bare participle** — dominated by what, over what? The field's own doc
one line up supplies both. *naming.*

**Mi8: "seed" carries the module's second-most-important idea, is never defined, and already means
three other things in this crate** — a PRNG seed (including inside `ng`), an iterative starting
value the EM then moves (`src/ssr/cohort/em_init.rs`), and initialising a structure. The second is
the trap: this module's own doc says in bold that it fits nothing. *naming.*

**Mi9: "at 2:1 at every realistic diversity" states as exact a ratio both sources hedge.** The
ratio is `2·α_ref : (1 + α_alt)`, so 2:1 only as `α_alt → 0`: 1.998:1 at a human θ, 1.98:1 at 1 in
100. *naming.*

**Mi10: `SpectrumSeed`'s two concentrations are `pub f64` with no constructor**, so `NaN` or a
negative is a legal value and the first check it meets is a `debug_assert` compiled out of release.
*idiomatic.* Two measured details: the `_private: ()` trick is rejected here by
`clippy::manual_non_exhaustive`, and `#[non_exhaustive]` would not block the in-crate caller that
matters.

**Mi11: the panic-message test drops a `&'static str` payload to the empty string and then reports
a false cause.** Measured: a check rewritten as `assert!(cond, "literal")` fires correctly and the
test fails with `said: ` — it cannot tell a wrong message from a non-`String` payload. *errors,
smells, naming (convergent).*

**Mi12: three value checks have no test** — a non-finite entry (dropping `is_finite()` only lets
`+∞` through, which nothing covers), the scratch length's exactness in the *over-long* direction
(the plausible one under a reused-buffer architecture), and an out-of-range homozygous allele id.
*errors, reliability (convergent).*

**Mi13: every test runs at ploidy 2 and none at a single-allele locus**, against a range commitment
that names polyploids and a doc that singles out the homozygous lookup as the thing deferred above
diploidy. Probed: none of the untested corners breaks, so this is a missing guard rather than a
live defect. *reliability.*

**Mi14: `pub fn assert_row_shapes` is `pub` to silence a dead-code warning.** `pub(crate)` plus a
scoped `#[allow(dead_code, reason = …)]` says the same without widening the API; `#[expect]` is the
wrong tool here and was measured to be — the test module calls the function, so the expectation is
unfulfilled under `--lib --tests`. *idiomatic.*

#### Nits

Four dead `let view = view;` rebinds. Two lines of the seed test assert what the derive already
states. `Concentration::get(&self)` and `allele_count(&self)` take `&self` on a `Copy` type where
`src/ng/types.rs` takes `self` throughout. `Concentration` lacks the `PartialEq` its two neighbours
in the same file derive. `{MIN_ALT_CONCENTRATION}` renders as `0.000000000001` where the doc and
spec both write `1e-12`. `genotype_count * allele_count` is an unchecked multiply inside a
release-held assertion, in a profile with `overflow-checks` off. The `#[allow]` carries no
`reason =` where the two nearest precedents in the crate both do. "on the GIAB trio" implies a
joint trio call where the measurement is three samples each called on its own. Test-local bindings
named `scratch`, `panicked`, `expected`, `run`, `fitted`, `copied`, `guessed`, `diploid`. Two test
names lead with "the seam", a label `cargo test` output never defines.

### 7. Out of scope observations

- `src/ng/mod.rs:29` and three sites in frozen production carry unresolved intra-doc links, so
  `cargo doc` is already failing on this branch independently of this diff.
- `examples/ng_cohort_merge_parallel_cost.rs` trips `clippy::manual_is_multiple_of` at the base
  commit, which is why `--all-targets` is red before this diff adds to it.
- Sealing the trait (`mod private { pub trait Sealed {} }`) is the checklist's answer to "the trait
  cannot make it do so", and cannot be written until an implementor type exists at B1.

### 8. Missing tests to add now

Grouped by what they guard, each checked against the mutation it is written for.

1. **The genotype order** — assert the whole triallelic row, which is not a palindrome, plus the
   two view contents the row's values never touch (M6).
2. **The release half** — a second run of this module's tests under `--release` (M5).
3. **An empty row** and **an over-long scratch** (M7, Mi12).
4. **A non-finite concentration entry** at `+∞` (Mi12).
5. **The shapes the caller commits to** — ploidy 1, ploidy 4 and 8, a single-allele locus, and a
   locus at the candidate cap (Mi13).
6. **The message names both buffers** (M8).
7. **`SeedRegime::NeutralShape` constructed at least once**, and `SpectrumSeed`'s refusals (Mi10).

### 9. What's good

- The pointer-identity assertion in `a_concentration_borrows_the_buffer_it_is_given` looks like
  ceremony and is the only thing catching a `get()` that copies — proved by a mutation returning a
  leaked `to_vec()`.
- The `assert_ne!` in the seed test is the only thing catching a hand-written `PartialEq` that
  drops the regime — proved by mutation.
- The four shape cases each differ from a passing call by exactly one truncated slice, and each was
  traced to reach its own check with the earlier ones passing.
- The no-`Result` ruling is honoured cleanly: no `unwrap`, `expect`, `panic!`, `todo!`, discarded
  `Result` or `Drop` in non-test code.
- `Concentration` earns its lifetime independently of the borrow: it is the only place the
  per-allele invariant is named, and without it three `f64` slices with three different meanings
  would sit side by side in one signature.

### 10. Commands to re-verify

- `<worktree>/scripts/dev.sh cargo fmt --check`
- `<worktree>/scripts/dev.sh cargo clippy --lib --tests --all-features -- -D warnings` — **the
  `--tests` is the correction M1 makes; `--lib` alone does not check the test module**
- `<worktree>/scripts/dev.sh cargo test --lib`
- `<worktree>/scripts/dev.sh cargo test --release --lib ng::calling::genotype_prior` — **new, and
  the only command that can fail on the "held in release" half**

Per-category files are left as an audit trail in `tmp/review_2026-08-21_ng-calling-prior-a2/`.
