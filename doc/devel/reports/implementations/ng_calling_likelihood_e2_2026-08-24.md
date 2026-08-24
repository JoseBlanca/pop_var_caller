# ng read likelihoods — E2: two cutoffs, and the mass they discard is now reported

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step E2 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone E,
on top of `44c39a83`.*

## 1. What it is

Two changes to [`alignment/stutter.rs`](../../../../src/ng/alignment/stutter.rs), which the plan
says to keep in their own commit because **an unreported loss compares candidates on different
scales silently**.

**One constant became two.** `MAX_SLIP = 10` applied a single number to the repeat count on the
whole-repeat branch and to the re-indexed base-pair count on the part-repeat branch. It is now
`MAX_WHOLE_REPEAT_SLIP = 10` (repeats) and `MAX_PART_REPEAT_SLIP = 10` (base pairs), each declared
inherited from production's provisional 10 rather than fitted. **Both are still 10, so no score
moves** — what changes is that the two are separately settable by whoever measures them.

**And the model now reports what it cannot place.** `StutterModel::truncated_mass_lost(period,
repeat_count)` returns one minus the total the distribution puts on everything a read of that
candidate could show.

## 2. Why the loss has to be reported, and how big it actually is

A model that quietly loses mass on some candidates and not others is comparing them on different
scales — and genotyping *is* a comparison between candidates. Three things go missing, and they
are not the same size:

| what is missing | how big, measured |
|---|---|
| **the two cutoffs' tails** | `(1 − one_step_share)^10` of a branch's mass. At HipSTR's shipped 0.95 that is `0.05^10`, about 1 part in 10¹³ — negligible. At a one-step share of 0.5 it is about 1 in a thousand of that branch |
| **contractions the tract is too short to reach** | **the term that varies per candidate.** On the shortest tract the copy floors admit — four repeats, at hexamers, slippage level 2 in 100 split 4:1 — **2.0 parts in a million** at a one-step share of 0.95 and **2.0 parts in a thousand** at 0.5 |
| **at period 1, the whole part-repeat branch** | every change is a whole number of repeats when the motif is one base, so that branch is unreachable and *all* its mass goes: **2 in 100** at HipSTR's shipped values |

The third is the largest by a factor of ten thousand and the one most easily mistaken for a defect.
The sweep in `the_reported_loss_equals_one_minus_the_reachable_sum` spans **2.06 in 100 down to
exactly zero** across parameters all inside §8's clamps.

## 3. The one place the specification says two things, and which one this follows

Spec §4.2 states the unreachable slips two ways, and they differ by one step of the geometric:

- **its prose** — "contracting away *more repeats than exist*" — would let a four-repeat tract lose
  all four;
- **its worked figures** — the unreachable tail at a four-repeat tract is "a contraction of **four
  repeats or more**" — would not.

**Only the second reproduces the two sizes the specification itself states**, and it does so
exactly: 2.0 parts in a million and 2.0 parts in a thousand, to the digits quoted. The first
reading gives 1.0 in ten million and 1.0 in a thousand — each one step away.

**So this follows the figures**: a read of a candidate must still show a repeat, so at most
`repeat_count − 1` of them can go. `the_shortest_admissible_tract_loses_the_size_the_specification_states`
is the record of that reading as much as a guard on the arithmetic, and flipping it is one
subtraction plus two numbers in that test.

**This is the deferred question for the owner** (§6): the prose and the figures should agree, and
the spec is not this step's to edit.

## 4. Closed form, checked by enumeration

`truncated_mass_lost` sums five geometric tails, because it is called once per candidate per read
group. The test enumerates instead — it walks every length a read could show, calling
`probability` for each, and requires the two to agree to 1e-12 across **six one-step shares × three
part-repeat one-step shares × six periods × eight repeat counts, 864 combinations**.

Two routes to one number is the point. A version that enumerated inside the model would make spec
§12's fifth test compare a thing with itself.

## 5. Validation

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,369 passed, 0 failed, 14 ignored** (4,358 before) |
| `./scripts/dev.sh cargo test --all-features` — every target | 4,463 / 0 / 18 |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0 |
| `cargo check --examples --all-features` | exit 0 |
| `cargo fmt --check` | exit 0 |

`ng::alignment::stutter::tests` holds **24 tests**, against 19 at E1's close. Five added: the
enumeration check, the specification's two sizes, the period-1 case, the short-tract case, and the
floored-model case.

**Mutation-tested**, each restored from a byte-compared copy:

| mutation | outcome |
|---|---|
| the contraction boundary moved to `repeat_count` | **killed** by 3 tests |
| the part-repeat expansion term dropped from the sum | **killed** by 4 tests |
| the part-repeat regime wired to `MAX_WHOLE_REPEAT_SLIP` | **survives, and must** — the two constants both hold 10, so which one feeds which branch is not observable. Recorded in the boundary test's own doc so the silence is not mistaken for coverage; it becomes checkable the moment a measurement moves either constant |

## 6. Deferred for the owner

1. **Spec §4.2's prose and its figures disagree about whether a tract can lose its last repeat**
   (§3). This step follows the figures. The prose should be corrected to match, or the figures
   recomputed — either way it is a one-line edit to a document this step does not own.
2. **`truncated_mass_lost` has no consumer yet.** It feeds `SsrScoringContext.truncated_mass_lost`,
   which F1 builds. Until then it is exercised only by its own tests.
