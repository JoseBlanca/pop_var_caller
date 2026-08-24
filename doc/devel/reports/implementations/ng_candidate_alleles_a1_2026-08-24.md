# ng candidate alleles — A1: the rule's constants, and a `const` path for the share

*Implementation report, 2026-08-24. Branch `ng-candidate-alleles`, worktree
`../pop_var_caller-candidate-alleles`. Step A1 of
[`candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md), Milestone A, on top of
`3edab4cd`.*

## 1. Plan

Three edits, no logic:

1. `MinAltReadShare` gains `new_const` — a `const` constructor that **panics** where the
   fallible `new` returns `None`, because a `const` declaration has nobody to hand an `Option` to
   and a bad value there should fail the compile.
2. A new module `src/ng/calling/allele_candidates/mod.rs` holding `CandidateSelectionConfig`,
   `DEFAULT_ALLELE_SUPPORT` (floor 2 reads, share 5 in 100) and
   `DEFAULT_MAX_CANDIDATE_ALLELES` (6), each constant carrying its source and its softness.
3. `MinAltReads::reached_by`'s doc comment widened: the numerator is the caller's, and the three
   callers count different reads.

## 2. Assumptions and departures

**The module declaration moved from A2 to A1**, and it is the only departure. The plan puts
`pub mod allele_candidates;` in A2's list. A `.rs` file no `mod` names is not compiled at all —
cargo neither builds it nor warns — so A1's constants and its tests would have been dead
text and the step's own gate (*"the module compiles"*) would have proved nothing. The declaration
is one line, A2 needs it too, and leaving it out would have made A1 unverifiable. Recorded here
rather than escalated: it changes no design and reaches into no other step.

Two smaller things, both inside the step:

- **`CandidateSelectionConfig::DEFAULT` and a `Default` impl were added** beside the two
  constants. The architecture (§2.1) declares the struct and the two constants and does not say
  how a run gets a config holding both. `MinAltReads`, `MinAltObs` and `MinAltReadShare` in the
  merge each carry exactly this pair, so the module follows the neighbour it reuses.
- **`calling/mod.rs`'s own module doc said one of the four sub-modules was present.** Two are
  now; the sentence was updated with the step rather than left to go stale.

**What A1 does not do, deliberately:** nothing reads these constants yet. No run is wired to the
config, and `select_generic` does not exist until Milestone C. The cap's value is declared here
and [`calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2.1 declares the same number for
the loop's post-discovery enforcement; making those one constant rather than two is that plan's
edit, not this one's (arch §4).

## 3. Changes made

- [`src/ng/run/cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs) —
  `MinAltReadShare::new_const`, and `MinAltReads::reached_by`'s widened doc comment. No behaviour
  changed: `new_const` has no caller in the merge, and a doc comment is not code.
- [`src/ng/calling/allele_candidates/mod.rs`](../../../../src/ng/calling/allele_candidates/mod.rs)
  — new. The module doc, `CandidateSelectionConfig`, the two constants, and the tests.
- [`src/ng/calling/mod.rs`](../../../../src/ng/calling/mod.rs) — `pub mod allele_candidates;`
  and the one-sentence doc correction.

**Why the range check is spelled out** rather than asking `(0.0..=1.0).contains(&share)`:
`RangeInclusive::contains` is not a `const fn`. After the review it lives in one private
`const fn is_a_fraction_of_one` that **both** constructors call, so the fallible one and the
`const` one cannot come to disagree about what a legal share is — the first draft wrote the test
out twice and a mutation of one copy left the other's tests green. `share >= 0.0 && share <= 1.0`
refuses everything the range plus an `is_finite()` call refuses: a negative share and `-∞` fail
the lower comparison, a share above one and `+∞` fail the upper, and `NaN` fails both. `-0.0` is
a share of nothing and both accept it. **No `clippy::manual_range_contains` allow is needed** —
clippy already exempts `const` contexts, measured by deleting the attribute (clippy stays clean)
and by adding a non-`const` twin (clippy fires).

## 4. Tests added

Ten, after the review — five in `allele_candidates::tests` and five moved to or added in
`cohort_merge::tests`, beside the constructor they guard.

In `allele_candidates::tests`:

| test | what it pins |
|---|---|
| `the_default_bar_is_two_reads_or_five_in_a_hundred` | the two numbers, **and the floor's coupling to `MinAltObs::DEFAULT`** rather than to the digit 2 |
| `the_floor_decides_at_three_reads_and_the_share_at_three_hundred` | spec §3's two worked examples — 2 at 3 compared reads, 15 at 300 — plus 16 at 301, because `0.05 × 300` is exactly 15 and the first two cannot see the rounding |
| `the_allele_share_binds_only_above_forty_compared_reads` | the claim the constant's doc makes, held against `MinAltReads::DEFAULT`: identical at 1, 3, 11, 20 and **40** compared reads, strictly larger at **41** and at 300 |
| `the_cap_default_is_six_and_the_config_carries_it` | 6, and that `Default` carries it |
| `the_default_config_is_the_two_announced_constants_and_not_the_merges_rule` | that `DEFAULT.support` is the allele rule and **not** the merge's — the assertion whose absence was the review's Blocker |

In `cohort_merge::tests`, beside `MinAltReadShare`:

| test | what it pins |
|---|---|
| `the_const_share_refuses_exactly_what_the_fallible_one_refuses` | ten values through both constructors under `catch_unwind`, requiring the same answer — the claim `is_a_fraction_of_one` exists to make true |
| `a_const_share_below_zero_fails_rather_than_clamping` | the lower bound, which nothing touched before |
| `a_const_share_above_one_fails_rather_than_clamping` | the upper bound |
| `a_const_share_that_is_not_a_number_fails` | `NaN` fails both comparisons |
| `a_const_share_may_sit_on_either_end_of_the_range` | 0.0 and 1.0 as *accepted*, which pins both boundaries as inclusive |

`a_share_outside_a_fraction_of_one_is_refused`, which already existed, gained `f64::NEG_INFINITY`.

**Every one of these was checked by mutation** — see the review report's §5.

## 5. What the review changed

Five agents, each in its own worktree, ran 18 mutations between them; 12 survived the six tests
this step originally shipped. Full account in
[the review report](../reviews/ng_candidate_alleles_a1_2026-08-24.md). In short:

- **A Blocker:** nothing pinned `CandidateSelectionConfig::DEFAULT.support`, so wiring it to the
  merge's rule instead of the allele rule passed every test while halving the share — 2,308
  alternatives kept against 5,596 on the GIAB trio at 300× (spec §3.3).
- **A Major:** `new_const`'s lower bound was untested, and a negative share does not crash — it
  saturates to 0 in `required_of`'s cast and deletes the share half of the bar at every depth. The
  range check is now one private `const fn` both constructors call.
- **A Major deferred to Checkpoint A:** a cap of 0 or 1 is representable and is refusal under
  another name; the fix changes a shape arch §2.1 declares, so what landed is the obligation in
  the field's doc — `select_generic` asserts a cap of at least 2 at step C2.
- **Four measured figures were quoted under conditions they were not measured under**, of 19 the
  naming agent checked. All four rewritten to name their baseline and their depth.
- **The test claimed to be discriminating was half vacuous**: its equality arm stopped at 20
  compared reads, so doubling the share to 0.10 passed it — one read short.
- **`cargo doc` was never in this step's gate**, and the module doc carried a link to a module
  that does not exist yet, against a `deny`-level lint.

## 6. Validation

All in the container, on the tree as committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo doc --lib --no-deps` — 23 unresolved intra-doc links, all pre-existing on `main`, none in
  this step's files (24 before the fix).
- `cargo test --lib allele_candidates` — 5 passed; `cargo test --lib const_share` — 5 passed.
- `cargo test --lib` — see the commit message for the run on the committed tree.

**`cargo clippy --all-targets --all-features` is red on this tree with 14 errors, none in
`src/`** — five benches and examples inherited from `main`, none touched by this plan. This step
is gated on `--lib --tests`, stated rather than quietly narrowed.

## 7. Tradeoffs and follow-ups

- **`new_const` is a second constructor for one type**, which is a cost. The alternative was to
  widen `MinAltReadShare`'s field to `pub(crate)` and let the constant build the struct directly,
  which would have removed the range check from the one place a wrong value is cheapest to catch.
- **Both constants are soft and say so.** The share is measured on one human trio over 572 kb
  (spec Q3) and the cap has never been measured at its own value (spec Q2). Neither has a code
  shape depending on the answer.
- **Nothing here reserves a field for a quality-sum bar** (spec §3.3's third bar). The merge
  already carries `q_sum` per row, so adding one later is a config field and a term in the fold.
