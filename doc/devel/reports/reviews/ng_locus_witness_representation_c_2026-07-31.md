# Code Review: ng locus witness representation — Milestone C (the fold)

**Date:** 2026-07-31
**Reviewers:** three category sub-agents, one isolated worktree each, plus an
orchestrator-run documentation pass
**Scope:** `ebe3685..82b13a0` on `ng-pileup-generator` — C1, C0, C2, C3
**Status:** Approve-with-changes — all applied, in `fc7d839`, `3fecdf6`, `4e0e02f`, `19e38f8`

---

## 1. Scope

- **What:** the four commits that make the fold speak in witnessed *sets* and stop discarding
  a read with a hole in its witness — `ebe3685` (C1, the fold hands over
  `WitnessedRefPositions`, still discarding a hole), `6805e42` (C0, added mid-milestone by
  owner decision: a read that never enters the tract is not an STR partial), `761d53e` (C2,
  `ReadWitness::Partial` carries a set; absorbed plan steps C4, D1, D2), `82b13a0` (C3, a
  holed witness is recorded).
- **In scope:** `src/ng/locus_generation/**`, `examples/ng_*_loci_dump.rs`,
  `examples/ng_ssr_aligner_bakeoff.rs`, and the milestone's docs.
- **Out of scope:** `src/pileup/` (production, frozen); Milestone D.
- **Categories dispatched:** `reliability` (mutation sweep — the invariants are the
  deliverable), `behaviour_safety` (did each step change only what it claimed?),
  `module_structure` + API-fitness (can the API carry Milestone D?).
- **Run by the orchestrator instead of an agent:** `naming` + documentation accuracy. Recorded
  as a departure, not a gap — see §3.

## 2. Verdict

**Approve-with-changes.** The milestone's headline claims survived independent checking, and
one came back stronger than it was stated: an eleven-fixture differential rendered through the
dump's own `report.render()` at five commits shows C1 and C2 **byte-identical**, C3 adding
holed rows and **modifying or deleting none**, with no `chain_ids` field on any surviving row
moving. `git bisect` has a clean single culprit.

What the review found instead: **five guards that survived their own mutation**, two of them
guarding this milestone's own decisions, and **two shapes that will drift silently** — a
hand-written sum over a struct's fields, and one fact stored in two places.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` | clean, at review time and after fixes |
| `cargo test --lib --bins --tests --examples --all-features` | 2,831 → **2,835** passed, 0 failed |
| `ng::locus_generation` | 300 → **304** passed |
| STR dump on the tomato CRAM | byte-identical to the C0 baseline throughout |

**A departure from the review protocol, recorded.** The four category agents were launched in
parallel and **all four stalled at their first command** — each tried to start a 16 GB
container VM and only one can run at a time on this machine. They produced nothing and were
re-run **sequentially**. The fourth category, naming and documentation accuracy, was then run
by the orchestrator by reading rather than by an agent, because it needs no container and the
sequential queue was the bottleneck. That is weaker than the other three: it is the one
category here whose findings were not adversarially produced. It found ten stale statements
(`fc7d839`) and should be re-run as an agent if the milestone is revisited.

## 4. Top 3 priorities (all applied)

1. **C3's one design decision was pinned by nothing** — reverting `read_agreed_with_reference`
   to its pre-C3 body left the whole suite green, including both fixtures that produce a holed
   witness.
2. **A four-term sum over `SsrGeneratorCounts`' fields, in two tools**, one of which had
   already been left behind once. Adding a fifth reason and changing nothing else compiled
   clean with both tools under-reporting.
3. **`apply_events_into` answered a `bool` *and* filled a buffer** — one fact twice, with the
   `bool` ignorable, which is the hazard `canonicalise_runs` had been reshaped to close one
   milestone earlier.

## 5. Findings

### Major

**M1: `open_record.rs` — C3's `read_agreed_with_reference` rule is unpinned**
- **Category:** reliability
- **Confidence:** High. Putting the pre-C3 enclosing-extent body back verbatim left
  `300 passed; 0 failed` and the full suite green. The two walk fixtures that produce a holed
  witness cannot discriminate: their concatenated bases are *shorter* than the enclosing
  slice, so both readings answer `false` for different reasons.
- **Applied:** `a_witness_with_a_hole_never_counts_as_agreeing_with_the_reference` — a read
  witnessing 100 and 102 of a three-base record, in a bucket holding the whole reference
  slice, so the two readings disagree. The old one withholds a chain id on the strength of the
  base at 101 the read never saw. Mutation now fails.

**M2: `ssr.rs` — C0's guard had become unfalsifiable, and C2 is what did it**
- **Category:** reliability
- **Confidence:** High. Guarding only the left border *and* deleting the guard outright both
  left everything green. C2 gave `from_left`/`from_right` an `Option`, so `partial()`'s
  handling of `None` already answered `OutsideTract` and no input could tell the paths apart.
  **The C0 commit's quoted mutation output no longer reproduced at HEAD.**
- **Applied:** the guard removed — one decision, one place. The test's claim that "a fix
  catching only one border would halve the population" was replaced by what is now true, and
  the discriminating mutation moved onto the constructor.

**M3: `SsrGeneratorCounts` — a hand-written reason sum in two tools**
- **Categories:** reliability (found the instance), module_structure (found the shape)
- **Confidence:** High. `ng_ssr_aligner_bakeoff` summed three of four, under-reporting by
  6,704 reads of ~9,265 on tomato chr01. Adding a *fifth* reason and changing nothing else
  then left `clippy … -D warnings` clean with **both** tools short.
- **Applied:** `SsrGeneratorCounts::reads_without_observation()`, an exhaustive destructure
  with no `..`. The same mutation now fails with `error[E0027]: pattern does not mention field`.

**M4: `open_record.rs` — `apply_events_into` stored one fact twice**
- **Category:** module_structure
- **Confidence:** High. Three sites reconciled the `bool` against the buffer by hand, and the
  `bool` had no `#[must_use]`.
- **Applied:** it returns nothing; the buffer is the only answer and both callers go through
  `take_from`/`refill_from`. Two `expect`s and a `debug_assert` gone; in `refold_live_reads`
  the refill moved ahead of the bucket work, so witness and `allele_index` have no window in
  which to disagree. **No mutation is quoted** — the hazard is removed rather than detected,
  and failing to consume the buffer is harmless because every call clears on entry.

### Minor (all applied)

- **`finalise`'s cap exclusion had become a guard that cannot fail** — C3 removed the state
  that made it fire. Now reached directly by a hand-built record with one read id in both
  lists.
- **`num_obs_along_locus`'s clamp**, which its own comment calls *the* guard, survived
  deletion: the only over-long run in a test came through a constructor that had already
  clamped it.
- **`is_flush_right`'s `>=` survived becoming `==`** — the over-long run its doc justifies was
  never built.
- **`refill_from` / `take_from` documented "the buffer untouched" on failure** while
  `mem::take` had emptied it, and a comment promised to "give it back". The emptied behaviour
  is the useful one; the docs now say so and a test asserts it.
- **`sort_key` published the encoding** the type documents as private behind `runs()`. It now
  returns the set, with `Ord` derived — sound for the same reason `Eq` is.
- **`witness_of` narrowed run coordinates through `LocusLen::from_positions`**, using the type
  as a saturating cast when its own doc says it exists to keep two same-shaped `u16`s apart,
  and clamping where spec §3.4 asks for an error.
- **`witness_of`'s panic message named only the invariant**, and an inverted footprint died
  inside `u32::clamp` because the width guard's `saturating_sub` reads 0 on an inversion.
- **`WitnessedLocusPositions` guaranteed non-empty and no accessor said so** — consumers paid
  an `expect` per border or hid it behind an `is_some_and` answering `false` for the impossible
  set. `first_run` / `last_run` / `span` are total.
- **The removed `tract.is_empty()` guard also covered an *inverted* span**, since
  `Range::is_empty` is `start >= end`; `partial()` underflows there. A `debug_assert` makes the
  failure legible.
- **Both `expect(dead_code)` attributes became `#[cfg(test)]`** — the `expect` fired only in
  the clippy step, since `cargo test --lib` cfg-es it away.
- **Ten stale doc statements** describing the shape the code had before the step that changed
  it (`fc7d839`), plus a wrong number in a C3 assertion message and an accounting assertion
  that cannot fail, relabelled to say what it proves.

### Recorded, not applied

- **`reads_without_observation` looks structurally unreachable on the generic path** after C3.
  Left standing: spec §6 owns the question and its answer needs the comparison against
  `reads_silent_over_footprint` that the spec asks for. The question is sharpened there,
  including that the STR path uses the same counter for four reasons.
- **Flush at both borders is not the same as pinned** (F10) — a hole in the middle, or an STR
  reach measured in read bases exceeding the reference tract. Only `Complete` says pinned.
  Recorded on `ReadWitness`, because it is also the real reason the constructors do not
  implement arch §1.1's `Complete` short-circuit: the implementation report justified that on
  the STR dump moving, and the stronger argument is correctness — `Complete` gates what a
  likelihood may score as an exact length. **For the owner at D3.**
- **Four `witness_label`s, three spellings**, two already inconsistent inside one function
  (`partial_left`/`partial_right` beside `partial:interior`). **D4 owns this** by the plan,
  which asks it to share the derivation and let each tool spell its own strings.

## 6. What's good

- **The differential was built rather than borrowed.** A reviewer refused to trust edited test
  expectations and rendered eleven fixtures through the tool's own output path at five
  commits. That is what makes "C3 is the only commit that moves generic output" a measurement.
- **The `sort_key` order claim was proven, not assumed** — algebraically and then exhaustively
  over a 132×132 grid.
- **The type boundary is real, compiler-verified.** Three bypass attempts, three errors;
  `WitnessedRefPositions` is tighter than arch §2 promises because the module itself is
  private to `pileup`.
- **Milestone D's call sites compile against the current types** — `from_run`, both dumps, and
  the census's per-run positions — which is the check that found the previous milestone's top
  defect.
- **C0's correction came from the owner questioning a number**, not from the code. "I don't
  believe half the reads are precisely anchored on the flank but cover 0 bases of the tract"
  was right, and the measurement that followed found a class of read the STR path should never
  have been minting.

## 7. Commands to re-verify

```
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --lib --bins --tests --examples --all-features
```

Each fix must fail under its own mutation, not merely pass:

| test | mutation |
|---|---|
| `a_witness_with_a_hole_never_counts_as_agreeing_with_the_reference` | revert `read_agreed_with_reference` to the enclosing-extent read |
| `a_read_covering_only_a_flank_is_outside_the_tract` | `from_left`: `covered.max(1)` |
| `depth_clamps_a_witness_that_reaches_past_the_locus_it_is_attached_to` | delete both `.min(len)` in `num_obs_along_locus` |
| `a_run_reaching_past_the_locus_is_flush_right` | `is_flush_right`: `>=` → `==` |
| `a_read_counted_out_and_also_named_to_the_cap_is_not_reported_as_capped` | delete the `!reads_without_observation.contains` clause |
| `taking_and_refilling_leave_the_callers_buffer_empty_and_reusable` | restore the buffer on a failed `refill_from` |
| *(compile-time)* `SsrGeneratorCounts::reads_without_observation` | add a fifth field to the struct |
