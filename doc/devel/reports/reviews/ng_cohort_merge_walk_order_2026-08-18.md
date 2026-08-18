# Code Review: the cohort merge's "which sample comes next"
**Date:** 2026-08-18
**Reviewer:** rust-code-review skill (orchestrator), two category sub-agents in isolated worktrees
**Scope:** replacing the scan over the cohort with an ordered structure, plus the probe that measured it
**Status:** Approve-with-changes

---

### 1. Scope

- **Reviewed:** the working-tree change as the stash commit `a58c5cfc` over `0166c70f` — at
  review time a `BinaryHeap`; the review's own Major finding replaced it with a tournament.
- **In scope:** `src/ng/run/cohort_merge/close.rs` (the ordering structure, `LocusCloser::over`,
  `take_head`, `sample_with_earliest_head` and the new randomised oracle) and
  `examples/ng_cohort_merge_walk_cost.rs`.
- **Out of scope:** the verdicts and judging in `close.rs`, and the rest of the module.
- **Categories dispatched:** reliability (the change must not move the output), and the hot-path
  items of extras plus smells and idiomatic (the owner asked for the fastest reasonable shape).

### 2. Verdict

**Approve-with-changes.** The output is unchanged and well-guarded; the structure was not the
fastest one available, and both reviewers found gaps in what the tests could see.

### 3. Execution status

- `cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean.
- `cargo test --lib ng::run::cohort_merge` → `165 passed; 0 failed` at review time.
- The probe, before and after, quoted in each finding.
- **Mutation totals:** 18 across the two agents plus the author's five; 3 survived, all three
  proven to change no behaviour.

### 4. Open questions and assumptions

1. **Is the fixture the right one to optimise against?** Every sample carrying a record at
   every position is what the generic mint produces where a sample has reads, but a region only
   part of the cohort covers is the ordinary case at three reads a position. Both were measured;
   the answer differed, and it decided the shape (§5, M1).

### 5. Top 3 priorities

1. **M1** — a tournament tree over the *covering* samples runs the same walk 20–45% faster than
   the heap, and finding the next sample is two thirds of the walk.
2. **M2** — the invariant that keeps the structure and the cursors in step was documented,
   guarded only in debug, and tested nowhere; violating it in release emits one wrong locus.
3. **M3** — the randomised oracle's generator reached none of the shapes that separate an
   ordered structure from a scan.

### 6. Findings

#### Majors

**M1: close.rs — the heap is the right shape but not the fastest one available.**
**Category:** performance. **Confidence:** High, measured.
The heap pays three times over: a 24-byte key where 16 will do, a pop followed by a push where
one replace-top would do, and a sift-down comparing *both* children a level. The merge has a
property none of that uses — **k is fixed for the whole walk** — so a tournament tree can be
built once and only have its values changed. Measured end to end through `build_region`, two
runs each: 63 samples 18.3/17.5 µs against 13.8/14.7; 250, 116/131 against 75.8/73.3; 1,000,
723/728 against 422/384; 3,000, 3342/3412 against 2222/2226. **The leaves must be the covering
samples**: over the whole cohort the tournament is slower than the heap where little of the
cohort covers the region — 74.4 µs against 37.7 at 3,000 samples and 1% covering — and keyed to
the covering samples it is the fastest candidate there, 27.7 µs.
*Also measured and rejected:* a 4-ary heap (worse at every k from 10 up), `peek_mut` alone (no
gain at 3,000), a small-k threshold keeping the scan (the scan never wins at any k — at 10
samples it is already twice the heap).

**M2: close.rs — the invariant that keeps the structure and the cursors in step is guarded only
in debug.** **Category:** reliability. **Confidence:** High, measured in release.
Consuming a sample other than the one the structure showed emits one wrong locus — holding one
sample's record while another's two observations vanish from the run — before the walk dies on
a key with no head. Loud rather than silent, but the doc claimed "silently" and no test drove
it.

**M3: close.rs — the randomised oracle's generator misses four of the six shapes it was asked
for.** **Category:** reliability. **Confidence:** High, measured by instrumenting the generator.
Over its 300 seeds: 0 cohorts with an empty sample, 0 where one sample's observations wholly
precede another's, none with a single observation, and a cohort size topping out at 8 against a
caller committed to thousands. Ties (267 of 300) and two contigs (150 of 300) were reached.
**No mutant escaped through the gap** — the finding is about what the check will cover for the
next change.

#### Minors

**Mi1 — the tie-break the doc fixes is pinned by no test.** Flipping it leaves every other test
green, which is precisely why the rule needs one: it is a promise to a later step.

**Mi2 — the performance property has no regression guard.** The change is
behaviour-preserving, so the whole suite passes either way; a scan-based rewrite would ship
green. A structural assertion, not a timing one, is the fix.

**Mi3 — the key is a tuple whose second field carries meaning only in prose**, and a third of
its 24 bytes is padding. Narrowing it alone took 1,000 samples from 913/884 µs to 645/675.

**Mi4 — the probe reports one mean and no spread**, on a machine whose own swing at 3,000
samples reached 30% between runs of the same binary; and its second column measured within that
swing of its first at every cohort size, so it implied a difference the measurement does not
carry.

**Mi5 — a wrong ratio in the doc comment.** "56.9 µs at 63 samples and 101 ms at 3,000 — four
times the cohort for fourteen times the time" does not follow from its own two numbers: 63 to
3,000 is 47.6 times the cohort and 1,783 times the time. The 4×/14× is the step between
adjacent rows.

#### Nits

The probe pins its substitution by two bare indices that must agree; `take_head`'s destructure
binds a value only a `debug_assert!` reads.

### 7. Out of scope observations

- The walk body is ~120 lines with three inline assertions and their justifications; whether
  that is one function or three is a structure question, not a performance one — the ordering,
  not the body, is what dominates.
- The coordinate-order assertion inside the walk compares positions only, and the loop breaks on
  a contig change before reaching it, so a sample descending by *contig* passes unchallenged.
  Unchanged by this work; the scan behaved the same.

### 8. Missing tests to add now

All added: the tie-break; the structural guard (one live leaf per unspent sample, one
consumption per observation); a sample wholly before another; and the identity check that every
sample is handed back its own observations, once each, in order, inside its locus's ground.

### 8a. The diff's own quantitative claims

| claim | verdict |
|---|---|
| the four timings in the doc comment | CHECKED-CORRECT against the probe |
| "four times the cohort for fourteen times the time" beside 63 and 3,000 | **WRONG** — that is the step between adjacent rows, not those two |
| the probe attributes its figure to the walk | **CHECKED-CORRECT** — measured separately, the walk is 97% of it at 3,000 samples and 90% at 63 |
| the two columns' difference | **not supported** — below the machine's own run-to-run swing at every cohort size |

### 9. What's good

- **The ordering is two thirds of the walk**, established by an oracle that already knows the
  merge order and so pays nothing to find the next sample — which is what makes the structure
  worth arguing about at all.
- **Every candidate drove the identical walk in one binary in one run**, with the loci asserted
  equal before anything was timed.
- **Construction was priced rather than assumed**: 1.1% of the walk at 63 samples, 0.35% at
  3,000, so making the closer reusable would buy at most that.
- **The desync class is structurally loud**: a wrongly consumed sample accumulates one entry
  more than it has observations left, so a key with no head is always reached.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh cargo run --release --example ng_cohort_merge_walk_cost
./scripts/dev.sh bash tmp/mutate_heap.sh
```
