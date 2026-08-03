# Fixes applied — ng read filtering in stages, C1

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `6e22718`
**Review:** [`ng_read_filtering_stages_c1_2026-08-03.md`](ng_read_filtering_stages_c1_2026-08-03.md)
**Impl report:** [`ng_read_filtering_stages_c1_2026-08-03.md`](../implementations/ng_read_filtering_stages_c1_2026-08-03.md)

Applied: 9 · Applied with adaptation: 2 · Already fixed: 0 · **Deferred: 2** · Disputed: 0

---

## 1. Findings table

| id | severity | subject | decision | status | files | validated |
|---|---|---|---|---|---|---|
| M1 | Major | a failed reposition can be swallowed | Apply (part) + Defer (part) | Applied with adaptation | `in_memory.rs`, `cursor.rs` | Pass — kills the surviving mutation |
| M2 | Major | `ReadFilterError::Source` names two conditions | Defer | Deferred | none | N/A |
| M3 | Major | the recorded justification is factually wrong | Apply | Applied | `in_memory.rs`, `cursor.rs`, `filtering.rs` | Pass |
| M4 | Minor | `with_failure_at_read` no-ops past the end | Apply | Applied with adaptation | `in_memory.rs` | Pass |
| M5 | Minor | the reference test near-duplicates its neighbour | Apply | Applied | `cursor.rs` | Pass |
| M6 | Minor | a strictly dominated test | Apply | Applied | `in_memory.rs` | Pass |
| M7 | Minor | `script_position` collides with ng's `Position` | Apply | Applied | `in_memory.rs`, `cursor.rs` | Pass |
| M8 | Minor | `read_error_at` is a prepositional fragment | Apply | Applied | `in_memory.rs` | Pass |
| M9 | Minor | contradictory ordinals | Apply | Applied | `cursor.rs`, `filtering.rs` | Pass |
| M10 | Minor | `step_one_failure` overclaims | Apply | Applied | `cursor.rs` | Pass |
| M11 | Minor | the decode note omits C2/C3's cliff | Apply | Applied | `filtering.rs` | Pass |
| nits | Nit | stale `expect`s, domain-term collision, test names | Apply | Applied | 4 files | Pass |

## 2. The two that were not simply applied

### M1 — applied in half, deferred in half, on the owner's ruling

**What was applied.** `InMemoryAlignedReadsReader::with_failing_seek()` — the second way a reader
can break — plus `a_reposition_that_fails_is_refused_rather_than_answered`. Together these kill
the surviving mutation: replacing the jump arm with `let _ = …jump_to(region); Ok(())` now fails
that test and only that test.

**What was deferred, and why it is not a punt.** Writing the test showed the defect is larger than
a missing assertion. `move_to_region` commits `region` and `last_region_start` *before* the
fallible `jump_to`, so a failed reposition leaves the cursor pointed at a region it never reached,
serving from wherever the seek abandoned the reader — and the next forward region then takes the
*reuse* path and reads on without jumping. The repair is the cursor's own `failed` flag covering
both routes into a stopped cursor, plus committing the state only after the jump succeeds.

That flag is **C2's named deliverable**. Building it in C1 would mean C2 merging with a half-built
version of its own central change. **Owner's ruling (2026-08-03): widen C2 to cover it.** Recorded
at the test, in C1's impl report, and carried into C2's brief.

The tree is green and honest: the assertion C1 makes is true, and the two it does not make are
named at the test rather than left to be discovered.

### M4 — applied with a different fix than the review proposed

The review proposed `is_some_and(|at| self.next_index >= at)`. **That does not work, and running
it is how we know:** `next_index` stops advancing once the script is exhausted, so it never
*reaches* a fault scripted beyond the last record. The first attempt at this fix left
`with_failure_at_read(5)` on a one-record script silently inert — the review's own diagnosis, one
step further along.

The fix that holds clamps to the end of the script:

```rust
.is_some_and(|failing| self.next_index >= failing.min(self.records.len()))
```

so every past-the-end fault fires where the file breaks. Pinned by three cases in one test — at
the end, well past it, and an empty script.

## 3. Deferred, with a home

- **M2 — `ReadFilterError::Source` names two unrelated fatal conditions.** Adding a variant is a
  design change, and arch §4 forbids it in the same sentence that mis-describes the enum. Both
  the code and the doc need a ruling. → **Checkpoint C**, on the owner's decision.
- **M1's second half** — the cursor must stop after a failed reposition. → **C2**, widened.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt` | clean |
| `cargo test --lib` | **2,847 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**Suite 2,842 → 2,847 (+5), fully accounted:** `in_memory.rs` +5 (four scripted-fault tests, one
seek-failure test, minus the dominated one), `cursor.rs` +3 (two re-pointed, one reposition),
`filtering.rs` −3.

### Mutations re-run after the fixes

| mutation | before the fixes | after |
|---|---|---|
| failed `jump_to` swallowed, `Ok(())` returned | **SURVIVED (2,845 / 0)** | killed by `a_reposition_that_fails_is_refused_rather_than_answered`, alone |
| fault scripted well past the end never fires | survived (silently inert) | killed by `a_fault_scripted_at_or_past_the_end_…` |
| fault fires on every read | killed | still killed |
| script consulted before the fault | killed | still killed |
| the fault is consumed | killed | still killed |
| source failure charged to `Decode` | killed by the new test alone | still killed |
| reference failure charged to `Source` | killed by the new test alone | still killed |
| the narrowing swallows the reader's failure | killed (new test + `t10`) | still killed |

## 5. Disputed findings

None. Two of the review's own briefing errors were corrected by its agents and are recorded in
the review report rather than disputed here.

## 6. Failed-validation findings

None — but M4's first attempt failed its own new test and was replaced, which is recorded above
rather than hidden: the review's suggested code was wrong in detail, and only running it said so.
