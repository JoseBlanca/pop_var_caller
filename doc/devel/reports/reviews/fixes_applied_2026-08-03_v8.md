# Fixes applied — ng read filtering in stages, C2

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `7e8cfce`
**Review:** [`ng_read_filtering_stages_c2_2026-08-03.md`](ng_read_filtering_stages_c2_2026-08-03.md)
**Impl report:** [`ng_read_filtering_stages_c2_2026-08-03.md`](../implementations/ng_read_filtering_stages_c2_2026-08-03.md)

Applied: 13 · Applied with adaptation: 2 · Already fixed: 0 · **Deferred: 3** · Disputed: 0

---

## 1. Findings table

| id | severity | subject | decision | status | validated |
|---|---|---|---|---|---|
| B1 | Blocker | a reposition may discard the tally | Apply | Applied | Pass — new test kills the mutation, alone |
| M1 | Major | the widening landed half-pinned; its comment was false | Apply | **Applied with adaptation** | Pass — §2 |
| M2 | Major | the checkpoint condition is not met | Apply | Applied (owner-ruled) | Pass — cycle closed, 961 lines |
| M3 | Major | `counts` / `read_group_counts` collision | Apply | Applied | field renamed `read_group_tally` |
| Mi1 | Minor | row 5 of the accounting overstates | Apply | Applied | comment rewritten |
| Mi2 | Minor | fusing no longer latches | Defer | Deferred | §3 |
| Mi3 | Minor | one-buffer-per-pass untested | Defer | Deferred | §3 |
| Mi4 | Minor | `"the walk stays stopped"` cannot fail | Defer | Deferred | §3 |
| Mi5 | Minor | three visibilities too wide | Apply | Applied | `pub(in crate::ng::read)`, clippy clean |
| Mi6 | Minor | `buffer` / "buffered read" vocabulary | Apply | Applied | doc states the relation |
| Mi7 | Minor | `verdict_on_raw_read` takes no raw read | Defer | Deferred | §3 |
| Mi8 | Minor | stale prose in three places | Apply | Applied | rewritten |
| nits | Nit | accounting comment, `decode_fails` doc | Apply | Applied | — |
| — | — | the timing challenge | Apply | **Applied with adaptation** | re-measured, six runs a side |

## 2. The two applied with adaptation

### M1 — fixed in the code, not the prose

The review's minimum was to correct the false comment. The code was fixed instead, because the
*claim* was the one worth keeping: the reposition now precedes **everything**, including the
eviction and the counters, so a failed jump leaves the kept set, the counters and the region state
exactly as they were. The promise is now true rather than narrowed.

That also made the ordering *testable*, which it was not: the flag masks it, so the ordering had no
observable consequence on its own. The prerequisite was making the scripted seek fault
**positional** — `with_failing_seek` → `with_failing_seek_at(n)` — because an all-or-nothing fault
cannot express "a seek that fails after a region has been served", and that is the only state in
which "left exactly as it was" differs from "was never anywhere".

### The timing challenge — re-measured rather than argued

The review refused a one-run-each-side comparison, citing B1's six-runs-a-side bar on the same
probe. Re-measured by rebuilding `6e22718` in a throwaway worktree:

| | six runs | mean |
|---|---|---|
| baseline `6e22718` | 1.820 1.821 1.825 1.829 1.830 1.844 | **1.828 s** |
| C2 | 1.824 1.824 1.829 1.830 1.834 1.834 | **1.829 s** |

Under 0.1 %. The single 1.782 s reading the first draft used as the baseline was an outlier. The
impl report's §6 is corrected.

## 3. Deferred, each with a reason

- **Mi2 — the fusing no longer latches.** Behaviour-preserving today, and changing it is a
  *contract* decision: either the cursor re-latches a clean end of input, or the readers' "latch
  `done`" behaviour becomes a written part of the `AlignedReadsReader` contract. Both are design
  moves, and the second is where it probably belongs. → recorded in the impl report; raise at
  Checkpoint C if the owner wants it settled before D.
- **Mi3 — the one-buffer-per-pass invariant is untested.** The same class C1b closed for
  `fill_raw_read`, one layer up, and the same fix shape (capacity assertions). Left because C3
  rewrites this seam when `RecordSource` goes, and the test wants writing against the final
  signature rather than twice.
- **Mi4 — `"the walk stays stopped"` cannot fail in the reference-failure test.** Real, and its
  sibling `Source` test does kill the fuse, so the property is covered. Fixing it well means a
  fixture whose fatal record is *not* last, which is a fixture change in a test the review
  otherwise endorsed. → folded into C3's pass over these tests.
- **Mi7 — `verdict_on_raw_read` takes `(flag, mapq)`, not a raw read.** Arch §7 left this open
  "for Milestone C", so it is legitimately in scope — but it is a signature change to a function
  C3 and C4 do not touch, and taking it here would put an unrelated API change into the step whose
  failure is silent. → **Checkpoint C**, with the reviewer's note that the whole-read form compiles
  and passes.

## 4. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,857 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**`cargo doc` reached 13 twice during this work and the count caught both** — a bad
`AlignmentCursor` intra-doc path, and `read/mod.rs` still pointing at `filtering::ReadFilterError`
after the type moved.

### Mutations re-run after the fixes

| mutation | before | after |
|---|---|---|
| clear the tally on the jump path | **survived 2,855** | killed by `the_tally_survives_a_reposition_that_drops_everything_and_jumps`, alone |
| revert the commit ordering | **survived 2,855** | killed by `a_failed_reposition_leaves_the_cursor_untouched`, alone |
| fold every read group into one entry | killed | still killed |
| `other_sample` rider on the last entry | killed | still killed |
| a failed reposition does not set `failed` | killed | still killed |

## 5. Disputed findings

None. The review contradicted the author on two points — the "out of scope" call for the re-homing,
and the timing evidence — and was right on both.

## 6. Failed-validation findings

None. One implementation misstep is recorded rather than hidden: the first `awk` cut that removed
the re-homed blocks also took the `RecordSource` trait, which the compiler caught immediately. Line
-range surgery on a file `cargo fmt` has reflowed is the hazard the plan's own notes name.
