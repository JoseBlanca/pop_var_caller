# Code review — ng read filtering in stages, C1

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `6e22718`
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **C1**.
**Impl report:** [`ng_read_filtering_stages_c1_2026-08-03.md`](../implementations/ng_read_filtering_stages_c1_2026-08-03.md).

---

## 1. Scope

The uncommitted working-tree diff for C1 — 3 files, 302 insertions / 64 deletions:

- `src/ng/read/input/aligned_reads_reader/in_memory.rs`
- `src/ng/read/input/cursor.rs`
- `src/ng/read/filtering.rs`

**Out of scope:** what C2/C3/C4 will change; production code outside `src/ng/`.

**Categories dispatched**, one `general-purpose` agent each, each in its own git worktree:

| category | why |
|---|---|
| `reliability` | always; and C1 is almost entirely tests |
| `refactor_safety` | C1 deletes three tests and adds six — the trade needed auditing |
| `errors` | C1's whole subject is the three fatal conditions |
| `naming` | C1 adds new API into the vocabulary Milestone A had just settled |

Because the change was uncommitted, each agent detached at `6e22718` and applied an exported
patch, rather than checking out a branch.

## 2. Verdict

**Approve with changes** — all applied; see the fixes report.

Nothing was wrong with what C1 *does*. Every finding was about what C1 *claims*, plus one live
defect the review's mutation pass uncovered in code C1 does not touch.

## 3. Execution status

Every agent independently reproduced the author's gate in its own worktree: `cargo fmt --check`
exit 0, `cargo test --lib` **2,845 passed / 0 failed**, `cargo clippy --all-targets
--all-features -- -D warnings` exit 0.

**Mutation coverage was the point of the fan-out, and it was used.** Across four agents,
**33 mutations** were applied, `grep -c`-confirmed present, run against the full suite and
reverted. Thirty-two were killed. **One survived.**

## 4. Open questions and assumptions

Both went to the owner and were ruled on the same day (2026-08-03).

1. **The surviving mutation's repair (M2 below) belongs to C2**, widened to cover both routes
   into a stopped cursor rather than only the one that replaces `FilterState`. C1 pins the
   refusal; C2 pins the stop.
2. **The two design-level Majors are carried to Checkpoint C** — `ReadFilterError::Source`
   conflating two fatal conditions, and the plan's refuted D2.

## 5. Top 3 priorities

1. **M1 — a surviving mutation on a fatal path.** `move_to_region`'s failed reposition could be
   swallowed entirely with the whole suite green.
2. **M3 — the diff's own justification was factually wrong**, in three comment blocks, in the
   direction that flatters the change.
3. **M4 — `with_failure_at_read` silently no-ops** for a fault scripted past the end of the
   script, while its doc claims the opposite.

## 6. Findings

### Major

#### M1: `src/ng/read/input/cursor.rs:474` — a failed reposition can be swallowed, and the whole suite stays green
**Categories:** reliability, errors (convergent — found independently by two agents)
**Confidence:** High (mutation-verified)

Replacing the jump arm with `let _ = self.filter.source_mut().jump_to(region); Ok(())` passes
**2,845 / 0**. Every other test touching that line asserts the move *succeeds*.

Writing the missing test then exposed more than an untested path. `move_to_region` sets
`self.region` and `self.last_region_start` **before** the fallible `jump_to`, so after a failure
the cursor is pointed at a region it never reached with the reader left wherever the failed seek
abandoned it — and `next_read` serves from there. Because `last_region_start` moved too, the
*next* forward region takes the reuse path and carries on reading without jumping at all. That is
exactly the "plausible, silently short answer" `CursorError::AfterFailure` exists to prevent; the
guard misses it because it asks `self.filter.has_failed()`, and a failed reposition never reaches
the filter.

**Latent, not live:** both real arms' `begin_region` are effectively infallible today (the BAM arm
queries an in-memory index, the CRAM arm resets state).

**Resolution:** C1 adds the fault-injection knob (`with_failing_seek`) and pins the refusal, which
kills the mutation. The stop is C2's, by the owner's ruling.

#### M2: `src/ng/read/filtering.rs:477` — `ReadFilterError::Source` names two unrelated fatal conditions
**Category:** errors · **Confidence:** High

`resolve_read_group`'s failure leaves `RegionRawAlignedReads::read_next` through the same `?` as
an I/O read failure, so an unresolvable `RG` tag renders as *"reading the next alignment record
failed"*. Rendered from a real BAM:

```
[0] reading alignment file '…/mixed.bam' failed
[1] reading the next alignment record failed
[2] read 'untagged' the record carries no RG tag, and its file declares several read groups
```

So there are **four** fatal conditions, not three, and **arch §4's table is wrong** where it says
the three variants name the three pieces. "Re-fetch the BAM" and "fix the `@RG` header" are
different operator actions.

**Resolution:** deferred to **Checkpoint C** — adding a variant is a design change, and arch §4
forbids it in the same breath that it mis-describes the enum.

#### M3: three comment blocks — the recorded justification for deleting the tests is factually wrong
**Categories:** reliability, refactor_safety, naming (convergent — three agents)
**Confidence:** High (mutation-verified in both directions)

All three blocks claimed these properties had been pinned only by doubles bypassing two layers.
Two pre-existing tests say otherwise, both present at `6e22718`:

- `open_bam.rs::t10_a_truncated_file_fails_once_and_then_refuses_later_regions` truncates a real
  indexed BAM mid-walk and asserts reads flow, exactly one error, the walk stops, and
  `AfterFailure`. It dies under "the narrowing swallows the reader's failure".
- `cursor.rs::a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short`
  already drove the reference failure through the whole chain, using the **identical** `overruns`
  fixture the new test uses, thirty lines above it.

**What C1 genuinely buys is narrower and was not stated:** no test anywhere pinned *which*
`ReadFilterError` a fault is charged to, because both of those match `Err(_)`. Charging a read
failure to `Decode` kills exactly one test in the tree — the new one; charging a reference failure
to `Source` kills exactly one — the other new one.

On a branch whose named recurring defect is a test claiming more than it proves, a *comment*
claiming a coverage gap that did not exist is the same failure of the record — and it is
load-bearing, because it is the argument the plan uses to order C1 before C3.

### Minor

- **M4 — `with_failure_at_read` silently no-ops past the end of the script** (errors). The doc
  said any past-the-end fault is the truncated file; the code compared for equality, and
  `next_index` stops advancing at the end, so only `== len` fired. A fault-injection knob that
  accepts an unreachable position in silence is the failure mode this plan is written against.
- **M5 — the new reference test near-duplicates its neighbour** (reliability). Same fixture, same
  region, same `AfterFailure`; one differing assertion. That is how the wrong one gets deleted
  later as "the duplicate", taking the assertion nobody noticed with it.
- **M6 — `a_reader_with_no_scripted_error_reads_its_whole_script_cleanly` is strictly dominated**
  (reliability, refactor_safety). Arming `new()` fails 60 tests; no mutation kills this one alone.
- **M7 — `script_position` is a third name for an index the type already calls `next_index`**
  (naming), and it borrows the word ng reserves for reference coordinates (`Position`). The
  `usize` does not disambiguate, because coordinates in these test modules are bare `usize` too.
- **M8 — `read_error_at` is a prepositional fragment**; the value is an index and the name never
  says so (naming).
- **M9 — two ordinals contradict each other and the enum** (naming). `cursor.rs` called the
  reference fetch "the second of the three"; `filtering.rs` called the decode "the third of the
  three"; `ReadFilterError` declares Source, Decode, Reference.
- **M10 — `step_one_failure`'s `.expect()` states a variant-wide invariant the variant lacks**
  (errors). `CursorError::ReadRecord`'s other construction site carries a raw `io::Error`.
- **M11 — the decode note omits that C2/C3 take the variant's last two constructions with them**
  (reliability, refactor_safety). `FakeRecord::decode_fails` is *not* dead.

### Nits

`.expect("an in-memory read cannot fail")` is no longer true of the type (4 sites); "short read"
collides with the domain term charged to `DropReason::TooShort`; `contains("read 1")` also matches
`read 10`; "rewinds into it rather than past it" reads backwards; three test names carry
"position" or do not name the alternative they rule out.

## 7. Out of scope observations

- `region_raw_aligned_reads.rs` propagates with a bare `?` and never names the **region** it was
  serving, though `self.region` is in hand. An operator gets the file and the reason but not
  where. A candidate for C2's rewrite of that loop.
- `read_named_with_length("overruns", 0, 95, 30)` passes three same-primitive `usize` domain
  scalars positionally — a transposition hazard, pre-existing in `test_fixtures.rs`.
- Whether ng's binaries walk `Error::source()` when reporting. If any top-level reporter prints
  only `{error}`, the whole chain collapses to "reading alignment file '…' failed". Not checked.

## 8. What the review confirmed rather than found

**All four agents were told to attack the author's central claim, and all four independently
confirmed it.** `ReadFilterError::Decode` is unreachable through the real chain:

- One agent asserted `decode(&*buf).is_ok()` at **both** of the narrowing's yield points and ran
  the whole suite — including the real BAM and CRAM walks, the SSR loci tests and the generic
  locus generator — without the assertion tripping once.
- One drove a **nine-record** adversarial script (no-start, no-contig, neither, unmapped-with and
  without a start, an out-of-range reference id, an empty-CIGAR/empty-sequence record, a
  non-string `RG` tag, one good read) through a cursor: no error of any kind.
- One checked the noodles side rather than assuming it — `RecordBuf::alignment_end` is
  `alignment_start().and_then(…)`, so `alignment_end().is_some()` *implies*
  `alignment_start().is_some()`, which is what makes the `overlaps` gate sufficient.

**And the D2 finding is worse than the author reported.** The author wrote that the plan's D2
would pass under the mutation it names. An agent hoisted the conversion above the pre-decode
filters and ran everything: **2,845 passed, 0 failed** — the hoist is invisible to the *entire
suite*, not merely to D2.

**Two corrections to the review's own brief**, both from agents contradicting the orchestrator:
`records_consumed` was cited as history when it is live in the file being edited, and the planted
objection to `step_one_failure`'s register was rejected with a reason that was accepted.
