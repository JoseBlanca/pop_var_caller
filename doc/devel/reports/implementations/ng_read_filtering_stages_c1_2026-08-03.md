# ng — read filtering in stages, C1: the in-memory reader can be scripted to fail

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `6e22718` (Checkpoint B)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **C1**.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §6, §8 ·
[arch](../../ng/arch/read_filtering_stages.md) §7.

---

## 1. Plan

C1's contract, from the plan:

> `InMemoryAlignedReadsReader` gains a scripted error — a read position at which it returns
> `Err` instead of a record. Re-point the three fatal-error tests
> (`read_filter_source_read_error_is_fatal`, `read_filter_decode_error_is_fatal`,
> `read_filter_reference_error_mid_stream_is_fatal`) through it, so they run the real chain
> rather than a double that bypasses two layers. **Before C3, which deletes the doubles.**

Built in three parts:

1. `InMemoryAlignedReadsReader::with_read_error_at(script_position)` and the
   `read_error_at: Option<usize>` field behind it, with four tests of the mechanism itself.
2. Two of the three fatal-error tests re-pointed into `input/cursor.rs`, driving the whole
   chain — `AlignedReadsReader` → `RegionRawAlignedReads` → `ReadFilter` → `AlignmentCursor`.
3. The third **deleted with no successor**, because the chain it would be re-pointed through
   cannot produce the error it tests. §2 has the measurement.

Arch §7 left the shape of the scripted error open — "an error at position *n*, or an arm that
always fails". The plan's own wording decides it (*"a read position at which it returns
`Err`"*), and position-*n* is the strictly more capable of the two: an always-failing arm
cannot produce a fault that arrives **mid-walk**, which is the interesting one and the one both
re-pointed tests use.

## 2. Assumptions, and one measured finding that changed the step

### 2.1 `ReadFilterError::Decode` is unreachable through the real chain — measured, not argued

`read_filter_decode_error_is_fatal` **cannot be re-pointed**, and this is the step's main
finding.

`RawAlignedRead::decode` refuses exactly three things
([`aligned_read.rs`](../../../../src/ng/read/aligned_read.rs)): a record with no reference
sequence id, one with no alignment start, and a buffer with no read group stamped.
`RegionRawAlignedReads::read_next` guarantees all three *before* it yields:

| what `decode` refuses | why the narrowing has already excluded it |
|---|---|
| no reference sequence id | it drops anything whose `reference_sequence_id` is not `Some(this contig)` |
| no alignment start | `overlaps` returns `false` unless both `alignment_start` and `alignment_end` are `Some` |
| no read group | the group is resolved and stamped on the record actually handed over |

So every buffer that reaches the conversion decodes.

**Driven as an experiment rather than reasoned about.** A cursor scripted with a
placed-but-unstarted record and a started-but-unplaced one, plus one clean read, yields:

```
outcome=["ok:clean"]
tally=[(Some(ReadGroupId(0)), ReadFilterCounts { kept: 1, .. all zero .. })]
```

Neither an error nor a read for the two broken records: both are discarded by the narrowing,
**below the filter and uncounted**.

What that leaves is an error variant that is defence in depth against a regression in the layer
below, not a response to any input. Keeping the variant is right; keeping a test double
substituted for the very chain that makes it unreachable is not — that double is what C1 exists
to remove. The finding is recorded as a comment block where the test was, so that "this error
cannot happen" is a dated claim with its evidence attached rather than folk knowledge.

**Spec §8 sanctions the reduction**: it lists *one* replacement test — "a scripted read error
still surfaces as fatal, through the real chain" — and says it "replaces the three test-double
tests". The plan's C1 wording is the more ambitious of the two. C1 delivers two.

### 2.2 It bears on the plan's D2, which is why it is carried to Checkpoint C

D2 proposes proving the conversion is not hoisted above the first filter by using *"a read that
would fail to convert: unmapped, with no alignment start"*, asserting *"a clean drop charged to
`unmapped` and no error"*, and mutation-verifying by hoisting the conversion.

No such read reaches filter #5. With no alignment start it has no footprint, so the narrowing
drops it first — uncounted, so not charged to `unmapped` either. An unmapped read that *does*
reach #5 has a start and a contig, so it would convert perfectly well if the conversion were
hoisted. **D2 as written would pass under the mutation it names**, which is the failure mode
this branch has hit eight times. Raised at Checkpoint C; Milestone D is a separate session and
nothing here pre-empts it.

### 2.3 Both re-pointed properties were already covered on real inputs — C1 buys something narrower

The plan's argument for C1 is that the doubles "bypass two layers". True, but it does **not**
follow that the chain was uncovered, and the first draft of this report said so anyway. The
review corrected it, from three agents independently:

- `open_bam.rs::t10_a_truncated_file_fails_once_and_then_refuses_later_regions` truncates a real
  indexed BAM mid-walk and asserts the whole contract — reads flow, one error, the walk stops,
  the next region is refused.
- `cursor.rs::a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short`
  already drove the reference failure through the whole chain, over the **identical** `overruns`
  fixture the new test uses, thirty lines above it.

Both die under "the layer in between swallows the fault", so that property had two real pins
before C1.

**What C1 actually buys, measured:** no test anywhere pinned *which* `ReadFilterError` a fault is
charged to, because both of those match `Err(_)`. Charging a read failure to `Decode` kills
exactly one test in the tree — the new one. Charging a reference failure to `Source` kills exactly
one — the other new one. A scripted fault is a fault whose kind the script chose, so the assertion
can name it; that is the whole return on the mechanism, and it is worth stating accurately rather
than claiming a gap that was not there.

### 2.4 A scope addition, on the owner's ruling

C1 gained a second fault-injection knob, `with_failing_seek`, and a test for it. The plan does not
mention one: it asks for "a read position at which it returns `Err` instead of a record". But a
reader can break in **two** places, and only one of them is a read — on a BAM, `begin_region` runs
an index query, so a corrupt index fails the *move* and no read is attempted. C1's review found
that route completely uncovered, by a mutation that survived the whole suite.

It is one field on the struct C1 was already editing, so it landed here rather than becoming its
own step. The behavioural repair that the missing test then exposed did **not** land here — see
§6.

### 2.5 Placement — a deviation, recorded

The plan does not say where the re-pointed tests live. They went to
[`input/cursor.rs`](../../../../src/ng/read/input/cursor.rs), not `filtering.rs`, for two
reasons: `filtering.rs`'s own module doc states that it no longer knows what a BAM is and the
drop-tally fixture was moved to `cursor.rs` at Milestone F for exactly this reason; and C2
deletes `ReadFilter`, so a test left in `filtering.rs` driving it would have to move again one
step later.

## 3. Changes made

### `src/ng/read/input/aligned_reads_reader/in_memory.rs`

- **New field `read_error_at: Option<usize>`** and **new method `with_read_error_at`**, a
  consuming builder so a scripted reader reads as one expression.
- `read_next` checks the fault **before consulting the script**, so a position past the end of
  the script is the *truncated* file rather than a clean stop, and does **not** advance
  `next_index`, so a file that cannot be read stays unreadable across both re-reads and
  repositions.
- A `# A read can be scripted to fail` section on the type's doc saying what the mechanism is
  for.

### `src/ng/read/input/cursor.rs`

- **`cursor_over_failing_at(records, script_position)`** — a cursor whose reader breaks at a
  chosen read.
- **`step_one_failure(&CursorError) -> &ReadFilterError`** — recovers the step-1 error the
  cursor wrapped, so a test can say *which* of the three fatal conditions fired rather than
  only that something failed.
- Two re-pointed tests, under a section header explaining what the doubles could not see.

### `src/ng/read/filtering.rs`

- `ErroringSource`, `read_filter_source_read_error_is_fatal`,
  `read_filter_reference_error_mid_stream_is_fatal` and `read_filter_decode_error_is_fatal`
  deleted, each replaced by a recorded note naming its successor — or, for the third, the
  measurement showing it has none.

## 4. Tests added / updated

| test | what it pins |
|---|---|
| `a_scripted_fault_fires_at_its_own_read_and_not_before` | the fault is at the scripted read **and not before** — the half that distinguishes this from an always-failing arm |
| `a_scripted_fault_is_not_consumed_by_reading_or_by_repositioning` | a fault does not heal: re-reading fails again, and `begin_region` does not clear it |
| `a_scripted_fault_survives_a_rewind_at_the_same_read` | at read 0 "rewound to the start" and "rewound into the fault" are the same observation; at read 1 they are not |
| `a_fault_scripted_at_or_past_the_end_of_the_script_fails_rather_than_ending_cleanly` | three cases — at the end, well past it, and an empty script — so the truncated file is never reported as a finished one |
| `a_scripted_seek_failure_breaks_the_reposition` | the second way a reader can break: the *move* fails and no read is attempted |
| `a_failure_reading_off_the_file_is_fatal_through_the_whole_chain` | a read failure is charged to **`Source`**, yielded once, fuses the walk, and refuses later regions |
| `a_reference_fetch_failure_mid_walk_is_fatal_through_the_whole_chain` | the same for filter #8's fetch running off the contig end, charged to **`Reference`** |
| `a_reposition_that_fails_is_refused_rather_than_answered` | a failed reposition is reported, not swallowed — the mutation that survived the whole suite |

Both fatal-path cursor tests deliver a clean read **first**, so a chain that never delivered
anything would not satisfy them. The variant each names is the assertion no other test in the
tree makes.

**Suite: 2,842 → 2,847 (+5)** — `in_memory.rs` +5 (five new, one strictly-dominated one deleted),
`cursor.rs` +3, `filtering.rs` −3. Fully accounted; no unexplained movement.

## 5. Validation

Run on the host, in debug, per the plan's precondition that the default gate is not the
right one here.

| command | result |
|---|---|
| `cargo fmt` | clean |
| `cargo test --lib` | **2,847 passed**, 0 failed, 5 ignored |
| `cargo test --examples` | ok (5 targets) |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0, no output |

**The four acceptance dumps, compared with `cmp` and not by line count:**

| dump | result |
|---|---|
| `ng_generic_loci_dump` HG002 chr21 (BAM), 251,792 lines | byte-identical |
| `ng_ssr_loci_dump` HG002 chr21 (BAM), 4,406 lines | byte-identical |
| `ng_generic_loci_dump` tomato SL4.0ch01 (CRAM), 1,718,914 lines | byte-identical |
| `ng_ssr_loci_dump` tomato SL4.0ch01 (CRAM), 11,945 lines | byte-identical |

`ng_generic_walk_probe` chr21: `loci=236081 observations=251786 reads_admitted=54709` — exact.
C1 touches no production path, so the `seconds` reading is run-to-run variation and not a
measurement of anything.

### Mutations run, and what each killed

Six by the author, **33 more across the review's four agents**, each in its own worktree. Every
mutation's marker was `grep -c`-confirmed present before the run, because `cargo fmt` reformats
between an edit and a substitution and a pattern that matched nothing looks exactly like a
surviving mutation.

| mutation | killed |
|---|---|
| fault fires on every read | `a_scripted_fault_fires_at_its_own_read_and_not_before`, `a_fault_scripted_at_or_past_…` |
| script consulted before the fault | `a_fault_scripted_at_or_past_the_end_…` |
| fault clamped by `==` rather than to the script's end | `a_fault_scripted_at_or_past_the_end_…` (well-past case) |
| fault consumed | all three scripted-fault tests |
| rewind shifts the fault relative to the script | `a_scripted_fault_survives_a_rewind_at_the_same_read` |
| **the middle layer swallows the fault** (`self.reader.read_next(buf).unwrap_or(false)`) | `a_failure_reading_off_the_file_…` **and** the pre-existing `t10_a_truncated_file_fails_once_and_then_refuses_later_regions` |
| source failure mis-charged as `Decode` | `a_failure_reading_off_the_file_…`, **alone in the tree** |
| reference failure mis-charged as `Source` | `a_reference_fetch_failure_mid_walk_…`, **alone in the tree** |
| reference failure becomes a silent drop | the new test **and** the pre-existing `a_cursor_whose_file_failed_refuses_later_regions_…` |
| **failed `jump_to` swallowed, `Ok(())` returned** | **survived the whole suite before C1's review**; now killed by `a_reposition_that_fails_is_refused_rather_than_answered`, alone |
| `new()` arms every reader | 60 tests |

## 6. What the review changed, and what it handed on

The full account is in the [review](../reviews/ng_read_filtering_stages_c1_2026-08-03.md) and the
[fixes report](../reviews/fixes_applied_2026-08-03_v6.md). The three that matter here:

- **A surviving mutation on a fatal path.** `move_to_region`'s failed reposition could be
  swallowed with the whole suite green. C1 now carries `with_failing_seek` and a test that kills
  it. Writing that test also exposed a **live defect**: the cursor commits `region` and
  `last_region_start` before the fallible jump, so a failed reposition leaves it serving from an
  unknown file position, and the next forward region then reuses rather than jumping. **Owner's
  ruling: C2 fixes it**, with its `failed` flag widened to cover both routes into a stopped
  cursor — see §2.4.
- **This report's own §2.3 was wrong** and has been rewritten. C1 buys variant discrimination,
  not chain coverage.
- **`with_failure_at_read` was silently inert** for a fault scripted past the end of the script.
  The review's suggested fix did not work either; the one that holds clamps to the script's end.

Still open, and not C1's:

- **`ReadFilterError::Decode` cannot be given a test** through the production chain (§2.1). The
  variant stays as defence in depth. **`FakeRecord::decode_fails` is not dead** — two surviving
  tests set it, and they are the only constructions of `Decode` in the tree. Both die at C2/C3,
  after which rewriting the decode arm as `Err(_) => continue` will survive the whole suite.
  Recorded at the code so it is a decision rather than a discovery.
- **`ReadFilterError::Source` names two unrelated fatal conditions** — an I/O read failure and an
  unresolvable `RG` tag — so arch §4's table is wrong. → Checkpoint C.
- **The plan's D2 does not survive §2.1**, and is worse than §2.2 stated: a reviewer hoisted the
  conversion above the pre-decode filters and the **entire suite** stayed green. → Checkpoint C.
- **Nothing in `filtering.rs` was re-homed**: the module still owns `ReadFilter` and its state
  machine, which is C2's subject.
