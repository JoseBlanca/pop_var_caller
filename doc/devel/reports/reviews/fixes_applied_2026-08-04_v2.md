# Fixes applied — D2 review (the conversion is asked for nothing when the first filter drops)

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `5320fe4`
**Review:** [`ng_read_filtering_stages_d2_2026-08-04.md`](ng_read_filtering_stages_d2_2026-08-04.md)
— three agents in isolated worktrees. Raw findings under
`tmp/review_2026-08-04_d2-conversion-not-asked/`.

**All four Majors applied.** The mechanism kept its shape — count the conversion inside the
conversion — and changed its carrier, its reach and its prose.

---

## 1. M3 — the counter became `CursorCounts::reads_converted`

The thread-local, `reset_conversions_attempted`, `conversions_attempted` and their 17-line
justification block are **deleted**. In their place, one field on the struct that already exists to
say what the cursor did:

```rust
pub reads_converted: u64,          // beside reads_decoded

fn convert_buffered_read(&mut self) -> io::Result<AlignedRead> {
    self.counts.reads_converted += 1;
    self.buffer.decode()
}
```

Two reviewers built this independently and measured identical detection power. It has no reset
protocol, no per-thread caveat, is per-cursor rather than thread-global, and is folded across a
sample's files for free by the exhaustive-destructure `AddAssign` — whose own doc says a new field
"is now folded everywhere by construction". It also removes what would have been the crate's first
`#[cfg(test)]` static and first `#[cfg(test)]` statement in a production function body.

**It ships, and that is a deliberate widening of the step.** Milestone D was framed as tests-only;
that framing was already broken by extracting the method, and this makes the counter a real
observable rather than test scaffolding. Re-measured for it — §5.

`reads_decoded` gained the sentence it has needed since it was written: it counts reads surviving
**both** filters, not decodes, and the difference from `reads_converted` is exactly the second
filter's drops.

## 2. M1 — the test now covers all six of the first filter's reasons

The submitted test scripted three. A reviewer converted the `Unmapped` and `LowMapq` drops before
rejecting them and **all 2,869 tests stayed green**.

The fixture is now eight records — two kept and one per reason: duplicate, low MAPQ (the one that
is not a flag bit, and the highest-volume drop on real data), supplementary, secondary, unmapped
(placed, so it reaches the filter at all) and QC-fail — asserted against a **whole-struct**
`ReadFilterCounts` compare rather than four fields, so a seventh reason appearing is a compile
error. New helper `read_at_mapq`.

**Re-measured:** the per-reason hoist now fails, `left: 4, right: 2`.

## 3. M2 and M4 — the two false claims

- *"an increment beside the call … reporting zero conversions"* → **it reports 2**, the number the
  test asserts, so the test passes while every read is converted. I reproduced this myself before
  writing it down: with the increment left at the old call site and the call hoisted, both D2 tests
  pass. The doc now says that, and adds the qualification a reviewer supplied — a hoist that
  *carried* the increment would still be caught, so the hazard is the ordinary shape where one
  statement moves and its neighbour does not.
- *"the only way this design's central ordering can be tested at all"* → **false**; a reviewer built
  a zero-sized witness that makes the hoist a compile error. The doc now records the alternative,
  why it was not taken (spec §1's "does not … add a type", and that it pins only half the property
  — no type can forbid a second copy of the length rule), and that it is a Checkpoint D question.

## 4. The Minors

| finding | applied |
|---|---|
| "2,866 passed" three times, stale and misattributed | all three copies gone; the claim is now the reproducible one — each test's unique detection power |
| the same argument written out three times | one full statement, on `convert_buffered_read`; the test block points at it |
| opposite directional verbs, and test 2's doc naming its mutation backwards | "hoist / above" throughout, matching spec §8 and the checkpoint |
| test 2 cited spec §2 for a claim §2 does not make about #7 | rewritten: #7 **could** move — a raw record carries a sequence — and the real reason it stays is that which side of the conversion a filter runs on is `read_filtering.md`'s to decide, not this plan's |
| the `an_unmapped_read_…` cross-reference cannot support the footprint-less claim | reworded; it now says explicitly that nothing pins the footprint-less case because nothing can observe it |
| the impossibility argument restates `ReadFilterError::Decode`'s doc | cites it instead |
| assertion messages misdiagnose the instrument-removed case | both now name all three causes and how to tell them apart |
| `reads_decoded` does not count decodes | stated on the field, beside `reads_converted` |
| `flagged_at` is a participle with no noun; `=` not `|=` | `flagged_read_at`, and `|=` |

## 5. Validation, re-run after the rework

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,869 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**Mutations, all re-run against the reworked version:**

| mutation | result |
|---|---|
| conversion hoisted above the first filter | test 1 alone — `2868 passed; 1 failed`, `left: 8, right: 2` |
| filter #7 hoisted above the conversion | test 2 alone — `2868 passed; 1 failed` |
| the increment removed from `convert_buffered_read` | **both** — `2867 passed; 2 failed` |
| two drop reasons convert before rejecting (M1's mutation) | test 1 alone — `2868 passed; 1 failed`, `left: 4, right: 2` |
| the increment left at the call site while the call is hoisted | **survives** — which is the point of M2, and why the count is in the callee |

**Speed, re-measured because the counter now ships** — `ng_generic_walk_probe` on HG002 chr21, six
runs a side, one machine, one session, nothing else running:

| | runs (s) | median |
|---|---|---|
| before | 1.839, 1.846, 1.826, 1.827, 1.819, 1.826 | 1.8265 |
| after | 1.853, 1.848, 1.834, 1.840, 1.830, 1.838 | 1.8390 |

Medians differ by **0.7 %** and the ranges overlap — the fastest *after* run (1.830) beats the
slowest *before* run (1.846) — so this is noise by the bar B1 set, not a cost. One `u64` increment
per conversion, ~55k of them on chromosome 21.

## 6. Not applied

- **The witness type.** Recorded on `convert_buffered_read` and carried to Checkpoint D — it needs
  a spec §1 amendment, and it would not remove the counter.
- **Printing `reads_converted` in `ng_generic_walk_probe`.** Now possible, since the counter ships,
  and `CursorCounts` does not currently reach the probe at all — that is plumbing through
  `SampleCursor` and the stream, and it belongs to whoever wants the number on real data.
