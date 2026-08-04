# Code Review: ng alignment cursor — Milestone B

**Date:** 2026-08-02
**Reviewer:** rust-code-review skill (orchestrator + isolated sub-agents)
**Scope:** the working-tree diff of each Milestone B step, reviewed before its commit
**Status:** Request-changes, then approve — all findings applied

---

## B1 — the cursor, without the rule

### 1. Scope

The uncommitted change for B1 against `0cb0b94` plus its patch. In scope:
`src/ng/read/input/cursor.rs`, `src/ng/read/input/region_records.rs` (new),
`src/ng/read/filtering.rs` (the three-state flag), `src/ng/read/input/mod.rs`. **Two agents**
over four checklists: `reliability`; `module_structure` + `naming` + `smells`.

### 2. Verdict

**Request-changes, then approve.** **Fifteen mutations run, seven survived** — and the three
Blockers all say the same thing: *the machinery B2 was to be built against was not, in fact,
exercised.*

### 5. Top 3 priorities

1. **B1** — the entire kept-read walk is unreachable.
2. **B2** — `restart_after_end_of_input` has no test; deleting the guard that is the sole
   reason for three states left the suite green.
3. **B3** — the two overlap tests genuinely disagree on a read that reaches both.

### 6. Findings

#### Blocker

**B1: the kept-read walk never executes.** `kept` is cleared on every `move_to_region` and
`examined` set to `kept.len()` after every push, so `kept.get(examined)` is always `None`.
Proved independently by both agents — one put `panic!` in the loop body, the other
`unreachable!()`; all 13 cursor tests passed either way. That made `read_start`, the free
`overlaps`, the `examined += 1` and the early `return None` all dead, and spec §11.6 names two
of them as things that *must* fail a test when removed. **The step whose purpose was "build
the machinery first, so the rule is a diff against something known correct" had shipped
scaffolding nothing could execute.** *Categories: reliability, smells — convergent.*

**B2: `restart_after_end_of_input` had no test.** Making it clear `Failed` too — the one thing
it must never do — left `1503 passed; 0 failed`. The two existing tests that sound relevant
never call it. *Category: reliability.*

**B3: the two overlap tests disagree on a read that reaches both.** noodles maps a zero CIGAR
span to *no* span, so `alignment_end()` gives an all-soft-clip read the one-base footprint
`start..=start`: `RegionRecords` **accepts** it, while the cursor's hand-written `overlaps`
**rejected** it. A mapped 30-base all-soft-clip read clears every step-1 filter — filter #7
uses SEQ length, #9 explicitly exempts all-clip. Measured on a `30S` read at position 40
against region 31..=60: the cursor yields it while `overlaps` says false. At B2 the same read
becomes order-dependent — yielded when read fresh, dropped when replayed. *Categories:
reliability, smells — convergent.*

#### Major / Minor

- The "no state moves before the chromosome check" obligation was only half pinned: hoisting
  `kept.clear()` and `region` above the check survived.
- `region_records.rs` had **no test module at all** — the other-sample skip and the
  never-pointed-anywhere path were both unreached.
- The `end >= region.start` boundary survived mutation against the hand-written six-region
  table. **A property test over random scripts and random region sequences kills it** — direct
  evidence that a table of examples was not enough.
- `restart_after_end_of_input` carried `source_mut`'s entire 25-line doc as a prefix, asserting
  three things false of the function it now documented — including that clearing the flag "is
  not an accessor's" business, on the accessor that clears it.
- `read_start` was a one-line wrapper for `read.pos`; folding it in was −5 lines, clippy clean.
- `CursorError`'s doc still named a variant `Io`, renamed to `ReadRecord` two commits earlier.
- `kept` / `examined` are the bare participles `naming.md` forbids; `kept_reads()` returns a
  count.

### 7. Out of scope observations

- **A design gap for B2, found by mutation rather than by reading.** The "don't clear `kept`"
  mutation fails not by losing reads but by **double-yielding**: `RegionRecords::move_to`
  repositions the reader unconditionally and the filter arm does not check `kept` before
  pushing. So spec §4's "partly held — hand over the kept reads, then carry on reading, no
  jump" has **no code path**; it needs the layers to agree on where reading resumes, not just
  a rule about what to keep.
- The arch's "Module home" file list contains no `region_records.rs` — the type is prescribed,
  its file is not.

### 9. What's good

- **Both agents found the dead walk independently, by executing rather than reading** —
  `panic!` and `unreachable!()` in the same loop body. No amount of careful reading was going
  to produce `1503 passed` as evidence.
- **The property test earned its place immediately**: 2,000 cases green, and it kills a
  boundary mutation the hand-written table misses.
