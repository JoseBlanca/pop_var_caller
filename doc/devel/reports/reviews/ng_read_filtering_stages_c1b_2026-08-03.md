# Code review — ng read filtering in stages, C1b (`container.rs` test module)

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `c718a1c` (C1)
**Impl report:** [`ng_read_filtering_stages_c1b_2026-08-03.md`](../implementations/ng_read_filtering_stages_c1b_2026-08-03.md).

---

## 1. Scope

The uncommitted working-tree diff for C1b — one file, +301 lines, entirely inside
`#[cfg(test)]`: a new test module for
`src/ng/read/input/aligned_reads_reader/container.rs`, which had none.

**Categories dispatched**, one `general-purpose` agent each, in its own git worktree:

| category | why |
|---|---|
| `reliability` | the change *is* tests; whether they can fail is the whole question |
| `naming` | new names entering the vocabulary Milestone A settled, plus substantial new prose |

Two rather than four, in proportion to a test-only diff. Each agent detached at `c718a1c` and
applied an exported patch, since the change was uncommitted.

## 2. Verdict

**Approve with changes** — all applied.

The tests were sound. Both findings that mattered were about what the module *claimed*, plus two
more surviving mutations the author had not looked for.

## 3. Execution status

Both agents reproduced the gate independently: `cargo fmt --check` clean, `cargo test --lib`
**2,856 passed / 0 failed**, `cargo clippy --all-targets --all-features -- -D warnings` exit 0.

**Mutation coverage.** The reliability agent ran **47 mutations** through a harness that verified
each marker was present exactly once, confirmed the file changed, ran the tests, and byte-restored
— asserting `restored == pristine` at the end. The naming agent ran instrumented probes rather
than reviewing by reading.

**The author's four reported mutations reproduced exactly** — kill set for kill set, not
approximately.

## 4. Top 3 priorities

1. **M1 — the module's central justification was false**, and correcting it makes the tests
   *more* important, not less.
2. **M2 — six of seven scalars could be read from the wrong entry** with the whole suite green.
3. **M3 — `Span::new`'s documented overflow refusal had zero coverage**, and is reachable in
   three lines.

## 5. Findings

### Major

#### M1: the module doc's "every caller passes a fresh buffer today" is false
**Category:** naming · **Confidence:** High (verified by running instrumented code)

Three places in the new prose claimed all four properties were latent, because every caller passes
a fresh buffer, and that a buffer *with a history* is a condition C2 would introduce. The claim was
inherited from B2's deferred finding and repeated without checking.

The agent instrumented `fill_raw_read` and ran the existing CRAM cursor walk:

```
PROBE i=0 incoming_seq_len=0  entry_seq_len=30
PROBE i=1 incoming_seq_len=30 entry_seq_len=30
…
PROBE i=39 incoming_seq_len=30 entry_seq_len=30
```

`ReadFilter::next` refills **one** buffer for a whole pass — as `ReadFilterBuffers::record_buf`'s
own doc already says — so every read after the first has arrived with a history since B2, not from
C2 onward. Asserting freshness instead of printing turns **three existing tests red**.

**The real reason the gap went unnoticed is different and worse:** no existing test compares a
served read's *content* against an independent expectation. With `sequence.clear()` deleted,
sequences grow past 300,000 bases and the pre-existing suite still reports `2,847 passed` — the
CRAM-versus-BAM oracle included. So a regression in these clears would corrupt production reads
**today**, silently, and only this module would catch it.

Two neighbouring claims were checked and held: deleting `out.data_mut().clear()` does leave the
pre-existing suite green, and a SAM `*` SEQ line does decode to an empty sequence.

#### M2: six of the seven scalars can be read from the wrong entry undetected
**Category:** reliability · **Confidence:** High (mutation-tested against the full suite)

Replacing `entry` with `self.index[0]` for `flags`, `reference_sequence_id`, `mapping_quality`,
`mate_reference_sequence_id`, `mate_alignment_start` and `template_length` each leaves the **entire
library suite green**. Only `alignment_start` dies — and only in the real-CRAM walks, not in the
new module.

The cause was the fixture: `full_record` varied only the name and the bases, so every record of
every multi-record fixture carried identical scalars. The module pinned *span* per-entry-ness and
not *scalar* per-entry-ness. This is the module's own failure mode one field along — a read
reporting another read's position, mate, mapping quality or template length.

#### M3: `Span::new`'s documented overflow refusal has zero coverage
**Category:** reliability · **Confidence:** High

Replacing the checked conversion with `start as u32` leaves the full suite green. The doc names the
failure it prevents — handing back another record's bytes — and nothing guarded it. `Span::new` is
private but the test module is in the same file, so it is directly callable.

### Minor

- **`shrink_to_fit` is unobservable and unpinned** (reliability). Deleting all three calls leaves
  the suite green, while `decode_container_at`'s comment quantifies what it buys (1.6 MiB of
  5.4 MiB per open file). To a reader who has not read that comment it looks like dead code.
- **The allocation-reuse property — the module's entire rationale — was untested** (reliability).
  The nine tests pinned that reuse is *safe*; none pinned that it *happens*. Replacing a field
  wholesale instead of clearing and refilling is behaviourally identical and leaves the suite
  green.
- **`filled` is a bare participle** (naming) at seven call sites: `filled(&container, 2)` does not
  say whether it returns a read, a buffer, a count or a bool. It also collides with
  `fill_raw_read`, which does the opposite — fills a caller's buffer rather than making one.
- **"span" collides with "reference interval"** (naming), which is what it means in the sibling
  module and in this file's own module doc. A test name is read in failure output, away from the
  declaration that disambiguates it. The name was also narrower than its assertions.
- **"short read" is a charged term** (naming): in this crate it means below `min_read_length`,
  charged to `DropReason::TooShort`. The fixture is 4 bases — genuinely under the default 30 — so
  the misreading is live, and points at a subsystem the test has nothing to do with.
- **"Four input classes" mislabelled two of its four items** (naming) — one is a property, one is a
  line of code.
- **`push`'s error path is unreachable** from any test (reliability): it needs a 4 GiB payload.
  Recorded rather than left looking uncovered.
- **The deferral's own text disagreed with itself** (reliability): §4 of the origin report lists
  "the clear-and-refill claim" as the third deferred item, §5 lists "the buffer-shrink claim".

### Nits

Two names understated their assertions; `expect` messages read as labels rather than statements;
one doc lead restated its own test name; plan step codes appear undefined; `fill_raw_read`'s panic
on an out-of-range index was neither documented nor pinned.

## 6. What the review confirmed rather than found

- **No test is strictly dominated — all nine earn their place.**
  `records_read_their_own_slices_of_the_shared_buffers` was the one the agent expected to delete:
  no unique kill, and a subset of another test's kill set under 30 mutations. It survives on a
  mutation written specifically to break that relation — `index[i]` → `index[(i+1).min(len-1)]`,
  serving every read as the *next* one — which the two-record test passes and the three-record
  fixture fails. Recorded because "three records instead of two" is exactly the kind of test that
  looks redundant and is not.
- **`quality_scores.clear()` and `cigar.clear()` do not ride along** on the sequence clear — the
  question the brief raised. Each of the four clears is independently load-bearing, verified by
  four separate single-line deletions. `name.clear()` is killed by one test alone.
- **The "four unreached classes" claim is correct.** Each of the four probing mutations is killed
  only by tests inside the new module, so the other 2,852 tests do not touch them — the same
  evidence the author had, obtained independently.
- **`a_container_counts_the_records_it_packed`'s second assertion is not a tautology**, though it
  looks like one: its real content is uniquely killed by adding `other_sample_records += 1` to
  `push`.

## 7. Out of scope observations

- A property test over the round trip (`push` + `fill_raw_read` is a round-trip law over a
  structured domain) was specified but not compiled by the agent. It would subsume the scalar gap
  automatically and reach orderings the fixed fixtures do not — long→short, short→long,
  unnamed→named, repeats of one index.
- The non-default validation gate is worth recording somewhere durable: a reviewer who runs the
  obvious `cargo test --all-targets` gets a red tree and no way to tell it is pre-existing.
