# window coverage — C3 review: the look-ahead, its fixtures, and the machinery that watches it

**Date:** 2026-09-06/07
**Reviewed:** the working tree of plan step C3, at `a5767706 wip: C3 for review`
**Implementation report:** [ng_window_coverage_c3_2026-09-06.md](../implementations/ng_window_coverage_c3_2026-09-06.md)
**Fixes applied in:** C3's own commit

Two agents on the step's diff, each in its own worktree detached at the step: one on reliability,
errors and every number the step's prose claims; one on naming, idiom, module structure and
refactor safety. **1 Major both found independently**, 4 further Major, 9 Minor. All applied
except the three recorded below as raised.

## Verdict

**The look-ahead's arithmetic is right and its two new tests hold it.** Five deliberate defects in
it — removed, doubled, clamp removed, clamp taken from the region's start, look-ahead added to the
region's start — were each caught, by between 3 and 30 tests. Its sufficiency argument holds for
the reason the code gives: the record that stops a draw is not merely held, it is *observed* by the
same cover, so a sample with any record past the chain has had its centres closed.

**What was wrong was around it.** A cover that draws a sample onto the next contig ended holding
*that* contig's bases, so the accessor a builder will read at C4 answers `None` at every position
of the region it owns — found independently by both agents and by C2's own re-run correctness
review, which is three times from three directions. A fixture this step changed stopped testing its
own name. And every field the recorder writes, and the whole of the probe's comparison, had no test
at all: four deliberate defects in them survived the entire suite.

## Findings

### 1. Major, CONFIRMED — a cover leaves the reference buffer on the *next* contig

`observation_cache.rs`, `read_the_ground_and_measure_coverage_over_it`.

The contigs a cover reads are sorted ascending, which the accumulators require, and the doc then
claimed the region's own comes last. Those are the same order only while no contig sorts after the
region's — and one regularly does: a sample is drawn one record past the reach, that record can be
the first of the next contig, and the next contig sorts after. The cover then ends with the next
contig's bases in the buffer.

**Measured**, on the tree as committed: one record at `contig 0:40-45`, one at `contig 1:10`, cover
`contig 0:40-50` — `reference_base_at(contig 0:45)` answers `None` where the reference has a base.
Without the next-contig record the same cover answers it.

**Latent, not live**: `reference_base_at` has no caller outside tests until C4's builder. But the
look-ahead turns this from a corner into **the ordinary case at every contig boundary**, which is
why it is fixed here rather than left to the step that would first suffer it.

**Fix.** The region's ground is read again at the end of a cover that ended elsewhere, observing
nothing — every record on that contig was observed on its own turn. Pinned by
`a_record_held_on_the_next_contig_does_not_take_the_buffer_with_it`, which fails when the refetch
is removed.

**The alternative was rejected on evidence.** Refusing to read contigs past the region's would
defer the accumulator's contig-change flush — and with the look-ahead clamped at a contig's end,
that flush is the *only* thing that can close the last 250 bases of a contig, since no later
position on it exists. That variant trades a latent defect for a live one.

### 2. Major, CONFIRMED — a fixture this step changed stopped testing its own name

`callers.rs`, `a_refused_record_is_the_runs_answer_even_when_the_source_then_fails`.

The test's subject is a precedence: when the output sink refuses a record *and* the source then
fails, the refusal is the run's answer. C3 added a record on the next contig to stage the two. But
that record stops the draw for the whole of `chr1`, and the run's analysed ground is `chr1` alone —
so **the source failure was never drawn**, and the test asserted a precedence between two failures
of which only one happened.

**Measured**: swapping the two arms so the source failure wins left all 515 `ng::run` tests green.
On the parent's fixture the same mutation was caught.

**Fix.** This one test is given a segmentation over both fixture contigs, so the merge reaches past
the record that shields `chr1`. The precedence mutation fails the test again, naming the failure it
got.

**The general fact, worth carrying:** `chr1` in that fixture is a hundred bases and half a window
is 250, so **every** cover on it reaches its end. A fixture that needs a cover to stop short of a
contig's end now needs a contig longer than 250 bases; a smaller building region will not do it.

### 3. Major, CONFIRMED — every field the recorder writes was untested

Its only test asserted the *off* path, and the row-building loop is unreachable under test because
no fixture names a file. Three mutations applied at once — sample index `+1`, position `+1`, absent
rows omitted — **survived 6,362 tests**. Each is silently wrong in its own way: an index off by one
compares every sample against its neighbour's windows; a position off by one finds no centre at
all; an omitted absent row turns "the run had nothing here" into "the run never reached this
locus", which is one of the two halves of the failure the facility exists to see.

**Fix.** Row-building split out of the I/O and tested directly. Re-running the mutations against it:
caught.

### 4. Major, CONFIRMED — the probe's comparison had no test, including its bit comparison

Replacing the bit comparison with field-wise `f32` equality **survived** the example's five tests.
It is a real defect: an absent window is a pair of `NaN`s, so a correct run would be reported as
disagreeing at every centre the position floor silenced — a large share of a store at three reads a
position.

**Fix.** The verdict is a named function delegating to `WindowCoverage`'s own `PartialEq`, which is
already bitwise, with a test over two absent windows and two a last bit apart. The mutation is
caught.

### 5. Major — the row format was written in the library and re-derived by column index in the example

Swap two columns on the writing side and every field still parses, and the comparison then reports
each locus as one the run had no window for — **which is the shape of a pass in this measurement**.
Both halves live in one file now, with a round-trip test over an ordinary window, one the floor
silenced, and one the run had none for.

### 6. Major — the recorder was in the wrong module, under a name promising a check it does not make

Live `pub` code beside `production_parity.rs`, which is `#[cfg(test)]`; naming `WindowedCohort`;
called from exactly one place in the cohort merge; and performing no comparison — the comparison is
the probe's. It is `cohort_merge/recorded_windows.rs` now, which also takes the dependency out of
the `window_coverage` ↔ `cohort_merge` cycle rather than adding a third file to it.

### 7. Major — two doc comments this step made false, and one it made incomplete

`where_the_chain_starts` was named for the *start* of the chain and returns its *reach*; both call
sites read `let mut chain_reach = self.where_the_chain_starts(…)`, name and variable disagreeing in
one line. It is `where_the_chain_reaches_before_any_sweep`. And `cover`'s own doc still said the
reach "starts at `region`'s last base" — the sentence this step falsified — while two further
paragraphs, on what bounds a cover and what one costs, were left describing only the chaining term.

### 8. Major — a totalling path a new counter falls out of without a compile error

The window comparison's totals were copied field by field, where the same file's other totalling
function destructures precisely so a count added later cannot be left out. It destructures now,
with the two example lists named and discarded.

### 9. Minor — the probe told the reader a locus written twice is routine

It is not: a builder skips a locus it does not own rather than building it, so a duplicate means
either the ownership rule broke or two runs wrote into one file. A reader told duplicates are normal
would explain away the one thing a duplicate can reveal.

### 10. Minor, applied — numbers and prose

- **"at ten bases apart the look-ahead alone holds twenty-five of the thirty"** — measured **26**,
  by restoring the pre-C3 fixture and reading the failure. Corrected.
- **"finalised only once the stream has passed `p + 250`"** — the accumulator closes a centre at
  `centre + half <= position`, so it finalises on *reaching* it. Reaching is what makes half a
  window exactly sufficient rather than one base short, so the stricter word understated the code.
  Corrected in both places.
- **"the gap between them is wider than the look-ahead's half window"**, of the records at 800 and
  900 — that gap is 100. The load-bearing fact is that 800 is past the 750 the region's own
  look-ahead reaches. Rewritten.
- **"waits for the contig to change or for the pass to end"** — nothing in the merge calls
  `finish`, so the second half described something that does not exist. It now says which step
  supplies it.
- Smaller: a flag whose doc described an `Option` that was not there; two allocations per record in
  the recomputation, one of them a copy of every record's reference bases; a per-record linear scan
  over the reference's contig list; `(u32, u64)` map keys where `GenomePosition` fits and already
  orders by genome order; a usage block that never said where the recorded file comes from; and a
  `# Panics` section the recorder had not.

## Raised, not fixed

- **A mis-ordered store list is caught by nothing but luck.** The recorded file keys its rows by
  sample *index*, and the probe's only defence is a doc comment. Listing the stores wrong compares
  one sample's walk against another's windows — usually loud, but loud in the wrong words, naming a
  cache defect. Writing the sample's **name** into each row and checking it against the store
  header's would close it; that changes the row format, and both agents left it as the author's
  call.
- **A mutation survivor in `ground_this_cover_holds_on`**: on the region's own contig, widening the
  ground by the last held summary's reach is dead — the chain already dominates it — so replacing
  `last.reach()` with its start survives the whole suite. It is *not* dead on the other contigs a
  cover reads, where a multi-base record's ground would be cut short and the measurement would
  panic on a missing base. **No fixture has a multi-base record on a contig a cover is leaving.**
  This is the same code C2's own review found wrong from the other side (its finding 3), and both
  are answered by the same fixture; it is written up there and carried with C2's fixes.
- **A region already past its contig's end** now has a fixture
  (`a_region_past_its_contigs_end_is_not_clamped_back_behind_itself`), which closes the one
  look-ahead mutation that had survived — the clamp's floor taken from the region's start rather
  than its end.

## Mutations run

| mutation | outcome |
|---|---|
| look-ahead removed | caught — 4 tests |
| look-ahead doubled | caught — 7 tests |
| clamp removed | caught — 30+ tests |
| clamp's floor from the region's start | **survived**; now caught by the fixture added for it |
| look-ahead added to the region's start | caught — 16 tests |
| ground follows the region's end, not the chain | caught — 30+ tests |
| ground's per-sample widening uses `start_position()` not `reach()` | **survived**; raised above |
| `draw_to` remembers where the last sweep stopped | caught — 5 tests |
| overshoot by two records instead of one | caught — 6 tests |
| `with_observations` checks the span's left edge | caught — 1 test |
| the serial driver covers the whole analysed region | caught — 1 test, exactly the one whose doc claims it |
| the source failure wins over the record refusal | **survived**; now caught, finding 2 |
| the recorder writes `sample + 1` and `position + 1` | **survived** 6,362 tests; now caught, finding 3 |
| the recorder omits the absent rows | **survived** 6,362 tests; now caught, finding 3 |
| the probe compares the two `f32` fields with `==` | **survived** 5 tests; now caught, finding 4 |
| the refetch of the region's ground removed | caught by the fixture added for it, finding 1 |

## Numbers checked

Every figure in the step's prose was run rather than recalled. All correct but one: `WINDOW_BP` is
500 so half a window is 250; the fixture reference's contigs are 2,000 bases and the unclamped
fetch would ask for 2,150; six records held at the moment of failure against the thirty a
whole-region cover holds; the fixture's `chr1` is a hundred bases; the failure names
`region_on(9, 40, 300)` because contig 9 is unknown and so unclamped; 540 + 250 = 790. **The
exception is "twenty-five of the thirty", measured 26**, corrected above.

## Looked at and found sound

The look-ahead's placement — one helper, both fixpoints, so the two cannot drift; `MergeReference`
gaining `ContigTable`, with the fixture that wraps a reference forwarding it; the recorder's
concurrency (one mutex for the process, one `write_all` of a whole locus's rows so rows cannot
interleave inside a line, poison recovered rather than propagated) and its inertness when off,
confirmed by mutation; the probe's refusal to report a pass having compared nothing, and its refusal
of `--windows-from-the-run` without `--reference`; the probe's contig-name check, which is reachable
and fires on the first record of each contig; that the recomputation has not drifted from the run —
the same depth rule, a summary built from the head exactly as the psp path builds it, the same
accumulator constants from the same module, the same canonicalising fetch; and every moved fixture
still catching the defect its own doc names.
