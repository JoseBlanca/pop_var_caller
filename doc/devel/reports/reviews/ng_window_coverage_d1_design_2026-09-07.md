# window coverage — D1 review (naming, structure, prose): the measurement is right, its shape and half its comments are not

**Date:** 2026-09-07
**Reviewed:** commit `444672b0 wip: D1 for review`, diffed against `cdfd7821`
**Remit:** naming, idiom, code smells, refactor safety, module structure, and the prose the step
ships. Correctness, error handling and the numbers are a separate reviewer's.
**Implementation report:** [ng_window_coverage_d1_2026-09-07.md](../implementations/ng_window_coverage_d1_2026-09-07.md)
**Fixes:** applied in this worktree, listed at the end.

## Verdict

**The measurement's central idea is sound and I could not find a way to make it lie.** Reading the
spread off the shipped accumulator — one instance per candidate floor, and the count it silences at
floor `F` *is* the number of windows holding fewer than `F` positions — means nothing in the probe
reimplements the window, which is the only way a probe of a rule can avoid checking one copy of the
rule against another. The four assertions in `close` are well chosen: they catch a sweep fed a
different stream, a sweep disagreeing with the walk where both answer the same question, an empty
store, and a non-monotonic distribution.

**What was weak was the shape around it and the comments on top of it.** One field of state —
"a total has no accumulator to build" — propagated into a second constructor, an `Option` around
every arm's accumulator, two `expect`s, a destructuring arm that had to explain why an accumulator
is not summed, and a comment describing a clone that did not happen. The file's own module doc,
which states how the file is organised, was left describing a file with two measurements after the
step added a third. A constant that the step spent 40 lines un-provisioning is still called
provisional 45 lines further down the same file. And the phrase "moves by 6 in 100" carries a
*relative change in a rate* in three documents whose every other "in 100" is a *share of windows*.

**8 findings applied, 1 left for its owner.** 4 Major, 5 Minor.

## Findings

### 1. Major — the module doc's "One walk, several measurements" describes a file that no longer exists

The section is the file's own statement of how it is organised, and it ended:

> `SingleBaseEquality` is the first question's and `WindowRecomputation` the second's; neither knows
> about the other, and a third would be a third implementor and one more entry in `main`'s slice.

The step adds a third question and does not add a third implementor. A reader who takes the section
at its word goes looking for a `FloorSweep: RecordMeasurement` that is not there.

**Riding inside `WindowRecomputation` is the right call, and the doc rather than the code is what
should have moved.** What the sweep counts is not a record: it is a covered position with its
reference base and its depth, and those exist only after the recomputation has looked the reference
bases up over the record's ground and run the depth rule over the record's body. A separate
`RecordMeasurement` would need its own `WindowedRefSeq` accessor, its own contig check against the
store's contig list, its own `bases`/`reported` scratch buffers and its own `for_each_reported_depth`
call — a second copy of the feeding logic, around 60 lines, plus a second reference accessor per
sample. Worse, `close`'s check of the arm at the shipped floor against the walk's own count of
absent windows would then compare two independently fed streams, so a disagreement would no longer
mean what the assertion message says it means. The cheap check exists precisely because the two are
fed from one place.

**Fixed** by extending the section: the third question rides inside the second, why the record is
the wrong granularity for it, and what a genuine third implementor would still be.

### 2. Major — `WindowCoverageConfig::min_window_positions` still calls the value provisional

`src/ng/window_coverage/mod.rs` line 171, untouched by the step:

> Fewest covered positions a window may be built from and still report a number; below this it comes
> back absent. **Provisionally** [`MIN_WINDOW_POSITIONS`] — see spec §3.3.

The step's whole point is that it is no longer provisional, and it rewrote the constant's own doc to
say so — but a reader arriving at the config type, which is where the field is actually set, is told
the opposite. Two doc comments 45 lines apart now contradict each other about the same number.

**Fixed**: the field points at the constant and says the value is measured, with the constant's doc
named as where the cost of each setting lives.

### 3. Major — "moves by 6 in 100" uses the document's own unit for a different quantity

The sentence appears in three places (the constant's doc, spec §3.3, the report's "The answer"):

> The last two rows are the same ground at a sixth of the depth and the answer moves by 6 in 100;
> between the third row and the fourth the intervals shorten 42-fold and it moves 6,000-fold.

Every other "in 100" and "in 10,000" in these documents is a **share of windows**. This one is a
**relative change in that share** — 360.13 per 10,000 against 381.69 is +6%. A reader who carries
the established unit across reads "the depth change silenced 6 more windows in every 100", which is
about 600 times the real effect and points the wrong way for the argument the sentence is making
(that depth barely matters). The same clause is quoted again in the report's "what this measurement
could not say", where it does load-bearing work.

**Fixed** in all three: the numbers are shown rather than the ratio between them — "the share
silenced at 50 barely moves: 360 windows in every 10,000 at 30 reads a position against 382 at 5" —
and the interval comparison is spelled the same way, "from 0.06 in 10,000 to 360, six thousand times
as many windows".

### 4. Major — `FloorSweep::empty()`, and the four other things it forced

`empty()` exists for one reason: `WindowTotals` has to sum stores' answers into something, and
`OneFloor` could not be built without an accumulator that a total has no use for. That one
requirement showed up five times over:

- a second constructor whose only job is to reproduce `new`'s floor list;
- `accumulator: Option<WindowCoverageAccumulator>` on every arm;
- `expect("positions are observed before the sweep is closed")` in `observe`;
- `expect("the sweep is closed exactly once")` in `close`;
- an `accumulator: _` arm in `add`'s destructure with a two-line comment explaining why an
  accumulator is not summed.

And the comment at the call site described an operation that does not happen:

> **The first store's arms are cloned into by summing an empty bank into them**

Nothing is cloned; an all-zero bank is built and added to.

**There is a form with one constructor, and it removes all five.** The live bank and its answer are
two different things with two different lifetimes: `FloorSweep` holds the accumulators and exists
only for the pass; `WindowsUnderEachFloor` — a `Vec` of `(floor, windows, silenced)` — is what gets
printed and summed. `close` now takes `self`, consumes the accumulators, and returns the answer, so
neither `expect` can be written and there is no half-closed sweep to describe. The total holds
`Option<WindowsUnderEachFloor>` and clones the first store's answer in, which is what the comment
already claimed.

`WindowRecomputation` grows a second field for the answer — the same `accumulator`/`histogram` shape
that struct already uses for exactly this reason, so the pattern is the file's own.

### 5. Minor — `FloorSweep::new` had no doc comment; `empty` had a paragraph

The constructor that builds 25 accumulators, and the one place the `..` over the run's configuration
is actually written, carried nothing — while the one that builds an empty list carried three lines. `silenced_at` in the tests
also carried none. **Fixed**: `new` documents what it builds and why the `..` is safe there;
`silenced_at` gets one line. I checked the whole diff for the defect this branch's review history
records — a new function landing inside the previous item's doc comment — and **found none**: every
other new item (`the_configuration_a_run_uses`, `FloorSweep`, `OneFloor`, `CANDIDATE_FLOORS`,
`the_floors_asked_about`, `observe`, `close`, `print`, `add`, and all six tests) carries its own.

### 6. Minor — `the_floors_asked_about`'s doc is in the future tense about the step that has happened

> …cannot quietly stop applying **when plan step D1 moves the default**.

D1 is this step, and it has ruled that the default does not move. The sentence also leans on an
internal step label to do work a plain phrase does. **Fixed**: "if `MIN_WINDOW_POSITIONS` is ever
changed", with an intra-doc link.

### 7. Minor — `CANDIDATE_FLOORS`'s doc: a subjectless first sentence, and an invariant nothing checked

> The floors the sweep asks about — closest together over the range spec §3.3 is choosing in.

*Closest together* has nothing to be close to; the sentence needs "spaced most closely". Two lines
later the doc states a property in the language of the assertions beside it:

> A floor of 1 silences nothing … so an arm there reporting anything above zero is a defect and not
> a reading.

Nothing checked it. Three of the four properties `close` cares about are asserted and this one, also
stated as a defect, was left to the reader. It is true — `finalise_centre` counts the centre itself
into `distinct_positions_summed`, so no window can hold fewer than one position — and checking it is
one comparison.

**Fixed**: sentence reworded, and a fifth assertion added at floor 1. `close`'s doc and the report's
list of checks both now say five.

### 8. Minor — the printed rows: a duplicate row, a duplicate name, and three untagged numbers

Two separate problems, and neither breaks a machine reader, because **there is none**: nothing under
`scripts/`, `tests/` or `benchmarks/` parses this probe's output, and its only references are the
implementation plan and the reports. Humans read it.

- `sweep-windows` printed the window count a third time, under a fourth name. `close` asserts every
  arm's count equal to the walk's, and both `WindowRecomputation::print` and `WindowTotals::print`
  already emit that number as `windows-finalised`. **Removed**; `print`'s doc now says the
  denominator is the row above and why that is safe.
- The floor rows were `label`, key, then three bare numbers — floor, silenced, per-10,000 — where
  every other row in the file is `label`, key, value. The file already has a precedent for a wider
  row: the counter-example rows put the identifier third and **tag** what follows
  (`head=`, `evidence=`, `partial-witness=`). **Fixed** to match:
  `…\twindows-under-a-floor-of\t50\twindows=867\tper-10000=1.13`.

`FloorSweep::print`'s doc also opened "One row a floor" while printing a header row that was not one
row a floor; that goes away with the header.

### 9. Minor — "arm" names two different things in the report; and `store` names two in `add`

The report uses *arm* for one accumulator at one candidate floor, and then twice for a missing
measurement case ("the one arm this measurement does not have"). One term, two concepts, in one
document. **Fixed** to "case".

`WindowTotals::add` wrote `if let Some(store) = &store.sweep`, shadowing a `&WindowRecomputation`
with a `&FloorSweep` inside a five-line body that uses both meanings. **Fixed**.

### 10. Not fixed — the report's **Review:** link points at a file that does not exist

`doc/devel/reports/reviews/ng_window_coverage_d1_2026-09-07.md` is not in the tree. That is the other
reviewer's filename, so I have left it alone rather than repoint it; this review is at
`ng_window_coverage_d1_design_2026-09-07.md` and should be linked beside it when the two are merged.

## Looked at and found sound

**`the_configuration_a_run_uses` keeps the guarantee the C5 review praised, and `..` does not
weaken it.** C5 praised the literal for spelling "the same six constants … with no `..` so it cannot
drift from the run's". The function still spells all six with no `..`, so a field added to
`WindowCoverageConfig` is a compile error there — the guarantee is intact, and now in one place
instead of two. The `..the_configuration_a_run_uses()` in `FloorSweep::new` is a different thing from
the hazard C5 was guarding against: the base is *the run's own configuration*, not a `Default`, so a
field added later reaches every arm carrying the run's value, which is exactly what an arm that
differs in one setting wants. The only sharp edge was that the function's doc said "no `..`" while
the code beneath it wrote one; I added a sentence saying the override is the intended use.

**The name says what the value is.** `the_configuration_a_run_uses` returns a configuration and reads
as a noun phrase, in keeping with this file's `the_floors_asked_about` and the module's
`histograms_beside` / `read_a_row`.

**Refactor safety, on the three cases asked about.**

| change | what happens |
|---|---|
| a field added to `WindowCoverageConfig` | compile error in `the_configuration_a_run_uses` — every field named, no `..` |
| a variant added to `SampleHistogram` | the sweep discards each arm's histogram (`let (tail, _) = accumulator.finish()`) and is unaffected; `report_the_histogram`'s match is exhaustive with no wildcard, so the compiler still catches it where it matters |
| `MIN_WINDOW_POSITIONS` renamed | compile error at all three uses |
| a floor added to `CANDIDATE_FLOORS` | compile error on the declared length `[u32; 25]` |
| a field added to `WindowRecomputation` | compile error in `RecordMeasurement::observe`'s `let Self { … }`, which names every field with no `..` |

The one place the step departs from the house style is `WindowTotals::add`, which reads
`store.positions_observed`, `store.windows` and `store.sweep` by field access while carrying a doc
comment advertising that it destructures. That predates this step — the doc's claim is about
`ComparisonWithTheRun`, which it does destructure — and the step only adds a third direct read, so I
left it. Worth folding into a later step: destructuring `WindowRecomputation` there would make a new
per-store count a compile error rather than a total that silently stays zero.

**The tests' arithmetic checks out.** In the 600-position test, a centre at `p ≤ 250` holds `p + 250`
positions and one at `p ≥ 351` holds `851 − p`; both fall under 300 for 49 centres at each end, so 98
— correct. The 100 centres at 251..=350 are the only ones holding the full 501, so a floor of 501
silences 500 — correct. The thinnest window is the one centred at base 1 holding 251, so a floor of
250 silences none — correct. Each of those doc comments states the arithmetic rather than the
conclusion, which is what they should do.

**Two test doc comments did not, and are fixed.**
`an_arms_silenced_count_is_the_windows_holding_fewer_positions_than_its_floor` opened "Ten covered
positions **inside one window**", while the assertion two lines down is `windows == 10` — there are
ten windows, each holding all ten positions, and the sentence's subject was wrong. It then said
"**This is the claim the whole measurement rests on**", which ranks the test instead of stating what
it pins; the fact survives, the ranking is gone.
`two_sweeps_sum_floor_by_floor_into_an_empty_one` named a constructor rather than the property, and
its doc asserted "its counts are the stores' summed" without the arithmetic a reader can check. It is
now `two_sweeps_sum_floor_by_floor`, and the doc gives the doubled numbers (1,200 windows against
600, 196 silenced at floor 300 against 98, 1,000 at 501 against 500).

The other four names each say what they pin, and the two `should_panic` tests' comments give the
arithmetic that makes the panic the expected one.

**The rest of the prose.** The module doc's "The third question" section leads with what a silenced
window is before the mechanism, which is the right order; I added a one-clause definition of *arm* at
its first use and dropped *bank*, which appeared undefined there and nowhere in the code. The report
is the strongest document in the step: it separates what the floor costs from what it buys and says
plainly that only the first was measured; it names the three stores it had to build and why the
six-accession slice settles nothing on its own; and its closing paragraph on the case it could not
build — a sample whose coverage is patchy over continuous analysed ground — states a real gap
instead of claiming coverage it lacks.

## What I changed

`examples/ng_window_coverage_probe.rs`

- Split the live bank from its answer: `FloorSweep`/`OneFloor` hold accumulators and exist only for
  the pass; new `WindowsUnderEachFloor`/`WindowsUnderOneFloor` carry the counts, `print` and `add`.
  `close` now takes `self` and returns the answer. Deleted `FloorSweep::empty`, the
  `Option` around each arm's accumulator, both `expect`s, and the `accumulator: _` destructure arm.
- `WindowRecomputation` gains `windows_under_each_floor: Option<WindowsUnderEachFloor>`, set by
  `close` out of `sweep`, mirroring the existing `accumulator`/`histogram` pair; both `print` sites
  and `WindowTotals` read it. `WindowTotals` clones the first store's answer instead of summing into
  an empty one, and no longer shadows `store`.
- Added the fifth assertion in `close`: the arm at a floor of 1 silenced nothing.
- Dropped the `sweep-windows` row; tagged the floor rows' extra columns
  (`windows=`, `per-10000=`).
- Module doc: "One walk, several measurements" now covers the third question and says why it rides
  inside the second; "The third question" defines *arm* and drops *bank*.
- Doc comments: `FloorSweep::new` and `silenced_at` given one; `CANDIDATE_FLOORS`'s first sentence
  fixed and its floor-1 claim now backed by the assertion; `the_floors_asked_about` in the present
  tense without the step label; `the_configuration_a_run_uses` says the `..` override is intended;
  the `sweep` field and the `--covered-positions-per-window` flag re-described.
- Tests: `two_sweeps_sum_floor_by_floor_into_an_empty_one` → `two_sweeps_sum_floor_by_floor`, with
  the doubled counts in its doc; the ten-position test's subject corrected and its ranking clause
  dropped.

`src/ng/window_coverage/mod.rs`

- `WindowCoverageConfig::min_window_positions` no longer calls the value provisional.
- `MIN_WINDOW_POSITIONS`: "measured, not inherited" instead of a step label; the "6 in 100"
  comparison replaced by the two shares; one short line rewrapped.

`doc/devel/ng/spec/window_coverage.md`

- §3.3: "no longer soft" pointed at a sentence the same edit deleted — now "measured rather than
  inherited"; the depth and interval comparisons show their numbers; split into two paragraphs.
- §9: "1 in 28" is the smaller of the measured pair — now "between 1 in 28 and 1 in 26"; one
  paragraph rewrapped. *The arithmetic behind those two figures is the other reviewer's to confirm.*

`doc/devel/reports/implementations/ng_window_coverage_d1_2026-09-07.md`

- The "6 in 100" sentences, in both places.
- *arm* used for the missing measurement case → *case*, in both places.
- "Four things are asserted" → five, with the new one listed; *bank* → *sweep*.
- The renamed test and the new type in the "Changes made" and "Tests added" entries.

## Validation

In the container, on this worktree:

- `cargo test --all-features --example ng_window_coverage_probe` — **13 passed**, 0 failed,
  unchanged from the step's own count.
- `cargo check --lib --tests --all-features` — clean.
- `cargo clippy --all-features --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in the library and all predating this branch; none in the example.
- `cargo fmt --check` — `examples/ng_window_coverage_probe.rs` clean; the other files it reports were
  already reported before these edits.
- **Not re-run: the five real-store walks.** No `.psp` exists in this worktree, so the printed-row
  change is unexercised on real data. It is a formatting change with no arithmetic in it, and the
  arithmetic behind the rows is what the six unit tests cover.
