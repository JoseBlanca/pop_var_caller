# window coverage — C3: the look-ahead, and the oracle that can see it

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone C, step C3
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.3, §6 trap 3, §10
**Review:** [ng_window_coverage_c3_2026-09-06.md](../reviews/ng_window_coverage_c3_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [B2](ng_window_coverage_b2_2026-09-06.md), [C2](ng_window_coverage_c2_2026-09-06.md)

## The answer

**Each cover now draws half a window past the region it was asked for, and that is what stops a
region's last centres being absent.** Measured on the six-accession tomato slice, over the 26,754
sample-loci a psp-mode run builds there:

| | run has no window, the walk does | run and walk agree bit for bit | disagree |
|---|---|---|---|
| without the look-ahead | **571** | 25,674 | 0 |
| with it | **12** | 26,233 | 0 |

The **12** that remain are two positions — `SL4.0ch01:13906474` and `SL4.0ch01:13906595`, the same
two in all six samples — and they are the last two loci of the last analysed interval, where the
sample has no later record for the accumulator to close the centre on. Nothing a cover can do
reaches them; `finish` will, at plan step C5.

**No VCF byte moved**: 2,311 records, sha256
`84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on both routes, before and
after this step, as at every step since A1.

## What was wrong, and how it was invisible

A window centred at `p` is finalised only once that sample's stream has passed `p + 250` (spec
§3.3). A cover drew each sample to its region's last base plus whatever chained past it, so the
centres in the region's last stretch were still open when the builder for that region ran: it read
`None` and the sample was absent there. **The measurement changes no VCF byte, so a run in which
every such centre was missing writes exactly what a correct run writes.** Nothing in the output,
and no existing test, could tell the two apart.

Without the look-ahead the absences are not scattered: on `SRS3394606` they arrive in clusters of
six to twelve loci spanning about 200 bases of reference, one cluster about every 16 kb — 96 of
that sample's 4,459 loci. With the look-ahead the same sample has 2, and they are the end of the
run's ground.

## What landed

### The look-ahead itself

`ObservationCache::where_the_chain_starts` seeds both covers' chain at the region's last base plus
`window_coverage::WINDOW_BP / 2` — **derived from the constant, not retyped** — clamped to the
contig's length. Both `cover` and `cover_in_parallel` call it, so the two fixpoints cannot start
from different places.

**The clamp is this step's own case and belongs to it.** C2 made a cover read the reference over
each contig it added records on, and the ground it reads runs to the chain; a chain pushed 250
bases past a region at the end of a contig asks the reference for bases past the contig's last one,
which is refused rather than truncated. `MergeReference` gains `ContigTable` for the lengths —
every reference the merge is given already carries them, and the fixture that wraps one now
forwards it.

**Drawing further is not the same as finalising, and the difference is why this is enough.** What
closes a centre is a *position* beyond it. The draw stops at the first record starting past the
chain and **holds it**, so a sample with any record after `p + 250` has had it observed by this
cover. What is left is a sample with no such record at all — the last half-window of its records —
which is exactly the 12 above.

### The recorder, and the whole-store recomputation

Two pieces, because the failure is invisible from inside the run:

- **[`src/ng/run/cohort_merge/recorded_windows.rs`](../../../../src/ng/run/cohort_merge/recorded_windows.rs)**
  — when `NG_WINDOW_COVERAGE_FILE` names a file, every sample's window at every built locus is
  written to it, as **bit patterns** (an absent window is a pair of `NaN`s, and a printed `NaN`
  cannot be told from another one). A sample with no window is written down rather than omitted,
  because "the run had nothing here" and "the run never reached this locus" are the two halves of
  what is being looked for. Off, it is one already-resolved `OnceLock` read per built locus.
  **Both halves of the row format live in this one file**, the writer beside the reader, because
  the reader is compiled into another crate target: a column reordered on one side alone still
  parses on the other, and the comparison would then report every locus as one the run had no
  window for — which is what a *pass* looks like in this measurement.
- **[`examples/ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)**
  gains `WindowRecomputation`: it walks a store end to end — one pass, no covers, no eviction —
  feeding the same accumulator the same positions and depths against the same reference, and
  compares its windows against the run's by bit pattern. **It differs from the run in how it reads
  the store and in nothing else**, which is what makes a difference attributable to the walk.

The comparison refuses an empty verdict: a dump naming no locus, or one in which nothing was ever
compared against a finalised window, fails rather than reporting zero disagreements.

`window_coverage::depth` became `pub` so the recomputation can call the very rule the cache calls.
An oracle that reimplemented the depth rule would be checking its own copy of it, where the only
thing that differs between the two sides is the walk.

## Assumptions and deviations, all minor and all recorded

1. **The plan says "a test-only hook"; this is an environment-variable-gated recorder in live
   code.** A `#[cfg(test)]` hook cannot run over the six-accession slice, which is where the plan
   asks for the comparison. The pattern is the repository's own (`NG_TRACT_DUMP`,
   `PVC_MINTED_ERROR_CENSUS`, `PVC_PARITY_*`).
2. **The probe takes `--reference` and `--windows-from-the-run` as flags**, and asks only B2's
   question without them — so B2's command still runs unchanged.
3. **The comparison is run through a script under `tmp/`**, not through
   `scripts/ng_mode_equivalence_oracle.sh`. `dev.sh` forwards only `CARGO_*` and `RUST*`, so the
   recorder's variable has to be set inside the container. Folding it into the standing oracle
   would make every step's oracle run twice as long for a check only this step and C5 need.
4. **Six existing fixtures moved their coordinates.** Each was written when a cover stopped at its
   region's end; a fixture packed inside 250 bases is now drawn whole by the first cover and no
   longer exercises what its name says. Their claims are unchanged and their prose says why the
   positions are far apart.
5. **One fixture in `callers.rs` gained a record on the next contig.** The test stages a refused
   record *before* a source failure, and the ten-base building regions it used to do that with are
   not enough now that a cover draws 250 bases further: the fixture's `chr1` is a hundred bases, so
   the first cover reached the end of the source. A record on the next contig is past every reach
   on this one, so it stops the draw there.

## Changes made

- **[`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)** —
  `where_the_chain_reaches_before_any_sweep`; `MergeReference` gains `ContigTable`; the cover leaves
  the buffer on the region's own contig; four new tests and five fixtures moved.
- **[`recorded_windows.rs`](../../../../src/ng/run/cohort_merge/recorded_windows.rs)** — new, with
  three tests.
- **[`window_coverage/mod.rs`](../../../../src/ng/window_coverage/mod.rs)** — `depth` is `pub`.
- **[`build.rs`](../../../../src/ng/run/cohort_merge/build.rs)** — the recorder's one call site,
  in the arm that keeps a built locus.
- **[`cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)** — the release-counting
  fixture reference forwards its contig table.
- **[`callers.rs`](../../../../src/ng/run/callers.rs)** — one fixture, above, plus the segmentation
  over both fixture contigs the review's finding needed.
- **[`ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** — the
  recomputation, the dump reader, the comparison and its verdict.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched.

## Tests added

Eight, six of which the reviews asked for.

| test | what it pins |
|---|---|
| `a_cover_reads_half_a_window_past_its_region_with_nothing_chaining_there` | the look-ahead: the ground reaches `region.end + WINDOW_BP / 2` where nothing chains there, and no further |
| `the_look_ahead_stops_at_the_contigs_last_base` | the clamp: a region within half a window of the fixture reference's 2,000-base contig covers rather than failing, and the ground stops at base 2,000 |
| `a_region_past_its_contigs_end_is_not_clamped_back_behind_itself` | the clamp's floor, which one review's mutation had found unguarded: a region overrunning its contig keeps its own ground rather than having the chain pulled back behind the builder it was made for |
| `a_record_held_on_the_next_contig_does_not_take_the_buffer_with_it` | the cover ends holding the region's ground, not the ground of the contig its overshoot draw reached |
| `recorded_windows::a_locus_writes_one_row_a_sample_naming_the_sample_and_the_locus` | the three fields the comparison keys on — the sample index, the position, and that a sample with no window is written down rather than omitted |
| `recorded_windows::a_row_survives_being_written_and_read_back` | the row format's two halves, which are compiled into two different crate targets — an ordinary window, one the floor silenced, and a sample the run had no window for |
| `recorded_windows::a_row_that_will_not_parse_says_which_column_failed` | a malformed row names its column |
| `recorded_windows::with_no_file_named_nothing_is_written_and_nothing_is_read` | off, the recorder reads nothing from the window it is handed |

And one in the probe, `absent_windows_agree_and_a_last_bit_apart_is_a_disagreement`: the comparison
is by bit pattern, so two absent windows agree and two values a last bit apart do not.

## What the review changed

Two agents, one on reliability and errors and the numbers this report claims, one on naming, idiom,
module structure and refactor safety; the full report is
[beside this one](../reviews/ng_window_coverage_c3_2026-09-06.md). Five Major findings and nine
Minor. The three that were defects rather than clarity:

- **A cover that draws a sample onto the next contig ended holding *that* contig's bases**, so the
  accessor a builder reads answers `None` at every position of the region it owns. **Both agents
  found it independently, and C2's own re-run correctness review found it a third time.** It is
  latent — that accessor has no caller outside tests until C4 — but the look-ahead turns it from a
  corner into the ordinary case at every contig boundary, which is why it is fixed here. The cover
  reads the region's ground again at the end, observing nothing.
- **The `callers.rs` fixture this step changed stopped testing its own name.** Its subject is that
  a refused record outranks a source failure behind it; the record added to stage the two stops the
  draw for the whole of `chr1`, and the run's analysed ground is `chr1` alone, so the source failure
  was never drawn. Swapping the two arms left all 515 `ng::run` tests green. The test now runs over
  a segmentation covering both fixture contigs, and the swap fails it again. **The general fact,
  which will bite the next fixture:** `chr1` is a hundred bases and half a window is 250, so every
  cover on it reaches that contig's end — a fixture needing a cover to stop short now needs a contig
  longer than 250 bases, not a narrower building region.
- **Every field the recorder writes, and the whole of the probe's comparison, had no test.** Four
  deliberate defects survived the suite: the sample index off by one, the position off by one, the
  absent rows omitted, and the bit comparison replaced by float equality — the last of which would
  report a *correct* run as disagreeing at every window the position floor silenced. All four are
  now caught.

The design half's findings, applied:

- **The look-ahead's helper was named for the wrong end of the value it returns.** It gives the
  chain's *reach*, and both call sites read `let mut chain_reach = self.where_the_chain_starts(…)`
  — the name and the variable it fills disagreeing in one line. It is
  `where_the_chain_reaches_before_any_sweep` now.
- **`cover`'s own doc comment said the reach "starts at `region`'s last base"**, which this step
  had just made false, and two more paragraphs of it were left incomplete: what bounds how far a
  cover goes, and what a cover costs. All three now carry the look-ahead.
- **The row format was written in the library and re-derived by column index in the example.**
  Swap two columns on the writing side and every field still parses, and the comparison then
  reports each locus as one the run had no window for — *which is the shape of a pass in this
  measurement*. Both halves live in one file now, with a round-trip test over an ordinary window,
  a silenced one and a missing one.
- **The recorder was in the wrong module.** It was live `pub` code beside `production_parity.rs`,
  which is `#[cfg(test)]`; it named `WindowedCohort`; its one caller is the cohort merge; and it
  performed no parity check — the comparison is the probe's. It is
  `cohort_merge/recorded_windows.rs`, which also takes the dependency out of the
  `window_coverage` ↔ `cohort_merge` cycle rather than adding a third file to it.
- **The probe told the reader that a locus written twice is routine.** It is not: a builder skips
  a locus it does not own rather than building it, so a duplicate means the ownership rule broke
  or two runs wrote into one file. A reader told duplicates are normal would explain away the one
  thing a duplicate can reveal.
- **The window comparison's totals were copied field by field** where the file's other totalling
  function destructures precisely so a count added later cannot be left out. It destructures now.
- Smaller: a flag whose doc described an `Option` that was not there; two allocations per record
  in the recomputation, including a copy of every record's reference bases; a per-record linear
  scan over the reference's contig list; `(u32, u64)` map keys where `GenomePosition` fits and
  already orders by genome order; and a usage block that never said where the recorded file comes
  from.

**One recommendation deferred to C4**, and it is right: the recorder is called where a locus is
built, which is correct today because the loci built are the loci written — but from C4 the pair
travels on the observation and the record sink is where "at every written record" is decided. The
call moves there then, and the comment says so.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,366 passed, 0 failed, 15 ignored** (6,358 before this
  step). The eight above are the whole increase.
- `cargo test --all-features --example ng_window_coverage_probe` — 6 passed, 0 failed.
- `cargo clippy --lib --all-features --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — every file this step touched is at the hunk count it had before it
  (`observation_cache.rs` 4, `serial.rs` 4, `build.rs` 1, `callers.rs` 0, `cohort_merge/mod.rs` 0
  of its own); `window_coverage/mod.rs`, `recorded_windows.rs` and the probe are clean.
- **The standing oracle either side of this commit**: 2,311 records, sha256 `84ad19c2…0590d` on
  both routes.
- **The whole-store recomputation**: 1,172,242 covered positions over six stores, the same number
  of windows finalised, **0 disagreements**, 26,233 loci agreeing bit for bit, 12 where the run has
  no window and the walk does.

## Tradeoffs and follow-ups

- **What the look-ahead costs is half a window more of held records per sample.** On this slice a
  sample has 198,710 covered positions over 200 kb of reference — 0.99 records a base — so half a
  window is about **248 more records held per sample**, at a mean 14.4 reads compared with the
  reference at each. Plan step D3 prices the milestone's memory.
- **Twelve loci a run still reads as absent**, at the end of the run's own ground. `finish` closes
  them and C5 is where `finish` is called; whether the pair at those loci can reach the record they
  belong to is C5's question, not this step's.
- **The recorder's cost when on** is one `Mutex` and one flush per built locus. It is a measuring
  tool and no run that is not measuring pays for it.
- **`ng_cohort_merge_real_cost` does not compile**, and has not since before this branch (its
  `MemberRange.observations` field does not exist on `main` either); C1 added two further errors to
  it by widening `ObservationCache::over` and `cover`. Fixing it needs a decision about a module
  another branch is working in, so it is raised rather than taken.
