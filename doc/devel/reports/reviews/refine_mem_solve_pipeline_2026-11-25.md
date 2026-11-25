# Code Review: refine_mem_solve_pipeline
**Date:** 2026-11-25
**Reviewer:** rust-code-review skill (orchestrator, single reviewer — see Scope)
**Scope:** correctness of the refinement changes merged to main between 01eb7714 and 1ee9488c
**Status:** Approve-with-changes

---

## 1. Scope

**What was reviewed.** The diff `01eb7714..1ee9488c` on `main`: three branches of work on
refinement, the neighbour-evidence pass that re-genotypes each sample from its neighbours'
haplotypes.

- **Memory cuts** (`refine_mem_build_2026-11-24.md`): the walks' per-thread arrays (paths handed to
  the phase as they are made, a path's states in two bytes, `Uniforms::skip`, arrays at exact
  length); a skipped record held in one form at a time, packed or whole, and packed again when no
  window walks it; one set of genotype sums per sample in a phased core; the budget counting the
  larger of the haplotype store and the written results; written lines counted at nothing.
- **Solve speed-ups** (`refine_solve_build_2026-11-24.md`): the parameter fit running the forward
  pass only, through a forward loop shared with the full walk; emissions from an allele-pair table
  where every template carries one allele for certain; genotype pairs indexed by allele count;
  template choice reading both copies of a sample in one pass; kept-out sites marked once a round.
- **Pipeline** (`refine_pipeline_build_2026-11-24.md`): arriving lines read into sites on the pool
  in batches of up to 1,024 lines or 8 MB of text, taken in arrival order; a line rewritten with
  borrowed fields. The bgzf-writer thread (c7082615) and its revert (1ee9488c).

**Reviewed against:** 1ee9488c (main).

**In-scope files:** `src/neighbour_evidence/{model,refine,stage,windows,templates,fit,random,line}.rs`
and their tests; `src/vcf/writer.rs` (revert only); the three implementation reports,
`PROJECT_STATUS.md` and `doc/devel/ng/spec/neighbour_evidence_pass.md` for claims.

**Out of scope:** the examples (`refine_mem_probe`, `solve_bench`, `refine_pipeline_probe`) —
measurement tools, not shipped code.

**Categories.** The caller asked for a correctness review with a named checklist, so this was run
as one reviewer working through that checklist, with probes in a private worktree, rather than the
skill's per-category fan-out. Style, naming and module-structure categories were not run.

## 2. Verdict

**Approve-with-changes.** Every path checked for output change keeps the old bytes: the random
streams, the emission table, the genotype-pair index, the shared forward loop, the template choice,
the read-in-batches queue and the borrowed rewrite. One regression was found and reproduced: a line
that zstd makes *longer* (any line under about 60 bytes) now panics with an arithmetic overflow in
test and `soak` builds. Release builds wrap and happen to compute the right total. The budget's new
"store or results, not both" rule leaves an uncounted result buffer uncovered, which matters most
below 256 samples. One implementation report still describes the reverted bgzf thread as live.

## 3. Execution status

- `cargo test --lib -- neighbour_evidence vcf::writer read_phase` at 1ee9488c, in a detached
  worktree, with two probe tests added: `test result: FAILED. 208 passed; 1 failed` — the one
  failure is the probe behind M1, every existing test passed.
- The same two probe tests at 01eb7714: `test result: ok. 2 passed; 0 failed`.
- Probe patch kept at `tmp/review_2026-11-25_refine/probe_tests.patch`; the worktree was removed.
- `zstd -1` on a 51-byte one-sample line gives a 60-byte frame (measured with the zstd CLI; the
  `bulk::Compressor` also writes the content size, so its overhead is the same order).
- Not run: `cargo fmt --check`, `cargo clippy`, `cargo doc`, `cargo audit`, the full test suite
  (out of the correctness scope asked for).
- Findings labelled "Needs verification": 1 (M2's size).

## 4. Open questions and assumptions

1. Is the uncounted 256-sample result chunk (M2) meant to be covered by the budget's headroom, as
   the walks' arrays are? If so, the headroom was measured with the store still counted beside the
   results, and needs measuring again. (M2)
2. Should `PROJECT_STATUS.md`'s "Current focus" block be updated by the author (it still says
   "Not merged" and "R2 … owner to decide")? This review did not edit it. (Mi1)

## 5. Top 3 priorities

1. **M1** — restore signed arithmetic in `compress_held_lines`; a short line panics test and soak
   builds.
2. **M2** — count, or shrink, the result chunk now that the store's slack no longer covers it.
3. **Mi1** — mark R2 as reverted in the pipeline report and the status file.

## 6. Findings

### Major

#### M1: src/neighbour_evidence/stage.rs:815-821 — a line zstd lengthens overflows `saved`

- [stage.rs](../../../../src/neighbour_evidence/stage.rs#L815-L821) — **Major**
- **Confidence:** High (verified by a failing probe test; passes on 01eb7714).
- **Problem:** e848066e changed the running saving from `i64` to `usize`:
  ```rust
  let mut saved = 0;
  for (index, held) in self.held.iter().enumerate().skip(from) {
      if let HeldText::Compressed { bytes, length } = &held.text {
          saved += *length - bytes.len();
      }
      self.held_bytes[index + 1] -= saved as u64;
  }
  ```
  zstd adds a frame header and block header (about 9 bytes) to input it cannot shrink, so a short
  line comes back longer than it went in and `*length - bytes.len()` goes below zero. The probe
  `review_probe_short_filtered_lines` pushes 3,000 one-sample lines of 44 bytes
  (`chr1\t…\tLowQual\tAF=1\tGT\t1/1`) and panics: `panicked at src/neighbour_evidence/stage.rs:818:26:
  attempt to subtract with overflow`. The same test passes on 01eb7714.
- **Why it matters:** the `test` profile and the `soak` profile (`overflow-checks = true`,
  Cargo.toml) abort on any run that holds a short line — one-sample cohorts and filtered lines with
  few columns are the likely cases. Release builds wrap; because every step wraps modulo 2⁶⁴ the
  final total is still right, so release output is correct, but only by that accident.
- **Suggested fix:** go back to a signed saving, or adjust each total by the line's own change:
  ```rust
  let mut saved: i64 = 0;
  for (index, held) in self.held.iter().enumerate().skip(from) {
      if let HeldText::Compressed { bytes, length } = &held.text {
          saved += *length as i64 - bytes.len() as i64;
      }
      let total = &mut self.held_bytes[index + 1];
      *total = (*total as i64 - saved) as u64;
  }
  ```
  Add the probe as a regression test (section 8).

#### M2: src/neighbour_evidence/refine.rs:1159-1205 — the result chunk is no longer covered by the store's count

- [refine.rs](../../../../src/neighbour_evidence/refine.rs#L1159-L1205),
  [windows.rs](../../../../src/neighbour_evidence/windows.rs#L307-L321) — **Major**
- **Confidence:** Medium. The mechanism is read from the code; the byte sizes are computed from the
  struct layouts, not measured. **Needs verification.**
- **Problem:** the claim that the store and the results are never held together is true:
  `drop(store)` precedes `RefinedRecord::unchanged` (refine.rs:1146-1155). But the results phase
  holds more than the results the count charges. For each chunk of up to `SOLUTIONS_AT_ONCE` = 256
  samples it builds a `SampleSolution` per sample, and that buffer is not counted. Per sample and
  core record (two alleles), a solution holds a `SolvedSite` (24 bytes), an `Option<WrittenPhase>`
  (32) and three `f32` posteriors (12): about 68 bytes, against 29 bytes for the counted written
  result. Each sample's state is freed as its solution is made: 16 bytes of sums per core record
  and 17 bytes of phase per window record. Where the window is about the size of its core — a
  contig shorter than a window, where margins are clipped — that frees about 33 bytes and adds
  about 68. Before this change, the store (20 bytes a sample a record) was counted on top of the
  results and absorbed most of the difference. Now the count takes the larger of the two, so
  nothing absorbs it.
- **Why it matters:** at 256 samples or fewer every sample's solution is held at once, together
  with every result. By the arithmetic above, that is about 29 + 35 = 64 bytes a sample a core record
  in the results phase, against 29 counted. Above 256 samples the gap is a fraction of the window
  (256 / samples). The budget is the run's memory guarantee, and this is the corner where it binds
  least often but is now furthest off.
- **Suggested fix:** one of the following.
  1. Count the chunk: add `min(samples, SOLUTIONS_AT_ONCE) × core records × solution bytes` to
     `held` in `plan_window`. It is a function of the records alone, so windows stay independent
     of the thread count.
  2. Shrink the solution so it fits inside the written result: `record: u32`,
     `posterior_start: u32`, and the phase as `(u8, u8, u64)` or an index into a side table.
  3. At least measure the results phase's peak against the count, using `refine_mem_probe` under
     dhat at 20 and 200 samples with a contig shorter than one window, and record the ratio next
     to the 1.09 figure the memory report quotes.

### Minor

#### Mi1: doc/devel/reports/implementations/refine_pipeline_build_2026-11-24.md:21-29, 42, 218-223 — R2 still described as live after the revert

- [refine_pipeline_build_2026-11-24.md](../implementations/refine_pipeline_build_2026-11-24.md#L21-L29) — **Minor**
- **Confidence:** High.
- **Problem:** the revert itself is clean: `git diff c7082615^ 1ee9488c -- src/vcf` is empty, so
  `writer.rs` and its tests are back to what R1 left. But the report still lists R2 in the
  changes table (line 42), says "R2 raises peak resident memory … **The owner should decide
  whether to keep R2**" (lines 25-29), and keeps "Keep R2 or revert it?" as an open decision (line
  218). `PROJECT_STATUS.md:29` says "R2 … owner to decide. Not merged", and lines 33 and 38 say
  "Not merged" for the two other branches, all three of which are on main.
- **Why it matters:** a reader of the report would think a bgzf hand-off thread is in the writer,
  and would look for a memory effect that no longer exists.
- **Suggested fix:** add a dated line at the top of the report: "R2 reverted on main in 1ee9488c
  (2026-11-25), as recommended; the output and memory figures for R1 and R3 stand." Mark decision 1
  as settled, and in `PROJECT_STATUS.md` replace "owner to decide. Not merged" with "R2 reverted
  (1ee9488c); merged to main", and "Not merged" with "merged to main" on lines 33 and 38.

#### Mi2: src/neighbour_evidence/stage.rs:689-709 — what the read-in-batches queue holds is not bounded by text alone

- [stage.rs](../../../../src/neighbour_evidence/stage.rs#L689-L709),
  [stage.rs](../../../../src/neighbour_evidence/stage.rs#L333-L339) — **Minor**
- **Confidence:** Medium (read from the code, not measured).
- **Problem:** `BYTES_READ_TOGETHER` bounds the waiting lines' *text* at 8 MB, but three other
  things wait with it, uncounted. (a) Each `ArrivedLine` carries its unpacked `RecordFragments`,
  which grow with read depth, not with line length. (b) `take_arrived_lines` turns the whole batch
  into `SiteRecord`s at once, before taking any of them. (c) A window completed by an early line is
  solved inside `take_line` while the rest of the batch's parsed records stay alive in `read`. So
  "held only until each line is taken, a moment later" (stage.rs:336) can mean the full length of
  a window solve, with up to 1,023 parsed records beside it.
- **Why it matters:** the amount is bounded — roughly the parsed size of 8 MB of text plus 1,024
  lines' fragments, tens of MB at worst. It does not grow without limit, but it sits outside the
  budget and outside the doc comment's description of it.
- **Suggested fix:** state the bound truthfully in the doc comment ("the parsed records of one
  batch, held through any window an early line completes"). Optionally count each line's
  fragment bytes in `arrived_bytes`.

#### Mi3: src/neighbour_evidence/stage.rs:665-667 — a stop drops the rest of a parsed batch silently

- [stage.rs](../../../../src/neighbour_evidence/stage.rs#L665-L667) — **Minor**
- **Confidence:** High (read from the code). There is no output defect today: the only caller,
  `VcfWriter::write_line`, turns any `StageStop` into a fatal `VcfWriteError`.
- **Problem:** when `take_line` returns an error part-way through a batch, the remaining lines in
  `read` are dropped, and so are lines already `mem::take`n from `arrived`. The doc says "A stopped
  stage is not fed again", but nothing enforces it. A later `push_writing` or `finish_writing`
  carries on with lines missing from the middle of the stream.
- **Why it matters:** a future caller that recovers from an over-budget stop, for example by
  writing the remaining lines unrefined, would lose up to 1,023 lines without an error.
- **Suggested fix:** add a `stopped: bool` that `push_writing` and `finish_writing` check first,
  returning the stop again (or panicking with a message naming the misuse).

## 7. Out of scope observations

- `src/neighbour_evidence/line.rs` — `rewrite_line` writes `PS` over a sample column's extra field
  when the column has more fields than `FORMAT`. The pipeline report already notes this; main and
  the branch behave alike.

## 8. Missing tests to add now

- `compress_held_lines_keeps_totals_when_a_line_grows` (stage): one-sample lines under 60 bytes,
  3,000 of them, through `push` and `finish` under the `test` profile. It catches M1. The probe in
  `tmp/review_2026-11-25_refine/probe_tests.patch` (`review_probe_short_filtered_lines`) is the
  test as written.
- `results_phase_peak_is_within_the_count_at_small_cohorts` (refine, measured under dhat in
  `refine_mem_probe` rather than as a unit test): 20 samples, a contig shorter than one window,
  peak requested bytes during the results phase against `plan.bytes`. It catches M2.
- `a_stopped_stage_refuses_further_lines` (stage): force an over-budget stop mid-batch, then call
  `push_writing` again and expect an error. It catches Mi3 once a guard exists.

## 9. What's good

- `Uniforms::skip` keeps every random stream aligned without drawing ahead. The walk advances the
  parent stream by exactly `sites.len()` per path, the backward sweep draws exactly one uniform per
  site per path, and `skipping_lands_where_drawing_does` pins the arithmetic
  ([random.rs](../../../../src/neighbour_evidence/random.rs#L67-L74)).
- `definite_emissions_into` recomputes each table entry with the general path's own operations in
  its own order, including the `allele == 0` assign-then-add, so the bits match by construction
  rather than by tolerance ([model.rs](../../../../src/neighbour_evidence/model.rs#L429-L533)).
- One `forward_pass` serves both the full walk and the fit's likelihood-only walk, so the fit cannot
  drift from the walk ([model.rs](../../../../src/neighbour_evidence/model.rs#L1126-L1158)).
- `check_the_kept_totals` now asserts that each site is held in exactly one form and at the
  packed length the budget counted. That turns the repack lifecycle into a checked invariant on
  every test stage ([stage.rs](../../../../src/neighbour_evidence/stage.rs#L1299-L1311)).

## 10. Commands to re-verify

- `git apply tmp/review_2026-11-25_refine/probe_tests.patch` in a scratch worktree at main, then
  `cargo test --lib review_probe` — expect `review_probe_short_filtered_lines` to fail until M1 is
  fixed.
- `cargo test --lib -- neighbour_evidence vcf::writer read_phase` — 208 existing tests pass at
  1ee9488c.
- `git diff c7082615^ 1ee9488c -- src/vcf` — empty, confirming the revert is clean.

### Checked and found unchanged

These are the checklist items from the request that were traced and found to keep the old output.

- **Random streams.** The round's own path and each extra path clone the stream and then skip the
  parent by `sites.len()`. The backward sweep calls `draw` at every site of every block, one uniform
  each, from the last site back, which is the order `stream[sites.len() - 1 - index]` read before.
  Multi-allelic records take one uniform, as before. Early exits do not exist inside the sweep.
- **Extra paths handed on as made.** `haplotypes_along` still runs own path first, then extras in
  draw order, with the same `reads` and `joints`. `on_path` records each path in that order, and
  with `extra_paths = 0` the own path is still recorded.
- **Two-byte states.** `K ≤ 255` is asserted only when a walk draws. `MAX_TEMPLATES` is 16, so
  `K + K²` is far below `u16::MAX`.
- **Emission table.** Used only where every template's emission shares are a single 1. Any other
  share (including values above 1, or an all-zero template) falls back to the general path.
- **Genotype pairs by count.** Same `GenotypeTable` order. `known(0)` now returns `None` where it
  returned an empty slice, but no record has zero alleles.
- **One set of sums.** `hidden` is fixed before the rounds, so a site sums the same distribution in
  every round. Hidden records are drawn only from informative core records, so the empty
  `posterior` placeholder is always replaced before it is read.
- **Repack lifecycle.** The solve reads records immutably (`self.sites.window`), so repacking
  stores the numbers the first packing stored. Whole sites held outside the window are counted
  whole (`resident_bytes_of` reads `sites.whole`). A site walked by two overlapping windows stays
  whole; one the next window skips is packed again in that window's `hold_whole`, before its solve.
- **Parse queue.** A contig change inside a batch goes through `take_line` → `finish_contig`, as
  before. End of input goes through `finish_writing` → `take_arrived_lines`. Parsing is per line
  and does not depend on the thread count. `lines_read_one_at_a_time_or_in_batches_write_the_same_bytes`
  passes.
- **Borrowed rewrite.** Columns are joined in the same order with the same separators. Short
  columns are padded with `.`, long columns are kept, and `:PS` is appended to `FORMAT` exactly
  where it was before.
