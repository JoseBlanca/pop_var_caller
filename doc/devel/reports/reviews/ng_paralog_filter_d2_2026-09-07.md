# Code Review: ng_paralog_filter_d2

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), two sub-agents in isolated worktrees
**Scope:** step D2 of the hidden-duplication filter plan — the filter at six accessions, the per-pass
timing and the spill's size that make it measurable
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the working tree on branch `ng-paralog-filter` over `c7ad76c3` — the timing
  and spill-size code in [finish.rs](../../../../src/ng/run/paralog_filter/finish.rs),
  [spill_file.rs](../../../../src/ng/run/paralog_filter/spill_file.rs) and their tests, and the
  step's implementation report with the run it reports.
- **Out of scope:** D1's `WhatTheFitCameTo` (reviewed at D1); production; the copied statistics under
  `src/ng/paralog/`; steps A1–C5.
- **Two sub-agents, in their own detached worktrees at `c7ad76c3` with the working tree's `src/`
  diff applied:** one re-measuring every figure in the report independently, one reviewing the code
  and mutation-testing it.

### 2. Verdict

**Request-changes.** Zero Blockers, eleven Major, thirteen Minor.

**Two themes, and they are the same theme.** Five of the report's Majors are claims about *why* a
number is what it is, and every one of them was wrong or unfounded; four of the code's Majors are
places where the code or its documentation said one thing and did another. The counts were right
throughout — the re-measurement reviewer confirmed twenty-two separate quantities exactly, to the
last digit — and what failed was every sentence explaining them.

**And five of nine mutations survived**, including the one that matters most: a build whose calling
pass reports zero shipped green, which would have inflated the scoring share this step exists to
measure.

### 3. Execution status

| command | with this step | at `c7ad76c3` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6643 passed; 0 failed; 15 ignored` | `6639 passed` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9 |
| `cargo fmt --check` | **10** unique files before the fix, 9 after | 9 |

Clippy counted by error kind and `fmt` by unique file, per the failure this plan logged on
2026-09-07. One integration test fails at both, which is `main`'s.

**Mutation totals: 9 run, 4 killed, 5 survived.** All five are closed; see the fix report.

### 4. Open questions and assumptions

1. **The plan's "six-accession slice" is five plants.** `SRS3394712` and `SRS3394712_SRR7279484` are
   two sequencing runs of one biosample, and
   [`rename_dup_samples.sh`](../../../../benchmarks/tomato1/scripts/rename_dup_samples.sh) exists to
   stop a cohort VCF seeing the duplicate column. **A fact about the plan's standing ground**, not
   about this step's code, and the ground cannot move without moving the sha256 every step checks
   against. Recorded in the report and in `PROJECT_STATUS`.
2. **The spill's cost at scale is a spec question.** About 350 GB at 3,000 samples and five million
   records, against about 14 hours of scoring. Spec §5 prices memory and time and not temporary
   disk. Raised for Checkpoint D.

### 5. Top 3 priorities

1. **M-code-3** — a build whose calling pass reports zero passes every test.
2. **M-report-3 / M-report-4** — the report's over-coverage figures divided by a fitted parameter
   and called the result copies, which produced a wrong mechanism.
3. **M-code-1 / M-code-2** — the code's two descriptions of what the timing clock measures were both
   wrong, and the one an operator reads was the wronger.

### 6. Findings

#### Major — the code

**M-code-1: `finish.rs` — the doc for `calling` names three things that happen before the clock starts**

**Confidence:** High.

The field's doc said the span was "the calling pass plus whatever the run does around it: opening
the alignments or the stored files, the merge, and the calls themselves". The reviewer traced
`SpillingSink::beside`'s call sites in both subcommands: before either, the run has already read and
verified the reference, built the segmentation, opened and index-checked every input file, fitted or
read the run parameters, and assembled the header's metadata. What is inside the span is the calling
loop and the spill's own flush — **a better span than the doc claimed**, and the right one, but a
reader believing the doc would read the file-opening cost into it and discount the scoring share.

**M-code-2: the operator-facing line said the clock starts at the first parked record; it starts before the calling loop**

**Confidence:** High.

The printed line said "time from the first record parked" and "the run's own startup before the
first record is not counted here". The clock is taken in `SpillFile::beside`, whose own doc two lines
above says it **creates nothing** — the file appears when the first entry is appended, somewhere
inside the calling loop. So the interval before the first record is *inside* the span, not outside
it, and on a cohort whose first called record is far into a chromosome that interval is large. The
struct field's doc had it right and the operator's line had it wrong.

**M-code-3: the test could not tell the calling clock from a stopwatch started and read on one line — mutation SURVIVED**

**Confidence:** High.

The only guard was `assert!(spent.calling > Duration::ZERO)`. Replacing the span with
`Instant::now().elapsed()` gives **83 nanoseconds**, which passes, and the whole suite went green
with the calling pass reporting 0% and the scoring share inflated by the whole of the calling pass.
That is exactly the number the plan turns into a decision.

**M-code-4: the printed size was not pinned — a wrong unit survived**

**Confidence:** High.

The test asserted the line contained the figure formatted from the same field the line was built
from — self-consistent by construction — and at the fixture's size every plausible divisor rendered
`0.0 MB`. Changing the divisor from 1024² to 1000² passed.

#### Major — the report

**M-report-1: the kept heterozygosity count missed every multiallelic heterozygote**

453 of 13,361 (3.4%) should be **474 of 13,361 (3.55%)**: the count matched the literal string `0/1`
and dropped the 21 heterozygotes spelled `0/2`, `0/4`, `1/2`. "Nine times more heterozygous" is
**8.7 times**.

**M-report-2: "9.7 times the shallowest" is a ratio of fitted parameters, presented as the cohort's depth range**

29.24/3.02 is arithmetically right but both are *fitted* values. The samples' measured median window
depths span **6.6×** (26.08 against 3.96), and reads kept span 6.1×.

**M-report-3: the fits do not sit at the samples' own typical depth, and the run's own dump says so**

Median copy number over all 2,311 records, by sample: 1.38, 0.89, 0.74, 1.31, 0.98, 1.13. The two
shallowest read about 35% over one copy everywhere and `SRS3394713` reads 26% under. The report said
the run "cannot say" whether the puzzle at the top-scoring records was the duplication or the fits;
it can, and the answer was in the file the report already had.

**M-report-4: "the heterozygous sample is at ordinary depth" is wrong, and the whole "does not fit" paragraph with it**

At the 14-record cluster the six window depths normalised to **each sample's own median** are 2.43,
1.22, 1.93, 2.75, 0.83, 0.99. **Three samples are locally over-covered, and the heterozygous one is
among them at 1.93×.** It read as 1.4 copies only because its fitted one-copy value is 35% above its
own median. The report's "the one thing in this run that does not fit" was an artefact of dividing by
a fitted parameter and calling the quotient copies — and the correction makes the model's own story
*hold* at that cluster rather than fail.

**M-report-5: the clustering null was never computed**

It holds, and harder than claimed: drawing 36 of the 2,311 written records at random 200,000 times,
the mean number landing in multi-member clusters is 12.6 and the highest 28, against 34 observed;
P(≥34) is 0 in 200,000. The report asserted "not what a statistical fluke looks like" without the
null that makes it a statement.

**M-report-6: the memory direction does not replicate, and its sign reversed**

Five replicates gave the filter-on runs 15 MB *lower*; eight gave them 10 MB *higher*. The report
offered a mechanism for the direction it happened to measure. Both are 2–3% of a 420 MB footprint;
the conclusion (below the noise) stands and the direction and the mechanism do not.

**M-report-7: the time is extrapolated and the disk is not, and the disk is the one that binds**

Running the report's own extrapolation on the spill gives about **350 GB** of scratch at 3,000
samples and five million records, against about 14 hours of scoring. The report called the spill "the
filter's real resource cost" and then priced only the time.

#### Minor

Code: everything under 51 kB printed `0.0 MB`, including an empty spill (**Mi1**); all four durations
printed `0.00s` while claiming to sum to `0.01s` (**Mi2**); the header assembly fell between two spans
and was charged to nobody while the line said "in all" (**Mi3**); `WhereTheTimeWent` was the one new
struct not destructured where it is read, so a fifth field would be silently dropped from the total —
**mutation SURVIVED** (**Mi4**); three new `rustfmt` violations, taking the file count from 9 to 10
(**Mi5**); `WhereTheTimeWent` missing from `mod.rs`'s re-exports (**Mi6**); a `None` size printed
nothing at all (**Mi7**); `spill_bytes` dropped the "on disk" that made the name say what the value is
(**Mi8**); `scoring` and `scoring_took` in one scope (**Mi9**); nothing pinned that the size survives a
second `finish_writing` (**Mi10**); the test's share parser split the whole report on `(` and found
four numbers only because no other line prints a percentage (**Mi11**).

Report: the 31 "deletions" are 22 deletions and 9 equal-length substitutions (**Mi12**); the
denominator is the 199,672 bases the run called, not 200,000, and the spill ratio depends on which
run's output it is against (**Mi13**).

### 7. What is right, and was checked

- **The size is read after the flush and after the handle is closed** — verified by the mutation that
  moves it before the flush, which fails because the fixture's records are still inside the 64 KiB
  buffer.
- **Nothing runs per record**: `Instant::now()` four times a run, `fs::metadata` once.
- **The clocks are monotonic and saturating**, so no duration can come back negative on a wall-clock
  step.
- **The filter-off path builds no `SpillFile`**, starts no clock, and never calls the report.
- **`FilteredRun` and `SpillFile` keep their exhaustive destructuring**, both updated for the new
  fields.
- **Pass two is O(records × samples)** — both grids are fixed constants — which is what the
  extrapolation assumes.
- **Twenty-two report quantities reproduced exactly**, including the six fitted depths and sigmas, all
  seven ratio percentiles, the six clusters and their spans, the flagged heterozygote count, the
  absent-window rate, the record kinds, and the six-accession oracle's sha256.
