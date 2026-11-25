# Refinement's pipeline around the solve: three speed-ups built (2026-11-24)

> **2026-11-25: R2 is reverted.** The bgzf writer's thread of its own (c7082615) was reverted on
> main in 1ee9488c, as this report recommends: it saved no time measurable at 2 or 6 threads and
> raised tomato2's peak resident memory by 30 to 70 MB. The writer is as R1 left it. R1 and R3 are
> on main, and their figures below stand. Where this report describes R2 below, it describes
> what was built and measured, not what main holds.

Branch `refine-pipeline`, from main bf887b45, which already holds the window-memory cuts
([refine_mem_build_2026-11-24.md](refine_mem_build_2026-11-24.md)) and the solve speed-ups
([refine_solve_build_2026-11-24.md](refine_solve_build_2026-11-24.md)). It builds items R1, R2 and
R3 of the review `doc/devel/reports/reviews/refine_pipeline_speed_opportunities_2026-11-24.md` —
a review on branch `pipeline-hunt` (head 5bbeac5c), not on this one. The review's prototypes were
written on 01eb7714, before both merged branches; each was written again here on the current code,
one commit each. R4 (skipping the window planner) and everything the review rejected are not built.

**What a reader should take away.**

- **R1 removes most of the serial work outside the solve.** The thread every line passes through
  spends 0.67 s of its own time on 300,000 lines of 10 drawn samples, and 0.15 s after R1; at 100
  samples 0.52 s becomes 0.08 s, at 500 samples 1.76 s becomes 0.32 s. At 6 threads that is 6.6% of
  the wall time at 10 samples (5.94 to 5.55 s) and 4.8% at 500 (33.1 to 31.5 s).
- **On tomato2 the three together save 1 s of a 97-s writing pass at 2 threads and 2 s of 43 s at
  6.** The solve is most of that pass, so a serial step of a second or two is all there was to win
  at these thread counts. The review's projection — a third of the writing pass at 48 threads — was
  not measured here: the laptop allows six.
- **R2 and R3 save nothing measurable on this laptop.** R2, the thread in front of the bgzf
  writer, did not change wall time beyond the run-to-run spread on any cohort; R3 cuts the
  allocations by 29% (10.9 to 7.8 million blocks at 10 samples) but the process's CPU by under 1%.
- **The output is byte for byte main's after every commit**: 12 drawn cohorts through the real
  `VcfWriter` (10, 100 and 500 samples, 3 and 30 reads a position, thinned margins and not), 12
  multi-allelic cohorts through the stage alone, and tomato2's VCF body.
- **R2 raises peak resident memory on tomato2 by 30 to 70 MB at 6 threads** (483 to 508 MB on
  main and with R1, 513 to 562 MB with R2), more than the 12 MB its buffers hold. The bytes the
  program asks of the allocator do not rise beyond the buffers, so the rest is most likely the
  allocator keeping memory for the threads that now work at once. **The owner should decide
  whether to keep R2** (see *Decisions for the owner*). It was reverted (note at the top).

Measured on a laptop (macOS, arm64) that was quiet for the timed runs (one-minute load 2.4 to 5);
at most six threads.

## What changed

| Review item | Commit | What |
|---|---|---|
| — | 25097147, f657c5c6 | `examples/refine_pipeline_probe.rs`, the review's timing tool, ported without its section timers; mimalloc as `pvc` uses it, and dhat when built with `dhat-heap`. |
| R1 | 4bd13d61 | Arriving lines read into sites on the pool, a batch at a time. |
| R2 | c7082615 | The multithreaded bgzf writer fed by a thread of its own. **Reverted in 1ee9488c (2026-11-25).** |
| R3 | 19f46c33 | A line rewritten with its fields borrowed, into one buffer. |

### R1 — arriving lines read into sites on the pool

Each line arriving at the refinement stage is read into a *site*: its `PL`, `GP` and `AF` parsed
into the numbers the solve works on. That was done on the writing thread — the one thread every
line passes through — one line at a time, while the pool of `--threads` workers waited.

Now lines wait in the stage (`NeighbourEvidenceStage::push_writing`) until 1,024 have arrived or
their text reaches 8 MB, are read into sites together on the pool (`take_arrived_lines`), and are
then taken one at a time in arrival order by the code that took each on arrival before
(`take_line`). Everything a window plan depends on — the running byte totals, the margin thinner's
token bucket, the planner's call after each line — sees the same lines in the same order, so the
output cannot change. `finish_writing` takes the lines still waiting; a contig change is handled
inside `take_line`, as before.

**What changes for a run that stops.** A window over the budget, or a write the file refuses, is
found when the line that reveals it is *taken*, which can be up to 1,023 lines after it arrived. On
tomato2 with `--refine-window-bytes 1000000` the error is the same window — "from contig 1
position 3406993, holds 2963 records of 59 samples … Raise --refine-window-bytes to at least
16576156" — but the record the writer names as the one it could not write moved from position
4,407,353 on main to 4,420,323. Nothing else in the message changed, and the run exits 1 as before
without hanging. That record is only the one being handed in when the stop was noticed; it is not
in the window.

**Memory.** At most 8 MB of line text wait, with the sites read from them until each is taken a
moment later. At 2,000 samples a line is about 90 kB, so about 90 lines wait. The window budget
does not count them, as it does not count the walks' arrays.

**Tests** (`stage/tests.rs`). Lines read 1, 7, 333 or 1,024 at a time plan the same windows, count
the same tally and write the same bytes (`lines_read_one_at_a_time_or_in_batches_write_the_same_bytes`);
a window over budget stops the stage with the same `WindowOverBudget`, and a write refused at the
300th line stops it with the writer's error after the same 299 lines, however many lines wait.
The batch size is a field of the stage, set to the constant outside tests.

### R2 — the bgzf writer fed by a thread of its own

`noodles`' multithreaded bgzf writer queues at most as many blocks as it has workers, so the thread
writing to it blocks while they are all compressing — and with refinement on, that thread is the
one that hands the pool its next chunk of lines to rewrite.

`BgzfHandOff` (in `vcf/writer.rs`) sits in front of it: the writing thread appends to a 4-MB buffer
and hands each full one over a channel of one to a thread that writes it into the bgzf writer. At
most three buffers exist, about 12 MB. bgzf cuts blocks by the bytes alone, so the file is the same;
a flush is passed on in order, so it would cut a block where it did before (nothing flushes
mid-file today). The thread is one beyond `--threads` (the owner's rule: `--threads` is the least a
run uses). One thread and standard output keep their own sinks, as before.

**Errors.** The thread stops at the first failed write; the next hand-over, flush or `finish` finds
the channel closed and returns the thread's error with its own message. Dropped without `finish`,
the hand-off waits for its thread, so no thread outlives the writer. On tomato2 with the output's
size limited (`ulimit -f`, the hidden-duplication filter off so that the VCF is written while
calling), both main and the branch stop with "could not write to the in-flight output
`…calls.vcf.gz.tmp`: File too large (os error 27)", exit 1, leave the `.tmp` file and no output;
the locus named moved (6:30508230 on main, 6:42862380 on the branch) because the bytes now reach the
file up to 12 MB and 1,023 lines later. With the filter on, the spill fills first and the message
is identical on both.

**The prototype's switch `PVC_NO_HANDOFF` is not kept**: an environment variable is a hidden option,
and an A/B is a build of R1's commit.

**Tests** (`writer/tests.rs`). The hand-off writes the same file as the bgzf writer fed directly,
over 9.9 MB of text and a flush in the middle; a destination that refuses writes after 100 kB is
reported with its own message; a hand-off dropped unfinished returns.

### R3 — a line rewritten with its fields borrowed

`rewrite_line` copied every field of every sample into a `String` of its own, then joined each
column and the line. Now each field is borrowed from the line and owned only where the pass writes
it (`Cow<str>`), and the line is written into one buffer. This is the prototype's change applied
unaltered: `line.rs` had not changed since the prototype's base.

**Test** (`line/tests.rs`): a column that drops trailing fields gets `.` up to its `PS`.

## Output identity

Every build was compared with main (bf887b45, with the probe added, which changes no library code):

- **Through the real `VcfWriter`** (`refine_pipeline_probe`, 6 threads, the `.vcf.gz` hashed
  whole): six drawn cohorts, each with margins thinned to 5 records a kb and not thinned — 10
  samples × 300,000 records at 3 and 30 reads, 100 × 30,000 at 3 and 30, 500 × 20,000 at 3 and 30.
  Biallelic records with read fragments. All 12 hashes equal after R1, R2 and R3. The paired timing
  runs at 2 threads gave the same hashes again.
- **Multi-allelic records** (`refine_mem_probe`, the stage alone, 2 to 6 alleles): 10, 100 and 500
  samples at 3 and 30 fragments a base, thinned and not. All 12 hashes equal after each commit.
- **tomato2** (`pvc call-from-psps`, 59 samples, 32 Mb, `--refine-genotypes on`, 6 threads; the
  inputs the review used, under the main checkout's `tmp/parallel`): the decompressed VCF body, `##`
  lines dropped, hashes `1ce10c62d9c48469` on main and after each commit, and in all 16 timing and
  memory runs at 2 and 6 threads. The `##` lines differ only by the output path in the command line.

## Time

**Drawn cohorts** (`refine_pipeline_probe`, margins thinned to 5 a kb, every build run one after
another in each repetition; median of three, of two at 500 samples). *Writing thread* is that
thread's own CPU time: the work the pool cannot share.

| Cohort | Threads | Wall: main → R1 → R2 → R3 | Writing thread: main → R1 → R2 → R3 |
|---|---|---|---|
| 10 samples, 300,000 records | 2 | 11.91 → 11.65 → 11.60 → 11.50 s | 0.67 → 0.15 → 0.15 → 0.15 s |
| 10 samples, 300,000 records | 6 | 5.94 → 5.55 → 5.57 → 5.52 s | 0.69 → 0.17 → 0.17 → 0.16 s |
| 100 samples, 30,000 records | 2 | 24.31 → 24.20 → 24.05 → 24.08 s | 0.52 → 0.08 → 0.07 → 0.07 s |
| 100 samples, 30,000 records | 6 | 10.51 → 10.39 → 10.26 → 10.29 s | 0.53 → 0.09 → 0.08 → 0.07 s |
| 500 samples, 20,000 records | 6 | 33.08 → 31.48 → 31.66 → 31.44 s | 1.76 → 0.32 → 0.29 → 0.29 s |

The repetitions of one build spread by up to 0.3 s at 10 samples and up to 0.7 s at 100, so R1's
saving is clear at 10 and 500 samples (main's slowest and R1's fastest runs do not overlap) and
within the spread at 100. R2's and R3's steps are within the spread everywhere.

**tomato2** (paired runs, main then the branch, twice at each thread count). `pvc` prints the
calling time and the whole run's time to the second; the writing pass — the hidden-duplication
filter's last pass, which holds the refinement stage — is the difference.

| Threads | Writing pass, main | Writing pass, R1–R3 | Whole run, main | Whole run, R1–R3 |
|---|---|---|---|---|
| 2 | 97, 97 s | 96, 96 s | 212.97, 212.50 s | 211.74, 212.34 s |
| 6 | 43, 43 s | 41, 41 s | 93.92, 93.94 s | 90.63, 90.95 s |

In ten more runs at 6 threads made for memory (below), main took 93.66 and 92.45 s, R1 91.27 and
90.67 s, R2 91.39 and 91.03 s, R3 90.78 and 90.50 s: R1 is the step that saves time; R2 and R3 add
nothing beyond the spread.

**Why less than the review expected.** The review projected R1's saving from where it measured
reading at 2.7 s of tomato2's writing pass, at a load of 13 to 30 on a shared laptop; on a quiet
machine the whole serial step was smaller, and at 2 and 6 threads the solve is 90% or more of the
pass. R1's saving grows with the thread count — the serial reading does not divide by threads and
the rest does — so the review's 48-thread case is where it matters; it was not measured here. R2's
wait on the bgzf workers was 1.4 s on tomato2 at 2 threads in the review; no step of that size
appears here.

## Memory

**Resident memory** (`/usr/bin/time -l`, maximum resident set size):

| Run | main | R1 | R2 | R3 |
|---|---|---|---|---|
| tomato2, 6 threads (paired, five runs of main and R3, three of R1 and R2) | 483–508 MB | 484–503 MB | 513–562 MB | 508–555 MB |
| tomato2, 2 threads (two runs each) | 313, 341 MB | — | — | 335, 345 MB |
| drawn, 10 samples, 6 threads (three runs; the probe holds its input too) | 852–880 MB | 848–891 MB | 879–891 MB | 882–901 MB |
| drawn, 100 samples, 6 threads | 684–699 MB | 671–695 MB | 659–687 MB | 658–710 MB |
| drawn, 500 samples, 6 threads (two runs) | 1,735–1,764 MB | 1,831–1,859 MB | 1,776–1,788 MB | 1,768–1,793 MB |

**Bytes asked of the allocator** (dhat, the probe built with `dhat-heap`, 10 samples × 30,000
records, 6 threads, one run each): the peak while writing was 132.0 MB on main, 128.1 with R1, 139.9
with R2 and 129.3 with R3. R2's run asked for 11.8 MB more than R1's, inside the 12 MB its buffers
hold; R1's and R3's asked for less than main's, so a single run moves by a few MB with the threads'
timing. The program holds no more than its bounded buffers.

So the resident memory on tomato2 at 6 threads rises with R2 by more than R2 holds — about 30 to
70 MB, 6 to 14% of a 500-MB run — while what the program asks for does not. The likely cause is
the allocator (mimalloc) keeping freed memory per thread: with R2 the bgzf workers compress while
the pool rewrites the next chunk, so more threads allocate at once. Handing each written 4-MB
buffer back to be filled again, rather than allocating a new one, did not change it (539 and 540
MB in two runs); that change is not committed. At 2 threads the rise is within the spread.

R3 asks for 29% fewer blocks (10.9 to 7.8 million on the same run) and the same peak.

## Validation

- Touched modules' tests after each commit (`cargo test --release --lib -- neighbour_evidence
  vcf::writer run::paralog_filter`): 359 passed before R1's new tests; 362 after R1, 365 after R2,
  366 after R3; none failed.
- `cargo test --release --no-fail-fast` on the branch head: 5,385 passed, 2 failed, 9 ignored —
  the library's 5,323 passed and 2 failed, and 62 in the integration tests and examples. The two
  failures are the tests that expect a panic only debug builds raise
  (`alignment::left_align_structured::tests::an_offset_past_the_reference_panics`,
  `genetics::tests::lgamma_panics_on_non_positive_in_debug`); they fail under `--release` on main
  too.
- `cargo clippy --release --all-targets`: no warning in a touched file. One it raised in R1's test
  helper (a complex return type) is fixed by a type alias in the report's commit; the 11 left are
  main's, in `parameter_estimation`, `run/census_fit.rs`, an example and the vendored
  `noodles-cram`.
- Stops: over budget and a refused write, both on tomato2, both as described under R1 and R2.

## Decisions for the owner

1. **Keep R2 or revert it?** **Settled 2026-11-25: reverted** on main in 1ee9488c. It saves no time measurable at 2 or 6 threads on this laptop, on
   drawn cohorts of 10 to 500 samples or on tomato2, and it raises tomato2's peak resident memory
   at 6 threads by 30 to 70 MB (the bytes the program asks for stay within its 12 MB). The review
   argued it matters for wide lines and many threads; neither was tested here. **Recommendation:
   revert it** (`git revert c7082615` applies cleanly: R3 does not touch `writer.rs`) unless a
   measurement at many threads on a wide cohort shows the wait it removes.
2. **The record named when a run stops** is now the one being handed in when the stop is noticed,
   up to 1,023 lines after the one main named; the window described, and the advice, are
   unchanged. Keeping the old record would mean carrying the place of the line that revealed the
   stop through `StageStop` to the writer. Recommendation: leave it; the record named was never in
   the window.

## Follow-ups not done

- A measurement at many threads (the review's 48-thread projection) on a wide cohort, which is
  where R1 and R2 were argued to matter most. It needs a machine with more than six cores to
  spare.
- An observation, not a change: `rewrite_line` writes `PS` over a column's seventh field when the
  column has more fields than `FORMAT` (`PS` goes at `FORMAT`'s length). Calling never writes such
  columns; main and the branch behave alike.
