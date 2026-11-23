# Thinned and packed margins, and a core floor (2026-11-23)

Branch `margin-density`, from main d7c85b45.

## Why

The second tomato25 run (2,176 samples) stopped in the writing pass at chromosome 1, 42.25 Mb: the
smallest window there (0.5-Mb core and margins) was over the 90-GB budget even after the exact
memory cuts of 2026-11-22. A margin is read only to tell which samples a record's haplotypes
follow, and once a stretch has enough records to tell that, more say little: an experiment that
walked every second or fourth far-margin record (`margin_thinning_experiment.patch`, never
committed) left the fill check unchanged and cut memory and time. The owner asked to apply it
genome-wide, with a density above which margins are thinned, and to measure accuracy; then, when
the thinned window still did not fit 70 GB, to hold the skipped records compressed and to let a
core shrink below the margin floor, keeping only what does not hurt accuracy as the default.

## Changes made

- [`windows.rs`](../../../../src/neighbour_evidence/windows.rs): **`MarginThinning`**, a density
  in thousandths of a record a kb (`--refine-margin-density`, `DEFAULT_MARGIN_DENSITY` = 5);
  `NEAR_MARGIN_BASES` (50 kb), walked in full; `WindowBytes::skipped`, what a skipped record saves,
  paid only within the near margins; `WindowPlan::near_start`/`near_end`. **`core_floor_bases`**
  (`--refine-core-floor`, the margin floor by default): below the margin floor the core shrinks
  alone, its margins at the floor.
- [`stage.rs`](../../../../src/neighbour_evidence/stage.rs): **`MarginThinner`**, a token bucket of
  `D` tokens a kb holding at most a kb's worth, decides at each record's arrival whether the far
  margins walk it, so the choice depends on neither the window nor the thread count; records a
  sample's last written block joins at are walked too. **`HeldSites::pack`** compresses a skipped
  record's numbers and reads with zstd when its line is compressed; **`hold_whole`** unpacks the
  records a window walks and drops the unpacked copy of the others. The header states
  `margin_density`, `margin_near` and `core_floor` only where they apply, so an unthinned run's
  header is as before.
- [`records.rs`](../../../../src/neighbour_evidence/records.rs),
  [`fragments.rs`](../../../../src/read_phase/fragments.rs): `SiteRecord::take_numbers`/
  `put_numbers` and `RecordNumbers`, and `PackedRecordFragments::write_to`/`read_from`: a record's
  numbers and reads to bytes and back.
- [`refine.rs`](../../../../src/neighbour_evidence/refine.rs),
  [`templates.rs`](../../../../src/neighbour_evidence/templates.rs): the solve takes the window's
  `walked` mask; walks, and now **template choice**, read only the walked records — the only ones
  whose numbers are held. Template choice already compared at most 500 records a sample a window.
- CLI ([`call_from_psps.rs`](../../../../src/cli/call_from_psps.rs),
  [`call_from_alignments.rs`](../../../../src/cli/call_from_alignments.rs),
  [`resume_calling.rs`](../../../../src/cli/resume_calling.rs),
  [`calling_run.rs`](../../../../src/cli/calling_run.rs)): `--refine-margin-density`,
  `--refine-core-floor`, and the hidden `--refine-window-max` and `--refine-core-max` for
  measurements on short simulated contigs; the window options travel as `RefineWindowArgs`. The
  resume manifest keeps the density and the core floor; a manifest from before them reads as
  unthinned with the core floor at the margin floor.
- [`refine_replay.rs`](../../../../examples/refine_replay.rs): `REPLAY_MARGIN_DENSITY`,
  `REPLAY_CORE_FLOOR`; `REPLAY_FIXED=<core>,<margin>` now plans a core shorter than its margins.
- Spec: [neighbour_evidence_pass.md](../../ng/spec/neighbour_evidence_pass.md) §2.2, §6 (thinned
  margins). User docs: [neighbour_evidence.md](../../../user/neighbour_evidence.md) (new section on
  memory), [calling.md](../../../user/calling.md) (resume options).

## Measurements

**Replays** of refinement over the run's kept spill, on kimura, 48 threads, `refine_replay`.
*Dense slice*: chromosome 1, 44.0–45.5 Mb, 130,402 records, every window fixed at a 500-kb core and
500-kb margins, so arms compare genotype for genotype. The fill check hides about 288,000
genotypes per arm.

| Arm | Fill check wrong | Genotypes differing from ref | Peak | Wall |
|---|---|---|---|---|
| ref, every record walked | 0.0649% | — | 52.1 GB | 1,680 s |
| 20 a kb, not packed | 0.0635% | 0.037% | 40.8 GB | 1,219 s |
| **5 a kb, packed (the default)** | **0.0648%** | **0.039%** | **30.7 GB** | **821 s** |
| every record, 125-kb core | 0.0618%* | 0.039% | 33.2 GB | 2,691 s |

\* 1.1 million hidden genotypes: four times the windows hide four times as many. Changing only the
draws' seed moves 0.026% of genotypes (2026-11-22).

*Failing stretch*: chromosome 1, 40–47 Mb, 825,782 records, windows sized by the budget.

| Arm | Budget | Result | Highest window peak | Wall |
|---|---|---|---|---|
| 20 a kb, not packed | 90 GB | through 42.25 Mb (stopped by hand) | 108.4 GB | — |
| 20 a kb, not packed | 70 GB | stopped at 42.05 Mb: smallest window 73.5 GB | 84.1 GB | — |
| 5 a kb, not packed | 70 GB | stopped at 42.10 Mb: smallest window 72.3 GB | 83.3 GB | — |
| 5 a kb, not packed | 80 GB | ran, fill check 0.0800% | 97.8 GB | 7,587 s |
| 5 a kb, packed, core floor 100 kb | 60 GB | ran, cores down to 400 kb, fill check 0.0766% | 69.5 GB | 7,250 s |
| 5 a kb, packed | 65 GB | stopped at 42.20 Mb: smallest window 65.9 GB | 70.9 GB | — |

Real memory ran 1.16–1.22 times the budget. Packing saved about 6.4 GB on the smallest window at
42.2 Mb, less than the 15–20 GB estimated from the breakdown: the core's results and the near
margins' solve state, which packing does not touch, are most of that window.

**Simulated cohorts** on rick (refinement benchmark, 909ea135's psps, parameters and truth; all
samples, the scorer's `all` group), with `--refine-window-max 250000` so that the 1-Mb contigs have
several windows. The simulations hold fewer than 5 records a kb, so 5 a kb thins nothing.

| Arm | Outbred wrong | Outbred no call | Selfing wrong | F2 wrong | F2 no call |
|---|---|---|---|---|---|
| default windows (one a contig) | 20,179 | 20,412 | 8,694 | 3,765 | 20,677 |
| 250-kb windows | 19,867 | 18,992 | 9,126 | 6,686 | 33,314 |
| packed build, no thinning | 19,867 | 18,992 | — | — | — |
| 5 a kb, packed | 19,867 | 18,992 | 9,126 | 6,686 | 33,314 |
| 1 a kb, not packed | 19,506 | 18,463 | 9,117 | 6,770 | 33,620 |
| 1 a kb, packed | 19,553 | 18,581 | 9,335 | 6,689 | 35,287 |
| 0.1 a kb, not packed | 19,879 | 18,668 | 9,149 | 7,577 | 36,129 |
| 0.1 a kb, packed | 20,187 | 18,305 | 9,178 | 7,228 | 35,248 |
| 125-kb core, 250-kb margins | 19,610 | 18,469 | 9,406 | 7,531 | 37,496 |
| 50-kb core, 250-kb margins | 19,672 | 18,519 | 9,352 | 7,976 | 39,466 |

The packed build without thinning reproduces the unpacked output exactly. Short cores cost CPU
(+18% at 125 kb, +74% at 50 kb in the outbred cohort) and, in the F2, accuracy: 13% and 19% more
wrong genotypes. The F2 loses accuracy to anything that cuts its context — 250-kb windows nearly
double its errors against one window a contig — since its haplotypes run the length of a contig.

## Decisions

- **`--refine-margin-density` defaults to 5** (owner, 2026-11-23): on the densest tomato stretch it
  kept the fill check's error and cut memory 41% and time by half; the simulations, sparser, are
  unchanged; at 1 a kb and below they begin to lose genotypes.
- **`--refine-core-floor` defaults to the margin floor, 500 kb** (owner): a short core did not hurt
  tomato, but cost the F2 up to a fifth more wrong genotypes; it stays an option for a window that
  would otherwise stop a run.
- **The tomato25 rerun's budget is 80 GB** (owner, `snp_calling/pipeline.toml`), for a peak of about
  93 GB at the measured 1.16 ratio, under 100 GB.

## Tests

`cargo test --release`: 5,363 passed; 2 failed, the two that expect a debug-only panic and fail
under `--release` on main too. Clippy: no new warning. New tests: the thinner's spacing and its
independence from the thread count, thinned margins sparser than the density writing what
unthinned ones write, the planner paying skipped records only near the core, a core floor fitting
a window the margin floor cannot and changing nothing else, the density's parsing and the header,
and records' numbers and reads written to bytes and read back in every form.

## Open

- Real memory runs 16–22% over the budget's count; what the count misses is not yet known.
- Not reviewed. Not yet run genome-wide.
