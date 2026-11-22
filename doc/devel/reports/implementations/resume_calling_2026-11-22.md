# `pvc resume-calling`: a run finished from the scratch directory it kept (2026-11-22)

Branch `resume-from-scratch`, from main c60a9bf6.

## Why

On 2026-11-22 the second tomato25 run (2,176 samples) stopped after 86 h: calling (65 h) and the
hidden-duplication filter's scoring (7 h 20 min) had finished, and the writing pass stopped where
genotype refinement's smallest window at chromosome 1, 42.25 Mb, needed 93.3 GB against a 90 GB
budget. `--keep-tmp-on-error` kept the 5.4 TB spill, but nothing could read it: the filter's
passes also need each sample's coverage histogram, which lived only in memory and cannot be rebuilt
without reading every sample's reads again. Every new attempt was another 65 h of calling.

## Plan

Keep beside the spill, as soon as calling ends, everything passes two and three take besides it;
add a command that reads it back and runs the two passes through the function a run calls.

## Assumptions

- **The resume point is kept by every run with the filter on**, without a flag: about 175 MB at
  2,176 samples (the histograms), against a spill of terabytes. A run that cannot keep it logs a
  warning and goes on.
- **A resumed run may change only what the writing does**: `--refine-window-bytes`,
  `--refine-window-floor`, `--phase-min-path-agreement`, `--threads`, the output and the log. What
  was called, and what the filter was asked, are the stopped run's.
- **The VCF's header names the stopped run's command line and log**, the run that called.
- **A resumed run that finishes removes the directory**, as a run would; with `--keep-tmp-on-error`
  it keeps it again if it fails.
- **`--explain-loci` is not carried over**: a resumed run explains nothing.
- **The spill's layout is numbered** (`spill::SPILL_FORMAT`, 1), since the spill itself has no
  version; a build whose number differs refuses the spill.
- **What cannot be reported**: calling's own report, which the stopped run logged.
- **The 2026-11-22 spill cannot be resumed**: it was written before the resume point existed.

## Changes made

- [`resume.rs`](../../../../src/run/paralog_filter/resume.rs) (new): `save` writes
  `coverage-histograms.tsv` (one line a sample, `recorded_windows::write_the_histogram`'s form:
  every cell, the bin width by its bit pattern) and then `resume.json` (spill entry count and layout
  number, ploidy, inbreeding by bit pattern, the filter's settings, the header's metadata, the
  refinement decision), each under a temporary name renamed into place; `load` reads them back
  through the types' own checks; `is_resumable`.
- [`finish.rs`](../../../../src/run/paralog_filter/finish.rs): `fit_score_and_write_the_calls_explaining`
  saves the resume point first — the one place both calling commands pass through after calling.
- [`recorded_windows.rs`](../../../../src/run/cohort_merge/recorded_windows.rs): `read_the_histogram`,
  the reader of `write_the_histogram`'s lines.
- [`spill_file.rs`](../../../../src/run/paralog_filter/spill_file.rs): `SpillFile::finished_in` (a
  kept spill, ready to read), `directory` (the run's scratch directory, `None` for a spill named
  beside an output), `path` no longer test-only.
- [`spill.rs`](../../../../src/run/paralog_filter/spill.rs): `SPILL_FORMAT`.
- [`calling_run.rs`](../../../../src/cli/calling_run.rs): `RunScratch` holds a path and removes the
  directory itself, so `kept_scratch_directory` can take over a kept one; the kept-directory warning
  names `pvc resume-calling` when the directory holds a resume point;
  `CallingRunError::ScratchDirectoryNotFound`.
- [`resume_calling.rs`](../../../../src/cli/resume_calling.rs) (new): the subcommand.
- [`header.rs`](../../../../src/vcf/header.rs): `low_base_quality_test()`.
- `doc/user/calling.md`, `doc/devel/ng/spec/command_line.md` §7.5.

## Tests added

- `recorded_windows::tests`: a histogram of each kind survives being written and read back, bit for
  bit; a line whose cells do not match its scheme is refused.
- `resume::tests`: a resume point reads back as written (histograms, inbreeding bits, the filter's
  settings, the header with its refinement decision, ploidy, entry count); one from another spill
  layout is refused; a directory without a manifest is not resumable.
- `resume_calling::tests`: a directory with no resume point, and one that is not there, are refused.
- `mode_equivalence::tests::a_run_resumed_from_its_kept_scratch_writes_what_an_unbroken_run_writes`:
  the 2026-11-22 failure in miniature — refinement forced on with a one-byte budget, the run fails
  and keeps its directory; resumed at the default budget, it writes exactly the lines of the same run
  made at that budget from the start, but for `##commandline` and `##pvcLog`; the directory is gone
  afterwards.
- The kept-directory test now expects the resume point beside the spill.

## Validation

- `cargo test --release --lib`: 5,287 passed, 2 failed (the two tests that expect a debug-only
  panic, as on main), 9 ignored. `cargo test --release --test '*'`: 33 passed.

## Tradeoffs and follow-ups

- **The scoring pass is run again** on resuming (7 h 20 min at 2,176 samples). Keeping its verdicts
  too would skip it when only the writing changes; it would also need the scoring model's report.
- The resume point is written by `call-from-alignments` too, untested there beyond the shared
  function.
