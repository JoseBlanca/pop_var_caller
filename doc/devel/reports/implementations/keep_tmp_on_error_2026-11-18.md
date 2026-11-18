# `--keep-tmp-on-error`: a failed run's scratch directory kept (2026-11-18)

Branch `keep-tmp-on-error`, from main 08959627.

## Why

On 2026-11-17 the tomato25 run (2,176 samples, `call-from-psps`, hidden-duplication filter and
genotype refinement on) stopped after 75 h 40 min in the filter's writing pass, when refinement's
smallest window at the start of chromosome 1 needed 88.3 GB against a 32 GB budget. The run's
scratch directory, holding the 5.4 TB spill of every called record, was removed on the way out,
so nothing of the 65 hours of calling could be looked at.

## Plan

A flag on both calling commands that keeps the run's scratch directory, with everything in it,
when the run stops on an error; a run that finishes removes it as before; without the flag nothing
changes.

## Assumptions

- **Only a returned error counts as a failure the flag acts on.** Release builds abort on a panic
  (`panic = "abort"`), so no `Drop` runs and the directory is left behind whatever the flag says —
  as `--tmp-dir`'s help already said. A debug build unwinds, and the unwinding counts as a failure.
- **A run succeeds at its last step**, after the report is logged; any `?` before that is a failure.
- **The spill inside the scratch directory no longer removes itself.** The directory owns its
  removal, so keeping the directory keeps the spill. A spill named beside an output
  (`SpillFile::beside`) is test-only and still removes itself.

## Changes made

- [`calling_run.rs`](../../../../src/cli/calling_run.rs): `scratch_directory_for` returns a
  `RunScratch` instead of a `tempfile::TempDir`. It removes the directory on drop unless
  `ScratchOnFailure::Kept` was asked for and `the_run_succeeded` was never called; a kept directory
  is logged at `WARN` with its path. A removal that fails is logged rather than silent.
- [`call_from_psps.rs`](../../../../src/cli/call_from_psps.rs),
  [`call_from_alignments.rs`](../../../../src/cli/call_from_alignments.rs): `--keep-tmp-on-error`;
  each run marks its scratch as succeeded as its last step. `--tmp-dir`'s help names the flag.
- [`spill_file.rs`](../../../../src/run/paralog_filter/spill_file.rs): a spill named with
  `in_directory` leaves its removal to the directory (`SpillRemoval::WithItsDirectory`).

## Tests added

- `calling_run::tests`: a failed run keeps its directory and the file in it when asked; a run that
  succeeded removes it even when asked to keep it; a failed run removes it by default; the flag's
  mapping.
- `spill_file::tests::a_spill_in_the_runs_scratch_directory_leaves_its_removal_to_the_directory`.
- `mode_equivalence::tests`, over the varying cohort walked into psps, filter on, refinement forced
  on with a one-byte window budget — the tomato failure in miniature: with the flag the run fails
  with "genotype refinement cannot go on" and leaves `pvc_tmp_*/hidden-dup-spill.tmp`, non-empty;
  without it `--tmp-dir` is empty; at the default budget the run finishes, writes at least three
  records, and leaves `--tmp-dir` empty though the flag was given.

## Validation

- `cargo test --release --lib`: 5,272 passed, 2 failed, 9 ignored. The two failures are the
  known tests that expect a debug-only panic.
- `cargo test --release --test '*'`: 33 passed.
- `cargo clippy --all-targets --all-features -- -D warnings`: the same 8 errors as main, all in
  files this branch does not touch (`parameter_estimation/joint/ssr_fit/*`, `run/census_fit.rs`).
- `cargo fmt --check`: clean on the files touched; main's two unformatted files
  (`calling/quality/mod.rs`, `calling/quality_parity.rs`) left as they are.
- `cargo test --all-targets --all-features` does not build on main either:
  `tests/ng_calling_loop_allocation.rs` and the untracked `examples/big_record_probe.rs`.

## Tradeoffs and follow-ups

- **A kept spill can only be looked at, not resumed from.** Re-running the filter's passes or
  refinement over a kept spill would need a command that reads one; none exists.
- A kept directory on a large cohort holds terabytes; the log line says to remove it.
