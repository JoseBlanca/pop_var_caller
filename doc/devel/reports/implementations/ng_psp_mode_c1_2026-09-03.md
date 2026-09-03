# ng psp mode — C1: the `generate-psps` subcommand

**Date:** 2026-09-03
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone C, step C1
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §2, §5.2
**Branch:** `ng-psp-mode`

## Plan

Give psp mode's walk stage a command line: reference, catalog, one `--alignment` per sample,
optional BED, `--output-dir`; assemble the reference and segmentation exactly as
`call-from-alignments` does; then one gatherer per sample, writing `<output-dir>/<sample>.psp`.

Two commits preceded it in Milestone C, both prerequisites rather than plan steps and both
recorded as such:

- `70385f5b` — the run-level shared test fixtures the B1 and B2 reviews both deferred. The
  catalog header claiming digest `[7;16]` had reached three copies and the not-default read
  filters two, because a `pub(super)` helper inside a `#[cfg(test)] mod tests` cannot be
  reached from a sibling file. Mechanical; no test changed what it asserts.
- `f00d56e9` — **the lift the plan's "reuse" requires.** C1's text says to reuse
  `segments_over`, `analysed_regions` and `build_read_groups`; the first two were private to
  `call_from_alignments.rs`, so reuse meant lifting them. They now live in
  `src/pop_var_caller_exp/run_ground.rs` behind a small `GroundRequest` each subcommand builds
  from its own flags, and the six refusals that belong to that assembly moved with them into
  `GroundError`, carried `#[error(transparent)]` by both commands so the rendered sentence is
  unchanged. Direct mode's 88 CLI tests pass untouched.

## Assumptions (choices the plan left open)

- **A sample's files are grouped by `@RG SM`, and the walk is per sample, not per flag.**
  The plan says "one `--alignment` per sample", but `build_read_groups` already establishes
  that files sharing one `SM` are one sample, and the gatherer refuses files naming two.
  So the command builds the run-wide table, then walks each `SampleReadGroups` entry with
  that sample's distinct files. `files_of` deduplicates because the table holds one row per
  read group, so a three-lane file appears three times in it.
- **The psp is named `<sample>.psp`, derived, with no per-sample output flag.** A cohort of
  psps is opened by naming files; a per-sample output flag would let two samples be written
  to one path with nothing noticing.
- **The output directory is created, not demanded** (`create_dir_all`). A walk is often the
  first thing run on a fresh machine.
- **Read filters and locus-generator settings are not flags**, matching direct mode, which
  hard-codes both defaults. They are recorded in every psp header regardless.
- **No thread flag**, per the owner's 2026-09-03 ruling: samples are walked one at a time and
  a cohort is parallelised by running invocations.
- **The census is not written.** Spec §2 gives the walk stage two files; Milestone G wires
  the second. Stated in the module doc so a reader is not left to infer it.

## Changes made

- `src/pop_var_caller_exp/generate_psps.rs` (new): `GeneratePspsArgs` (10 flags),
  `GeneratePspsCliError` (6 variants, one of them the transparent `Ground`),
  `run_generate_psps`, plus `ground_request`, `files_of`, `provenance`, `psp_path_for`.
- `src/pop_var_caller_exp/generate_psps/tests.rs` (new): 12 tests.
- Wiring: `mod.rs` (module + re-exports), `cli.rs` (the `GeneratePsps` subcommand with its
  help text), `src/main_exp.rs` (the dispatch arm).

## Tests added

Twelve, in three groups. **The command surface:** the subcommand spelling and its flags; the
repeating `--alignment` keeping its order; a walk with no alignment refused by clap; the
catalog and ground defaults; and — the one that guards the two modes against drift — the
routing defaults parsed through clap and compared against `call-from-alignments`' own, since a
psp written under different criteria is refused by a later calling run.
**The refusals:** a missing catalog rendering the shared message that names the file and the
command that builds it. **The files it writes** (driving `run_generate_psps` over a real
on-disk reference, catalog and two BAMs): one psp per sample named for the sample and each
opening; both psps recording the same analysed ground and the criteria actually routed with;
the provenance naming tool, subcommand and the file opened, with the applied read filters;
the output directory created when absent; two files naming one sample becoming one psp with
both read groups in its walk-local table; and a file holding several of one sample's read
groups being opened once.

## Validation results

In the container, on the tree reported:

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'pop_var_caller_exp'` — **100 passed, 0 failed, 1 ignored** (88
  pre-existing direct-mode tests plus the 12 new).
- **Run for real**, release build: one tomato accession
  (`tmp/tomato_slice/SRR7279481.p1.bench.cram`) over the first two 100 kb intervals of
  `benchmarks/tomato1/regions.bed`, through the real catalog — wrote
  `SRS3394712.psp`, **914,715 bytes, in 3.0 s**, exit 0. A run with a missing BED printed the
  shared refusal naming the file, exit 1.
- **Why that file is not the 948,689 bytes the B2 oracle wrote from the same slice —
  measured, not inferred.** The oracle harness asks the catalog with
  `StrRepeatCriteria::default()`, the file's own storage floors, where the command asks with
  ng's calling floors, so the two route different amounts of ground to the tract path. Running
  the command twice over the same slice, changing only the routing flags:
  **914,716 bytes at the calling floors against 949,322 with the floors relaxed toward the
  catalog's** (`--min-copies 5,5,5,5,5,5 --min-purity 0.8 --max-str-len 100`) — which lands on
  the oracle's figure, so the criteria account for the difference. The 1-byte gap from the
  914,715 quoted above is the recorded command line: that run wrote to a directory whose name
  is one character shorter.

## Tradeoffs and follow-ups

- **C2 and C3 are deliberately absent**: no per-sample report yet (a finished sample printing
  its `WriteStats`), and no `--force` or interrupted-rerun refusal. A rerun today overwrites
  a finished psp silently, which is `PspWriter::create`'s documented truncation — C3 is where
  that becomes a refusal.
- The five repeat-routing `#[arg]` blocks are near-copies of the sibling command's. Kept
  separate deliberately — the crate's convention is one module per subcommand owning its own
  `Args` — with the *values* pinned equal by a test rather than by sharing a struct.
