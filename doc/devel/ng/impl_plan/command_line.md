# The `pvc` command line — implementation plan

**Status:** draft, 2026-10-16. The build order for [`command_line.md`](../spec/command_line.md),
whose design the owner settled on 2026-10-15 and 2026-10-16. **This plan is build order, not new
design**: where a step needs a choice the spec did not make, the step stops and asks.

**Branch:** `command-line`, worktree `../pop_var_caller-command-line`, from `main` at 98230b01.

---

## The constraint the order comes from

Another branch, `str-penalised-fit`, is about to merge into `main`. It rewrites the files the
parameter fit lives in:

- `src/cli/estimate_parameters.rs` and its tests;
- the parameters-file module, `src/calling/parameters_file/`, which moves to format 3;
- `src/run/census_fit.rs`, where the fit prints its progress;
- `src/cli/cross_platform_digests.rs` and `src/cli/test_fixtures.rs`.

So the steps that stay out of those files come first and can start now, and the steps that need
them wait for the merge. Before step 3 this branch is rebased onto the merged `main`.

## How each step is checked

The light development loop: each step is built, the tests it touches are run, and the step is
committed. The owner looks at the result before the next step starts.

**Steps 1 to 5 must not change a call.** The VCF's records stay byte-identical; what changes is two
header lines (§7.6 of the spec) and the filter's name, `hiddenParalog` → `hiddenDup`. That is
checked once, with one oracle run at the end of step 5, rather than after every step. Step 6
changes calls and has its own check.

---

## Step 1 — the logging foundation, piloted on `call-from-alignments` (now)

Spec §7.1, §7.2, §7.3 (the header part) and §7.7.

- Add `tracing` and `tracing-subscriber` to `Cargo.toml`.
- One module in `src/cli/` that installs two destinations: the terminal (stderr, `info` and above,
  `warn` and above under `--quiet`) and the log file (`debug` and above, written line by line). It
  knows whether stderr is a terminal, and draws the progress line in place only when it is.
- `--quiet` on every subcommand (a global option); `--log` on `call-from-alignments`; the default
  log `pvc_calls.log` in the VCF's directory, refused if it exists unless `--force` — which this
  step adds to `call-from-alignments`, for the VCF as well. `call-from-psps` gets both in step 2.
- `call-from-alignments` says what it is doing (the first line of spec §7.2), shows progress
  through the new module, and ends with the `done:` line. Its log gets the header block: version,
  command line, working directory, start time, host, threads; the inputs by name; every setting's
  value.
- The run report moves from stdout to the log.

**Done (2026-10-16).** Four tomato accessions over three regions of `benchmarks/tomato1`: 3,119
records, VCF identical to `main`'s apart from `##commandline`, parameters file identical. The
progress percentage line and the report's wording are unchanged; rewording them is step 5's.

**Done when:** a direct-mode run on the small test fixture prints the §7.2 shape on the terminal and
writes a log with the header, and its VCF is byte-identical to `main`'s.

## Step 2 — reshaping the command line (calling and psp parts now; `estimate-parameters` after the merge)

Spec §4, §5.3, §5.4 and §6.

**Done before the merge (2026-10-16):**

- `pvc dev` with `inspect-psp`, `type-regions`, `estimate-contamination`; subcommands listed in
  workflow order with one-line summaries; `pvc --help` describes the two workflows and what a psp
  is.
- Renames on both calling commands, flag and Rust field alike: `--max-locus-span`,
  `--locus-builder-regions-len`, `--max-generic-alleles` (its help says repeat loci are not
  affected), `--hidden-dup-fdr`, `--hidden-dup-tag`; the VCF filter `hiddenDup`. The fit branch
  builds none of these structs, so renaming the fields costs no conflict beyond
  `cross_platform_digests.rs`, which builds `CallFromPspsArgs` and had to change anyway.
- `--max-depth` replaces `--max-read-group-depth` on both calling commands and sets the read cap
  too; `--max-reads-per-position` is a *Developer* option that defaults to it. On `generate-psps`,
  which calls nothing and so has only the cap, `--max-depth` is the cap. **Its Rust field keeps the
  name `max_reads_per_position` until step 5**: `census_fit.rs` and the `estimate-parameters`
  tests build that struct.
- The repeat settings removed from `call-from-psps`, which reads them from the psp headers. The
  test that a run routing differently from its psps is refused became a test that psps written
  under non-default settings are called under them.
- `call-from-psps` gets step 1's treatment: `--force`, `--log`, its first line and `done:` line, the
  report in the log.
- Help rewritten for `call-from-alignments`, `call-from-psps`, `generate-psps` and
  `repeat-catalog`: one or two sentences an option, no spec paths or internal step numbers;
  *Advanced*, *What counts as a repeat* and *Developer* headings; `repeat-catalog`'s scan settings
  under *Advanced*, `--threads`' default of 1 explained.
- The three VCF checksums in `cross_platform_digests.rs` re-recorded: only the `##FILTER` line
  moved, and putting `hiddenParalog` back reproduces the old checksums.

**Checked:** four tomato1 accessions over three regions, both modes, 3,119 records each: the VCFs
are identical to `main`'s apart from `##commandline` and the filter's name
(`tmp/step2/run_both_modes.sh`).

**After the merge:**

- `regenerate-census` deleted, with its CLI module and tests, and `estimate-parameters`' refusal
  rewritten to tell the user to write old psps again with `generate-psps`. **Moved here from the
  first half** because that refusal builds a `regenerate-census` command line
  (`the_command_that_regenerates`, `THE_COMMAND_THAT_REBUILDS_A_CENSUS`), and both live in files
  the fit branch rewrites. The library's `census_from_psp` stays: three tests and an example use it
  to check the walk's census.
- `--threads` on `estimate-parameters`; `--str-strata-at-once`; its description's inbreeding
  sentence corrected; its help rewritten like the others'.

## Step 3 — what a run writes beside its output (after the merge)

Spec §7 opening, §7.3 and §7.6.

- Calling runs stop writing `calls.parameters.toml`. The refusal of a `--parameters` file that is
  the one the run would write goes with it. `parameters_file.md` §7 is marked superseded.
- The log gets the parameters the run scored with: their source, and the full values as a TOML block
  in the parameters-file format (format 3), copyable as a `--parameters` file — a test checks that
  the block read back gives the same parameters.
- The log gets the report and the end block: finish time, elapsed, peak memory, exit status, and on
  failure the full error.
- The VCF header gains `##pvcLog` and `##pvcParameters`; `##source` becomes `pvc <version>
  (<commit>)`, which needs the commit at build time.
- Logs for the other user subcommands: `pvc_call_parameters.log`, `<sample>.psp.log`,
  `<catalogue>.log`.

## Step 4 — `--output -` and `--tmp-dir` (now, after step 1)

Spec §7.4 and §7.5.

- `--output -` on both calling commands writes uncompressed VCF to stdout; the log then defaults to
  `pvc_calls.log` in the working directory.
- `--tmp-dir`, defaulting to the output's directory; the run makes `pvc_tmp_<random>/` inside it,
  puts the hidden-dup spill file there, and removes it however the run ends. Outputs in the making
  (`<name>.tmp`) stay beside their output.

**Done (2026-10-16).** The VCF writer gained a third destination beside plain text and bgzip,
standard output, which it only flushes at the end; the scratch directory is a `tempfile::TempDir`
made only when the filter is on, so a run with `--hidden-dup-fdr 0` makes none. Until step 3
removes it, a VCF on standard output puts its parameters file in the working directory as
`pvc_calls.parameters.toml`. Checked on the tomato1 psps of step 2: the same run to a file and to
stdout give the same VCF (3,119 records); both scratch directories were made inside `--tmp-dir`
and gone afterwards; and both modes still match `main`'s VCFs.

## Step 5 — every message through the logging framework (after the merge)

Spec §7.7.

- The remaining `println!` and `eprintln!` in the subcommands and the library become `tracing`
  events at the levels of spec §7.7, including the fit's progress lines in `census_fit.rs`.
- A test that fails on a `println!` or `eprintln!` in `src/` outside `src/main.rs` and `src/cli/`.
- The Rust field names follow the step 2 renames.
- `--quiet`, `--log` and the §7.2 first line on every user subcommand.
- **The oracle run**: the tomato and GIAB oracles, confirming the records are unchanged.

## Step 6 — `--call-satellites` and `--generic-path-only` (last)

Spec §4.4.1. New behaviour, not interface, so it gets its own short plan before code, written once
steps 1 to 5 are merged. It touches how ground is sorted, adds a field to the psp header, and needs
a parameters file without a stutter section. It also touches the STR-bundle code, which the bundle
performance work is changing on `main`; that work should have settled first.

**Its check is an accuracy one**: `--call-satellites` against GATK on GIAB over satellite ground,
before it is considered as a default.

## After the plan

The user documentation, written from the spec once step 5 is merged: the interface it describes
stops changing there.
