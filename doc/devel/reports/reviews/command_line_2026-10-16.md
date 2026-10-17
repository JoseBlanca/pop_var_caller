# Review: the command-line branch (2026-10-16)

**Scope.** Branch `command-line` against `main` (`git diff main...command-line`, eleven commits,
`01e017e2` to `ff1d2538`): the revision of `pvc`'s command line and of what a run prints and
writes, built to `doc/devel/ng/spec/command_line.md` steps 1 to 5. Three reviewers read the diff,
none of them the author: one general correctness pass, one on the log's global slot, the stdout
writer and the parameters-block extractor, and one on help text against behaviour and on the
scripts that now read logs. Every finding below was checked against the code before it was
accepted.

**Verdict.** No defect moves a call: every finding is in what a run writes beside its calls (the
log, the terminal, the help) or in the scripts around the binary. Two kinds mattered:

- **Scripts that no longer ran.** Thirteen scripts either read the `<output>.parameters.toml`
  that calling runs no longer write, or made two calling runs into one directory, where the
  second is refused because the first left `pvc_calls.log`. One more — the paralog-filter
  comparison — ran but failed its first check on every run, because the two VCFs it compares
  named different logs in their headers.
- **Logs that misdescribed their run.** A failed run's log refused the corrected retry; every
  per-sample log of `generate-psps` but the last lacked its "finished" line; a panic left a log
  that ended like a killed run's; `repeat-catalog` recorded the wrong thread count.

All of these are fixed. Eleven findings are left as they are, each with its reason in §2.

## 1. Fixed

| # | Finding | Fix |
|---|---|---|
| 1 | **Scripts read the removed parameters file, or are refused on their second calling run.** `ng_psp_concurrency_invariance.sh` moved `run.parameters.toml`; `run_tract_qual_experiment.sh` and `ng_str_stutter_benchmark.sh` (with `ng_str_stutter_true_parameters.py`) started from the defaults run's parameters file; `run_ng_coverages.sh` no longer left the `HG002_<cov>.raw.parameters.toml` that `sweep_tract_parameters.sh` requires, so the sweep always stopped. Each of those, and `portable_float_baseline.sh`, `ng_str_accuracy/{call_and_score,hg002_two_arms,hg002_two_binaries}.sh` and `benchmarks/lib/run_ng.sh`, makes several calling runs into one directory without `--log`, so every run after the first was refused on `pvc_calls.log`. | Each calling run names its own `--log`; where a parameters file was read, it is copied out of the log with the same `sed` the adapted oracle uses. Scripts that re-run into an existing VCF by design (`call_and_score.sh`, `sweep_tract_parameters.sh`, `run_ng_coverages.sh`, `run_ng.sh`) pass `--force`, which keeps the overwrite they had on `main`, where the calling commands did not refuse an existing VCF. `hg002_two_binaries.sh` passes `--log` only to a build whose help lists it, since its base build may predate the option. |
| 2 | **`ng_paralog_filter_runs.sh` failed its first relation on every run.** Each side wrote `--log $side.pvc.log`, the header records the log's name (`##pvcLog=`), and the check that the filter changed no record it did not flag compares the filter-off and filter-on headers. The script's own comment still described the deleted `##parametersFile=`. Separately, its "what the filter told the operator" section dropped the filter's warnings: a report line beginning `warning:` is logged as a `WARN` event, without the `report:` prefix the script searches for. | Every side writes `run.pvc.log`, moved after the run like `run.vcf`; the `WARN` lines are read too. |
| 3 | **A failed run's log refused the corrected command.** Every subcommand opened its log with "refuse if it exists" after checking only some of its arguments. A run that failed afterwards — a reference mismatch, `estimate-parameters --ploidy 0`, a bad parameters file — left the log, and the fixed command was refused with "a log is already at …" although no output existed. | `open_log` replaces, without `--force`, a log whose last line records a non-zero exit status (`run_output.rs`, `left_by_a_failed_run`). A log that ends any other way — a run that succeeded, one still running, one killed — is still refused, so two runs in one directory still cannot overwrite each other's record. Checked by hand: a failed `repeat-catalog`, run twice, fails the same way twice instead of being refused the second time. |
| 4 | **Only the last per-sample log of `generate-psps` said how its run ended.** Opening the next sample's log replaced the first in the global slot without closing it, so `a.psp.log` and `b.psp.log` of a three-sample run ended with no exit status, elapsed time or peak memory, and read as killed runs (spec §7.3 item 7). | Each sample's log is closed, with exit status 0, when the next sample's opens; the last is closed by the binary after the cohort's report, which stays in the last sample's log as before. |
| 5 | **A panic left the log ending at its last event.** Release builds set `panic = "abort"`, so nothing after a panic runs: the log got neither the error nor the "finished" line, and could not be told from a kill. | The binary installs a panic hook that keeps Rust's terminal message and then writes the panic and a "finished: exit status 101" line to the log (`run_output::close_log_after_a_panic`). |
| 6 | **`--tmp-dir` was first checked when the scratch directory was made**, after the reference had been read and verified and every input opened — minutes on a genome. | Both calling commands refuse a `--tmp-dir` that is not a directory with the rest of the typed arguments (`calling_run::refuse_a_tmp_dir_that_is_not_a_directory`). |
| 7 | **`call-from-psps` judged `--max-locus-span` and `--locus-builder-regions-len` after its log opened**, against its own comment that a typing mistake leaves no log. | Judged before the log opens; the psp list it needs is already expanded there. |
| 8 | **`repeat-catalog` stated the wrong thread count.** The terminal printed `--threads 0` as "0 threads" although the scan uses 1, and the log header's `threads:` line read rayon's global pool, which nothing had sized, so it recorded every core. | The count is taken as `max(1, --threads)` once, printed, and used to size the pool the scan runs on, which the header reads. |
| 9 | **`--quiet` was missing from the log's settings** (spec §7.3 item 3: every setting). `settings_of` read the subcommand's options from a command clap had not built, and clap attaches global options only when it builds. | The command is built first. Checked: `setting: quiet = true` now appears. |
| 10 | **Help and documentation that described old or unbuilt behaviour.** `--tmp-dir` promised the scratch directory is removed "however it ends"; a kill or a panic leaves it. `--hidden-dup-fdr` placed the filter's scratch "beside the output" although it goes to `--tmp-dir`. `call_from_alignments.rs`'s module documentation still described `<output stem>.parameters.toml`; the README listed `regenerate-census`; `run_ng.sh` and `benchmarks/lib/README.md` said the run report is teed to the script's log. | Reworded to what the code does. |

## 2. Not fixed, and why

| # | Finding | Why not |
|---|---|---|
| 11 | **One default log name a directory**: two calls writing `chr1.vcf.gz` and `chr2.vcf.gz` into one directory collide on `pvc_calls.log`. | The owner's ruling (spec §7.3): the second run is refused, and runs that share a directory each name their own `--log`. The refusal is made by the open itself (`create_new`), so two runs started at once cannot both take the file. |
| 12 | **`--force` replaces the log when the run starts, but the VCF only when it succeeds.** A forced re-run that fails leaves the earlier VCF beside a log of the failed run; the VCF's `##pvcLog` then names a log that does not describe it. | The alternatives each cost something the project has chosen against: deleting the earlier output when the run starts gives up the file a failed re-run should not destroy (the reason `generate-psps` writes a `.partial` and renames it), and keeping the earlier log under a second name adds a file every forced run leaves behind. The failed run's log ends with its error and "exit status 1", so it does not pass for the earlier run's. **Recommendation: leave it.** |
| 13 | **`generate-psps` keeps no log of a failure before its first sample** (an unreadable reference or catalogue). | Its logs are per sample, and there is no sample yet to name one after. The error is on the terminal. |
| 14 | **`estimate-parameters` judges `--ploidy` and `--str-param-estimates-at-once` after its log opens**, so a refused value leaves a log. | After fix 3 that log no longer blocks the corrected command; it records the refusal, which is harmless. Moving the checks means lifting them out of `fit_and_assemble`, which its other callers rely on. |
| 15 | **The scripts' `sed` extraction of the parameters block, given a log whose block was cut short, prints to the end of the file and drops the last line**, where `run_output::parameters_in` returns nothing. | A block is written in one write, so only a log whose run never reached it can lack the end marker, and every script extracts only after a calling run that exited 0. |
| 16 | **`--output - \| head` ends the run with exit status 1** and a message worded for a file ("could not write to the in-flight output `standard output`"). | Spec §7.4 accepts exit status 1 on a stream that stops; only the wording is off. |
| 17 | **Library modules import from the CLI module** (`run::calling_progress` and `parameter_estimation::progress` use `cli::run_output`'s `PROGRESS` and formatters), and `run_output` re-implements `progress.rs`'s peak-memory reading and size formatter. | No behaviour follows: the two size formatters give the same text at every size (whole MiB below a GiB, GiB to one decimal above). A move into a module of their own is a refactor for when step 6 reworks these files. |
| 18 | **The log does not record how much scratch space the run used** (spec §7.5); it records only where the scratch directory was. | Not built in steps 1–5. Left for step 6. |
| 19 | **`--hidden-dup-tag` with `--hidden-dup-fdr 0` is accepted and does nothing.** | Unchanged from `main`; the options were renamed on this branch, not given new behaviour. |
| 20 | **`repeat-catalog` leaves `<catalogue>.tmp` behind when it fails** (seen while checking fix 3). | Unchanged from `main`: this branch does not touch the catalogue writer. |
| 21 | **`benchmarks/giab/src/ng_missed_sites_probe.sh`**: under `set -e`, a run that fails before opening its log makes the `sed` on that log stop the script, where the old `\|\| true` carried on; and its pattern `repeat tracts this caller has not built yet:` matches nothing since the repeat path was built, so `tract_bp` is always 0. | The first only changes how a broken run stops; the second predates this branch. |

## 3. Areas checked and found sound

- **The log's global slot.** No test installs the logging subscriber (only `main` calls `install`), so tests running in parallel cannot put events into each other's logs. The parameters block goes through the per-run `RunLog` handle, never through the slot, so the tests that read it back are deterministic. In a real run, one process has one log open at a time (generate-psps after fix 4), and every worker thread's events reach it. The cost is a gap in coverage: no test checks that the header, events, report and "finished" line reach the right log.
- **The stdout writer.** `--output -` is always uncompressed; there is no `.tmp`, rename or index; nothing else in the calling commands writes to stdout; the default log is `./pvc_calls.log` and the header says `##pvcLog=pvc_calls.log`, both as spec §7.3 and §7.6 say.
- **The parameters block.** Written in one write directly to the file, so no timestamp or level can land inside it; markers matched as whole lines, tolerant of CRLF; at most one block a log; a cut-short block yields nothing. The text read back is accepted by `--parameters` in the calling commands' tests.
- **Help against behaviour** for `--max-depth`, `--max-reads-per-position`, `--output -`, `--threads`, `--force`, `--log` and the renamed options: defaults and scope match the code.

## 4. Verification

- `cargo test --release --no-fail-fast -- --skip an_offset_past_the_reference_panics --skip lgamma_panics_on_non_positive_in_debug`: 5,210 passed, 0 failed, 7 ignored.
- `cargo clippy --release --lib --tests --bins`: no warning in the crate's own code.
- Every changed script passes `bash -n`.
- **The identity oracle was not re-run**: no fix reaches calling. They change what is written to the log and the terminal, which arguments are checked before the log opens, the help text, and the scripts; `repeat-catalog`'s thread pool changes only how many contigs are scanned at once, and the catalogue is the same at any number.
