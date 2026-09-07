# Code Review: ng_paralog_filter_c4

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), three sub-agents in isolated worktrees
**Scope:** step C4 of the hidden-duplication filter plan — pass three, and the run that carries the filter end to end
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `627104f4` on branch `ng-paralog-filter`. Merge base `420b041f`.
- **In-scope files:** [pass_three.rs](../../../../src/ng/run/paralog_filter/pass_three.rs) and its
  tests, [finish.rs](../../../../src/ng/run/paralog_filter/finish.rs) and its tests,
  [vcf/header.rs](../../../../src/ng/vcf/header.rs), the two subcommands and their tests,
  [calling_run.rs](../../../../src/pop_var_caller_exp/calling_run.rs), the step's implementation
  report and the `PROJECT_STATUS.md` entry.
- **Out of scope:** `src/ng/paralog/*.rs` (guarded copies and the differential); production;
  `src/ng/window_coverage/`; steps A1–C3; **step C5's four end-to-end oracles**, which are the next
  step and not this one's.
- **Categories dispatched:** reliability + extras (mutation testing, named as the main task since
  this step's code had had none); errors + defaults + refactor_safety + smells; naming + idiomatic
  + module_structure + the diff's own numbers. `unsafe_concurrency` and `tooling` skipped.

### 2. Verdict

**Request-changes.** One Blocker, seventeen Major, twenty-eight Minor.

**The theme is that almost nothing this step wrote was pinned by a test.** The reliability reviewer
ran **31 mutations and 19 survived** — five in pass three, eight in the run's finish, four in the
header, two at the call sites. Every survivor was proved to change behaviour before being recorded.
They are not exotic: the report's four counts swapped, the header saying `3/1` where it means
`1/3`, `em_converged` negated, four decimals of the fitted rate lost, the calibration line moved
below the declarations, and both error messages naming a record as contig 0 position 0.

**And the Blocker is a fixture, not an oracle.** The only command-line test that runs the filter
runs it over a cohort that writes **zero records** — nothing is scored, dropped, tagged or written
through a subcommand, and no sample gets a coverage model. Direct mode has no such test at all.

### 3. Execution status

| command | at `627104f4` | at the merge base `420b041f` |
|---|---|---|
| `cargo test --all-features --lib` | `ok. 6621 passed; 0 failed; 15 ignored` | `6605 passed` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9 |
| `cargo fmt --check` | 9 unique files, none this step's | the same 9 |

One integration test fails at both — `main`'s. `--all-targets`, `cargo doc` and `cargo audit` not
run.

**Findings labelled "Needs verification": 0.** Every finding came from a mutation or a probe run in
a reviewer's own worktree.

**Mutation totals: 31 run, 19 survived, 3 changed no behaviour, 9 killed.**

### 4. Open questions and assumptions

1. **⛦ Spec §10's second oracle cannot pass as written, and C5 is the step that runs it.** It says:
   filter on at an unreachable target, strip the two `INFO` keys, and the result must equal the
   filter-off file exactly. It cannot — the on-run's header carries four lines the off-run does
   not. That follows from deviation 1, and the alternative (declaring them always) breaks the
   *first* oracle instead. **The two as written are incompatible and the first is the one that
   matters.** Recommendation: amend §10's second oracle to strip the filter's four header lines as
   well. **A spec edit, so it is the owner's, and C5 needs it settled.**
2. **Is a record-producing command-line fixture worth building now?** The current cohort — two
   samples, three reads — cannot fit a coverage model or write a record. Bears on the Blocker; see
   its entry for what was done instead.

### 5. Top 3 priorities

1. **B1** — nothing puts a record through a subcommand's filter path.
2. **M1** — a target of zero would have removed records had it reached the scoring; the off switch
   lived in two copied comparisons.
3. **M2–M13** — nineteen surviving mutations across the report, the header and the errors.

### 6. Findings

#### Blocker

**B1: src/pop_var_caller_exp/call_from_psps/tests.rs — the only command-line test that runs the filter runs it over a cohort that writes no records, and direct mode has none**

**Confidence:** High. **Categories:** reliability.

Measured by probe on the shipped fixture:

```
PROBE3 psps: dropping wrote 0 record(s), tagging wrote 0, of which 0 carry the filter
hidden-duplication filter: 0 record(s) dropped, 0 tagged, 0 written; 0 were scored, 0 could not be
  duplication rate 0.030000 (…), cut at 8.1500, over 0 of 2 sample(s) with a coverage model
```

Negating the tag flag at either call site left the output byte-identical — one of the three
mutations that "changed no behaviour", and the reason it changed none. Direct mode's path works
(probed) but nothing would notice if it stopped.

**Partly fixed, and the remainder is recorded.** The filter's whole path — fit, score, calibrate,
write, report — is now exercised with real records that are actually removed and tagged, through
the one function both subcommands call
(`the_report_counts_what_went_and_what_stayed`, on an asymmetric three-record fixture). What is
still untested is the ~45 lines of per-subcommand wiring around that call. **Building a
record-producing command-line cohort is deferred**: the shared fixture is two samples and three
reads, and giving it enough coverage to fit a model and write records changes the ground every
command-line test in three files runs on. It is filed for C5, which builds real-data runs anyway.

#### Major

**M1: pass_two.rs — a target of zero would have removed records, and the off switch was two copied comparisons**

**Confidence:** High. **Categories:** reliability, defaults, refactor_safety.

Measured: one record at ratio 120 with `target_fdr = 0.0` gives `dropped: 1`. The tail
false-discovery value underflows to exactly zero above about ratio 40, and zero is not above zero.
Pass two's own doc records a ratio of 1,663 at 63 samples, so this is the ordinary case at cohort
size, not an edge. `--paralog-fdr 0` was "off" only because two call sites each wrote
`if target_fdr.get() > 0.0`.

**Fixed by making it unspellable**: `TargetFdr::try_new` refuses zero, and
`WhatTheOperatorAskedFor::from_the_flags` is the one place where zero becomes *no filter*,
returning `Option<Self>`. A target in hand now always means the filter runs.
`a_target_of_zero_is_no_filter_rather_than_a_filter_that_removes_nothing` holds it.

**M2–M13: nineteen surviving mutations**

**Confidence:** High. **Categories:** reliability.

| what a mutation changed | now killed by |
|---|---|
| the report swaps `dropped` and `tagged`, or `written` and `unscored` | `the_report_counts_what_went_and_what_stayed` |
| tag mode drops anyway (the flag hardcoded at the call to pass three) | same |
| the header swaps the fitted count and the cohort size (`3/1` for `1/3`) | `the_header_says_how_much_of_the_cohort_was_fitted_and_whether_the_rate_was` |
| the header negates `em_converged` | same |
| the fitted rate or the cut loses four decimals, in the header or the report | `the_provenance_writes_the_cut_and_the_rate_to_the_records_own_precisions`, `the_report_writes_the_rate_and_the_cut_to_their_precisions` |
| `em_converged=` deleted from the calibration line | `the_provenance_writes_the_cut_and_the_rate_to_the_records_own_precisions` |
| the calibration line moves below the `INFO` declarations | `the_calibration_line_comes_before_the_info_declarations` |
| the saturation line is never printed | `the_report_says_when_ratios_sat_past_the_ends_of_the_range` |
| `PARALOG_POST` loses four decimals | `the_two_info_fields_are_written_to_the_headers_precisions` |
| the patch failure's ordinal counts from zero, or names contig 0 position 0 | `a_line_that_cannot_be_patched_names_its_record_and_its_ordinal` |
| the write failure names contig 0 position 0 | `a_write_failure_names_the_record` |
| a ratio vector **longer** than the spill is accepted and truncated | `a_ratio_vector_longer_than_the_spill_is_refused` |

The two remaining survivors are stdout-only — the run report's lines, and their order — and
nothing in the suite captures stdout. `calling_run.rs`'s own doc already records that class.

**M14: finish.rs — the histograms were cloned, so spec §5's "dropped once the models are fitted" did not hold**

**Confidence:** High. 80.2 kB a sample — 5 MB at tomato's 63, **241 MB held to the end of the run**
at spec §4's three thousand. **Fixed**: `std::mem::take` at both call sites; nothing else reads the
field.

**M15: finish.rs — the report named a rejected sample by its index, three lines under a comment headed "Named, not counted"**

**Confidence:** High. An index is a fact about the order the run's arguments were typed in.
**Fixed**: `FilteredRun` carries the sample names and the report prints them.

**M16: pass_two.rs / finish.rs — C3's fallback warning was never printed**

**Confidence:** High. `why_the_paralog_rate_is_not_fitted` says of itself that the run report prints
it, and C3's report listed printing it among C4's jobs; it was called only from tests. What C4
shipped instead dropped the actionable half — that the removed records were calibrated against a
constant rather than against this run. **Fixed**: printed, and the rate line's parenthesis
shortened so the two do not say the same thing twice.

**M17: vcf/header.rs — `header_text` read the metadata field by field**

**Confidence:** High. So the next field added is silently missing from the header — which is what
nearly happened to this step's own provenance field. **Fixed**: destructured exhaustively.

#### Minor

Applied: the filter's id and the two rounding precisions had two spellings each and now have one,
in `vcf/header.rs`, held to the declaration by
`the_declared_filter_id_is_the_one_a_record_carries` (a record carrying an id its header does not
declare is invalid VCF); `paralog_filter/mod.rs`'s inventory omitted `finish`, the one module both
subcommands call; a stale doc paragraph in `call_from_psps/tests.rs` still said non-zero targets
were refused; the flag's help said nothing about the disk the filter needs.

Deferred, with reasons in the fix-application report: the ~73 lines of wiring duplicated between
the two subcommands; three configurations chosen inside `fit_score_and_write_the_calls` that cannot
be overridden or reported, and the error variant that is unreachable because of it; the renames
(`WhatTheFilterDid`, `WhatTheOperatorAskedFor`, `what_to_tell_the_operator` say what a value is
*about* rather than what it *is*); `write_the_records_the_filter_kept` also writing the records the
filter did not keep; and `PARALOG_POST` printing `0.000000` for a confidently-clean record.

### 6a. The diff's own numbers

**One wrong, repeated in three places; everything else correct.**

| claim | verdict | correct value |
|---|---|---|
| the gate's library count, in the report's table, its prose and the commit message | **WRONG** | **6,621**, not 6,619 — and the delta is **16**, not 14. Confirmed by re-running at `627104f4` and at the merge base |
| "6,619 is 6,605 plus the 14 library tests this step adds; the other two are the command-line ones, which are in the same lib target" | **WRONG, and self-contradicting** | the sentence concedes the pair is in the same target and then omits it |
| "clippy first went to 13, in four kinds" | **UNVERIFIABLE** | describes a working tree no commit preserves. Kept, now marked as such |
| 16 tests; the 9 / 5 / 2 split; 6,605 at the merge base; 9 clippy errors in three kinds; 9 fmt files; the four-and-six decimals; the three tail values `0.667`, `1.17e-5`, `0.500` and which records they remove at which target; the filter-off header carrying none of the four lines; the oracle test comparing against the parked entry's own bytes; a fitted rate of certainty being reachable | **CHECKED-CORRECT** | — |

**The mechanism, because it is the transferable part:** the gate was run, then two more tests were
added, then the report was written from the earlier run. That is the first entry in this project's
own writing log — prose written from the memory of a measurement rather than after the run it
reports — and the tell was that the sentence had to explain away a discrepancy of exactly two.

### 7. Out of scope observations

- **Nothing removes `<output>.tmp` when a run fails between the writer's creation and its finish.**
  Pre-existing, but this step lengthens that window by a whole pass over every record; on the
  filter path the spill's `Drop` removes the parked records while the partial `.tmp` stays.
- **`did.unscored` is reached only because the verdict screens on finiteness first.** If that ever
  moved, an unscored record would be dropped and never counted. Same shape as C3's Blocker, on the
  other side of the boundary.
- Two mutations at the call sites are invisible because nothing captures stdout — a class
  `calling_run.rs` already documents.

### 8. Missing tests to add now

All applied except the record-producing command-line cohort (B1's remainder), which is filed for
C5.

### 9. What's good

- **The spill's `Drop` removes the parked records on every unwind path traced, including a panic**,
  and `format_error_chain` renders the whole chain, so no failure loses its context.
- **`HiddenParalogProvenance::as_header_value` destructures its own fields exhaustively**, so a
  field added to the provenance cannot be silently left out of the line — which is exactly the
  defect found four lines away in `header_text`.
- **The oracle test compares a written line against the parked entry's own bytes**, not against a
  re-encoding.
- **The three mechanism claims in the step's report all hold**, checked by measurement rather than
  by reading: the conditional header lines, the byte comparison, and the reachability of a fitted
  rate of certainty.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo test --all-features --lib --bins --tests
./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings 2>&1 \
  | grep -E "^error: " | grep -v "could not compile" | sort | uniq -c
./scripts/dev.sh cargo fmt --check          # counted by unique file
./scripts/dev.sh cargo test --lib --all-features ng::run::paralog_filter
./scripts/dev.sh cargo test --lib --all-features ng::vcf::header
```
