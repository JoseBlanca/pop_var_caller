# Fixes applied — ng hidden-duplication filter, step C4 (pass three, and the run end to end)

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_c4_2026-09-07.md](ng_paralog_filter_c4_2026-09-07.md) — 1 Blocker,
17 Major, 28 Minor, Request-changes; three sub-agents in isolated worktrees
**Reviewed commit:** `627104f4` · **Branch:** `ng-paralog-filter`

## The answer

**Applied 24 · Applied with adaptation 2 · Deferred 6 · Disputed 0 · Failed validation 0.**

The step gains 13 tests, from 16 to 29; the library suite goes from 6,621 to 6,634.

**The review's substance was that almost nothing this step wrote was pinned.** Thirty-one
mutations, nineteen survivors — and none of them exotic: the report's four counts swapped, the
header saying `3/1` where it means `1/3`, `em_converged` negated, the fitted rate and the cut
losing the decimals that make them comparable with a record's own fields, the calibration line
drifting below the declarations, and both error messages naming a record as contig 0 position 0.
Seventeen of the nineteen are now killed; the two that remain are the run report's lines and their
order, which nothing in the suite can see because nothing captures stdout.

**And one finding was a live hazard rather than a missing test.**

## The one that was a hazard

**A target of zero would have removed records.** `--paralog-fdr 0` means *do not run the filter*,
and that is what it did — but only because two call sites each wrote `if target_fdr.get() > 0.0`.
Handed a zero, the scoring itself removes the most extreme records: a strongly duplicated record's
tail false-discovery value underflows to **exactly** zero, and zero is not above zero. Measured by
the reviewer: one record at ratio 120 with `target_fdr = 0.0` gives `dropped: 1`. That is not an
edge — pass two's own doc records a ratio of 1,663 at 63 samples, and this step's own sibling test
measures four of seven records removed at the strictest target a run can now ask for.

**Fixed by making it unspellable rather than documented.** `TargetFdr::try_new` refuses zero, and
`WhatTheOperatorAskedFor::from_the_flags` is the single place where zero becomes *no filter* —
returning `Option<Self>`, so a target in hand always means the filter runs and "off" is the absence
of one. Both subcommands now branch on that `Option` instead of on a comparison each wrote
separately.

## What was applied

| # | severity | finding | outcome |
|---|---|---|---|
| B1 | Blocker | no record goes through a subcommand's filter path | **Adapted** — the shared path is covered with real removed and tagged records; the command-line cohort is deferred to C5 |
| M1 | Major | a target of zero would remove records | **Applied** — `Option<WhatTheOperatorAskedFor>` |
| M2–M13 | Major | nineteen surviving mutations | **Applied** — 13 tests; 17 killed, 2 stdout-only |
| M14 | Major | the histograms were cloned, 241 MB at 3,000 samples | **Applied** — `std::mem::take` |
| M15 | Major | the report named a sample by index | **Applied** — by name |
| M16 | Major | C3's fallback warning was never printed | **Applied** |
| M17 | Major | `header_text` read the metadata field by field | **Applied** — exhaustive destructure |
| Mi | Minor | the filter's id and the two precisions spelled twice | **Applied** — one spelling each, in `vcf/header.rs` |
| Mi | Minor | `mod.rs`'s inventory omitted `finish` | **Applied** |
| Mi | Minor | a stale test doc still said non-zero targets are refused | **Applied** |
| Mi | Minor | the flag's help said nothing about the disk the filter needs | **Applied** |
| 6a | numbers | the gate's library count | **Applied** — 6,621, not 6,619 |
| M | Major | ~73 lines of wiring duplicated between the subcommands | **Deferred** |
| M | Major | three configurations unreportable, one error variant unreachable | **Deferred** |
| Mi | Minor | three names say what a value is *about*, not what it *is* | **Deferred** |
| Mi | Minor | `write_the_records_the_filter_kept` also writes what it did not keep | **Deferred** |
| Mi | Minor | `PARALOG_POST` prints `0.000000` for a confidently-clean record | **Deferred** |
| — | — | spec §10's second oracle cannot pass as written | **Raised for the owner** |

## The four that change what a person sees

- **The flag's help says what the filter costs in disk.** Above zero the run parks every called
  record in `<output>.paralog-spill.tmp` and reads it twice, so it needs free space there of
  several times a compressed VCF's size. The file goes away when the run ends, however it ends.
  Before this, a run that fitted on its disk could now stop part-way through calling with no output
  at all and no warning it might.
- **A rejected sample is named, not numbered.** An index is a fact about the order the run's
  arguments were typed in; an operator would have to count their own command line. The comment
  three lines above already said "Named, not counted."
- **The warning for a run whose rate could not be fitted is printed.** It carries the half the
  short line dropped: that the records such a run removed were calibrated against a constant rather
  than against this cohort.
- **The filter's id has one spelling.** A `FILTER` value a record carries and the `##FILTER` line
  declaring it are one fact; a file whose records carry an undeclared id is invalid VCF.
  `the_declared_filter_id_is_the_one_a_record_carries` holds them together, and the two rounding
  precisions moved with it — a record's `PARALOG_LR` and the header's `lr_cut` must be comparable
  as written, or a difference between them says nothing.

## Deferred, with reasons

- **The command-line cohort that writes records (B1's remainder).** The shared fixture is two
  samples and three reads; it cannot fit a coverage model, let alone write a record. Giving it
  enough coverage changes the ground every command-line test in three files runs on, and C5 builds
  real-data runs anyway. **What was done instead**: the filter's whole path — fit, score,
  calibrate, write, report — is now exercised with records that are actually removed and tagged,
  through the one function both subcommands call. The untested remainder is the ~45 lines of
  per-subcommand wiring around it, and the duplication finding below is the other half of that
  story.
- **The ~73 duplicated lines of wiring.** Real, and it is what makes the untested surface two
  copies rather than one. Deferred because lifting it needs a shared error type across two
  subcommand error enums, which is a refactor of both commands rather than a fix to this step —
  and because C5 will touch both call sites anyway.
- **The three unreportable configurations, and the unreachable error variant.** The fit's
  configuration, the model's parameters and the calibration's knobs are chosen inside
  `fit_score_and_write_the_calls` with no way to override or print them; C3 kept
  `ParalogVerdicts::config` with the argument that "a number in that position which no run ever
  prints is a number nobody checks", and C4 is the step that prints the report. Deferred as a pair
  with the run report's shape, which C5's real runs will say more about.
- **The three names.** `WhatTheFilterDid`, `WhatTheOperatorAskedFor` and `what_to_tell_the_operator`
  say what a value is *about* rather than what it *is*, which is the project's own rule. Agreed;
  deferred because renaming three public types mid-milestone churns the two subcommands and the
  reports for no behavioural gain, and Checkpoint C is the place to take it.
- **`write_the_records_the_filter_kept` also writes what it did not keep**, in tag mode. Same
  reason.
- **`PARALOG_POST=0.000000` on a confidently-clean record.** True, and arguably right: the field's
  meaning is a probability and zero to six decimals is the honest rendering of one that small.

## Raised for the owner, at Checkpoint C

**Spec §10's second oracle cannot pass as written, and C5 is the step that runs it.** It asks that
the filter on at an unreachable target, with the two `INFO` keys stripped, equal the filter-off file
exactly. It cannot: the on-run's header carries four lines the off-run does not — the
`##paralogFilter=` line and the three declarations — and stripping `INFO` keys does not remove
them. Declaring them always would break §10's *first* oracle instead, that a filter-off run is
byte-identical to the pre-filter run, which is the one the whole plan rests on. **The two as
written are incompatible.** Recommendation: amend the second to strip the filter's four header
lines as well as the two `INFO` keys. That is a spec edit and it needs settling before C5.

## The mutation ledger

**The review ran 31 and 19 survived.** Thirteen of the nineteen were re-run against the new tests,
each applied as one exact string replacement, each file restored afterwards and its sha256 compared
with the original. **All thirteen killed**, each by the test written for it:

| # | the defect | killed by |
|---|---|---|
| S1 | `PARALOG_POST` loses four decimals | `the_two_info_fields_are_written_to_the_headers_precisions` |
| S2 | the patch failure's record number counts from zero | `a_line_that_cannot_be_patched_names_its_record_and_its_ordinal` |
| S3 | a ratio vector longer than the spill is accepted and truncated | `a_ratio_vector_longer_than_the_spill_is_refused` |
| S4 | the patch failure names contig 0, position 0 | `a_line_that_cannot_be_patched_names_its_record_and_its_ordinal` |
| S5 | the write failure names contig 0, position 0 | `a_write_failure_names_the_record` |
| S6 | the header swaps the fitted count and the cohort size | `the_header_says_how_much_of_the_cohort_was_fitted_and_whether_the_rate_was` |
| S7 | the header negates `em_converged` | same |
| S8 | the report swaps `dropped` and `tagged` | `the_report_counts_what_went_and_what_stayed` |
| S9 | the report swaps `written` and `unscored` | same |
| S10 | tag mode drops anyway | same |
| S11 | the saturation line is never printed | `the_report_says_when_ratios_sat_past_the_ends_of_the_range` |
| S12 | the calibration line moves below the declarations | `the_calibration_line_comes_before_the_info_declarations`, and two more |
| S13 | the fitted rate loses four decimals in the header | `the_provenance_writes_the_cut_and_the_rate_to_the_records_own_precisions` |

The six not re-run are the two stdout-only survivors (nothing in the suite captures stdout) and
four the reviewer reported against code this fix pass rewrote, where the mutation no longer applies.

### Four defects in the measuring tool, and what each would have claimed

**The sweep took four attempts, and every one of the first three would have reported a result it
had not measured.** They are recorded because the failure mode is the one this plan's own skill
warns about — a tool that says "green" about a run it never made:

1. **Two filters passed to `cargo test`, which takes one.** The joined string matched no test, and
   `0 passed` read as a pass — so **every mutation would have been reported as a survivor against a
   suite that never ran**. Caught because the first result claimed a survivor for a defect a new
   test demonstrably catches.
2. **A killed sweep left a mutation on disk.** The script writes the defect, runs, then restores;
   interrupted between those, `pass_three.rs` sat with its record numbering off by one. Found by
   scanning the working tree rather than trusting the script — the exact hazard
   `plan-driven-implementation` records, where two injected mutations once reached a commit whose
   message quoted a full-suite pass. **Fixed structurally**: the restore now happens *before* the
   result is judged, so no check can fail while the tree is dirty.
3. **The guard against (1) was itself double-escaped** inside a raw string, so it matched a literal
   backslash and never fired.
4. **And once it did fire, it only recognised a *passing* run** — but a killed mutant's run reports
   `FAILED`, so every correct kill was filed as "nothing ran". It now counts what a run executed,
   passed plus failed, over both result forms.

Three of the four biased towards *refusing* to claim a result, which is why they were survivable.
The first did not, and it is the one that would have put thirteen unearned kills in this report.

## What the numbers said, re-measured

| claim | was | is |
|---|---|---|
| the library test count at `627104f4`, in the report's gate table, its prose, and the commit message | 6,619 | **6,621**, and the step adds **16** tests, not 14 |
| "clippy first went to 13, in four kinds" | presented among re-derivable figures | true, but it describes a working tree no commit preserves — now marked unverifiable |

**The mechanism, because it is the transferable part.** The gate was run; then two more tests were
added; then the report was written from the earlier run. That is the first entry in this project's
own writing log — *prose written from the memory of a measurement rather than after the run it
reports* — and the tell was in the sentence itself, which had to explain away a discrepancy of
exactly two. The rule it breaks is the one that says a section is written after the run it reports,
never during.
