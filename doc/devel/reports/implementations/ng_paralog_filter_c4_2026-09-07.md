# The hidden-duplication filter — C4: the verdict on the file, and the run that carries it

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone C, step C4
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.5, §3.6
**Branch:** `ng-paralog-filter`

## The answer

**The filter runs.** A run reads its sources, parks every finished record, fits each sample's
coverage model, scores every parked record, resolves the operator's target false-discovery rate to
a cut, and writes the VCF — leaving out the records the cut removes, or writing them on the
`hiddenParalog` filter when the operator asked to see them. Both subcommands do this through one
call, so direct mode and psp mode cannot drift into filtering differently.

**And the two departures C2 recorded are reversed, which is what this step was told to do.**
`--paralog-fdr` defaults to spec §3.6's `0.01` again, and the refusal of any non-zero target is
deleted. A run that says nothing about the filter now filters; `--paralog-fdr 0` turns it off and
writes byte for byte what the run wrote before this work.

## What was added

| file | what it is |
|---|---|
| [`pass_three.rs`](../../../src/ng/run/paralog_filter/pass_three.rs) | the verdict on each parked record and the file that comes out; `WhatTheFilterDid`, `PassThreeError` |
| [`pass_three/tests.rs`](../../../src/ng/run/paralog_filter/pass_three/tests.rs) | 9 tests |
| [`finish.rs`](../../../src/ng/run/paralog_filter/finish.rs) | what a run does after its calling pass: the fit, both passes, the header's provenance, and the report's words |
| [`finish/tests.rs`](../../../src/ng/run/paralog_filter/finish/tests.rs) | 5 tests |
| [`vcf/header.rs`](../../../src/ng/vcf/header.rs) | `HiddenParalogProvenance`, the three declarations, and the `##paralogFilter=` line — all conditional |
| both subcommands | the sink chosen by the target, the filter run after the calling pass, the report's extra lines; the refusal and its error variant deleted |

## Assumptions and recorded deviations

### 1. The filter's header lines are written only when it ran

Spec §3.5 gives the three declarations and the `##paralogFilter=` line but does not say whether a
run with the filter off carries them. **It cannot**: the standing oracle for this whole plan is
that `--paralog-fdr 0` reproduces the pre-filter run byte for byte, and a declaration emitted
always would break that on the header alone. There is a second reason that would hold even without
the oracle — a file declaring a filter it never applied cannot be told apart, by a reader, from one
that applied it and found nothing.

So `VcfHeaderMetadata` gains one field, `None` by default, and
`the_hidden_duplication_filter_ran` is the only way to set it. It is consuming, so a header cannot
be built and then quietly amended.

**`VcfHeaderMetadata` lost its `Eq` derive** as a consequence: it now holds the fitted rate and the
cut, which are floating-point, and `Eq` promises a reflexivity floats do not have. Nothing used it.

### 2. `WhatTheOperatorAskedFor` bundles the two knobs

`fit_score_and_write_the_calls` reached eight arguments, and clippy said so. The two that belong
together are the operator's: a target and a bare `bool`, adjacent in the list and swappable past
the type checker. They are one decision — at a target of zero the filter does not run and the tag
flag has nothing to act on — so they travel as one value.

### 3. The two `INFO` fields go on a record with a finite ratio; the posterior can still be absent

Spec §3.5 says both fields ride "on every record with a finite ratio". `PARALOG_POST` needs the
run's fitted rate to be strictly between none and all, and the copied calibration answers `None`
where it is not, rather than a saturated `0` or `1`. So a scored record always carries the ratio
and *usually* the probability. **An earlier draft of this step keyed both fields on the posterior
being present**, which would have counted a scored record as unscored on a run whose rate fitted to
a certainty; caught before the tests were written.

### 4. The ratio is written to four decimals, the posterior to six

The same precisions the header's `lr_cut` and `pi` carry, so an operator auditing why a record went
can compare a record's `PARALOG_LR` against the header's `lr_cut` as written, without wondering
whether a difference is real or a rounding.

### 5. The sink is taken apart rather than consumed

`CalledRecordSink::finish` drops the spill, and passes two and three need it. Each subcommand
matches the sink out instead: the off variant finishes its writer, the on variant flushes its spill
and hands it back for the next two passes.

### 6. `print_report_with_the_filter_s_lines`, with an empty slice on an off run

One function rather than two, and the empty slice is what says an off run prints what it printed
before the filter existed.

## Tests

**16** — 9 of pass three, 5 of the run's finish, and 2 of the command line.

**The one that carries the oracle** is `a_record_the_verdict_does_not_touch_is_written_unchanged`:
the written line is compared against the parked entry's **own bytes**, not against a re-encoding,
because a re-encoding would agree with a writer that rebuilt the line wrongly in the same way. And
`tagging_writes_what_dropping_leaves_out_on_the_filter_s_id` is spec §10's third relation at the
level of two records: the tag-mode file with its tagged lines removed *is* the drop-mode file.

**Which dimensions the fixtures vary.** Whether a record is removed, kept or tagged; whether it was
scored at all; whether the operator asked to drop or to tag; the contig, at two of them, so a place
built from the wrong entry shows up as an ordering refusal rather than passing; the `FILTER`
column's prior contents, since a record already carrying one is joined and not replaced; and the
record count, at zero, one and three.

**A test's premise was wrong and the measurement said so.** `the_verdict_follows_the_file_s_order…`
first ran at a target of one in two, where two of its three records are removed rather than the one
the test is about. Measured, the three records' tail false-discovery values are `0.667`,
`1.17 × 10⁻⁵` and `0.500`; at a target of three in ten only the middle one is reached, and the third
is close enough that at one in two it goes too — so the fixture is not passing because the other
two are far away.

**And one gap the suite could not have seen.** Every existing command-line test sets
`--paralog-fdr` explicitly, so nothing would have noticed the default going back to `0.01` — or
failing to. Both subcommands now have `a_run_that_says_nothing_about_the_filter_takes_the_specs_target`,
which reads it off the parsed command line.

## The gate

| command | at this step | at C3's last commit `420b041f` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6621 passed; 0 failed; 15 ignored` | `6605 passed; 0 failed; 15 ignored` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9, same kinds |
| `cargo fmt --check` | dirty on 9 files, none this step's | the same 9 |

6,621 is 6,605 plus this step's 16 tests. The one integration failure is `main`'s.

**This figure was wrong in the first version of this report, and the way it went wrong is the one
this project's writing log names first.** It said 6,619 — the number from a gate run made *before*
the last two tests were added, and then written up from memory of that run rather than after a new
one. The two missing were the command-line tests two paragraphs above, which is why the sentence
also had to explain away a discrepancy of exactly two. Caught by the step's own review, re-measured
at `627104f4`, and corrected forward.

**Clippy first went to 13, in four kinds the baseline does not have** — a loop counter, two unused
imports, and the eight-argument function that became deviation 2. Counting error *kinds* is what
saw them; a check greping for the baseline's three could not have. **That figure describes a
working tree no commit preserves**, so unlike every other number here a reader cannot re-derive it;
it is kept because the lesson is the counting method, and marked because it sits among numbers that
can be checked.

## What this leaves for C5

The four oracles, on real reads: off is byte-identical; on at an unreachable target with the two
`INFO` keys stripped equals off; the drop-mode file equals the tag-mode file minus its
`hiddenParalog` lines; and both modes produce the same file with the filter on.
