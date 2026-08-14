# B2 — a census that is a file, read one section at a time

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step B2 — implementation plan 2.
**Design authority:** [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§1.1a; [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§6.2.
**Date:** 2026-08-14.

---

## 1. What landed

`Sections` has its second state. `open_census(path)` checks the header and reads the directory —
**and decodes nothing else**; a scoped call then fills the sections it was asked for, hands the
closure borrows of them, and drops them when it returns. Nothing is retained between calls,
because there is no field to retain it in.

**The estimator did not change.** `fit_jointly`, `fit_contamination` and `gather_strata` take the
cohort and call `with_generic` / `with_strata`, which now answer from memory or from a file
without the caller knowing which. That is what milestone A was for.

## 2. The one deviation, and the number behind it

**`Sections::Backed` holds a path, where the architecture holds a `Box<dyn ReadSeek>`.**

An open reader means **one open file descriptor a sample, held for the length of the fit**. This
caller commits to cohorts of several thousand (`design_principles.md` §0), and the default soft
descriptor limit is 256 on macOS and 1,024 on Linux — before the pileups beside them are counted.
A thousand-sample cohort would fail to open its own evidence.

A path opens the file **once a call**, seeks to each section inside it, and closes it again. The
contract §1.1a states is unchanged in every particular — one read a section, decoded from a slice,
nothing retained between calls — and what it buys back is that a census can be cloned and compared,
which the estimator's own tests and three examples need. The door §1.1a wanted left open, memory
mapping, is if anything wider from a path than from a reader.

## 3. The scoped calls now return a `Result`

Recorded at A2 as the churn B2 would land: nothing can fail while every section is resident, so
the calls returned their value directly until there was a file to fail on. They now return
`Result<R, CensusError>`, and `JointFitError` gains a `Census` variant so the estimator's own
error type carries it.

Every call site inside the estimator propagates with `?`. The four examples all hold a census a
walk just built, so they `.expect` on one named constant that says exactly that.

## 4. Tests

| test | what it pins |
|---|---|
| `a_census_read_from_a_file_answers_what_the_one_in_memory_answers` | **the oracle.** Every section of the corner fixture, asked for through the same scoped calls off a file and off the value in memory, plus the terms, the read groups and the strata |
| `a_call_for_one_stratum_reads_that_stratum_and_leaves_the_rest_alone` | a band of one stratum out of two comes back with one section, and that section's extent is smaller than the file |

**What is not asserted yet, and it is the half worth writing** (spec §7.15): that the *bytes read*
are the section's own, which needs a counting reader. That is B4, where the same cohort is fitted
from memory and from files.

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo clippy --lib --all-features -- -D warnings` | clean; the examples this step touched are clean too |
| `cargo test --lib ng::parameter_estimation::joint::census_file` | `10 passed; 0 failed` (8 before) |
| `cargo test --lib` | `3,602 passed; 0 failed; 11 ignored` (3,600 before) |
| the 88-second tomato oracle | byte-identical to B1's |
| the 74-second trio oracle | byte-identical to B1's |
