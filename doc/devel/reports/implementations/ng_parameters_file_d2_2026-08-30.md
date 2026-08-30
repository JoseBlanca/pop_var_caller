# ng parameters file — D2: the three refusals

**Date:** 2026-08-30
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone D, step D2
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §6, §13 test 4
**Code:** `refuse_if_not_this_runs_inputs` in
[bindings.rs](../../../../src/ng/calling/parameters_file/bindings.rs), and
`ParametersFileError::FittedFromOtherInputs` in
[mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. What it does

Spec §6 binds the file to four things. Three of them **refuse** a run they do not match, because a
file that fails one cannot be *interpreted* against that run: one fitted against another assembly
has its repeat strata cut at other tract lengths, one listing other samples has its inbreeding
coefficients against other plants, and one whose read-group table does not cover the run's leaves a
library with no calibration and no contamination row — which surfaces as a panic at whichever locus
first carries one of that library's reads. The fourth, the census, demotes instead, and is D3's.

`refuse_if_not_this_runs_inputs(&self, reference: &ReferenceDigest, read_groups: &ReadGroups)`
takes **the same two arguments `of_run` writes from**, so the writer and the check read one pair of
inputs and a file this run wrote is a file this run accepts. The sample list is not a third
argument: it is derived from the read-group table exactly as `of_run` derives it.

## 2. The refusal's shape — the owner's ruling, and what it buys

**Owner's ruling of 2026-08-30: exceed the census.** The census's own refusal is
`Freshness::{Rebuild, Refused}(&'static str)` — a field name, with the two values `freshness`
compared already discarded. Spec §9 says this refusal is "in the shape the census's own refusal
already uses", and that clause is what is wrong: §6 and §13's fourth test both ask for the field
*and* the two values, and a bare field name cannot answer *which of my three copies of this
assembly was the file fitted on*.

So `FittedFromOtherInputs { field, in_the_file, in_the_run }`, rendering as

> the parameters file was fitted from other inputs: `fitted_from.samples[1]` is `"Ailsa ‘Craig’
> \"×2\""` in the file and `"TS-9"` in this run

**`field` is a key path and not prose**, the same vocabulary `Meaningless` already uses, so a
reader meets one way of being pointed at a line and can find it by searching the file in front of
them. `every_refusal_names_a_key_the_file_contains` holds it, over all four shapes.

## 3. The two axes are compared differently, and neither choice is free

**Samples by position.** Every per-sample row is read *by name* into the file's list and handed to
calling as a *position*, so a list holding the right names in another order is a run that gives
each plant its neighbour's inbreeding coefficient and batch, with nothing downstream able to see
it. `the_same_samples_in_another_order_are_refused`.

**Read groups by their number.** Row order in `fitted_from.read_groups` carries no meaning anywhere
in the module — `validate` sorts the ids before checking they are dense, the projection reads the
table only for its length, and every other section joins on the `read_group` key. **The first draft
compared these two tables positionally and refused a legitimate file** (§6).

## 4. What it does not check, and the two holes that leaves

It checks nothing about the file's agreement with itself — that is `validate`'s, and this does
**not** call it. A file that has not been through one is still compared soundly: every join here is
a lookup rather than a position, and there is no indexing and no panic path in the function.

**⚑ Two things do follow from that, and both are why D3 runs the two together.** A file whose own
sample list disagrees with its own read-group table is refused here *as though the run were wrong*,
where `validate` names the file precisely; and a file with no samples and no read groups matches a
run with none, which `validate` refuses and this cannot see. Neither is reachable once `validate`
has run first, and **D3's own test is where that ordering is pinned**.

## 5. Tests

Thirteen added. Beyond one per binding:

| test | what it holds |
|---|---|
| `the_run_that_wrote_the_file_is_not_refused` | the case every refusal is a departure from — a suite of refusal tests alone passes on a function that always refuses |
| `the_same_samples_in_another_order_are_refused` | the positional axis, with both names asserted |
| `a_cohort_with_one_more_sample_is_refused_by_its_count` | the half `zip` cannot reach: a cohort agreeing on every shared position |
| `a_read_group_row_that_differs_in_any_name_is_refused` | all three names, **and both values on each** — see §6 |
| `a_read_group_table_numbered_otherwise_is_refused` | the same lanes filed under numbers the run does not use |
| `a_file_whose_rows_are_written_in_another_order_is_the_same_file` | the Blocker's regression: it validates, and it is accepted |
| `the_reference_is_the_first_thing_checked` | §6's order, which decides what a run mismatched twice hears |
| `every_refusal_names_a_key_the_file_contains` | all four shapes, every dotted segment, against the produced text |
| `a_refusal_over_a_cohort_of_thousands_still_fits_on_a_line` | 3,000 samples — §9's top of range — in **67 characters**, asserted rather than described |
| `a_binding_refusal_has_no_line_in_the_file_to_send_anyone_to` | `line()` is `None`, and both values reach the message |

## 6. Validation

All in the container, by absolute path:

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --lib ng::calling::parameters_file`: **158 passed, 0 failed** (145 before).
- `cargo test --lib`: **5,529 passed, 0 failed** (5,516 before).
- `cargo doc --no-deps`: **25 unresolved-link errors**, the baseline.

## 7. Mutation testing

**Thirteen mutants, every one fails a test.** Each spliced from a pristine copy and the restore
`diff`-verified before the next.

| mutant | tests it fails |
|---|---|
| print no names at all, only the tally | 4 |
| delete the reference refusal | 3 |
| delete the sample-position loop | 3 |
| delete the read-group id-set check | 3 |
| check the reference after the samples | 3 |
| delete the sample-count check | 2 |
| drop the `declared_id` column | 2 |
| name read groups without their number | 2 |
| swap the two values at the reference site | 1 |
| swap the two values on a read-group row | 1 |
| drop the library column | 1 |
| drop the sample column | 1 |
| **join the two read-group tables by position again** | 1 |

## 8. What the review found

Two agents in isolated worktrees: correctness, and design fidelity plus the quality of the messages.

### The Blocker: a legitimate file was refused

The first draft joined the file's read-group table to the run's **by position**. Row order in that
table carries no meaning anywhere else in the module, so a file with two rows swapped is the same
file — and where the two rows are two lanes of one plant, the file's own first-seen sample order
does not move either, so it validates and projects. The reviewer ran it: `validate` said `Ok`, both
projections succeeded, and the binding check refused with *the dense index of read group 1 is 1 in
the file and 0 in this run*.

Joined on the read group's own number now, which also makes the gap check exactly *equal id sets*.
`a_file_whose_rows_are_written_in_another_order_is_the_same_file` is the regression, and the
mutant that restores the positional join fails it.

### Three messages a geneticist could not act on

- **`the dense index of read group 3 is 3 in the file and 2 in this run`** — circular: the field
  names the row by the number it then reports as wrong, there is no read group 3 in the run at all,
  and *dense index* is a phrase that appears nowhere in a produced file. Gone with the key join.
- **`the number of samples is 2 in the file and 3 in this run`** — §6 says *refuse, **naming the
  samples***, and this is the commonest real refusal there is: a cohort that gained an accession.
  A reader was handed two cardinalities and left to diff two lists by eye, one of which is not
  written down anywhere they can look. It names them now, capped at five with a tally.
- **Two identical lists beside "these differ"** — found by printing every message this code can
  produce, after the key join was in. A file numbered `0, 1, 3` for this run's three lanes has the
  same three `@RG ID`s as the run, so naming lanes by `@RG ID` alone printed
  `3 read groups ("HWI.3", "HWI.4", "HWI.5")` on both sides. Each lane is now its number **and**
  its name.

### Two doc comments that said something untrue

- **`field` was documented as "named so that a reader can find it in the file"** and two of its own
  three examples appeared nowhere in a produced file. Rather than weaken the promise, the fields
  became key paths, which makes it true and brings this error under the same guarantee
  `validate`'s refusals carry.
- **The check order was defended with a claim the file's comment does not make.** That comment
  contrasts refuse-with-demote, not walk order, and the property it names holds whatever the order
  is. What the order actually decides is which refusal a run mismatched in more than one place
  hears first, and the reference leads because its consequence is the worst of the three.

### The value-swap gap

The three name columns asserted only the *field*. The reviewer swapped `in_the_file` and
`in_the_run` in the library arm and the module stayed green — a message telling a geneticist the
file says `lib9` where the file says `lib4`. All three arms assert both values now, and the mutant
fails.

### Recorded and not fixed

- `validate` never checks that `reference_digest` is 32 lower-case hex characters, so a file with
  upper-case hex refuses with two 32-character strings differing only in case. A later `validate`
  rung, recorded in `PROJECT_STATUS.md`.
- No message names the *file*. A run is handed a path; the wrapping belongs at F1's call site.
- Names are rendered through Rust's `Debug` escaping, which diverges from TOML's for control
  characters — a name a reader then could not find by searching the file. Unlikely in an `@RG SM`
  and cheap to note.
