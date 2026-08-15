# A2 — `records.rs` → `census.rs`, and the type names

**Plan:** [census_rename_and_encoding.md](../../ng/impl_plan/census_rename_and_encoding.md), milestone A,
step A2.
**Design authority:** [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§4's rename table.
**Date:** 2026-08-14.

---

## 1. Plan

One mechanical sweep, no arithmetic touched. The module moves with `git mv`; the type names are
substituted with word-boundary `sed` inside the dev container; the prose the substitution damaged is
repaired by hand; the compiler and the library suite are the oracle.

## 2. Assumptions and absorbed deviations

Four, each recorded because the plan's step text does not name it.

1. **`RecordError` → `CensusError` is a no-op.** The plan lists it and the arch table lists it, but no
   such type exists in the code: `rg -w RecordError src examples tests benches` returned nothing
   before the sweep. The error enum is specified (arch §2.3) and belongs to the census *file*, which
   is plan 2. Nothing was created.

2. **`KeptLociDigester` → `CensusLociDigester` was renamed too**, although neither the plan's list nor
   the arch table names it. It is the builder of `CensusLociDigest`; leaving it behind would leave the
   module saying both vocabularies, which is the state the arch says this rename exists to end.

3. **The field `SampleCensusEvidence::identity` was renamed to `terms`**, and with it the bindings and
   test helpers that carried a `SelectionTerms` or a `RecordingTerms` value under the name `identity`.
   The arch's own type sketch (§1.1) names the field `terms`; `RecordingTerms` stored in a field called
   `identity` says the thing the rename removed. The compiler checks every site.

4. **Function names holding the old vocabulary were left alone** — `select_kept_loci`, and
   `CensusLoci::from_parts`'s callers. *"Select the kept loci"* is still true English about what the
   function does, and renaming functions is not what the step names.

## 3. Changes made

| what | where |
|---|---|
| module moved | `src/ng/parameter_estimation/joint/records.rs` → `census.rs` (`git mv`), `mod.rs` updated |
| `KeptLoci` → `CensusLoci` | `loci.rs` and its users |
| `KeptLociDigest` → `CensusLociDigest`, `KeptLociDigester` → `CensusLociDigester` | `loci.rs`, `census.rs`, three examples |
| `SampleRecords` → `SampleCensusEvidence` | `census.rs`, `fit.rs`, `ssr_fit.rs`, `contamination.rs`, four examples |
| `GenericRecords` → `GenericEvidence`, `SsrRecords` → `SsrEvidence` | `census.rs` and its users |
| `RecordWriter` → `CensusWriter` | `census.rs`, `examples/ng_joint_records_walk.rs` |
| `RecordIdentity` → `RecordingTerms`, `SelectionIdentity` → `SelectionTerms` | `census.rs`, `loci.rs`, `fit.rs`, `contamination.rs`, four examples |
| field `identity` → `terms` | `SampleCensusEvidence`, `CensusWriter`, and every reader |
| module doc rewritten to say what a census is | `census.rs` header |

**One line of harness output changes text:** `examples/ng_joint_records_walk.rs` printed
`--- the identity check, which is what lets these be pooled ---` and now prints
`--- the recording-terms check, … ---`. It is named here because milestone A's oracle is that the two
cohorts' output does not move, and this is one of the places it does.

## 4. Tests

No test was added and none was removed: an identifier substitution has no new behaviour to assert, and
the suite's job here is to prove that none of the old behaviour moved. Four test helpers were renamed
with their types (`identity()` → `terms()`, `selection_identity()` → `selection_terms()`).

## 5. Validation

Run in the dev container.

| command | result |
|---|---|
| `cargo fmt --check` | clean (after `cargo fmt` reordered five import blocks the substitution had left out of alphabetical order) |
| `cargo check --all-targets` | clean |
| `cargo test --lib` | `3581 passed; 0 failed; 11 ignored`, 469.33 s — the same 3,581 the plan's preconditions record |
| `cargo test --all-targets` | every test target green — 3,581 in the library plus 61 across the integration targets — and then **one benchmark target panics**: `benches/psp_writer_perf.rs:386` — `index out of bounds: the len is 3300000 but the index is 3300000`. **The same panic, same file, same line, on `ng-joint-fit` before this branch's first commit.** It is in the `.psp` writer's benchmark, which this plan does not touch. Not chased. |
| `cargo clippy --all-targets --all-features -- -D warnings` | red, **and red identically on `ng-joint-fit` before this branch's first commit** — three errors in `lib test` (`this function has too many arguments (9/7)`, `useless use of vec!` twice) plus errors in `ng_duplicated_class_harness` / `ng_joint_fit_harness`. None is in a line this step changed. Not chased. |

`rg -w` over `src`, `examples`, `tests` and `benches` for the ten old names returns nothing.

## 6. Tradeoffs and follow-ups

- **The coverage-by-window summary is still here**, `CoverageByWindow` and `coverage.rs` and the
  `coverage` field and `RecordingTerms::coverage_window` with it. It is A3's deletion, deliberately
  not bundled: A2 must be provably answer-preserving and a deletion that removes printed lines is not.
- **`SampleCensusEvidence` still holds its sections as two public `BTreeMap`s** rather than the
  private `Sections` of arch §1.1. That shape belongs to the census file (plan 2) and nothing here
  touches it.
