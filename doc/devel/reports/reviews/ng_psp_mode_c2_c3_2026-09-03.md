# Code Review: ng_psp_mode_c2_c3
**Date:** 2026-09-03
**Reviewer:** rust-code-review skill (orchestrator; three sub-agents in isolated worktrees, seven category checklists between them)
**Scope:** steps C2+C3's uncommitted diff — the walk report and the overwrite refusal, on 7a164243
**Status:** Request-changes (all applied — see the fix report)

---

### 1. Scope

- Reviewed: the working-tree diff of plan steps C2 and C3, exported as `tmp/review_2026-09-03_c2_c3_report/c2c3.patch`.
- Against: commit `7a164243` + the patch, branch `ng-psp-mode`.
- In-scope: the C2/C3 additions to [generate_psps.rs](../../../src/pop_var_caller_exp/generate_psps.rs) and its [tests](../../../src/pop_var_caller_exp/generate_psps/tests.rs) — `WalkReport`, `SampleWalkOutcome`, `walk_every_sample`, `--force`, `PspAlreadyThere`, and seven new tests.
- Out of scope: the rest of the command (reviewed at C1), `gatherer.rs`, `run_ground.rs`.
- Categories (7): reliability, errors, extras, defaults, naming, idiomatic, smells. Skipped: module_structure (no module change), refactor_safety and unsafe_concurrency (no new type crossing a boundary, no concurrency primitive), tooling (no manifest change).

### 2. Verdict

**Request-changes.** The two features work — the refusal-before-any-walking property is genuinely pinned, and `--force` is well covered — but **20 mutations were run and 15 survived**, and the review found three claims in this diff that are wrong rather than merely untested: the report calls loci "observations", divides by a total that is not its numerator's whole, and its own doc says the report is a value a test can hold while a bare `println!` prints half of it where no test can see.

### 3. Execution status

- All three agents detached at `7a164243`, applied the patch, verified the marker symbols, and restored their trees (verified by diffing against the patch).
- Mutation totals: **20 run, 15 survived, 0 changed-no-behaviour** (reliability/errors 7/5, extras/defaults 6/5, naming/idiomatic/smells: static, with two rewrites verified by build).
- Orchestrator-side: fmt clean, clippy `-D warnings` exit 0, `cargo test --lib 'pop_var_caller_exp'` 119 passed; and a real three-run sequence on a tomato slice (write, refuse, force).

### 4. Open questions and assumptions

1. **Should the per-sample progress line exist at all, given the report repeats it?** Kept — the plan asks for both — but moved to stderr and made to share one formatter with the report, so the two cannot drift and a shell capturing the report gets only the report.
2. **How much of `LocusCounts` does a walk owe a reader?** Five of eight fields, unchanged: the region counts and the three base counters. `loci_emitted` equals the record count the line already prints, and the two unhandled *region* counts say less than their base counterparts.

### 5. Top 3 priorities

1. **M1/M2** — the two wrong claims a person reads: loci reported as "observations", and shares taken over the ground *asked for* rather than the ground *walked*, which is the exact arithmetic the sibling report records printing 200.0% once.
2. **M3** — C3's headline property is not pinned: the truncation test tests `PspReader`, and the stopped-walk test's new assertion cannot fail.
3. **M4** — C3's own guard silently defanged a C1 test, which now stops at the door and never reaches the walk it claims to exercise.

### 6. Findings

#### Major

- `generate_psps.rs:302,345` — **M1: "observations" names the wrong thing, by a level of nesting.** `WriteStats::records` counts psp records, and one record *is* a `SampleLocusObservations` — a locus, which itself holds a field named `observations: Vec<SequenceObservation>` documented as "Observations, not alleles". So "193,603 observations" is really 193,603 loci holding several times that many observations, in a report whose reader is a geneticist. **Categories:** naming, extras.
- `generate_psps.rs:345-350` — **M2: the share's denominator is not its numerator's whole.** A repeat tract is typed and walked whole even where a BED cuts one, so `regions_handled_bp` and its two siblings can sum past the BED's ask. `src/ng/run/report.rs:160-212` documents this same bug being found and fixed in the sibling command — "a BED of 120 bases inside two tracts charged 240 to *not built yet*, and dividing by the 120 printed **200.0%**" — with a committed regression test. Both commands reach these counters the same way. **Categories:** naming, reliability, extras (three agents).
- `generate_psps.rs:539-546` — **M3: the C2 split does not do what its own doc claims.** `walk_every_sample` exists so the report is "a value a test can hold, rather than something only a terminal ever sees" — and then prints a bare `println!` inside itself, in different words from `lines()`, which no test covers and whose deletion survives. **Categories:** naming, idiomatic, smells, reliability, extras (all five agents that could see it).
- `tests.rs` — **M4: C3's guard silently defanged a C1 test.** `a_stopped_rewalk_does_not_destroy_the_psp_it_was_replacing` now stops at the new pre-flight check (`PspAlreadyThere`) and never opens an alignment file. It still passes because `expect_err` does not check the variant, and its own message — "the re-walk's file has no index" — is now false. **Categories:** extras.
- `tests.rs` — **M5: C3's headline property is not pinned.** "An interrupted run leaves nothing that reads as whole" is tested by truncating a psp by hand and asserting `PspReader` refuses it — a property of the reader, which the plan says the format already guarantees, and a state this command cannot produce now that it renames. The stopped-walk test's `.partial` assertion cannot fail either: the sample stops inside `open`, before any `.partial` exists, which is why deleting the cleanup left the suite green. **Categories:** reliability, extras.
- `generate_psps.rs` — **M6: the overwrite guard is advisory, and the scratch name makes the race worse.** Two invocations naming one sample — the parallelism this module's own doc recommends — both pass the existence check, both write the same fixed `<sample>.psp.partial`, and both rename over the final path. **Categories:** errors, extras.
- `generate_psps.rs` — **M7: the report's content is barely pinned.** Of the eight numbers `lines()` prints, three are asserted; mutations that swapped the two region counts, blanked the denominator, or deleted the uncovered-ground clause entirely all survived. **Categories:** reliability, extras.

#### Minor

`psp.exists()` answering "no" both when there is no file and when the process cannot tell (errors); the refusal naming the first sample rather than the blocked one, mutation-proven (reliability); `--force`'s flag name untied to the `--force` in its own message, so renaming the argument leaves the message wrong and the suite green (reliability); the two unhandled kinds rendered as their internal names ("not built yet", "out of scope") where the sibling report deliberately renders the same fields in plain English (defaults); nothing saying a psp was *replaced* when `--force` acts (defaults); the report never naming the ground, which `run/report.rs` records being fixed once already (smells); `WalkReport`/`SampleWalkOutcome` `pub` with a private constructor and no `Debug` (smells); `SampleWalkOutcome`'s field names diverging from the crate's existing `SampleWalkTallies` for the same concepts (naming); the `--force` test building a second fixture cohort and dropping its `TempDir` guards on the same line (smells); `push_str(&format!(…))` where `write!` is the idiom, with a comment prescribing the wrong cure (idiomatic).

### 7. Out of scope observations

- Nothing from Milestone D or later crept in: no non-test code in the diff reads a psp.
- The command's `--help` and the refusal message were both verified against the rendered binary.

### 8. Missing tests to add now

The report's numbers pinned one by one; the two uncovered-ground clauses present when non-zero and absent when zero; the refusal naming the blocked sample rather than the first; the `--force` flag tied to clap; the stale-`.partial` cleanup made observable by planting one; and the defanged C1 test given `force: true` so it reaches the walk again.

### 9. What's good

- Two properties are genuinely pinned by mutation: moving the overwrite check into the walk loop fails `the_refusal_comes_before_the_first_sample_is_walked`, and ignoring `--force` fails its own test.
- `--force` passes almost all of the defaults checklist: read in the open at the call site, its consequence in `--help`, a refusal naming sample, path and remedy, and the invocation recorded in the psp header.
- The diff is C2+C3's scope with nothing from later milestones.

### 10. Commands to re-verify

`scripts/dev.sh cargo fmt --check`; `… clippy --all-targets --all-features -- -D warnings`; `… cargo test --lib 'pop_var_caller_exp'`; and the three-run real sequence: write, refuse without `--force`, replace with it.

### Author response convention
Address findings by identifier (M1–M7) in the fix-application report.
