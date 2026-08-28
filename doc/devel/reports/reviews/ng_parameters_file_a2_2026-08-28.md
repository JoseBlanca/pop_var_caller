# Code Review: ng parameters file — A2, a value and its warrant as one shape

**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator) over two category agents in isolated worktrees
**Scope:** the uncommitted working-tree change of step A2 of [parameters_file.md](../../ng/impl_plan/parameters_file.md)
**Status:** Request-changes (all applied — see [the fixes report](fixes_applied_ng_parameters_file_a2_2026-08-28.md))

---

## 1. Scope

- **What was reviewed:** the working-tree diff of step A2 — one new type, `WarrantedValue`, and
  four numbers moved into it.
- **Reviewed against:** branch `ng-parameters-file`, uncommitted, on top of A1 at `2ac36b53`. Both
  agents detached to that commit and applied the diff as a patch.
- **In-scope files:** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)
  and its `testdata/every_shape.toml`, **and only the parts the patch touches** — A1's other
  thirty types were reviewed and settled in the previous step.
- **Categories dispatched:** `reliability` and `naming`. Two rather than A1's five, because A2 is
  one type and one idea on a module the previous step swept in full. `idiomatic`, `smells` and
  `module_structure` were not dispatched: A2 adds no module, no error path, no concurrency and no
  duplication, and the shape question it does raise — one shape against per-row spellings — is
  what `naming` and `reliability` were pointed at directly.

## 2. Verdict

**Request-changes.** Two Majors from each agent, and the three that matter are all about claims
the new code makes rather than about what it does.

## 3. Execution status

Run by the orchestrator in the container, on the reviewed tree, and passed to both agents:
`cargo fmt -- --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib ng::calling::parameters_file` 11 passed, 0 failed, 1 ignored; `cargo doc
--no-deps` zero unresolved links in this module.

**Both agents kept their worktrees this time** — the fix for A1's three lost reports was to tell
each to write its findings file early and update it as it went, and both did. `reliability` ran
**4 mutations, 0 survived, 0 changed no behaviour**, and verified its restore by content: the
working-tree diff over the two in-scope files came back byte-identical to the patch it was given.

Findings labelled "Needs verification": none.

## 4. Open questions and assumptions

1. **Whether editing `repeat_tract_outlier_weight` in the file makes it `Supplied`.** The fix
   applied assumes it does, on spec §2's ladder and spec §3.8's statement that a person editing it
   is why it is written down. If the owner's view is that this number is *always* `Defaulted`, the
   field's doc should say so and the number should join the list of quantities that carry a
   different kind of warrant. Recorded rather than assumed silently.

## 5. Top 3 priorities

1. **A number with a four-state warrant available was left out** (`reliability`, Major). Spec §8
   names three parameters with an honest default that must all be marked `Defaulted`; A2 moved
   two of the three and left `repeat_tract_outlier_weight` a bare `f64`.
2. **Both new statements of what the evidence count counts are wrong** (`reliability`, Major), and
   A2 had just removed the unit from the key, making those comments the only place a reader could
   learn it. An inbreeding coefficient counts covered positions, not windows — and a window is
   100,000 bases. A repeat-tract substitution rate counts bases compared, not reads.
3. **A test assertion pinned a claim the spec contradicts** (`reliability`, Major). `assert!(
   calibration.get("observations").is_none())` fixes "the calibration multiplier never carries a
   count"; spec §3.3 asks for that multiplier "with its warrant and the number of observations
   behind it", and the count exists upstream on the `Estimate<ErrorRate>` the seam drops. When B1
   supplies it, the test goes red and reads as *put it back to None*.

## 6. Findings

**4 Majors and 4 Minors**, plus two nits and four cross-category notes. The per-finding text is in
`tmp/review_2026-08-28_parameters-file-a2/`, which is gitignored;
[the fixes report](fixes_applied_ng_parameters_file_a2_2026-08-28.md) accounts for every one.

### Major

- **M1** `repeat_tract_outlier_weight` is one of spec §8's three honest defaults and is the only
  one that did not move into the shape. *(reliability)*
- **M2** both new statements of what `observations` counts name the wrong unit, and the key is
  deliberately unnamed so those comments are the only place a reader can learn it.
  *(reliability, convergent: `naming`'s Minor on the same key)*
- **M3** `every_warranted_number_is_written_the_same_way` pins "the calibration multiplier never
  has an observation count", which spec §3.3 contradicts. *(reliability)*
- **M4** the new doc's heading and its slippage bullet call four non-`Warrant` mechanisms "a
  different kind of warrant", contradicting the reservation `Warrant`'s own doc makes
  thirty-five lines above — the same overload A1's review filed as its M13, re-opened in prose
  rather than in identifiers. *(naming)*

### Minor

- **m1** the wrong seam is cited: `RunParameters::assemble` strips the inbreeding warrant, not
  `joint/census_moments.rs` — which is in fact the one cited file that *does* keep a count.
  *(reliability)*
- **m2** `observations` names a count with neither its unit nor its subject, and the shared-shape
  argument does not carry: two of the four quantities have no count, and `ContaminationRow` one
  section away spells its evidence counts with their units. *(naming)*
- **m3** the produced file no longer contains a `[repeat_tracts]` heading at all, because every
  field of that section became a table. *(naming)*
- **m4** three test-quality points: the `7` counts fixture *rows* rather than warranted *fields*,
  so a balanced change leaves it satisfied; the test's doc credits the "written flat" mechanism
  with a case the count is actually catching; and
  `a_stated_concentration_of_one_says_whether_it_was_fitted` is near-unfalsifiable after A2,
  because `warrant` is a required field and two files whose warrants differ differ by
  construction. *(reliability)*

### Nits

Two from `naming`: the ordinary-site prior's "**pair**" uses a term before anything says it is two
concentrations; and "one **read** off a curve through four cells" puts the sequencing sense of
*read* and the verb sense in one clause, so it parses as one sequencing read on the first pass.

## 7. Out of scope observations

`naming` confirmed A1's M13 is **not** re-introduced at the identifier level: `warrant` appears as
an identifier only on `Warrant`, `WarrantedValue` and its field, and the three structures M13 named
still carry `origin`.

## 8. Missing tests to add now

Three, from `reliability`, all applied: the evidence count's unit by name, the two absences, and
the reader-side recovery of a stated concentration's warrant.

## 9. What's good

- `every_warranted_number_is_written_the_same_way` **is discriminating**: all four mutations the
  brief named were killed by it, including flattening a row and adding a fifth warranted number.
- The list of quantities that carry a different kind of warrant was checked against spec §2, §3 and
  §5 and against `run_parameters.rs`, and found correct and complete for the four it named.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file -- --ignored regenerate`
