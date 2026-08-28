# Code Review: ng parameters file — A3, absence is a missing key

**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator) over two agents in isolated worktrees
**Scope:** the uncommitted working-tree change of step A3, and — for the design-fidelity pass — the whole of Milestone A against spec §3 and §5
**Status:** Request-changes (all applied except two design questions raised at Checkpoint A — see [the fixes report](fixes_applied_ng_parameters_file_a3_2026-08-28.md))

---

## 1. Scope

- **What was reviewed:** the working-tree diff of A3 — `contamination` becoming `Option`, and a
  row's measurement becoming a nested `Option`. **The second agent's scope was wider on purpose**:
  Checkpoint A's own question is whether the types cover §3 and whether §5's five states are
  expressible and distinct, and that is a question about Milestone A as a whole.
- **Reviewed against:** branch `ng-parameters-file`, uncommitted, on top of A2 at `7ab5cf4f`.
- **Categories dispatched:** `reliability`, and a **design-fidelity** pass (the `smells` checklist
  pointed at spec conformance). Two rather than five: A3 adds no module, no error path and no
  concurrency, and the names it introduces are two, both settled by §5's own vocabulary.

## 2. Verdict

**Request-changes.** One Blocker, four Majors, four Minors. The Blocker is that the step's own
test could not fail on the hazard its own message names, and the reviewer proved it by writing
that hazard.

## 3. Execution status

Run by the orchestrator in the container and passed to both agents: `cargo fmt -- --check` clean;
`cargo clippy --all-targets --all-features -- -D warnings` clean; `cargo test --lib
ng::calling::parameters_file` 13 passed, 0 failed, 1 ignored; `cargo doc --no-deps` zero
unresolved links in this module.

`reliability` ran **5 mutations, 1 survived, 0 changed no behaviour**, restoring between each and
verifying its final tree by content against the patch it was given. The design-fidelity agent ran
five throwaway probes in the module's own test block — 18 passed, 0 failed — and wrote out the two
Rust values for each of §5's five rows.

Both agents kept their worktrees. Findings labelled "Needs verification": none.

## 4. Open questions and assumptions

Two, both raised at Checkpoint A rather than settled, because both are the spec and the code
disagreeing rather than the code disagreeing with itself:

1. **Spec §2.1's wholesale demotion has nowhere to write itself.** Only five numbers carry a
   `Warrant`; the slippage numbers, the contamination fraction and the ordinary-site prior's two
   concentrations carry a different kind, and none of those vocabularies has a word for
   `Supplied`. §13's fifth test — "same genotypes, every warrant `Supplied`" — cannot be written
   as stated. The in-memory `LevelProvenance` has the same gap.
2. **The repeat-tract substitution rate's axis is not priced.** It is keyed by (read group ×
   stratum × ploidy); §3.7 says "per stratum" and §9's three axes do not include it.

## 5. Top 3 priorities

1. **Part 5 of the new test cannot fail on the hazard it names** (`reliability`, Blocker). It
   asserted a property of `Vec::pop` over data it had just edited and never read the emitted
   document. A writer densifying the (stratum × slippage group) axis with zero rows — §5's fifth
   row's exact prohibition — added three rows to the file and left the test green.
2. **An all-unmeasured contamination table is a second, untested spelling of §5's first state**
   (`reliability`, Major). Read literally it takes the mixture path where absence takes the plain
   formula.
3. **A `ContaminationMeasurement` with both counts zero is writable** (`reliability`, Major) — the
   in-memory `UNMEASURED_READ_GROUP` shape, which a projection written from the view rather than
   from the estimate would produce.

## 6. Findings

**1 Blocker, 4 Majors, 4 Minors, 2 nits.** Per-finding text is in
`tmp/review_2026-08-28_parameters-file-a3/`, which is gitignored.

### Blocker

- **B1** part 5 of `each_of_the_five_states_is_a_missing_key_and_not_a_value` cannot fail on the
  zero-row writer. *(reliability, proved by mutation)*

### Major

- **M1** an all-unmeasured contamination table says §5's first state a second way, untested.
  *(reliability)*
- **M2** `ContaminationMeasurement` admits both evidence counts as zero. *(reliability)*
- **M3** §2.1's wholesale demotion has nowhere to write itself for the four quantities that carry
  a different kind of warrant. *(design fidelity — open question 1)*
- **M4** the substitution-rate axis grows with the cohort and §9 does not price it. *(design
  fidelity — open question 2)*

### Minor

- **m1** the empty-list refusal is handed to step C1, which the plan gives parsing only; the
  semantic checks are C2's. *(reliability)*
- **m2** part 4 clears the whole spectrum table when the fixture already contains the real gap —
  "an empty vector writes an empty array" is weaker than "the furnished stratum writes no row
  while its neighbour writes one". *(reliability)*
- **m3** `a_mistyped_key_is_refused_rather_than_absorbed` uses `replace`, inserting the key twice,
  because the fixture emits `ploidy = 2` at two levels. *(reliability)*
- **m4** `Inbreeding::by_sample`'s "at least one is required" is enforced nowhere and owned by
  nobody, unlike its sibling. *(design fidelity)*

### Nits

Two: a trailing `assert_ne!` in parts 1 and 4 that adds nothing the assertion above it did not
prove; and "rides on two counts being zero", where the predicate is *either* count being zero —
`likelihood/mod.rs` tests both one-sided cases.

## 7. Out of scope observations

The design-fidelity pass reports that **the types cover §3 in full** — every subsection has a
home, in the right unit, at the right grain — and that **all five states of §5 are expressible and
distinct**, which is Checkpoint A's criterion. On spec §12 question 2 (whether the prior's moments
are written beside the seed) it notes the shape can take them later, **but in one direction only**:
`deny_unknown_fields` means a file written after the addition is *refused* by an older build rather
than ignored, so "marked informational and ignored on read" needs the field declared now if a
version-tolerant read is wanted.

## 8. Missing tests to add now

Two, both applied: the all-unmeasured table and the evidenceless measurement, recorded as accepted
today so C2's refusal has a failing test to flip.

## 9. What's good

- Parts 1, 2 and 3 each killed their own targeted mutant, including the one that drops the warrant
  — so the answer to "could part 3 pass for an implementation that ignores the warrant" is no.
- **The golden file is the module's strongest test**: it killed four of the five mutants,
  including the one the five-states test missed.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file -- --ignored regenerate`
