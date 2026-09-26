# Code Review: estimation_memory_issue1
**Date:** 2026-09-26
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** issue 1 of `doc/devel/implementation_plans/estimation_memory.md` — `CohortFit` keeps each stratum's substitution counts instead of the whole repeat-tract evidence
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** PR diff, commit `a75855a9` against `1d5ca86b`, branch `evidence-release`.
- **In-scope files:** [census_fit.rs](../../../../src/run/census_fit.rs),
  [ssr_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit.rs) (the new
  `StratumSubstitutionCounts`, `StratumEvidence::substitution_counts`, `substitution_rate`, the new
  test), the implementation report, the `PROJECT_STATUS.md` block.
- **Out of scope:** the rest of `ssr_fit.rs`, which issue 2 changes.
- **Categories dispatched:** three agents, each in its own worktree, each covering a group of
  categories — reliability + errors + refactor_safety; naming + idiomatic + smells + defaults;
  module_structure + float_portability + tooling + extras. `unsafe_concurrency` was not dispatched:
  the diff adds no concurrency. **Deviation from the skill's one-agent-a-category fan-out**,
  recorded here: the diff is about 60 lines in two files, and grouping kept the run to three agents.

## 2. Verdict

Approve-with-changes. The change does what plan §4 asks and moves no number; the findings are test
coverage and wording.

## 3. Execution status

- `cargo fmt --check` → clean. `cargo clippy --all-targets --all-features -- -D warnings` → clean.
- `cargo test --all-targets --all-features --no-fail-fast` → 4,865 passed, 3 failed, all three in
  `examples/ng_generic_loci_dump.rs` and `examples/ng_ssr_loci_dump.rs`, which reference none of the
  changed code and whose failures `PROJECT_STATUS.md` already records.
- Reliability agent, targeted baseline: `cargo test --lib -- substitution_counts census_fit
  cross_platform` → 11 passed. Mutation testing: 8 mutations run, 0 survived that changed
  behaviour, 1 changed no behaviour (it zeroed a count the fixture already had at zero).
- Structure agent: 2 mutations, both caught (one only by the checksum test).
- `cargo doc`, `cargo audit`: not run (no public-API or dependency change).
- Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. `fit_strata` returns one outcome a stratum, in the evidence's order — assumed by Mi1's test.
   Checked: both of its arms map over `strata` in order, and the test passes.

## 5. Top 3 priorities

1. **Mi1** — nothing but the checksum sees the counts handed from the fit to the file.
2. **Mi2** — no end-to-end fixture has a stratum with nothing compared, or a nonzero rate.
3. **Mi3/Mi4** — two doc comments still say the rates are "read off the evidence".

## 6. Findings

### Minor

**Mi1: src/run/census_fit.rs:70 — the handover of counts has no test of its own.**
**Categories:** reliability. Confidence: High. Emptying the vector, doubling `bases_compared`, or
shifting a stratum's repeat count by one each changed the checksum test's MD5, and all 8 of
`census_fit`'s own tests still passed. When the checksum is re-recorded for an intended change, a
handover defect landing in the same commit would be recorded with it. Fix: assert in
`a_cohort_of_censuses_is_fitted_both_halves` one count a stratum, in the fit's order, and at least
one with bases compared.

**Mi2: src/run/census_fit.rs:310 — no end-to-end test of a stratum with nothing compared, or of a
nonzero rate.** **Categories:** reliability. Confidence: High. The fixture's one stratum compares
460 bases with 0 mismatching, so the loop never skips a stratum and never writes a nonzero rate. A
later `unwrap_or(0.0)` — writing a fitted zero where nothing was measured — would pass every test.
Fix: a test that adds two strata's counts to a real fit before assembling the file.

**Mi3: src/run/census_fit.rs:306 — the loop variable is `stratum` but holds counts**, so the code
reads `stratum.stratum.period`; the comment above still says "read off the evidence". **Categories:**
naming, extras. Fix: rename to `counts`, reword.

**Mi4: src/run/census_fit.rs:70 — the field doc says the rates are "read off the evidence"** and
mixes what the field holds with the history of the field it replaced. **Categories:** naming,
extras. Fix: say what it holds, then why it survives the evidence.

**Mi5: implementation report — the pre-existing failures are placed in the wrong status block.**
**Categories:** extras. They are recorded in "Step 5/6 — STR observations through a run", which
also records the `ng_ssr_loci_dump` failure the report presents as unrecorded. The claim that
neither example references the changed code is true.

### Nits

- The `StratumSubstitutionCounts` doc's "two numbers a stratum where the evidence is…" does not
  parse on first reading, and the 0.7 GiB figure has no source.
- The new type's `stratum` field has no doc.
- The new test's last assertion is true by construction: `StratumEvidence::substitution_rate` now
  calls the counts' method.
- The `map` closure's argument is named `stratum` for a `StratumEvidence`.
- The three count fields are declared on both types; `StratumEvidence` could hold the new type.
  Not taken: `StratumEvidence` is built field by field in several places, and that change reaches
  past this step.
- `drop(evidence)` frees nothing earlier than the function's end would; kept for the intent.

## 7. Out of scope observations

- Plan §4 says the smaller-evidence follow-up is "listed in §7"; it is in §8.
- Whether freeing the evidence lowers resident memory depends on the allocator returning pages;
  only the kimura rerun (plan §8) shows it.

## 8. Missing tests to add now

- `a_cohort_of_censuses_is_fitted_both_halves` — extend with Mi1's assertions.
- `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own` — Mi2.

## 9. What's good

- The rate is one expression in one place ([ssr_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit.rs)),
  so the parameters file and the harness that prints from the evidence cannot disagree.
- The new struct literal spells out every field, so a new count added to the evidence is a compile
  error there.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib -- census_fit substitution_counts cross_platform`
