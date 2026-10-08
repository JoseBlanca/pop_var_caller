# Code Review: fit_precision_f1
**Date:** 2026-10-08
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step F1 — own-fit standard errors kept beside blended slippage numbers; a log line counting every stratum
**Status:** Request-changes (all findings since applied — see [fixes_applied_fit_precision_f1_2026-10-08.md](fixes_applied_fit_precision_f1_2026-10-08.md))

---

## 1. Scope

- **What was reviewed:** the diff of commit `5e1d43ec` (a temporary commit, folded into the step's commit) against
  `950909f9`, branch `fit-precision-followup`.
- **In-scope files:** `src/parameter_estimation/joint/stratum_fits.rs`, `src/calling/parameters_file/validate.rs`,
  `src/calling/parameters_file/mod.rs` (fixture), `src/calling/parameters_file/to_toml.rs` (note),
  `src/calling/parameters_file/testdata/every_shape*.toml`, `src/parameter_estimation/joint/ssr_fit.rs`
  (`strata_outcomes_summary`, its call site and test); the spec amendment and plan section for consistency.
- **Out of scope:** the rest of `ssr_fit.rs`.
- **Categories dispatched:** two reviewers in their own worktrees. One ran reliability, errors, refactor_safety and
  float_portability, with mutation testing; the other ran naming, idiomatic, smells, defaults and extras, and
  checked every number the diff states (step 8a). Per-category files: `tmp/review_2026-10-08_f1_own_fit_errors/`.

## 2. Verdict

Request-changes. The change does what the amendment says, but a blended fall-off was not tested at all, and the
reader's rule for the two shares was not pinned in either direction.

## 3. Execution status

- `cargo test --release --lib -- parameters_file stratum_fits ssr_fit::tests::the_strata_outcomes ssr_fit::tests::the_subsets`:
  290 passed, 0 failed.
- Mutation testing: 19 mutations, 6 survived, 0 changed no behaviour (each survivor shown to change the function's
  answer by a probe test).
- fmt, clippy and the full suite: not run at review time; run before the commit.
- Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

None needed an answer from the owner.

## 5. Top 3 priorities

1. B1 — a blended fall-off is untested; the old rule for it survives the suite.
2. M1 — the reader's share rule is unpinned; reverting it would refuse every file a real fit now writes.
3. M3 — four figures or mechanisms in the new prose are wrong.

## 6. Findings

### Blocker

**B1: src/parameter_estimation/joint/stratum_fits.rs:1216 — no test has a blended fall-off.** Putting the old rule
back for the fall-off only (`matches!(…, ShareSource::Stratum)`) passes the whole suite, and on a blended fall-off
the error is dropped again — the defect the amendment exists to fix. *Fix:* a fitted group with a blended fall-off,
asserting its error is kept. **Categories:** reliability.

### Major

**M1: src/calling/parameters_file/validate.rs:1741, mod.rs:1780 — the reader's rule for the two shares is not
pinned.** The fixture's blended shorter share carries no error, and it has no blended fall-off, so reverting either
share's rule survives — and every file a fit now writes with a blended share would then be refused on reading. An
error beside a shorter share taken whole from the curve is also accepted unnoticed. *Fix:* an error beside the
fixture's blended shorter share; accept and refuse cases for each share. **Categories:** reliability, extras.

**M2: four doc comments outside the diff still state the old rule** — `LevelOrigin::own_fit_standard_error`
(mod.rs:1132), `SharesOrigin`'s share errors (mod.rs:1157), `StratumFit::standard_errors` (ssr_fit.rs:442), and the
validator test's doc (validate.rs:1735). **Categories:** smells, reliability (cross-category).

**M3: wrong figures and one wrong mechanism in the new prose.** (a) The plan says the furthest SNP/indel parameters
were "400 or more errors away"; start 2's furthest was 44.27. (b) The weight range 2 × 10⁻⁵ to 1 × 10⁻⁴ was measured
on the 22 one-base strata, not all 31 blends. (c) "`curve_weight` says how far the written number departs from the
own fit" is wrong: the level is blended on the logit scale, and the weight says how much of the number is the
curve's. (d) "usually tiny" rests on one large cohort; on a thin one the curve is most of the number. (e) The file
note said a `this_stratum` or `blend` origin "carries" an error; it carries one only where the own fit gave one.
**Categories:** claims, extras.

### Minor

**Mi1: stratum_fits.rs:124/127/130 and validate.rs:815/859/869 — the rule is a deny-list (`!matches!(…, Curve)`)**,
so a source added later would keep the own-fit error without the compiler asking. *Fix:* one exhaustive match per
enum. **Categories:** refactor_safety, idiomatic.

**Mi2: validate.rs:1765 — fixture rows picked by bare index** where named constants exist for exactly that.
**Categories:** naming.

**Mi3: spec fit_precision.md:607 — the table row still states the old rule in bold**, the amendment appended below
it. **Categories:** extras.

**Mi4: no fitted stratum in the tests takes its level whole from the curve**, so deleting the level filter passes.
**Categories:** reliability.

**Mi5: ssr_fit.rs — deleting the call that writes the new log line passes every test.** The function itself is
well tested (all five of its mutations caught). **Categories:** reliability.

### Nits

`let below` names a clause, not a count; "1 were furnished"; the two below-the-floor phrasings differ; "holding 3
to 3" for one stratum; the summary prints the configured floor rather than each refusal's own (always equal today);
the "Its fall-off" comment in the validator test refers to the wrong row; the bool parameter of `an_own_fit_error`
reads as a double negative.

## 7. Out of scope observations

- No cross-platform checksum fixture writes a slippage row, so the slippage rows' bytes are pinned across platforms
  only through the four-accession oracle. Pre-existing.

## 8. Missing tests to add now

- `own_fit_errors_are_kept_for_every_number_the_own_fit_went_into`: a third group — level from the curve whole,
  shorter share the stratum's, fall-off a blend (B1, Mi4).
- In `every_key_version_two_added_is_refused_where_it_means_nothing`: a shorter share taken whole from the curve
  refused; a blended fall-off accepted (M1).

## 9. What's good

- `strata_outcomes_summary` matches `StratumOutcome` and `StratumRefusal` with no catch-all, so a new outcome is
  flagged where it is counted (ssr_fit.rs).
- The error kept is taken from the own fit computed before the curves are drawn, so it cannot be the blend's by
  accident (stratum_fits.rs).

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --release --lib -- parameters_file stratum_fits the_strata_outcomes`
- `./scripts/dev.sh cargo test --release --lib cross_platform_digests`
