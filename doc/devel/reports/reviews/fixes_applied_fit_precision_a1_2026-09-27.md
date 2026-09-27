# Fix Application Report: fit_precision_a1_2026-09-27.md

**Date:** 2026-09-27
**Source review:** `doc/devel/reports/reviews/fit_precision_a1_2026-09-27.md`
**Source state reviewed against:** `190c549e` (review-only commit), branch `fit-precision`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 2
- Minors: 14
- Nits: 9

### Outcome totals
- Applied: 11
- Applied with adaptation: 2
- Already fixed: 0
- Deferred: 4
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0 (two owner decisions are carried to checkpoint A, §3)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features --no-fail-fast` → counts in the commit message, from
  `tmp/fit_precision/suite_a1c.log` on the committed tree; the three pre-existing failures only
  (`examples/ng_generic_loci_dump.rs` ×2, `examples/ng_ssr_loci_dump.rs` ×1); both
  `cli::cross_platform_digests` tests pass unchanged
- `cargo test --release --lib parameter_estimation::joint::fit::information` → 0, 10 passed
- `cargo doc --no-deps`, `cargo audit` → not run (no public API or dependency change)
- Performance check → not applicable (nothing calls the scorer outside its tests)

### Unresolved high-priority findings
- None. Both Majors applied.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | held-out reads and coverage odds untested | Apply | Applied | No | `information.rs` | Pass | No |
| M2 | Major | duplicated class off untested | Apply | Applied | No | `information.rs` | Pass | No |
| Mi1 | Minor | "10⁻⁷ or better" | Apply | Applied | No | `PROJECT_STATUS.md` | N/A | No |
| Mi2 | Minor | disagreement not relative; tolerance loose | Apply | Applied | No | `information.rs`, report | Pass | No |
| Mi3 | Minor | digamma quoted at 12 nodes only | Apply | Applied | No | `information.rs`, report, PROJECT_STATUS | Pass | No |
| Mi4 | Minor | later read groups: no recommendation, no test | Apply | Applied with adaptation | No | `information.rs`, report | Pass | Yes — owner |
| Mi5 | Minor | `cohort` index doc wrong without the class | Apply | Applied | No | `information.rs` | Pass | No |
| Mi6 | Minor | index constants rather than named enums | Defer | Deferred | No | None | N/A | A2 |
| Mi7 | Minor | scorer's inputs unchecked against the pass's | Defer | Deferred | No | None | N/A | A2 |
| Mi8 | Minor | model re-derived inline | Apply | Applied | No | `fit.rs`, `information.rs` | Pass | No |
| Mi9 | Minor | `score_position` 230 lines | Defer | Deferred | No | None | N/A | A2 |
| Mi10 | Minor | §5 evidence pointed at a commit message | Apply | Applied | No | report | N/A | No |
| Mi11 | Minor | status line read as all built | Apply | Applied | No | `PROJECT_STATUS.md` | N/A | No |
| Mi12 | Minor | undefined terms | Apply | Applied | No | report, `information.rs` | N/A | No |
| Mi13 | Minor | no-drop rule untestable here | Apply | Applied with adaptation | No | `information.rs`, report | Pass | No |
| Mi14 | Minor | test docs' slope sizes wrong | Apply | Applied | No | `information.rs` | Pass | No |
| Nits | Nit | nine items | Apply | Applied (8), Deferred (1) | No | `information.rs`, `fit.rs`, report | Pass | module name, at A2 |

## 3. Questions asked and answers

None asked in this run. Two owner decisions are carried to checkpoint A, each with a recommendation
in the A1 implementation report: what the parameters file writes for a sample's later read groups
(Mi4, review open question 1), and whether to change the maximisation's Beta-shape update (review
open question 2).

## 4. Per-finding log

### M1 — held-out reads and coverage odds never reach the scorer
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation summary:** the reviewer's test, adapted to the module's helpers, as
  `the_slopes_hold_with_held_out_reads_and_coverage_odds`: one reference read of every fourth
  sample-position becomes a held-out read (583 sample-positions), coverage odds 0.25, 3, 1, 8; the
  carrier shapes, the duplicated share and all eight rates against a central difference.
- **Review suggestion used verbatim?:** No — crate paths instead of `super::super`, the relative
  measure of Mi2, tolerance 10⁻⁶.
- **Verification performed:** passes at 1.38 × 10⁻⁷; the reviewer ran it against both mutants
  (fails at 292 and 0.269).
- **Files changed:** `src/parameter_estimation/joint/fit/information.rs`
- **Validation:** module tests 10 passed; full suite §1.
- **Residual risk:** None.

### M2 — no test with the duplicated class off
- **Severity:** Major · **Initial decision:** Apply · **Final status:** Applied
- **Implementation summary:** `the_slopes_hold_with_the_duplicated_class_off`: the five cohort
  coordinates, eight rates and four excesses of the layout without the class, against a central
  difference, and the class's three slots asserted exactly zero. It kills M9 (the assertion on the
  slots).
- **Verification performed:** passes, worst 1.27 × 10⁻⁷ (sample 2's excess).
- **Files changed:** `information.rs` · **Residual risk:** None.

### Mi1 — "10⁻⁷ or better"
- **Final status:** Applied. Current focus now reads "to within 2 parts in 10 million (1.7 × 10⁻⁷ at
  worst)".

### Mi2 — the disagreement measure and the tolerance
- **Final status:** Applied. `relative_disagreement` divides by `max(|numeric|, 10⁻³)`; `TOLERANCE`
  is 10⁻⁶ for every slope test; the report's table row for the Beta shapes reads 1.70 × 10⁻⁷
  (carrier `a`). Kills M39 (a uniform 5 × 10⁻⁶ error).

### Mi3 — digamma figure at the shipped node count
- **Final status:** Applied. The test asserts the digamma form's miss at both 12 (> 3%) and 16
  (> 2%) nodes; measured 4.71 × 10⁻² and 2.88 × 10⁻², this module's 9.00 × 10⁻⁹ and 9.56 × 10⁻⁹.
  The module doc, the report and PROJECT_STATUS quote 2.9% at 16 beside 4.7% at 12.

### Mi4 — later read groups
- **Final status:** Applied with adaptation.
- **Implementation summary:** the report now recommends writing a later read group's rate with its
  first group's value and error, marked as shared, rather than `defaulted`, and says the 482 on
  kimura is a difference of two counts. New test `a_second_read_groups_rates_have_no_slope`: sample
  0 gains a second read group at other rates; moving either rate leaves the log-likelihood
  bit-for-bit unchanged, and the sample's clean-rate slope still matches its first group's.
- **Adaptation:** the spec amendment the reviewer proposed is not made here — the plan-driven skill
  edits no design text, and the spec note depends on the owner's ruling. Carried to checkpoint A.

### Mi5 — `cohort` index doc
- **Final status:** Applied. The doc says the orders match only with the class fitted.

### Mi6 — named parameter slots
- **Final status:** Deferred to A2, where the information blocks are indexed and A3/A5 first print
  per-parameter names; changing the indices now would be reworked there.

### Mi7 — the scorer's inputs
- **Final status:** Deferred to A2, which adds the one production call site; the bundle is designed
  there against that caller.

### Mi8 — the model re-derived inline
- **Final status:** Applied. `candidate_read_probability` (now also what `ReadLogs::of` calls — same
  expression, same bits), `candidate_read_probability_slope`, `reference_read_probability_slope`,
  `genotype_frequencies_slope_in_frequency` and `genotype_frequencies_slope_in_excess` sit in
  `fit.rs` beside the functions they differentiate; `the_models_slopes_are_the_derivatives_of_its_functions`
  checks each against a finite difference over a grid.

### Mi9 — the function's length
- **Final status:** Deferred to A2, together with precomputing each read's slope outside the node
  loop (review §7): both change the function's shape, and A2 is where it gains its caller.

### Mi10 — the report's evidence
- **Final status:** Applied. §5 quotes the suite's counts, names the three pre-existing failures and
  the two checksum tests, and cites the logs.

### Mi11 — the status line
- **Final status:** Applied. "`fixes-applied` through step A1 of the plan's 17 … The programme: …".

### Mi12 — undefined terms
- **Final status:** Applied. The report defines the digamma function, the Gauss–Jacobi rule, SQUAREM
  and the pass's skipping of small posteriors, and replaces "the E2 ruling" with the owner's ruling
  in words; the module doc defines `ψ` where it first appears and no longer says "attribution".

### Mi13 — the no-drop rule
- **Final status:** Applied with adaptation: the limitation is recorded on the shares test's doc and
  in the report (§2 item 5); the per-position test near a bound was not built — the reviewer's own
  proposal was unrun, and at a share of 10⁻⁹ a step that is both precise and unbiased is not
  established.

### Mi14 — the test docs' sizes
- **Final status:** Applied. The fixture doc gives the measured range, 0.35 to about 36,000
  log-likelihood units a unit of the parameter over every test, except the carrier shapes' −0.20 and
  0.063.

### Nits
- Applied: "four shares"; the range test asserts more than 1,200 of 2,400 sample-positions ranged
  (1,729 measured); the report's `fit.rs` length (4,880 before, 4,944 after); one doc comment per
  `Scratch` field; `node_slope` / `ln_weight_slope`; crate paths in the tests; `2 * samples` named;
  `per_share` renamed `ln_posterior_per_share`.
- Deferred: the module's name, until A2 decides what else lives in it.

## 5. Deferred findings to carry forward
- Mi6 — named enums for the parameter slots (A2).
- Mi7 — one bundle of the scorer's per-pass inputs (A2).
- Mi9 — split `score_position`, and compute each read's slope once per (class, candidate, sample,
  genotype) (A2).
- Nit — the module name (A2).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — the scorer has no caller outside its tests until A2.
- **Outcome:** skipped.

## 10. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit::information -- --nocapture`
- `scripts/dev.sh bash -c 'cargo fmt && cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- module tests → 0, 10 passed (`tmp/fit_precision/module_a1c.log`)
- fmt, clippy, full suite → see §1 and the commit message (`tmp/fit_precision/suite_a1c.log`)

## 12. Notes
- For A2/A3, from the correctness review: at a homozygote excess of exactly 1 a sample's per-position
  slope can reach order 10⁶ (analytic −2.86 × 10⁶ against a finite difference of −6.44 × 10⁵, which
  the reviewer judges unable to resolve it); such a parameter's information would be dominated by a
  few positions.
