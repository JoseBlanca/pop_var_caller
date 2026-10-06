# Fix Application Report: fit_precision_a7_2026-09-28.md

**Date:** 2026-09-28
**Source review:** `doc/devel/reports/reviews/fit_precision_a7_2026-09-28.md`
**Source state reviewed against:** `8be6683f` (review-only commit on no branch; parent `41e0f9bf`)
**Execution mode:** non-interactive
**Overall status:** Completed (two items carried to the owner at checkpoint A′)

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 5
- Minors: 10
- Nits: 5

### Outcome totals
- Applied: 15
- Applied with adaptation: 4
- Already fixed: 0
- Deferred: 3
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0 (M2's cost and the "never lowering" ruling go to the owner at the checkpoint; the code needs no answer)

### Validation summary
- `cargo fmt --check` → 0, clean (after one `cargo fmt` of `information.rs`)
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (the vendored `noodles-cram` warnings, as before)
- `cargo test --all-targets --all-features --no-fail-fast` → 101: 4,930 passed, 3 failed, 5 ignored; the three are the pre-existing failures (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1). `tmp/fit_precision/suite_a7b.log`.
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 72 passed, 1 ignored (last run before commit)
- `cargo doc --no-deps` → not run
- `cargo audit` → not run
- Performance check → not run as a criterion comparison (see §9)

### Unresolved high-priority findings
- M2 — the update takes more passes than the one it replaces; measured and stated, raised with the owner (checkpoint A′), no code change.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | the step at a bound | Apply | Applied with adaptation | No | `fit.rs` | Pass | No |
| B2 | Blocker | the pass's slope sums untested | Apply | Applied | No | `fit/information.rs` | Pass | No |
| M1 | Major | monotonicity measured on interior cohorts only | Apply | Applied | No | `fit.rs` | Pass | Owner rules at A′ |
| M2 | Major | more passes than the update it replaces | Apply | Applied (doc and report) | No | `fit.rs`, report | Pass | Owner at A′ |
| M3 | Major | the monotonicity test passes wrong updates | Apply | Applied | No | `fit.rs` | Pass | No |
| M4 | Major | the flat point checked only by the ignored test; no carrier | Apply | Applied with adaptation | No | `fit/information.rs` | Pass | No |
| M5 | Major | "wider than their range" is false | Apply | Applied | No | `cli/cross_platform_digests.rs` | Pass | No |
| m1 | Minor | the report does not say how the oracle's calls moved | Apply | Applied | No | report | N/A | No |
| m2 | Minor | "0.4 to 3.1 beside the maximum" | Apply | Applied | No | `fit.rs` | N/A | No |
| m3 | Minor | module doc's guarantee | Apply | Applied | No | `fit.rs` | N/A | No |
| m4 | Minor | spec §3.2 names `fit_beta_shapes` | Defer | Deferred | No | None | N/A | Owner at A′ |
| m5 | Minor | the carrier's rule and slopes as two `Option`s | Apply | Applied | No | `fit.rs`, `fit/information.rs` | Pass | No |
| m6 | Minor | the slope formula twice | Apply | Applied with adaptation | No | `fit/information.rs` | Pass | No |
| m7 | Minor | parameter `weight` | Apply | Applied | No | `fit.rs` | Pass | No |
| m8 | Minor | a non-finite slope becomes NaN shapes | Apply | Applied | No | `fit.rs` | Pass | No |
| m9 | Minor | rest-point wording without the bound exception | Apply | Applied | No | `fit.rs`, report | N/A | No |
| m10 | Minor | no limit on the step's length | Defer | Deferred | No | None | N/A | Yes |
| n1 | Nit | the determinant guard is silent | Apply | Applied with adaptation | No | `fit.rs` | Pass | No |
| n2 | Nit | `RuleSlopes`' fields `pub(super)` | Defer | Deferred | No | None | N/A | No |
| n3 | Nit | the report's test table | Apply | Applied | No | report | N/A | No |
| n4 | Nit | the step's test shares the code's formula | Apply | Applied | No | `fit.rs` | Pass | No |
| n5 | Nit | the seed-search claim has no log | Apply | Applied | No | report | N/A | No |

## 3. Questions asked and answers

None asked during the run. Two items go to the owner at checkpoint A′ with a recommendation: whether
"never lowering the log-likelihood" may stay measured rather than enforced (M1), and the pass count
(M2).

## 4. Per-finding log

### B1 — the step at a bound
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** confirmed by two reviewers independently, with measurements; one correct path.
- **Implementation summary:** `step_beta_shapes` takes the joint step; if a shape's step leaves the
  interval, that shape is held at the bound it crosses and the other takes the conditional step
  `(g_j − h_ji·Δi)/h_jj` given the held shape's actual move; if that too leaves the interval, both
  are clamped. A non-finite slope moves nothing (m8).
- **Review suggestion used verbatim?:** No. **Adaptation:** the design reviewer's version drops the
  held shape's row whenever it sits on its bound and pushes out; the correctness reviewer's takes
  the conditional step given the held shape's move. This follows the second, which also covers a
  step that *crosses* a bound from inside. On all 36 of the correctness reviewer's fits (12 cohorts ×
  3 starts) this code gives that reviewer's passes and log-likelihoods to the printed digits
  (`tmp/fit_precision/a7b_passes.log` against its `passes_fix_full.log`).
- **Verification performed:** the reviewer's bound probe re-run on this code
  (`tmp/fit_precision/a7b_bounds.log`): at 20 samples × 3 reads drawn at (15, 200) the first shape
  rests at 3.44–3.47 with its slope 0.007 errors from zero, 341 to 355 units above the clamped fit;
  at 20 × 8 drawn at (2, 120), 0.46–0.47 and 331 to 372 units; at 4 × 3, 28 to 31 units.
- **Files changed:** `src/parameter_estimation/joint/fit.rs`
- **Tests added or modified:** `the_shapes_rest_where_the_likelihoods_slope_is_zero_at_a_bound_too`
  (new); projected cases in `the_shapes_step_follows_the_slope_by_the_complete_data_curvature`; a
  bound cohort in `plain_passes_never_lower_the_log_likelihood`. Returning to per-shape clamping fails
  three tests (`tmp/fit_precision/a7fix_mutations.log`).
- **Validation:** `cargo test --release --lib parameter_estimation::joint::fit` → 0, 72 passed.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** the oracle and the fixture moved again; both measured (report §2).

### B2 — the pass's slope sums untested
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `the_passs_shape_slopes_are_the_likelihoods` — the pass's summed
  density and carrier slopes equal the scorer's (finite-differenced elsewhere) within 10⁻⁷ relative
  over a cohort of more than two chunks, with and without the duplicated class;
  `the_passs_carrier_slopes_follow_the_coverage_odds` — the same with each sample's coverage odds set.
  Measured agreement at most 1.8 × 10⁻⁹.
- **Review suggestion used verbatim?:** Yes, adapted to the module's helper names.
- **Files changed:** `src/parameter_estimation/joint/fit/information.rs`
- **Validation:** fit module tests → 0.
- **Residual risk:** None

### M1 — monotonicity measured on interior cohorts only
- **Final status:** Applied. The bound fixed (B1); the measurement re-run with nine bound cohorts
  added (`tmp/fit_precision/a7b_mono.log`: 10,020 plain passes, 161 fall, by at most 2.9 × 10⁻¹¹
  units — the rounding of a log-likelihood of about 10⁵ that has stopped moving; the eight interior
  cohorts none of 1,920, as before); a bound cohort (4 samples, b drawn at 150, 150 passes) is in the
  test and fails the clamped step. The correctness reviewer's duplicated 20-sample cohort is in the
  measurement, not the test: the test's bound cohort already fails the defect,, and adding the
  other would lengthen the fit module's tests by a 20-sample fit run three times over. **Whether "never lowering" may stay measured is the
  owner's to rule** (checkpoint A′).

### M2 — more passes than the update it replaces
- **Final status:** Applied (doc and report); the cost itself is for the owner.
- **Implementation summary:** `step_beta_shapes`' doc now says it takes more passes, with the
  measured range; the report states the pass counts on the projected code
  (`tmp/fit_precision/a7b_passes.log`). With the bound fixed, the slowdown holds where the shapes are
  inside their bounds (20 samples × 3 reads: 18–30 accelerated passes against 12–18 before) and
  reverses where one sits on a bound (20 × 8 duplicated, carrier at 50: 28–47 against 32–78).

### M3 — the monotonicity test passes wrong updates
- **Final status:** Applied. Its doc says what it cannot see (a frozen or timid update) and names the
  tests that hold the slopes and the rest point; the one-sample cohort, which could not fail, is
  replaced by the bound cohort.

### M4 — the flat point checked only by the ignored test; no carrier
- **Final status:** Applied with adaptation. The fast guards are B2's two tests and the bound
  rest-point test; no duplicated regime was added to the coverage test (it would redraw every cohort
  and double its 35 minutes); its doc now says the carrier's rest point is not checked there and names
  the tests that hold its slopes.

### M5 — "wider than their range"
- **Final status:** Applied. The note now says the fit reports neither shape as identified, and
  describes the projected step's move (report §2).

### m1 — the oracle's calls — Applied (report §2: records, genotypes, inbreeding coefficients, the
prior's concentrations, the duplicated share, and that every start stops at the pass limit).

### m2, m3, m9 — doc wording — Applied (`step_beta_shapes`' doc: 0.65 to 2.3 from the flat point; the
module doc: one step, measured not guaranteed; the bound exception stated in the step's doc).

### m4 — spec §3.2 — Deferred: the spec is the owner's; offered as an amendment at checkpoint A′.

### m5 — the carrier's two `Option`s — Applied. `PassModel.carrier: Option<CarrierRule>`, a struct
holding the rule and its slopes, the slopes sized from the carrier rule's own node count;
`one_position` and the scorer read both from it.

### m6 — the slope formula twice — Applied with adaptation: not merged into one function; B2's tests
hold the two copies together within 10⁻⁷, the alternative the design reviewer offered. The skip
thresholds differ by design (the pass skips nodes whose share is at most 10⁻¹², for speed); the
measured difference is at most 1.8 × 10⁻⁹ relative.

### m7 — `weight` — Applied: renamed `positions` (the posterior count the Beta covers). The shapes stay
two positional arguments.

### m8 — non-finite slope — Applied: returns the shapes unmoved; asserted in the step's unit test.

### m10 — no step-length limit — Deferred. The step divides by the complete-data curvature, which is
never flatter than the likelihood's own (the missing information is subtracted from it, not added),
so it is shorter than a Newton step on the likelihood; no measurement here showed a jump to a bound
that the likelihood did not ask for. Recorded in case plan step B's stopping rule meets one.

### n1 — the determinant guard — Applied with adaptation: kept (the reliability reviewer shows it is
reachable at a vanishing count of positions) and asserted: a count of 10⁻¹⁷⁰ leaves the shapes
where they were.

### n2 — `RuleSlopes`' fields — Deferred: a read-only accessor changes nothing a test can see; left for
a cleanup.

### n3, n5 — report — Applied: the test table describes the projected cases; the seed-search sentence
is gone (the trace test now reads the winner from its rows).

### n4 — hand-worked value — Applied: at (1, 1) one position's curvature is [[1, −c], [−c, 1]] with
c = π²/6 − 1, so a slope of (1, 0) moves the shapes by (1, c)/(1 − c²); held to 10⁻⁶.

## 5. Deferred findings to carry forward
- m4 — spec §3.2 amendment, the owner's.
- m10 — no step-length limit; revisit if step B's stopping rule meets a jump.
- n2 — `RuleSlopes` accessor.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check

- **Triggered:** no criterion comparison. The fixes add no work to a pass (the projection is in the
  once-a-pass maximisation; the carrier's rule and slopes are regrouped, not recomputed). The step's
  own per-pass cost was measured before review (+4 to 7%, `tmp/fit_precision/a7_timing.log`) and is
  unchanged by the fixes; its pass counts are measured (M2).

## 10. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`
- `scripts/dev.sh cargo test --release --lib plain_steps_never_lower -- --ignored --nocapture` (probe appended, then removed)
- `scripts/dev.sh cargo test --release --lib review_probe::probe_passes_and_chain -- --ignored --nocapture` (reviewer's probe, then removed)
- `scripts/dev.sh cargo test --release --lib review_probe::probe_bounds -- --ignored --nocapture`
- `scripts/dev.sh cargo test --release --lib cli::cross_platform_digests`
- `DEV_EXTRA_MOUNT=… scripts/dev.sh sh tmp/fit_precision/run_oracle.sh …/oracle_a7b`
- `scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture`
- `scripts/dev.sh cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`; `cargo test --all-targets --all-features --no-fail-fast`

## 11. Command results
- fit module tests → 0, 72 passed, 1 ignored
- monotonicity probe → 0; 161 of 10,020 passes fall, by at most 2.9 × 10⁻¹¹
- pass-count probe → 0; 36 fits, the reviewer's fixed code's figures
- bound probe → 0; see B1
- cross-platform digests → 101 before re-recording (the expected move); passed after it, in the full suite
- oracle → 1 (the stale baseline's two fit lines, as since step A1); two checksums moved, report §2
- coverage run → see the implementation report §3
- fmt → 0; clippy → 0; suite → 101, 4,930 passed, 3 failed (pre-existing), 5 ignored

## 12. Notes
- The design and reliability reviewers' worktrees, and the correctness reviewer's, were removed after
  their findings and evidence were copied to `tmp/review_2026-09-28_fit_precision_a7/`.
