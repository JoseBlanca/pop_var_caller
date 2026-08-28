# Fix Application Report: ng_parameters_file_b1_2026-08-28.md

**Date:** 2026-08-28
**Source review:** [ng_parameters_file_b1_2026-08-28.md](ng_parameters_file_b1_2026-08-28.md)
**Source state reviewed against:** the uncommitted B1 diff over `90885e48`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 9
- Minors: 24
- Nits: 7 classes

### Outcome totals
- Applied: 31
- Applied with adaptation: 3
- Deferred: 2
- Disputed: 2
- Awaiting the owner: 2 (M9 and Mi2, both raised at Checkpoint B)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::calling::parameters_file` → 0, **45 passed, 0 failed, 1 ignored**
- `cargo test --lib ng::parameter_estimation::joint::stratum_fits` → 0, **30 passed, 0 failed**
- `cargo test --lib` → 0, **4,968 passed, 0 failed, 12 ignored**
- `cargo test --all-targets --all-features` → not run as a gate; pre-existing panic in
  `benches/psp_writer_perf.rs:386`, verified on clean `main`
- Performance check → skipped; nothing in the diff is reachable from `benches/`

### The six mutations the review proved survived were re-run against the fixed tree, and **all six
now fail a test**

Each was applied to the fixed tree, the module suite run, and the tree restored from a pristine
copy and verified byte-identical by `diff`:

| mutation | test that now fails |
|---|---|
| the library-count guard deleted | `a_read_group_table_with_another_runs_library_count_is_refused` |
| `intercept`/`slope` swapped in **both** curve conversions | `every_field_of_both_curves_reaches_the_file_under_its_own_name` |
| `shares: row.shares[0]` instead of `[index]` | `a_number_off_a_curve_carries_the_curve_and_a_stratums_own_carries_neither` |
| `was_declared_by_the_run: true` | `a_run_that_declared_no_batching_writes_the_flag_false` |
| the stated concentration given an observation count | `the_stated_concentration_says_whether_the_run_fitted_it` |
| `FittedShape::Flat => ShareShape::Sloping` | `every_pre_pass_word_maps_to_its_own_word_in_the_file` |

### Unresolved, and both are the owner's
- **M9** — a slippage number carries an origin and no warrant, so spec §2.1's wholesale demotion
  has nowhere to write itself. Convergent with the item Checkpoint A already raised.
- **Mi2** — recorded as decided by the coder (below) and flagged at Checkpoint B in case the owner
  reads it the other way.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| B1 | Blocker | the refusal test cannot fail | Apply | Applied |
| B2 | Blocker | the calibration count's join is unchecked and untested | Apply | Applied |
| M1 | Major | seven panics abort a finished run | Ask → coder's | Applied with adaptation |
| M2 | Major | eight of the two curves' sixteen fields untested | Apply | Applied |
| M3 | Major | `was_declared_by_the_run` only tested `true` | Apply | Applied |
| M4 | Major | eight enum-mirror arms entered by no test | Apply | Applied |
| M5 | Major | a wider rate set accepted in silence | Apply | Applied |
| M6 | Major | three `expect()`s, one naming nothing | Apply | Applied with adaptation |
| M7 | Major | the three new iterators untested; a rung nobody observes | Apply | Applied |
| M8 | Major | five documented guards untested | Apply | Applied |
| M9 | Major | a slippage number has an origin and no warrant | Ask | **Awaiting the owner** |
| Mi1 | Minor | a substitution rate keyed past the axis | Apply | Applied |
| Mi2 | Minor | the median concentration's warrant | Ask → coder's | Applied with adaptation |
| Mi3 | Minor | an all-unmeasured contamination table is writable | Defer | Deferred to C2 |
| Mi4 | Minor | an unmeasured row drops a fraction | Apply | Applied |
| Mi5 | Minor | `ploidy` from a fixture where four sources agree | Apply | Applied |
| Mi6 | Minor | the shares index is never exercised | Apply | Applied |
| Mi7 | Minor | rebuilding a `ReadGroupId` from an index | Apply | Applied |
| Mi8 | Minor | `warranted`'s function-pointer parameter | Apply | Applied |
| Mi9 | Minor | `&BTreeMap` in a `pub` signature | Apply | Applied |
| Mi10 | Minor | helpers taking the whole `&RunParameters` | Apply | Applied |
| Mi11 | Minor | the `SeedRung` doc opens with a bug's history | Apply | Applied |
| Mi12 | Minor | rationale in an unrendered `//!` | Apply | Applied |
| Mi13 | Minor | the iterator re-implements `at`'s tail | Apply | Applied |
| Mi14 | Minor | the six-argument call spelled out five times | Apply | Applied |
| Mi15 | Minor | `of_each_…` for two filtering iterators | Apply | Applied |
| Mi16 | Minor | `parameters` as an argument name | Apply | Applied |
| Mi17–Mi20 | Minor | `\|of\|`, `held`, `views`, `batching`, `what` | Apply | Applied |
| Mi21 | Minor | `warranted` a bare participle | Apply | Applied |
| Mi22–Mi24 | Minor | four wrong numbers in changed prose | Apply | Applied |
| Nits | Nit | seven classes | Apply (5) / Dispute (2) | Applied / Disputed |
| — | — | `reference_digest` newtype | Defer | Deferred to D1 |

## 3. Questions asked and answers

None asked of the owner mid-run. **Two are carried to Checkpoint B**: M9, and Mi2 in case the
coder's ruling is not the owner's.

## 4. Per-finding log — the ones that changed a decision

### B1 — the refusal test cannot fail
The two guards' messages now share no opening: the library check opens *"this run's read-group
table covers … libraries"* and the sample check keeps *"the read-group table names … samples"*. The
test split in two, each with a fixture that trips only its own guard — four lanes over the run's
own two plants, and three lanes over three plants. **Verified by re-running the review's own
mutation**: with the library guard removed, `a_read_group_table_with_another_runs_library_count_is_refused`
fails where the whole suite previously stayed green.

### B2 — the calibration count's join
`calibration_rows` now asserts that the calibration's warrant is the rate's, allowing the one
legitimate disagreement: assembly substitutes `Defaulted` when `from_fitted_rate` refuses a rate,
so a defaulted calibration may sit beside a fitted rate and nothing else may. Two tests:
`a_rate_set_from_another_fit_is_refused` and `a_rate_set_missing_one_of_the_runs_read_groups_is_refused`.

### M1 — panic or `Result`: **the panics stay, and the doc now argues this function's case**
*Applied with adaptation.* The reviewer's asymmetry is real and is now written into the `# Panics`
section in as many words: this runs after the last locus, so a panic discards a cohort's calling
work where `assemble`'s equivalent checks discard a startup. Kept as panics because **every one of
these conditions is a wiring bug in whoever assembles the six arguments, not a state any input data
can reach**, and the alternative to refusing is writing a false provenance into the run's own
record. **What the reviewer is right about is that the decision needs a call site, and there is
none** — `of_run` has no caller yet. The doc says so and names step F1, which decides the order of
the two writes and therefore whether a failed projection can lose a VCF.

### M6 — the three `expect()`s
*Applied with adaptation.* Two of the three are **removed rather than commented**: the length-spectrum
iterators now hand out `(key, &[f64], f64)` instead of a `LengthSpectrum` whose `fitted_weights` is
an `Option` the producer can never leave empty, so the consumer has nothing to unwrap. The third is
now a `panic!` naming the stratum and the slippage group, and it lives in one shared helper —
`fitted_at` — that `StratumFits::at` and the new iterator both use, which also closes Mi13.

### Mi2 — the median concentration's warrant: **`FittedHere`, and the reasoning is now in the field's doc**
*Applied with adaptation.* The reviewer is right that the spec settles neither, and right that §2's
`Borrowed` describes an average over neighbours. The ruling taken: **this warrant is about the
number, not about the locus that ends up using it.** The number in the file did come out of this
run's own data — it is the median of the strata this run fitted — and which rung a *tract* landed on
is `LengthSpectrumRung`, carried per locus. Marking it `borrowed` would put a second answer to that
per-locus question in a run-level field. The field's doc now says this. **Raised at Checkpoint B**
because it is one line for the owner if they read it the other way.

### Mi3 — an all-unmeasured contamination table
*Deferred to C2.* The reviewer's fix collapses it to absence at the writer. Not taken: the writer's
job is to record what the run used, and in that state the run took the *mixture* path with every
fraction zero, where absence asks the reader for the plain formula. Collapsing would make the file
say something the run did not do. The state is not reachable from the estimator, which refuses
below 100 markers — the reviewer could construct it only by hand. C2's `validate` already owns
refusing it.

### Nits disputed
Two, both about the `as` casts. The idiomatic agent verified they cannot truncate; three of the
five are gone anyway because two helpers now zip `ReadGroups::iter()`, and the two that remain use
`u32::try_from(..).expect(..)`, the crate's own idiom. The `cells as u64` / `strata as u64` widenings
stay `as` because they sit inside `From` impls under `fallible_impl_from = "warn"`.

## 5. Deferred findings to carry forward
- **Mi3** → step C2's `validate`.
- **`reference_digest` as a newtype** → step D1, which is where the digest's producer lands. The
  expected form is now named in `of_run`'s doc.

## 6. Disputed findings to return to reviewer
- The two cast nits above.

## 7–8. Failed validation / blocked by context mismatch
None.

## 9. Performance check
Skipped — nothing in the diff is reachable from any harness in `benches/`.

## 10–11. Commands run

| command | exit | result |
|---|---|---|
| `./scripts/dev.sh cargo fmt --check` | 0 | clean |
| `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `./scripts/dev.sh cargo test --lib ng::calling::parameters_file` | 0 | 45 passed, 1 ignored |
| `./scripts/dev.sh cargo test --lib ng::parameter_estimation::joint::stratum_fits` | 0 | 30 passed |
| `./scripts/dev.sh cargo test --lib` | 0 | 4,968 passed, 12 ignored |

## 12. Notes

- **Test counts:** the projection module went 16 → **31** tests and `stratum_fits.rs` 27 → **30**;
  the whole library suite 4,950 → **4,968**.
- **The fixture grew a fourth stratum** — period 3 at 9 repeats, derived, whose slippage group 1
  carries a blended shares provenance where its group 0 carries none. That one change reaches the
  `ShareSource::Blend` arm nothing entered, and kills the `row.shares[0]` mutation, which was
  byte-identical on the old fixture.
- **The second substitution-rate key is now tetraploid**, so a row that wrote the run's own ploidy
  (2) or its period in place of the key's is visible.
