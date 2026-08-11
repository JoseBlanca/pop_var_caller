# Fix Application Report: ng_parameter_prepass_generic_a1a2a3_2026-08-06.md

**Date:** 2026-08-06
**Source review:** `doc/devel/reports/reviews/ng_parameter_prepass_generic_a1a2a3_2026-08-06.md`
**Source state reviewed against:** `5b3a646` (ng-parameter-estimation)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 4
- Minors: 16
- Nits: 10

### Outcome totals
- Applied: 26
- Applied with adaptation: 2
- Already fixed: 0
- Deferred: 4
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib ng::types::` → 0, **28 passed** (was 26)
- `cargo test --lib ng::parameter_estimation` → 0, **6 passed** (was 3)
- `cargo test --all-targets --all-features` → 101, **2,906 passed, 1 failed, 5 ignored**. The one failure is `ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes`, pre-existing at `HEAD` and out of scope. 2,906 = 2,901 + the five tests added here.
- `cargo doc --no-deps --lib` → 101, **12 unresolved links, all pre-existing**; the three this commit introduced are gone and no file under `parameter_estimation/` contributes one.
- Performance check → not applicable: no `Apply` touched code reachable from `benches/`.

### Unresolved high-priority findings
None. The four deferrals are all Minor or Nit and each reaches outside the step.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | Untested rejection direction on two rates | Apply | Applied | `types.rs` | Mutation-verified |
| M2 | Major | `Ploidy` doc contrasts with newtypes that are checked | Apply | Applied | `types.rs` | Read back |
| M3 | Major | `is_finite` dead, and its stated reason false | Apply | Applied | `types.rs`, impl report | 28 tests pass without it |
| M4 | Major | `histogram.rs` "nothing is lost" contradicts the design | Apply | Applied | `generic/histogram.rs` | Checked against arch §2.2 |
| Mi1 | Minor | `WindowIndex` + ladder one level above their consumer | Apply | Applied | module move | 6 tests pass |
| Mi2 | Minor | Banner states absent consumers; miscounts what step 4 fits | Apply | Applied | `types.rs` | Read back |
| Mi3 | Minor | Saturating cast: two silent ladder failures | Apply | Applied with adaptation | `generic/mod.rs` | Mutation-verified |
| Mi4 | Minor | Ladder-edge constants restate their value | Apply | Applied | `generic/mod.rs` | Read back |
| Mi5 | Minor | Five copies of one `[0, 1]` predicate | Apply | Applied | `types.rs` | 28 tests pass |
| Mi6 | Minor | `PartialEq` over `f64`: `NaN` error ≠ itself | Apply | Applied | `types.rs` | Read back |
| Mi7 | Minor | `DomainError::ErrorRate` has three producers | Defer | Deferred | none | — |
| Mi8 | Minor | `Ploidy` never offered a large copy number | Apply | Applied | `types.rs` | Mutation-verified |
| Mi9 | Minor | Ladder pinned by literals, not by its constants | Apply | Applied | `generic/mod.rs` | Mutation-verified |
| Mi10 | Minor | Impl report's test count and milestone wrong | Apply | Applied | impl report | Re-ran the suite |
| Mi11 | Minor | "Phred appears here and nowhere else" false at crate scope | Apply | Applied | `generic/mod.rs` | Read back |
| Mi12 | Minor | "coarsest"/"finest" for the ladder's endpoints | Apply | Applied | `generic/mod.rs` | Read back |
| Mi13 | Minor | "Finer than a caller can feel" asserted where the spec marks it soft | Apply | Applied | `generic/mod.rs` | Checked against spec §3 |
| Mi14 | Minor | `INBREEDING_WINDOW_BP` contrast has no referent | Apply | Applied | `generic/mod.rs` | Read back |
| Mi15 | Minor | Pure-reference skip misattributed to production's caller | Apply | Applied | `parameter_estimation/mod.rs` | Checked against spec §2.1 |
| Mi16 | Minor | Two files never define their filename's term | Apply | Applied | `mixture_weights.rs`, `runs.rs` | Read back |
| N1 | Nit | Three unresolved intra-doc links | Apply | Applied | `parameter_estimation/mod.rs` | `cargo doc` |
| N2 | Nit | No `#[must_use]` on `error_rate_ladder()` | Apply | Applied | `generic/mod.rs` | clippy |
| N3 | Nit | `WindowIndex`/`ContigId` comparison supports the opposite | Apply | Applied | `generic/mod.rs` | Read back |
| N4 | Nit | "unchecked" against "unconstrained" | Apply | Applied | `types.rs` | Read back |
| N5 | Nit | `top_rung` holds an index | Apply | Applied with adaptation | `generic/mod.rs` | Variable removed entirely |
| N6 | Nit | `histogram.rs` worked example does not match the design table | Apply | Applied | `generic/histogram.rs` | Example removed |
| N7 | Nit | `fitting/mod.rs` does not name its milestone | Apply | Applied | `fitting/mod.rs` | Read back |
| N8 | Nit | `.expect()` renders `Debug`, not the `Display` message | Apply | Applied | `generic/mod.rs` | Mutation-verified |
| N9 | Nit | `WindowIndex::get()` untested | Apply | Applied | `generic/mod.rs` | New test |
| N10 | Nit | `ng/mod.rs` re-export list has no criterion | Defer | Deferred | none | — |
| Q1 | — | Open question: split `DomainError::ErrorRate`? | Ask | Deferred | none | — |
| Q2 | — | Open question: `DomainError` shape at six variants | Ask | Deferred | none | — |

## 3. Questions asked and answers

None asked of the user. Q1 and Q2 are recorded as open items in `PROJECT_STATUS.md` rather than blocking the step, because both are pre-existing ng-wide conventions whose change would reach outside step 4.

## 4. Per-finding log

Grouped, because the fixes cluster into six coherent patches rather than 30 independent ones. Each group was validated before the next began.

### Group 1 — M1, M5, M8, N4, Mi2, Mi6, M2, M3: `src/ng/types.rs`

- **Implementation.** Added a private `checked_probability(x, reject)` taking the `DomainError` tuple-variant constructor as a `fn(f64) -> DomainError` value, and routed all three rate constructors through it. The `is_finite` clause is gone — it rejected nothing the range check did not already reject. Rewrote the banner (three fitted, one handed; the consumers named as anticipated rather than present). Fixed `Ploidy`'s contrast to "the unconstrained newtypes elsewhere in this file", naming them, and adopted the file's existing word. Documented on `DomainError` that `PartialEq` is IEEE equality so a `NaN` error is not equal to itself, and on `DomainError::ErrorRate` that three constructors raise it deliberately. Added an explicit no-upper-bound paragraph to `Ploidy` — polyploids are in scope.
- **Tests.** `each_constrained_rate_rejects_out_of_range_with_its_own_variant` renamed to `..._in_both_directions` and given the two missing directions. Two proptests added: the rates accept exactly the probabilities and round-trip bit for bit, and ploidy accepts every non-zero `u8`.
- **Verification.** Mutating `checked_probability` to `(0.0..=2.0)` fails the round-trip proptest at `minimal failing input: x = 1.118791511682117`. Before the fix the equivalent mutation on `InbreedingF` alone left the suite green.
- **Adaptation:** the review offered "keep the guard and correct the comment" as the smaller change. Dropped it instead, because the shared helper makes the three constructors one predicate — and a single predicate is what stops a sixth copy drifting to `0.0..1.0` and rejecting a genotype frequency of exactly one.

### Group 2 — Mi1: the module move

`WindowIndex`, `INBREEDING_WINDOW_BP`, the three Phred constants and `error_rate_ladder()` moved from `parameter_estimation/mod.rs` into `generic/mod.rs`. No `use` line anywhere else in the crate changed, because nothing imports them yet. `parameter_estimation/mod.rs` now holds the step's surface and says why it holds no vocabulary: the STR sub-unit joins at that level and has no use for a window size or a ladder of per-base error rates.

### Group 3 — Mi3, Mi9, N2, N5, N8: the ladder builder

- **Implementation.** `ERROR_RATE_LADDER_RUNGS: usize = 161` is now stated, and `error_rate_ladder()` builds that many rungs from `MIN` and `STEP`. The saturating `as u32` cast is gone entirely — there is no rung count to compute. A `const _: () = assert!(…)` pins the three constants to "upward from a non-negative Phred in positive steps", and `.expect()` became an `unwrap_or_else` naming the offending rung and its Phred.
- **Verification, all three failure modes:**
  - `MAX = 50.1` → two tests fail, reporting `the ladder's span must be a whole number of steps, got 160.39999389648438` and `last rung 0.00001 vs 0.00000977237564304849`. Before the fix this was silent.
  - `MIN`/`MAX` swapped → **`error[E0080]: evaluation panicked: the error-rate ladder runs upward from a non-negative Phred in positive steps`**. A build failure, where before it was a one-rung ladder that ran.
  - The ratio test now asserts `windows(2).count() == RUNGS - 1` first, so it cannot pass by iterating nothing.
- **Adaptation:** the review's first choice was a fixed-size array `[ErrorRate; 161]`. Kept `Vec`, which was the review's own stated fallback, because the plan's out-of-scope list anticipates the reference-bias term *lengthening* the ladder later — a length baked into the return type would make that a signature change. The stated constant plus the const assertion closes both silent failures without it.

### Group 4 — Mi4, Mi11–Mi14, N3: the constants' documentation

Both ladder edges now carry what `INBREEDING_WINDOW_BP` already carried: what the rate means in plain terms, "fixed, not a knob", the DRAGstr provenance, and the statement that the remedy for a read group outside the range is the endpoint-argmax flag rather than a wider ladder. "Coarsest/finest" became "noisiest/cleanest", and the test's `coarser`/`finer` bindings became `higher_rate`/`lower_rate`. The step's doc now says the spec *argues* the spacing is below what a caller can feel and *marks the argument soft*. "Phred appears here and nowhere else" became "in step 4 only here", with `BaseQual` and `MapQual` named. `INBREEDING_WINDOW_BP`'s dangling contrast was replaced with the shortest run worth resolving. `WindowIndex`'s `ContigId` comparison now compares on the one axis meant.

### Group 5 — M4, Mi15, Mi16, N1, N6, N7: the placeholder documentation

`histogram.rs` rewritten: the key is three things, not two, and the paragraph that follows says why the library attribution is the part that is easy to drop and cannot be. `parameter_estimation/mod.rs` now says production writes the pure-reference columns and it is the *heterozygosity accumulator* that never looks at them. `mixture_weights.rs` defines "mixture weights" where the filename's term first does work; `runs.rs` defines "run of homozygosity". `fitting/mod.rs` names its milestone. The three intra-doc links became plain code spans. The worked example whose numbers disagreed with the design table was dropped rather than restated.

### Group 6 — Mi10: the implementation report

Corrected the full-suite row to `2,901 passed, 1 failed (pre-existing), 5 ignored`, and "the first thing that reads a number is Milestone D" to name Milestone B's depth means. Corrected assumption 3, which asserted the false `is_finite` rationale, to record what actually happened: the guard was dropped and one predicate now lives in one place.

### Deferred

- **Mi7 / Q1 — splitting `DomainError::ErrorRate`.** Real, and the fix touches `alignment/emission.rs` and `alignment/ssr_marginal_sequence.rs`. That is blast radius beyond this step; it goes on the plan's open list.
- **Q2 — `DomainError`'s shape at six variants.** An ng-wide convention. A decision, not a fix.
- **N10 — the `ng/mod.rs` re-export list.** Pre-existing and partial before this commit; picking a criterion is a separate change.
- **The 1-based window formula.** `naming` noted as a cross-category observation that `start / INBREEDING_WINDOW_BP` on a 1-based `Position` puts 99,999 bases in window 0 and 100,000 in every other. No code computes it yet — it lands in Milestone C3. Recorded in `WindowIndex`'s own doc so the implementer meets it there, and on the plan's open list.

## 5. Deferred findings to carry forward
Mi7 (Q1), Q2, N10, and the 1-based window formula — see above.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no `Apply` touched perf-sensitive code. Nothing under `src/ng/parameter_estimation/` or the four newtypes in `types.rs` is reachable from any harness in `benches/`; the module has no consumers yet.

## 10. Commands run
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::types::`
- `./scripts/dev.sh cargo test --lib ng::parameter_estimation`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh cargo doc --no-deps --lib`
- `./scripts/dev.sh cargo build --lib` (against three deliberate mutations)

## 11. Command results
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, no output
- `cargo test --lib ng::types::` → 0, 28 passed
- `cargo test --lib ng::parameter_estimation` → 0, 6 passed
- `cargo test --all-targets --all-features` → 101, 2,906 passed / 1 failed (pre-existing) / 5 ignored
- `cargo doc --no-deps --lib` → 101, 12 unresolved links, all pre-existing, none in this step's files
- `cargo build --lib` with `MIN`/`MAX` swapped → 101, `error[E0080]` as intended

## 12. Notes

- **The highest-value finding came from mutation, not from reading.** Four of six agents deleted the `is_finite` clause and re-ran rather than reasoning about `RangeInclusive::contains`; the one thing none of them could have got wrong is that the suite stayed green. The same method produced M1, Mi3, Mi8 and Mi9.
- **`cargo doc --no-deps --lib` is now part of this feature's validation set.** `Cargo.toml` denies broken intra-doc links and `fmt`/`clippy`/`test` do not cover rustdoc, so three deny-level errors passed the author's gate.
- **The five documentation-only files got Major-severity findings.** On a milestone whose deliverable is types and prose, the prose *is* the surface, and `generic/histogram.rs` was describing a data structure whose whole subtlety it had dropped.
