# A4 — the estimator reads the cohort by section, and not one number moves

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step A4 — implementation plan 2,
the last step of milestone A.
**Design authority:** [arch/parameter_prepass_joint_fit.md](../../ng/arch/parameter_prepass_joint_fit.md)
§2.1; [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§2.2.
**Date:** 2026-08-14.

---

## 1. What this step is for

**The estimator stops holding whole record sets.** `fit_jointly` takes
`&mut CohortCensusEvidence` and reads what it needs through one scoped call: the sections are
lent for the length of that call and taken back when it returns, so a file-backed census will
have nothing decoded afterwards. Nothing about the arithmetic changes, which is what §3 measures.

## 2. What changed

| | before | after |
|---|---|---|
| `fit_jointly` | `&[SampleCensusEvidence]` | `&mut CohortCensusEvidence` |
| `fit_contamination` | `&[SampleCensusEvidence]` | `&mut CohortCensusEvidence` |
| `ssr_fit::gather_strata` | `&[SampleCensusEvidence]` | `&mut CohortCensusEvidence` |
| the recording-terms refusal | made inside `fit_jointly` | made when the cohort is built |

**The refusal moved to the door and kept its name.** `CohortCensusEvidence::new` compares the
twelve terms before a section is decoded (A3); `impl From<TermsDisagreement> for JointFitError`
turns its refusal into `IdentityMismatch`, whose name
[`parameter_prepass_joint_loci.md`](../../ng/spec/parameter_prepass_joint_loci.md) specifies.

**Contamination is fitted inside the same scoped call as the alternation.** It reads the same
ordinary-position sections, so lending them twice would decode the generic half twice on a
file-backed census. Nothing else about it changes: it still runs *after* the alternation, on the
converged error rates and homozygote excess (spec §3.4).

### One thing this step does not do, and it belongs to a later measurement

**`gather_strata` is handed every stratum in one band.** The slippage fit borrows a thin stratum
from its neighbours across the whole set (`fit_strata`), so it needs them together as it is
written. How many strata may be resident at once — and how the fit runs inside a memory ceiling —
is [`parameter_prepass_joint_fit.md`](../../ng/spec/parameter_prepass_joint_fit.md) §11 questions
8 and 10, which the plan puts out of scope. What A4 delivers is that the access is *through* the
cohort's scoped call, so the band is a parameter rather than a property of the code.

## 3. What moved on real reads

**Nothing. Both oracles are byte-identical to A2's output**, wall-clock times masked
(`diff` of `tmp/a2_tomato.txt` against `tmp/a4_tomato.txt`, and the same for the trio, prints
nothing). That is the whole assertion this step needs: it is a change of access, and a change of
access that moved a fitted number would be a defect in the plumbing.

## 4. Tests

No test was added: this step moves signatures rather than behaviour, and the tests that already
exist are what say the behaviour is unchanged. Two were rewritten to reach the new door:

| test | what changed |
|---|---|
| `samples_that_disagree_on_the_ladder_are_refused_and_the_field_is_named` | the refusal now comes from `CohortCensusEvidence::new`, and the test converts it with `JointFitError::from` — so it still pins the variant *and* the field |
| the drawn-cohort fits in `fit.rs` and `contamination.rs` | wrap their samples with a one-line `as_cohort` helper |

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo check --all-targets` | 0 errors |
| `cargo clippy --lib --all-features -- -D warnings` | clean; the three examples this step touched are clean too |
| `cargo test --lib` | `3,592 passed; 0 failed; 11 ignored` — unchanged from A3 |
| the 88-second tomato oracle | §3 — byte-identical |
| the 74-second trio oracle | §3 — byte-identical |

**The two red gates are the two that were red before this branch's first commit**, neither in
code this plan touches: `cargo clippy --all-targets --all-features -- -D warnings`, and
`cargo test --all-targets`, which panics in `benches/psp_writer_perf.rs:386`.
