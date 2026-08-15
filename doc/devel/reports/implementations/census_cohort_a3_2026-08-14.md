# A3 — the whole cohort, refused at the door and lent a band at a time

**Plan:** [census_file.md](../../ng/impl_plan/census_file.md) step A3 — implementation plan 2.
**Design authority:** [arch/parameter_prepass_joint_records.md](../../ng/arch/parameter_prepass_joint_records.md)
§2.2; [spec/parameter_prepass_joint_records.md](../../ng/spec/parameter_prepass_joint_records.md)
§5, §6.2; [spec/parameter_prepass_joint_loci.md](../../ng/spec/parameter_prepass_joint_loci.md) §3.6.
**Date:** 2026-08-14.

---

## 1. What this step is for

**A tract's length frequencies are fitted from every sample with reads there**, so the unit the
fit borrows is one band of strata across the *whole cohort* — not one sample's stratum. That puts
the scoped calls on a type that owns every sample, which is what `CohortCensusEvidence` is.

**And the twelve recording terms are checked before a single section is read.** They say which loci
were asked for, which came back, and in what units the evidence was written down; two samples that
disagree on any of them hold rows meaning different things, and every one of them fails silently.
The refusal is at the door, where it costs nothing, rather than after a pass over two million
positions.

## 2. What landed

| item | what it is |
|---|---|
| `CohortCensusEvidence::new` | adopts a `Vec<SampleCensusEvidence>`, refusing on the first disagreement |
| `TermsDisagreement` | which two samples, and the value they first differ on, by name |
| `read_groups`, `strata` | the union across samples, in the order the fit visits them |
| `with_generic(groups, f)` | every sample's ordinary-position sections for a band of read groups |
| `with_strata(strata, f)` | every sample's tracts for a band of strata |
| `SampleTractSections<'a>` | one sample's `(read group, stratum, section)` triples, as a call lends them |

**A sample answers with the sections it holds**, not with one per group asked for: a cohort's
samples need not have been sequenced the same way, and `fit.rs` has always read each sample's own
read groups. A band naming a stratum the cohort does not hold lends nothing, rather than a zero
that would read as a sample with no reads there.

**`TermsDisagreement` rather than a variant of `CensusError`.** `CensusError`'s two variants
(arch §2.3) both need the file, which does not exist until milestone B, and the refusal this
carries is spec §5's, whose name in the fit is fixed by
[`parameter_prepass_joint_loci.md`](../../ng/spec/parameter_prepass_joint_loci.md) —
`JointFitError::IdentityMismatch`. A4 builds that variant from this one, so the specified name
survives and the census does not depend on the module that does the mathematics.

## 3. Nothing calls it yet, so nothing can move

**This step is additive: no oracle run was made, because no code path reaches the new type.**
`fit_jointly`, `fit_contamination` and `gather_strata` still take `&[SampleCensusEvidence]`; A4 is
where they move onto the cohort and where both oracles are run again. `grep -rn
CohortCensusEvidence src/ examples/` finds it only in `census.rs`.

## 4. Tests

| test | what it pins |
|---|---|
| `a_cohort_refuses_two_samples_that_did_not_record_the_same_thing` | the refusal names both samples and the value — two depth ladders, every other term agreeing |
| `a_cohort_lends_a_band_across_every_sample_and_nothing_outside_it` | one row a sample rather than one value; the union of read groups and strata; a band naming a stratum the cohort does not hold lends nothing |

## 5. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::parameter_estimation::joint::census` | `43 passed; 0 failed` (41 before) |
| `cargo test --lib` | `3,592 passed; 0 failed; 11 ignored` (3,590 before) |
