# Fix Application Report: ng_calling_loop_a1_2026-08-25.md

**Date:** 2026-08-25
**Source review:** [ng_calling_loop_a1_2026-08-25.md](ng_calling_loop_a1_2026-08-25.md)
**Source state reviewed against:** `5843f60a` on branch `ng-calling-loop`, over branch point `bbcf2165`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 14
- Minors: 13
- Nits: 9 (grouped)

### Outcome totals
- Applied: 24
- Applied with adaptation: 3
- Already fixed: 0
- Deferred: 1
- Disputed: 1
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary

Every command in the container, from this worktree's own `scripts/dev.sh`.

- `cargo fmt --check` → **0**, no output
- `cargo clippy --lib --all-features --tests -- -D warnings` → **0**, no warnings
- `cargo test --lib` → **0**, `4517 passed; 0 failed; 14 ignored`
- `cargo test --release --lib ng::calling::tests` → **0**, `50 passed; 0 failed`
- `cargo doc --no-deps` → not run — no public item outside `src/ng/calling/mod.rs` changed, and the module's own doc links are exercised by the build's `rustdoc::broken_intra_doc_links` handling at lint level.
- `cargo audit` → not run — no dependency changed.
- Performance check → **not applicable**: nothing changed is reachable from any harness in `benches/`. `grep -rl 'ng::calling' benches/` returns nothing.

**The counts, against what:** the branch point had `4488 passed`; step A1 as reviewed had `4502`; after these fixes it is **`4517`**. So the review's fixes added **15** tests to A1's original 14.

**The release run is the one that matters here and it is new.** `cargo test --release --lib ng::calling::tests` is the only command that can tell `assert!` from `debug_assert!`, and it passes all 50 of the module's tests — including every one of the eleven `#[should_panic]` cases, which is what proves each new check is held in release rather than compiled out.

### Unresolved high-priority findings
- **M14** — no CI gate holds this module's assertions to release. Deferred; it cannot be added from this branch. See §5.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | repeat-bundle path cells untested | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M1 | Major | six accessors silent on an unprepared scratch | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M2 | Major | genotype table vs allele count asserted nowhere | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M3 | Major | `advance_expected_copies` returns stale values | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M4 | Major | empty contamination slice is a second spelling | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M5 | Major | cohort and sample copies share a name stem | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M6 | Major | `Missing` at a tract has no assertion | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M7 | Major | empty read-group list accepted | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M8 | Major | empty-evidence messages name no locus | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M9 | Major | `LocusEvidence::ssr`'s guard untested | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M10 | Major | infinite copy count untested | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M11 | Major | the two axes indistinguishable in the fixture | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M12 | Major | four buffers unreachable | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M13 | Major | shipped emission scratch never instantiated | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| M14 | Major | no release gate in CI | Defer | **Deferred** | No | None | N/A | Yes — §5 |
| Mi1 | Minor | `lg_table` naming, and the wrong quantity in its doc | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi2 | Minor | `FrozenParameters` fields name their topic | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi3 | Minor | three undefined words | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi4 | Minor | `fill_poisoned`, and the unnamed `NaN` | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi5 | Minor | `selection_mut` is half a name | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi6 | Minor | `Missing` names the output symbol | Ask→resolve | **Applied with adaptation** | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi7 | Minor | default type parameter picks an arm invisibly | Apply | **Applied with adaptation** | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi8 | Minor | empty gather has no expression | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi9 | Minor | overflow messages name no operand | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi10 | Minor | called allele ids never range-checked | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Mi11 | Minor | report uses bare step codes | Apply | Applied | No | impl report | N/A | No |
| Mi12 | Minor | report's "both paths" overstates coverage | Apply | Applied | No | impl report | N/A | No |
| Mi13 | Minor | three small coverage gaps | Apply | Applied | No | `src/ng/calling/mod.rs` | Pass | No |
| Nits | Nit | nine grouped | Apply | **Applied with adaptation** (8 of 9) | No | `src/ng/calling/mod.rs` | Pass | No |
| Nit (`GenericLocusSample`) | Nit | name is the pair's coordinates | Dispute | **Disputed** | No | None | N/A | No |

## 3. Questions asked and answers

None. The review's one open question (§4.1, whether `SampleGenotypeCall::Missing` is the right name) was resolved against the spec's own vocabulary rather than referred — see Mi6.

## 4. Per-finding log

Grouped where the fix is one edit serving several findings; every finding appears.

### B1 — repeat-bundle cells of the path matrix untested
- **Final status:** Applied. **Review suggestion used verbatim?** No — adapted.
- **Implementation:** two tests, `ssr_evidence_against_a_bundle_allele_table_is_accepted` and `generic_evidence_against_a_bundle_allele_table_is_refused`. The adaptation: the panic message now says *"its allele table is a repeat-bundle locus"* rather than printing the kind with `{:?}`, so the refusing test's message names the bundle in words (see Nits).
- **Verification:** both pass; the second is `#[should_panic(expected = "belong to different loci")]` and passes under `--release` too, so the assertion it rests on is release-held.
- **Residual risk:** the matrix is now 5 of 6 cells. The sixth — repeat evidence at a SNP/indel allele table — is the symmetric refusal and is not covered; it shares its one `matches!` arm with the covered `(Generic, Ssr)` case, so no mutation distinguishes them.

### M1, M12, M13 — the scratch's unreachable and unguarded surface
- **Final status:** Applied.
- **Implementation:** `assert_prepared()` is called first in every accessor, including through the two row-range helpers, so all fourteen refuse an unsized scratch by name rather than handing back an empty slice. The four buffers that had no accessor — the prior row, the posterior row, the sample concentration and the prior's per-allele workspace — gained read and write pairs. `CallingScratch::default()` remains the only constructor, which is right: a worker allocates its scratch before it has met a locus.
- **Tests added:** `an_unprepared_scratch_is_refused`; `the_prior_and_posterior_buffers_are_sized_on_their_own_axes`, which asserts the two row buffers are 3 long and the two per-allele buffers 2 at a fixture where `assert_ne!(genotype_count, allele_count)` holds, so a buffer sized on the wrong axis fails; `the_shipped_emission_scratch_builds_and_can_cross_a_worker_boundary`, which builds `CallingScratch<StutterSubstitutionScratch>` — the configuration a run uses and no test built — touches all three sub-scratches, and pins `Send`.

### M2 — the genotype table and the allele count
- **Final status:** Applied. **Adaptation:** none; the reviewer's signature was taken as proposed.
- **Implementation:** `prepare_for_locus(sample_count, alleles, genotypes)` asserts `genotypes.allele_count() == alleles.len()` before anything else. Placed here rather than at A2's `call_locus` because this is the one point where the locus's shape is fixed; the review's open question §4.3 records the alternative.
- **Test:** `a_genotype_table_built_for_another_allele_set_is_refused` admits a third allele to the table and hands in the two-allele genotype table — the discovery-round shape, not a hypothetical.

### M3 — `advance_cohort_expected_copies`
- **Final status:** Applied.
- **Implementation:** fills the returned buffer with `UNWRITTEN_SCRATCH_VALUE` on the way out, and the doc states a fact rather than a requirement. The cost the original comment gave for not filling does not survive the arithmetic, as the review says: the fill is one write per allele against a pass costing samples × genotypes evaluations.
- **Test:** the existing advance test now asserts the buffer arrives unwritten; `preparing_a_locus_poisons_the_previous_passs_copies_as_well` covers the other half — the previous-pass buffer surviving into the next locus.

### M4, M7 — the two axes of `FrozenParameters`, and how *absent* is spelled
- **Final status:** Applied.
- **Implementation:** `new` now refuses an empty contamination list, naming `FrozenParameters::uncontaminated` — which is added, alongside `contamination_is_absent()`. Both constructors funnel through one private `gather`, which is where the two non-empty-axis assertions live, so neither door can drop one. This is the spelling `ContaminationMixture::uncontaminated` / `is_absent` already carries for the per-locus half of the same mixture.
- **Tests:** `a_run_with_no_contamination_fitted_says_so_by_name`, `an_empty_contamination_list_is_refused_in_favour_of_the_named_constructor`, `frozen_parameters_refuse_a_run_with_no_read_groups`.
- **Note:** the test helper now routes to whichever constructor the contamination list calls for, so the existing fixtures did not have to choose.

### M5, Mi1, Mi2, Mi3, Mi4, Mi5 — the naming pass
- **Final status:** Applied.
- **Implementation, the renames that carry weight:**
  - the cohort's and one sample's expected copies now say whose they are — `cohort_expected_copies`, `sample_expected_copies`, `previous_cohort_expected_copies`, `advance_cohort_expected_copies`, `per_sample_expected_copies` — reusing the word `LocusInference::cohort_expected_copies` already uses;
  - `lg_table` → `genotype_likelihoods`, `lg_row` → `sample_genotype_likelihoods`. **And the doc comments were wrong, which is the half that mattered**: they called the buffer "the read likelihood", where the read-likelihood spec reserves that for one read against one allele. What it holds is the *genotype* likelihood. `Lg` survives only as a cross-reference in prose;
  - `FrozenParameters`' fields carry their axis: `calibration_by_read_group`, `contamination_by_read_group`, `inbreeding_coefficient_by_sample`, `prior_seed`, `ssr_slippage_fits`;
  - `fill_poisoned` → `resize_and_fill`, and the value it fills with is now the named `UNWRITTEN_SCRATCH_VALUE`, referenced from every accessor that can hand a caller an unwritten buffer;
  - `selection_mut` → `candidate_selection_mut`; `prior_allele_scratch` → `prior_row_workspace`.
- **Doc fixes:** *concentration* is defined at first use (chromosomes the prior behaves as though it had already seen), *stratum* likewise (tracts sharing a motif length and a repeat count), and *fold* is replaced by "per-allele running totals".

### M6, Mi10 — the two rulings `LocusInference::new` did not enforce
- **Final status:** Applied.
- **Implementation:** two assertions beside the existing gene-diversity one. A `Missing` call at a repeat tract or bundle is refused — the mirror of the marker that is refused at a SNP/indel locus, one ruling per path. And every called genotype's allele ids are checked against the table it was called over, which catches the stale-after-prune half.
- **Tests:** `a_repeat_tract_locus_cannot_carry_a_missing_call`, `a_call_naming_an_allele_the_locus_lost_is_refused`.
- **Residual risk, stated by the review and not closed:** an id that stays in range after a renumber names a *different* allele silently. The fix for that is the prune returning its remapping, which is the step that builds the prune, not this one.

### M8, Mi9 — panic messages that named nothing
- **Final status:** Applied. Both empty-evidence messages now interpolate the locus's region and say which path they are on; the two overflow messages name both operands and what the table would have been.

### M9, M10, M11, Mi13 — the coverage gaps
- **Final status:** Applied.
- **Tests:** `ssr_evidence_naming_no_sample_at_all_is_refused` (the one assertion in the module that the 16-way downgrade experiment found untested); `expected_copies_reject_an_infinite_count`; `read_group_count_and_sample_count_are_two_different_axes`, built at one library and three samples so the two axes separate; `the_allele_table_and_its_copies_are_never_empty`; the repeat arm's `region()` asserted alongside the SNP/indel arm's; and `Send` pinned in the shipped-scratch test.

### Mi6 — `SampleGenotypeCall::Missing`
- **Final status:** **Applied with adaptation.** The name is kept; the hazard is closed in the doc.
- **Reasoning:** the reviewer's substantive point is right — a geneticist may read `Missing` as *missing data*, and a zero-coverage sample is `Called`, so the natural reading is backwards. But the enum qualifies the *call*, not the data, and *missing* is the spec's own word for what emission writes (§9). Renaming to `Uncallable` would move the code's vocabulary away from the design's for a hazard that a sentence closes. The doc now says outright that a sample with no reads is `Called`, and that `Missing` is not missing data but *the caller declined to invent a genotype over a set that cannot hold this sample's allele*.
- **Verification:** the review's own verification step — confirm zero-coverage samples really are `Called` — checks out: `GenericSampleEvidence::empty()` scores every genotype alike and the prior decides alone (spec §7).

### Mi7 — the default type parameter
- **Final status:** **Applied with adaptation** — the reviewer offered two fixes and preferred the cheaper one; the stricter one was taken.
- **Implementation:** the default is **dropped**. `CallingScratch<SsrEmissionScratch>` now matches `SsrRowScratch<ModelScratch>`, the type it wraps, so the seam has one convention rather than two, and every construction names the emission model it is for. The cost the reviewer priced — four test call sites — was nil, since all of them already wrote the parameter.

### Mi8 — the empty slippage gather
- **Final status:** Applied. The constructor's doc bullet now names the expression, `StratumFits::over(&[], BTreeMap::new())`, rather than describing it. The better home is a named constructor in `parameter_estimation/joint/`, which is out of this session's scope and is routed in the review's §7.

### Mi11, Mi12 — the implementation report
- **Final status:** Applied. Three bare plan step codes are glossed, the two fan-out plans are named rather than counted, and the test-table row that read as full coverage of the discriminant now says which two of the three kinds it covers and points at the Blocker for the third. A note at the top of the report sends the reader here for what the review changed.

### Nits
- **Final status:** Applied with adaptation — eight of nine.
- Applied: `prepare_for` → `prepare_for_locus`; `assert_agrees_with` → `assert_matches_locus_and_run`; `prior_allele_scratch` → `prior_row_workspace`; `row_range` is private behind two typed wrappers, `genotype_row_range` and `allele_row_range`, so no caller picks the width, and its message names which table was indexed; `LocusKind` is rendered as a word rather than with `{:?}`, which for a tract printed both flanks as decimal byte arrays; the three test fixtures are `diploid_ploidy`, `outbred_samples`, `frozen_parameters`; the advance test's name loses its two dropped apostrophes; and *seam* is replaced by "the boundary where calling begins".
- **Disputed — `GenericLocusSample`.** The nit reads the name as the pair's coordinates where the value is evidence plus a ruling. It is *a sample at this generic locus*, which is what the value is, and the doc's first sentence names both halves. The proposed alternatives (`RuledSampleEvidence`, `GenericSampleEvidenceAndRuling`) are longer without saying more. Kept.

### M14 — no release gate — **Deferred**
- **Reasoning:** the fix is a CI step, `cargo test --release --lib ng::calling`, and it cannot be added from this branch. At `5843f60a` that command is `461 passed; 4 failed`, and all four failures are `#[should_panic]` tests in `src/ng/calling/likelihood/` — the `ng-calling-likelihoods` branch's file, which this session must not edit. Adding the gate now would land a red CI step.
- **What this session did instead:** ran the narrower `cargo test --release --lib ng::calling::tests` as part of this fix run's validation, where all 50 pass. That establishes the property for today's tree and is recorded in §1; it does not stop a later edit from downgrading a check.
- **Follow-up:** raise the four `--release` failures with the `ng-calling-likelihoods` branch, then add the gate.

## 5. Deferred findings to carry forward
- **M14** — no CI gate holds this module's assertions to release. Blocked on four pre-existing `--release` failures in `src/ng/calling/likelihood/`.

## 6. Disputed findings to return to reviewer
- **Nit (`GenericLocusSample`)** — the name says what the value is; the alternatives are longer without being clearer.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing changed is reachable from any harness in `benches/`; `grep -rl 'ng::calling' benches/` returns nothing.
- **Baseline saved:** not applicable.
- **Outcome:** skipped.

## 10. Commands run
- `./scripts/dev.sh cargo fmt`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features --tests -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::tests`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --release --lib ng::calling::tests`

## 11. Command results
- `cargo fmt --check` → 0, no output
- `cargo clippy --lib --all-features --tests -- -D warnings` → 0, no warnings
- `cargo test --lib ng::calling::tests` → 0, `50 passed; 0 failed`
- `cargo test --lib` → 0, `4517 passed; 0 failed; 14 ignored`
- `cargo test --release --lib ng::calling::tests` → 0, `50 passed; 0 failed`

## 12. Notes
- **The branch is one compiler pin behind `main`.** `main` moved to 1.98 on 2026-08-25 (`54a0fd96`) and made `cargo clippy --all-targets --all-features -- -D warnings` exit 0 (`f3c8c797`); this branch still pins 1.97.1, where that wider scope is red for reasons under `examples/` and `benches/`. Taking main's pin, and then running the wider clippy scope, is the next thing to do on this branch and is deliberately not bundled into the step's own commit.
- **One review worktree's directory survived cleanup** at `.claude/worktrees/agent-a96585d42bd9ce617` — git no longer tracks it, but about 3.6 GB of build output is still on disk. Removing it needs a `rm -rf`, which this session is not permitted to run.
