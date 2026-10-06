# Code Review: fit_precision_e2
**Date:** 2026-10-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step E2 — the parameters file, version 2: each number's standard error, the own-fit errors, the samples a stratum was fitted on, the SNP/indel fit's starts, and numbers with no information written `defaulted`
**Status:** Request-changes

---

## 1. Scope

- **Reviewed:** review object `c2f32f0c` against `25af661c` (branch `fit-precision`), plan step E2 of
  [fit_precision.md](../../implementation_plans/fit_precision.md); spec [fit_precision.md](../../ng/spec/fit_precision.md)
  §5.2–5.3, §8 question 1.
- **In scope:** every file of the diff (25), chiefly `src/calling/parameters_file/`, `src/parameter_estimation/joint/{fit,stratum_fits}.rs`,
  `src/run/census_fit.rs`, `src/calling/likelihood/mod.rs`, `src/cli/{calling_run,call_from_psps,call_from_alignments,estimate_parameters}.rs`.
- **Dispatched**, two reviewers, each in its own worktree at `c2f32f0c`:
  - reliability, float portability, refactor safety, and the parser/validator and stable-output items of extras;
  - errors, defaults, naming, idiomatic Rust, smells, module structure, and the prose the file carries.

## 2. Verdict

Request-changes: four Majors, one of them a confident-and-wrong number reachable at one sample.

## 3. Execution status

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean at `c2f32f0c`.
- Release lib tests of the touched modules: pass at `c2f32f0c` but for the cross-platform checksum, re-recorded with
  its proof.
- Mutations: 5 by the author before review, all killed; 5 by the reliability reviewer — 2 killed, 2 survived (M4 and
  Mi1 below), 1 changed no behaviour; 1 attempted by the prose reviewer and refused by the permission system.
- Every figure in the oracle baseline note and the checksum note was recomputed by both reviewers and holds.

## 4. Open questions and assumptions

1. Whether the cross-platform digest should pin the new errors (M3): needs a fixture that produces them and a
   re-record on macOS and Linux — a decision for the owner.

## 5. Top 3 priorities

1. M1 — one sample with reads beside a tracts-only sample gets its excess written `fitted_here`.
2. M2 — empty claimed-error totals are given to every read group, which turns a lost read-group axis into a silent
   multiplier of one.
3. M5 — two refusal messages carry a run of 26 spaces mid-sentence.

## 6. Findings

### Major

**M1: [fit.rs:2203](../../../../src/parameter_estimation/joint/fit.rs#L2203) — the excess is fitted across samples
no read informs.** Confidence High. *Reliability.* The rule "fit the excess at two samples or more" counts every sample,
at four places (the maximisation, both inversions, the warrant). One sample with reads beside one whose census holds
repeat tracts only is a one-sample cohort, and its excess came back `fitted_here` with an error. Before this step that
cohort could not reach a file (M2's panic). **Fix:** count only samples that inform their excess, at every place.

**M2: [census_fit.rs:415](../../../../src/run/census_fit.rs#L415) — empty totals for every read group.** Confidence
High. *Reliability.* The step needs them only for a read group nothing scored; given to all, a fitted rate that lost its
totals — the lost axis `checked_read_group_count_of` exists to catch — becomes the defaulted multiplier silently, and
the new run-level test's assertion on the tracts-only read group could not fail, since its drawn samples had no totals.
**Fix:** empty totals only for a defaulted rate from no observation; real totals for the test's drawn samples.

**M3: [cross_platform_digests.rs:107](../../../../src/cli/cross_platform_digests.rs#L107) — the portability digest
pins none of the new numbers.** Confidence High on what is pinned, Low on whether any of it differs between platforms.
*Reliability.* The fixture's file carries no `standard_error`: its rates sit on their bounds and its substitution
rates are zero. **Fix (owner's decision):** a fixture that produces errors, re-recorded on both platforms.

**M4: [calling_run.rs:705](../../../../src/cli/calling_run.rs#L705) — nothing tests that a calling run writes back the
starts it read.** Confidence High. *Reliability, convergent with the prose reviewer's builder finding.* Dropping them
survived every test. **Fix:** a test through the whole command, and one door both commands write through.

**M5: [validate.rs:151](../../../../src/calling/parameters_file/validate.rs#L151), `:170` — two refusals carry 26
spaces mid-sentence.** Confidence High. *Errors, convergent with the reliability reviewer.* A wrapped string lost its
`\`. **Fix:** the continuation, and `refused` asserting no message carries a run of spaces.

### Minor

- **Mi1** `validate.rs:157` — `agreed_with_start = 0` untested (mutant survived).
- **Mi2** `validate.rs` — accepts: a start agreeing with one that never converged; an empty starts list;
  `samples_fitted_on` above the cohort's size; a `standard_error` on the outlier weight or the fallback concentration,
  then dropped on write-back. (A version-1 file carrying version-2 keys is also accepted.)
- **Mi3** `to_toml.rs:341` — the subset sentence also fires where the subset grew to every sample, and gives no size.
- **Mi4** `to_toml.rs:106` — the header's account of a missing error omits the commonest case (a rate of exactly
  zero: 52 of 188 rows on the oracle) and whose error a `supplied` number carries.
- **Mi5** `validate.rs:163`, `:176`, `:1152`, `:803` — refusals that name the key but give no fix, or a wrong one.
- **Mi6** two new tests spliced between an existing test's doc comment and its `fn`
  (`likelihood/mod.rs:2654`, `from_run_parameters.rs:1631`).
- **Mi7** `StartEnded` (file) beside `StartEnding` (fit): two names a tense apart; not in either spelling test.
- **Mi8** `ReadsBehindEachCalibration::nothing_was_fitted` used by the fit's own writer, whose doc says every
  calibration is `Defaulted`.
- **Mi9** stale docs: the version check's "until there is one", `of_run`'s "six arguments", the module's example tree
  at version 1.
- **Mi10** `fit.rs:2191` and two test docs — "the start the fit never moved it from" is wrong: the accelerated steps
  move an uninformed excess, so writing the default value is load-bearing.
- **Mi11** `parameter_estimation/mod.rs:137` — `Estimate::standard_error` says a supplied value has none; a supplied
  number from a file keeps its error.

### Nits

Internal step labels in comments; the baseline citing an untracked scratch script; `_ => (None, None)` over
`StratumOutcome`; `expect`s without a note; "a tenth" written where `SETTLED_FRACTION` is meant; long inline paths in
tests; shorthand in `OwnFitStandardErrors`'s doc.

## 7. Out of scope observations

- `run_parameters.rs:210` — the "both halves by construction" comment is incomplete on the census path (goes with M2).
- A no-information coordinate drifts during the fit under acceleration (Mi10); pinning it during the fit rather than
  overwriting it afterwards is a design question for a later plan.

## 8. Missing tests to add now

- `one_sample_with_reads_beside_a_tracts_only_sample_holds_its_excess` (M1).
- `a_library_with_no_reads_beside_one_with_reads_comes_out_defaulted` (the A6 case beside a sibling with reads).
- `each_ending_is_written_with_its_number_and_passes` (the starts' mapping, `Agreed` included).
- `a_run_writes_back_the_starts_its_file_carried` (M4, through the whole command).

## 9. What's good

- The proof that no number moved: every version-2 key and note removed, both oracle files are the old ones line for
  line, and the fixture's file hashes to its old checksum — re-run by both reviewers.
- `version_1_as_written.toml` kept unregenerated, so a version-1 file is read as a real file, not a hand-made one.
- Own-fit errors filtered by the provenance source at the one place they are gathered (`StratumFits::over`).

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --release --lib -- calling::parameters_file run:: cli:: parameter_estimation::joint calling::likelihood`
- `scripts/promote_ng_oracle.sh` per `tmp/fit_precision/run_oracle_e2b.sh`

Per-category files: `tmp/review_2026-10-06_fit_precision_e2/`.
