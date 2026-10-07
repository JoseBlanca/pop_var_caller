# Code Review: fit_precision_zero_substitution_rate
**Date:** 2026-10-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** a repeat-tract substitution rate is never exactly zero or one (owner, checkpoint E of `fit_precision.md`)
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** review object `d696f374` against `479e5ddb`, branch `fit-precision`: `StratumSubstitutionCounts::substitution_rate`,
  its error in `parameters_from_the_fit`, the docs and note it changes, the golden files and the re-recorded checksums.
- **Dispatched:** one reviewer, in its own worktree, over reliability, float portability, errors, naming, smells and the
  prose of the comments and notes — the change is one function and its readers.

## 2. Verdict

Approve-with-changes; one Major for the owner.

## 3. Execution status

- The five targeted tests pass at `d696f374` (`--release`, Linux container).
- Mutations: 4 run, 4 killed.
- The claim that the old rule reproduces the old calls checksum was re-checked independently: `2c0a8946…` exactly.
- Not re-run by the reviewer: the oracle figures (taken from the author's run, `tmp/fit_precision/oracle_e3/`).

## 4. Open questions and assumptions

1. M1: refuse a rate of exactly zero or one in a file, or convert it on reading — the owner's decision.

## 5. Top 3 priorities

1. M1 — a zero still reaches calling from a version-1 or hand-edited file.
2. Mi1 — the rule's own doc mis-states its reason and an example.
3. Mi3 — `>=` turns an impossible count into a plausible rate.

## 6. Findings

### Major

**M1: [validate.rs:884](../../../../src/calling/parameters_file/validate.rs#L884) — a zero rate still reaches calling
from a file.** Confidence High. `validate` checks a substitution rate only for `[0, 1]`, and every version-1 file holds
the zeros the old rule wrote (52 of 188 on the oracle); calling then gives a read with one mismatched base about e⁻⁷⁰⁸
under every length (`FlatEmission` floors at `f64::MIN_POSITIVE`), so it falls to the outlier term. **Options:** refuse
0 and 1 in `validate`, the message giving the replacement from the row's `bases_compared`; or rewrite a zero with a count
to `0.5 / (n + 1)` on reading, which is exact; a zero with no count must be refused either way.

### Minor

- **Mi1** [ssr_fit.rs:236](../../../../src/parameter_estimation/joint/ssr_fit.rs#L236) — the doc calls the rate a
  "prior" and its example "one over 34" reads as one mismatch in 34 bases, which the rule leaves alone.
- **Mi2** [ssr_fit.rs:336](../../../../src/parameter_estimation/joint/ssr_fit.rs#L336) and spec
  `parameter_prepass_ssr.md` §4.2 — still "one division".
- **Mi3** [ssr_fit.rs:252](../../../../src/parameter_estimation/joint/ssr_fit.rs#L252) — `>=` gives an impossible count
  (more mismatching bases than compared) a fitted rate near one, where it used to be refused as above one.
- **Mi4** [promote_ng_oracle.baseline:23](../../../../scripts/promote_ng_oracle.baseline#L23) — the 52 rows also gain a
  standard error, which the note omits.
- **Mi5** the rule is tested inside two tests named for other behaviours.

### Nits

The locals `compared` and `compared_f`.

## 7. Out of scope observations

- `types.rs:1197` — `DomainError::SsrBaseComparison` cites a type that no longer exists and nothing constructs it.
- `fit_precision.md` spec §5.2's table cites "§3.7", which is the parameters-file spec's section, not this spec's.

## 8. Missing tests to add now

`a_count_with_one_outcome_only_takes_half_a_count_of_the_other` (Mi5); for M1, per the owner's choice,
`a_substitution_rate_of_zero_or_one_is_refused_with_its_replacement` or its conversion counterpart.

## 9. What's good

- The fixture's calls change was attributed by putting the old rule back and reproducing the old checksum, not argued.

## 10. Commands to re-verify

`./scripts/dev.sh cargo test --release --lib -- substitution cross_platform_digests run::census_fit`
