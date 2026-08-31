# ng parameters file — F1 design-fidelity review

**Date:** 2026-08-31
**Step:** [parameters_file.md](../../ng/impl_plan/parameters_file.md) F1
**Design authority:** [parameters_file.md](../../ng/spec/parameters_file.md) §1.2, §2, §2.1, §3, §5,
§6, §7, §8, §10; [run_streaming.md](../../ng/spec/run_streaming.md) §2; `CLAUDE.md` §"the range,
not the example"
**How it was run:** one agent, in a worktree detached at `ede29317` with the step applied,
building through the worktree's own `scripts/dev.sh` (220 passed, 0 failed).

**Verdict: mostly faithful, with one design-level defect and two seam problems in the new optional
census.**

---

## Blocker — the file's opening sentence contradicts the warrants under it

Only two of the seven `was_fitted` predicates read a warrant; four read row presence and one reads
a rung. `demoted_to_no_better_than_supplied` moves warrants and nothing else, by design. So a file
demoted under §2.1 opened with *"5 of the 7 groups … were fitted from your data"* over rows whose
every warrant is `supplied`, while `WhatTheRunFitted`'s own doc promised the summary *"cannot come
to disagree with the warrants beside the numbers"*.

**The reviewer split it correctly, and the split is what the fix follows.** The substitution-rate
rows *do* carry a warrant (`SubstitutionRateRow.rate: WarrantedValue`, demoted at
`bindings.rs`), so that predicate was simply reading the wrong thing — a one-line fix. The other
four carry no *handed-over* state at all, which is the D3 gap already recorded in
`PROJECT_STATUS.md` for `SeedRung::FittedCurve`; F1 had built a headline on top of it and asserted
it did not exist.

**Both other reviewers found the same thing independently**, from the code's side and from the
artefact's. See the fixes report.

## Major — the `Option<&CensusIdentity>` ruling holds, but the spec does not state it

The reviewer checked §2.1 and §6 and agrees the reading holds: §2.1's *rule* sentence and §6's
fourth bullet are both about a **mismatch** and neither covers *there is nothing to compare
against*, but §2.1's **rationale** names and rejects "demote on every read" precisely because it
breaks the two-mode oracle — and direct mode is the mode with no census, so demoting on absence is
that rejected alternative under another name.

**It offers a second argument, stronger than the one the code gave, and the code now carries it:**
under the alternative, **warrants decay on every round trip.** A psp run fits and writes
`fitted_here`; a direct run reads it, reports everything `supplied`, and §7 makes it write that
back. §7's *a run is reproducible from its own output* would be true of the genotypes and false of
the warrants after one hop.

**And it raises something larger, which is the owner's.** `run_streaming.md` §2 gives psp mode's
**calling** stage the psps only — census files belong to the walk and the fit — so unless a driver
deliberately opens them, a separately-invoked psp-mode call has no `CensusIdentity` either. On
today's spec the fourth binding may be reachable from **no real run at all**, and §13 test 5 asks
for the demotion to be visible in a run's report. Recorded for Checkpoint F; the recommendation is
one sentence in §2.1 or §6.

## Major — `None` conflated "no census" with "the censuses agreed"

`ParametersForThisRun.fitted_under_another_census: Option<String>`, `None` documented as *they
agree*. After F1 `None` also meant *no comparison was made*, so a driver or F2's report reading
`is_none()` would tell the user the file matched this run's census when no census existed. Three
states in, two out — and the step's own test pinned the conflation rather than catching it.
**Fixed** with a three-state `CensusAgreement`.

## Major — `Some(of_a_run_with_no_census())` and `None` gave opposite outcomes

The writer mints an empty-term identity for *this run had no census* and takes it by value; the
reader took `Option<&CensusIdentity>` where `None` keeps warrants. A driver holding one identity —
the natural shape, since `of_run` takes one — would hand the reader the writer's spelling and get
everything demoted. **Fixed**: an identity naming no terms is treated as no identity, with a test.

## Major — "not fitted" was the wrong thing to say about an absent contamination table

The unfitted-groups sentence ended *"so every number in those is a constant compiled into this
caller or a value somebody handed it"*. For an absent `[contamination]` section there are no
numbers and no constant: §3.4, §5's first row and §8 all say absence is a **real model state**.
**Fixed** — the blanket clause is gone and each section says what it holds instead.

## Major — the header's measured claim was wrong

Refuted the "eleven lines" half correctly: `to_toml`'s opening note is unconditional and renders
**39 lines** for every file. Its replacement figure — 106 identical lines, first difference at 107
— was measured against `ng_parameters_file_e2_defaults_run_as_written_2026-08-31.toml`, an artefact
written at **E2** (`6e434561`) whose preamble E3 later changed; that comparison actually first
differs at line 23. **Measured directly on the base commit's own golden**
(`git show ede29317:…/every_shape_as_written.toml`): the first `warrant` value is on **line 105**,
which is what the source already said. 39 stands; 105 stands; 107 does not.

## Minor

- *"a run scoring from a file could not write the file it used at all"* overclaims — it could, by
  inventing a placeholder error rate, which is what the retired helper did. **Fixed.**
- `RepeatTractLengthSpectra::key()` named one of the group's two tables, so a run whose spectra
  came off its periods' curves is pointed at the empty one. **Fixed.**
- `write_beside_the_vcf` overwrites silently, and §7 tells users to edit the file their last run
  wrote — a re-run whose supplied input and whose VCF share a stem overwrites its own input.
  **Documented**, and the write is atomic now so the replacement is at least whole; whether a
  driver should refuse it is the driver's.
- **Range.** `to_toml` builds the whole file as one `String`, and §9 prices the substitution-rate
  axis at up to 62 MB at 3,000 samples — C4 re-measured 185 bytes a row, nearer 79 MB. So this is
  a single allocation of that size taken *after* the last locus. **Recorded** in the module header.
- Every `was_fitted` predicate is `any()` over a cohort-sized axis, so at 3,000 samples one fitted
  row in 3,000 makes a group "fitted". The sentence that keeps it honest is on all three arms.
- `to_string_lossy` and paths with no file name — the correctness review's, and fixed there.

## Scope

**Clean on F2**: `run_report.rs` is untouched. **Nothing left for F2 that is F1's**, given no run
driver exists. **The naming rule is a fair reading of "beside its VCF"** and F1 is its right home,
since the driver plan defers to §7; the reservation is that it lives only in a module header, and
it is a format decision a user will depend on.

**The prose work is authorised by E3's handback rather than by F1's plan line** — correct, but the
step's largest new user-facing surface comes from a deferral, and the plan's landed note says so.

**`ReadsBehindEachCalibration` is a real seam, not a rename.** It moves the two fit-specific
invariants into the fit's own constructor and retires an in-band sentinel §5's own rule forbids —
before this, an empty `BTreeMap` *meant* the run had no fit. No fourth source was found that the
three constructors cannot express.

**The seven groups are defensible against §3.** The two exclusions rest on the spec's own words;
`RepeatTractLengthSpectra` as one group is right, since §3.7's top and middle rungs are two answers
to one question; `OrdinarySitePrior`'s predicate matches `SeedRung`'s own documentation; `Borrowed`
counting as measured is well argued from §2's ladder.

**Still not met, and correctly recorded:** §7's *writing is unconditional* fails for a run whose
reference came from a `.fai` alone. The step names it, argues why the digest must not be made
optional (§6's strongest binding refuses rather than demotes), and pushes it to the run driver.
