# Code Review: ng_paralog_filter_c3

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), three sub-agents in isolated worktrees
**Scope:** step C3 of the hidden-duplication filter plan — pass two: scoring the parked records and resolving the cut
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `14c64563` on branch `ng-paralog-filter` — one new module and its
  tests, one method added to an existing type with four accessors narrowed, and the module's
  declarations. Merge base `546845b2`.
- **In-scope files:** [pass_two.rs](../../../../src/ng/run/paralog_filter/pass_two.rs),
  [pass_two/tests.rs](../../../../src/ng/run/paralog_filter/pass_two/tests.rs),
  [scoring_context.rs](../../../../src/ng/run/paralog_filter/scoring_context.rs),
  [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs), the step's implementation report, and
  the `PROJECT_STATUS.md` entry.
- **Out of scope:** `src/ng/paralog/*.rs` (guarded copies — findings raised, never edited);
  production; `src/ng/window_coverage/`; steps A1–C2; pass three, the CLI flags' behaviour, the run
  report and the header line, all C4's.
- **Categories dispatched:** reliability + extras (mutation testing and the test challenge);
  errors + defaults + refactor_safety + smells; naming + idiomatic + module_structure + the diff's
  own numbers. `unsafe_concurrency` skipped — no `unsafe`, no threads, no shared state.
  `tooling` skipped — `Cargo.toml` untouched.

### 2. Verdict

**Request-changes.** One Blocker, eight Major, sixteen Minor.

**The theme is that the step guards the ratio's *value* thoroughly and never checks what the
verdict does with it.** Thirteen tests cover two thirds of the step's own stated contract — an
unscored record is *kept* and *out of the fit* — and none covers the third, that it is *never
flagged*. Measured, that third is closer to failing than any of them: the false-discovery curve's
own answer for a value that is not a number is `0.4999999999999852`, which sits **inside** a target
of one in two. The only thing holding an unscorable record out of the removed set is the finiteness
screen in the frozen `ParalogCalibration::flags`, and this suite could not see it working.

**And three fixtures were uniform in dimensions the pass reads.** Every record sat on contig 0, so
a failure that reported a hardcoded first contig passed all thirteen — which is C2's Blocker
repeating, in a step whose own module doc claims to enumerate its constant dimensions and names
only one of four.

### 3. Execution status

| command | at `14c64563` | at the merge base `546845b2` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6595 passed; 0 failed; 15 ignored` | `6582 passed; 0 failed; 15 ignored` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9, same kinds |
| `cargo fmt --check` | dirty on 9 files, none this step's | the same 9 |
| `cargo test --lib --all-features ng::run::paralog_filter` | `ok. 121 passed; 0 failed` | — |

One integration test fails at both — `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`,
which is `main`'s. `--all-targets`, `cargo doc` and `cargo audit` not run, for B1's reasons
(`main` is red; `Cargo.toml` untouched).

**Findings labelled "Needs verification": 0.** Every finding came from a mutation or a probe run in
a reviewer's own worktree.

**Mutation totals across the fan-out: 6 run by the reviewers, 4 survived, 0 changed no behaviour**,
on top of the author's 8 run and 8 killed. One of the reviewers' six re-verified the author's
central deviation by asking it backwards — installing *production's* screening rule in place of
ng's — and one test killed it.

### 4. Open questions and assumptions

1. **Does pass three need random access to the ratios, or only one per record as it streams?**
   Bears on M6. Answered by the author during synthesis: spec §3.5 walks the spill "per entry, in
   spill order, with its ratio from pass two", so it streams, and a consuming cursor is the right
   shape. Deferred to C4 rather than applied, because spec §3.7's type block declares the field.
2. **May `src/ng/paralog/` gain a new file?** Bears on M7. Yes — only the five copied files and the
   span in `calibration.rs` are byte-guarded, and the module already holds three files of ng's own.
   The author's stated reason for placing the calibration elsewhere was therefore wrong.
3. **Is the histogram's fixed ±100 range right at a large cohort?** Bears on Mi8. Unresolved and
   not this step's to resolve; raised for Checkpoint C and made measurable here.

### 5. Top 3 priorities

1. **B1** — assert the "never flagged" third of the contract; it is inside the cut by 1.5 × 10⁻¹⁴.
2. **M2** — the fallback rate and the iteration's seed are the same number, so two mutations that
   read the wrong one survive the whole suite.
3. **M4** — every fixture on contig 0; the failure that names a record can name the wrong contig
   and pass. C2's Blocker, repeating.

### 6. Findings

#### Blocker

**B1: src/ng/run/paralog_filter/pass_two/tests.rs — nothing asserts that an unscored record is never flagged, and it is inside the cut by 1.5 × 10⁻¹⁴**

**Confidence:** High. **Categories:** reliability.

The module doc, the commit message and the report all state the contract three ways: an unscored
record is *kept*, *never flagged*, and *out of the fit*. Thirteen tests cover the first and third.
**No test calls `flags` on an unscored ratio.** Measured on the mixed fixture at a target of one in
two:

```
PROBE ratios=[-31.149604210372544, NaN, 156.2578875006895] flagged=[true, false, true] q_of_lr(NaN)=0.4999999999999852 target=0.5
```

`curve.q_of_lr(NaN) <= target_fdr` is **true**. The only thing keeping that record out of the
removed set is the `lr.is_finite() &&` conjunct in `ParalogCalibration::flags`
([calibration.rs:92](../../../../src/ng/paralog/calibration.rs#L92)), which is in a frozen file this
pass does not own. If that screen ever moves, a record no sample could speak for is dropped from
the VCF as a hidden duplication, nothing panics, and `records_in_the_fit` still says it was never
scored — spec §6 trap 4 arriving through the other door.

**Fix:** `an_unscored_record_is_never_flagged_however_loose_the_target`, asserting across four
targets up to the loosest the type admits, plus a companion assertion that both *scored* records
are removed at that target so the test cannot pass by removing nothing.

#### Major

**M1: pass_two.rs — `target_fdr` was an unvalidated `f64`, and both wrong ends failed silently**

**Confidence:** High. **Categories:** errors, defaults.

Measured on a three-record spill over six samples:

```
PROBE target_fdr NaN => threshold None, flagged 0/3
PROBE target_fdr -1  => threshold None, flagged 0/3
PROBE target_fdr 5   => threshold Some(-99.95), flagged 3/3
```

A target of `5` — what someone typing five for "five percent" gets — removes every scored record. A
negative or `NaN` one removes none and reports no cut, which is the same thing the calibration
reports when a target is simply unreachable, so a run report cannot tell a wrong knob from a clean
cohort.

**Fix:** `TargetFdr`, a newtype admitting `[0, 1)` — the range `--paralog-fdr` already parses —
mirroring `InbreedingF`, with `NotATargetFdr` naming the value.

**M2: pass_two.rs — no test could tell the configured fallback rate from the iteration's seed, because both defaults are `0.03`**

**Confidence:** High. **Categories:** reliability.

`DEFAULT_EM_START` and `DEFAULT_FALLBACK_PARALOG_PRIOR` are the same number, and an empty
histogram's estimate *is* the seed. Two mutations survived all thirteen tests: reading
`config.em.start` instead of `config.fallback_prior`, and skipping the substitution entirely when
the histogram is empty — the exact case `an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing`
exists to cover. Both give `0.17` where `0.41` was configured.

**Fix:** the two tests that assert the fallback now run on a configuration whose fallback (`0.41`)
is not its seed (`0.17`); a new `the_shipped_fallback_rate_is_the_documented_one` pins the default.

**M3: pass_two.rs — the warning's record count was only checked where the two counts are equal**

**Confidence:** High. **Categories:** reliability.

`records_in_the_fit` earns its place precisely because it differs from the parked count when
records could not be scored. The only test reading the warning's number used a spill in which both
records score, so replacing `records_in_the_fit` with `ratios.len()` survived all thirteen.

**Fix:** `the_warning_counts_the_records_in_the_fit_and_not_the_records_parked`, on a three-record
spill of which one cannot be scored, asserting the count **and** the rate the sentence quotes.

**M4: pass_two/tests.rs — every fixture sits on contig 0, so the failure that names a record can name the wrong contig and pass**

**Confidence:** High. **Categories:** reliability, extras.

`contig: entry.contig.get()` → `contig: 0` survived all thirteen. `CohortSizeMismatch` names the
record because "pass two walks millions of them and a count alone would say nothing about which" —
and on a real run, where almost nothing is on the first contig, a hardcoded one is wrong almost
every time. **This is C2's Blocker in the same words.**

**Fix:** the mismatch fixture moved to contig 7 and the field asserted; plus
`a_record_wider_than_the_runs_cohort_stops_the_pass`, since a wide record and a narrow one are two
different wrong answers the check exists to prevent.

**M5: pass_two.rs — `ParalogVerdicts` recorded the answer but none of the six constants that produced it**

**Confidence:** High. **Categories:** defaults.

The EM's start, tolerance and cap, the fallback rate, and the histogram's range and bin count. The
curve's binning is private with no accessor and the configuration was dropped after use, so spec
§3.5's run report *could not* print them even if C4 wanted to — and spec §3.1 says these are
prototype-tuned constants nothing has re-measured, which is exactly when an unreportable number
never gets checked.

**Fix:** `config: CalibrationConfig` (it is `Copy`) and `lr_histogram: LrHistogramShape` on the
verdicts, with the histogram built from the three named constants so the recorded shape and the
folded-into shape are the same one —
`the_recorded_shape_is_the_histogram_the_ratios_were_folded_into` holds that.

**M6: pass_two.rs — "in spill order" is the only thing tying a ratio to its record, and it is prose**

**Confidence:** High for the coupling; Medium that it will be broken.
**Categories:** refactor_safety.

`pub ratios: Vec<f64>` is indexable, carries no key, and has no length assertion. A pass three that
read the spill filtered, chunked or sharded — which is how the rest of this run is built — would
give every record its neighbour's verdict, with nothing about the VCF looking wrong.

**Applied in part, deferred in part.** A `debug_assert_eq!` against `entries_written()` and a
sharpened field doc now stand; the reviewer's consuming `RatiosInSpillOrder` cursor is the right
answer and **is deferred to C4**, because spec §3.7's type block declares `ratios: Vec<f64>` and
changing it is a design change rather than a fix. Question 1 above is answered: pass three streams,
so the cursor costs nothing.

**M7: pass_two.rs — the calibration and its differential were both in the run stage, and the stated reason for it was wrong**

**Confidence:** Medium. **Categories:** module_structure.

`calibrate_from_the_ratio_histogram` reads a histogram and a configuration and returns a
calibration; it knows nothing about the spill or the passes. The reason given for placing it beside
its caller was that `src/ng/paralog/calibration.rs` "may not gain a line" — true of that *file* and
not of the *module*, which already holds three files of ng's own. The same mismatch repeated in the
test: every other bit-for-bit differential in the crate lives in `production_parity.rs`, including
the one covering the three pieces this function composes.

**Fix:** new [src/ng/paralog/calibrate.rs](../../../../src/ng/paralog/calibrate.rs), declared to
the copy guard as ng's own; the differential moved to
[production_parity.rs](../../../../src/ng/paralog/production_parity.rs) beside its nine siblings,
and widened from two configurations to three so it also covers M2's case.

**M8: the implementation report and the module doc name one held-constant dimension where there are four**

**Confidence:** High. **Categories:** extras (diff matches stated intent).

The list names GC and argues it harmless. Contig, the inbreeding coefficient (always `0.0`) and the
coverage model (every sample fitted from the same histogram) are constant too and unmentioned. The
GC argument covers two of the three; it does not cover contig, which this pass reads and reports —
and the list is what a reader trusts instead of re-deriving the fixture matrix, so naming one
constant reads as "we checked, there is one".

**Fix:** both paragraphs now name all four and say which are argued and which are tested.

#### Minor

**Mi1** — the `spill.read()` failure path was unreachable from any test; a spill pass one has not
finished writing is now refused by `a_spill_pass_one_has_not_finished_writing_is_refused`.
**Applied.**

**Mi2** — the spill was only tested short, never long. The reader compares its count at end of file
only, so surplus records are decoded, scored and folded *before* the mismatch is seen.
`a_spill_holding_more_records_than_pass_one_counted_is_refused`. **Applied.**

**Mi3** — the single-sample test asserted only finiteness, where spec §4 says four things. Measured:
one sample over two records fits π = 0.50 and reports it converged. The test now says what it does
and does not settle, and D1 is named as where the real answer comes from. **Applied.**

**Mi4** — `observation_of`'s comment promised that a corrupt read-pair row is "counted by the
caller"; this step made `score` the only caller and it counts nothing. Comment corrected.
**Applied.**

**Mi5** — `the_ratios_are_in_spill_order_…` parks its records with the evidence rising, so spill
order and sorted order coincide there and a sorted return would satisfy it.
`the_ratios_follow_the_spill_and_not_the_evidence` parks them out of order. **Applied.**

**Mi6** — `records_scored` named the wrong set: `score` runs on every record and the field counts
only what was folded. Renamed `records_in_the_fit`. **Applied.**

**Mi7** — two private one-line getters left over from the visibility narrowing, adding indirection
inside the type that owns the fields. Inlined. **Applied.**

**Mi8** — the histogram's range is fixed at ±100 while the ratio grows with the cohort: measured on
one duplication-shaped record, 24.2 at one sample, 156.3 at six, **1,663 at 63**. Harmless while it
happens to one class, whose probability is saturated anyway; **not** harmless if a cohort is large
enough to push both classes past the same edge, where they share one bin and the target has nothing
to move. **Made measurable, not fixed**: `ratios_outside_the_histogram` on the verdicts, with
`a_ratio_past_the_end_of_the_histogram_is_counted`. **Raised for Checkpoint C.**

**Mi9** — `a_looser_target_removes_more_records`'s `strict > 0` passes only because four tail
q-values underflow to *exactly* zero (measured `[0.2747, 0.1538, 2.23e-12, 0, 0, 0, 0]`). Commented.
**Applied.**

**Mi10** — a wildcard `other =>` arm in a test match over a `#[non_exhaustive]` enum, which has no
effect inside the defining crate. Replaced with `let … else`. **Applied.**

**Mi11** — the capacity reservation discarded its conversion failure without a word. Commented.
**Applied.**

**Mi12** — `mod.rs`'s "What exists so far" had become six clauses joined by "and". Now a list, one
line a module. **Applied.**

**Mi13** — `#[derive(Clone)]` on `ParalogVerdicts`, whose main field is one `f64` a record.
Dropped. **Applied.**

**Mi14** — `pass_one.rs`'s comment claimed C3 catches a short entry "only on the generic branch";
the check is on the entry's row count, which answers for both variants. Corrected. **Applied**
(C2 code, one comment).

**Mi15** — the report's "three lines" against the commit message's "four lines" against a five-
statement body. Now consistent. **Applied.**

**Mi16** — `ParalogVerdicts` holds no verdict, and `PassTwoError`'s variants are named after the
thing rather than the act. **Not applied**: spec §3.7 names the type, and renaming it is a spec
change. Recorded.

#### Nits

Import ordering in the test module (`super` before `crate`, the reverse of `pass_one.rs`); the
doubled `pub mod` + `pub use` surface giving every type two rustdoc paths, which is the module's
existing convention and wants one decision for the whole module rather than six; the
`records_in_the_fit` binding hoisted above a struct literal it need not be hoisted above; and the
run-stage re-export of `WindowCoverage`, now pure convenience since `ng-window-coverage` merged.

### 6a. The diff's own numbers

Twenty-three claims re-derived by running them. **One wrong, one mechanism incomplete, one
inconsistent, one not reproducible; nineteen correct.**

| claim | verdict | correct value |
|---|---|---|
| mutation 7 is invisible to "the **eleven** behavioural tests" | **WRONG** | **twelve** — 13 tests, exactly one of which bypasses the pass. Re-running the mutation gives `120 passed; 1 failed` |
| the twenty-record fixture gives one cut at every target because the classes "separate completely" | **MECHANISM INCOMPLETE** | separation is real (two ratios, `+156.26` and `−31.15`) but is not the reason. Measured `q(+156.26) = 0` exactly and `q(−31.15) = 0.7499999999999927`: a target of zero was reachable because the duplicated class's q **underflows to exactly zero**, and the cut was constant because the other class sits **above every target tried**. The same fixture with ten duplicated records removes all twenty at 0.5 |
| the fallback is "ng's own **three** lines" (`pass_two.rs`) against "**four** lines" (report, commit message) | **INCONSISTENT** | the body is five statements; all three now say four |
| "12 files fmt-dirty before, 9 after" | **NOT REPRODUCIBLE** | a pre-commit working-tree state git does not hold. Self-consistent, but the step touches four `.rs` files, so one was already clean. Reworded |
| 13 tests; 6595 = 6582 + 13; clippy 3+3+3; fmt 9 files; mutation 4 killed by six tests; mutation 7 by the parity test alone and *because* `flags` reads the curve; the seven-record ratios; 4/5/7 removed at 0, 0.01, 0.5; the parity test's three histograms × five targets; `SampleObservation`'s four fields; the CLI refusal | **CHECKED-CORRECT** | — |

The `PROJECT_STATUS.md` entry repeats no wrong number, notably not the "eleven".

### 7. Out of scope observations

- **Spec §3.2's "a record with no alternative allele is not scored" is implemented in neither
  pass.** Pass one spills such a record — a refused tract written `ALT .` with every sample
  no-called — and pass two scores it on coverage like any tract, giving it a finite ratio that
  enters the fit and can be flagged. Plan step C3 does not ask for it and no deviation records it.
  **For the plan owner; C4 or a spec note.**
- **`SampleObservation::inbreeding_coefficient` is filled and never read by the copied scorer** —
  per-sample `F` reaches the score only through `ParalogScorePrecompute::new`. C1 code, unchanged
  here; an altitude question for whoever revisits the observation type.
- `spill_file.rs`'s `Drop` uses the only `eprintln!` on a non-test path in this module.
- `mod.rs`'s re-export of `WindowCoverage` from a run-stage module is now pure convenience.

### 8. Missing tests to add now

All applied. `an_unscored_record_is_never_flagged_however_loose_the_target`;
`the_warning_counts_the_records_in_the_fit_and_not_the_records_parked`;
`the_ratios_follow_the_spill_and_not_the_evidence`;
`a_record_wider_than_the_runs_cohort_stops_the_pass`;
`a_spill_pass_one_has_not_finished_writing_is_refused`;
`a_spill_holding_more_records_than_pass_one_counted_is_refused`;
`a_target_that_is_not_a_fraction_of_one_is_refused`;
`the_recorded_shape_is_the_histogram_the_ratios_were_folded_into`;
`a_ratio_past_the_end_of_the_histogram_is_counted`;
`the_shipped_fallback_rate_is_the_documented_one`.

### 9. What's good

- **The error chains were checked by rendering them, not by reasoning about them.** All four
  failure paths read well through the crate's own `format_error_chain`: the truncated spill gives
  *"…could not be read from the paralog spill `<path>`: …could not be read at sample 2 of a record:
  …ends inside a record, while reading its mean_depth"*. The outer message adds the stage rather
  than repeating the inner one.
- **Both `SpilledSamples` variants and both row structs are destructured by name with no `_` arm**,
  and `ParalogCalibration`, `ParalogPrior` and `ParalogVerdicts` are built by full struct literals —
  so a field added to any of them fails the build rather than being dropped.
- **Trap 4 is defended on both sides of the copy boundary**: the frozen histogram's own `is_finite`
  guard is pinned from outside by two of this step's tests asserting the folded count.
- **The deviation from production is real, reachable and separated by one test** — installing
  production's screening rule in place of ng's is killed by
  `a_record_whose_samples_have_models_but_no_usable_spread_is_unscored` alone.
- **Memory is flat in the cohort**, as spec §1.1 goal 6 asks: the ratios are 8 bytes a record
  whatever the sample count, the observation buffer is one cohort-length allocation cleared and
  refilled, and each entry is dropped at the end of its iteration.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo test --all-features --lib --bins --tests
./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings 2>&1 \
  | grep -E "^error: " | grep -v "could not compile" | sort | uniq -c
./scripts/dev.sh cargo fmt --check          # counted by unique file
./scripts/dev.sh cargo test --lib --all-features ng::run::paralog_filter
./scripts/dev.sh cargo test --lib --all-features ng::paralog
```
