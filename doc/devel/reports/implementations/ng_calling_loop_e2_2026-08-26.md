# ng calling loop — E2: the pre-pass's outputs, gathered once

**Step:** E2 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — `FrozenParameters`
assembly.
**Design authority:** [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2;
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2, §3.6, §4.3, §4.4;
[`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §2.3.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed

**`src/ng/calling/run_parameters.rs`** — `RunParameters`, which owns what `FrozenParameters`
borrows, with one constructor per run from the pre-pass's outputs and a `view()` that hands the
borrowed form to calling.

**And one field on `FrozenParameters`: the STR substitution rate.** This is what D1's repeat-tract
refusal names. The slippage lookup carries three of a scoring context's four fitted numbers — the
level, the direction split and the fall-off; the fourth, the per-base substitution rate, is fitted
alongside them per `(read group, stratum)` and the pre-pass emits it as a map of its own. It is now
on the borrowed view, with a lookup keyed **by the candidate's repeat count, never the reference
tract's** — the same rule and the same argument naming as `StratumFits::at`, because a read's
chance of mismatching is a property of the tract it was copied from.

## 2. The four rules, each with a failure attached

**A read group with no usable rate gets scale one and says `Defaulted`.** The scale is the fitted
rate over the reads' own mean reported error, so both halves must be there. *No fit* and *a fitted
zero* are different claims and only one is safe to multiply by: a zero rate gives a zero scale,
which charges every read of the library the floor — maximal confidence about every base, from a
number saying the fit found no errors at all. **And the calibration keeps the rate's own warrant**:
a rate borrowed from a sibling read group makes a *borrowed* calibration.

**Contamination is absent or measured, never a fitted zero.** Where no read group identified a
fraction the run is `uncontaminated` and the read likelihood computes its plain formula. Where some
did, every read group needs an entry, and one that identified nothing gets a view of zero fraction
and **zero evidence counts** — because there is no correction to make for what could not be
measured, and `ContaminationView::was_measured` is the only thing that tells that apart from a
library measured and found clean. Both come back near zero.

**The inbreeding coefficients arrive in the run's sample order and are stored as they arrive.** The
pre-pass keys them by sample *name*; the run owns the order, and mapping one onto the other belongs
where the run is assembled. What this checks is that there is at least one.

**The prior's seed is projected once**, by a named step of its own, and not per locus: what varies
per locus is how the seed is spread across that locus's alleles.

## 3. The two axes, and what a gap in either actually costs

`ReadGroupParameters::calibration_of` indexes by `read_group.get() as usize`, so the run's
read-group ids are `0..n` with nothing missing. The pre-pass's maps are keyed by id and could carry
any set.

**⚑ The first draft said a gap "slides every later group's calibration onto its neighbour", and
that is wrong** — the review measured it. The dense vectors are built over `0..count` by *keyed
lookup*, so nothing slides: the missing id's slot takes `defaulted()`, and the vectors come out
shorter than the highest id, so the **highest read group is dropped entirely**. Its symptom is a
panic in `calibration_of` at whichever locus first carries one of that library's reads — a message
about a locus, arriving after the whole pre-pass is finished. The check is worth keeping for that
second reason rather than the first: failing at assembly is the difference between naming the run
and naming a locus.

**And the review found a gap the module did not check at all.** Nothing looked at the
*contamination* map's keys, so an estimate for a read group past the axis was dropped in silence —
a contaminated library left uncorrected. Worse, the walk that builds the views was pinned by a
fixture that could not tell "one view per read group" from "one view per estimate", because every
read group in it had an estimate; a per-*estimate* walk passed the whole suite while charging one
library's 3% to another. Both are fixed: the axis is checked, and the fixture has a fourth read
group with no entry.

The pairing between a fitted rate and its accumulator total is checked in the same pass — the two
come from one pass over one set of reads, so one without the other means they saw different data.
**Only one direction was tested**, and dropping the other half of the key union left the suite
green; the reverse now has its own test, and it is the likelier direction, since the accumulator
runs over every read of every library while a fit can decline to model one.

## 4. Tests

**Seventeen**, 4,716 → 4,733 on the library target. **Seven of them are the review's**, and each
closes something measured to survive the first suite.

| test | what it pins |
|---|---|
| `a_measured_read_group_gets_a_scale_and_keeps_the_rates_warrant` | 0.002 over a reads' mean of 0.004 is a scale of 0.5 — and a borrowed rate makes a borrowed calibration |
| `a_read_group_the_fit_could_not_measure_gets_scale_one_and_says_so` | a zero rate is refused rather than trusted |
| `a_run_where_nothing_was_identified_is_uncontaminated` | absent, not a fitted zero |
| `an_unmeasured_read_group_is_told_apart_from_one_measured_and_clean` | three read groups: 0.03 fitted, 0.00001 measured, and one not measured at all — only the counts separate the last two |
| `the_inbreeding_coefficients_keep_the_order_they_arrive_in` | the run owns the order |
| `the_substitution_rate_is_keyed_by_the_candidates_repeat_count` | 6 repeats and 12 are different strata, at 0.001 and 0.004; 30 is an ordinary absence |
| `a_read_group_whose_accumulator_saw_no_read_gets_scale_one` | the third of the three routes to `Defaulted` |
| `a_fit_noisier_than_the_reads_reported_scales_above_one` | 0.008 over 0.004 is 2.0 — every other fixture puts the fit *below* the reported mean, where a transposed division gives the same 0.5 the test asserts |
| `the_substitution_rate_is_keyed_by_the_runs_ploidy` | a haploid run finds its own strata: hard-coding diploid in the key survived every other fixture |
| `the_substitution_rate_is_keyed_by_the_read_group` | two libraries at one stratum, at 0.001 and 0.006 |
| `the_read_group_count_is_the_calibration_axis` | not the contamination vector's, which is empty on an uncontaminated run |
| six `#[should_panic]` | a gap in the read-group ids, a rate without its accumulator total **and the reverse**, a contamination estimate off the axis, a run with no samples, a run with no read groups |

**Two fixtures were rebuilt because they could not fail.** The contamination fixture gave every
read group an estimate, so "one view per read group" and "one view per estimate" were the same
number — a per-estimate walk passed. And the inbreeding fixture was `[0.0, 0.9, 0.0]`, a
palindrome: reversing the stored order left a bit-identical vector, and reversal is the wrong
implementation a map-keyed source most easily produces. It is now `[0.1, 0.9, 0.5]`, with every
slot asserted.

**One measured detail worth recording**: the calibration scale comes back `0.500000000282`, not
0.5. The accumulator sums each read's log error in **fixed point**, precisely so that merging
shards in different orders gives the same denominator, and the price of that determinism is a
quantised mean.

## 5. Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4733 passed; 0 failed; 14 ignored`. Before E2: **4,716**.
- `cargo test --release --lib ng::calling --all-features` — `687 passed; 0 failed; 3 ignored`.
- `cargo test --test ng_calling_loop_allocation --features dhat-heap` — `1 passed`.
- **`cargo doc --no-deps --lib`** — the review ran it and this diff had added a broken intra-doc
  link, which the crate denies. Fixed; the module now contributes none of the build's remaining
  (pre-existing) 28.
- **The release-held checks: E2 adds five.** Downgraded all five to `debug_assert` together and
  re-ran under `--release`: **6 failed**, every one reached.

## 6. What the review found

**One agent, carrying reliability, step 8a and a lighter naming/errors/smells pass** — the diff is
one type, one constructor, one lookup. **19 mutations run, 4 survived, 3 changed no behaviour**;
**1 Blocker, 5 Majors, 5 Minors**, and **6 of 24 claims wrong**. All applied.

**The Blocker was a mechanism stated in six places and true in none** (§3). **Three of the five
Majors were tests that could not fail**: the contamination walk, the reverse pairing direction, and
the ploidy in the substitution lookup's key — each mutation left the suite green.

**A fourth was a field this module fabricates.** An unmeasured read group's `ContaminationView`
carries a `source` of `TheWholeSamplesReads`, whose own documentation defines it as *"fitted from
every read of the sample and copied onto this read group"* — a positive claim about a number that
was never fitted. The fraction and the counts honour *absent is not a fitted zero*; `source` cannot,
because the type has no variant that says *not measured*. It is now documented at the constant and
pinned by a test, and the type's owner is the one who could make it unrepresentable.

**The fifth was `cargo doc`**, which the branch's gate does not run and which the crate denies
broken links on.

**One correction worth keeping**: the `.zip()` that pairs a rate with its accumulator total has an
unreachable empty arm — the pairing check upstream guarantees both — and its comment claimed that
arm was the `Defaulted` route. An `expect` planted there never fired across the whole suite. What
`from_fitted_rate` actually refuses is a *zero* rate.

## 7. What this step owes, and what it does not unblock

**⚑ A repeat tract is still refused by the driver, and E2 was only half of what it needed.** The
rate D1's message names is now reachable through `FrozenParameters`. What is still unbuilt is the
**assembly** of a tract's scoring contexts: one `SsrScoringContext` per `(read group, candidate)`,
each holding a `StutterModel` built from the slippage fit, that rate, the unreachable mass and a
warrant — plus the reachable-length support and the outlier weight the row also takes. That is a
per-locus build with a borrowing shape of its own (the contexts borrow the stutter models, so the
two cannot live in one struct), and it belongs with the table build in D1's
`build_genotype_likelihood_table`.

**Recommendation: give it its own step before E3**, rather than folding it into the integration
fixture. E3's job is to show that genotypes come out of real evidence; a step that first has to
invent a borrowing shape for the tract's parameters is a second job, and the plan's own rule is one
step, one commit.

**⚑ The substitution lookup collapses three absences into one `None`.** Its sibling
`StratumFits::at` returns a `Result` that separates *unknown read group*, *no such stratum* and
*group not in the fit* — and that last one, its own comment says, "is not a quiet library": it means
two things were assembled from different runs. The new lookup says only `None`, which its doc calls
ordinary. The two fill the same scoring context, so the assembled context will carry the weaker
answer. **Recommendation: mirror the error type when the tract's context assembly is built**, which
is the step that first has both in hand.

**⚑ And E2's constructor takes the pre-pass's outputs as arguments rather than finding them.** The
name-to-run-order mapping for the inbreeding coefficients, and the choice of which fitting route's
outputs feed the calibration, belong where a run is assembled — which the plan puts out of scope
("wiring `call_locus` into the merge's builder for real runs"). This step's contract is *one
constructor per run from the pre-pass's outputs*, and that is what it is.
