# ng step 4, the STR path — D1: how far a read slips from the allele it came off

*Implementation report, 2026-08-12. Step D1 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — three agents; the reliability agent ran 21 mutations of which 6
survived and 3 named real defects, and I ran 8 more, four before the review and four after it to
confirm the fixes. Design authority: [`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md)
§3 and §4.1; ported from [`examples/shared/stutter_model.rs`](../../../../examples/shared/stutter_model.rs)
`Slip::p`, as the plan requires.*

## What the step is

`SlippageModel::probability_of_slipping_by(whole_repeats)` — how often a read shows exactly that
many whole motif copies more than the allele it came off, negative for fewer. Three questions asked
in order, one per parameter: did the read slip at all, which way, and how far. The distance is
geometric, truncated at `MAX_SLIP_STEP = 8` copies and renormalised over what is left, so the kernel
is a proper distribution rather than one quietly missing its tail.

**Nothing reads a locus yet**, which is Milestone D's own rule: the mathematics is built and proven
against algebraic identities and the exact-bias harness before a single locus is walked.

## The harness was confirmed green first

The plan makes the exact-bias harness
([`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs)) D's oracle and
requires it green before D starts. Its `gates` section, which is the one D1 and D2 answer to, reads:

| check | result |
|---|---|
| the rule sums to one over the cell space, marginal end buckets, at ±1 / ±2 / ±3 | 1.000000000 — PASS |
| the same rule scored by plugging in the edge offset | 0.948834392 / 0.995446812 / 0.999603090 — FAIL, as designed |
| no bucket charged a negative number of reads, 3,002 cells | 0 violations — PASS |
| a silent kernel puts mass anywhere but the allele's own bucket | 0.000e0 — PASS |
| the control: generate and fit under the same key, reference origin | level +0.000%, up-share +0.0000, fall-off +0.0000, spread across starts 1.000× |

So the 0.9488 figure the plan quotes for the un-rescaled plug-in at ±1 reproduces exactly.

## Recorded deviations from the ported reference

1. **The truncation's renormaliser is the geometric's partial sum, `1 + f + … + f⁷`, not the
   reference's closed form `(1 − f)/(1 − f⁸)`.** Same algebra, different arithmetic: the closed form
   subtracts two nearly equal numbers as the fall-off approaches one, and its relative error runs at
   about `1e-16/(8(1 − f))`. At `1 − f = 1e-9` and a level of one the kernel sums to 0.999999996500,
   which is 3,500 times the tolerance the sums-to-one gate is asserted at — so the gate stops being
   an identity over a band of legal fall-offs. The partial sum's worst deviation over the same rows
   is 2.2e-16. It also removes a special case rather than adding one: at `f = 1` the sum is 8, so the
   uniform limit falls out instead of being branched to.
2. **A fall-off of exactly one returns the uniform distance where the reference returns `NaN`.**
   `SlipStepDecay` accepts 1.0, and `(1 − f)·f^{s−1}/(1 − f⁸)` is `0/0` there. A `NaN` would not fail
   loudly: it reaches the likelihood, and the searches in this crate pick their maximum with
   `total_cmp` ([`coupled_fit.rs:1492`](../../../../src/ng/parameter_estimation/generic/coupled_fit.rs)),
   which ranks `NaN` above every finite score — so a fall-off of one would be *selected* rather than
   skipped, and reported with a `NaN` likelihood beside it.
3. **The offset is range-checked before its absolute value is taken.** The reference takes
   `d.abs()` first, which overflows at `i32::MIN`. Faithful to port, wrong for library code.

## What the review changed

**Blocker — the sums-to-one gate was false over a band of legal inputs (deviation 1).** The reviewer
measured the closed form summing to 0.999999996500 at `1 − f = 1e-9`, and showed the shipped sweep
could not see it because it stepped 0.95 → 1.0 straight over the band. Fixed as above; the sweep now
carries `1 − 1e-8` and `1 − 1e-9`, and reverting to the closed form fails on them.

**Blocker — `i32::MIN` breaks the function two ways, and the test stepped around it.** The truncation
test asserted on `i32::MIN + 1`, which is the tell. At `i32::MIN` itself, `whole_repeats.abs()` panics
in a debug build; in a release build it wraps to `i32::MIN`, which is negative, so a
`copies > MAX_SLIP_STEP` guard never fires — measured, a model of `(0.15, 0.5, 1.0)` charged
**0.009375 to a slip of −2,147,483,648 copies, exactly what it charges to a slip of one**. The check
now runs on the signed offset, before the absolute value, and the test asserts `i32::MIN`.

**Major — the fall-off was pinned on one rung of eight.** The ratio test compared only two copies
against one. The sums-to-one gate is invariant under any *permutation* of the eight distances and the
direction test looks only at the one-copy arms, so exchanging the masses of seven and eight copies
passed all fourteen tests while putting an eight-copy slip at `1/f` times its true probability —
2.5-fold at a fall-off of 0.4. The test now walks every rung; the worst residual on correct code is
5.6e-17 against a 1e-12 tolerance.

**Major — `MAX_SLIP_STEP`'s value was pinned nowhere.** Every assertion spelled the constant
symbolically, so setting it to 2, 4, 5 or 40 left all fourteen tests green. Two is the damaging one:
the kernel then returns exactly zero for a read that slipped three copies. Four is the confusion the
constant's own doc comment warns about in prose — it is then the same width as `OFFSET_HALF_RANGE`,
which is a distance between a different pair of things. `assert_eq!(MAX_SLIP_STEP, 8)` now sits in the
truncation test, matching what `offset_bucket.rs:118` already does for `OFFSET_BUCKETS`.

**I re-ran all four decisive mutants after applying the fixes** rather than taking the agents' word:
the closed form fails the sums-to-one test, `abs()`-before-range-check fails the truncation test,
`MAX_SLIP_STEP = 4` fails it too, and the seven/eight swap fails the ladder test.

## ⚠ Six wrong claims of mine, all about my own arithmetic

Every figure quoted from the design documents was correct; every wrong one was mine.

1. **"the mass beyond eight copies is under 1e-8 of the slipped reads"** — written twice. The mass is
   `f⁸`, which is 3.2e-10 at tomato's fall-off of 0.065 but **5.2e-8** at spec §3's largest measured
   value, 0.123. Now written as the range it is.
2. **"under 1e-10 of all reads at the highest slippage level any stratum reaches"** — at the level of
   0.150 that is **7.9e-9**, eighty times the stated bound. The sentence did not survive its own
   chain either: 1e-8 of slipped reads at a 15% level is 1.5e-9.
3. **"a sweep over realistic values alone would pass with the renormalisation deleted"** — false,
   measured. What the gate sees is `level · f⁸`, not `f⁸`, and the realistic rows fail by 39 times the
   tolerance. The hostile fall-offs earn their place by making the failure unmissable, not by being
   the only thing that catches it.
4. **"leaves this at 0.57 on the worst row"** — 0.57 is the 0.9 row; the worst is 0.95, at 0.34.
5. **"and leaves it green on every row that looks like real data"** — see 3.
6. **A wrong mechanism: a `NaN` "makes every candidate in a search score alike, reading from outside
   as a flat surface."** A `NaN` compares false against everything, so candidates are incomparable
   rather than equal, and what happens depends on the argmax idiom. This crate's is
   `max_by(total_cmp)`, under which `NaN` ranks **above** every finite score — verified by running it
   — so the `NaN` candidate wins rather than being skipped. The corrected sentence is on the test.

## Tests

Six new, fourteen in the module.

| test | what it pins |
|---|---|
| `the_slip_kernel_is_a_distribution_at_every_parameter_setting` | sums to one over 225 parameter triples, and every probability in `[0, 1]` |
| `a_stratum_where_nothing_slips_puts_every_read_on_its_own_allele` | the third algebraic gate, at the kernel |
| `each_further_copy_costs_the_fall_off_whichever_way_the_read_slipped` | every rung of the ladder, both directions |
| `a_read_loses_copies_far_more_often_than_it_gains_them` | the direction is the way round the data — the one thing the three gates cannot see |
| `no_read_slips_further_than_the_truncation_and_the_last_copy_inside_it_is_kept` | the boundary, `i32::MIN`, and the constant's value |
| `a_fall_off_of_one_spreads_the_distance_evenly_instead_of_returning_a_nan` | the `0/0` the closed form has there |

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` and
`cargo test --lib --bins --tests --all-features` in the container. Suite 3,499 → 3,505.
