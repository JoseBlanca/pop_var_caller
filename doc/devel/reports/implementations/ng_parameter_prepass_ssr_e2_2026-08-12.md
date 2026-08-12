# ng step 4, the STR path — E2: the slippage search, per stratum

*Implementation report, 2026-08-12. Step E2 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), with the review that
followed and the fixes applied — three agents; the reliability agent ran 32 mutations of which 15
survived, and I re-ran eight decisive ones after the fixes. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.2,
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §3, §4.1.*

## What the step is

The three slippage parameters — how often a read slips, which way, how far — searched from four
starting points per stratum, with the genotype frequencies climbed at every candidate.
`fit_slippage` fits one table; `fit_slippage_by_stratum` walks the accumulator;
`starts_must_agree` refuses a fit whose starts landed further apart than
`START_AGREEMENT_LIMIT`.

The mathematics was proven in Milestone D. What E2 adds is the wiring: a stratum's entries become
cells at the key's ploidy, the model is built from the stratum's own repeat count so the low-end
clip happens, the starts are placed around the stratum's own share of moved reads, and everything
the search reports comes back with the answer.

## The two decisions this step had to make

**A search whose inner genotype climbs ran out of passes is reported, not refused.** That surface
is concave, so a climb that ran out did not find a wrong summit — it ran out of time on the right
one. A stratum here has 66 to 91 genotypes against the sibling path's three, and its pass cap was
measured on those three, so making this fatal would end most real runs. It rides on every fit,
the summary will count it, and the walk does not act on it. It is not decoration: a candidate
scored short of its own summit scores low, so a neighbouring candidate can win on that alone.

**A stratum no read of which sits on the whole-repeat grid is not fitted at all.** This was a
review finding and it is the sharpest thing the step turned up. The model scores only whole-repeat
movement, so on such a stratum every genotype and every candidate score alike — and the search
does not then return something obviously wrong. Measured on 500 loci whose every read carried a
one-base deletion inside the tract, it returned a slippage level of **0.5976**, the top of its
range, with the four starts agreeing to 1.00. The one diagnostic guarding this fit could not see
it; the borrowing step's floor could not either, because such a stratum has as many loci as any
other. `fit_slippage` now returns `None` there and the walk leaves the key out, which is what
makes the stratum borrow.

## What else the review changed

**`slipped_reads` was counting reads that did not slip.** It was filled from every read at a
length other than the reference's, which includes the reads that moved by something that is not a
whole number of copies — an ordinary indel inside the tract, an interruption. Those carry nothing
about *which way* a read slipped or *how far*, and those two are exactly what the count gates
(`MIN_SLIPPED_READS_TO_FIT_SHARES`). At the guard shares spec §5 measures, up to 58.5% of the
moved reads, a stratum reporting the floor's 4,000 had put about 1,660 reads behind its fall-off.
The fit now reports three counts — reads that moved by whole copies, reads that moved off the
grid, and reads the model scored — and the level's starting point is the first over the third.

**The outer search's termination was dropped.** `arch` §3 asks for it by name, and what was being
carried instead was the *inner* climb's flag. Both now travel.

**The frequencies had nothing saying which allele lengths index them.** After a merge that is not
recoverable from the key: the monotonicity walk refits with the model of the **lower** repeat
count, so a consumer rebuilding the support from the stratum would index them by a support two
lengths too wide. The support is now emitted beside them.

**The "search did not settle" error now names the ploidy too**, since one stratum is fitted once
per ploidy its loci sat on and the other three fields do not name one fit on a genome with more
than one.

## A defect in shared code, found by a probe rather than by a test

Looking for a fixture whose starts disagree, I hit a `debug_assert` in the shared genotype climb:
*a pass lost ground*. It was not a defect in the climb. The assertion allows a relative slack of
`1e-9 × |score|`, and a table the model explains almost perfectly scores just below zero — one
locus whose single read is where the fitted genotype puts it gives −3.7e-8, where the slack is
3.7e-17, finer than the rounding of the sum that produced it. The measured step was
−3.748510590817489e-8 → −3.748510613021949e-8, a loss of 2.2e-16. The slack gained an absolute
floor of 1e-12, a thousandth of the score's own convergence tolerance.

**It fires only in a debug build, which is what the test suite is**, so on real data this would
have arrived as a panic at the first stratum with a near-perfect fit rather than as a wrong
number. Recorded here because it is a change in `fitting/`, shared with the SNP/indel path.

## ⚠ Five wrong claims of mine

Every figure quoted from the design documents held; all five of these were mine, and two were
wrong *mechanisms*.

1. **"a search that ignored its cells entirely would return whatever its starting points averaged
   to"** — nothing averages the starts (the best start's endpoint is returned), and on a flat
   surface the golden section's ties collapse toward the **top** of the range, so such a search
   returns about 0.6 rather than about 1e-4. The sentence also defeated its own test: all four
   starts sit below the bound it was defending.
2. **"a search that returned its starting points unchanged has not searched"** — an assertion
   that cannot fail. Every axis is line-searched over its whole range at every sweep, so a start's
   own value is overwritten before it is read.
3. **"a genome a quarter covered by runs of homozygosity"** — 29%, in `slippage.rs` and spec §4.2.
4. A documented parameter `starting_level` that the signature does not have.
5. **"66 to 91 genotypes"** is diploid-only — the same two supports of eleven and thirteen allele
   lengths give 1,001 and 1,820 at four copies.

## Tests

Five new beyond the four written before the review, ten in all for this step.

| test | what it pins |
|---|---|
| `a_stratum_where_nothing_slipped_is_fitted_at_almost_no_slippage` | the level runs down to its floor where nothing moved |
| `a_stratum_that_loses_a_repeat_in_one_read_of_ten_is_fitted_at_about_that` | 0.09986 against a truth of exactly 0.1, and gains at 0.005 of the moved reads |
| `the_search_reports_all_four_starts_best_scoring_first` | four starts, distinct on all three axes, placed at three times **this stratum's** share |
| `a_fit_whose_starts_landed_far_apart_is_refused_and_says_whose_it_was` | the limit, its boundary, and that an unsettled climb is not refused |
| `the_walk_fits_every_stratum_under_its_own_key` | one fit per stratum, under the accumulator's keys |
| `a_stratum_whose_starts_land_far_apart_is_refused_by_the_walk` | **the only fixture where the four starts disagree** — spread 67, climbs unsettled, and the walk refuses it |
| `a_short_tract_is_fitted_over_a_clipped_support` | eleven allele lengths and 66 genotypes at four reference copies |
| `a_haploid_stratum_is_fitted_over_single_alleles` | thirteen frequencies where the key says one genome copy |
| `the_reads_that_moved_are_counted_apart_from_the_reads_that_moved_off_the_grid` | 200, 200 and 9,800 over two entries of unequal weight |
| `a_stratum_no_read_of_which_sits_on_the_grid_is_not_fitted` | the 0.5976 case, absent rather than fitted |

**The mutation run's finding was one thing said four ways: on every settled fixture the four
starts reach a bit-identical point and score a bit-identical likelihood.** So "best-scoring
first", "the answer is the best start's", the spread, and the reason there are four starts at all
were each comparing a number with itself — five of the fifteen survivors. The fixture whose starts
genuinely disagree is what gives those assertions content, and it is also the only one that
reaches the walk's error path and the only one where the climbs fail to settle.

Re-run after the fixes, on my own rather than on the agent's word: four identical starts, a fixed
starting level, a constant spread, a constant settled-flag, cells built at a fixed ploidy, a model
built from a fixed repeat count, the agreement check dropped from the walk, and the read counts
losing their per-locus weighting — **all eight now fail, at seven different tests**.

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` clean, and
`cargo test --lib --bins --tests --all-features` in the container: 3,448 → **3,458** lib tests, 0
failed, 11 ignored. The module's tests take 9.6 s, of which 4 s is the one fixture whose climbs
never settle — run at a deliberately coarse search setting for that reason, which changes the
spread from 193 to 67 and not the conclusion.
