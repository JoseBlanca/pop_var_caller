# ng read likelihoods — B2: the closed form, and a claim the specification could not keep

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step B2 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone B, on
top of `f9e38349`.*

## 1. What it is

`genotype_log_likelihood_row` — **the first thing in this plan that computes a likelihood.** One
call, one sample, one number per candidate genotype:

```text
log Lg(g)  =   Σ         n_o · log( k_{a(o)} / P )              ← the reads it explains
            o : k_{a(o)} > 0

            +  Σ         [ q_sum_o + n_o · ( log scale − log m ) ]
            o : k_{a(o)} = 0                                     ← the reads it calls errors

            +  q_sum_other                                       ← the pooled leftover
```

A read either shows something the genotype can produce, in which case it is charged only for which
copy it came from, or it does not, in which case it is charged for being wrong. That is the whole
of it.

## 2. The decision B1's review left to this step

**The error-spread table now stores `log m` rather than `m`**, and is renamed accordingly.

The argument is not the 2.5× the lookup measured. It is that **spec §3.3's closed form takes no
logarithm at all inside its loop** — every logarithm it needs is a property of something that does
not vary per term: `log(k/P)` of the copy count and the ploidy, `log scale` of the read group,
`log m` of the `(allele, genotype)` pair. All three are computed before the observation walk, and
what is left inside is a multiply and an add. A table holding `m` would have put an `ln` back into
the one place the formula was shaped to keep clear of.

It also settles the naming question the review raised: **nobody divides by 1.0986**, so *divisor*
stops being the word. `fill_log_error_spreads` fills, `LogErrorSpreadTable` reads.

## 3. The production differential — this milestone's whole claim

ng's row and production's `standard_log_likelihood` are the same closed form with two recorded
changes: production carries a multinomial coefficient that ng drops (spec §3.4), and ng divides a
wrong read's error mass by three where production divides by nothing (spec §3.5).

**Add the coefficient back, take the spread out, and the two agree — for all six genotypes of a
triallelic diploid, to better than 10⁻⁹.** Every difference attributed; none unexplained.

**Three differences, not two, and the third is what review found.** Production has no calibration
at all, so ng's `n·log scale` on the error side has to come back out too. At a defaulted
calibration that term is exactly zero — which is how **a mutation deleting `log scale` outright
survived every test in this file.** Measured cost of that mutation on a four-allele fixture at
scale 0.37: **372.84 nats, about 1,620 Phred, with the whole suite green.** The differential now
runs at a scale of 2.5 and reconciles the term explicitly, and the hand-computed row is checked
twice, once at each.

**Neither shortcut that would let the comparison agree with itself is taken.** The spread comes out
of the table rather than a literal `n · ln 3` — the gate B1's review asked this step to pass, since
with the literal it would go on passing with the whole error-spread step deleted. And the
coefficient's `ln(n!)` is written locally rather than borrowed from production's lookup table,
which would cancel a wrong entry on both sides.

This needed one change to frozen production — `standard_log_likelihood` widened from a private
`fn` to `pub(crate)`. That is the single exception the freeze allows: visibility so a parity test
can name the oracle, and nothing else. Nothing shipped in `src/ng/` calls it.

**B1 is now load-bearing, and at the library level rather than only in tests.** The row's own
signature names `LogErrorSpreadTable`, so deleting B1's contents is `error[E0425]` before any test
runs, and the re-export makes it `error[E0432]` as well. That is the gate the plan set at this
step. What has *not* happened is the row itself becoming load-bearing — nothing outside its own
tests calls it, and nothing will until the calling loop does.

## 4. A specification claim that was false, found by refusing to write the easy test

**Spec §2.3 said the aggregation identity is bitwise.** The identity is that the likelihood of a
list of individual reads and the likelihood of the merge's fold of those same reads are the same
number — the property that makes it safe for the merge to throw the individual reads away.

The first fixture written for it **did** agree to the last bit. That is what the plan asked for and
it would have been easy to stop there.

**It was luck.** A 144-combination sweep — six quality sets from two to seven reads and Phred 45
down to 3, three read-group scales, ploidy 2 and 4 — finds the two forms up to **two units in the
last place apart**. They are not the same sequence of additions: five reads accumulate
`q_r + log scale − log m` one at a time, where the fold accumulates
`Σ q_r + n·log scale − n·log m` once, and floating-point addition is not associative.

**Spec §2.3 and §12's ninth test are corrected**, and the correction keeps two claims apart that
the old sentence ran together:

- **The model requirement is untouched and is the one that matters.** No term may be a non-linear
  function of a per-read quality, because `q_sum` recovers only the geometric mean of the reads'
  error probabilities. §3.3's formula satisfies that *exactly*, and that is what the test pins.
- **The arithmetic claim was too strong.** Summation order is not exact.

**The single-fixture test now asserts the measured bound too**, with a comment saying its own
exactness is a property of its fixture rather than of the formula. A test that passes for a reason
its fixture does not establish is a test that breaks on a change costing nothing.

### The review broke this correction in four places, and it was right to

**A unit-in-the-last-place count was the wrong bound, because it grows with depth.** The first
sweep stopped at seven reads. Widened to the depths this caller commits to, the same axes give six
ulps at twenty reads and **past a hundred at 300** — it is repeated summation, so it grows with how
much is summed. The bound is now **relative, 2 × 10⁻¹⁴ over 864 combinations** with read counts
from 2 to 300, four quality profiles including all-equal at Phred 93, three read-group scales and
two ploidies. At the depths where it is largest that is about 10⁻¹⁰ Phred.

**"About 10⁻¹⁵ nats" understated it fourteenfold** — measured 1.4 × 10⁻¹⁴ at the worst point, and
nats are the wrong unit anyway, because the row's magnitude scales with depth and so does the
absolute error.

**"No formula of this shape could make it exact" was too strong, and the true statement is
sharper.** A shape *does* exist for the single-observation case: accumulate the read *counts* —
integers, so their sums are exact — per `(read group, spread class)` and per copy count, and
multiply each logarithm in once at the end. A reviewer implemented it and measured **zero
disagreement over four million comparisons**. It is refused for two stated reasons rather than
because it is impossible: it buys nothing past that corner, since with two or more observations the
per-read form sums flat where the fold sums a tree and the row cannot recover a grouping the fold
has already consumed — the merge's own `q_sum` is itself a tree sum — and it costs
`genotypes × read groups × 2` accumulators held across the observation walk, which spec §8's
no-allocation contract does not provide.

**And §12's *eighth* test carried the same false claim, eight sections from the one that was
corrected.** It said permuting the observations changes no genotype's likelihood "by a single
bit"; permuting the observations *is* changing the summation order, and the row's own fixture is
one ulp apart. Corrected, along with the test's name and a comment that claimed its tolerance was
"stronger" than the aggregation bound when it was about twice as loose.

**The survivor the review found has no fix here, and that is the point of the disagreement
check.** Rewriting the row's accumulation as the exact shape above leaves all 32 tests green — so
nothing pinned the summation shape that §2.3 now makes a claim about. The sweep now asserts that
**some** comparison disagrees, which the exact shape would fail.

## 5. What the tests pin

| test | the defect it fails on |
|---|---|
| `a_biallelic_diploid_row_is_what_the_formula_says_term_by_term` | any term of the closed form, at three genotypes hand-computed independently — **and again at a read-group scale that is not one**, because `log scale` is exactly zero at a defaulted calibration and a mutation dropping it slipped past the first half |
| `ng_and_production_agree_once_the_two_recorded_changes_are_undone` | the closed form diverging from production by anything other than the two recorded changes |
| `how_far_the_bitwise_aggregation_claim_reaches` | the aggregation identity widening past what was measured — and it is the test that found the specification wrong |
| `pooling_an_observations_reads_does_not_change_the_answer` | a term that is a non-linear function of a per-read quality, which is the contract §2.3 actually asks for |
| `permuting_the_observations_moves_a_row_only_by_what_the_order_costs` | a row that sorted, bucketed or re-grouped internally |
| `the_pooled_leftover_shifts_every_genotype_by_the_same_amount` | the leftover added per observation rather than once, which moves no genotype and every QUAL |
| `a_sample_with_no_evidence_scores_every_genotype_at_zero` | a branch where an empty sum should do the work |
| `the_spread_is_not_zero_in_the_production_differentials_fixture` | the differential reconciling a spread that was zero anyway |

## 5a. Four more corrections the reviews earned

**The architecture document had not been updated at all** — it still declared the names B1's rename
retired, still posed the `log m` question as open, still showed a row signature with a contamination
slice and a scratch the row does not take, and still said the aggregation identity is bitwise.
Corrected, including its own version of the row contract: *reproducible* at any thread count with
the observations summed in the merge's order, rather than "bit-identical at any observation order",
which is the same overreach in a third place.

**Spec §3.3's bound on the dropped `log(1 − ε)` ran the wrong way.** The omission is `−n·ln(1 − ε)`,
which is `n·ε` **and a little more**, not `n·ε` at most; and the 0.75 nats quoted at 300 reads and a
poor library is the linearisation — the exact value is **0.77 nats, 3.3 Phred**, 2.6% larger.

**`MAX_PLOIDY_COPIES`'s doc claimed an assertion that did not exist.** `Ploidy::try_new` rejects
only zero, so seventeen copies is constructible and used to panic with `index out of bounds` — the
one caller bug in this file whose panic said nothing about what went wrong. It asserts now.

**And the row's own per-term argument was not applied to its own lookup.** It reached the spread
through a checked accessor per `(observation, genotype)` — a bounds check, a `checked_mul` and a
formatting closure inside a loop the doc describes as a multiply and an add. For a fixed
observation the allele is fixed, so the loop now walks that allele's *column*: checked once,
strided after.

## 5b. The gate this step set was still open one level down

**The mutation review's headline: replacing the row's table lookup with the constant
`LOG_ERROR_SPREAD` passed every test.** That is exactly the shape B1's review made a gate for the
production differential — take the spread from a literal and the test agrees with itself — and it
was still open in the row, for a reason no reviewer had to guess at: **every row fixture used
single-base alleles**, where every allele pair is one substitution apart, so the `m = 1` arm never
appeared on the error side of any row assertion. A locus carrying a deletion now does, and against
its deletion homozygote *neither* wrong read is one substitution away, so the charge is the two
`q_sum`s alone.

Three more terms turned out to be pinned by nothing:

- **The ploidy in the copy share.** `ln(k/P)` at `k = P` is zero at every ploidy, and a diploid
  heterozygote's `ln(1/2)` was all any fixture constrained — so a hardcoded ploidy of two survived,
  returning **positive log-probabilities** at a tetraploid and reordering its genotypes. The sweep
  could not see it: it compares a tetraploid row against another tetraploid row.
- **`LOG_ERROR_SPREAD` itself**, pinned only to four decimals by the size test, and tied to nothing
  — `ERROR_SPREAD_BASES` is read by no code at all. The two are now asserted equal.
- **Three assertions with no test**: the row's allele range, the row's spread-table stride, and
  `LogErrorSpreadTable::over`'s own length check, which is the one the type exists for and which
  only the *fill*'s equivalent had been tested.

And **two doc comments claimed guarantees they did not provide**: the calibrated-scale helper said
it would catch a dropped `log scale` (every fixture using it compared the row against itself), and
the spread-is-not-zero test said the differential would fail without it (both sides read the same
table, so a uniformly zero table cancels).

## 6. Validation

In the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,268 passed / 0 failed / 14 ignored**; the likelihood module holds 81,
  of which 18 are this step's.
- `cargo doc --no-deps` — 23 unresolved links, the same 23 that are on `main`.
