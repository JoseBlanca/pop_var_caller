# ng genotype prior — B1: the ported primitive

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step B1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone B.*

> **This report describes the code as submitted for review.** The review found a Blocker and four
> Majors and changed a good deal of it — most consequentially, the fixtures all ran at a reference
> concentration of 1, and the provenance prose named the wrong production function. What landed is
> in [the review](../reviews/ng_calling_prior_b1_2026-08-22.md) and
> [the fixes](ng_calling_prior_b1_fixes_2026-08-22.md).

## 1. Plan

Port production's Dirichlet-multinomial log-prior primitive into
[`genotype_prior/dirichlet_multinomial.rs`](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs),
filling the caller's row instead of returning a `Vec`, and carry across the independent
rising-factorial oracle that checks it.

This is the **random-mating half** of the row: what a genotype would be worth if the sample's
copies were independent draws from the population. Step B2 wraps the inbreeding mixture around it.

Design authority: [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §3.1 (the formula and
why the genotype-independent term is dropped), §8 (no allocation, assertions not `Result`s), §9
(the reuse map and the parity oracle); [`arch/calling_priors.md`](../../ng/arch/calling_priors.md)
§6.

## 2. Assumptions and deviations

**The function is named `fill_random_mating_log_priors`, not `dirichlet_multinomial_log_priors`.**
Two reasons, and the second is the one that decided it. It fills a buffer, and this project names
functions with verbs — the A2 review filed the same complaint against a noun-named method on the
trait. And the name says *which half of the mixture* the values are, which is what step B2 needs to
know when it reads them: the file is named for the distribution, the function for the branch. The
doc comment names production's function as its source in the first line of its provenance section.

**Both oracles are used, not one.** The plan names the rising-factorial computation (spec §12
test 4). This adds a second: bit-for-bit equality against production's own
`dirichlet_multinomial_log_priors`. They check different claims — the first that the *mathematics*
is right, the second that the *transcription* is, which is the claim a port actually owes and the
one the foundations plan's C2 settled by the same route. Neither replaces the other, and §5 shows a
mutation that only one of them catches.

**Two properties from spec §12 landed here rather than at B2** — the 2:1 tripwire (test 1) and the
invariant mass tracking θ (test 2). The plan lists both under B2. They hold at `F = 0`, where the
mixture is the identity on this function's output, so testing them against the primitive tests them
against less machinery; B2 re-checks the ratio through the mixture. Recorded rather than absorbed
silently, because it moves work between two steps of the same milestone.

## 3. Changes made

One file, `src/ng/calling/genotype_prior/dirichlet_multinomial.rs`, **+350 / −5**
(`git diff --numstat`), of which the test module is about two thirds.

`fill_random_mating_log_priors(row: &mut PriorRow<'_>)`:

- fills the caller's per-allele scratch with `lgamma(α_a)`, the baseline every genotype's term
  subtracts — the one thing that changed shape from production, which allocated a `Vec` for it;
- then, per genotype, folds `log C(m; k) + Σ_a [lgamma(α_a + k_a) − lgamma(α_a)]` over the
  alleles the genotype carries a copy of, in the same order as production, and writes one
  `LogProb`.

Every precondition is already held: the structural ones by `PriorRow::new` in release, the value
one by `Concentration::new` in debug.

### The correction the mutation pass forced

The doc comment first said that skipping a zero-count allele was a cost saving, because the term
`lgamma(α_a + 0) − lgamma(α_a)` is zero. **That is true of the arithmetic and false of the code.**
The fold associates as `(acc + lgamma(α_a + k_a)) − lgamma(α_a)`, so a skipped allele's two large,
nearly equal logarithms would enter and leave the accumulator with a rounding in between. Measured:
removing the branch moves a diploid biallelic hom-ref row from `0.6931471805599453` to
`0.6931471805599454`, one unit in the last place, and the parity test fails on it. The doc now says
that, and so does the test that pins it.

## 4. Tests added

Five.

| test | what it pins |
|---|---|
| `every_log_prior_matches_the_rising_factorial_oracle` | **The independent oracle.** Every row matches a rising-factorial computation using no `lgamma` at all, over 348 genotype rows — the ten shapes' genotype counts (1, 4, 1, 3, 6, 21, 20, 15, 7, 9) at four diversities each, from a haploid monomorphic locus to an octoploid one. The count is asserted so a shape dropped from the grid fails here rather than quietly narrowing the check. |
| `the_port_matches_production_bit_for_bit` | **The transcription.** Bit equality with production's own primitive across six shapes and three diversities. Bit rather than tolerance, because the only way to reach it is the same operations in the same order — which is what a port promises and what a re-associated fold quietly breaks. |
| `an_allele_a_genotype_does_not_carry_cannot_move_its_prior` | Moving `α₂` from `1e-4` to `0.75` leaves every genotype without a copy of allele 2 bit-identical, and moves `2/2` — the second assertion is what stops the test passing on a function that ignores the concentration entirely. |
| `the_heterozygote_is_twice_the_homozygous_alternative_at_every_realistic_diversity` | **Spec §12 test 1, the §2.3 tripwire.** The ratio is exactly `2·α_ref : (1 + α_alt)` and within 1% of 2:1 at four diversities. It fails the moment anyone raises `α_ref`: at production's old plug-in value of 10 it is 20:1. |
| `the_homozygous_reference_weight_tracks_the_diversity` | **Spec §12 test 2.** The hom-ref weight is `1 − 3θ/2` to a tolerance that scales with `θ²` — the neglected term — so the assertion tightens as the approximation improves rather than passing on slack, and it falls monotonically as θ rises. |

## 5. Mutation results

**Five run, five killed, none surviving and none inert.**

| mutation | outcome |
|---|---|
| the zero-count branch removed (always compute the term) | killed — **only** by the bit-parity test, at one ulp; the rising-factorial oracle's `1e-12` sails past it |
| the `− lgamma(α_a)` baseline subtraction dropped | killed — three tests |
| the fold re-associated to sum the terms then add the coefficient | killed — only by the bit-parity test |
| the per-allele scratch never filled (stale baselines) | killed — three tests |
| *(the grid count asserted at a guessed 1,264)* | **caught by the test itself before review** — the real figure is 348, and the assertion existed precisely so a made-up number could not survive |

The two mutations that only the parity test catches are the reason both oracles are kept: an
independent re-derivation checks the mathematics and cannot see a floating-point re-association,
which is exactly the class of change a later refactor makes without meaning to.

## 6. Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | 0 | `Finished dev profile … in 6.48s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `test result: ok. 25 passed; 0 failed` |
| `cargo test --release --lib ng::calling::genotype_prior` | 0 | `test result: ok. 22 passed; 0 failed` |
| `cargo test --lib` | 0 | `test result: ok. 4032 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 646.86s` |

## 7. Trade-offs and follow-ups

- **Nothing implements `GenotypePriorModel` yet.** This is a free function; `MarginalizedDirichletPrior`
  arrives at B2 and calls it, then rewrites the row in place with the inbreeding branch.
- **`PROBABILITY_FLOOR` is still unimported.** Nothing here reaches a logarithm of a probability;
  its consumers are B2's Wright test oracle and the comparator at F1.
- **The scratch size question from arch §8 is settled by this file**: the primitive hoists exactly
  one `f64` per allele, which is what the seam carries.
