# ng calling loop — E2f: the ordinary-site prior's seed, from the fit

**Step:** E2f of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md).
**Design authority:** [`spec/population_diversity.md`](../../ng/spec/population_diversity.md) §3.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

The SNP/indel genotype prior takes two numbers about how variable the population is, and **both
were already fitted and neither reached the caller**. The joint fit produces the population's
allele-frequency density — two point masses and a Beta over what segregates — and the expected
heterozygosity read off it; the prior's projection takes an `ExpectedHeterozygosity` and a
`FittedSpectrum` of `2N + 1` allele-count class weights, and `RunParameters::project_seed` has taken
both as arguments since step E2 with nothing supplying them. This is the adapter: `FrequencyDensity`
gains the projection into allele-count classes, `JointFit` gains the wrap, and `run_parameters` gains
the type that owns the class weights for as long as the seed is being projected. **No estimator, and
nothing about what the caller does with either number afterwards.**

## 2. The projection

```text
class k  =  p_segregating · BetaBinomial(k | 2N, a, b)      for every k
class 0     also carries p_invariant
class 2N    also carries p_fixed_alt
```

**Where each point mass lands is the whole of why they are separate parameters.** A position where
the population carries only the reference base is one where no chromosome of any panel carries the
alternative — class 0, at every panel size — and a position fixed for the reference accession's own
private allele is the mirror. Folding either into the Beta would make it a *frequency* a panel could
sample away from, which is exactly what it is not.

**The panel size is the run's fact, not the density's.** The density describes a population and has
no panel in it; which panel it is projected at is the run's sample count. The same density at one
individual and at 63 gives 3 class weights and 127, and a different seed.

**The two bookkeeping numbers `FittedSpectrum` carries**: the count of variable census sites is the
density's segregating share times the positions it was fitted on, and **the regulariser weight is
zero — a measurement rather than a placeholder**. That field says how many sites' worth of
pseudo-counts held the spectrum at the neutral shape; the joint route holds it at nothing, because it
fits the density by maximum likelihood over the census positions with no pseudo-counts anywhere. So
`census_sites_outweigh_regularizer` is true for this route whenever any position segregates — the
honest answer, not a flattering one: the flag exists to catch a *regularised* spectrum whose real
sites lost, and this route has no regulariser to lose to.

## 3. The decision this step had to take: the panel-size floor

**Not ruled — and what this step contributes is a measurement that changes the question.**
`population_diversity.md` §9's third question asks where a panel-size floor belongs, names the
statistic to set it from, and ends *"Confirm before code"*. This step does not set one, and
`FittedFrequencySpectrum::of` does not apply one: it projects at whatever panel it is given.

**The statistic §9 names cannot locate a floor.** It says to sweep how far the fit's answer sits
from the measurement it was fitted to and *"put the floor where it stops falling"*. Swept across
five densities in [`examples/ng_spectrum_panel_floor.rs`](../../../../examples/ng_spectrum_panel_floor.rs),
**it does not fall — it is smallest at the smallest panel**: 1.5 × 10⁻⁹ nats at one individual
against 6.4 × 10⁻⁴ at two hundred, on a tomato-like shape, rising monotonically on all five. That
is structural rather than a property of any cohort. At one diploid individual the spectrum's three
classes are two free numbers after normalisation, against the two-parameter family's two — and
more directly, **the Beta-binomial at two draws *is* a Dirichlet-multinomial**, so the family lands
on the measurement exactly and would do so whatever the measurement said.

**And the sweep is the wrong experiment for the question, which is a sampling one.** It projects one
*exact* density at many panel sizes; nothing in it is about a small panel's estimate being noisy,
which is what a floor is for. The experiment that would settle it is already named, in
[`parameter_prepass_cohort.md`](../../ng/spec/parameter_prepass_cohort.md) §10's third question —
subsample the tomato cohort and watch where the spectrum stops being stable. **That question is
still open and this step did not run it.**

**Who applies a floor, if one is set: not this function.** Whether the run's assembly hands the
prior `Some(spectrum)` or `None` at a small panel is the assembly's decision, and the assembly is
unbuilt. So the ordinary-site ladder's middle rung stays reachable — by that decision, by the
cohort gather below its own designed floor, and by the per-sample histogram route, which supplies a
diversity and no density at all (`population_diversity.md` §3.5).

**One thing worth stating because a first draft of this report got it wrong.** At one individual
the two rungs a floor would switch between agree on the heterozygosity to within 0.1% — but they do
**not** agree on the seed. On the tomato-like density the projection returns a total of 0.223
chromosomes where the neutral rung's is 1.001: the same expected frequency, held with **four and a
half times less conviction**, which is the half of a Dirichlet seed this project has already traced
a genotype-concordance defect to. Nothing here has measured which is better. The draft said a floor
would "replace an exact match with an exact match", and that was true of one of the seed's two
numbers.

## 4. ⚑ What the sweep found that is not this step's to fix

**The projected pair's implied heterozygosity falls away from the density's as the panel grows, and
at real cohort sizes the gap is large.** The projection into classes is exact — the classes carry
the density's own heterozygosity at every panel size, pinned by
`fit::tests::the_classes_carry_the_densitys_heterozygosity_at_every_panel`. What loses it is the
two-parameter fit to those classes, which cannot represent a point mass and trades the ends against
the middle as the class count grows.

| density | at 1 individual | at 63 | at 200 |
|---|---:|---:|---:|
| tomato-like, strong rare-allele pile-up | −0.1% | **−9.9%** | −15.4% |
| human-like, moderate pile-up | −0.1% | **−18.6%** | −25.8% |
| flat over what segregates | +0.1% | −40.9% | −49.4% |
| the unit tests' own lopsided fixture | −0.0% | −15.9% | −24.3% |
| middling frequencies | +0.0% | −53.9% | −61.0% |

**Read that way the sweep argues for a ceiling rather than against a floor**, which is the opposite
of what a first draft of §3 concluded from the same data.

**⚠ Every density above is illustrative, not fitted**, and the caveat belongs with the numbers: no
cohort's fitted `FrequencyDensity` is recorded in this repository, so the grid is a set of Beta
shapes chosen to span this project's two benchmark cohorts' diversities. All the panels are outbred
(`F = 0`); `project_spectrum_seed`'s own documentation records the reference concentration moving
8.6% to 14.0% across `F = 0.6` to `0.9`, so these are outbred-panel figures.

**What it does not do is make this step a regression, and it is worth being precise about why,
because a first draft of this paragraph was not.** **Nothing runs yet** — `RunParameters::assemble`,
`project_seed`, `FittedFrequencySpectrum::of` and `JointFit::fitted_diversity` have no callers
outside tests and examples, and `ng::calling` has no driver. So the comparison is between two
designs, not two runs. Under the design this replaces the prior takes
`ExpectedHeterozygosity::SPECIES_FALLBACK`, a human figure of 1 difference per 1,000 bases; the
tomato-like density's own diversity is 6.06 per 10,000, which the constant overstates by **65%**.
Seeded from that density at 63 individuals the pair implies 9.9% below it. **Closer by a factor of
about seven** — 65.0 against 9.9.

## 5. Two sentences retired

One was flagged by `population_diversity.md` §3.4's ⚠ (`seed_generic.rs:604`); the other two are the
same sentence copied into `SeedRegime`'s own documentation and into a test's doc comment in the same
file, and nothing had flagged either. All three said: *the pre-pass emits the spectrum as absent
below a panel-size floor, so a cohort of five arrives here without one while a single sample arrives
with one.*

**⚑ Whether that sentence is wrong is contested, and this step does not settle it.** §3.4 reads it as
having the two cohort sizes the wrong way round. But `calling_priors.md` §4.1, which owns this
consumer, has a route behind each half — the floor is the *cohort gather's*, and the single sample's
spectrum was to come from the *per-sample histogram* route — so under §4.1 the sizes are not
transposed at all. **§4.1 then contradicts itself two sentences earlier**, saying the single-sample
case *"rests on the per-sample windowed histogram"* and yields `(1, θ)`, which is the *absent*
branch. One of the two spec passages is wrong and it is a spec matter, not this step's.

What all three now say is what is true either way and what this code actually does: **a run arrives
without a spectrum for one of three reasons — the histogram route, the cohort gather below its floor,
or an assembly that chose not to project one — and nothing in this module branches on cohort size.**

## 6. Verification

> The numbers below are the state after the reviews' fixes were applied; §7 says what they found.

| gate | before (E2e, `efd5e9af`) | after |
|---|---|---|
| `cargo test --lib` | 4,842 passed | **4,858** passed, 0 failed, 14 ignored |
| `cargo fmt --all -- --check` | 0 | 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | 0 |
| `cargo doc --no-deps --lib` | 28 unresolved links | **28**, unchanged |

**One release-held check outside a test module** — the panel bound in `allele_count_classes`.
Downgraded to `debug_assert!` and re-run under `cargo test --release --lib --all-features
ng::parameter_estimation::joint::fit`, **two tests fail**, so both of its ends are reached. The
density's two-mass check beside it is `debug_assert!` deliberately: no converged fit can reach it,
the M-step clamps both masses off a normalised branch total, and what it guards is the type's public
fields.

**Nine mutations, nine caught**, each applied alone to the shipped source with the three modules'
tests run against it: the Beta's two shape parameters swapped above one individual; the finished
vector reversed at the seam; both point masses dropped at the seam; both masses put into class 0;
the two masses swapped; the ceiling made exclusive; the variable-site count taken from the
population's segregating share; the same taken as every position fitted; and the hoisted `lgamma`
taken at the wrong argument. **The first three survived the reviews and now fail** — see §7.

## 7. What the three reviews found

Three agents in worktrees from `efd5e9af`, one brief each: arithmetic and numerical behaviour;
tests and mutation; design conformance and claim-checking.

**The arithmetic is clean and was checked against an independent integration.** No Blocker, no
Major: the Beta-binomial formula agrees with a tanh-sinh quadrature of
`∫ Beta(f; a, b) · C(2N, k) f^k (1−f)^(2N−k) df` to within 9.2 × 10⁻¹⁴ over nine `(a, b, N)`
combinations; no `NaN`, infinity or negative weight anywhere from `N = 1` to the ceiling with `a`
and `b` from 0.01 to 1,000; the division by the total is unreachable at zero. Cost is `O(N)` — 0.161
ms at 3,000 individuals, negligible against the 3.3 ms the projection it feeds costs at 200.

**The Blocker was a test blindness, and it is the seventh consecutive step of this plan whose
largest review category was that.** Two mutations survived all 175 tests: **swapping the Beta's two
shape parameters at every panel above one individual** — projecting a rare-allele pile-up as a
common-allele one — and **the seam handing the prior the finished vector reversed**, which on the
fixture density puts 0.995 of the population's weight into the last class, a population read as
nearly fixed for the alternative allele. Four fixture accidents lined up, one per test: the
heterozygosity check weights class `k` by `2k(2N − k)`, which is unchanged under `k → 2N − k`; the
point-mass check ran at `a = b = 2`, where the Beta half of the vector is its own mirror; the only
class-by-class check ran at **one** individual with both masses zero; and every density in the set
had `b ≥ a`, so which end the pile-up sits at was never varied. There is now a class-by-class check
at five individuals with `a ≠ b` and both masses non-zero, against a Beta-binomial computed from the
ratio recurrence — which touches no `lgamma` and shares no code with the implementation — plus an
assertion that the fixture is not its own mirror. The point-mass check moved to `a = 3, b = 1.5`.

**A second Blocker was a test that could not fail, found independently by two of the three
agents.** `the_class_weights_return_the_densitys_own_shares` asserted that the classes total 1 —
which the function's own final divide guarantees — and then that `total − p_invariant −
p_fixed_alt` equals `p_segregating()`, which is that value's definition. Both halves were
identities; measured, the test passed with the body replaced by a flat vector. It now takes each
mass out of **its own end class**, which reads the weights rather than their total.

**`variable_census_sites` was the wrong quantity.** It was the population's segregating share times
the positions fitted, where both the producer's and the consumer's documentation ask for a count of
positions **variable across this panel** — a position can segregate in the population and show one
allele in every chromosome a small panel holds. It is now one minus the two end classes, times the
positions. On the fixture density that is 6.1 in 10,000 at one individual against the 4.0 in 1,000
the old code reported — a 6.6-fold over-report at the small end, shrinking as the panel grows — and
the test that asserted the old value locked the mismatch in.

**The floor ruling was wrong-headed and §3 above is the rewrite.** The first draft ruled *no floor*
on this sweep, and the design review's three objections all stand: the sweep is a fixed-density
experiment where the floor question is a sampling one; the statistic §9 names is degenerate at
`N = 1` for a structural reason, which is a reason to change statistic rather than a measurement of
absence; and the draft's supporting claim — that a floor would "replace an exact match with an exact
match" — was true of the heterozygosity and false of the concentration, which the two rungs
disagree on 4.5-fold. It also erased a second designed producer, the cohort gather, from the
consumer's documentation while that route's spec still stands.

**Nine further wrong claims, all mechanisms or locations, none a number.** Every one of the twenty
measured figures the report and the doc comments quote was reproduced by two agents independently.
The nine: a citation of `population_diversity.md` §3.3 for a claim that section contradicts;
`SeedRegime::NeutralShape` called unreachable when `project_spectrum_seed` reaches it on
`(None, Some(θ))`; "the fitted pair never settles, so a floor would not create one", a non sequitur;
"the family reproduces it whatever it says", false as stated — the Dirichlet-multinomial at two
draws covers only the overdispersed part of the simplex, and the real reason is that it *is* the
Beta-binomial there; "the prior today runs on `SPECIES_FALLBACK`", where nothing runs at all;
a −26% attributed to a fixture whose own figure is −24.3%; `fitted_diversity`'s `None` described as
"a run to refuse" where the ladder below it seeds from the fallback; a citation of `CLAUDE.md` §0,
which has no numbered sections; and a rustdoc link written as a filesystem-relative path.

**And the "two sentences retired" claim was itself half wrong** — see §5, which now records that
whether the sentence was ever wrong is contested between two spec passages, and that a third copy
of it survived in the same file.
