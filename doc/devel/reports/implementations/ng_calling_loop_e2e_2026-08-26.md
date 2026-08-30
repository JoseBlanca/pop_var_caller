# ng calling loop — E2e: the repeat tract's prior seed, from the fit rather than from a construction

**Step:** E2e of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md).
**Design authority:** [`spec/population_diversity.md`](../../ng/spec/population_diversity.md) §4,
§5, which supersedes [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §5's constructed
seed.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed, in one paragraph

A repeat tract's prior belief about which lengths are plausible was built here: mass falling off
geometrically from the cohort's commonest length at the tract, scaled so that the prior's own
implied gene diversity reproduced a cohort-wide repeat gene diversity the pre-pass was supposed to
measure. **Neither of its two run-level inputs had a producer, and the joint repeat fit was already
estimating the object the prior actually wants** — per stratum, a *length spectrum* (how that
stratum's chromosomes are spread over tract lengths, indexed in whole repeat units from the
reference tract length) and a *concentration* (how monomorphic those tracts are). `StratumFits`,
the one lookup that crosses into calling, gathered each stratum's slippage numbers and **dropped
both**. They are now carried across that seam with the rung they came from, `fill_ssr_seed` is
rebuilt on them, and the refusal that fired at every tract at one outbred sample is gone — not
handled, gone: the fitted pair asserts no scaling and has nothing to fail at.

## 2. The three pieces

### 2.1 The seam carries two more values

`StratumFits::over` already reads every fitted stratum's outcome, and `StratumFit` already carries
`length_spectrum` and `concentration`, so the top rung cost no new argument: the gather harvests
them where it harvests the slippage vectors. A stratum furnished from its period's slippage curves
(`DerivedStratum`) carries **no** spectrum by construction — nothing about it was estimated — and
that absence is what the middle rung exists to answer.

**The lookup is keyed by the tract, where its two neighbours are keyed by the candidate.**
`StratumFits::at` and `FrozenParameters::ssr_substitution_rate_at` answer *how does a read of **this
candidate** go wrong*, which is a property of the tract a read was copied from, so they take a
candidate's repeat count. The two of them already sit in one loop in
`inference/repeat_tract_parameters.rs`; the third has no caller yet — E3b is what puts it beside
them, and until then the guard against confusing the two keys is the argument names and this
paragraph. `length_spectrum_at` answers *which lengths can **this tract** be*, which is one question per
locus: the spectrum is indexed by offset from the reference tract length, so passing a candidate's
count re-centres the shape on that candidate. Measured, on the fixture in
`seed_ssr::tests::the_reference_repeat_count_is_the_tracts_and_nothing_here_can_check_it`: with a
candidate's count passed instead of the tract's, the tract's own reference length falls from
**0.595 of the prior's mass to 0.091** and the candidate whose count was passed rises from **0.048
to 0.909**.

### 2.2 The seed is the fitted pair, conditioned onto the locus's candidates

```text
offset:  d_j = repeat count of candidate j − the tract's reference repeat count
shape:   w_j = max(spectrum(d_j), SHAPE_FLOOR)            (1/K on the stated-flat rung)
seed:    α_j = concentration · w_j,  floored at MIN_ALT_CONCENTRATION
```

**Nothing divides by the retained mass, and that is the Dirichlet's own arithmetic rather than a
choice.** A Dirichlet over the fit's `2·span + 1` length classes, conditioned on the tract carrying
one of this locus's candidate lengths, is exactly the Dirichlet over those candidates with each
class's own `α` kept. Two consequences, both asserted:

- normalised, the seed **is** the fitted spectrum restricted to the candidate lengths — which is
  `population_diversity.md` §8's first check, and it holds by construction rather than by
  measurement;
- the total is `concentration × (mass the candidates cover)`, so a locus whose candidates cover a
  tenth of what the stratum spreads over is held with a tenth of the conviction. On the fixture:
  four candidates covering 0.84 of the fitted mass are held with 10.08 chromosomes and one
  candidate covering 0.50 with 6.0, at a fitted concentration of 12.

**A candidate outside the fit's reach takes a floor, not the end class's weight.** The end class
holds the mass the fit put at exactly `±span` repeats; handing it to a candidate five repeats out
would claim the fit measured something its own span is the statement that it declined to measure.
The floor is production's `G0_FLOOR`, kept for its original reason — a masked long heterozygous copy
the candidate set nearly missed has to stay recoverable rather than fall into an absorbing zero —
and here it does a second job: it makes the shape's total strictly positive, so the shares are
always a distribution.

**Two inputs disappeared and one failure with them.** The cohort's modal repeat count at the tract
had no source, because repeat-tract candidate selection is unwritten; the run-wide repeat gene
diversity has no producer at all. `SsrSeedOutcome::DiversityUnreachable` went with them: it existed
only because a constructed shape had to be scaled to a measurement, which is possible only below a
ceiling the shape itself sets — at most 0.625 over the three lengths a single diploid can show,
against the ~0.72 HG002 actually has, so **at one outbred sample it fired at every tract**
(`population_diversity.md` §4.2). `fill_ssr_seed` now returns a `Concentration` rather than an
outcome enum, and there is no refusal left in the module.

### 2.3 The tract ladder, three rungs, each reported

| rung | when | what it is |
|---|---|---|
| `StratumsOwnFit` | the stratum was fitted on its own tracts | its own length spectrum and concentration |
| `PeriodsPooledTracts` | it was not, or it is not in the fit at all | one fit over every tract of its motif period |
| `StatedFlat` | the period has no pool either | a flat shape at a stated concentration |

The rung travels on the value the lookup returns, and `LocusInference` now has a field for it —
but **nothing carries it there yet**, because the driver still refuses every repeat tract. So
`population_diversity.md` §1's third goal, that a call resting on a measurement and one resting on a
stated constant be distinguishable without re-running anything, is **prepared and not met**; E3b
meets it.

**The middle rung is new fitting work, and this is the one place the spec contradicts itself.**
§6 says *"Nothing here is fitted; everything is already computed once per run"*; §4.4 decides the
rung is *"one fit over every tract of that period"*. §6's sentence is true of the seam and of the
ordinary-site side and false of this rung, and §4.4 is the one that decides. So
`fit_period_length_spectra` is a **second, opt-in call** rather than a widening of `fit_strata`: a
run that does not ask pays nothing and still gets an answer at every tract, one rung lower, and says
which. `StratumFits::with_period_length_spectra` is how a run that does ask hands the pools over.

**The bottom rung's concentration is §9's open question 4, taken at its leaning**: the run's own
median fitted concentration where any stratum was fitted, and a stated constant —
`STATED_FLAT_CONCENTRATION`, one chromosome's worth of belief — only where none was. A period's pool
does not enter that median, because it is fitted from the very same tracts the strata's own fits
read.

## 3. Deviations, recorded rather than done quietly

**It edits `genotype_prior/seed_ssr.rs`, which is [`calling_prior.md`](../../ng/impl_plan/calling_prior.md)'s
module**, and rewrites the function that plan's E1 shipped. The plan's own E2e entry says to record
this, and this is the record.

**It deletes `SeedDecayPerRepeat` and its `DomainError` variant** (`src/ng/types.rs`). That type was
the decay of the constructed geometric shape and had no other consumer; with the shape gone it names
nothing. **`RepeatGeneDiversity` stays**, with its documentation corrected: nothing reads it today,
but [`parameter_prepass_cohort.md`](../../ng/spec/parameter_prepass_cohort.md) §3 still specifies the
STR gene diversity as one of the two diversities the pre-pass emits, separately from
`ExpectedHeterozygosity` and precisely so a consumer cannot confuse them.

**It formats `tests/ng_calling_loop_calls_genotypes.rs`**, which was committed unformatted at
`424a0808` and made `cargo fmt --all -- --check` exit non-zero before this step touched anything.
Two hunks, no behaviour.

**It takes `population_diversity.md` §9's question 4 at its leaning, where that question says
"Confirm before code".** The leaning is the run's own median fitted concentration where any stratum
was fitted, and a stated constant only where none was; that is what shipped, and four tests pin its
four properties (the median at an odd count and at an even one, the pools excluded from it, and the
constant reached only by a run that fitted nothing). Taken rather than asked because the standing
instruction for this run is to bank questions that can wait, and this one can: it is one number, it
is named, and it is reached only on the ladder's bottom rung.

**The bottom rung spreads its stated total over the locus's *candidate* lengths where §4.4's table
says the *reachable* lengths** — every length the stutter model can produce from a candidate, a
strictly larger support. The candidate set is what the seed builder is handed; the reachable lengths
are built by the read likelihood (`likelihood::ssr`'s `fill_reachable_lengths`) and are not in the
prior's hands. On that rung the shape is flat either way, so the two differ only in how the total is
divided, and nothing has measured what that costs. Recorded on the constant itself, not only here.

**It replaces `LocusInference::seed_diversity_unreachable` with
`length_spectrum_rung: Option<LengthSpectrumRung>`.** The old field's entire subject was the
refusal this step deleted; it is `None` everywhere today, exactly as the flag was `false`
everywhere, and it is the carrier `population_diversity.md` §8's third check needs. The
"never on the SNP/indel path" check it inherited is kept and re-aimed: a frequency spectrum has a
ladder of its own with different rungs, so a tract rung at an ordinary site is one path's ladder
wired onto the other.

## 4. What it does not do

- **The ordinary-site side**, which is E2f and independent.
- **Wiring the seed into the driver.** `call_locus` still refuses every repeat tract at its front
  door; the branch that would call `fill_ssr_seed` at a tract is E3b's.
- **Reporting the rung in the run's output**, which needs the output stage E3b reaches.
- **Making the tract parameters representable as absent** (`population_diversity.md` §5:
  *absent-or-present as a whole*, and §6: *the only refusal in this document is a repeat tract in a
  run with no repeat-tract parameters at all*). **This is a real gap, not a deferral with a
  reason.** `FrozenParameters` holds a bare `&StratumFits`, so absence has no spelling, and
  `length_spectrum_at` on a gather over no outcomes answers `StatedFlat` at one chromosome — which
  `a_run_that_fitted_no_stratum_states_the_constant` pins. That is §4.4's ladder doing exactly what
  §4.4 says and §5 saying it should have been refused instead, and the two sections are not
  reconciled here.

  **What a later step needs, and the predicate is not the obvious one:** `StratumFits::strata()`
  counts strata carrying *slippage*, and a run whose strata were all furnished from curves has
  `strata() > 0` with no length spectrum anywhere. `strata_with_a_length_spectrum()` is the count
  that answers this question, and it exists for it.

## 5. Verification

> **The numbers below are the state after the review's fixes were applied**; §6 says what the
> reviews found and what moved.

### 5.1 The gates

All in the container, `./scripts/dev.sh`:

| gate | before (424a0808) | after |
|---|---|---|
| `cargo test --lib` | 4,815 passed | **4,842** passed, 0 failed, 14 ignored |
| `cargo test --release --lib ng::calling --all-features` | 754 passed | 751 passed |
| `cargo test --test ng_calling_loop_calls_genotypes` | 10 passed | 10 passed |
| `cargo test --test ng_calling_loop_allocation --features dhat-heap` | 1 passed | 1 passed |
| `cargo fmt --all -- --check` | **non-zero** (see §3) | 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | 0 |
| `cargo doc --no-deps --lib` | 28 unresolved links, exits 101 | **28**, unchanged |

The per-module counts, measured rather than derived: `seed_ssr` 24 → **20** tests (the rewrite
retired the refusal and everything about the decay, and the review added three), `stratum_fits`
6 → **27**, `ssr_fit` 14 → **23**, `calling::tests` +1. Net **+27**, which is the library total's
move.

### 5.2 The release-held checks are all reached

**Ten new `assert!`s live outside a test module**, plus one `panic!` that cannot be downgraded:
two in `pool_a_period` (the read span and the slippage-group count, both properties of the run),
one in `fit_pooled` (an allele span of zero), three in `checked_length_spectrum` (the class count,
the total, the concentration), three in `LengthSpectrum`'s two constructors, and one in
`with_period_length_spectra` (a pool filed under another period's key).

Downgraded to `debug_assert!` in one run and re-run under
`cargo test --release --lib --all-features ng::parameter_estimation::joint`, **13 tests fail** —
every one of the ten checks is reached by at least one test that dies without it. The `panic!` on a
negative length share still fires in release, which is why its test is the one that stays green, and
it is reached by `a_negative_length_share_is_refused`.

### 5.3 Fifteen mutations, fifteen caught

Each applied alone to the shipped source, with the three modules' tests run against it and the
source restored afterwards. The harness was scratch and is not in the tree.

| mutation | caught by |
|---|---|
| the offset's sign flipped | `seed_ssr` |
| the class index off by one | `seed_ssr` |
| a candidate past the reach takes the end class instead of the floor | `seed_ssr` |
| the flat rung writes 1 per candidate instead of `1/K` | `seed_ssr` |
| the seed normalised before scaling | `seed_ssr` |
| the concentration floor dropped | `seed_ssr` |
| the shared export left unnormalised | `seed_ssr` |
| the period's pool consulted before the stratum's own fit | `stratum_fits` |
| the stated concentration taken as the mean | `stratum_fits` |
| an even-count median taking the lower middle | `stratum_fits` |
| the strata's own spectra not harvested at all | `stratum_fits` |
| a pool dragging the stated concentration | `stratum_fits` |
| the pool reading only its first stratum | `ssr_fit` |
| every period pooled together | `ssr_fit` |
| the refusal floor made inclusive | `ssr_fit` |

**Three of these are the fixtures' shape rather than the assertions'**, and are the reason the
fixtures look the way they do. The spectrum's two ends differ (`0.04` against `0.05`) and no two of
its five classes share a weight, so a flipped sign and a shifted index are different numbers rather
than the same one by symmetry. The two strata pooled in `a_periods_pool_reads_every_stratum_of_it`
are drawn from spectra tilted opposite ways, so a pool that read one of them would recover that
one's tilt exactly — drawn from one truth, that mutation survives. And the two periods in
`each_motif_period_is_pooled_apart` are likewise opposed, so pooling them together is one number
twice.

### 5.4 What the middle rung costs, measured

`cargo test --release --lib …::timing_scratch -- --nocapture`, on two strata of 300 tracts each, 8
samples, span 1, in the container: **`fit_strata` 2.68 s, `fit_period_length_spectra` 1.68–1.72 s
on top of it** — three runs, 2.677/2.681/2.690 against 1.668/1.683/1.718. So asking for the middle
rung costs the repeat-tract half of a fit about **60% more**, not double: the pooled climb reads the
same tracts the per-stratum climbs read, and comes out cheaper than their sum because it runs one
climb where they run two. The harness was scratch and is not in the tree; the numbers are quoted
with the fixture that produced them so they can be reproduced.

**This is why the call is opt-in.** A run that does not want to pay still answers at every tract,
from the bottom rung, and says which rung it used.
