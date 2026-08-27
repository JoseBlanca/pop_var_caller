# The ordinary-site seed: the diversity is pinned, and the shape's ramp does not point where it was thought to

**Branch** `ng-seed-shrinkage`, worktree `../pop_var_caller-seed-shrinkage`, cut from
`ng-calling-loop` at `1474c5cc`.
**Spec** [`ordinary_site_seed.md`](../../ng/spec/ordinary_site_seed.md) — §3 and §4, built here.
**Date** 2026-08-26.

---

## 1. What this is, in one paragraph

The SNP/indel genotype prior starts from two numbers: how many chromosomes' worth of prior belief
sit on the reference allele and how many on the alternatives. Until now both came from one place —
a two-parameter fit to the panel's allele-count classes. **That fit loses the population's
diversity, and loses more the larger the panel**: 9.9% at 63 individuals on a tomato-like shape,
18.6% on a human-like one (spec §1.2). The change here separates the pair into the two things it
really is — an expected allele frequency and a total conviction — and takes each from where it is
measured well. The total is now solved from the run's own fitted heterozygosity, so the seed
reproduces that measurement exactly at every panel size and every shape. The frequency is blended
between the neutral shape and the panel's own.

**And the measurement that was supposed to set the blend's one constant says the blend points the
wrong way.** §4 assumed a panel's own shape is a poor guess when the panel is small and a good one
when it is large, so the weight should rise with the panel. Measured, the opposite holds: **the
panel's own shape is exact at one individual and degrades monotonically from there.** §5 below
gives the evidence and the mechanism. The ramp is built as specified and its constant measured;
**the measurement's own answer is zero — the ramp deleted — and what ships is a quarter of an
individual, a hedge whose price §5.2 states and whose fate §8 puts to the owner.**

---

## 2. Plan, as carried out

1. **§3 first**, because it is what makes both ends of §4's ramp imply the same diversity, so the
   ramp then interpolates shape alone.
2. Split what the search returns into the part the seed keeps and the part it replaces
   (`FittedShape`), so the discarded number stays reachable for the two programs that measure what
   discarding it costs.
3. Build the ramp with a placeholder constant, and the tests that do not depend on its value.
4. **Fit the constant** on drawn cohorts across both axes the caller commits to — cohort size and
   read depth — in a new example.
5. Set the constant, retune the tests that depend on it, and re-run the sweep so the example's
   output is the shipped one.

---

## 3. Assumptions, where the spec left something open

Each of these was a choice; each is stated so it can be overruled.

- **A spectrum that arrives with no diversity beside it is refused rather than projected.** After
  §3 the pair's total comes from the measurement, so a run with no measurement has nothing to pin a
  shape to; the run falls to the species-range guess and reports `FallbackDiversity`, discarding
  the shape. In practice the two arrive together — the joint route reads its heterozygosity off the
  same density it projects — so this is the degenerate-fit path. **Before this change that case
  produced a fitted seed**, because the diversity was not read at all.
- **A diversity of exactly zero short-circuits before the search runs.** It costs nothing to fit a
  shape that will scale nothing, so the 399 predictions are skipped and the seed is
  `(1, MIN_ALT_CONCENTRATION)` with its own regime.
- **The zero-diversity regime carries no shape weight.** With no diversity there is nothing for a
  shape to scale, so a run that had a spectrum and one that did not receive the same seed and the
  same record.
- **`SpectrumMatch` is carried on the unreachable-diversity regime as well as on the fitted one.**
  The spec says the two things already carried stay carried; on this path the pair the search
  produced is discarded, and how good that pair was is the difference between a contradictory
  measurement and a search that never got near one.
- **No floor is applied to the solved pair.** An earlier draft floored the alternative
  concentration at `MIN_ALT_CONCENTRATION`, and that floor binds below a measured diversity of
  about `2e-12` — where it breaks the one thing §3 guarantees. The floor belongs to the
  zero-diversity case, which is taken separately, and to the per-locus expansion, which already
  applies it.
- **The comparison against the shape's ceiling is `≥`, not `>`.** At exactly the ceiling the solved
  total is infinite and `SpectrumSeed` refuses a non-finite concentration, so a strictly-greater
  comparison would turn a reported fall-back into a panic at the run's assembly.

---

## 4. Changes made

### `src/ng/calling/genotype_prior/seed_generic.rs`

- **`FittedShape` and `fit_spectrum_shape`** — what the search returns, split into the expected
  frequency the seed keeps and the pair the seed replaces. Public, because the sweep that set the
  blend's constant has to measure the shipped search rather than a copy of it.
- **`HALF_WEIGHT_PANEL_SIZE` and `panel_shape_weight`** — the ramp, `w = N / (N + N₀)`.
- **`blend_expected_frequency`** — the geometric blend, in log space because the two ends can be
  orders of magnitude apart.
- **`total_for_diversity` and `PinnedTotal`** — §3's identity, solved, with the one case that has
  no answer named rather than clamped.
- **`project_spectrum_seed`** — rewritten around the two: blend the shape, solve the total, and
  branch on the three ways the solve can fail.

### `src/ng/calling/genotype_prior/mod.rs`

`SeedRegime` gains the weight on its fitted variant and two variants for §3.1's refusals:
`DiversityUnreachable` and `ZeroDiversity`. The two runs a reader must be able to tell apart — one
that fell to the neutral rung because no spectrum arrived, and one that fell there because no total
could reach its measured diversity — are now different values.

### `src/ng/calling/run_parameters.rs`

Its three projection tests handed the projection a hard-coded `1e-3` while the seed did not read
it. They now pass the density's own heterozygosity. **This was a fixture accident of exactly the
kind this plan keeps finding**: the argument was there, and nothing read it, so nothing checked it.

### `examples/ng_seed_shape_weight_sweep.rs` — new

The measurement that set the constant. §5 is its result.

### Two existing examples

`ng_spectrum_panel_floor.rs` now reports the *search's* pair rather than the seed's, because that
is what it was measuring — the loss the pin exists to remove. `ng_inbreeding_sensitivity.rs`
reports **both**, in four columns, and §5.6 is what that showed. Each carries a note saying the
shipped seed no longer takes its total from the search.

### Five documents, because the sentences this change retires were standing in five places

The first sweep found three and the design review found two more, which is the same failure this
plan keeps repeating — the sixth step running where a retired sentence outlived the behaviour it
described.

| document | what it said | what it says now |
|---|---|---|
| `seed_generic.rs`'s own module header and `project_spectrum_seed` | the seed is *"read off the pre-pass's fitted frequency spectrum"*, and *"a spectrum makes the diversity moot: it carries its own scale"* | rewritten; the diversity is read on every path |
| `population_diversity.md` §3.4 | the top rung takes *"shape and scale both from it; the diversity is not read"* | the scale is the measurement's, the shape is blended, and the two rungs are the ends of one ramp |
| `population_diversity.md` §9, question 3 | where to put a panel-size floor — *"Confirm before code"* | **answered**: no floor, because there is no switch left to place one at, and the statistic that question named is smallest where a floor would fire |
| `calling_priors.md` §4.1 | the seed is the pair whose predicted spectrum best matches | a superseded note at the head; the search now supplies the shape alone |
| `impl_plan/calling_loop.md`, step E2f | *"⚑ One decision this step must take rather than inherit: the panel-size floor"* | retired, with a pointer — it was a live instruction to take a decision that no longer exists |
| `arch/calling_priors.md` §2.3 | `SeedRegime` drawn with three variants and a `data_dominated` field | a ⚠ naming the five variants that shipped and the weight the fitted one carries |

And `run_parameters.rs`'s own `FittedFrequencySpectrum::of` carried four claims about behaviour this
change alters, one of them — *"a prior four and a half times more easily moved by the reads"* —
false of the shipped seed. It is now a fact about the search, with the shipped pairs beside it.

---

## 5. The measurement, and what it says

`examples/ng_seed_shape_weight_sweep.rs`. Drawn cohorts at known parameters: 1 to 63 diploid
individuals, at 3, 8 and 20 reads a sample, two population shapes, 3,000 positions a cohort, six
drawn cohorts a cell for the fit and four held out. **The same drawn positions are refitted at
every cohort size**, which is what stops the answer moving for a reason that has nothing to do with
panel size.

### 5.1 The control that settles it: the projection alone

The density is handed straight to `allele_count_classes` and the shipped search is run over the
result. **No cohort is drawn and nothing is estimated**, so the only thing that moves between rows
is how many allele-count classes the two-parameter family is asked to fit at once. The number is
the expected alternative-allele frequency the search reads back, over the density's own.

| individuals | 1 | 2 | 5 | 10 | 25 | 63 | 200 |
|---|---:|---:|---:|---:|---:|---:|---:|
| tomato-like, Beta(0.20, 1.00) | 0.999× | 1.021× | 1.063× | 1.097× | 1.141× | 1.177× | 1.217× |
| human-like, Beta(0.35, 1.20) | 1.000× | 1.019× | 1.054× | 1.079× | 1.114× | 1.141× | 1.164× |
| flat, Beta(1.00, 1.00) | 1.000× | 0.979× | 0.942× | 0.912× | 0.882× | 0.862× | 0.843× |
| lopsided, Beta(0.50, 2.00) | 1.000× | 1.031× | 1.091× | 1.138× | 1.184× | 1.217× | 1.242× |
| middling, Beta(4.00, 4.00) | 1.000× | 0.968× | 0.905× | 0.872× | 0.845× | 0.831× | 0.818× |

**Exact at one individual on all five shapes, and monotonically worse from there.** The direction
of the error depends on the shape; the fact that it grows with the panel does not.

**The mechanism is one individual's arithmetic.** A panel of one has three allele-count classes,
which after normalisation are two free numbers, against the family's two parameters — so the fit
reproduces the panel's first two moments exactly, point masses included, and those are the
population's own. At 63 individuals the same two parameters are fitted over 127 classes and can no
longer absorb the mass piled at *invariant*, so they compromise. **This is spec §1.2's mechanism,
showing up on the ratio of the pair rather than on its total.**

### 5.2 What that does to §4.1's criterion, and what the constant was fitted to instead

§4.1 sets the constant at the panel size where the panel's own shape and the neutral shape are
equally trustworthy, read off where their errors cross. **They cross the other way round, on both
shapes and at all three depths**: the panel's own shape is the nearer one at the smallest panels
and loses ground as the panel grows. At a single genome it lands within 0.4% to 5.7% of the
frequency that genome was drawn with, in all six cells, against 7% to 21% for the neutral shape.

So the constant was fitted the other way it can be: **the value that puts the blended shape
nearest the truth, averaged over every panel size, depth and population.**

| half-weight panel size | 0 | 0.25 | 1 | 12 | 200 |
|---|---:|---:|---:|---:|---:|
| drawn cohorts, held out | **0.1486** | 0.1526 | 0.1616 | 0.1709 | 0.1604 |
| drawn cohorts, fitted on | **0.1543** | 0.1573 | 0.1653 | 0.1859 | 0.1854 |
| no drawing, 4 in 1,000 segregating | **0.0994** | 0.1255 | 0.1869 | 0.4046 | 0.6599 |

**All three arms put the minimum at zero** and rise monotonically over the first decade. A
half-weight panel size of zero is a weight of exactly one at every panel — the blend inert, and
§4's ramp a no-op.

**So `HALF_WEIGHT_PANEL_SIZE = 0.25` is a hedge and not a fit, and the report says so rather than
dressing it up.** It costs 2.7% against the minimum on the held-out cohorts, 1.9% on the ones the
constant was fitted on and **26% on the arm with a real cohort's share of variable positions**.
What it buys is a fifth of the shape kept on the neutral side at a single genome — a guard against
a pathological fit at a small panel, which is `ordinary_site_seed.md` §1.1's concern and which this
sweep **cannot see**, because every density it draws from is well behaved. **⚑ Whether to pay that
is the owner's call**, and §8 states the alternative.

### 5.3 The two open questions §4.1 asked

- **Does `N₀` depend on depth?** **No, and cleanly.** The best value is **0** on the strong
  rare-allele pile-up at 3, 8 and 20 reads a sample, and **200** on the moderate one at all three.
  Depth moves it not at all; the two population shapes disagree by at least two hundred-fold. The
  spec's leaning was that depth would matter. **An earlier version of this sweep said depth moved
  it** — that was the draw-to-draw scatter §5.4 describes, and correcting the scoring removed it.
- **Is the blend ever worse than both ends?** **It cannot be, and an earlier version of this
  program said otherwise.** In log space the blend is a convex combination of the two errors, so
  per drawn cohort it is never outside them. Over 168 held-out cohort-and-panel rows at the shipped
  constant: at least as good as both in 8, between the two in 160, worse than both in 0. The
  earlier "worse than both in 1" compared three separately-taken medians, which can produce that
  and means nothing.

### 5.4 What each cohort is scored against, and why it is not the density

**Each drawn cohort is scored against the alternative-allele frequency it was itself drawn with** —
the share of that panel's chromosomes that ended up carrying the alternative allele — and not
against its density's expectation.

The difference is larger than anything the tables above resolve. Over 3,000 positions the realised
frequency of a Beta(0.70, 2.50) cohort has a standard deviation of **7.7%** of its own mean, and a
Beta(0.20, 1.00) one **10.3%**; the candidate half-weight panel sizes are separated by less than a
tenth of a per cent. Scoring against the expectation puts that scatter into every error column,
where it belongs to the draw rather than to either guess — and taking the median over six cohorts
does not remove it. **An earlier version of this sweep did exactly that**, and the two things it
changed are recorded above: the best value appeared to move with depth, and it appeared to be flat
below a quarter of an individual rather than rising from zero.

### 5.5 What this is not

**Drawn cohorts, not a real one.** This checkout cannot rebuild the tomato census — the read files
are not in the repository — so nothing here is a confirmation on real data. Spec §7's second open
question keeps that open, and it names the experiment: subsample the tomato panel and refit.

**The drawn cohorts segregate at 10% of positions**, well above either benchmark cohort, because a
tomato-like 4 in 1,000 would put about a dozen variable positions in a cohort of this size and
would measure the draw rather than the panel. §5.1's control is at the realistic 4 in 1,000, and it
is the arm with no sampling in it.

### 5.6 A side effect worth naming: the seed stops caring what the inbreeding coefficient is

`examples/ng_inbreeding_sensitivity.rs` holds a panel at a known inbreeding coefficient and tells
the fit a wrong one. Before this change the seed was the search's own pair, and the error went
straight through. Now it does not.

| panel | true F | F used | the search's `α_ref` | the shipped `α_ref` |
|---:|---:|---:|---:|---:|
| 1 | 0.60 | 0.50 | 0.667 | 1.0027 |
| 1 | 0.60 | 0.60 | 1.001 | 1.0027 |
| 1 | 0.60 | 0.70 | 2.005 | 1.0027 |
| 63 | 0.85 | 0.75 | 0.985 | 1.0172 |
| 63 | 0.85 | 0.85 | 1.001 | 1.0029 |
| 63 | 0.85 | 0.95 | 1.018 | 0.9866 |

**At one individual the seed does not move at all**, to five digits, while the search's answer moves
by a factor of three. The reason is that at one individual a wrong `F` rescales the search's pair
without moving its *ratio*, and the ratio is the only thing the seed keeps — the total it would
have taken is the one the pin replaces. At 63 individuals a `±0.10` error in `F` moves the shipped
reference concentration by at most 1.4%, against 4% for the search's own.

**This does not close `ordinary_site_seed.md` §7's third open question** — whether inbreeding
should enter the *weight* — but it changes what that question is worth: on the ordinary-site path
`F` now reaches the seed only through the shape, and at a single genome it does not reach it at all.

---

## 6. Tests

Sixteen new tests, and two rewritten. What each would take to be false is stated in its own doc
comment; the ones worth naming here:

| test | what breaks it |
|---|---|
| `the_seeds_implied_diversity_is_the_measured_one_at_every_panel_and_shape` | any error in the solved total, over five densities × five panel sizes, checked against the module's own spectrum machinery rather than against the identity that produced it |
| `the_solved_total_is_what_the_beta_binomial_needs` | the same, against a bisection oracle over a grid of frequencies and ceiling-shares |
| `the_blend_is_geometric_and_reaches_both_ends_exactly` | an arithmetic blend, a swapped weight, or either end taken alone |
| `the_bigger_the_panel_the_more_of_its_own_shape_the_seed_takes` | a weight that does not rise with the panel — on a spectrum whose shape is deliberately ten times the neutral one, because **on a neutral fixture the two ends of the blend are the same number and the weight could be anything** |
| `a_neutral_panel_projects_to_one_and_theta` | the search failing to read the neutral shape back; the pin failing to reproduce the measurement; the pinned pair drifting from `(1, θ)` by more than the `3θ` the pin costs |
| `a_fully_invariant_cohort_at_a_measured_diversity_falls_to_the_neutral_rung_and_says_so` | §3.1's first refusal going silent |
| `a_heterozygosity_above_a_half_refuses_the_run` | §3.1's second |
| `a_cohort_with_no_variation_is_floored_and_says_the_diversity_was_zero` | §3.1's third |
| `a_measurement_exactly_at_the_shapes_ceiling_has_no_total` | the ceiling comparison written as strictly-greater, which panics instead of reporting |
| `two_panels_that_leaned_differently_emit_different_records` | the weight not reaching the run's output |

**The fixture accident this step found in its own inheritance.** `projection_tests::project` handed
the projection a hard-coded diversity of `1e-3` whatever spectrum it was given, and three tests in
`run_parameters` passed `None`. Both were invisible while the diversity was not read. Making it
read changed six tests, and every one of the six was a fixture saying something it did not mean.

### 6.1 The mutation battery

Fifteen mutations of the shipped code, each run against the tests under
`ng::calling::genotype_prior`. **Fourteen were
caught on the first pass and one survived**, and the survivor is the one worth recording.

| mutation | tests that failed |
|---|---:|
| the blend's two weights swapped | 3 |
| an arithmetic blend instead of a geometric one | 2 |
| the weight is `N₀ / (N + N₀)` rather than `N / (N + N₀)` | 4 |
| the solved total is `t` rather than `t / (1 − t)` | 5 |
| the seed's two concentrations swapped | 5 |
| the panel size is a constant, not the spectrum's | 3 |
| the shape is the panel's own, never blended | 1 |
| the shape is the neutral one, never blended | 2 |
| the ceiling comparison written strictly-greater | 1 |
| the zero-diversity branch removed | 1 |
| the impossible-diversity refusal removed | 1 |
| the search's shape is `α_ref / total` rather than `α_alt / total` | 4 |
| the search's total is its reference concentration alone | 2 |
| the zero-diversity seed floors its reference concentration too | 1 |
| the unreachable-diversity fall-back keeps the fitted shape | 1 |

**The survivor: writing the neutral shape's expected frequency as `θ` rather than `θ / (1 + θ)`.**
All 143 tests passed. The two differ by a factor of `1 + θ` — 1 part in 1,000 at a human
diversity, 0.06% at tomato's — which sat inside every tolerance in the module, because **every
diversity in every fixture was `10⁻⁴` to `10⁻²`**. That is the shared accident: the tests covered
four decades of diversity and all four were small.

It matters despite being small: the neutral end of the ramp has to be *exactly* the pair the
no-spectrum branch returns, or the two rungs `population_diversity.md` §3.4 switched between are
not the two ends of one ramp. The fix is a named `neutral_expected_frequency` and a test that
compares it against the ratio of the seed that branch actually builds, at diversities up to 0.4 —
where the mutation is wrong by 40%.


### 6.2 What the review's own battery found, and the two Blockers in it

A review agent ran **42 further mutations** against the same tests. **Ten survived**, two of them
Blocker-class, and every one is closed now — each with a mutation re-run to prove it.

**Blocker: squaring the blend's weight passed all 144 tests.** The seed would take 0.64 of its
shape from the panel at one individual while `SeedRegime::FittedSpectrum` reported 0.80 — the seed
and its own record disagreeing, which is exactly what the spec's goal 3 forbids. It survived
because `the_bigger_the_panel_the_more_of_its_own_shape_the_seed_takes` bounded the share only
loosely: below 0.85 at one individual and above 0.99 at 63, and 0.64 and 0.992 are both inside.
**A band wide enough to hold both ends of a ramp is wide enough to hold a different ramp.** The
test now asserts the share *equals* `panel_shape_weight(N)` to 5 parts in 1,000 at every panel
size; the shipped code achieves 8 parts in 10,000, so that leaves six-fold headroom and refuses the
square by a factor of forty.

**Blocker: the survivor §6.1 records came back one line away from where it was fixed.** Writing
`measured_diversity` where `neutral_expected_frequency(measured_diversity)` belongs — at the *call
site*, leaving the function itself correct — passed all 144. `the_ramps_neutral_end_is_the_pair_the_neutral_rung_returns`
tests the function, and the branch it compares against never calls it, so nothing touched the line
where the value is consumed. **No test of the pin can ever see this**, because the total is
re-solved from whatever the shape turns out to be — the implied heterozygosity comes back exactly
`θ` under the mutation too. `a_neutral_panel_projects_to_one_and_theta` now has an arm at `θ = 0.1`
and `0.4`, where the mutation lands 68% from the pinned neutral pair against the shipped code's
0.4%.

**The other eight, each now caught:**

| what survived | what closes it |
|---|---|
| `FittedShape::concentrations` swapped, or dropping `(1 − f)` on the reference side — **no test touched it at all**, and its two readers produce §1.2's headline figures | a test that the pair it rebuilds is the one the search returned |
| `its_own_diversity` put back to the hard-coded `1e-3` — **the fixture fix §6 above boasts about was not load-bearing**; it changed which constant was wrong | `the_same_density_at_two_panels_is_two_different_spectra` now asserts both seeds imply the density's own heterozygosity, which its doc already claimed |
| `DiversityUnreachable` reporting the search's frequency instead of the blended one | asserted strictly between the two ends; the search's `1.0e-9` satisfies the ceiling inequality just as well as the blend's `6.7e-8` |
| `census_sites_outweigh_regularizer` written `≥` | a fixture row where the two are equal |
| `MAX_IMPLIED_DIVERSITY` set to 0.8 — it was pinned only to the interval `[0.5, 0.9)` | a refusal test at `0.500001` |
| `fit_spectrum_shape` switched to `SearchPrecision::fine()`, tripling every run's cost | the prediction-budget test now checks the wrapper against `fit_pair` at `fast()` |
| clamping the blended shape at a half | a panel whose alternative allele is the **common** one — the only fixture in the module on that side of a half, and a real case on a crop reference |
| `its_own_diversity`'s `.ok()` swallowing an out-of-range value into the fallback path | `.expect()` |

**And a recorded count was wrong**: §6.1's "the seed's two concentrations swapped" fails 4 tests,
not 5.

**The shape all the fixtures shared, in three sentences.** Every diversity was between `10⁻⁴` and
`10⁻²`, so a factor of `1 + θ` hid. Every expected frequency was below 0.29, so which side of a
half the shape sat on was never varied. And every spectrum built from `(1, θ)` has the two ends of
the blend at the same number, so on those the weight could be anything at all — which is why the
one fixture built from a shape ten times the neutral one is carrying the whole of the ramp's
correctness.

---

## 7. Validation

All run through `./scripts/dev.sh` from this worktree.

| command | result |
|---|---|
| `cargo test --lib` | **4,874 passed / 0 failed / 14 ignored** (4,858 at `1474c5cc`) |
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo doc --no-deps --lib` | **27** unresolved intra-doc links, against the branch's 28. **None added**; the one that went was `projection_tests::the_projection_returns_one_pair_at_every_inbreeding_coefficient`, in a doc comment this change rewrote |

**The `--release --lib` suite was run too**: 4,850 passed, **8 failed**, 14 ignored at the time it
was run. All eight are
pre-existing and all eight are the same thing — a `#[should_panic]` test on a check that is a
`debug_assert!`, which release compiles out — in `genetics`, `ng::alignment`, `sample_summary`,
`ssr::cohort` and `var_calling::dust_filter`, none of which this work touches. **None of this
change's own release-held checks is among them**, which is what that run is for: the one
release-held assertion added here is the refusal of a heterozygosity above a half, and demoting it
to `debug_assert!` makes `a_heterozygosity_above_a_half_refuses_the_run` fail in release and pass in
debug, which is how it was confirmed reachable.

The release total is below the debug one because twenty items in the tree are behind
`#[cfg(debug_assertions)]`; none of them is new here.

`cargo test --all-targets` fails in `benches/psp_writer_perf.rs`, also pre-existing.

---

## 8. Trade-offs and follow-ups

**⚑ Two decisions for the owner, and the first is the one that matters.**

**1. The ramp has no measured support, and what ships is a hedge.** Every arm of the measurement
puts the best half-weight panel size at **zero** — which is a weight of exactly one at every panel,
the blend inert and §4's ramp a no-op. What ships is `0.25`, the smallest value that keeps the
mechanism alive; it costs 2.7% on the held-out drawn cohorts, 1.9% on the ones it was fitted on and
26% on the arm with a real cohort's share of variable positions (§5.2). **What it buys is a guard
this sweep cannot price**: a fifth of the shape kept on the neutral side at a single genome,
against a pathological fit at a small panel — `ordinary_site_seed.md` §1.1's concern, which does
not reach the quantity the seed reads and which every drawn density here is too well behaved to
produce. **Recommendation: keep 0.25 until a single real genome has been fitted end to end**, then
either confirm the guard is doing something or set the constant to zero and delete the ramp. The
alternative — take the fit's answer now — is defensible and is a one-line change.

**2. §4's premise is contradicted, and the fix is not in this file.** The panel's own shape
degrades with the panel because the shape is read back *out of* allele-count classes rather than
off the density that produced them. **On the joint route the density is right there** —
`FrequencyDensity` carries `a` and `b`, and the expected frequency is
`p_segregating · a / (a + b) + p_fixed_alt`, exactly the quantity §5.1 shows the refit losing.
Taking it from the density instead of refitting would make the panel-size dependence disappear, and
with it the reason for a ramp at all. What stops it here is that the seam takes class weights
rather than a density, and spec §2 lists changing the search as a non-goal. **Recommendation: raise
it against `calling_priors.md` §4.1 as its own step** — it is a change to what the seam carries,
and it would retire this ramp rather than retune it.

**Smaller things, none of them blocking.**

- **A measurement within one part in `10¹⁶` of the shape's ceiling gives a total of `9.0e15`** — a
  prior no depth of reads could move, and nine orders above `MAX_PROJECTION_CONCENTRATION`, which
  is the bound this module states for its own spectrum machinery. Nothing predicts a spectrum from
  a seed today, so it is latent. It is recorded rather than guarded, because clamping it would
  break the pin for a case no fit can produce.
- **The pin is the seed's guarantee and `fill_locus_concentration` does not extend it.** That
  function floors the *alternative* entry alone, so below a measured diversity of about `1e-13` it
  moves the locus row's ratio rather than its scale — the homozygous-alternative weight goes from 5
  in 10,000 at `θ = 2e-12` to 0.48 at `1e-15`. The smallest diversity a real cohort could fit is
  one variable position in the census, `1e-7` at the two-million-position budget, so this is
  unreachable; the comment at the seed says where the guarantee stops rather than implying it
  carries.
- **Spec §6's checklist needed three amendments** and they are recorded in the spec itself: item 2's
  0.15% is the *diversity's* agreement and the pair moves about `3θ`; item 3's "better than either
  end" is not testable in the direction it states, because per cohort the blend is a convex
  combination of the two errors; item 6's second half is arithmetic rather than a test.
- **Nothing here has been run on a real cohort.** Spec §7's second open question, unchanged.
- **Whether inbreeding belongs in the weight** is spec §7's third question, untouched: every
  measurement here is at `F = 0`. §5.6 changes what that question is worth rather than answering it.
