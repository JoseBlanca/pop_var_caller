# ng genotype prior — D2: the run's two starting numbers, read off the fitted spectrum

*Implementation report, 2026-08-22. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step D2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone D.*

## 1. What it is

`project_spectrum_seed` turns the pre-pass's fitted frequency spectrum — one weight per
allele-count class, *what share of sites carry the alternative allele on exactly `j` of the
panel's chromosomes* — into the two numbers the genotype prior starts every locus from: the
reference allele's concentration and the total shared out across whatever alternative alleles a
locus turns out to carry.

It is a **change of representation, not a second estimate**: the neutral `1/p` density and the
neutral frequency spectrum are the same statement written twice, once at a locus and once across a
panel, so nothing is fitted here that the pre-pass has not already fitted
([`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §4.1).

The fit maximises the likelihood of the fitted spectrum's class weights under a candidate pair's
predicted spectrum — equivalently minimises the Kullback–Leibler divergence — **over every class
including the monomorphic one**, predicting with step D1's `fill_expected_spectrum`, which carries
the panel's inbreeding.

Design authority: spec §4, §4.1, §4.2, §12 tests 5–7;
[`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §2.3, §4.

## 2. The plan's impl-time confirmation, answered

The plan and arch §8 left open **which concrete type the projection consumes**, offering
`FrequencyDensity` ([`joint/fit.rs:87`](../../../../src/ng/parameter_estimation/joint/fit.rs)) or
the pre-pass cohort gather's own wrapper.

**It is neither, and the reason is that they are different objects.** `FrequencyDensity` is four
numbers describing how the **population's** allele frequency is distributed — two point masses and
a Beta over what segregates. The projection matches **class weights**: how a panel's `2N`
chromosomes came out. The cohort gather that will produce those does not exist in the code yet
(nothing under `src/ng/parameter_estimation` names a spectrum).

So `FittedSpectrum<'_>` is **a borrowed view, not a copy**, as the plan asks: it wraps
`&[f64]` class weights plus the two numbers spec §4.1 requires the run to report — how many sites'
worth of pseudo-counts held the estimate at the neutral shape, and how many census sites came out
variable. When the gather lands it owns the weights and this borrows them.

**The panel size is read off the class count** rather than taken as an argument. Arch §4's sketch
has a `chromosomes: u32` beside the spectrum; `2N + 1` classes already fix `N`, and a second
argument is a second place for it to disagree. The constructor refuses an even class count.

## 3. Two departures from arch §4's signature, both to make the plan's own rule expressible

- **`diversity: Option<ExpectedHeterozygosity>`**, not a bare value. The plan's D2 line requires
  *no fitted θ → `FallbackDiversity`*, and a function handed a value cannot tell a fitted `1e-3`
  from the fallback `1e-3`. With the option, the regime is derived inside rather than asserted by
  a caller.
- **No `chromosomes` argument**, for the reason above.

## 4. The optimiser is not `fit_by_multistart`, and that is a measurement not a preference

Arch §4 says the projection "reuses the pre-pass's fitting machinery
([`fitting/multistart.rs`](../../../../src/ng/parameter_estimation/fitting/multistart.rs)) — a
two-parameter fit needs nothing new."

`fit_by_multistart` scores **one cell at a time** through
`NoiseModel::append_genotype_likelihoods`, which takes `&self` and so cannot cache. The natural
cell here is one allele-count class, so every class would rebuild the entire spectrum: at 3,200
individuals that is **6,401 predictions where one is needed** — about 1.7 hours per candidate
against 0.96 seconds. It also has no notion of a search direction that is not an axis, which §5
below needs.

What is reused is the **shape** — several starts, coordinate line searches by golden section on
each axis's own scale, the best point evaluated kept rather than the last reached, a capped search
reported rather than asserted — and `SearchPrecision` itself, whose `fast()`/`fine()` reasoning
transfers unchanged: 1% of a concentration is finer than a genotype likelihood registers.

## 5. The search's coordinates, which are the one thing that was not obvious

**Coordinate descent on `(ln α_ref, ln α_alt)` finds the answer at 26 and 63 individuals and fails
at one.** Two independent implementations agree on that and disagree on the size — 0.844 in mine
and 0.244 in the review's, with the starts spread 3,206-fold and 7,421-fold — because the answer
depends on the search box each used and neither recorded it. **Take the shape, not the digits:**
in those coordinates one individual comes back somewhere between a quarter and five-sixths of the
answer, and 26 and 63 come back right.

Searching the **total and the ratio** instead does the opposite: 7 parts in ten million high at
one individual, 1.3% high at 63.

The reason is that these are the same two quantities twice, rotated 45° in log space. Writing
`t = ln(α_ref + α_alt)` and `r = ln(α_alt/α_ref)`, and with the ratio small:

```text
  [1,  0]      the total, both concentrations moving together
  [0,  1]      α_alt alone — α_ref does not move
  [√½, −√½]    α_ref alone — α_alt does not move
```

So **three directions, one for each quantity the two parametrisations name between them**, and
neither panel size is a special case. With them the fit returns `α_ref` within 0.25% of 1 at every
panel size from 1 to 150, at three diversities and two inbreeding coefficients.

**A fourth direction, `[√½, √½]`, was in the first version and the review removed it.** It is a
coordinate of neither parametrisation — along it `ln α_alt` moves twice as fast as `ln α_ref` —
and it did harm rather than nothing: on a spectrum with all its mass at intermediate frequency,
the shape spec §4.1 says two parameters cannot hold, the sweep including it ended at
`α_ref = 3.59` where the best point in the box is 498 (log-likelihood −3.410 against −2.232). It
also cost up to 314 predictions a fit. The mutation that deleted it was one of only two survivors
of thirteen.

## 6. What a fit costs

**399 predictions, the same at every panel size and every inbreeding coefficient measured**,
asserted in `a_fit_costs_at_most_450_predictions`. Measured wall clock in release, one fit, on the
exact spectrum of `(1, θ)` at `θ = 6 in 10,000` and `F = 0.8`
(`the_cost_of_one_fit_by_panel_size`, `#[ignore]`d):

| individuals | 400 | 800 | 1,600 | 3,200 |
|---|---|---|---|---|
| one fit | 3.8 s | 22 s | 2.2 min | **11.8 min** |

About `N^2.5` against D1's `N^2.45` for one prediction — inside the spread of a four-point fit, so
the search adds a constant factor rather than a power.

**This is larger than the 2.6 minutes the plan's D1 line and
[`spectrum_projection_cost_2026-08-22.md`](../../ng/reports/spectrum_projection_cost_2026-08-22.md)
predict for 3,200 samples, and the per-prediction measurement in that report is not what moved.**
Two things did. The prediction count is 399 rather than the "about 160 objective evaluations" that
report assumes, and **a prediction averages 1.78 s inside a fit against 0.96 s at the neutral
pair** — the search spends most of its predictions away from that pair, where the branch-tail trim
drops fewer splits. Raised for the milestone review at the end of step D; the cheapest lever if it
matters is fewer starts.

## 7. Tests

Twenty-one new, one of them `#[ignore]`d. Every projection target is **the exact expected spectrum
built by D1, never `θ/k` and never sampled** — `θ/k`'s own error is 0.272% at tomato's diversity
and 4.4% at a `θ` of 1 in 100, larger than what these tests measure (spec §12 test 5).

| test | what it pins | spec §12 |
|---|---|---|
| `a_neutral_panel_projects_to_one_and_theta` | the pair comes back at `(1, θ)` over 6 panel sizes, 3 diversities and 2 inbreeding coefficients — worst 0.25% on `α_ref`, 0.31% on `α_alt`, inside the 1% the search resolves; and 2 parts in 100,000 at a thousand-fold finer resolution, which is what says that residue is the resolution and not a bias | test 5 |
| `the_projection_returns_one_pair_at_every_inbreeding_coefficient` | one density gives one pair at `F = 0, 0.6, 0.8, 0.9`, and the independent-chromosome projection still returns its own biased answer — with `F = 0` as the comparison's zero, where the two predictions are the same model | test 6 |
| `at_one_individual_the_projection_is_still_the_neutral_pair` | the low end of the committed cohort range, at `F = 0` and `F = 0.9` | test 7 |
| `an_absent_spectrum_is_the_neutral_pair_at_the_fitted_diversity` | `(1, θ)` exactly, regime `NeutralShape` | §4.1 |
| `no_fitted_diversity_falls_back_to_the_species_value_and_says_so` | regime `FallbackDiversity` | §4 |
| `the_regularizer_weight_and_whether_the_census_sites_outweighed_it_travel_with_the_seed` | both, in both directions | §4.1 |
| `both_variant_classes_get_the_same_seed_today` | Q1's seam: the argument exists, the behaviour does not yet split | §4.2 |
| `every_start_reaches_the_same_pair_within_one_percent` | the four starts swept individually to convergence land within 1%, the widest disagreement 1.0039-fold at 26 individuals — what makes "only the best start continues" safe | — |
| `the_objective_is_maximised_by_the_truth` | Gibbs' inequality, an independent check of the objective that uses no spectrum at all | §4.1 |
| `the_winning_score_is_the_spectrums_own_entropy` | that the fit reaches the maximum, not merely a point near the right pair — at the true pair the divergence is zero, so the score is the entropy | §4.1 |
| `an_impossible_class_is_floored_rather_than_sent_to_negative_infinity` | `0 × ln 0` is skipped, not a `NaN`; an occupied class predicted at zero pays the floor, stays finite, and still ranks below every candidate that can produce it | — |
| `a_fully_inbred_panel_whose_spectrum_holds_heterozygotes_still_returns_a_pair` | the blocker of §8 below, pinned | — |
| six `#[should_panic]` tests | every guard on `FittedSpectrum::new`: an even class count, one class, a panel past the projection's range, a negative weight that still sums to one, and each of the two site counts | — |
| `a_fit_costs_at_most_450_predictions` | the count above, asserted so a search that gets more expensive has to say so | — |
| `the_cost_of_one_fit_by_panel_size` | the wall clock in the table above, `#[ignore]`d | — |

**The independent-chromosome numbers this step measures, on a panel of 26 individuals at
`θ = 6 in 10,000`** — spec §4.1 states the first and third of these and this reproduces both:

| `F` | 0 | 0.6 | 0.8 | 0.9 |
|---|---|---|---|---|
| two-branch `α_ref` | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| independent-chromosome `α_ref` | 1.0000 | 0.9144 | 0.8793 | 0.8599 |
| how far low | — | 8.6% | 12.1% | 14.0% |
| independent-chromosome `α_alt` | 1.000 θ | 0.893 θ | 0.848 θ | 0.824 θ |

The last row settles a reading spec §4.1 leaves open: a fit against an independently-called VCF of
18 tomato accessions returned `α_alt = 0.81 θ`, and it was worth asking whether a domesticated
selfer stretches the two-parameter family. **A perfectly neutral panel run through an
independent-chromosome projection returns 0.824 θ by itself at `F = 0.9`**, so that number is
consistent with being nothing but the bias this step removes, and says nothing about tomato either
way.

## 8. What the review found, and what is still open

Four agents in isolated worktrees — mathematics and numerics, reliability and contracts,
conformance to the settled design, and naming with prose and smells. Twenty-seven mutations
between the first two, of which six survived and were each proved to change behaviour.

### The blocker: a fully inbred panel aborted the run, blaming the wrong number

Both agents found it independently. At `F = 1` the prediction puts **exactly zero** in every odd
allele-count class, so a fitted spectrum carrying any weight in one of them scores `−∞` at every
candidate pair — 441 of 441 points across the search box, measured. No start could then beat the
sentinel the search began from, the sentinel was `NaN`, and the run died three frames later
saying *the reference concentration must be finite … got NaN*.

Both inputs are legal today. `InbreedingF::try_new(1.0)` returns `Ok` — the `[0, 1)` tightening
spec §7 requires is [`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md)
Milestone A, which has not been started — and a census spectrum with singletons is ordinary. The
two arrive from different pre-pass accumulators, so nothing makes them consistent.

**Fixed by flooring the prediction before the logarithm** at `PROBABILITY_FLOOR`, which is spec
§8's own rule for every logarithm in this module. That is not only a `NaN` repair: `−∞` is not an
ordering, so golden section compares `−∞ > −∞` as false and walks blind. **That region is not
hypothetical at ordinary inbreeding either** — on the exact neutral spectrum at `F = 0.8`, none of
441 points across the box scores `−∞` at 150 individuals, 17 do at 400, and 28 of 225 do at
1,600, which is the middle of the committed cohort range.

### The best point evaluated was not the one returned, and the doc said it was

A line search ended at its bracket's midpoint without ever comparing that against where it
started. Measured: 5 of 80 line searches on the exact neutral spectrum at 26 individuals ended
below their start, and 31 of 80 on a flat spectrum. Keeping the better of the two costs no
prediction, because both scores are already in hand. The mutation that made a whole sweep return
its last point rather than its best changed the answer on 6 of 20 spectra, worst
`α_ref` 64.22 → 3.25, and it survived the test set.

### Four guards had no test, and each was proved to matter

Every one survived a mutation that deleted it, and each was shown to return a plausible answer
rather than fail: a spectrum of one class (a panel of no individuals) returns `α_ref = 996.87`,
the box's top corner; weights of `[1.5, −0.5, 0]`, which sum to one, return the bottom corner; a
`NaN` regularizer reaches the run's output *and* reports that the prior dominated, because
`NaN > x` is false. Six `#[should_panic]` tests now hold them, and the class count gained a
ceiling — at 10,000 individuals a fit is already about two hours, and nothing on the way would say
the run had stopped being a run.

### Still open, and it is the one thing this step does not do

**A fit that could not match the spectrum is emitted as though it had.** Two cases reach it: a
spectrum no pair in the family can produce (the `F = 1` case above), and one whose best pair lies
outside the search box — a fully invariant cohort lands on the box's floor at `α_alt = 1e-12`, and
a spectrum with its mass at intermediate frequency lands on the ceiling at `α_ref = α_alt = 500`.
In every case `SeedRegime::FittedSpectrum` comes back looking like an ordinary fit.

That is the same complaint spec §12 test 11 makes about the STR seed — *what it must not do is
return the closest total it can reach as though it had met the target* — and the same complaint
spec §4 makes about a run on the fallback diversity. **The fix is a field on
`SeedRegime::FittedSpectrum`, which is arch §2.3's type**, so it reaches past this step and is
raised for the milestone review rather than taken here.

Worth knowing about the scope: for spectra genuinely drawn from the two-parameter family, `α_ref`
comes back within 0.6% of 1 across `θ` from 1e-9 to 6e-4 at 1 and 26 individuals, and the fits
this module's tests run sit far from both walls.

## 9. Corrections owed to the design documents — raised, not applied

1. **The independent-chromosome bias, in `arch/calling_priors.md` §4 and this plan's D2 line**,
   both of which say 9–14% at tomato's fitted `F`. Measured here: **8.6% at `F = 0.6`, 12.1% at
   0.8, 14.0% at 0.9**. Spec §4.1's "12 to 14% at tomato's fitted `F` of 0.8 to 0.9" is right.
   (Owed since D1; now measured by a test rather than quoted.)
2. **Spec §4.1 line 329 contradicts spec §12 test 5.** It says the gap between `θ/k` and the exact
   spectrum "is why §12's test builds its target by *sampling* the Dirichlet rather than by writing
   `θ/k`". Test 5 forbids sampling in as many words — "**And not by drawing sites at random
   either**" — and requires the closed form, which is what D1 built and what these tests use. One
   of the two sentences has to go, and it is the §4.1 one.
3. **Two sentences ask for an exactness no bounded-resolution search can deliver**, and they are
   not the three the first draft of this report named. Spec §12 **test 5** — "the projection must
   return `α_ref = 1` and `α_alt = θ` to floating-point tolerance at several panel sizes" — and
   spec **§4.1** — "At one sample the projection returns `(1, θ)` **exactly**". Tests 6 and 7 ask
   only for one pair from three densities and for `(1, θ)` at one sample, which are met. What is
   achievable, and what the tests assert: the answer inside the resolution asked for (0.45% at the
   shipped setting) and shrinking with it (2 parts in 100,000 at a thousand-fold finer), which is
   the statement that separates a resolution from a bias.
4. **Arch §4's "the optimiser reuses `fitting/multistart.rs` — a two-parameter fit needs nothing
   new."** Measured false, §4 above. `SearchPrecision` is reused; the driver is not.
5. **The cost of a fit at the top of the cohort range**, §6 above — 12.7 minutes at 3,200
   individuals against the 2.6 minutes recorded in the plan's D1 line and in
   `spectrum_projection_cost_2026-08-22.md`, which assumed "about 160 objective evaluations".

Four more the review added, none of them this step's to fix:

6. **No document says where the projection's single panel `F` comes from.** Spec §4.1 requires the
   prediction "at the panel's `F`", and
   [`parameter_prepass_cohort.md`](../../ng/spec/parameter_prepass_cohort.md) §4.1 says "at the
   panel's inbreeding from §3" — but §3 is per sample, and `JointFit` carries one
   `HomozygoteExcess` per sample. Nothing defines the aggregation. This is not cosmetic: this
   step's own test measures `α_ref` moving from 1.000 to 0.860 as `F` goes 0 to 0.9, so the wrong
   panel `F` is a 14% error in the number §4.1 is about, and a cohort mixing selfers with
   outcrossers has no single `F` to take. Owed by `parameter_prepass_cohort.md` §3 or §4.1.
7. **Arch §4's `seed_for_locus` signature omits `VariantClass`; plan D3 requires it.** Arch writes
   `seed_for_locus(seed, allele_count, out)`; the plan says "`VariantClass` stays an argument even
   while both classes pass one θ". The plan is right — `seed_for_locus` is the port of
   `alpha_from_diversity` ([`genetics.rs:214`](../../../../src/genetics.rs)), which is where
   production's 8:1 split actually lived — so arch §4's signature should gain the argument. **Step
   D3 hits this**, and takes the plan's side.
8. **`parameter_prepass_cohort.md` §4.1 should name the exact expected spectrum as the
   regularizer's shape, not `θ/k`.** It says the neutral shape is "the neutral frequency density
   put through the same two-branch sampling", which is `fill_expected_spectrum` — but it also
   writes the shape as `θ/k` repeatedly, and spec §4.1 says the one-sample "spectrum is its `θ/k`
   prior untouched". A pre-pass that implements `θ/k` literally stops the one-sample projection
   returning `(1, θ)` by 0.272% at tomato's diversity and 4.4% at a `θ` of 1 in 100 — the gap this
   module's own `the_neutral_shape_appears_in_the_small_diversity_limit` measures.
9. **Arch §2.3 should say at the field that `data_dominated` is the aggregate**, and that the
   per-class ratio spec §4.1 asks for is the pre-pass's to emit. Arch §4 does assign it there, but
   in a subordinate clause of a different section, and spec §4.1 is explicit that the panel-wide
   ratio is the wrong number to quote as reassurance.
