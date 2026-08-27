# ng — the ordinary-site prior's two numbers: implementation report

**Branch** `ng-prior-moments`, worktree `../pop_var_caller-prior-moments`, cut from `main` at
`9f15f5e5`.

**Plan** [`../../ng/impl_plan/ordinary_site_prior_moments.md`](../../ng/impl_plan/ordinary_site_prior_moments.md).
**Design authority** [`../../ng/spec/ordinary_site_prior_moments.md`](../../ng/spec/ordinary_site_prior_moments.md),
with [`../../ng/spec/ordinary_site_seed.md`](../../ng/spec/ordinary_site_seed.md) §3 for the
identity that turns two moments into a concentration pair.
**The measurements this work rests on** were made before it started and are in
[`../ng_ordinary_site_prior_moments_2026-08-27.md`](../ng_ordinary_site_prior_moments_2026-08-27.md);
no step here is a sweep.

**One report, one section a step.** The plan's steps are small and several are deletions, so a
file apiece would be a file of two paragraphs; each section below carries the step's own contract,
what shipped, and what was measured about it.

---

## A1 — the population's mean alternative-allele frequency, in closed form

**Contract (plan A1).** Add `p_fixed_alt + p_segregating · a/(a+b)` beside the heterozygosity
integral that already exists, with the same population-not-panel framing. Its own commit, because
a wrong mean frequency is a plausible number at every panel size and nothing downstream refuses it.

### What shipped

`FrequencyDensity::expected_alternative_frequency`, in
[`src/ng/parameter_estimation/joint/fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs),
one line of arithmetic beside `expected_heterozygosity`. Two of the density's three parts
contribute: positions where the population carries only a non-reference base are at frequency one
and carry their whole share, positions that segregate contribute their Beta's mean, and positions
carrying only the reference base contribute nothing.

Nothing calls it yet. A2 is where the seed starts reading it.

### The oracle: the search's own answer, where the search is exact

**The cheapest available proof that the replacement is right is the thing it replaces.** At one
diploid individual the panel has three allele-count classes, two of them free once normalised,
against the two-parameter family's two parameters — so the search reproduces those classes exactly
and its mean frequency is the density's own (spec §9, report §9's third reading). At larger panels
it is fitted over more classes than it has parameters and drifts, up to 1.22× the truth at 200
individuals, so the comparison is available only here.

Measured, on six densities at one individual and no inbreeding, the search's mean frequency over
the closed form's:

| density | search | closed form | ratio |
|---|---:|---:|---:|
| tomato-like, `Beta(0.20, 1.00)` | 1.664544e-3 | 1.666667e-3 | 0.9987 |
| human-like, `Beta(0.35, 1.20)` | 1.459505e-3 | 1.461290e-3 | 0.9988 |
| flat, `Beta(1.00, 1.00)` | 2.999262e-3 | 3.000000e-3 | 0.9998 |
| the unit tests' lopsided fixture, `Beta(0.50, 2.00)` | 2.798057e-2 | 2.800000e-2 | 0.9993 |
| middling, `Beta(4.00, 4.00)` | 2.999262e-3 | 3.000000e-3 | 0.9998 |
| reference base rare, `Beta(3.00, 0.60)` | 4.337823e-3 | 4.333333e-3 | 1.0010 |

**The band is 0.9987× to 1.0010×**, an order of magnitude inside the search's own 1% resolution
(`SearchPrecision::fast`). The two routes share no algebra: one maximises a log-likelihood over
Beta-binomial class weights, the other is one line of Beta moments.

**This test dies at A5**, with the search it is measured against. What survives as the permanent
check is the hand-computed one below.

### What the fixtures share, and what was added because of it

**The five densities this repository already sweeps all have `a ≤ b`**, and reading `b/(a+b)` for
`a/(a+b)` does one of two things on them. On the three where `a < b` it returns a number too
*high*. On the two symmetric ones, `Beta(1, 1)` and `Beta(4, 4)`, it returns the **identical**
number — the swap is not merely consistent there, it is invisible. `Beta(3, 0.6)` — the population
where the reference base is the rare one at the positions that vary (report §2) — was added
because it is the only fixture on which the swap points the other way.

**⛦ An earlier version of this paragraph said the swap "returns a number too high on every one of
them"**, which is wrong on two of the five for the same reason the paragraph below gives: those
two are symmetric. Corrected after review.

**Six rows, five distinct answers.** `Beta(1, 1)` and `Beta(4, 4)` are both symmetric, so both have
mean a half and both give 3.000 in 1,000. They differ in spread and not in mean, so for this
quantity they are one fixture. The set is kept whole because it is the set the earlier measurements
used; the duplication is recorded rather than left to be found.

### Mutations run

Four, on the formula itself, each against all three of the step's tests
(`the_expected_alternative_frequency_is_the_densitys_own`,
`the_two_point_masses_carry_their_own_ends`,
`the_closed_form_frequency_is_the_searchs_own_answer_at_one_individual`):

| mutation | tests failing |
|---|---|
| `a` and `b` swapped | 2 of 3 |
| the fixed-non-reference share dropped | 3 of 3 |
| the invariant share read in its place | 3 of 3 |
| the segregating share dropped | 3 of 3 |

The swap leaves `the_two_point_masses_carry_their_own_ends` green, which is correct: at
`p_fixed_alt = 1` and at `p_invariant = 1` the Beta contributes nothing either way, so that test
pins the two ends and says nothing about the shape between them.

The source was restored from a backup and the restore checked with `git diff` before anything
else ran.

### Validation

All in the container, from this worktree:

- `cargo fmt --all -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,877 passed, 0 failed, 14 ignored**, against the branch's 4,874 at
  `9f15f5e5`. Three tests added, none removed.
- `cargo doc --no-deps --lib` — **27 unresolved links, the same 27 as at `9f15f5e5`**. The crate
  denies broken intra-doc links, so this command exits 101 on the pre-existing set; what matters is
  that the count did not move.

---

## A2 — the seed takes the two moments, and nothing else

**Contract (plan A2).** The seed builder's signature becomes the two moments — no spectrum, no
panel size, no inbreeding coefficient. Its body is A1's frequency and the existing
`total_for_diversity`. The three regimes stay what they are: a fitted curve, a fitted diversity
with no curve, neither.

### What shipped

- **`ExpectedAlternativeFrequency`**, a new type in
  [`src/ng/types.rs`](../../../../src/ng/types.rs) beside `ExpectedHeterozygosity`, with its own
  `DomainError` variant. **It is a type and not an `f64` for the reason every other newtype in that
  section is one**: both moments are probabilities in `[0, 1]`, both are population quantities with
  no panel in them, and both reach the seed builder in the same call — so as bare floats a swapped
  pair compiles and returns a seed no downstream check refuses. Their sizes do not separate them
  either: where the alternative allele is rare the heterozygosity is about twice the frequency, and
  where the reference base is the rare one it is far smaller.
- **`seed_from_population_moments`** replaces `project_spectrum_seed` in
  [`seed_generic.rs`](../../../../src/ng/calling/genotype_prior/seed_generic.rs). It takes
  `Option<ExpectedAlternativeFrequency>` and `Option<ExpectedHeterozygosity>` and returns the seed.
  **The rename is the implementer's**, not the plan's: the function projects nothing any more, and
  a name that says it would be the sort of retired sentence this area keeps finding six copies of.
- **`RunParameters::seed_from_moments`** replaces `RunParameters::project_seed` at the seam in
  [`run_parameters.rs`](../../../../src/ng/calling/run_parameters.rs).

### Three deviations from the plan's step boundaries, and why each was forced

**1. `SeedRegime::FittedSpectrum` became `SeedRegime::FittedCurve`, a variant with no fields —
in A2 rather than in A3 and A5.** All four of its fields were computed from the two arguments A2
removes: `shape_from_panel` from the panel size, `spectrum_match` from the search, and
`regularizer_site_weight` and `census_sites_outweigh_regularizer` from the spectrum. **Keeping any
of them for one commit would have meant committing a fabricated value into a field a run reports.**
The end state is the plan's; only which commit removes which field moved. The two regulariser
fields have no future producer at all — they described a spectrum emission the cohort gather no
longer makes (owner's ruling of 2026-08-27, `parameter_prepass_cohort.md` §4) — so they are gone
rather than owed.

**2. `SeedRegime::DiversityUnreachable` kept only `expected_frequency`**, for the same reason. The
variant itself goes at A4.

**3. Two tests died here rather than at A3**, because the mechanism they pinned left the seed at
A2: `the_bigger_the_panel_the_more_of_its_own_shape_the_seed_takes` and
`two_panels_that_leaned_differently_emit_different_records`. Three more went with the regime's
deleted fields: `the_regularizer_weight_and_whether_the_census_sites_outweighed_it_travel_with_the_seed`,
`census_sites_equal_to_the_regulariser_do_not_outweigh_it`, and
`a_spectrum_with_no_diversity_falls_to_the_species_range_guess` — the last replaced by
`a_frequency_with_no_diversity_falls_back_and_says_so`, which asks the same question of the new
signature.

**The ramp's own three functions are still standing**, uncalled, under a one-commit
`#[allow(dead_code)]` that names A3 — so that the ramp's deletion is its own commit and a bisect
can find it, which is what the plan asked for.

### What the existing tests could not have caught, and what was added

**Every projection test in this module builds its spectrum from `(1, θ)`.** On such a spectrum the
two ends of the blend are the same number, so the blend's removal moves nothing — which is why
`a_neutral_panel_projects_to_one_and_theta` and
`at_one_individual_the_projection_is_still_the_neutral_pair` passed unchanged through a commit that
deleted the blend. **This is the same fixture weakness the previous branch's review recorded**, and
it is recorded again here rather than repaired: those tests are about the search and go at A5.

Two tests were added at the seam, in `run_parameters`:

- **`a_fitted_density_seeds_the_run_from_its_own_moments`** asserts the seed's implied
  heterozygosity is the density's own to within 1 part in 10¹² **and that its expected frequency is
  the density's own** — the second is the one that matters, because the implied heterozygosity
  `2 f (1 − f) · A/(A+1)` is symmetric under `f → 1 − f`, so it cannot see a swapped pair.
- **`the_seed_is_what_the_identity_gives_from_the_four_fitted_numbers`** writes the whole chain out
  from the fixture density's four fitted numbers — both moments, then
  `t = θ/(2f(1−f))`, `A = t/(1−t)`, `(A(1−f), A f)` — and compares. It is the only check here that
  does not go through the library's own accessors.

### Mutations run

Two, on the seed builder, each against the whole library suite:

| mutation | tests failing |
|---|---|
| the two concentrations swapped — `α_ref` given the frequency's share | **4** |
| the regime reported as `NeutralShape` instead of `FittedCurve` | **3** |

The swap is killed by the two frequency assertions and by the two neutral-panel tests; it is **not**
killed by `the_seeds_implied_diversity_is_the_measured_one_at_every_shape`, for the symmetry
reason above.

The source was restored from a backup after each and the restore checked with `git diff` before
anything else ran.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,874 passed, 0 failed, 14 ignored**, from 4,877 after A1. **Five tests
  removed and two added**, all named above; the fall is the deletion the plan expects and not a
  regression.

  **⛦ That counts three renames among the "removed", and a review caught it.** Only
  `the_seed_is_what_the_identity_gives_from_the_four_fitted_numbers` is net-new;
  `a_fitted_density_seeds_the_run_from_its_own_moments`,
  `no_fitted_frequency_is_the_neutral_pair_at_the_fitted_diversity` and
  `the_seeds_implied_diversity_is_the_measured_one_at_every_shape` are the same tests under new
  names, two of them carrying new assertions. The arithmetic nets out to the same −3 and the suite
  really did move 4,877 → 4,874; what was wrong was the composition.
- `cargo doc --no-deps --lib` — **27 unresolved links, the same count as at `9f15f5e5`**. One new
  break was introduced and fixed: `HALF_WEIGHT_PANEL_SIZE`'s documentation linked to
  `SeedRegime::FittedSpectrum`, which no longer exists.

---

## A3 — the blend goes

**Contract (plan A3).** Delete `HALF_WEIGHT_PANEL_SIZE`, `panel_shape_weight`, the log-space blend
and `shape_from_panel`, and the three tests that pin the ramp; keep
`examples/ng_seed_shape_weight_sweep.rs` as the record of why, with its head rewritten to say it
measured a mechanism that is now deleted.

### What shipped

`shape_from_panel` went at A2, for the reason that section gives. This step deletes the rest:

- **`HALF_WEIGHT_PANEL_SIZE`**, **`panel_shape_weight`**, **`blend_expected_frequency`** and
  **`neutral_expected_frequency`** — the last because it was the blend's lower end and had no other
  caller. The `#[allow(dead_code)]` markers A2 left on the two private ones are gone with them.
- **The three tests**: `the_blend_is_geometric_and_reaches_both_ends_exactly`,
  `the_ramps_neutral_end_is_the_pair_the_neutral_rung_returns` and
  `the_weight_rises_with_the_panel_and_stays_inside_zero_and_one`. What the second protected — that
  the no-frequency branch returns exactly `(1, θ)` — is still asserted, by
  `no_fitted_frequency_is_the_neutral_pair_at_the_fitted_diversity`.
- **The re-export** of the constant and the weight function from `genotype_prior`.

**A comment stands where the ramp was**, carrying what it did, the one constant it had, and the two
measurements that retired it: **all three arms of the sweep's headline table — the one averaged
over panel size, depth and population — put the best half-weight panel size at zero**, and the
blended seed came back at **0.62× to 0.92× of the truth at one individual** across four
populations.

**⛦ An earlier version of that comment and of this paragraph said "every arm" put it at zero**,
which the project's own record contradicts: the sweep's *depth-crossed* arms put it at 0 on a
strong rare-allele pile-up and at **200** on a moderate one, and that two-hundred-fold
disagreement is the reason no single constant is right
(`ng_seed_shrinkage_2026-08-26.md` §5.2, `PROJECT_STATUS.md`). Corrected after review, in the code
comment and here.

### The example was kept here, and A5 deleted it

> **⛦ Superseded by A5 below.** `examples/ng_seed_shape_weight_sweep.rs` also reads
> `fit_spectrum_shape`, which A5 deletes, so keeping it was possible for exactly one commit. What
> follows is what this step did; A5's table says what became of it.

`examples/ng_seed_shape_weight_sweep.rs` reads `HALF_WEIGHT_PANEL_SIZE` in four places, so keeping
it meant giving it a local `const HALF_WEIGHT_PANEL_SIZE: f64 = 0.25;`. **That is a copy of a
deleted constant and it is labelled as one** — the doc comment on it says the library no longer
holds it and why. Its head now opens with what the sweep found and what became of the mechanism:
the panel's own fitted shape is at its best at *one* individual and degrades as the panel grows, so
the sweep answered its own question the wrong way round.

Three sentences in the program's own printed output claimed the library *ships* the constant. All
three were corrected — this is the "grep for the retired sentence" check applied to that one file.
**Applied to the whole tree it finds more, and a review found them**: see the Milestone A review
section at the end.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,871 passed, 0 failed, 14 ignored**, from 4,874. Three tests removed,
  none added: the expected fall.
- `cargo doc --no-deps --lib` — **27 unresolved links**, unchanged.

---

## A4 — `DiversityUnreachable` goes, and the reason is a proof

**Contract (plan A4).** A curve's own two moments always satisfy `E[2f(1−f)] ≤ 2 E[f](1 − E[f])` by
Jensen, so no total is ever out of reach on this route. Delete the variant and its test, and record
the inequality where the variant was. `ZeroDiversity` stays.

### What shipped

- **`SeedRegime::DiversityUnreachable` is gone**, and **a comment stands exactly where it was**,
  carrying what it did, why it existed — this is the failure the repeat-tract seed used to have and
  it must not return silently — and the argument that makes it unreachable now.
- **`PinnedTotal` is gone with it.** `total_for_diversity` returns an `f64` and holds the
  inequality as a **release assertion**. The enum existed only to carry the "no total reaches it"
  answer to a regime that no longer exists.
- **The argument, written out where the code is**: `θ = E[2f(1−f)] = 2E[f] − 2E[f²]` against a
  ceiling of `2E[f](1 − E[f]) = 2E[f] − 2E[f]²`, and `E[f²] ≥ E[f]²` — Jensen, whose slack is
  exactly the spread of the population's frequencies, over the whole density and not over the Beta
  alone. So the measurement sits below the ceiling by twice that spread, with equality only where
  the whole population sits at one frequency, which no density the fit can produce does.
- **How much room that leaves has a closed form, and it is now pinned by a test.** Where every
  position segregates and neither point mass carries anything, the density is a bare `Beta(a, b)`
  and the share of the ceiling is exactly `(a + b)/(a + b + 1)`, so the solved total is exactly
  `a + b`. The fit clamps `a` and `b` to `[0.02, 50]` independently, so the tightest case in the
  box is `a = b = 50`: **one part in 101 of the ceiling left unused, and a total of 100
  chromosomes.** Adding either point mass widens it — swept over the whole box, the largest share
  any density asks for is that same 100 in 101, and it is asked for only where the masses carry
  nothing (`seed_tests::no_density_the_fit_can_produce_comes_within_one_part_in_a_hundred_of_its_ceiling`).

  **⛦ An earlier version of this argument said the tightest density was `Beta(50, 50)` "the
  narrowest of those, still with a spread of 2.5 in 1,000".** Two things were wrong. `Beta(50, 50)`
  is not the narrowest the clamp admits — `a` and `b` are clamped independently, so `Beta(0.02, 50)`
  is reachable and its spread is 7.8 in a *million*, 316 times narrower. And absolute spread was
  the wrong quantity: what bounds the assertion is spread **relative to the ceiling**, which is
  where the closed form above comes from. The conclusion held; the reason given for it did not.
  Corrected after review, and replaced by a swept test rather than a second argument.

### Why it is a release assertion and not a `debug_assert!`

**Because without it the run dies three frames later naming the wrong thing.** Measured, by
downgrading the check and running the module in release:

| what a caller passed | what it panics with instead |
|---|---|
| frequency `1e-9`, heterozygosity `6e-4` | `the reference concentration must be finite and strictly positive, got -1.0000033323444377` |
| frequency `0.5`, heterozygosity `0.5` | `the reference concentration must be finite and strictly positive, got inf` |

Both come from `SpectrumSeed::new` in a different module. Neither says that the two numbers handed
over cannot both describe one population, which is the only thing a reader needs to know.

### Tests

Three went with the mechanism —
`a_fully_invariant_cohort_at_a_measured_diversity_falls_to_the_neutral_rung_and_says_so`,
`a_heterozygosity_of_exactly_a_half_is_not_refused`, and
`a_measurement_exactly_at_the_shapes_ceiling_has_no_total` — and four arrived:

- **`two_moments_that_cannot_belong_to_one_curve_are_refused`**, at the state a fully invariant
  panel used to produce: a mean frequency of 1 in a thousand million against a heterozygosity of 6
  in 10,000, five orders of magnitude past the ceiling.
- **`a_heterozygosity_exactly_at_the_ceiling_is_refused`**, at frequency and heterozygosity both a
  half — the `≥` boundary rather than the `>` one, which is where the solved total goes infinite.
- **`just_below_the_ceiling_still_has_a_total`**, at 999 parts in a thousand of the ceiling, where
  the total is near a thousand. Without it the two above would pass with the comparison written as
  `share_of_ceiling > 0.0`.
- **`one_bit_below_the_ceiling_still_has_a_total`**, keeping the 9.0 × 10¹⁵ figure the deleted test
  carried.

**The release-held assertion battery was run**: the new `assert!` downgraded to `debug_assert!`,
`cargo test --release --lib ng::calling::genotype_prior --all-features` — **2 tests failed and both
are the two `#[should_panic]` tests above**, so the check is reached and nothing else depends on it.
The source was restored and the restore checked with `git diff` before anything else ran.

**The two `#[should_panic]` strings are the same, and deliberately**: they exercise one check from
two sides. The module's other refusal — a heterozygosity above a half, which is caught before the
frequency is looked at — keeps its own distinct string, *"is not a thin estimate"*.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,872 passed, 0 failed, 14 ignored**, from 4,871. Three tests removed,
  four added.
- `cargo doc --no-deps --lib` — **27 unresolved links**, unchanged.

---

## A5 — the projection and the search are deleted

**Contract (plan A5).** Delete `FittedSpectrum`, `fit_spectrum_shape`, `fit_pair`,
`fill_expected_spectrum`, `SpectrumMatch`, `MAX_PROJECTION_INDIVIDUALS`,
`FrequencyDensity::allele_count_classes` and `FittedFrequencySpectrum` — **and not the
variable-census-site count**, which Milestone B re-sources. `examples/ng_inbreeding_sensitivity.rs`
and `ng_spectrum_panel_floor.rs` consume the deleted machinery; each is either retired with a note
saying what it measured, or kept against a local copy.

### What was deleted

`seed_generic.rs` falls from **3,440 lines to 958**. Gone: the two-branch spectrum prediction
(`fill_expected_spectrum`, `fill_expected_spectrum_at`, `log_branch_split`,
`MAX_PROJECTION_CONCENTRATION`, `BRANCH_TAIL_TOLERANCE`, `NEGLIGIBLE_BRANCH_WEIGHT`), the search
(`fit_pair`, `SpectrumScorer`, `ScoredPoint`, `sweep_from`, `sweep_once`, `line_search`,
`bounds_along`, `concentrations_at`, `at_search_limit`, `spectrum_log_likelihood`,
`spectrum_entropy`, `ProjectionFit`, `SEARCH_STARTS`, `SEARCH_DIRECTIONS`, both search ranges,
`MAX_PROJECTION_INDIVIDUALS`, `SPECTRUM_NORMALISATION_TOLERANCE`), the two wrapper types
(`FittedSpectrum`, `FittedShape`) and `fit_spectrum_shape`.

Elsewhere: `SpectrumMatch` from `genotype_prior/mod.rs`;
`FrequencyDensity::allele_count_classes` and `MAX_PROJECTED_PANEL` from
`parameter_estimation/joint/fit.rs`; `FittedFrequencySpectrum` from `calling/run_parameters.rs`.

**Three comments stand where the machinery was** — in `seed_generic`'s module header, at
`allele_count_classes`, and at `FittedFrequencySpectrum` — each saying what the thing did, what
consumed it, and where the numbers it produced now live.

### The debt A5 must not lose, recorded where the code was

Spec §5 marks **the count of census positions that came out variable across the panel** as the one
thing that must survive. Its only producer was `FittedFrequencySpectrum::of`, which computed it as
one minus the two end classes — so the code could not survive the deletion, and what survives is
the requirement. The note left in `run_parameters.rs` states it: the quantity is not the share that
segregates in the population (the two differed 6.6-fold at one individual on a tomato-like
density), spec §6.2 re-sources it from the fit's own per-position posteriors as a **soft count**,
and **nothing computes it today** — it is step C2.

### The oracle the tests needed, rebuilt

`implied_heterozygosity` — the oracle for the pin, and the reason
`the_seeds_implied_diversity_is_the_measured_one_at_every_shape` is a check rather than a
restatement — read the module's own spectrum machinery at one individual. That machinery is gone.
It is rewritten as `2 · B(1 + α_alt, 1 + α_ref) / B(α_alt, α_ref)` through `lgamma`: the same
Beta-binomial by the same route, four lines, and it still **shares no line of arithmetic** with the
pin's `t / (1 − t)`.

### The examples: four retired, two cut down

The plan named two programs. **Four more consume the machinery**, and the same rule was applied to
all six.

| program | what happened | why |
|---|---|---|
| `ng_spectrum_projection_cost.rs` | **deleted** | it measured what one spectrum prediction costs; there are no predictions |
| `ng_spectrum_panel_floor.rs` | **deleted** | it measured what the search's pair loses against the density it was fitted to |
| `ng_seed_shape_weight_sweep.rs` | **deleted** | it swept the blend's constant *through the search*; A3 kept it, and A5 makes that impossible |
| `ng_inbreeding_sensitivity.rs` | **cut to its live half** | one of its two routes was the search; the other — the diversity divided by `1 − F` — needs none of it and still runs |
| `ng_prior_moment_estimators.rs` | **one arm retired** | its "what the caller gets today" table ran the search |
| `ng_prior_moments_from_reads.rs` | **one arm retired** | same, and its `TodaysPath` became `FromTheCurve`: the two integrals and the seed they imply |

**Why deleted rather than kept against a local copy**, which the plan offers as the alternative:
the machinery is about 1,300 lines and its correctness rested on roughly fifty tests that go with
it. A copy in `examples/` cannot carry those tests — `cargo test` does not run an example's test
module — so what would be kept is an untested copy of deleted code, whose figures would then be
facts about the copy. Every retired program's numbers are already in the reports it was written
for, and each retirement note says which report.

**⚑ Three citations now point at files that do not exist**, all in documents this step does not
edit: `research/ordinary_site_prior_moments.md` line 49, `spec/ordinary_site_seed.md` line 93 and
`spec/population_diversity.md` line 459 cite `examples/ng_spectrum_panel_floor.rs`, and
`ng/reports/spectrum_projection_cost_2026-08-22.md` cites `ng_spectrum_projection_cost.rs`. The
first three are live specs, and the second and third of them already carry supersession banners
naming the moments spec. **Repointing them at the reports that hold the same figures is owed and
not done here**, because these are the design documents the implementation skill does not edit.

### Two `#[should_panic]` messages that moved

`a_regularizer_weight_that_is_not_a_count_of_sites_is_refused` and
`a_variable_site_count_that_is_not_a_count_is_refused` tested `FittedSpectrum::new` and went with
it, as did every other refusal that constructor held.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,822 passed, 0 failed, 11 ignored**, from 4,872 and 14. **Fifty tests and
  three ignored ones deleted with the code they pinned**; the plan expects this fall and it is not
  a regression. The three ignored were the search's own wall-clock measurements.
- `cargo doc --no-deps --lib` — **25 unresolved links, down from 27**. The two that went were
  pre-existing breaks inside the deleted machinery
  (`projection_tests::a_fit_costs_at_most_450_predictions` and
  `::the_cost_of_one_fit_by_panel_size`). One new break was introduced and fixed: a paragraph on
  `fill_locus_concentration`'s floor argued from the search box's bottom corner, and now argues
  from the identity instead.
- **Both surviving harnesses were run**, not merely compiled:
  `cargo run --release --example ng_prior_moment_estimators -- 500 2` and
  `cargo run --release --example ng_prior_moments_from_reads -- 400 0.15 0.0 1 0`, both to
  completion with the expected number of columns.
- **One stale sentence in a harness's printed output was found by running it**: the estimators
  program told its reader that the mean frequency column decides *whether the seed still needs its
  blend*. It does not still need one.

---

## B1 — the two estimators, over the fit's own posteriors

**Contract (plan B1).** A function in `parameter_estimation::joint` taking the converged
posteriors, the sample count and the panel's `F`, returning the two numbers. Types and the
reduction only; Milestone C wires it to a run.

### What shipped

A new module,
[`parameter_estimation/joint/census_moments.rs`](../../../../src/ng/parameter_estimation/joint/census_moments.rs),
holding `CensusMoments` — the two moments averaged over the census positions — and
`CensusMoments::from_posteriors(genotype_posterior, samples, positions)`.

Per position it forms the panel's expected alternative-copy count `E[k]` and the sum of the
samples' own posterior variances, then averages `k/2N` and `2k(2N − k)/(2N(2N − 1))` over
positions. **The fit's third number a sample — the posterior that the sample carries an extra copy
of the position — takes no part**: that is a mapping fact rather than an allele count, and the fit
scores it as its own class precisely so it need not be read as a heterozygote.

Nothing calls it. C1 is where it meets a run.

### One deviation: the inbreeding coefficient arrives at B3, not here

**The plan's contract lists `F` among B1's arguments and B3 is what applies it.** An argument that
nothing reads is a `clippy -D warnings` failure and, worse, a signature that promises a correction
the body does not make. So the function takes the posteriors and the two counts, and B3 adds `F`
along with the division that uses it. The end state is the plan's.

### What the heterozygosity owes as of this step, said in its own documentation

It substitutes `E[k]` into a formula quadratic in `k`, so it comes back **high by the variance** —
2.538 ± 0.165 times the truth at one sample and three reads a position, against 1.219 ± 0.152 with
the term (report §4.1) — and it applies no inbreeding correction, a further 80% at one individual
at `F = 0.8`. Both are named on the field itself, with the plan steps that close them.

### The fixtures, and what each is for

Nine tests. Two are the ones that could have been vacuous and are not:

- **`point_mass_posteriors_return_the_census_s_own_moments`** runs at **one individual and at a
  thousand**, which is spec §9's first test minus its inbreeding half. The two sizes are not a
  formality: at a thousand individuals writing `2N` where `2N − 1` belongs is a 0.05% error and
  sits inside any tolerance a test would set, and at one individual it is 50%.
- **`the_frequency_and_the_heterozygosity_are_not_the_same_number`** runs at **two** individuals
  with one alternative copy in four chromosomes, where the frequency is 1 in 4 and the
  heterozygosity 1 in 2. At one individual the two agree, so no single-individual fixture can tell
  a reduction that returned one for the other apart.

The remaining seven pin the ends (a panel fixed for either allele has no heterozygosity), the
carrier posterior's non-participation, the variance being summed per sample rather than squared as
a whole, a certain genotype having no variance, and the two refusals.

### Mutations run

Five, each against the module's nine tests:

| mutation | tests failing |
|---|---|
| the finite-panel correction dropped — `2N` for `2N − 1` | 2 of 9 |
| the `(2N − k)` factor dropped, so the heterozygosity becomes the frequency | 3 of 9 |
| a homozygous-alternative sample counted as one copy rather than two | 4 of 9 |
| `E[k²]` written `P(het) + 2·P(both)` instead of `+ 4·P(both)` | 2 of 9 |
| the carrier posterior read where the both-non-reference one belongs | 4 of 9 |

The source was restored from a backup after each and the restore checked with `git diff`.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings`, and a
  plain `cargo build --lib` — all clean.
- `cargo test --lib` — **4,831 passed, 0 failed, 11 ignored**, from 4,822. Nine tests added.
- `cargo doc --no-deps --lib` — **25 unresolved links**, unchanged.

---

## The Milestone A review, and what it changed

Three agents in worktrees at `02735054`, one brief each: arithmetic and numerics; tests and
mutation; design conformance and claim-checking. **Every finding below is the author's own claim
about the author's own work** — the figures quoted from the design documents and the reports came
back clean, which is the split this project keeps measuring.

### Wrong mechanism claims — two, and both were the reason given for a decision

**1. `Beta(50, 50)` is not the narrowest density the fit can produce.** The doc comment justifying
the release assertion that replaced `DiversityUnreachable` said the Beta's shapes are clamped to
`[0.02, 50]` "and the narrowest of those, `Beta(50, 50)`, still has a spread of 2.5 in 1,000". The
clamps are **independent**, so `Beta(0.02, 50)` is reachable and its spread is 7.8 in a million —
316 times narrower. Worse, absolute spread was the wrong quantity: what bounds the assertion is
spread *relative to the ceiling*. Replaced by a closed form — the bare Beta's share of its ceiling
is exactly `(a + b)/(a + b + 1)`, so the solved total is exactly `a + b` and the tightest case in
the box is 100 chromosomes — **and by a test that sweeps the box and finds 100 in 101**, so the
claim is now measured rather than argued. The conclusion never moved; the reason did.

**2. "every arm of the sweep put the best half-weight panel size at zero" is contradicted by this
project's own record.** All three arms of that sweep's *headline* table did. Its depth-crossed arms
put it at 0 on a strong rare-allele pile-up and at **200** on a moderate one — the two-hundred-fold
disagreement that is the actual reason no single constant is right. Corrected in the code comment
where the ramp was, here, and in `PROJECT_STATUS.md`.

### Wrong numbers — one

**"reading `b/(a+b)` for `a/(a+b)` returns a number too high on every one of them" is false on two
of the five.** `Beta(1, 1)` and `Beta(4, 4)` are symmetric, so the swap returns the *identical*
number there — invisible rather than consistent. The argument for adding `Beta(3, 0.6)` survives
and is now stated correctly. The same report's next paragraph already said those two are symmetric,
so the two paragraphs contradicted each other.

### A defect the deletion introduced

`total_for_diversity`'s rustdoc summary line was the **first line of the deleted `PinnedTotal`'s
doc comment**, orphaned when the enum went and rendering as a sentence that stops mid-clause.
Removed.

### An over-claimed oracle, and the assertion that closes it

`implied_heterozygosity` is `2 α_ref α_alt / (A(A + 1))`, **symmetric in the two concentrations**,
so it returns the same number for a pair and for its mirror. Its doc claimed only that it shares no
arithmetic with the pin, which is true and reads as more than it is. Two changes: the doc now says
what the oracle is blind to, and
`the_seeds_implied_diversity_is_the_measured_one_at_every_shape` now asserts the seed's own
expected frequency beside the heterozygosity. **Before that, the only assertion in the tree that
could see the pair swapped lived in another module.**

### Retired sentences the step's own grep missed

The A3 section above records finding three, all in one file. A tree-wide grep found more:

| where | what it still said |
|---|---|
| `examples/ng_prior_moment_estimators.rs`, twice | that the mixtures let the comparison against the current path go through the caller's own projection — both the arm and the projection are deleted |
| `src/ng/calling/mod.rs`, a doc comment and **a release panic message** | "a SNP/indel locus seeds from a frequency spectrum" |
| `doc/devel/ng/arch/calling_priors.md` §4 and its reuse map | the whole SNP/indel architecture section, stating `project_spectrum_seed(spectrum, diversity, panel_inbreeding)` and "399 predictions and 11.8 minutes at 3,200 individuals" as current |
| `doc/devel/ng/spec/ordinary_site_seed.md` §2, §6, §7 | the projection and the search as non-goals that "stay"; a seven-item checklist five of whose items rest on deleted machinery, including one asking that `DiversityUnreachable` **remain reachable**; two open questions about a constant that no longer exists |
| `doc/devel/ng/spec/population_diversity.md` §3.2, §8, §9 | the `FittedSpectrum` adapter as the current seam, and the class-weight projection as a live acceptance criterion |

The two source files were fixed outright. **The three design documents were given supersession
banners rather than rewritten** — a banner saying the code no longer matches is a factual note and
follows the owner's own precedent of 2026-08-27; rewriting the design is not this loop's to do.

**⚑ Still owed, and deliberately not touched:** `doc/devel/ng/impl_plan/calling_loop.md` step E2
still names `project_spectrum_seed` in its contract. That plan is being executed on a sibling
branch right now, and editing it here would collide.

### Two A5 changes the report did not mention, now named

- **`ng_prior_moments_from_reads.rs`'s "calls moved, trebled" column changed what it measures.** It
  was the control on the *search-versus-census* comparison and is now the control on the
  *curve-versus-census* one, because the first comparison is gone. Report §9.2's published trebling
  figures were computed the other way.
- **A5's verification line in the plan asks for the implied-heterozygosity pin "at every panel size
  and shape"; what shipped runs over six shapes and no panel sizes.** That is correct — there is no
  panel size left in the seed to vary, which is the point of the milestone — but it is a change to
  a stated acceptance criterion and belongs in the record.

### What came back clean

Every test count, line count and broken-link count in this report reconciles against the tree. The
A1 closed form and its six hand-computed values, the identity in `total_for_diversity` and its
boundary handling, the rebuilt oracle's algebra, A4's two measured panic messages, and every figure
quoted from the design documents and the measurement report — 34 of 36, 0.749×, 0.62×–0.92×,
1.22×, 399 predictions, 11.8 minutes, 6.6-fold — all verified correct. The arithmetic review also
swept the box the fit clamps to and found no input on which the new release assertion can fire from
one fitted curve, with about 1% of headroom.

---

## B2 — the variance term

**Contract (plan B2).** `E[k(2N − k)] = 2N·E[k] − E[k]² − Var(k)`. Own commit, do not bundle, and
the fixture is named in the spec: posteriors midway between genotypes at one sample, where dropping
the term returns 2.5× the truth. **A cohort test cannot catch this.**

### What shipped

`nei_heterozygosity` takes the position's copy-count variance and subtracts it. That is the whole
change to the arithmetic; the variance was already being formed per position at B1 and was not
being read.

### The oracle, and it is exact

**At one individual the heterozygosity is exactly the posterior that the individual is
heterozygous.** That is what the question means — an individual's two chromosomes differ at a
position exactly when it is heterozygous there — and the algebra collapses to it. With `h` the
heterozygous posterior and `d` the both-non-reference one, `E[k] = h + 2d` and
`Var(k) = h + 4d − (h + 2d)²`, so

```text
2 E[k] − E[k]² − Var(k)  =  2(h + 2d) − (h + 2d)² − h − 4d + (h + 2d)²  =  h
```

— every term but `h` cancels. **Drop `Var(k)` and it does not**: what is left is
`2(h + 2d) − (h + 2d)²`, a different number at every `d` above zero. Pinned over five posterior
shapes.

### The fixture the spec names, and the number it gives

Posteriors `(0.3, 0.4, 0.3)` over reference, heterozygous and both non-reference at **one sample**
— reads that have barely decided, which is the shape three reads a position produces.

| | value |
|---|---:|
| the truth (the heterozygous posterior) | 0.400 |
| substituting `E[k]` and stopping — `E[k]` is exactly 1 here, so `2·1 − 1²` | 1.000 |
| **the ratio** | **2.5** |

Against the **2.538 ± 0.165** the report measures through a whole fit at one sample and three
reads. The fixture is not that measurement — it is a hand-chosen posterior — and it lands within
1.5% of it.

### Why a cohort test cannot catch it, shown rather than asserted

The same posteriors at **63 samples**: `E[k] = 63`, `Var(k) = 63 × 0.6 = 37.8`, and the term moves
the heterozygosity from 0.5040 to 0.4992 — an inflation of **0.96%**. `Var(k)` grows with the panel
while `E[k]²` grows with its square, so the term's share falls like `1/N`; and a real 63-sample fit
is far more certain per sample than this fixture, which is where the report's *agree to three
decimals* comes from.

**That figure is pinned, not bounded.** A test asserting only "under 1%" is also satisfied by zero,
and zero is exactly what a deleted term gives.

**And a note against a false comfort**: the point-mass tests B1 added say nothing about this term.
A certain genotype has no variance, so `E[k(2N − k)]` and `2N·E[k] − E[k]²` are the same number
there — stated in a test of its own so nobody reads those as covering it.

### Mutation run

Setting the variance passed to `nei_heterozygosity` to zero — the exact defect this step
prevents — fails **2 tests** of the whole 4,835-test library suite, and both are the two written
here for it. The source was restored from a backup and the restore checked with `git diff`.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,835 passed, 0 failed, 11 ignored**, from 4,832. Three tests added.

### The tests-and-mutation review, and the three fixtures it closed

Run at `02735054`, forty-three tool uses, nine container runs, and it found **one surviving
mutation on shipped arithmetic and two fixture sets that could not fail.**

**1. `FrequencyDensity::p_segregating`'s clamp had no test.** Removing `.max(0.0)` from
`1 − p_invariant − p_fixed_alt` left the whole 4,822-test suite green. It is the shared input to
**both** integrals, so a density whose two masses total above one sends both moments negative at
once — and every fixture in the tree has the two masses totalling at most 0.991, so nothing sat
near the saturation point. A doc comment elsewhere already leaned on the clamp by name. Closed by
`two_point_masses_totalling_above_one_leave_nothing_segregating`; the mutation now fails 1 test.

**2. Every heterozygosity fixture had `a · b = 1`.** `Beta(0.5, 2.0)` is the density in
`the_expected_heterozygosity_is_the_densitys_own`, in `a_lopsided_density`, and in row 4 of the
shape list — and in all three the product in the numerator is exactly one. So **deleting `a * b`
from the formula altogether**, which removes the Beta's whole shape dependence, passed the test
written to pin that formula. Closed by adding `Beta(0.35, 1.20)` beside it, where `a · b = 0.42`;
the mutation now fails 5 tests where it failed 2.

**3. The six-shape list was decoration.** Its twenty-five lines of justification claim it spans
`a ≤ b` and `a > b` so that an `a`-for-`b` swap shows. Its only consumer asserted the seed's implied
heterozygosity against `density.expected_heterozygosity()` — **an identity with respect to the
moment functions**, since the same call supplied the seed's input. The frequency assertion added
earlier in this review round was an identity in the same way. Measured: the exact `a`↔`b` swap left
that test green.

Closed by making the list carry **both of each density's moments as hand-computed literals**, in a
type of its own, and adding `both_closed_forms_are_what_a_hand_calculation_gives` to compare them.
The shape test now seeds from the literals rather than from the accessors, so it is no longer an
identity either. The swap mutation fails 3 tests where it failed 2, and — the point — one of the
three is now the test whose doc comment claims to catch it.

**And one thing the list still cannot do, now said on the list itself**: the heterozygosity column
cannot see an `a`-for-`b` swap on any row, because `2ab/((a+b)(a+b+1))` is symmetric in its two
arguments. That is a fact about the quantity, not about the fixtures, and the swap is visible only
in the mean frequency.

### What the review found sound, and what it recorded as owed

The `should_panic` battery came back clean: all four release-held checks in `seed_generic`
demoted at once, and **all seven `should_panic` tests failed**, none satisfied by a different panic.
`fill_locus_concentration`, `total_for_diversity`'s solved total, and
`ExpectedAlternativeFrequency`'s constructor are each pinned by exact-value assertions.

Recorded as gaps rather than defects, because Milestone A is not the step that closes them:

- **Nothing in `src/` reads `SeedRegime`**, and the test that pinned two runs emitting different
  records went with the variant's payload fields. What a run reports is spec §7 and Milestone C's
  step C3.
- **`RunParameters::seed_from_moments` has no production caller.** Expected at Milestone A — C1 is
  where the seam meets a run — but it means two tests are the entire integration story.
- **No test pins an absolute seed value.** The change moves `α_ref` off 1.0 for the first time: on
  the fixture density the pair is `(0.2223, 3.711e-4)` where production's was `(1.0, 6.06e-4)`. A
  reader cannot see the size of that from the suite.

---

## B3 — the inbreeding correction

**Contract (plan B3).** Divide the heterozygosity by `1 − F/(2N − 1)`. Own commit, do not bundle.
The fixture is one individual at `F = 0.8`, where the factor is `1 − F` and its absence is an 80%
error; the companion is a panel of a thousand.

### What shipped

`CensusMoments::from_posteriors` takes the panel's inbreeding coefficient and divides the
heterozygosity by `inbreeding_factor(F, 2N) = 1 − F/(2N − 1)`. **The frequency is untouched**, and
that is not an omission: inbreeding rearranges copies between an individual's two chromosomes
without changing how many the panel holds, and the frequency is linear in that count.

**Why the factor is what it is, in one sentence**: a pair of chromosomes drawn at random from the
panel comes from the *same individual* with probability `1/(2N − 1)`, and with probability `F` such
a pair is one ancestral copy counted twice and cannot differ.

**It never divides by zero**, because `InbreedingF` admits `[0, 1)` and the factor is smallest at
one individual, where it is `1 − F`.

### The fixture the spec names, and its companion

| | one individual | a thousand individuals |
|---|---:|---:|
| the factor `1 − F/(2N − 1)` at `F = 0.8` | 0.200 | 0.99960 |
| what the panel shows, against the population | 20% of it | 99.96% of it |
| what putting it back lifts the answer by | ×5 | 4.004 parts in 10,000 |

**A test written at a thousand would pass with the correction deleted** at any tolerance loose
enough to survive a real census's own scatter. That is the whole reason the fixture runs at one
individual, and the companion is there to say so with a number rather than a claim.

**⛦ And a slip caught by the test itself, of exactly the kind the arithmetic review caught earlier
in this branch.** The first version asserted the lift at a thousand was `0.8/1999`. It is not: that
is the *shortfall* — what the panel is missing — and the lift that puts it back is
`(0.8/1999)/(1 − 0.8/1999)`, 4.004 parts in 10,000 against 4.002. The two are the same fact from
opposite sides and they are not the same number; the test now asserts both and a comment says which
is which.

### Mutations run

Three, each against the whole 4,840-test library suite:

| mutation | tests failing |
|---|---|
| the correction removed entirely | 3 |
| `2N` written where `2N − 1` belongs | 2 |
| the sign flipped — multiplying by `1 + F/(2N − 1)` | 3 |

The source was restored from a backup and the restore checked with `git diff`.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,840 passed, 0 failed, 11 ignored**, from 4,837. Three tests added.
- `cargo doc --no-deps --lib` — **25 unresolved links**, unchanged.

---

## C1 — where the moments are accumulated

**Contract (plan C1).** Spec §5 leaves the implementer a choice: accumulate inside the expectation
step's per-position loop with `genotype_posteriors` off, or turn the flag on and reduce afterwards.
**Whichever ships states its memory cost in the module's own documentation** — the flag is 12 bytes
a position a sample, 1.5 GB for fifty samples over two million positions.

### What shipped: the running sums, which is what the design wants

`CensusMomentSums` — four numbers, whatever the census and whatever the cohort — lives on the fit's
per-chunk `Statistics`, is fed one position at a time from the expectation step's own scratch
buffer, and is merged chunk into chunk the way every other sum there is. **Nothing is copied and
nothing is stored.** `JointFit` carries it.

**The memory the other route would have cost is stated on the field itself**: the stored array is
three `f32`s a sample a position, 1.5 GB at fifty samples over the shipped two-million-position
census, held only to be summed once.

**It is collected on every pass and only the last pass's value is read.** A sum costs less than the
branch that would skip it, and a pass that skipped it would leave the field describing whichever
earlier pass last filled it.

### Two decisions inside that, both about not diverging from the array route

- **A position whose likelihood underflowed counts, carrying nothing.** That path returns early and
  pushes zeros into the stored array so later positions are not attributed to their neighbours; the
  sums now record a position of no alternative copies there too. Otherwise the two routes would
  divide by different position counts.
- **The sums go on `JointFit`, not the finished moments.** The heterozygosity needs the panel's
  inbreeding coefficient, and §4.1 prefers one this route does not produce — the runs estimator
  walks genome windows in the per-sample histogram route and this one walks census positions.
  `CensusMomentSums::finish` takes the coefficient, so choosing where it comes from stays with
  whoever assembles the run. **That choice is step C4.**

### The test that makes the memory decision safe

**`the_summed_moments_and_the_stored_array_agree`** runs a real fit on a drawn cohort of six samples
over 2,000 positions with the flag on, and compares the running sums against
`CensusMoments::from_posteriors` over the stored array. **The stored array is what every
measurement behind this design was made on**, so if the two ever part company the sums are the
wrong ones.

**They agree to `f32`, which is where the array rounds and the sums do not.** Measured: **2.9e-11 on
the frequency and 8.2e-10 on the heterozygosity**, four to five orders inside the `1e-6` asserted.
The tolerance is deliberately loose against the measurement, because it bounds a rounding rather
than an arithmetic, and a figure that tight would fail on a different cohort. The test also pins
that both routes saw the same number of positions, and that both moments are above `1e-4` — so the
agreement is not two zeros matching.

### Mutations run

Two, on the chunk merge, each against the whole library suite: dropping the position count fails
1 test, dropping the frequency sum fails 1 — the agreement test in both cases, which is the only
test in the tree that runs a fit across more than one chunk.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,841 passed, 0 failed, 11 ignored**, from 4,840. One test added.
- `cargo doc --no-deps --lib` — **25 unresolved links**, unchanged.

---

## C2 — the segregating-position count and the spreads

**Contract (plan C2).** A **soft** count — `Σ P(the position segregates)` — not a count of
positions with a non-zero expected copy count, which spec §6.2 records an earlier draft getting
wrong and which returns 100% of positions. The spread of each moment travels with it and is
**labelled a floor, not an interval**, because linked positions make an independence assumption too
narrow by 3 to 16 times. **No floor is applied.**

### The soft count, and the trap it avoids

```text
P(the position segregates)  =  1  −  Π over samples of P(no alternative copy)
                                  −  Π over samples of P(both copies non-reference)
```

formed in the same loop over samples that `E[k]` is already made in, so the hot path does not gain
a pass.

**The test for it is the trap itself.** Five samples, each with a 1 in 100 posterior of being
heterozygous and no more — the shape a real census has. Every such position's expected
alternative-copy count is 0.05, **above zero**, so the hard version calls every one of them
segregating and a run over two million positions reports 100% segregating. The soft count asks what
the words mean: the panel is all-reference at `0.99⁵ = 0.951`, so the position segregates at
**0.049**. Over a hundred such positions the run reports **4.9** where the hard version reports 100.

**Two things it assumes, both stated where it is computed.** The samples are treated as independent
given the position's posteriors, and they are not — they are coupled through the frequency they
share, exactly as `Var(k)` is. Positive coupling makes both ends more likely than the products, so
the true probability of segregating is a little lower and this count runs a little **high**; **its
size has not been measured and nothing here claims one**. And the carrier posterior takes no part,
for the same reason it takes no part in `E[k]`.

### The spreads, and why they are floors

The plain standard error of the mean across positions, on both moments. **What makes it a floor is
that census positions are linked** — a run of homozygosity or a shared haplotype makes neighbours
carry the same evidence twice — so a spread computed as though they were independent counts that
evidence more than once and comes back too narrow, by a factor
`parameter_prepass_census_sites.md` §5 puts between **3 and 16**. The field names say `floor` so a
reader cannot take one for an interval.

The heterozygosity's spread carries the inbreeding correction and the frequency's does not, because
the correction scales the one and leaves the other alone — a spread that did not travel with its own
number would describe a different quantity. There is a test for that.

Below two positions the spread is zero rather than `NaN`: one position has nothing to disagree with.

### No floor is applied, and no threshold

**Nothing branches on the segregating count.** Spec §6.2 forbids picking a floor until it is
measured and the measurement needs a real census; the run reports the count and takes no action,
which is distinguishable in the output from a floor that never fires. Nothing here computes the gap
between the two heterozygosity estimates either — that is C3's to print and §7's fourth open
question to threshold.

### Mutations run

Two, on the soft count, against the module's twenty-one tests: dropping the all-alternative end
fails 1, and **replacing the soft count with the hard one the spec forbids fails 2**.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,847 passed, 0 failed, 11 ignored**, from 4,841. Six tests added.
- `cargo doc --no-deps --lib` — **25 unresolved links**, unchanged. One new break was introduced
  and fixed: a `#[cfg(test)]` function cannot be the target of a doc link from code that is always
  compiled.

---

## C3 — what the run reports

**Contract (plan C3).** The two measured moments; the fit's own `expected_heterozygosity` beside
them, because two routes to one quantity is a diagnostic; the segregating count; the spreads; and
where the inbreeding coefficient came from. **The last one carries a circularity the output must
not hide.**

### What shipped

`CensusMomentsReport`, with a `Display` that prints spec §7's whole list, and
`InbreedingSource` — three variants, because the three sources are not interchangeable:

- **`RunsOfHomozygosity { windows }`** — the source §4.1 prefers, on three reasons of which the
  first is decisive: it reads the *distribution* of heterozygosity along a genome and needs no
  population expectation, so nothing about it depends on the diversity this correction is
  computing. The window count travels with it because its estimator's own floor is 3,000, below
  which what it returns is its noise.
- **`JointFitHomozygoteExcess`** — **circular here, and the report says so unconditionally**. That
  excess is measured against a population expectation the same fit produced, and the correction it
  feeds divides a diversity by `1 − F`. `parameter_prepass_generic.md` §6.3 states the rule in as
  many words, and the warning quotes the mechanism rather than the rule.
- **`User`** — per sample or one value for the whole panel, including zero. Not a single-sample
  feature: a user who knows how their material was bred knows it whatever the cohort size.

### Two thresholds not applied, and both are on the type's own documentation

**Nothing branches on the segregating count** (spec §6.2's floor, unmeasured) **and nothing
branches on the gap between the two heterozygosity estimates** (§8's fourth open question). The
second has a test of its own: a report whose two routes disagree **by a factor of two** produces no
warning. **A threshold at a tenth would fire on good runs** — a converged, healthy fit already
shows the curve's number 10.7% above the census average's on one of three populations measured, and
that population is the one whose shape the curve can hold exactly.

### The one-sample warning, and why it stops at one sample

Where the coefficient is the fit's own homozygote excess **and the panel holds one sample**, the
report adds a second warning: that excess is 0.000 whatever the truth, so the coefficient is a
floor and the heterozygosity is a floor with it — and it says what to multiply by. **That the
warning stops at one sample is a measurement rather than a taste**: the fit's coefficient goes from
0.000 at one sample to 0.833 at two against a truth of 0.8, and stays within 0.03 of the truth from
three samples to sixty-three (report §3.5). Pinned by a test at one sample and at two.

### What the tests pin

Six. Beyond the two above: that both spreads print the words *a floor, not an interval* — checked
by counting two occurrences, so dropping the label on one of them fails; that the runs estimator
warns at 1,200 windows, does not at 8,004, and **does not at exactly 3,000**, which is the boundary;
and that the three sources print differently, which is the whole requirement §7 restates from
`calling_priors.md` §4.

**⛦ A test caught a wrong claim of mine while it was being written.** The fixture's four positions
carry 1, 0, 0 and 2 alternative copies in a panel of two identical samples, and I asserted that two
of them segregate. One does: the position where both samples carry two alternative copies is fixed
*for the alternative* and does not segregate any more than the all-reference ones do — which is the
same fact the soft count's own end-case test pins, read from the other side.

### Validation

- `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib` — **4,853 passed, 0 failed, 11 ignored**, from 4,847. Six tests added.
- `cargo doc --no-deps --lib` — **25 unresolved links**, unchanged.

---

## The merge with `main`, and what it brought back

`main` at `283c5a28` — the calling loop through step F2 — merged in after C3, at the owner's
request. Two conflicts, both where the two branches edited the same lines for unrelated reasons:
`run_parameters.rs`'s import block, and the `calling/mod.rs` field whose doc comment this branch had
corrected and which main replaced outright.

**⛦ It also brought back two copies of a retired sentence and one `compile_fail` doctest that now
passes for the wrong reason** — which is what a tree-wide grep is for after every merge, not only
after every deletion.

- *"a SNP/indel locus seeds from a frequency spectrum"*, in a **release panic message** in
  `calling/mod.rs` and in a test assertion in `summarise_condition.rs`. Both corrected.
- `stratum_fits.rs`'s doctest handed a `FittedSpectrum` where a `LengthSpectrum` belonged, to prove
  the two spectra cannot be crossed. **That type is deleted, so the test still fails to compile —
  for the wrong reason.** A `compile_fail` test that passes because the type it names is gone
  proves nothing about the type it was written to protect. Removed, with what it protected and why
  it has no subject any more recorded in its place: this caller fits one spectrum now, not two.

Library target 4,853 → **4,904**; `cargo test --tests` green across every integration target;
clippy, fmt and the 25 broken intra-doc links all unchanged.
