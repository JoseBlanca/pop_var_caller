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

**The five densities this repository already sweeps all have `a ≤ b`** — the alternative allele
rare, or the two shapes balanced. On every one of them, reading `b/(a+b)` for `a/(a+b)` returns a
number too *high*, so a reader checking the sign of an error would see a consistent story and
conclude the formula was merely mis-scaled. `Beta(3, 0.6)` — the population where the reference
base is the rare one at the positions that vary (report §2) — was added for that reason, and the
swap moves the answer **down** there and up on the other five.

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

**A comment stands where the ramp was**, carrying what it did, the one constant it had, and the
two measurements that retired it: every arm of the sweep put the best half-weight panel size at
zero, and the blended seed came back at **0.62× to 0.92× of the truth at one individual** across
four populations.

### The example is kept, and it needed its own copy of the constant

`examples/ng_seed_shape_weight_sweep.rs` reads `HALF_WEIGHT_PANEL_SIZE` in four places, so keeping
it meant giving it a local `const HALF_WEIGHT_PANEL_SIZE: f64 = 0.25;`. **That is a copy of a
deleted constant and it is labelled as one** — the doc comment on it says the library no longer
holds it and why. Its head now opens with what the sweep found and what became of the mechanism:
the panel's own fitted shape is at its best at *one* individual and degrades as the panel grows, so
the sweep answered its own question the wrong way round.

Three sentences in the program's own printed output claimed the library *ships* the constant. All
three were corrected — this is the "grep for the retired sentence" check, and it found three.

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
  exactly the spread of the population's frequencies. So the measurement sits below the ceiling by
  twice that spread, with equality only where the whole population sits at one frequency. **No
  density the fit can produce does**: its Beta's shape parameters are clamped to `[0.02, 50]`, and
  the narrowest of those, `Beta(50, 50)`, still has a spread of 2.5 in 1,000.

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
