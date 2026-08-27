# ng — the ordinary-site prior's two numbers (implementation plan)

**Status:** draft, 2026-08-27. The build order for
[`../spec/ordinary_site_prior_moments.md`](../spec/ordinary_site_prior_moments.md), whose §2 splits
the work into two changes and says the first does not depend on the second. **This plan turns that
settled design into order; it is not a place for new design.** Where a step below meets a decision
the spec did not take, it says so and stops.

**Milestone A is step one** — the two numbers integrated off the fitted curve, and the deletion of
the projection, the search and the blend. It is a complete repair on its own: a run that takes only
Milestone A is better off than today whether or not B and C are ever built.

**Milestones B and C are step two** — the same two numbers averaged over the census positions, which
buys an error that goes to zero as the cohort grows and costs an expectation-step change, a variance
term and an inbreeding coefficient.

**The measurements are done.** [`../../reports/ng_ordinary_site_prior_moments_2026-08-27.md`](../../reports/ng_ordinary_site_prior_moments_2026-08-27.md)
answers the six research questions and the three harnesses behind it are in the tree. **No step
below is a sweep.** A step that finds itself wanting one has met a design gap, not a measurement gap.

---

## Scope

**In:** the closed-form moments on `FrequencyDensity`; the deletion of the class projection, the
two-parameter search, `FittedSpectrum`, `FittedFrequencySpectrum`, the blend and the inbreeding
coefficient at that seam; the census-average estimators on the joint fit, with their variance term
and their inbreeding correction; the segregating-position count and the spreads; and what the run
reports about all of it.

**Out:**

- **The runs-of-homozygosity estimator itself** — it exists
  (`parameter_estimation::generic::runs`) and this plan calls it; making it reachable from the joint
  route is spec §8's fifth open question and **Milestone C stops at the point that decision is
  needed.**
- **A floor on the segregating-position count.** Spec §6.2 forbids picking one until it is
  measured, and the measurement needs a real census.
- **A threshold on the two heterozygosity estimates' disagreement** — spec §7's fourth open
  question; the run prints both and judges neither.
- **The cohort gather's own emission.** The owner's ruling of 2026-08-27 is recorded at
  [`../spec/parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md) §4; building that
  step is its own plan.
- **The repeat-tract prior**, which takes no such detour.

---

## Principles (how the order was chosen)

- **The heart before the plumbing.** The two closed forms are three lines and the deletion is
  large; the arithmetic lands and is proven first, and only then does anything get removed.
- **A deletion is not a step to bundle.** Removing the search changes what every existing seed test
  is testing. It is its own commit so that a bisect can find it.
- **Isolate the steps whose failure is silent.** Three of them are: the expected-frequency formula,
  the variance term, and the inbreeding correction. Each returns a plausible number when it is
  wrong, each is marked below, and each names the fixture that catches it.
- **Types first, then implementation**, within every milestone (project rule).
- **Verify against ground truth.** The oracle is not self-consistency: it is the three harnesses
  already in the tree, and — for Milestone A — the fitted curve's own moments, which have closed
  forms.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

---

## Preconditions (already in place)

- **The spec is settled** and its §2 splits the work; §5 says where each piece lives.
- **The identity that turns two moments into a pair is built and shipped** —
  `total_for_diversity` in `seed_generic`, from [`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md)
  §3, which that document's supersession note leaves standing.
- **`FrequencyDensity::expected_heterozygosity` exists** and is the first of the two closed forms.
- **The three measurement harnesses are in the tree**: `examples/ng_prior_moment_estimators.rs`,
  `ng_prior_moments_from_reads.rs`, `ng_prior_moment_one_sample_inbreeding.rs`.
- **`parameter_estimation::generic::runs` exists** — Milestone C calls it and does not build it.
- **Not in place:** any consumer of the class projection other than the seed. `FittedSpectrum`'s
  only other intended consumer was the cohort gather, and the ruling above removed it.

---

## The steps

### Milestone A — step one: the moments off the curve, and the deletion

**A1. `FrequencyDensity::expected_alternative_frequency`.**  ✅
`p_fixed_alt + p_segregating · a/(a+b)`, beside `expected_heterozygosity`, with the same
population-not-panel framing in its documentation. **Own commit, do not bundle** — a wrong
expected frequency is a plausible number at every panel size and nothing downstream refuses it. The
oracle is the closed form checked against the search's own answer at **one individual**, where spec
§9 records the search is exact; the two agree to the search's 1% resolution on all five densities.
*Depends:* none. *Source:* spec §2, step one.

**A2. `project_spectrum_seed` takes the two numbers.**  ✅
Its signature becomes the two moments and nothing else — no spectrum, no panel size, no inbreeding
coefficient. The body is A1's frequency and the existing `total_for_diversity`. The three regimes
`SeedRegime` distinguishes stay what they are: a fitted curve, a fitted diversity with no curve,
neither. *Depends:* A1. *Source:* spec §6; [`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md) §3.

**A3. The blend goes.**  ✅
`HALF_WEIGHT_PANEL_SIZE`, `panel_shape_weight`, the log-space blend and `shape_from_panel` on
`SeedRegime::FittedSpectrum`. **The three tests that pin the ramp go with it**, and
`examples/ng_seed_shape_weight_sweep.rs` stays as the record of why — with its head rewritten to
say it measured a mechanism that is now deleted. *Depends:* A2. *Source:* spec §6.1.

**A4. `DiversityUnreachable` goes, and the reason is a proof.**  ✅
A curve's own two moments always satisfy `E[2f(1−f)] ≤ 2 E[f](1 − E[f])` by Jensen, so no total is
ever out of reach on this route. **Delete the variant and its test, and record the inequality where
the variant was** — a later reader must find why the refusal is absent rather than wonder.
`ZeroDiversity` stays: it is a real cohort state. *Depends:* A2. *Source:* spec §6, and
[`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md) §3.1.

**A5. The projection and the search are deleted.**  ✅
`FittedSpectrum`, `fit_spectrum_shape`, `fit_pair`, `fill_expected_spectrum`, `SpectrumMatch`,
`MAX_PROJECTION_INDIVIDUALS`, `FrequencyDensity::allele_count_classes`,
`FittedFrequencySpectrum` — **and not the variable-census-site count**, which spec §5 marks as the
one thing that must survive the deletion and which Milestone B re-sources. `examples/ng_inbreeding_sensitivity.rs`
and `ng_spectrum_panel_floor.rs` consume the deleted machinery; each is either retired with a note
saying what it measured, or kept against a local copy. **Own commit, do not bundle.**
*Depends:* A3, A4. *Source:* spec §5.

> **Checkpoint A:** the seed is two closed forms and an identity, the panel size appears nowhere in
> it, and the search is gone. **Milestone A is a complete repair — a run stopping here is better off
> than today.** Pause for review.

### Milestone B — step two's estimators, on known genotypes

**B1. The two estimators, over posteriors.**  ✅
A function in `parameter_estimation::joint` taking the converged posteriors, the sample count and
the panel's `F`, returning the two numbers. **Types and the reduction only** — Milestone C wires it
to a run. *Depends:* A5. *Source:* spec §3, §5.

**B2. The variance term.**  ✅
`E[k(2N − k)] = 2N·E[k] − E[k]² − Var(k)`. **Own commit, do not bundle, and the fixture is named
in the spec**: posteriors midway between genotypes at **one sample and three reads**, where
dropping the term returns 2.5× the truth. **A cohort test cannot catch this** — at 63 samples the
two agree to three decimals — so a test written on a panel is not a test of it.
*Depends:* B1. *Source:* spec §3.1, §9 test 2.

**B3. The inbreeding correction.**  ✅
Divide the heterozygosity by `1 − F/(2N − 1)`. **Own commit, do not bundle.** The fixture is one
individual at `F = 0.8`, where the factor is `1 − F` and its absence is an 80% error; the companion
is a panel of a thousand, where a `2N` written for `2N − 1` is 0.05% and only the single-individual
case would see it. *Depends:* B1. *Source:* spec §4, §9 test 1.

> **Checkpoint B:** both moments return the census's own exact values from point-mass posteriors, at
> one individual and at a thousand, and the two silent-failure terms each have a fixture that dies
> without them. Pause for review.

### Milestone C — step two wired to a run, and what it reports

**C1. Where the moments are accumulated.**  ✅
Spec §5 leaves the implementer a choice: accumulate inside the expectation step's per-position loop
with `genotype_posteriors` off, or turn the flag on and reduce afterwards. **Whichever ships states
its memory cost in the module's own documentation** — the flag is 12 bytes a position a sample, 1.5
GB for fifty samples over two million positions. *Depends:* B2, B3. *Source:* spec §5.

**C2. The segregating-position count and the spreads.**  ✅
A **soft** count — `Σ P(the position segregates)` — not a count of positions with a non-zero
expected copy count, which spec §6.2 records an earlier draft getting wrong and which returns 100%
of positions. The spread of each moment travels with it and is **labelled a floor, not an
interval**, because linked positions make an independence assumption too narrow by 3 to 16 times.
**No floor is applied** (Out, above). *Depends:* C1. *Source:* spec §6.2.

**C3. What the run reports.**  ☐
The two measured moments; the fit's own `expected_heterozygosity` beside them, because two routes to
one quantity is a diagnostic; the segregating count; the spreads; and where the inbreeding
coefficient came from. **The last one carries a circularity the output must not hide** — the joint
fit's homozygote excess is measured against a diversity the same fit produced. *Depends:* C2.
*Source:* spec §7.

**C4. Where the inbreeding coefficient comes from — stop and ask.**  ☐
Spec §4.1 prefers the runs estimator and §8's fifth open question records that **the joint route
does not produce one**: the runs estimator walks genome windows in the per-sample histogram route,
and this route walks census positions. Three shapes are possible and the spec does not choose.
**This step is a stop-and-ask, not an implementation** — bring the three options and a
recommendation, and do not pick one in the plan. *Depends:* C3. *Source:* spec §4.1, §8.5.

> **Checkpoint C:** a run estimates both moments from its own census, reports what they rest on, and
> says where its inbreeding coefficient came from. Pause for review — and C4 is a pause whether or
> not the checkpoint is honoured.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the closed forms against the search's own answer at one individual, where the search is exact on all five densities; the seed's implied heterozygosity equal to the measured one at every panel size and shape (spec §9 test 3); the three harnesses unchanged in what they report |
| B | point-mass posteriors returning the census's exact values at **one individual and at a thousand** (§9 test 1); midway posteriors at one sample and three reads returning less than the mean-substituted value, by the variance (§9 test 2) |
| C | `ng_prior_moments_from_reads.rs` reproducing the report's own figures against the wired-in estimators — a regression check on numbers that already exist, not a new sweep |

**And the harnesses stay.** A change that moves these numbers and leaves the tests green is a change
whose effect nobody has looked at.

---

## Out of scope (next plans)

- **The cohort gather's emission** — [`../spec/parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md) §4,
  under the owner's ruling of 2026-08-27; its own plan.
- **Reaching the runs estimator from the joint route** — spec §8's fifth question, unblocked by C4.
- **The thin-census floor** and **the two-estimate disagreement threshold** — spec §8's second and
  fourth, both needing a real census.
- **Confirmation on a real panel**, which stays open for everything in this area.
