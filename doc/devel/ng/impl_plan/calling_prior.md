# ng genotype prior (step 8) — implementation plan

**Status:** draft, 2026-08-21. The build order for **step 8 whole**: the
`calling/genotype_prior/` module — the ported Dirichlet-multinomial primitive and its inbreeding
mixture, the per-sample (leave-one-out) concentration, the SNP/indel seed read off the fitted
spectrum, the STR seed with its measured total, and the plug-in comparator behind the same seam.
Design is settled in [`calling_priors.md`](../spec/calling_priors.md) (spec) and
[`../arch/calling_priors.md`](../arch/calling_priors.md) §2–§5 (types & interfaces). This plan
turns that design into build order; it is **not** a place for new design — the open questions
(Q1–Q6) are all carried by the spec, and none blocks a step below.

**Where this sits.** Six plans build calling:
[`calling_prerequisites`](calling_prerequisites.md) ∥
[`calling_foundations`](calling_foundations.md) → `calling_prior` ∥
[`calling_read_likelihoods`](calling_read_likelihoods.md) → [`calling_loop`](calling_loop.md) →
[`calling_bakeoffs`](calling_bakeoffs.md). **This plan needs nothing from the prerequisites
plan** — its upstream inputs all exist already: the fitted spectrum
(`FrequencyDensity`, [`joint/fit.rs:87`](../../../../src/ng/parameter_estimation/joint/fit.rs)),
θ (`expected_heterozygosity`, [`joint/fit.rs:109`](../../../../src/ng/parameter_estimation/joint/fit.rs)),
and the per-sample `F`. It needs only the foundations plan (the genotype table's flat views,
`AlleleId`). **That asymmetry is deliberate and worth acting on: this plan can branch the day
foundations merges, while the read-likelihoods plan waits for prerequisites — the two fan-out
plans start on different days and still run in parallel.**

---

## Scope

**In:** `src/ng/calling/genotype_prior/` — `mod.rs` (the `GenotypePriorModel` trait,
`Concentration`, `sample_concentration`), `dirichlet_multinomial.rs` (the ported primitive +
`MarginalizedDirichletPrior`), `seed_generic.rs` (`project_spectrum_seed`, `seed_for_locus`),
`seed_ssr.rs` (`ssr_seed`, `seed_length_distribution`), `hardy_weinberg.rs` (`PlugInWrightPrior`);
`types.rs` gains `ExpectedHeterozygosity` and `DEFAULT_SPECIES_DIVERSITY_FALLBACK`.

**Out (later plans or upstream):**

- **When the prior is called, and with which expected copies** — the loop's
  ([`calling_loop.md`](calling_loop.md)); this plan builds pure functions.
- **The GIAB single-sample 5× regression under both seam impls** — spec §12's definition of done
  needs genotypes, which need the loop; recorded there and in
  [`calling_bakeoffs.md`](calling_bakeoffs.md)'s out-of-scope as an unscheduled measurement.
- **`InbreedingF`'s `[0, 1)` tightening** — [`calling_prerequisites.md`](calling_prerequisites.md)
  Milestone A. This plan's mathematics is indifferent to it: the `F = 1` limit test drives the
  mixture on a raw value through a test-only path, not through the newtype (arch §2.1).
- **The rung-weight division between same-length spellings** (spec Q3) and **the
  `DiversityUnreachable` policy** (spec Q2) — open by design; the builder's signature and the
  outcome type isolate both, and the provisional behaviour (ceiling total + marker on the output)
  is what E1 builds.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The row function (Milestone B) is built and
  proven against two independent oracles before any seed exists to feed it; the seeds (D, E) are
  built against closed-form targets before anything composes them.
- **Reuse over rewrite.** The primitive, the floors, the mixture, the leave-one-out arithmetic
  and the STR shape are **ports** — arch §6's reconciliation table is the map, and no formula is
  re-derived. The projection is the one genuinely new function, and its optimiser reuses
  [`fitting/multistart.rs`](../../../../src/ng/parameter_estimation/fitting/multistart.rs).
- **Verify against ground truth.** The primitive against the rising-factorial oracle
  (`pochhammer_ln` / `dm_log_prior_oracle`,
  [`genetics.rs:240`](../../../../src/genetics.rs)); the mixture against the Wright biallelic
  formulas ([`genetics.rs:66`](../../../../src/genetics.rs)) in the concentrated limit; the
  projection against **exact expected spectra built in closed form — never `θ/k`, never
  sampled** (spec §12 tests 5–7 spell out why both shortcuts destroy the test).
- **Isolate the silent steps.** A biased projection, a wrong leave-one-out subtraction, or the
  STR total's units error each move genotypes without crashing; those land as their own commits,
  marked below.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- **Foundations merged:** `GenotypeTable` + `GenotypeTableView` + `homozygous_allele_for`,
  `AlleleId`, `CandidateAlleles` ([`calling_foundations.md`](calling_foundations.md) Milestones
  A–C). The prerequisites plan is **not** a precondition (see above).
- The production mathematics to port: `dirichlet_multinomial_log_priors`
  ([`genetics.rs:127`](../../../../src/genetics.rs)), `PROBABILITY_FLOOR` ([`:18`](../../../../src/genetics.rs)),
  `MIN_ALT_CONCENTRATION` ([`:187`](../../../../src/genetics.rs)), `alpha_from_diversity`
  ([`:214`](../../../../src/genetics.rs)), the mixture
  ([`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs)), the
  leave-one-out twins ([`em.rs:278`](../../../../src/ssr/cohort/em.rs),
  [`posterior_engine.rs:3183`](../../../../src/var_calling/posterior_engine.rs)), the STR shape
  (`g0_pseudocounts`, [`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs))
  and its fallback decay
  ([`param_estimation.rs:167`](../../../../src/ssr/cohort/param_estimation.rs)).
- The spectrum inputs: `FrequencyDensity` and `expected_heterozygosity` on the joint fit (above).
  **Impl-time confirmation, not design:** the concrete type `project_spectrum_seed` consumes is
  the pre-pass cohort gather's to pin (arch §6); confirm at D2 and adapt with a view, not a copy.

## Worktree, branch, merge

- **Worktree** `../pop_var_caller-calling-prior`, **branch** `ng-calling-prior`, from `main`
  **after `ng-calling-foundations` has merged** (nothing else need have).
- **Runs in parallel with** `ng-calling-read-likelihoods`. Conflict surface:
  `src/ng/calling/mod.rs` — each branch adds one `pub mod` line and its re-exports. Keep the
  additions on separate, alphabetically-placed lines. `types.rs` is this branch's alone in phase
  2 (`ExpectedHeterozygosity` appended to the population-genetics section; the likelihood plan
  adds nothing to `types.rs`).
- **Merge order back: this branch merges first** — it is the smaller surface — and
  `ng-calling-read-likelihoods` merges second, resolving any adjacent-line `mod.rs` conflict.
  The loop plan branches only after both are in.

---

## The steps

### Milestone A — scaffold + types (no logic)

**A1. Scaffold `genotype_prior/` and seed `types.rs`.**  ✅
`calling/genotype_prior/mod.rs` (declares the four files) wired into `calling/mod.rs`.
`types.rs` gains `ExpectedHeterozygosity` (constrained to `[0, 1]`, `try_new`/`get`; **the
cohort's expected heterozygosity at ordinary sites — not the non-reference rate**, and the doc
comment says so) and `DEFAULT_SPECIES_DIVERSITY_FALLBACK = 1e-3` (port of
`DEFAULT_DIVERSITY_PRIOR`, [`diversity.rs:78`](../../../../src/var_calling/diversity.rs), with
its "weakly informative, overridable, must be visible in output" reasoning). *Source:* arch §2.1;
spec §4.

**A2. The local types and the seam.**  ✅
In `mod.rs`: `Concentration` (borrow of caller scratch; invariant — every entry
`≥ MIN_ALT_CONCENTRATION`, length = allele count, checked in debug), `SeedRegime`
(`FittedSpectrum` / `NeutralShape` / `FallbackDiversity` — **a branch on absence, never on cohort
size**), `SpectrumSeed`, and the `GenotypePriorModel` trait taking the flat views (concentration,
`genotype_allele_counts`, `log_multinomial_coeffs`, `homozygous_allele_for`, `InbreedingF`,
`&mut [LogProb]`). No `Result` anywhere in the module — mis-shaped input is a caller bug →
assertion, structural ones held in release. *Source:* arch §2.2, §2.3, §3.2, §1.1.

> **Checkpoint A:** types compile; the trait's contract is documented (bit-identical rows at any
> thread count; no allocation). Pause for review.

### Milestone B — the row: primitive + mixture (the decision, pure)

**B1. Port the primitive.**  ✅
`dirichlet_multinomial.rs`: `fill_random_mating_log_priors` — `dirichlet_multinomial_log_priors` ported as-is from
[`genetics.rs:127`](../../../../src/genetics.rs) with one change — fill a caller slice instead of
returning `Vec` (the no-alloc contract, spec §8). Import `PROBABILITY_FLOOR` and
`MIN_ALT_CONCENTRATION` with their reasons. Carry the **independent parity oracle** across:
`pochhammer_ln` / `dm_log_prior_oracle` ([`genetics.rs:240`](../../../../src/genetics.rs)) — an
independent implementation, not golden values, so it keeps checking after constants move (spec
§12 test 4). Keep production's release-mode structural assertions (a short coefficient array
silently truncates the iteration otherwise). *Source:* spec §3.1, §8, §9; arch §6.

**B2. `MarginalizedDirichletPrior` — the two-branch inbreeding mixture.**  ✅
§3.2's mixture over the primitive: `logsumexp` on rows where `homozygous_allele_for` names an
allele, `log(1 − F)` alone elsewhere. **The homozygous test is the table's precomputed lookup,
consumed — never an inline comparison** (the one function the above-diploidy spec will change).
Tests, each pinning a spec §12 property: the **2:1 tripwire** (test 1 — het:hom-alt stays 2:1 at
`F = 0` across realistic θ; fails the moment anyone raises `α_ref`); invariant mass tracks θ
(test 2); the `F = 1` limit via a test-only raw-value path (test 3 — heterozygotes at the floor,
homozygotes at `α_ref : α_alt`); and the **Wright oracle** — at a concentration scaled to
dominance (`α × 10⁶` at fixed ratio), the row converges to
`wright_genotype_log_priors` ([`genetics.rs:66`](../../../../src/genetics.rs)) at `F = 0` and
`F = 0.5`, biallelic diploid. *Depends:* B1. *Source:* spec §3.2, §3.3, §12; arch §3.2.

> **Checkpoint B:** the row function is proven against two independent oracles and the four
> ported property tests. Pause for review.

### Milestone C — the per-sample concentration

**C1. `sample_concentration`.**  ✅ *(shipped as `fill_sample_concentration`, taking two checked
copy-count types rather than bare slices — both departures recorded in the step's report and owed
to arch §3.1)*
`α'_s(a) = seed(a) + max(0, cohort − own)` filling caller scratch — the port of
`leave_one_out_alpha` ([`em.rs:278`](../../../../src/ssr/cohort/em.rs); SNP twin
[`posterior_engine.rs:3183`](../../../../src/var_calling/posterior_engine.rs), identical
arithmetic). The `max(0, ·)` guards float noise only. Tests: **one sample ⇒ `out == seed` bit for
bit, no tolerance, no branch** (spec §12 test 8); **monotone in cohort evidence** — raising an
allele's cohort copies never lowers its weight for a sample that did not contribute the rise
(test 9). **Own commit, do not bundle** — a wrong subtraction double-counts a sample's own reads
into its prior and nothing crashes; tests 8 and 9 are the oracle, green before and after.
*Depends:* A2. *Source:* spec §6; arch §3.1.

> **Checkpoint C:** the leave-one-out builder is exact at both ends of the cohort range. Pause
> for review.

### Milestone D — the SNP/indel seed: the projection

**D1. The exact expected spectrum, in closed form.**  ✅ *(shipped as `fill_expected_spectrum`.
Its cost is per objective evaluation rather than per run — about a minute per fit at 800
individuals and hours by several thousand — which is D2's to answer; see the step's report.)*
The function that predicts a candidate `(α_ref, α_alt)`'s allele-count class probabilities at
`2N` chromosomes under **§3.2's two-branch sampling at the panel's `F`** — used twice: inside the
projection's objective, and to build the tests' targets. Closed form; nothing simulated.
Property test: the class probabilities sum to 1 across `N`, `θ`, `F` grids. *Depends:* B2 (the
two-branch sampling is the mixture's). *Source:* spec §4.1; arch §4.

**D2. `project_spectrum_seed`.**  ☐
Maximum-likelihood fit of D1's predicted class probabilities to the fitted spectrum's class
weights, over **all** classes including monomorphic, via
[`fitting/multistart.rs`](../../../../src/ng/parameter_estimation/fitting/multistart.rs) — a
two-parameter fit needs nothing new. `None` spectrum → `NeutralShape` at the fitted θ; no fitted
θ → `FallbackDiversity`. `SeedRegime` (with the regularizer weight and prior/data-dominated flag)
carried to the output. Confirm the concrete spectrum type here (`FrequencyDensity`, or the
cohort gather's wrapper — impl-time confirmation, arch §8). Tests, targets **built by D1 in
closed form**: a neutral spectrum projects to `(1, θ)` at several panel sizes (spec §12 test 5);
**invariance to inbreeding** — one density at `F = 0 / 0.6 / 0.9` returns one pair, where an
independent-chromosome projection returns `α_ref ≈ 0.91 / 0.86` (test 6, the test that holds the
two-branch requirement in place); at `n = 1` the pair is `(1, θ)` with no test of `n` — the only
branch is on the spectrum being absent (test 7). **Own commit, do not bundle** — a projection
biased by independent-chromosome sampling is 9–14% off on `α_ref` at tomato's `F` and nothing
crashes; tests 5–7 are the oracle. *Depends:* D1. *Source:* spec §4.1; arch §4.

**D3. `seed_for_locus`.**  ☐
Expand the run's pair onto one locus's table: `α_ref` first, the ALT total split evenly across
the locus's alternative alleles, floored at `MIN_ALT_CONCENTRATION` — the shape of
`alpha_from_diversity` ([`genetics.rs:214`](../../../../src/genetics.rs)) with the pair as input
instead of hard-coded. `VariantClass` stays an argument even while both classes pass one θ
(spec Q1: splitting later must not touch call sites). Test: a triallelic locus carries the same
total polymorphism as a biallelic one. *Depends:* D2. *Source:* spec §4; arch §4.

> **Checkpoint D:** the projection returns the neutral pair on neutral input, is invariant to
> `F`, and collapses correctly at one sample and on absence. Pause for review.

### Milestone E — the STR seed

**E1. `ssr_seed` — the shape ported, the total new.**  ☐
Geometric decay from the cohort's modal repeat count (shape of `g0_pseudocounts`,
[`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs), floored so a far
allele stays recoverable), scaled so the prior's own implied gene diversity equals the measured
`D`: `Σα = D / (1 − c − D)`, `c` the shape's Simpson index. `DEFAULT_G0_FALLBACK_DECAY = 0.5`
imported and **renamed for what it decays** (the genotype prior's pseudocount decay — not the
stutter one-step share; the sibling spec's §4.2 trap). Where `D ≥ 1 − c`, return
`SsrSeedOutcome::DiversityUnreachable { measured, ceiling }` — **reported, never silently
rescaled**; the provisional consumer behaviour (ceiling total, marker on the output) is the
loop's to wire and Q2's to settle. Tests: **the implied-diversity identity** —
`A(1 − c)/(A + 1)` recovers `D` to floating-point tolerance, whatever the decay and allele count
(spec §12 test 10 — *not* `Σα = D`, which is the units error this section exists to kill); the
refusal fires exactly at the bound (test 11). **Own commit, do not bundle** — the historical
defect here was a prior asserting two-fifths of the measurement, silently; test 10 is the oracle.
*Depends:* A2. *Source:* spec §5.1; arch §5.

**E2. `seed_length_distribution`.**  ☐
The seed shape normalised to a distribution over tract lengths — the one export the likelihood's
STR contamination stand-in composes, defined here so the prior's shape has one spelling. Test:
sums to 1; proportional to E1's weights. *Depends:* E1. *Source:* arch §5;
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §4.1.

> **Checkpoint E:** the STR seed reproduces the measured diversity or refuses aloud; the shared
> length distribution exists. Pause for review.

### Milestone F — the comparator

**F1. `PlugInWrightPrior`.**  ☐
Hardy–Weinberg at the plug-in frequency `α'_s(a)/Σα'_s` with the same `F` mixture, behind the
same trait — kept **only** for the spec's change measurements and the production differential,
never a shipping default. It runs on **the same seed** as the marginalized prior: a test pins
that no code path supplies production's `DEFAULT_REF_PSEUDOCOUNT = 10`
([`posterior_engine.rs:107`](../../../../src/var_calling/posterior_engine.rs)) — that constant is
the §2.3 trap, not a config. Test: at a concentrated input the two impls agree (Var(p) → 0,
spec §6); at one diploid sample's thin input they differ in the recorded direction (plug-in
under-counts homozygotes). *Depends:* B2, C1. *Source:* spec §2, §5.3; arch §3.2, §7.

> **Checkpoint F:** both seam impls run under one trait; the recipe can select either. **Step 8
> is complete as a set of pure functions.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | type-level tests; `ExpectedHeterozygosity` range |
| B | **two independent oracles** — the rising-factorial primitive oracle, the Wright concentrated-limit oracle — plus the 2:1 tripwire, θ-mass, and `F = 1` property tests |
| C | bit-equality of seed and leave-one-out at one sample; cohort-evidence monotonicity |
| D | **closed-form exact-spectrum targets** — neutral → `(1, θ)`, `F`-invariance, one-sample collapse — never `θ/k`, never sampled |
| E | the implied-diversity identity (`A(1 − c)/(A + 1) = D`); the refusal at the bound; the normalised export |
| F | comparator-vs-default agreement at concentration, divergence at thin input; the no-`α_ref = 10` pin |

## Out of scope (next plans)

- **Calling the prior per pass, the scratch that backs `Concentration`, the first-pass rule** —
  [`calling_loop.md`](calling_loop.md); the exact scratch slot sizes shared with `CallingScratch`
  are that plan's impl-time confirmation.
- **The GIAB single-sample 5× regression under both seam impls** — spec §12's definition of done;
  needs the loop, and its scheduling is recorded in
  [`calling_bakeoffs.md`](calling_bakeoffs.md)'s out-of-scope.
- **Q2's `DiversityUnreachable` policy and Q3's rung-weight division** — spec-owned open
  questions; E1's outcome type and counts-not-sequences signature isolate both.
- **A shared-prior fast path for large cohorts** — the loop document's deferred perf item
  (spec §10).
