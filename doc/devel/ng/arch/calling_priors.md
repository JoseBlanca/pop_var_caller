# ng — the genotype prior: types & interfaces

*Status: architecture draft (2026-08-21), companion to the spec
[`../spec/calling_priors.md`](../spec/calling_priors.md) (the design and its rationale) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step traits) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). It is one of three coordinated calling
arch docs — the siblings are [`read_likelihoods.md`](read_likelihoods.md) and
[`calling_em_loop.md`](calling_em_loop.md) — and §0 says which of the three owns each shared type.
`ng_step_interfaces.md` §3 sketches step 8 as `genotype_log_prior(&self, genotype, freq, f)`; that
sketch is superseded here (§3.3). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types, verbs
for functions, newtypes for domain scalars, **STR** in prose ↔ `ssr` in code. Signatures are
illustrative; the **contract** is the deliverable. See the spec for the why behind every decision.*

## 0. Who owns what — the three calling docs

One owner per shared type; the other two documents point here or at the sibling, never restate.

| shared thing | owner |
|---|---|
| `Concentration` semantics, the two seed builders, the per-sample (leave-one-out) builder | **this doc** (§2, §3) |
| `ExpectedHeterozygosity` scalar (the spec's `θ`) — seeds `types.rs` | **this doc** (§2.1) |
| `CandidateAlleles` (the allele table), `GenotypeTable` + `AlleleId` (the genotypes and their indexing), `ExpectedAlleleCopies`, `CallingScratch`, `Phred` | [`calling_em_loop.md`](calling_em_loop.md) §2 |
| the per-sample evidence views, the `Lg` row contract, the stutter distribution and its parameter context | [`read_likelihoods.md`](read_likelihoods.md) §2, §4 |
| `LogProb`, `InbreedingF`, `Ploidy`, `ErrorRate` | exist in [`src/ng/types.rs`](../../../../src/ng/types.rs) — reconciled in §6, not re-minted |

The prior takes **flat slices**, never the loop's types: the leave-one-out builder takes `&[f64]`
copy counts, the row function takes the genotype table's flat views (§3). **The reason is not the
module tree** — the shared vocabulary lives in `calling/mod.rs` and every sub-module imports
downward, so nothing forbids the import. It is the no-allocation contract (spec §8) plus the reason
production's own primitive takes flat arrays: it avoids a back-reference into its caller
([`genetics.rs:127`](../../../../src/genetics.rs)).

## Module home

`src/ng/calling/genotype_prior/` — step 8, inside the `calling/` folder that holds steps 6–9
([`module_layout.md`](module_layout.md); the four are one folder because keeping them apart forced a
no-import rule between them). A folder of its own because the step has a real bake-off
(marginalized against plug-in, spec §2):

```
src/ng/calling/genotype_prior/
├── mod.rs                   – GenotypePriorModel trait, Concentration, the per-sample builder
├── dirichlet_multinomial.rs – the ported log-prior primitive + MarginalizedDirichletPrior
├── hardy_weinberg.rs        – PlugInWrightPrior, the comparator impl
├── seed_generic.rs          – SNP/indel starting point: the spectrum projection (spec §4.1)
└── seed_ssr.rs              – STR starting point: geometric shape, total from gene diversity (spec §5.1)
```

New scalar newtypes seed `src/ng/types.rs` beside the ones already there (its line 447 already
reserves the diversity scalar's slot by name).

## 1. General — what every path shares

The prior is two functions with one contract between them:

1. **Build a concentration** — per sample, per pass: the run's seed plus what the *other* samples
   showed at this locus (spec §6). Cheap arithmetic, no `lgamma`.
2. **Turn a concentration into a log-prior row** — one `LogProb` per candidate genotype, up to a
   shared additive constant (softmax-ready). This is where the marginalized/plug-in fork lives, and
   it is the expensive call the EM loop multiplies by the cohort size
   ([`calling_em_loop.md`](calling_em_loop.md) §3).

Both are pure functions of their arguments (spec §1.1 goal 5); neither allocates — the caller
hands in scratch (spec §8; the buffers live in the loop's `CallingScratch`,
[`calling_em_loop.md`](calling_em_loop.md) §2).

### 1.1 Errors, numerics, determinism

Contract only; the reasoning is spec §8. A mis-shaped input (concentration length ≠ allele count,
row length ≠ genotype count, a non-positive concentration entry) is a **caller bug → assertion**,
the structural ones held in release exactly as production's primitive holds them
([`genetics.rs:127`](../../../../src/genetics.rs)) — no `Result` anywhere in this module.
Probabilities are floored before any logarithm (`PROBABILITY_FLOOR`,
[`genetics.rs:18`](../../../../src/genetics.rs), imported); every concentration entry is floored
strictly positive (`MIN_ALT_CONCENTRATION`, [`genetics.rs:187`](../../../../src/genetics.rs)). No
RNG, no clock; the one order-sensitive sum (cohort expected copies) is the loop's to fix, not
this module's.

## 2. Types

### 2.1 Scalars — seed `types.rs`

```rust
/// The cohort's expected heterozygosity at ordinary sites — the spec's θ (spec §4).
/// NOT the non-reference rate. Constrained to [0, 1]; checked constructor per the
/// types.rs convention. Source: the pre-pass (JointFit::expected_heterozygosity, or the
/// histogram route's mean); fallback DEFAULT_SPECIES_DIVERSITY_FALLBACK below.
pub struct ExpectedHeterozygosity(f64);   // try_new / get

impl ExpectedHeterozygosity {
    /// Species-range fallback for a run with no fitted diversity, ~human nucleotide
    /// diversity. Weakly informative, overridable; a run on it must say so in its output
    /// (spec §4). Port of production's DEFAULT_DIVERSITY_PRIOR
    /// (src/var_calling/diversity.rs:78).
    pub const SPECIES_FALLBACK: Self = Self(1e-3);   // unitless heterozygosity
}
```

**A value of the type, not a free `f64` — decided at implementation, 2026-08-21.** As a loose
float the constant constructs an `ErrorRate`, an `InbreedingF` and a `GenotypeFrequency` just as
happily, because `1e-3` is a legal value of each — which is the confusion the four separate
`[0, 1]` scalars in `types.rs` exist to prevent, and it would be the only *diversity* constant in
the shared vocabulary for an STR-path author to reach for. Precedent in the same file:
`AlleleId::REFERENCE`. Report:
[`ng_calling_prior_a1_fixes_2026-08-21.md`](../../reports/implementations/ng_calling_prior_a1_fixes_2026-08-21.md).

`InbreedingF` already exists ([`types.rs:388`](../../../../src/ng/types.rs)) **but its constructor
accepts `1.0` and the spec requires `[0, 1)`** — the ceiling is a property of the type, not of one
estimator (spec §7). Tighten the constructor; the `F = 1` mathematical-limit test (spec §12 test 3)
drives the mixture function on a raw value through a test-only path, not through the newtype.

**Three things the tightening touches, and the third is the one that bites.**

- **It cannot be done in `checked_probability`** ([`types.rs:326`](../../../../src/ng/types.rs)),
  which is shared with the other fraction newtypes and is right to admit `1.0` for them. Give
  `InbreedingF` its own half-open check and a `DomainError` variant that says so.
- **An existing test asserts the current behaviour** —
  [`types.rs:862`](../../../../src/ng/types.rs) requires `try_new(1.0)` to succeed. It moves to the
  rejection list beside `1.5`.
- **The pre-pass builds one from a fitted value with `.expect(…)`**
  ([`runs.rs:634`](../../../../src/ng/parameter_estimation/generic/runs.rs)) — a
  coverage-weighted posterior occupancy, which can in principle reach exactly `1.0` on a fully
  homozygous sample. **Tightening the newtype alone turns a legitimate fit into a panic.** The fitted
  path must clamp below the ceiling before constructing, not assert; a caller cannot recover from a
  panic and the estimator's own §7 reasoning is that `F = 1` is unreachable *in practice*, which is
  not the same as unreachable in arithmetic.

### 2.2 The concentration

```rust
/// One strictly-positive number per allele of the locus's table, parallel to
/// CandidateAlleles — read as "chromosomes the prior behaves as though it had already
/// seen" (spec §1). Never constructed free-standing: the builders below fill a
/// caller-owned buffer, so nothing allocates per sample per pass (spec §8).
/// Invariant: every entry ≥ MIN_ALT_CONCENTRATION; length == allele count.
pub struct Concentration<'a>(&'a [f64]);   // borrow of scratch, checked in debug
```

**It, `PriorRow` (§3.2) and `SpectrumSeed` (§2.3) live in a private `mod checked` inside
`genotype_prior/mod.rs`, re-exported — decided at implementation, 2026-08-21, and the nesting is
load-bearing.** A private field is visible to a module's *descendants*, and the four sub-modules
here are descendants: with these types declared directly in `genotype_prior`, a struct literal in
`dirichlet_multinomial.rs` builds one field by field and skips the constructor. Measured — it
compiled and ran. One level of nesting makes those four siblings instead, and the literal fails
with `error[E0451]`.

### 2.3 The seed — per run (SNP/indel) or per locus (STR), with its provenance

```rust
/// How the run's starting concentration was obtained — must reach the run's output,
/// because a run on the fallback and a run on a fitted spectrum are otherwise
/// indistinguishable (spec §4, §4.1).
pub enum SeedRegime {
    /// Read off the pre-pass's fitted spectrum by the §4.1 projection.
    FittedSpectrum { regularizer_site_weight: f64, data_dominated: bool },
    /// No spectrum emitted (absent below the panel-size floor, or one sample):
    /// the neutral pair (1, θ) at the fitted θ. A branch on ABSENCE, never on cohort
    /// size (spec §4.1).
    NeutralShape,
    /// No fitted θ at all: the species-range fallback.
    FallbackDiversity,
}

/// The SNP/indel seed: two numbers for the whole run (per variant class — Q1 keeps the
/// class argument even while both classes pass the same θ, spec §4.2).
pub struct SpectrumSeed {
    pub alpha_ref: f64,            // 1.0 on a neutral panel — the fit's landing point, not a knob
    pub alpha_alt_total: f64,      // θ on a neutral panel; shared across the ALT alleles per locus
    pub regime: SeedRegime,
}
```

## 3. Interfaces

### 3.1 The per-sample concentration (general — both paths, identical)

```rust
/// α'_s(a) = seed(a) + max(0, cohort expected copies of a − this sample's own).
/// The max(0,·) guards float noise only (spec §6). Fills `out`; allocates nothing.
/// At one sample the two slices are equal and out == seed, bit for bit — no branch
/// (spec §6; pinned by test 8 of spec §12).
pub fn sample_concentration(
    seed: &[f64],
    cohort_expected_copies: &[f64],   // the loop's ExpectedAlleleCopies, as a slice
    own_expected_copies: &[f64],
    out: &mut [f64],
);
```

### 3.2 The row — the step-8 seam (the marginalized/plug-in bake-off)

```rust
/// One sample's log-prior over every candidate genotype, up to a shared additive
/// constant. Flat views (counts, coefficients, homozygous lookup) come from the
/// loop's GenotypeTable (calling_em_loop.md §2); taking them flat is what keeps this
/// module free of a back-reference into inference/ (spec §9).
/// The six buffers one row call reads and writes, with every shape check already run —
/// `new` is the only constructor and it runs them, so no implementation is reachable
/// with mis-matched buffers. Borrow-only: nothing allocates.
pub struct PriorRow<'a> { /* private */ }

impl<'a> PriorRow<'a> {
    pub fn new(
        concentration: Concentration<'a>,           // α'_s, parallel to the allele table
        genotype_allele_counts: &'a [u32],          // n_genotypes × n_alleles, row-major
        log_multinomial_coeffs: &'a [f64],          // ln C(ploidy; counts), per genotype
        homozygous_allele_for: &'a [Option<AlleleId>], // the ONE homozygous test (spec §3.3)
        per_allele_scratch: &'a mut [f64],          // working space, one entry per allele
        out: &'a mut [LogProb],
    ) -> Self;
    // + by-value accessors for the four views, and scratch_and_out(&mut self)
}

pub trait GenotypePriorModel {
    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF);
}

/// Default: §3.1's Dirichlet-multinomial with §3.2's two-branch inbreeding mixture
/// (logsumexp on the homozygous rows). Ports genetics.rs:127 + the engine mixture.
pub struct MarginalizedDirichletPrior;

/// Comparator: Hardy–Weinberg at the plug-in frequency α'_s(a)/Σα'_s, same F mixture.
/// Kept ONLY for the spec's change measurements (§2.2's GIAB regression, §12) and the
/// production differential — never a shipping default.
pub struct PlugInWrightPrior;
```

**A checked bundle rather than eight flat parameters — decided at implementation, 2026-08-21,
owner-authorised 2026-08-22.** The sketch this replaces passed the six buffers directly and put
the shape checks in a helper the trait's prose asked every implementation to call. Three
independent review agents reached the same conclusion: those are one defect, because a trait
cannot require that a method body opens with a particular call — measured, deleting the call left
every test passing. **The contract is unchanged** and that is what §7's decision is about: the
same six caller-owned buffers, nothing allocated, the same flat views, no back-reference into the
loop. **The per-allele scratch is the one genuinely new argument** — spec §8 requires scratch
"sized by allele count and genotype count", and without it the ported primitive allocates a `Vec`
of `lgamma(α_a)` per call. At a diploid six-allele locus (21 genotypes, 36 non-zero allele slots)
keeping it costs 42 `lgamma` calls per sample per pass against 72.

**Contract.** Same inputs → bit-identical rows at any thread count. The homozygous branch reads
`homozygous_allele_for` and nothing else decides homozygosity — that lookup is the “one function”
the above-diploidy spec will later change (spec §3.3). Cost: one `lgamma` per (allele, non-zero
count) pair plus one `logsumexp` per homozygous genotype (spec §8).

### 3.3 Correction to `ng_step_interfaces.md` §3 step 8 — recorded, not applied there

The sketch `genotype_log_prior(&self, genotype: &Genotype, freq: &[AlleleFreq], f) -> LogProb` is
superseded twice over. **Input:** the prior is built from a *concentration* — expected allele
copies with the sample's own contribution subtracted, added onto the seed — not from a frequency
vector; a frequency vector cannot carry the conviction (`Σα`) that §2.2 of the spec shows is the
whole repair (why: spec §6, §2). **Output:** one genotype at a time re-derives the per-allele
`lgamma` baseline per call; the primitive returns a whole row from flat arrays and ng ports that
shape (why: spec §3.1, [`genetics.rs:127`](../../../../src/genetics.rs)). An `AlleleFreq` newtype
is therefore not minted — nothing in the three calling docs consumes one.

## 4. The SNP/indel path

The seed comes from the pre-pass's fitted spectrum, projected onto `(α_ref, α_alt)`:

```rust
/// Project the fitted spectrum onto the two-parameter family — maximum-likelihood fit
/// of predicted class probabilities to the fitted spectrum's class weights, over ALL
/// classes including monomorphic, predicting with §3.2's two-branch sampling at the
/// panel's F (independent chromosomes bias α_ref down 9–14% at tomato's F; spec §4.1).
/// A change of representation, not a second estimate. `None` spectrum → NeutralShape.
pub fn project_spectrum_seed(
    spectrum: Option<&FittedSpectrum>,   // absent below the pre-pass's panel-size floor
    diversity: ExpectedHeterozygosity,
    panel_inbreeding: InbreedingF,
    chromosomes: u32,                    // 2N of the panel the spectrum was fitted at
    class: VariantClass,                 // Substitution | InsertionOrDeletion — Q1's argument
) -> SpectrumSeed;

/// Expand the run's two numbers to one locus's table: α_ref first, the ALT total split
/// evenly across the locus's alternative alleles, floored (spec §4). Port of
/// alpha_from_diversity (genetics.rs:214) with the pair as input instead of hard-coded.
pub fn seed_for_locus(seed: &SpectrumSeed, allele_count: usize, out: &mut [f64]);
```

The projection's optimiser reuses the pre-pass's fitting machinery
([`fitting/multistart.rs`](../../../../src/ng/parameter_estimation/fitting/multistart.rs)) — a
two-parameter fit needs nothing new. Census-site exclusion on depth, the regularizer sweep, and
the per-class reporting are the **pre-pass's** obligations (spec §4.1's traps); this module only
carries `SeedRegime` through to the output.

## 5. The STR path

The projection does not reach here (spec §5); the seed is per locus, from three inputs the STR
side already has:

```rust
/// The STR seed: geometric decay from the cohort's modal repeat count (the shape,
/// production's G₀), scaled so the prior's own implied gene diversity equals the
/// measured D:  Σα = D / (1 − c − D),  c the shape's Simpson index (spec §5.1).
pub fn ssr_seed(
    candidate_repeat_counts: &[u32],     // parallel to the locus's CandidateAlleles
    modal_repeat_count: u32,             // the cohort's mode at this locus
    decay: f64,                          // fitted per group of loci; fallback DEFAULT_G0_FALLBACK_DECAY
    gene_diversity: f64,                 // D, the pre-pass's STR diversity — never the SNP θ
    out: &mut [f64],
) -> SsrSeedOutcome;

pub enum SsrSeedOutcome {
    Seeded,
    /// D ≥ 1 − c: no total reproduces the measurement — the geometry cannot hold it.
    /// REPORTED, never silently rescaled (spec §5.1, test 11). Until Q2 settles the
    /// policy, the loop uses the ceiling total and this marker travels onto the locus's
    /// output through the provenance channel (read_likelihoods.md §8).
    DiversityUnreachable { measured: f64, ceiling: f64 },
}

/// Coded fallback decay for a period with no fitted value — the genotype prior's
/// pseudocount decay, NOT the stutter one-step share (spec sibling §4.2's trap).
/// Source: production's DEFAULT_G0_FALLBACK_P (ssr/cohort/param_estimation.rs:167).
pub const DEFAULT_G0_FALLBACK_DECAY: f64 = 0.5;   // unitless, per repeat unit of offset
```

**Two alleles on one rung** (spec §5.2): v1 gives each same-length spelling the rung's full
weight — production's behaviour — because the division needs the interrupted-repeat work to say
how to weight it. `OPEN:` spec Q3; the builder takes the counts, not the sequences, precisely so
the change lands in one function.

**One export the likelihood composes** (its §4.5.1 contamination stand-in): the seed shape
normalised to a distribution over tract lengths —
`pub fn seed_length_distribution(candidate_repeat_counts, modal, decay, out: &mut [f64])`.
Computed once per locus by the loop and handed into the STR scoring context
([`read_likelihoods.md`](read_likelihoods.md) §4); defined here so the prior's shape has one
spelling.

## 6. Reconciliation with existing code

Every row read on 2026-08-21.

| ng name | existing code | action |
|---|---|---|
| `fill_random_mating_log_priors` (the port) | [`genetics.rs:127`](../../../../src/genetics.rs) | **port as-is**, one change: fill a caller slice, and put the per-allele `lgamma(α_a)` baseline in caller scratch, instead of returning `Vec` and allocating a second one (spec §8 no-alloc). Renamed from `dirichlet_multinomial_log_priors` at implementation: the file carries the distribution, the function carries which half of §3.2's mixture the values are, and it fills rather than returns. **Production has a second, private copy of the same mathematics** — `fill_log_indep_per_g_from` ([`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs)), which the SNP/indel engine runs and which already takes a caller's `out` and `lgamma_alpha`. It associates differently and disagrees with the shared primitive on 112 of 492 measured genotype values, by at most one unit in the last place; ng's port is bit-identical to the shared one. |
| `PROBABILITY_FLOOR`, `MIN_ALT_CONCENTRATION` | [`genetics.rs:18`](../../../../src/genetics.rs), [`genetics.rs:187`](../../../../src/genetics.rs) | import with their reasoning |
| `seed_for_locus` | `alpha_from_diversity`, [`genetics.rs:214`](../../../../src/genetics.rs); `ALPHA_REF` [`:179`](../../../../src/genetics.rs) | **shape ported, source not** — the pair comes from the projection; `(1, θ)` is where a neutral panel lands |
| the inbreeding mixture | [`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs) (`fill_log_prior_per_g_homogeneous`), STR port [`em.rs:290`](../../../../src/ssr/cohort/em.rs) | port as §3.2's two-branch form, inside `MarginalizedDirichletPrior` |
| `sample_concentration` | `leave_one_out_alpha`, [`em.rs:278`](../../../../src/ssr/cohort/em.rs); SNP twin at [`posterior_engine.rs:3183`](../../../../src/var_calling/posterior_engine.rs) | port (identical arithmetic in both) |
| `ssr_seed` shape | `g0_pseudocounts`, [`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs) | **shape ported, total mass new** (spec §5.1) |
| `DEFAULT_G0_FALLBACK_DECAY` | `DEFAULT_G0_FALLBACK_P`, [`param_estimation.rs:167`](../../../../src/ssr/cohort/param_estimation.rs) | import, renamed for what it decays |
| `PlugInWrightPrior`'s pseudocounts — what it must NOT inherit | `DEFAULT_REF_PSEUDOCOUNT = 10`, [`posterior_engine.rs:107`](../../../../src/var_calling/posterior_engine.rs) | the comparator runs on the same seed as the marginalized prior; `α_ref = 10` is the §2.3 trap, not a config |
| `project_spectrum_seed` | — (production never fitted a spectrum) | **new**; optimiser reuses [`fitting/multistart.rs`](../../../../src/ng/parameter_estimation/fitting/multistart.rs) |
| `FittedSpectrum` input | joint route's `FrequencyDensity`, [`joint/fit.rs:87`](../../../../src/ng/parameter_estimation/joint/fit.rs); `expected_heterozygosity` on `JointFit`, [`joint/fit.rs:199`](../../../../src/ng/parameter_estimation/joint/fit.rs) | consume; the concrete spectrum type is the pre-pass cohort-gather's to pin (impl-time confirmation) |
| `DEFAULT_SPECIES_DIVERSITY_FALLBACK` | `DEFAULT_DIVERSITY_PRIOR`, [`diversity.rs:78`](../../../../src/var_calling/diversity.rs) | import value + reasoning; carried as overridable, regime-reported |
| `InbreedingF` | [`types.rs:388`](../../../../src/ng/types.rs) | **reuse, tighten range to `[0, 1)`** (today it admits `1.0`; spec §7) |
| `SeedRegime` reporting | `Provenance` / `Estimate<T>`, [`parameter_estimation/mod.rs:60`](../../../../src/ng/parameter_estimation/mod.rs) | same idea, prior-specific variants; do not force-fit the four-variant enum |
| Wright biallelic formulas | [`genetics.rs:66`](../../../../src/genetics.rs) | **test oracle only** (spec §3.2) — plus the row basis of `PlugInWrightPrior` |
| independent parity oracle | `pochhammer_ln` / `dm_log_prior_oracle`, [`genetics.rs:240`](../../../../src/genetics.rs) | carry the test across (spec §9) |

## 7. Design decisions — decided

- **Marginalized is the default and plug-in is a comparator behind the same trait — decided.**
  One seam (`GenotypePriorModel`), two impls; the recipe selects. Production ships plug-in default
  behind an env toggle ([`driver.rs:287`](../../../../src/ssr/cohort/driver.rs)); ng inverts that,
  and the comparator exists so spec §2.2's measurement stays re-runnable. Why: spec §2, §5.3.
- **The prior takes flat slices, not the loop's types — decided.** Nothing allocates per sample per
  pass, and the primitive keeps no back-reference into its caller — production's own reason. **Not**
  a module-tree constraint: under `calling/` the shared types sit one level up and importing them
  would be legal. Why: spec §8, §9.
- **The homozygous test is the table's precomputed lookup, consumed — never an inline comparison —
  decided.** One place for the above-diploidy spec to change. Why: spec §3.3.
- **Seed regime and STR seed refusal travel to the output — decided.** Two runs on different
  information must be distinguishable. Why: spec §4, §4.1, §5.1.
- **`VariantClass` is an argument today even though both classes pass one θ — decided.** Splitting
  later must not touch call sites. Why: spec Q1 (“confirm before code”).
- **No `AlleleFreq` newtype — decided.** Nothing in the three calling docs consumes a frequency
  vector; minting one would re-open the §2 plug-in shape. Why: §3.3 above. **Note for whoever reads
  `types.rs` first:** its own comment names `AlleleFreq` among the constrained types still to come
  ([`types.rs:447`](../../../../src/ng/types.rs)). That was written before this seam existed and is
  an example, not an obligation — calling does not owe it, and a later step that genuinely needs a
  frequency should mint it for that use rather than to satisfy the comment.

## 8. Open items

- `OPEN:` rung-weight division between same-length spellings (spec Q3) — builder signature already
  isolates it.
- `OPEN:` the policy at `DiversityUnreachable` (spec Q2) — the outcome type is settled, the
  consumer's choice among the three candidate answers is not; provisional behaviour in §5.
- Impl-time confirmations: the concrete `FittedSpectrum` type from the pre-pass cohort gather.
  **The scratch slot the row function needs is settled: one `f64` per allele**, which is what the
  ported primitive hoists (`lgamma(α_a)`, one per allele); `CallingScratch` owns it.

## Test & bench shape

Unit tests beside each file, pinning spec §12's eleven properties — notably the 2:1 ratio tripwire
(test 1), the exact-expected-spectrum projection targets built in closed form, never `θ/k` and never
sampled (tests 5–7), bit-equality of seed and leave-one-out at one sample (test 8), the STR seed's
implied-diversity identity and the refusal (tests 10–11), and the pochhammer oracle (test 4). The
end-to-end anchor is the GIAB single-sample 5× regression (genotype accuracy at true variants; count
of true hom-alt called het) run under both impls of the seam — spec §12's definition of done.
