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
    /// `data_dominated` is the PANEL-WIDE comparison, and spec §4.1 is explicit that a
    /// panel-wide ratio is the wrong number to quote as reassurance — the tail is where
    /// the regularizer binds. The per-class ratio is the pre-pass's to emit beside its
    /// spectrum (§4); this carries the aggregate and claims nothing more.
    /// `spectrum_match` says HOW FAR the pair is from what was measured, and whether the
    /// search ran out of range before it got there — decided by the owner at Checkpoint D,
    /// 2026-08-22, because a run on a compromised starting point and one that matched were
    /// otherwise identical in the output.
    FittedSpectrum {
        regularizer_site_weight: f64,
        data_dominated: bool,
        spectrum_match: SpectrumMatch,
    },
    /// No spectrum emitted (absent below the panel-size floor, or one sample):
    /// the neutral pair (1, θ) at the fitted θ. A branch on ABSENCE, never on cohort
    /// size (spec §4.1).
    NeutralShape,
    /// No fitted θ at all: the species-range fallback.
    FallbackDiversity,
}

/// How far the fit's pair is from the measured spectrum. REPORTED, never returned as though it
/// had matched — the rule spec §12 test 11 sets for the STR seed, applied here.
pub struct SpectrumMatch { divergence_nats: f64, at_search_limit: bool }   // both by accessor

impl SpectrumMatch {
    /// The Kullback-Leibler divergence of the measurement from the fitted pair's prediction,
    /// in nats. ZERO means the family reproduced the measurement exactly. Free: the fit's
    /// objective is already the measurement's own entropy minus this, so it is the winning
    /// score subtracted from that entropy and costs no prediction.
    fn divergence_nats(self) -> f64;
    /// The pair sits on the edge of the range searched, so a better one may lie outside it.
    /// Carried separately because it is not derivable from the divergence — a pair pinned
    /// against a bound can still predict the measurement well, and a fully invariant cohort
    /// reaches this legitimately.
    fn at_search_limit(self) -> bool;
}

/// The SNP/indel seed: two numbers for the whole run. When Q1 splits the estimate by
/// variant class this carries both alternative totals and `fill_locus_concentration`
/// picks between them, so the run still holds one seed (spec §4.2, Q1).
pub struct SpectrumSeed {
    pub alpha_ref: f64,            // 1.0 on a neutral panel — the fit's landing point, not a knob
    pub alpha_alt_total: f64,      // θ on a neutral panel; shared across the ALT alleles per locus
    pub regime: SeedRegime,
}
```

## 3. Interfaces

### 3.1 The per-sample concentration (general — both paths, identical)

```rust
/// The two copy-count arrays, each checked when it is built. They borrow; they own
/// nothing, so wrapping the loop's buffers costs no allocation.
pub struct CohortAlleleCopies<'a>(/* private */);
pub struct SampleAlleleCopies<'a>(/* private */);

/// α'_s(a) = seed(a) + max(0, cohort expected copies of a − this sample's own).
/// The max(0,·) guards float noise only (spec §6). Fills `out`; allocates nothing.
/// At one sample the two arrays are equal and out == seed, bit for bit — no branch
/// (spec §6; pinned by test 8 of spec §12).
pub fn fill_sample_concentration(
    seed: &[f64],
    cohort_copies: CohortAlleleCopies<'_>,   // wraps the loop's ExpectedAlleleCopies
    own_copies: SampleAlleleCopies<'_>,
    out: &mut [f64],
);
```

**Two changes from the sketch, decided at implementation and owner-authorised 2026-08-22.**

**The two copy-count arguments are checked types, not bare slices**, for the reason §3.2's checked
bundle exists: the compiler could not otherwise tell them apart. They are the same shape and the
same unit, their difference is the whole of the leave-one-out term, and **passed the wrong way
round the function returns the bare seed at every allele** — the cohort's evidence gone, nothing
raised in release, and the debug guard blind to it at one sample where the two arrays are equal by
construction. Measured on the flat-slice version; with the types it is `error[E0308]`. The types
also carry a debug check that the entries are finite and non-negative, which the flat version had
nowhere: `f64::max` returns the other operand on a `NaN`, so a `NaN` copy count silently became a
zero difference and left the allele carrying nothing but its seed.

**`fill_sample_concentration`, not `sample_concentration`** — it fills a caller buffer and returns
nothing, as the module's three other buffer-fillers do, and *sample* reads as a verb in a
statistics context.

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
    /// Which prior this is, for a run to record beside the genotypes it produced — the seam
    /// exists to compare two, and an arm without a label is a number nobody can act on.
    /// Added at F1. NOT `Debug`: B2 recorded that deriving `Debug` would do in the meantime,
    /// and it never worked — the trait has no `Debug` supertrait, so `Box<dyn
    /// GenotypePriorModel>` does not implement it. Not the seed's provenance either:
    /// `SeedRegime` and `SpectrumMatch` describe the input the two impls SHARE.
    fn name(&self) -> &'static str;
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
/// panel's F. Independent chromosomes bias α_ref down 8.6% at F = 0.6, 12.1% at 0.8 and
/// 14.0% at 0.9, measured on 26 individuals at tomato's diversity (spec §4.1).
/// A change of representation, not a second estimate. `None` spectrum → NeutralShape.
pub fn project_spectrum_seed(
    spectrum: Option<FittedSpectrum<'_>>, // absent below the pre-pass's panel-size floor
    diversity: Option<ExpectedHeterozygosity>, // None → FallbackDiversity; the regime is
                                         // derived here rather than asserted by a caller
    panel_inbreeding: InbreedingF,
) -> SpectrumSeed;                       // no VariantClass: the split lands one step later

/// Expand the run's two numbers over one locus's alleles: α_ref first, the ALT total
/// split evenly across the locus's alternative alleles, floored (spec §4). Port of
/// alpha_from_diversity (genetics.rs:214) with the pair as input instead of hard-coded.
/// Named for what it does, like the module's other two fillers.
pub fn fill_locus_concentration(
    seed: SpectrumSeed,                  // Copy, 24 bytes — by value, like its accessors
    class: VariantClass,                 // Q1's argument here too; see below
    allele_count: usize,                 // checked against out.len(): `out` is scratch the
    out: &mut [f64],                     // loop slices, so its length is a slicing decision
);
```

**The panel size is not an argument to the projection**: `2N + 1` class weights already fix `N`,
and a second argument is a second place for it to disagree. `FittedSpectrum` is a **borrowed view
over the class weights** plus the regularizer's weight and the variable-site count — decided at
implementation, 2026-08-22, and the answer to §8's open item. Neither `FrequencyDensity`
([`joint/fit.rs`](../../../../src/ng/parameter_estimation/joint/fit.rs)) nor a gather wrapper is
it: the first is a density over the *population's* allele frequency, and the projection matches a
*panel's* allele counts.

**The projection's optimiser is not `fit_by_multistart`, and that is a measurement.** That driver
scores one cell at a time through `NoiseModel::append_genotype_likelihoods`, which takes `&self`
and so cannot cache; the natural cell here is one allele-count class, and every class would
rebuild the whole spectrum — 6,401 predictions where one is needed at 3,200 individuals, about 1.7
hours per candidate against 0.96 seconds. It also has no notion of a search direction that is not
an axis, and this surface is a ridge whose direction depends on the panel size, so the search
sweeps three: the total concentration, `α_ref` alone and `α_alt` alone. What is reused is the
driver's *shape* and its `SearchPrecision`. A fit is 399 predictions and 11.8 minutes at 3,200
individuals, once per run.

**Which end owns a SNP-versus-indel split — settled by the owner, 2026-08-23: the per-locus
expansion.** `project_spectrum_seed` reads the *shape* of variation off the panel's allele counts,
which the pre-pass fits without separating the two classes; a class-specific *scale* belongs where
the run's total is shared out over a locus's alleles. So the class argument sits on
`fill_locus_concentration` alone, and when the pre-pass measures two diversities `SpectrumSeed`
carries both totals — the run still holds one seed, which is what the calling loop's frozen
parameters assume.

**Still `OPEN:` (spec Q3's sibling under Q1): a locus carrying one alternative of each kind.** One
class per locus cannot express it. The classes are derivable from `CandidateAlleles`, which hold
the bases, so the change lands in one function when Q1 is settled.

Census-site exclusion on depth, the regularizer sweep, and the per-class reporting are the
**pre-pass's** obligations (spec §4.1's traps); this module only carries `SeedRegime` through to
the output.

## 5. The STR path

The projection does not reach here (spec §5); the seed is per locus, from three inputs the STR
side already has:

```rust
/// The STR seed: geometric decay from the cohort's modal repeat count (the shape,
/// production's G₀), scaled so the prior's own implied gene diversity equals the
/// measured D:  Σα = D / (1 − c − D),  c the shape's Simpson index (spec §5.1).
pub fn fill_ssr_seed<'a>(
    candidate_repeat_counts: &[u32],     // parallel to the locus's CandidateAlleles
    modal_repeat_count: u32,             // the cohort's mode at this locus
    decay: SeedDecayPerRepeat,           // fitted per group of loci; fallback ::FALLBACK
    gene_diversity: RepeatGeneDiversity, // the pre-pass's STR diversity — never the SNP θ
    out: &'a mut [f64],
) -> SsrSeedOutcome<'a>;

#[must_use]
pub enum SsrSeedOutcome<'a> {
    Seeded(Concentration<'a>),
    /// D ≥ 1 − c: no total reproduces the measurement — the geometry cannot hold it.
    /// REPORTED, never silently rescaled (spec §5.1, test 11). No Concentration comes
    /// back; what does is the buffer holding the normalised shape, so a caller with a
    /// policy for these loci scales it in place and wraps the result itself. Until Q2
    /// settles that policy the loop's provisional one is a ceiling total, and this
    /// marker travels onto the locus's output through the provenance channel
    /// (read_likelihoods.md §1.4).
    DiversityUnreachable { measured: f64, ceiling: f64, shape: &'a mut [f64] },
}
```

**Two departures from the sketch above, both shipped at E1 and both recorded here rather than in
the code alone.**

- **The two scalars are checked types in `ng/types.rs`, not bare `f64`** — `RepeatGeneDiversity`
  and `SeedDecayPerRepeat`, beside `ExpectedHeterozygosity`, with `try_new` returning
  `DomainError` like every other measured scalar the caller consumes. They are the pre-pass's
  outputs, so a degenerate fit returning a `NaN` is a run to refuse with a message rather than a
  process to abort — which is why they are *not* in `genotype_prior`'s `checked` module with its
  panicking constructors. The coded fallback decay is `SeedDecayPerRepeat::FALLBACK`, an
  associated constant following `ExpectedHeterozygosity::SPECIES_FALLBACK`, rather than a
  free-standing `DEFAULT_G0_FALLBACK_DECAY`: as a loose `f64` it is exactly as constructible into
  a stutter one-step share as into this, which is the trap the rename was for. Its doc carries
  that trap verbatim.
**A distance and not a verdict — revised 2026-08-23, after review.** The first version of
`SpectrumMatch` was an enum whose `Reproduced` variant claimed something it never checked: it was
set whenever the search finished inside its range and no allele-count class came back at exactly
zero, neither of which measures how close the answer is. Measured on a panel of 26 individuals
whose alleles sit at two middling frequencies, the fitted pair and the measurement shared **4
parts in 100** of their mass and the marker said they matched. It now reports the distance, in
nats, and names no threshold: nobody has measured how far off the pair has to be before a genotype
moves, so classifying was the part that had to go rather than the checking. Reference values from
the module's own tests: **1.1e-9 nats** where the family can hold the shape, **0.481** and
**3.153** on two it cannot. The `F = 1` case that the old `Unreproducible` variant caught now
appears as a divergence above 10 nats, because the objective charges the impossible classes
`ln(PROBABILITY_FLOOR)`.

- **The refusal withholds the concentration where `SpectrumMatch` marks a returned value**, and
  the difference is deliberate. The spectrum fit runs **once per run** and the run cannot start
  without a seed, so withholding would leave the caller nothing to do but invent one; a marker on
  a returned value is the right shape there. The STR refusal is **per locus** and the loop can
  pick a policy for that locus and carry on, so withholding costs nothing and stops a caller
  falling into a buffer it was never handed. Neither shape makes the mistake unrepresentable —
  `Concentration::new` is public and has to be, because the provisional ceiling-total policy needs
  it.

**One case the rule above does not cover: a locus with one candidate length.** Its shape has a
Simpson index of exactly 1 and therefore a ceiling of 0, so `D ≥ 1 − c` would refuse every
monomorphic tract whatever the measurement, including a measurement of zero. E1 seeds it at
`ALPHA_REF` instead: one length is one genotype, whose prior probability is 1 at any positive
concentration, so no total can be wrong there.

**Two alleles on one rung** (spec §5.2): v1 gives each same-length spelling the rung's full
weight — production's behaviour — because the division needs the interrupted-repeat work to say
how to weight it. `OPEN:` spec Q3; the builder takes the counts, not the sequences, precisely so
the change lands in one function.

**One export the likelihood composes** (its §4.5.1 contamination stand-in): the seed's shape
before it is scaled, normalised to sum to one —

```rust
pub fn fill_seed_share_per_candidate(
    candidate_repeat_counts: &[u32], modal_repeat_count: u32,
    decay: SeedDecayPerRepeat, out: &mut [f64],
);
```

Computed once per locus by the loop and handed into the STR scoring context
([`read_likelihoods.md`](read_likelihoods.md) §4.1); defined here so the prior's shape has one
implementation behind both consumers.

**`OPEN:` it is per candidate, and the term it feeds is per observed length.** This sketch called
it `seed_length_distribution`, and E2 shipped it under a name that says what the buffer actually
holds, because the two are not the same thing wherever a locus has more candidates than lengths.
The mixture's third term is `c · seed(o)` with `o` an observation, and three cases separate the
two supports:

- **two candidates of one length each take the rung's full share** — deliberate as a
  concentration and open as spec Q3, but read as a claim about lengths it double-counts:
  measured at `0.8` for the modal length against the geometry's own `0.667`, on a tract with the
  mode spelled twice and one length above it, at the fallback decay;
- **the candidate set is post-prune**, while the mixture's sibling uniform term is spread over
  every length the stutter model can reach from a candidate — a strictly larger support;
- **a censored read carries no length at all**, only a lower bound.

None of the three is the prior's to settle, and the likelihood step should meet them as a
decision rather than at the point of use. `SsrContamination::length_distribution`
([`read_likelihoods.md`](read_likelihoods.md) §4.1) is the field that has to say which support it
means.

**A caller that has just built the seed already holds this shape** and should not rebuild it: on
a seed it is the concentration divided by its own total, and on a refusal it is exactly the
buffer `SsrSeedOutcome::DiversityUnreachable` hands back. The export is for a caller that wants
the shape without the seed.

**It does not survive a candidate being added.** A discovery round appends candidates mid-locus
([`calling_em_loop.md`](calling_em_loop.md) §5), and a frozen candidate-parallel buffer is then
one entry short with nothing to raise, because by construction it is not refilled. Discovery is
off by default; a loop that turns it on has to rebuild this.

## 6. Reconciliation with existing code

Every row read on 2026-08-21.

| ng name | existing code | action |
|---|---|---|
| `fill_random_mating_log_priors` (the port) | [`genetics.rs:127`](../../../../src/genetics.rs) | **port as-is**, one change: fill a caller slice, and put the per-allele `lgamma(α_a)` baseline in caller scratch, instead of returning `Vec` and allocating a second one (spec §8 no-alloc). Renamed from `dirichlet_multinomial_log_priors` at implementation: the file carries the distribution, the function carries which half of §3.2's mixture the values are, and it fills rather than returns. **Production has a second, private copy of the same mathematics** — `fill_log_indep_per_g_from` ([`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs)), which the SNP/indel engine runs and which already takes a caller's `out` and `lgamma_alpha`. It associates differently and disagrees with the shared primitive on 112 of 492 measured genotype values, by at most one unit in the last place; ng's port is bit-identical to the shared one. |
| `PROBABILITY_FLOOR`, `MIN_ALT_CONCENTRATION` | [`genetics.rs:18`](../../../../src/genetics.rs), [`genetics.rs:187`](../../../../src/genetics.rs) | import with their reasoning |
| `seed_for_locus` | `alpha_from_diversity`, [`genetics.rs:214`](../../../../src/genetics.rs); `ALPHA_REF` [`:179`](../../../../src/genetics.rs) | **shape ported, source not** — the pair comes from the projection; `(1, θ)` is where a neutral panel lands |
| the inbreeding mixture | [`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs) (`fill_log_prior_per_g_homogeneous`), STR port [`em.rs:290`](../../../../src/ssr/cohort/em.rs) | port as §3.2's two-branch form, inside `MarginalizedDirichletPrior` |
| `sample_concentration` | `leave_one_out_alpha`, [`em.rs:278`](../../../../src/ssr/cohort/em.rs); SNP twin at [`posterior_engine.rs:3183`](../../../../src/var_calling/posterior_engine.rs) | port (identical arithmetic in both) |
| `fill_ssr_seed`'s shape (`fill_seed_shape`) | `g0_pseudocounts`, [`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs) | **shape ported, total mass new** (spec §5.1) |
| `SeedDecayPerRepeat::FALLBACK` | `DEFAULT_G0_FALLBACK_P`, [`param_estimation.rs:167`](../../../../src/ssr/cohort/param_estimation.rs) | import, **retyped** and renamed for what it decays |
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
  One seam (`GenotypePriorModel`), two impls; the recipe selects. **The trait is object-safe and
  the selection was written and run at F1** — `Box<dyn GenotypePriorModel>` and
  `Arc<dyn GenotypePriorModel + Send + Sync>` both work, the second across a worker boundary, so
  the loop's own field should carry `+ Send + Sync` rather than the bare `Box<dyn …>`
  [`ng_step_interfaces.md`](ng_step_interfaces.md) sketches. **What does not exist is the recipe**:
  no plan builds one, and the measurement the comparator is kept for also needs the loop and a
  specification for candidate selection. Production ships plug-in default
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
- ~~Impl-time confirmation: the concrete `FittedSpectrum` type from the pre-pass cohort gather.~~
  **Settled 2026-08-22 at step D2** — a borrowed view over the class weights; see §4.
- ~~`OPEN:` which end of the SNP/indel path owns a class split (spec Q1).~~ **Settled 2026-08-23:
  the per-locus expansion — §4.** What remains open under Q1 is a locus carrying one alternative
  of each kind.
  **The scratch slot the row function needs is settled: one `f64` per allele**, which is what the
  ported primitive hoists (`lgamma(α_a)`, one per allele); `CallingScratch` owns it.

## Test & bench shape

Unit tests beside each file, pinning spec §12's eleven properties — notably the 2:1 ratio tripwire
(test 1), the exact-expected-spectrum projection targets built in closed form, never `θ/k` and never
sampled (tests 5–7), bit-equality of seed and leave-one-out at one sample (test 8), the STR seed's
implied-diversity identity and the refusal (tests 10–11), and the pochhammer oracle (test 4). The
end-to-end anchor is the GIAB single-sample 5× regression (genotype accuracy at true variants; count
of true hom-alt called het) run under both impls of the seam — spec §12's definition of done.
