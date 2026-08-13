# ng — the joint parameters fit, the estimator: types & interfaces

*Status: architecture draft (2026-08-12), companion to the spec
[`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) (the design and its
rationale) and to the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary)
and [`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types, verbs
for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the **contract** is the
deliverable. Every "why" lives in the spec — this doc does not re-argue one.*

*It reads the records of [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md)
at the loci of [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md), and nothing
else.*

## Module home

`src/ng/parameter_estimation/joint/fit.rs`, with the machinery it shares taken from `fitting/`
(`module_layout.md`: *parameter_estimation/* owns `fitting/`, `generic/`, `ssr/` — this route is a
fourth sub-unit and not a top-level split).

**It runs once, after every sample's genome walk.** That is a scheduling fact the driver owns, and the only
shape it imposes here is that the entry point takes *every* sample's records at once and returns one
value — there is no per-sample entry point, because there is no per-sample answer (spec §1).

---

## 1. Types

### 1.1 `JointFit` — what the route emits

```rust
/// Every parameter this route produces, for the whole cohort, in one value.
///
/// **Per-sample entries are maps keyed by sample**, not a parallel vector: the parameters fit's
/// output outlives the order the samples were walked in, and a caller looks a sample
/// up by name.
pub struct JointFit {
    /// Chemistry: three numbers per read group (spec §3.1).
    pub noise: BTreeMap<ReadGroupId, Estimate<SiteClassNoise>>,
    /// The population's allele-frequency density — four fitted numbers, this route's own
    /// parameter and not a separate step (spec §2.1.2).
    pub density: Estimate<FrequencyDensity>,
    /// The panel allele-count spectrum the density implies. **Derived from `density`,
    /// never fitted beside it** — a consumer asking for "the spectrum" means this, and
    /// emitting both under one name is how the two get confused (spec §2.1.2).
    pub spectrum: AlleleCountSpectrum,
    /// Per sample, the departure from the Hardy–Weinberg proportions the spectrum
    /// predicts. **Not the autozygosity coefficient** (§1.4).
    pub hom_excess: BTreeMap<SampleName, Estimate<HomozygoteExcess>>,
    /// Per sample, the fraction of reads from another individual (spec §3.4).
    pub contamination: BTreeMap<SampleName, ContaminationEstimate>,
    /// Derived from the converged posteriors rather than fitted (spec §3.2).
    pub rates: BTreeMap<SampleName, Estimate<SampleRates>>,
    /// The STR path's four slippage numbers, per read group × stratum (spec §4).
    pub slippage: BTreeMap<(ReadGroupId, Stratum), Estimate<SlippageParams>>,
    /// Per locus, the posterior that it belongs to each site class — the thing no
    /// per-sample marginal can produce (spec §2.2). **Diagnostic**, and the input to
    /// the duplicated-locus report.
    pub site_class_posteriors: SiteClassPosteriors,
}
```

**Contract.** Every entry carries `Estimate<T>`'s provenance and observation count, as the per-sample
route's do — a parameter that fell back rather than being fitted here says so
([`parameter_estimation/mod.rs:48,72`](../../../../src/ng/parameter_estimation/mod.rs)).

### 1.2 `SiteClassNoise` — three classes, not two

```rust
/// A read group's noise, as the three classes a locus can be drawn from (spec §2.2).
///
/// **`clean` and `noisy` are the pair the histogram route already fits**; `duplicated`
/// is this route's addition, and it is the one the error-rate ladder cannot reach —
/// its alternative-read fraction sits near a half rather than on a rung.
/// **Measured 2026-08-12: the noisy class needs about twenty-five samples** (spec §12.5).
/// Below ten its share comes back four to five times the truth, in the same direction at
/// every panel size. So this field's `Estimate`'s provenance is load-bearing rather than
/// informational — a `w` five times too large is a background subtraction, and on a sample
/// at tomato's heterozygosity floor the background is the whole of the answer.
pub struct SiteClassNoise {
    pub clean: ErrorRate,
    pub noisy: SiteNoise,
    /// `None` where the run has no coverage-by-window summary to condition it on, or
    /// where the class was fitted at zero weight.
    pub duplicated: Option<DuplicatedSiteClass>,
}

/// The class for a locus **a given sample** carries more copies of than the caller
/// assumes.
///
/// **Its grain is the (locus, sample) pair, where `clean` and `noisy` are properties of
/// the locus alone** — a collapsed paralog mismaps in every sample, and a duplication is
/// carried by an individual (spec §2.2). So only these two numbers are cohort-level;
/// **membership is decided per sample**, from that sample's window coverage.
///
/// **Conditioned on the window's relative copy number, never on the site's own depth**:
/// per-base coverage at 6× cannot tell a two-copy carrier from a sample reading high,
/// so the discriminator is the window (spec §2.2). The window summary it reads is
/// `records.rs`'s third object, not something this parameters fit derives.
pub struct DuplicatedSiteClass {
    pub weight: f64,
    pub alternative_fraction: f64,
}
```

**OPEN:** `DuplicatedSiteClass` is `PROPOSED` in the spec and gated on one measurement (spec §2.2's
closing paragraph). The type is stated so the parameters fit's shape does not have to change when it lands; an
implementation may ship with it always `None`.

### 1.3 `FrequencyDensity` — four numbers, and the spectrum it implies

**Changed 2026-08-12, and the change is the estimator's, not the type's.** An earlier version of this
section held one weight per allele count, `2N + 1` of them, on the strength of a `realSFS` reference
the spec has since withdrawn: that estimator conditions on the panel's allele *count*, and the
cancellation which makes it work fails as soon as a per-sample inbreeding coefficient enters — which
this route requires (spec §2.1.1). What is fitted is a **population frequency density**, and a free
weight per grid point would be an unregularised deconvolution (spec §2.1.2).

```rust
/// How the population's allele frequency is distributed at an ordinary position:
/// a mass on invariant, a mass on fixed-non-reference, and a Beta over what segregates.
///
/// **Four fitted numbers, and the quadrature over `f` is accuracy rather than freedom** —
/// doubling the nodes costs time and adds no parameter (spec §2.1.2).
pub struct FrequencyDensity {
    p_invariant: f64,
    p_fixed_alt: f64,
    /// The Beta's shape over the segregating sites. `a < 1` is the rare-allele pile-up.
    a: f64,
    b: f64,
}

impl FrequencyDensity {
    /// The only constructor: refuses masses that do not leave room for the Beta, and
    /// shapes outside `(0, ∞)`.
    pub fn try_new(p_invariant: f64, p_fixed_alt: f64, a: f64, b: f64) -> Result<Self, DomainError>;
    /// `∫ π(f) · 2 f (1 − f) df` — the population's expected heterozygosity, with no
    /// finite-sample correction because there is no panel in it (spec §5.3).
    pub fn expected_heterozygosity(&self) -> f64;
    /// The panel allele-count spectrum this density implies, given each sample's
    /// homozygote excess. **The inbreeding belongs in the signature**: it changes how a
    /// frequency turns into counts, and a spectrum computed without it is a different
    /// panel's.
    pub fn implied_spectrum(
        &self,
        ploidy: Ploidy,
        hom_excess: &[HomozygoteExcess],
    ) -> AlleleCountSpectrum;
}

/// How common each allele count is in the panel: `2N + 1` weights for `N` diploid
/// samples, summing to one. **Derived, never fitted** (spec §2.1.2).
pub struct AlleleCountSpectrum {
    weights: Vec<f64>,
}
```

**Measured 2026-08-12: one Beta is enough, and a second is not worth its three numbers.** Against a
drawn cohort whose true density has two bumps, four numbers return `Hexp` 4.9% high and seven return
5.8% high — no better on the truth the second bump exists for (spec §11 question 6). **So this struct
stays at four fields**, and an implementation that finds itself wanting a fifth should re-read that
measurement before adding one.

### 1.4 `HomozygoteExcess` — a second inbreeding coefficient, and a second type

```rust
/// How much less heterozygous an individual is than random mating in the panel would
/// predict: `1 − Hobs/Hexp`, with `Hexp` from the fitted spectrum.
///
/// **A different quantity from [`InbreedingF`]**, which is the fraction of a genome
/// lying in runs of homozygosity — the one the caller's genotype prior multiplies.
/// Two types rather than two values of one type, because a consumer handed the wrong
/// one gets a plausible answer (spec §5).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct HomozygoteExcess(f64);

impl HomozygoteExcess {
    /// **The only constructor, and it refuses anything outside `[0, 1]`.** `F_IS` is
    /// negative under a heterozygote excess, and an unconstrained fit will go there —
    /// which is exactly the escape route spec §2.2 argues the duplicated loci do not
    /// have. A negative value books one sample's mismapping as biology, with a
    /// plausible number and no error (spec §5.1).
    pub fn try_new(value: f64) -> Result<Self, DomainError>;
}
```

**Decided — a new newtype, and `InbreedingF` is left alone.** `InbreedingF`
([`src/ng/types.rs:355`](../../../../src/ng/types.rs)) is already the autozygosity coefficient
everywhere it is used, and the per-sample route is built on it; renaming it would touch shipped code
for a documentation gain. **What the compiler must refuse is the mix**, and two types do that (spec
§5.1). *Its doc comment carries the mapping to the literature: `F_ROH` for `InbreedingF`, `F_IS` for
this.*

### 1.5 `ContaminationEstimate` — a number, or a reason there isn't one

```rust
/// The fraction of a sample's reads that came from another individual (spec §3.4).
pub enum ContaminationEstimate {
    Estimated {
        alpha: Contamination,
        /// How many segregating markers stood behind it. **Reported beside the value**,
        /// because that is what says how far to trust it (spec §3.4.4).
        segregating_markers: u64,
    },
    /// The depth or the marker count could not identify it. **Emitted rather than
    /// substituted by a zero** — a caller told "no contamination" would act on it
    /// (spec §3.4.4).
    NotIdentified { reason: NotIdentifiedReason },
}
```

**Contract.** The search is restricted to `α ≤ ½`: the sequence-only likelihood is symmetric, so a
sample swap is invisible by construction and an estimate above a half is not a stronger claim but a
mirror image (spec §3.4.1).

### 1.6 `SequencingBatches` — who was sequenced beside whom, and it comes from the user

**The contaminant is a neighbouring library, not a random member of the species** (spec §3.4.3), so
the population the contaminant's genotype is drawn against is *the samples that ran together*. That
grouping is **stated, never inferred**: it is absent from both cohorts' alignments — the tomato
archive's `@RG` lines carry no `PU` and SRA rewrote the read names to
`SRR7279481.37559618:TTAGGC:37559618`, keeping the barcode and losing the flowcell — and a pipeline
that guessed it from what survives would be wrong silently.

```rust
/// Which read groups were sequenced together, as the run was told.
pub struct SequencingBatches {
    /// Every read group of one batch. **Read groups and not samples**: one sample's
    /// libraries may have run on different flowcells, and the read group is the grain
    /// the header gives.
    batches: Vec<BTreeSet<ReadGroupId>>,
}

impl SequencingBatches {
    /// Every read group in one batch — **the default**, and the honest statement of what
    /// a run knows when nobody has said otherwise.
    pub fn all_together(groups: &ReadGroups) -> Self;

    /// The batching the CLI was given.
    ///
    /// # Errors
    ///
    /// [`JointFitError::ReadGroupNotBatched`] naming every read group the batching left
    /// out. **A run that names any batch must name them all**: a user who lists three
    /// plates and forgets four samples would otherwise get a wrong contaminant prior for
    /// those four with nothing said.
    pub fn from_groups(
        groups: &ReadGroups,
        batches: Vec<BTreeSet<ReadGroupId>>,
    ) -> Result<Self, JointFitError>;

    /// The batch this read group belongs to.
    pub fn batch_of(&self, group: ReadGroupId) -> &BTreeSet<ReadGroupId>;

    /// Whether this is the default. **Travels with `α`** (spec §3.4.3): two runs under
    /// different batchings produce different numbers and neither is comparable to the
    /// other, and a number fitted under one batch for the whole cohort is the weaker
    /// kind.
    pub fn is_default(&self) -> bool;
}
```

**Contract.** A partition: every read group in exactly one batch, checked at construction. The default
is one batch holding everything, so the type is never optional and no consumer branches on its
absence.

**Where `PU` fits, and it is a default rather than an answer.** `ReadGroup`
([`read/input/read_groups.rs`](../../../../src/ng/read/input/read_groups.rs)) reads `ID`, `SM`, `LB`
and `PL` and has no platform-unit field. Adding one, and seeding the batching from it where a file
declares it, is worth doing — GIAB's read names carry `HISEQ1:23:H9UD5ADXX:2:…` even though its `PU`
says `unknown` — but **the user's `--sequenced-together` always wins**, because a declared `PU` is as
untrustworthy as the `PL` this module already refuses to group by.

---

## 2. Interfaces

### 2.1 The entry point

```rust
/// Fit every parameter from every sample's records, once.
///
/// # Errors
///
/// [`JointFitError::IdentityMismatch`] before any arithmetic, when two samples did not
/// keep the same loci — the refusal `loci.rs` defines and this call enforces.
pub fn fit_jointly(
    samples: &[SampleRecords],
    config: &JointFitConfig,
) -> Result<JointFit, JointFitError>;

/// What the run was asked for, beside the records.
pub struct JointFitConfig {
    /// Who was sequenced beside whom (§1.6). `SequencingBatches::all_together` by default.
    pub batches: SequencingBatches,
    /// How many principal components stand behind the individual-specific allele
    /// frequencies — spec §11 question 4, still open.
    pub components: usize,
    // … the starting points, the quadrature, the ladder
}
```

**Contract.** Deterministic: loci in `KeptLoci` order and samples in name order, so no parameter
varies with thread count and multiple starting points are enumerated rather than sampled (spec §7).
The identity check runs first and completely — a run that would fail on the fiftieth sample fails
before the first likelihood evaluation.

**`&[SampleRecords]` was the fifty-sample signature and it does not survive a thousand — CHANGED
2026-08-13.** A slice of whole record sets requires every sample's whole evidence resident, which is
6 GB at a thousand samples and 30 GB at five thousand (spec §7, §11 question 10). What the estimator
actually consumes is narrower: every sample's generic records, and then — after those are dropped —
every sample's tracts for one band of strata, a band at a time (spec §4, §4.1). **So the parameter is
the cohort, which lends sections for the length of a call and cannot be made to hand one over**
(records arch §2.2):

```rust
pub fn fit_jointly(
    records: &mut CohortRecords,
    config: &JointFitConfig,
) -> Result<JointFit, JointFitError>;
```

`CohortRecords` is built from every sample's `SampleRecords` — resident ones in the run that never
writes a file, file-backed ones otherwise — and this signature does not distinguish the two, which is
the point. The identity check across samples happens when it is built, before any section is decoded.

**What is still open** is whether the fit also reads locus-major within the generic half, which
question 10 measures. That would add a call to `CohortRecords`, not change the ones here: the unit of
lending would become a range of loci across every sample instead of a whole section, and the scoped
shape is what makes that additive rather than a redesign.

### 2.2 The likelihood, as the shared seam already shapes it

**The innermost sum is `NoiseModel::append_genotype_likelihoods`, unchanged**
([`fitting/mod.rs:64`](../../../../src/ng/parameter_estimation/fitting/mod.rs)). What this route adds
is the two sums outside it — over the locus's allele frequency, and over the site class — and those
are its own code rather than the trait's:

```rust
/// One locus's contribution to the log-likelihood, at these parameters.
///
/// Read outward, this is spec §3.1's expression: a sample's reads scored under each
/// genotype and summed against what the allele frequency implies; the samples
/// multiplied; the site class summed over; the frequency summed over the spectrum.
fn locus_log_likelihood(
    locus: &LocusEvidence<'_>,
    noise: &BTreeMap<ReadGroupId, SiteClassNoise>,
    density: &FrequencyDensity,
    quadrature: &FrequencyQuadrature,
    hom_excess: &BTreeMap<SampleName, HomozygoteExcess>,
) -> f64;
```

**Contract — which base is the alternative one is summed over inside this call, never chosen by the
caller.** Three terms with an equal prior, and only on the segregating branch: the invariant term
carries no allele at all (spec §3.1.1). **There is no `alt: ObservedAllele` parameter and there must
not be one** — a signature that took the allele would push the choice up to whoever builds
`LocusEvidence`, which is where the maximum-of-three bias enters.

**Contract.** `LocusEvidence` borrows every sample's record at one locus and owns nothing; the parameters fit
holds one locus at a time, so the working set is one number per frequency quadrature node rather than the cohort (spec §7).

### 2.3 How the maximum is found

Four blocks, alternating (spec §3.3), each of which is a call into machinery that exists:

| block | what moves | what it calls |
|---|---|---|
| 1 | each read group's three noise numbers | a **climb**, not `fit_by_profile_scan` — see the decision below |
| 2 | the density's four numbers | `climb_mixture_weights` for the two masses; the Beta's `a` and `b` are **not** weights and need their own climb (spec §2.1.2) |
| 3 | each sample's `HomozygoteExcess` | this route's own, one scalar per sample, constrained to `[0, 1]` (§1.4) |
| 4 | each sample's contamination | this route's own (spec §3.4), or after the loop |

**Blocks 2 and 3 can trade against one another** — both control how often a genotype comes out
heterozygous, and they are separated only by *where* the variation sits, across loci for one and
within a locus for the other. Not a shape change, but the alternation's termination has to be judged
on both together rather than block by block, and spec §12.9 measures how strong the trade is at 5, 10
and 50 samples.

**Decided — climb rather than scan, from several starting points.** The per-sample route's ladder
scan prices 161 rungs over a few hundred binned cells; here one score is a pass over two million loci
× fifty samples, so a scan is a pass over the data rather than over a table (spec §3.3). **The starting
points must span the separation between the clean and noisy classes**, for the reason the per-sample
route's two-state model records: a start that puts them close together empties one and reports
convergence.

**Measured 2026-08-12, and it settles what they must span.** A start that puts the clean and the noisy
class *close together* collapses them into one and reports convergence, costing **46% of the clean
error rate** and putting `Hexp` 10.6% high; a start that puts them far apart returns the clean rate to
−1.3% on the identical data (spec §11 question 3). **So the starting points are part of the
estimator, not a tuning detail**, and `fit_jointly` takes several and returns the best-scoring fit
rather than the last one.

**OPEN:** how many. Three was enough for one of them to land well; whether nine is needed, as the
per-sample route settled for its own two-state model, wants the profile curve of spec §11 question 2.
The signature takes a `&[StartingPoint]` so the answer is a constant and not a shape change.

### 2.4 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum JointFitError {
    /// Samples that did not keep the same loci. **Refuses rather than averaging**
    /// (spec §7).
    #[error("samples disagree on {field}; they did not keep the same loci")]
    IdentityMismatch { field: &'static str },
    /// The panel is too small for the parameter to be identified at all. At one sample
    /// this is `HomozygoteExcess` and contamination — the route still fits everything
    /// else and degenerates to the per-sample estimator (spec §6.1, §12.5).
    #[error("{parameter} is not identifiable with {samples} samples")]
    NotIdentifiable { parameter: &'static str, samples: usize },
    /// The alternation ran out of passes. **Never reported as convergence**, which is
    /// the failure mode the per-sample route's own termination handling exists for.
    #[error("the joint parameters fit did not converge in {passes} passes")]
    DidNotConverge { passes: u32 },
    /// A batching was given and left read groups out (§1.6). **Naming them is the
    /// point**: the user forgot a plate, and the alternative is a wrong contaminant
    /// prior for exactly those samples with nothing said.
    #[error("{} read groups were not assigned to a sequencing batch: {}", .groups.len(), .groups.join(", "))]
    ReadGroupNotBatched { groups: Vec<String> },
}
```

---

## 3. Design decisions — decided

- **One entry point taking every sample, no per-sample call.** Fitting against a per-locus frequency
  means a sample's evidence cannot be reduced alone — spec §1.
- **`HomozygoteExcess` is a new type beside `InbreedingF`.** §1.4; spec §5.1.
- **Who was sequenced beside whom is an input with a default, never an inference** — §1.6, spec
  §3.4.3. One batch holding every read group unless the CLI says otherwise; a partial batching is
  refused rather than completed.
- **Contamination is an enum, not an `Option<f64>`.** *Not identified* and *zero* are different
  answers and a caller must not be able to read one as the other — spec §3.4.4.
- **The site class is a per-locus latent variable with cohort-level parameters.** `w` and `ε_noisy`
  stay cohort-level; what changes is that each locus's posterior is computed from every sample's
  evidence — spec §2.2.
- **`Hobs` and `π_hom_alt` are derived, not fitted.** They come out of the converged posteriors, which
  is why they share `SampleRates` with the other route but not its provenance — spec §3.2.
- **A population frequency density with four numbers, not one weight per allele count.** The count
  form cannot carry per-sample inbreeding or contamination, and a free weight per grid point is an
  unregularised deconvolution — spec §2.1.1, §2.1.2. The panel allele-count spectrum is derived from
  it and emitted beside it under a different name (§1.3).
- **Which non-reference base segregates is summed over, never chosen from the data** — spec §3.1.1,
  and §2.2 above states it as a signature constraint.
- **`HomozygoteExcess` is constrained to `[0, 1]` by its constructor** — §1.4; spec §5.1.
- **The duplicated class is per (locus, sample); the other two are per locus** — §1.2; spec §2.2.
- **The route runs at one sample and says what it cannot fit there.** `HomozygoteExcess` and
  contamination come back `NotIdentifiable`; everything else is fitted, and the likelihood is the
  per-sample estimator — spec §6.1.
- **The spectrum is a parameter of this parameters fit, not a later step** — spec §2.1.
- **Nothing here reads the windowed histogram.** The autozygosity coefficient is unreachable from
  scattered loci and stays in the genome walk — spec §6.
- **OPEN:** whether contamination is a fourth block inside the alternation or a step after it
  (spec §3.4); nothing in these types turns on it. **2026-08-13 tilts it towards "after"**: spec
  §3.4.4 now re-fits `α` alone over a census far larger than the other blocks see, streamed a locus at
  a time across every sample's file, with each sample's principal-component coordinates held at what
  the small census gave. A block that reads a different set of loci from its neighbours is easier to
  reason about outside the alternation than inside it, and the cost is one extra sweep.

---

## 4. Reconciliation with existing code

| this doc | existing code | how they meet |
|---|---|---|
| the sum over genotypes at one sample's locus | `NoiseModel` + `SubstitutionNoiseModel` ([`fitting/mod.rs:64`](../../../../src/ng/parameter_estimation/fitting/mod.rs), [`generic/noise_model.rs:257`](../../../../src/ng/parameter_estimation/generic/noise_model.rs)) | used unchanged as the innermost term; `NoiseParams` is an associated type, which is the seam that lets this route carry a third site class without touching the trait |
| the two site classes | `SiteNoise` ([`generic/mod.rs:90`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | `SiteClassNoise` **holds** one rather than restating its two fields, so the pair keeps one home |
| the climb over mixture weights | `climb_mixture_weights` + `MAX_CLIMB_PASSES`, `CLIMB_STILLNESS` ([`fitting/mixture_weights.rs:50,66`](../../../../src/ng/parameter_estimation/fitting/mixture_weights.rs)) | the density's two masses are mixture weights and reuse it. **The Beta's `a` and `b` are not**, and reaching for this function for them is the mistake this row exists to head off |
| the alternation's shape and its termination | `CoupledFit` ([`generic/mod.rs:474`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | same loop shape and same rule that running out of passes is not convergence, with four blocks instead of two |
| provenance and evidence count | `Provenance`, `Estimate<T>` ([`parameter_estimation/mod.rs:48,72`](../../../../src/ng/parameter_estimation/mod.rs)) | used as-is on every emitted parameter |
| the per-sample rates | `SampleRates` ([`generic/mod.rs:342`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | reused so the two routes' `Hobs`/`π_hom_alt` are the same type and directly comparable — the comparison is the point |
| the autozygosity coefficient | `InbreedingF` ([`src/ng/types.rs:355`](../../../../src/ng/types.rs)) | **not reused and not renamed**; this route emits `HomozygoteExcess` and never writes to a field of that type |
| `ErrorRate`, `Ploidy`, `GenotypeFrequency`, `ReadGroupId` | [`src/ng/types.rs:311,386,333,199`](../../../../src/ng/types.rs) | used as-is |
| the ladder scan | `fit_by_profile_scan`, `fit_by_fixed_frequency_scan` ([`fitting/ladder_scan.rs:131,248`](../../../../src/ng/parameter_estimation/fitting/ladder_scan.rs)) | **not used** — a rung here is a pass over the data (§2.3). Recorded so an implementer does not reach for them by analogy |

---

## 5. Open items

- **`OPEN:`** the duplicated-locus class (§1.2), gated on the measurement in spec §2.2 — and on the
  window summary it reads, which is `records.rs`'s third object rather than this unit's.
- **Decided by measurement:** one Beta, not two (§1.3) — spec §11 question 6.
- **Impl-time:** the frequency quadrature's node count. A `FrequencyQuadrature` value threaded
  through `locus_log_likelihood`, so the accuracy knob is visible in the signature and cannot be
  mistaken for a parameter count.
- **`OPEN:`** the number of principal components behind contamination's individual-specific allele
  frequencies — spec §11 question 4. The type is a `Vec<f64>` of coordinates either way.
- **`OPEN:`** starting points for block 1 — spec §11 question 3.
- **`OPEN:`** whether the short-tract STR strata contribute to contamination — spec §11 question 5.
  Until it is settled, `fit_jointly` reads contamination from the generic records only.
- **Impl-time:** `Stratum` as a key type. The catalog keys strata by `(u8, u64)`
  ([`repeat_catalog/strata.rs:21`](../../../../src/ng/repeat_catalog/strata.rs)); pin whether this
  route takes that tuple or a newtype over it, and do not mint a second stratum notion.
- **Impl-time:** whether `SiteClassPosteriors` is materialised for every locus or streamed. Two
  million floats per class is small; the contract is only that a consumer can ask about a named locus.

---

## 6. Test and bench shape

Tests live in `joint/fit.rs`'s `#[cfg(test)] mod tests`, and — as with the per-sample route — **the
oracle is a world whose truth is known rather than a fixture**: fill records directly from drawn
parameters, with no reads and no alignments, and require the parameters fit to return what was drawn (spec
§12.1). Three assertions are the ones a plausible-but-wrong implementation would fail: the derived
rates matching the drawn genotypes, which is what catches a biased kept set (spec §12.2); the two
inbreeding coefficients moving in **opposite** directions under a false-heterozygote floor (§12.4);
and an uncontaminated but structured panel returning `α ≈ 0`, which a fit using the pooled spectrum
passes the mixture test and fails this one (§12.6). **Measured 2026-08-12 and the assertion needs
turning round**: the pooled spectrum does not inflate `α` on a structured panel, it *deflates* it — a
sample truly at 3% comes back at 0.5% at an `F_st` of 0.20 — so the test that catches it is a **spiked**
structured panel whose contaminated sample must be found, not a clean one whose samples must all read
zero. A clean structured panel is passed by the broken fit
(`../reports/joint_contamination_2026-08-12.md` §3). **A fourth assertion belongs beside them**: a
batching that leaves a read group out is refused and names it (§1.6). The panel-size sweep at 2, 5, 10 and 50 samples
reports where each parameter stops being estimable — a number the user needs and nothing currently
states (§12.5).

**No parity oracle, and no bench yet.** Neither route is a port of the other and agreement is the
thing being measured (spec §9). The cost measurements belong to the comparison, at the shape of a
real cohort, and are reported in `doc/devel/ng/reports/` (spec §8).
