# ng — the joint fit, the estimator: types & interfaces

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

**It runs once, after every sample's walk.** That is a scheduling fact the driver owns, and the only
shape it imposes here is that the entry point takes *every* sample's records at once and returns one
value — there is no per-sample entry point, because there is no per-sample answer (spec §1).

---

## 1. Types

### 1.1 `JointFit` — what the route emits

```rust
/// Every parameter this route produces, for the whole cohort, in one value.
///
/// **Per-sample entries are maps keyed by sample**, not a parallel vector: the fit's
/// output outlives the order the samples were walked in, and a caller looks a sample
/// up by name.
pub struct JointFit {
    /// Chemistry: three numbers per read group (spec §3.1).
    pub noise: BTreeMap<ReadGroupId, Estimate<SiteClassNoise>>,
    /// The panel's allele-frequency spectrum — this route's own parameter, not a
    /// separate step (spec §2.1).
    pub spectrum: Estimate<AlleleFrequencySpectrum>,
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
pub struct SiteClassNoise {
    pub clean: ErrorRate,
    pub noisy: SiteNoise,
    /// `None` where the run has no coverage-by-window summary to condition it on, or
    /// where the class was fitted at zero weight.
    pub duplicated: Option<DuplicatedSiteClass>,
}

/// The class for a locus the sample carries more copies of than the caller assumes.
///
/// **Conditioned on the window's relative copy number, never on the site's own depth**:
/// per-base coverage at 6× cannot tell a two-copy carrier from a sample reading high,
/// so the discriminator is the window (spec §2.2).
pub struct DuplicatedSiteClass {
    pub weight: f64,
    pub alternative_fraction: f64,
}
```

**OPEN:** `DuplicatedSiteClass` is `PROPOSED` in the spec and gated on one measurement (spec §2.2's
closing paragraph). The type is stated so the fit's shape does not have to change when it lands; an
implementation may ship with it always `None`.

### 1.3 `AlleleFrequencySpectrum` — a mixture over allele counts

```rust
/// How common each allele count is in the panel: `2N + 1` weights for `N` diploid
/// samples, summing to one.
///
/// **A fitted parameter of this route**, not a post-hoc computation — it is the outer
/// sum of §2.1's likelihood, and the diversity `Hexp` is read off it directly with no
/// division by an inbreeding coefficient (spec §5.3).
pub struct AlleleFrequencySpectrum {
    /// Index `c` is the weight on "c copies of the alternative allele in the panel".
    weights: Vec<f64>,
}

impl AlleleFrequencySpectrum {
    /// The only constructor: rejects weights that do not sum to one.
    pub fn try_new(weights: Vec<f64>) -> Result<Self, DomainError>;
    /// The panel's expected heterozygosity — spec §5.3.
    pub fn expected_heterozygosity(&self) -> f64;
}
```

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
        /// because that is what says how far to trust it (spec §3.4.3).
        segregating_markers: u64,
    },
    /// The depth or the marker count could not identify it. **Emitted rather than
    /// substituted by a zero** — a caller told "no contamination" would act on it
    /// (spec §3.4.3).
    NotIdentified { reason: NotIdentifiedReason },
}
```

**Contract.** The search is restricted to `α ≤ ½`: the sequence-only likelihood is symmetric, so a
sample swap is invisible by construction and an estimate above a half is not a stronger claim but a
mirror image (spec §3.4.1).

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
```

**Contract.** Deterministic: loci in `KeptLoci` order and samples in name order, so no parameter
varies with thread count and multiple starting points are enumerated rather than sampled (spec §7).
The identity check runs first and completely — a run that would fail on the fiftieth sample fails
before the first likelihood evaluation.

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
    spectrum: &AlleleFrequencySpectrum,
    hom_excess: &BTreeMap<SampleName, HomozygoteExcess>,
) -> f64;
```

**Contract.** `LocusEvidence` borrows every sample's record at one locus and owns nothing; the fit
holds one locus at a time, so the working set is `2N + 1` numbers rather than the cohort (spec §7).

### 2.3 How the maximum is found

Four blocks, alternating (spec §3.3), each of which is a call into machinery that exists:

| block | what moves | what it calls |
|---|---|---|
| 1 | each read group's three noise numbers | a **climb**, not `fit_by_profile_scan` — see the decision below |
| 2 | the spectrum | `climb_mixture_weights`, the spectrum being a mixture over allele counts |
| 3 | each sample's `HomozygoteExcess` | this route's own, one scalar per sample |
| 4 | each sample's contamination | this route's own (spec §3.4), or after the loop |

**Decided — climb rather than scan, from several starting points.** The per-sample route's ladder
scan prices 161 rungs over a few hundred binned cells; here one score is a pass over two million loci
× fifty samples, so a scan is a pass over the data rather than over a table (spec §3.3). **The starting
points must span the separation between the clean and noisy classes**, for the reason the per-sample
route's two-state model records: a start that puts them close together empties one and reports
convergence.

**OPEN:** how many starting points, and spanning what — spec §11 question 3, which needs the profile
curve of question 2. The signature takes a `&[StartingPoint]` so the answer is a constant and not a
shape change.

### 2.4 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum JointFitError {
    /// Samples that did not keep the same loci. **Refuses rather than averaging**
    /// (spec §7).
    #[error("samples disagree on {field}; they did not keep the same loci")]
    IdentityMismatch { field: &'static str },
    /// The panel is too small for the parameter to be identified at all — the spectrum
    /// has `2N + 1` weights, so at two samples it has five (spec §12.5).
    #[error("{parameter} is not identifiable with {samples} samples")]
    NotIdentifiable { parameter: &'static str, samples: usize },
    /// The alternation ran out of passes. **Never reported as convergence**, which is
    /// the failure mode the per-sample route's own termination handling exists for.
    #[error("the joint fit did not converge in {passes} passes")]
    DidNotConverge { passes: u32 },
}
```

---

## 3. Design decisions — decided

- **One entry point taking every sample, no per-sample call.** Fitting against a per-locus frequency
  means a sample's evidence cannot be reduced alone — spec §1.
- **`HomozygoteExcess` is a new type beside `InbreedingF`.** §1.4; spec §5.1.
- **Contamination is an enum, not an `Option<f64>`.** *Not identified* and *zero* are different
  answers and a caller must not be able to read one as the other — spec §3.4.3.
- **The site class is a per-locus latent variable with cohort-level parameters.** `w` and `ε_noisy`
  stay cohort-level; what changes is that each locus's posterior is computed from every sample's
  evidence — spec §2.2.
- **`Hobs` and `π_hom_alt` are derived, not fitted.** They come out of the converged posteriors, which
  is why they share `SampleRates` with the other route but not its provenance — spec §3.2.
- **The spectrum is a parameter of this fit, not a later step.** The `realSFS`-style
  expectation-maximization is the inner loop, not a separate pass — spec §2.1.
- **Nothing here reads the windowed histogram.** The autozygosity coefficient is unreachable from
  scattered loci and stays in the walk — spec §6.
- **OPEN:** whether contamination is a fourth block inside the alternation or a step after it
  (spec §3.4); nothing in these types turns on it.

---

## 4. Reconciliation with existing code

| this doc | existing code | how they meet |
|---|---|---|
| the sum over genotypes at one sample's locus | `NoiseModel` + `SubstitutionNoiseModel` ([`fitting/mod.rs:64`](../../../../src/ng/parameter_estimation/fitting/mod.rs), [`generic/noise_model.rs:257`](../../../../src/ng/parameter_estimation/generic/noise_model.rs)) | used unchanged as the innermost term; `NoiseParams` is an associated type, which is the seam that lets this route carry a third site class without touching the trait |
| the two site classes | `SiteNoise` ([`generic/mod.rs:90`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | `SiteClassNoise` **holds** one rather than restating its two fields, so the pair keeps one home |
| the climb over mixture weights | `climb_mixture_weights` + `MAX_CLIMB_PASSES`, `CLIMB_STILLNESS` ([`fitting/mixture_weights.rs:50,66`](../../../../src/ng/parameter_estimation/fitting/mixture_weights.rs)) | the spectrum is a mixture over allele counts, so the same climb applies with a different component set |
| the alternation's shape and its termination | `CoupledFit` ([`generic/mod.rs:474`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | same loop shape and same rule that running out of passes is not convergence, with four blocks instead of two |
| provenance and evidence count | `Provenance`, `Estimate<T>` ([`parameter_estimation/mod.rs:48,72`](../../../../src/ng/parameter_estimation/mod.rs)) | used as-is on every emitted parameter |
| the per-sample rates | `SampleRates` ([`generic/mod.rs:342`](../../../../src/ng/parameter_estimation/generic/mod.rs)) | reused so the two routes' `Hobs`/`π_hom_alt` are the same type and directly comparable — the comparison is the point |
| the autozygosity coefficient | `InbreedingF` ([`src/ng/types.rs:355`](../../../../src/ng/types.rs)) | **not reused and not renamed**; this route emits `HomozygoteExcess` and never writes to a field of that type |
| `ErrorRate`, `Ploidy`, `GenotypeFrequency`, `ReadGroupId` | [`src/ng/types.rs:311,386,333,199`](../../../../src/ng/types.rs) | used as-is |
| the ladder scan | `fit_by_profile_scan`, `fit_by_fixed_frequency_scan` ([`fitting/ladder_scan.rs:131,248`](../../../../src/ng/parameter_estimation/fitting/ladder_scan.rs)) | **not used** — a rung here is a pass over the data (§2.3). Recorded so an implementer does not reach for them by analogy |

---

## 5. Open items

- **`OPEN:`** the duplicated-locus class (§1.2), gated on the measurement in spec §2.2.
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
parameters, with no reads and no alignments, and require the fit to return what was drawn (spec
§12.1). Three assertions are the ones a plausible-but-wrong implementation would fail: the derived
rates matching the drawn genotypes, which is what catches a biased kept set (spec §12.2); the two
inbreeding coefficients moving in **opposite** directions under a false-heterozygote floor (§12.4);
and an uncontaminated but structured panel returning `α ≈ 0`, which a fit using the pooled spectrum
passes the mixture test and fails this one (§12.6). The panel-size sweep at 2, 5, 10 and 50 samples
reports where each parameter stops being estimable — a number the user needs and nothing currently
states (§12.5).

**No parity oracle, and no bench yet.** Neither route is a port of the other and agreement is the
thing being measured (spec §9). The cost measurements belong to the comparison, at the shape of a
real cohort, and are reported in `doc/devel/ng/reports/` (spec §8).
