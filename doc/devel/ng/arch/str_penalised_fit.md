# The repeat-tract fit by penalised EM: types & interfaces

**Architecture draft, 2026-10-14**, companion to two approved specs:
[`../spec/str_penalised_fit.md`](../spec/str_penalised_fit.md) (the fit) and
[`../spec/str_fit_to_calling.md`](../spec/str_fit_to_calling.md) (the hand-over to calling). Every
*why* is there; this doc gives the code shape. Shared conventions:
[`module_layout.md`](module_layout.md). Line numbers are on branch `str-penalised-fit` after merging
main (`e0b41ad4`).

---

## Module home

- **`parameter_estimation/joint/ssr_fit/penalised_fit.rs` becomes the repeat-tract fit.** It already
  holds the loop, the curves, the stiffness search and the errors; it gains the production entry
  point (§2.1) and the conversion into `StratumOutcome`.
- **`ssr_fit/composite_fit.rs` is renamed `ssr_fit/expected_counts.rs`**: once the climb goes, what it
  owns for production is step 1 of the loop (`expected_counts`, `ClassCounts`) and the pair likelihood
  (`SamplePairs`). The counting fit's own drivers (`fit_class_by_counting*`) stay in it as the
  zero-penalty oracle the tests compare against (spec §4.1).
- **`ssr_fit/counting_fit_errors.rs`** keeps the sandwich pieces the penalised errors share.
- **`ssr_fit.rs` keeps the shared vocabulary** (`Slippage`, `StratumEvidence`, `gather_strata`,
  `StratumOutcome`, the provenance types). The climb's code in it is deleted at plan milestone E.
- **`calling/parameters_file/`** owns format 3 and the migration; **`stratum_fits.rs`** owns the new
  ladder step.

One implementation, so **no trait**: the research bench's `PeriodFit` trait stays in
`bench_fixtures` for comparing methods and is not production surface.

---

## 1. Types

### 1.1 The production settings

The research defaults of `PenalisedFitConfig` (`penalised_fit.rs:130-155`) are not what production
runs (spec §2.7). Production gets one named constructor; the research fields that only measured
alternatives are removed (`allele_prior_fixed_draws`, `allele_prior_skips_within`,
`weigh_allele_draws_by_locus`, `Smoothing::Fixed` stays for the zero-penalty test only).

```rust
impl PenalisedFitConfig {
    /// **What `estimate-parameters` runs** (spec §2.7): start 0.05 / 0.5 / 0.3, at most 2,000
    /// rounds, stop below 1e-7, one locus's worth of allele draws plus one draw a class as the
    /// floor, F from pairs at most 10 times, standard errors on.
    pub fn production() -> Self;
}

/// One tract's worth of allele draws a class, spread over the lengths (spec §2.4). Measured.
pub const ALLELE_PRIOR_LOCI: f64 = 1.0;
/// Readings of F off pairs of samples, at most (spec §2.5). The bench's value; soft.
pub const PAIR_ALTERNATIONS: usize = 10;
/// A share's error is withheld when the class's loci with reads times the share is below this
/// (spec §7, owner 2026-10-14). Rests on two settings; soft.
pub const LOCI_FOR_A_SHARE_ERROR: f64 = 3.0;
```

### 1.2 What the fit returns per class: the existing `StratumOutcome`

**No new outcome type.** The fit produces `StratumOutcome::Fitted` for every class with reads and
`Derived` for a class without (values read off the curves), as the climb does
(`ssr_fit.rs:854-865`). Two changes to `StratumFit` (`ssr_fit.rs:537`):

```rust
pub struct StratumFit {
    pub stratum: Stratum,
    /// One entry, slippage group 0 (spec §1.4: one stutter for all read groups).
    pub slippage: Vec<Option<Slippage>>,
    pub allele_length_variation: Vec<f64>,
    /// **`None` where no locus of the motif length has two samples with reads** (spec §2.5).
    /// Was `f64`.
    pub concentration: Option<f64>,
    pub standard_errors: Option<StratumErrors>,
    pub level_provenance: Vec<Option<LevelProvenance>>,
    pub shares_provenance: Vec<Option<SharesProvenance>>,
    pub samples_fitted_on: Option<usize>,
    pub tracts_of_its_own: usize,
    pub reads_crossing: u64,
    /// **The climb's own record** (`ending`, `walks`, `log_likelihood_a_tract`, `borrowed`,
    /// `tracts_fitted`), read only by the climb's log summaries. `None` from the penalised fit;
    /// the field is deleted with the climb (plan milestone E).
    pub climb: Option<ClimbRecord>,
}
```

*The `climb` wrapper is the seam that lets both fits build one struct until the climb is deleted*;
nothing outside `ssr_fit.rs` reads those fields (spec §2.1; `StratumFits::over` reads only the fields
kept above, `stratum_fits.rs:504-607`).

### 1.3 Where a number came from

A value on a penalised curve is neither a cell's own fit nor main's parametric curve (spec
`str_fit_to_calling.md` §2.3). One new variant in each source enum:

```rust
pub enum LevelSource {          // slippage_curve.rs:519
    Cell,
    Curve,
    Blend { curve_weight: f64 },
    /// A point on the motif length's penalised curve, fitted with this class's reads.
    PenalisedCurve(PenalisedPoint),
}

pub enum ShareSource {          // share_curve.rs:580
    Stratum,
    Curve,
    Blend { curve_weight: f64 },
    PenalisedCurve(PenalisedPoint),
}

/// What the file says about a point on a penalised curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PenalisedPoint {
    /// The curve's stiffness strength (larger is stiffer).
    pub strength: f64,
    /// Whether a bend was eased at this repeat count: the reads insisted on a kink.
    pub eased_here: bool,
}
```

`LevelProvenance::reach` (`ssr_fit.rs:733-743`) carries whether the class lies inside the range of
classes with reads; `slipped_reads` stays the class's expected slipped reads. `curve` is `None` for a
penalised point.

### 1.4 Errors for a withheld share

```rust
pub enum StratumError {         // ssr_fit.rs:632
    Estimated(f64),
    NotIdentified,
    NotPlaced,
    NoShare,
    /// **Withheld: too few of the class's loci carry this length** for its error to be honest
    /// (`LOCI_FOR_A_SHARE_ERROR`; spec §7).
    TooFewLoci,
}
```

### 1.5 Calling's tables (companion spec §3)

```rust
pub struct FittedAlleleLengthVariation {      // stratum_fits.rs:379
    pub weights: Vec<f64>,
    /// **`None` where the fit reported no F** (companion spec §3.2). Was `f64`.
    pub concentration: Option<f64>,
}

pub enum AlleleLengthVariationRung {          // stratum_fits.rs:359
    StratumsOwnFit,
    /// The class's own allele length variation, held at the stated concentration because its F
    /// is absent (companion spec §3.2, owner 2026-10-14).
    StratumsOwnSharesStatedSpread,
    PeriodsPooledTracts,
    StatedFlat,
}
```

`checked_allele_length_variation` (`stratum_fits.rs:1022`) checks the concentration only when it is
`Some`; `median_concentration` (`:1100`) takes the median over the classes that have one.

### 1.6 The parameters file, format 3

```rust
pub const FORMAT_VERSION: u32 = 3;                     // parameters_file/mod.rs:192

pub struct StratumAlleleLengthVariationRow {           // mod.rs:1354
    pub period: u8,
    pub reference_repeats: u64,
    /// F = 1/(1+c). Absent where the fit's F is.
    pub between_locus_share: Option<f64>,
    pub between_locus_share_standard_error: Option<f64>,
    pub share_of_chromosomes_by_repeat_offset: Vec<f64>,
}
// PeriodAlleleLengthVariationRow (mod.rs:1374): the same three renames.

pub struct RepeatTracts {
    /// Was `fallback_length_spectrum_concentration` (mod.rs:1065); now F.
    pub fallback_between_locus_share: WarrantedValue,
    // … unchanged …
}

pub enum LevelSmoothing {                               // mod.rs:1199
    // … existing variants, still readable …
    /// A point on a penalised curve (§1.3).
    PenalisedCurve { strength: f64, eased_here: bool, reach: CurveReach },
}
// ShareSmoothing (mod.rs:1226): the same variant.
```

---

## 2. Interfaces

### 2.1 The fit

```rust
/// **Fit every repeat-tract class of a cohort**, one motif length at a time, by penalised EM
/// (spec §2). Returns one `StratumOutcome` per class handed in, in the same order.
pub fn fit_repeat_tracts(
    strata: &[StratumEvidence],
    homozygote_excess: &[f64],
    config: &PenalisedFitConfig,
) -> Vec<StratumOutcome>;
```

**Contract.**
- Groups `strata` by period and calls `fit_period_by_penalised_counting` (`penalised_fit.rs:328`)
  once per period; periods run in parallel, and within one the E-step runs over classes and loci in
  parallel.
- **Deterministic:** the same bits on every run and at every thread count. Every sum over loci,
  samples or classes is taken in index order (spec §3, first trap).
- **Memory:** at most one class's table of distinct rows alive per worker at once (spec §3).
- A class with reads → `Fitted`; with none → `Derived` from the curves inside their range, held
  flat outside (`CurveReach`), or `Refused { NoSpanningReads }` when its period has no curve at all.
  No refusal floor (spec §2.2).
- `concentration` is `None` when no locus of the period has two samples with reads;
  `working_concentration` is never written into a `StratumFit`.
- Infallible: no input it accepts makes it fail; malformed evidence is a caller bug and panics as
  `gather_strata`'s outputs already do.

### 2.2 Where it is called

`fit_a_cohort` (`run/census_fit.rs:220`) replaces
`ssr_fit::fit_strata_on_sample_subsets(&evidence, &homozygote_excess, &order, tracts)` (`:262`) with
`ssr_fit::penalised_fit::fit_repeat_tracts(&evidence, &homozygote_excess, &PenalisedFitConfig::production())`.
`sample_order` (`:259-260`) is no longer needed. The `CohortFit` it returns is unchanged.

### 2.3 The file's edge

```rust
/// F = 1 / (1 + concentration). The only place the formula appears (companion spec §2.1).
pub fn between_locus_share_of(concentration: f64) -> f64;
/// concentration = (1 − F) / F.
pub fn concentration_of(between_locus_share: f64) -> f64;

/// **Rewrite a format-1 or format-2 file's parsed TOML into format 3** before it is
/// deserialised (companion spec §2.2): the renamed keys, each `concentration` converted to
/// `between_locus_share`, the fallback's value converted with its warrant kept. A
/// non-positive or non-finite concentration is left as is, so that validation refuses it.
fn migrate_to_format_3(file: &mut toml::Value);   // from_toml.rs, called from from_toml (:90)
```

Both conversions live in `calling/parameters_file/`; nothing else in the crate computes F.

### 2.4 The CLI

`--str-param-estimates-at-once` (`cli/estimate_parameters.rs:145`) becomes a hidden argument that,
if given, exits with an error saying the repeat-tract fit no longer holds per-stratum tables (owner,
spec §7). Its doc and the `strata_at_once` plumbing (`:361-364`) go.

---

## 3. Design decisions — decided

- **One implementation, no trait.** Nothing competes with it in production (spec §1).
- **`StratumOutcome` is reused, not replaced**, because `StratumFits::over` already reads it and
  is unchanged except for the optional concentration (spec §2.1).
- **The climb's record is wrapped, then deleted**, so both fits share `StratumFit` until milestone E.
- **New provenance variants, not reuse of `ThisStratum`**: a pulled value is not the stratum's own
  (companion spec §2.3).
- **Migration on the parsed TOML, not serde aliases**, because the concentration's value changes,
  not just its key (companion spec §2.2).
- **The bundle seed uses the edge rule run by run** (`run_lengths_at_edge`, `seed_ssr.rs`; companion
  spec §3.3).

## 4. Reconciliation with existing code

| new / changed | existing code | how |
|---|---|---|
| `fit_repeat_tracts` | `fit_period_by_penalised_counting`, `penalised_fit.rs:328` | wraps it per period; adds the `StratumOutcome` conversion |
| production settings | `PenalisedFitConfig::default`, `penalised_fit.rs:130-155`; bench overrides `fit_bench.rs:497-515` | `production()` holds the bench's values plus `ALLELE_PRIOR_LOCI` |
| step 1 and pairs | `expected_counts`, `ClassCounts`, `SamplePairs` in `composite_fit.rs` | moved to `expected_counts.rs` unchanged |
| outcome | `StratumOutcome`, `ssr_fit.rs:854`; `StratumFit`, `:537` | reused; `concentration: Option`, `climb: Option<ClimbRecord>` |
| derived class | `DerivedStratum`, `ssr_fit.rs:826` | reused for classes with no reads |
| provenance | `LevelProvenance`, `ssr_fit.rs:733`; `LevelSource`, `slippage_curve.rs:519`; `ShareSource`, `share_curve.rs:580`; `CurveReach`, `slippage_curve.rs:654` | one variant each; `reach` reused |
| errors | `StratumErrors`, `ssr_fit.rs:611`; `StratumError`, `:632` | one variant |
| calling's tables | `FittedAlleleLengthVariation`, `stratum_fits.rs:379`; rung `:359`; `allele_length_variation_at` `:781`; `checked_allele_length_variation` `:1022`; `median_concentration` `:1100` | optional concentration; new rung |
| file | `FORMAT_VERSION`, `mod.rs:192`; rows `:1354`, `:1374`; fallback `:1065`; `LevelSmoothing` `:1199`; `ShareSmoothing` `:1226`; `from_toml`, `from_toml.rs:90` | format 3 + migration |
| caller | `fit_a_cohort`, `run/census_fit.rs:220`, the call at `:262` | one call replaced |
| CLI | `str_param_estimates_at_once`, `cli/estimate_parameters.rs:145` | removed with an error |

## 5. Open items

None of design. Impl-time confirmations:
- the exact place the periods' parallel loop lives (inside `fit_repeat_tracts` or `fit_a_cohort`);
- whether `DerivedStratum` needs its provenance filled with `PenalisedCurve` points or `Curve` ones
  (a class with no reads has no reads of its own; `PenalisedCurve` with `reach` says it best).

## 6. Test & bench shape

- Unit tests next to the code, carried from the prototype (spec §4.1), plus: determinism (two runs,
  1 and 18 threads, same bits); `StratumOutcome` conversion for a class with and without reads;
  migration of a format-1 and a format-2 fixture; the new rung.
- The bench (`bench_fixtures`, `examples/ng_str_fit_bench.rs`) is the accuracy oracle: production
  settings reproduce `D_tract_prior` on the cells of `compare_D_prior_2026-10-14.md`.
- The oracle script and `cli::cross_platform_digests` are re-recorded once the move is explained.
