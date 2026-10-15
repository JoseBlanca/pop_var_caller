# The repeat-tract fit by penalised EM: types & interfaces

**Architecture draft, 2026-10-14**, companion to two approved specs:
[`../spec/str_penalised_fit.md`](../spec/str_penalised_fit.md) (the fit) and
[`../spec/str_fit_to_calling.md`](../spec/str_fit_to_calling.md) (the hand-over to calling). Every
*why* is there; this doc gives the code shape. Shared conventions:
[`module_layout.md`](module_layout.md). Line numbers are on branch `str-penalised-fit` after merging
main (`e0b41ad4`). **Amended 2026-10-14 after an adversarial review**; records that carry the lead's
decisions say so and point to the spec's §7 or the companion's §6, where the owner confirms them.
**Amended 2026-10-15, approved by the owner the same day:** the stutter prior of spec §2.10, in §1.1a
and one row of §1.4; and the stated F's median over classes with at least 8 loci with reads
(companion §3.5), in §1.5. **Amended again 2026-10-15, owner, after step C8's review:** the
slip rate joins §1.4's prior row (§1.1a); and §2.4a, the cohort put in one fixed order before the
fit.

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
  `StratumOutcome`, the provenance types). The climb's code in it is deleted at plan milestone E (§4.1
  lists what goes).
- **`calling/parameters_file/`** owns format 3, the migration and `BetweenLocusShare`;
  **`stratum_fits.rs`** owns the new ladder step.

One implementation, so **no trait**: the research bench's `PeriodFit` trait stays in
`bench_fixtures` for comparing methods and is not production surface.

---

## 1. Types

### 1.1 The production settings

The research defaults of `PenalisedFitConfig` (`penalised_fit.rs:130-155`) are not what production
runs (spec §2.7). Production gets one named constructor. The research fields that only measured
alternatives are removed (`allele_prior_fixed_draws`, `allele_prior_skips_within`,
`weigh_allele_draws_by_locus`), and with them the bench's methods that set them (`D_prior_*`,
`examples/ng_str_fit_bench.rs:191-209`) and the clustering argument of
`examples/str_fit_penalised_em.rs:226-234`; `Smoothing::Fixed` stays for the zero-penalty test only.

**The allele window moves in** *(lead, 2026-10-14; spec §7)*: the prototype reads `allele_span` from
`SsrFitConfig` (`penalised_fit.rs:331`, `:497`; `composite_fit.rs:546`, `:840`; `allele_classes_of`,
`ssr_fit.rs:2664`), the only field of it the new fit needs. It becomes `PenalisedFitConfig::allele_span:
i32`, and every function on the fit's path that took `&SsrFitConfig` for the window takes it from
there. `SsrFitConfig` is deleted at milestone E.

```rust
impl PenalisedFitConfig {
    /// **What `estimate-parameters` runs** (spec §2.7): window ±9 (`ALLELE_SPAN`, `ssr_fit.rs:950`),
    /// start 0.05 / 0.5 / 0.3, at most 2,000 rounds, stop below 1e-7, one locus's worth of allele
    /// draws plus one draw a class as the floor, F from pairs at most 10 times, standard errors on.
    pub fn production() -> Self;
}

/// One tract's worth of allele draws a class, spread over the lengths (spec §2.4). Measured.
pub const ALLELE_PRIOR_LOCI: f64 = 1.0;
/// Readings of F off pairs of samples, at most (spec §2.5). The bench's value; soft.
pub const PAIR_ALTERNATIONS: usize = 10;
/// A share's error is withheld when the class's loci with reads times the share is below this
/// (spec §7, owner 2026-10-14). Rests on two settings; soft.
pub const LOCI_FOR_A_SHARE_ERROR: f64 = 3.0;
/// **Every error of a class is withheld below this many loci with reads** (spec §2.6; confirmed
/// by the owner, 2026-10-14). Errors on thinner classes were never measured; soft.
pub const LOCI_FOR_AN_ERROR: usize = 8;
```

### 1.1a The stutter prior — approved by the owner 2026-10-15

Spec §2.10; the owner's answers are in spec §2.10.9 (weight one read, centre calling's default, no new
file key, the `NotPlaced` rule below, `RIDGE` kept). Built in plan step C8. Line numbers at commit
`47e21b8c` (step C7).

```rust
pub struct PenalisedFitConfig {
    // … unchanged …
    /// **Pseudo-reads added to each class's three stutter binomials before the curves are drawn**
    /// (spec §2.10): `stutter_prior_reads × centre` successes in `stutter_prior_reads` trials, at
    /// every repeat count with a class. Zero turns it off, which is the fit before step C8.
    pub stutter_prior_reads: f64,
    /// Where the prior sits: calling's default stutter in the fit's three numbers.
    pub stutter_prior_centre: Slippage,
}

/// One read's worth a class, each curve (spec §2.10.2; owner, 2026-10-15). Soft: C8's sweep of
/// 0, 1 and 5 on held-out loci checks it, and a different winner goes back to the owner.
pub const STUTTER_PRIOR_READS: f64 = 1.0;
/// `StutterModel::hipstr_shipped`'s whole-repeat numbers in the fit's terms: 0.05 + 0.05 of reads
/// slip, half of them shorter, one-step share 0.95 (spec §2.10.2). Not a measurement.
pub const STUTTER_PRIOR_CENTRE: Slippage = Slippage {
    level: 0.10,
    shorter_share: 0.5,
    fall_off: 0.05,
};
```

- `production()` sets `STUTTER_PRIOR_READS` and `STUTTER_PRIOR_CENTRE`; `prototype()` sets weight 0,
  so the research examples keep their numbers.
- **The centre is a copy, pinned by a test.** `parameter_estimation` imports nothing from `calling`
  (`calling/parameters_file/mod.rs:37-38`), so the constant is written out here. A test on the
  calling side checks that `stutter_rates_for(&STUTTER_PRIOR_CENTRE)` (`stutter_rates.rs:86`) gives
  `hipstr_shipped`'s whole-repeat shorter, longer and one-step shares (`alignment/stutter.rs:310-319`).
- **Where it enters the loop:** the `binomial` closure in the M-step (`penalised_fit.rs:1071-1084`)
  takes the prior for the curve it builds and adds it to `successes` and `trials` at points whose
  counts are `Some`; the three calls at `:1086-1102` pass the centre's slip rate, shorter share and
  fall-off. F's call (`:1106`) passes none. Points with no class get nothing, so the curve still
  interpolates them.
- **`refit`'s "two points or fewer" test** (`:1921`, `:1931`) counts points with trials, which the
  prior makes every class. Spec §2.10.3 says this is intended; `bends()` (classes minus two) is
  unchanged.
- **Where it enters the errors:** in `standard_errors`, after the bend penalties (`:1384-1408`), add
  `a·p(1 − p)` at each class's diagonal entry of `penalty` for components 0, 1 and 2, with `p` that
  class's slip rate, shorter share and fall-off and `a` the weight. `penalised_sandwich` (`:2301`)
  then carries it into both H and the middle term.
- **`PenalisedClass::expected_slipped_reads`** (`:319-328`) stays the reads' own count, without the
  prior. Two siblings beside it carry the trials of the direction numbers, for the error rule:
  `expected_slips_with_a_direction` (`shorter + longer`) and `expected_repeats_slipped` (`steps`),
  both from the class's last counts.
- **The error mapping** (§1.4) gains one row, after the `TooFewLoci` rows: a slip rate whose
  `expected_slipped_reads`, a shorter share whose `expected_slips_with_a_direction`, or a fall-off
  whose `expected_repeats_slipped`, is below `stutter_prior_reads` → `NotPlaced`. The slip rate's
  part is the owner's, 2026-10-15, after C8's review (spec §2.10.4): with no slipped read the
  sandwich gave a logit error of 0.65, which the `NOT_PLACED` test does not withhold. In
  `errors_of_a_class` (`penalised_fit.rs`) the slip rate goes through the same check as
  `of_a_direction`, with `fitted.expected_slipped_reads` as its count.
- **The run log** (`penalised_fit/run_log.rs`): `opening_line` (`:67`) states the weight and centre;
  `how_the_fit_ran` (`:156`) counts the motif length's classes whose direction numbers rest on the
  prior. **The file:** no new key (owner, 2026-10-15); one sentence in the stutter-provenance note
  (`to_toml.rs:343`), which is the lead's reading of that answer and which the owner may drop.

### 1.2 What the fit returns per class: the existing `StratumOutcome`

**No new outcome type.** The fit produces `StratumOutcome::Fitted` for every class with reads and
`Refused { NoSpanningReads }` for a class with none, as the climb does (`ssr_fit.rs:2691-2697`). **It
produces no `Derived`** *(lead, 2026-10-14; spec §2.2)*. Two changes to `StratumFit`
(`ssr_fit.rs:537`):

```rust
pub struct StratumFit {
    /// **`None` where F is absent** (spec §2.5). Was `f64`.
    pub concentration: Option<f64>,
    // … `slippage` (one entry, group 0: spec §1.4), the shares, errors and provenance, unchanged;
    // `samples_fitted_on` is `None` from the penalised fit, which reads every sample …
    /// **The climb's own record** (`ending`, `walks`, `log_likelihood_a_tract`, `borrowed`,
    /// `tracts_fitted`). `None` from the penalised fit; deleted with the climb (plan milestone E).
    pub climb: Option<ClimbRecord>,
}
```

*The `climb` wrapper is the seam that lets both fits build one struct until the climb is deleted.*
`StratumFits::over` reads only the fields kept above (`stratum_fits.rs:504-607`). **Outside
`ssr_fit.rs`, four examples read the climb's fields, `StratumFit { … }` literals sit at 16 sites in 9
files, and a bench calls the climb**; plan step A2 lists them.

### 1.3 Where a number came from

A value on a penalised curve is neither a cell's own fit nor main's parametric curve (companion §2.3).
One new variant in each source enum:

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

`LevelProvenance` (`ssr_fit.rs:733-743`): `curve` is `None` for a penalised point; `reach` is
`Some(CurveReach::Inside)`, since the fit reports only classes with reads and they set the range;
`slipped_reads` stays the class's expected slipped reads.

**`curve_weight()` returns `Option<f64>`, `None` for a penalised point** *(lead, 2026-10-14)*
(`slippage_curve.rs:539-545`, `share_curve.rs:592-598`). **Calling's `warrant_of`**
(`calling/inference/repeat_tract_parameters.rs:1154`) **tests the source by an exhaustive match**:
`Cell`, `Stratum` and `PenalisedCurve` are `FittedHere`; `Curve`, and `Blend` at any weight above 0,
are `Borrowed`, as today — a `Blend` at weight 0 carries only the class's own fit and stays
`FittedHere` (companion §2.3; owner, 2026-10-14).

**Errors beside a penalised point are kept and accepted** (companion §2.3):
`OwnFitStandardErrors::of_the_numbers_emitted` (`stratum_fits.rs:107-126`) and `an_own_fit_error`'s
callers (`validate.rs:810-817`, `:851-865`; the function at `:1184-1200`) treat `PenalisedCurve` as
they treat `Cell` / `Stratum` / `ThisStratum`, each by an exhaustive match so that the next variant
cannot slip past them.

### 1.4 Errors: the prototype's into `StratumError`

```rust
pub enum StratumError {         // ssr_fit.rs:632
    Estimated(f64),
    NotIdentified,
    NotPlaced,
    NoShare,
    /// **Withheld: too few loci for the error to be honest** — a class below
    /// `LOCI_FOR_AN_ERROR` loci with reads, or a share below `LOCI_FOR_A_SHARE_ERROR` loci's
    /// worth (spec §2.6, §7).
    TooFewLoci,
}
```

`StratumErrors` (`ssr_fit.rs:611-621`) is unchanged: **F's error travels as the concentration's**
(`concentration`, from `PenalisedClassErrors::concentration`, `penalised_fit.rs:236`) and is
converted only at the file's edge (§2.3) *(lead, 2026-10-14)*. The prototype's F-scale error
(`PenalisedClassErrors::between_locus_share`, `:1046`) is not used.

**The mapping** *(lead, 2026-10-14)*, applied in this order to each number of a class:

| case | `StratumError` |
|---|---|
| F absent (its error, `concentration`) | `NotPlaced` |
| class with fewer than `LOCI_FOR_AN_ERROR` loci with reads | `TooFewLoci` |
| a share with share × loci with reads < `LOCI_FOR_A_SHARE_ERROR` | `TooFewLoci` |
| *(owner, 2026-10-15, §1.1a)* a slip rate whose class's expected slipped reads, a shorter share whose slipped reads with a direction, or a fall-off whose repeats slipped, are below `stutter_prior_reads` | `NotPlaced` |
| an error wider than `NOT_PLACED` (3, `ssr_fit.rs:683`) on the scale the number is fitted on — logit, or log for the concentration — as the climb judges it (`ssr_fit.rs:2172`) | `NotPlaced` |
| a positive, finite error | `Estimated` |
| zero, absent or not finite | `NotPlaced` |

**Never `Estimated(0)`**: the prototype clamps negative variances to zero (`penalised_fit.rs:1002`,
`:1015`, `:1037`), and the validator refuses an error of zero (`validate.rs:1154-1166`).

### 1.5 Calling's tables (companion §3)

```rust
pub struct FittedAlleleLengthVariation {      // stratum_fits.rs:379
    pub weights: Vec<f64>,
    /// **`None` where the fit reported no F** (companion §3.2). Was `f64`.
    pub concentration: Option<f64>,
    /// The concentration's standard error; `None` where it has none (companion §2.1).
    pub concentration_standard_error: Option<f64>,
}

/// **Declaration order is rank, weakest last**: the enum derives `Ord` (`stratum_fits.rs:358`)
/// and `BundlePrior::rung` reports a bundle's `.max()` as its weakest (`seed_ssr.rs:374-380`).
pub enum AlleleLengthVariationRung {          // stratum_fits.rs:359
    StratumsOwnFit,
    /// The class's own allele length variation, held at the stated concentration because its F
    /// is absent (companion §3.2, owner 2026-10-14).
    StratumsOwnSharesStatedSpread,
    PeriodsPooledTracts,
    StatedFlat,
}
```

- `FittedFrom` (`:190`) gains `ThisStratumsSharesStatedSpread`, between its two variants, which
  `AlleleLengthVariation::fitted` (`:257-279`) maps to the new rung.
- `allele_length_variation_at` (`:781-805`): a class with shares and no concentration gets
  `fitted(weights, stated, ThisStratumsSharesStatedSpread)`. `checked_allele_length_variation`
  (`:1022`) checks the concentration only when it is `Some`.
- `fitted_allele_length_variation_of_each_stratum` (`:954-960`) yields the optional concentration and
  its error to the writer (`from_run_parameters.rs:696-713`).
- `median_concentration` (`:1100-1112`) becomes the median taken in F, converted once; its warrant
  (`:597-601`) is `Defaulted` when no class has an F (companion §3.5).
- **The median counts only classes with at least `LOCI_FOR_AN_ERROR` loci with reads** (owner,
  2026-10-15; companion §3.5; plan step C9). Line numbers here are at the working tree of
  2026-10-15. The median is computed by `median_concentration` (`stratum_fits.rs:1250`), called once
  from `StratumFits::over` (`:657`). Today it reads only the map of `FittedAlleleLengthVariation`s
  (`:431-443`), which holds no locus count, so it cannot apply the rule by itself.
  - **The count is `StratumFit::tracts_of_its_own`** (`ssr_fit.rs:577`): the class's loci with at
    least one spanning read. The penalised fit sets it from `evidence.tracts_with_reads()`
    (`penalised_fit.rs:796`, stored at `:810`; the method at `ssr_fit.rs:405-407`), the same count
    the error rule compares with `LOCI_FOR_AN_ERROR` (`run_log.rs:473`). It is not in the fit's
    errors or in `FittedAlleleLengthVariation`.
  - **Where to read it:** `over` already matches `StratumOutcome::Fitted(fit)` to harvest each
    class's allele length variation (`:639-655`); it collects `fit.tracts_of_its_own` beside each
    concentration there and hands both to `median_concentration`, which takes the classes with an
    F and at least `LOCI_FOR_AN_ERROR` loci, or every class with an F if there are none. The
    constant is `pub` in `ssr_fit/penalised_fit.rs:260`, reached as
    `super::ssr_fit::penalised_fit::LOCI_FOR_AN_ERROR` (`pub mod penalised_fit`, `ssr_fit.rs:78`).
    `FittedAlleleLengthVariation` does not gain the count: the path that reads a file
    (`StratumFits::of_gathered_rows`) takes the median as written and never recomputes it.
  - **The warrant code does not change** (`:664-674`): `FittedHere` when any class has an F, which
    is exactly when the median is taken from some class's F under the new rule.

### 1.6 The parameters file, format 3

```rust
pub const FORMAT_VERSION: u32 = 3;                     // parameters_file/mod.rs:192

/// **F, checked: strictly between 0 and 1, and its concentration (1 − F)/F finite.** The file's
/// rows and its edge hold this, never a bare f64, so F and a concentration cannot be swapped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetweenLocusShare(f64);
impl BetweenLocusShare {
    pub fn new(share: f64) -> Option<Self>;
    pub fn get(self) -> f64;
}

pub struct StratumAlleleLengthVariationRow {           // mod.rs:1354
    pub period: u8,
    pub reference_repeats: u64,
    /// Absent where the fit's F is.
    pub between_locus_share: Option<BetweenLocusShare>,
    pub between_locus_share_standard_error: Option<f64>,
    pub share_of_chromosomes_by_repeat_offset: Vec<f64>,
}
// PeriodAlleleLengthVariationRow (mod.rs:1374): the same renames.

pub struct RepeatTracts {
    /// Was `fallback_length_spectrum_concentration` (mod.rs:1065); now F, checked with
    /// `BetweenLocusShare::new` by validation.
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

Validation (`validate.rs:1364-1396`): a `PenalisedCurve`'s `strength` is finite and positive. How a
row's `BetweenLocusShare` is deserialised through its checked constructor (`serde(try_from = "f64")`
or a check in `validate`) is an impl-time choice; either way no unchecked F reaches calling.

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
  once per period, **periods one at a time** *(lead, 2026-10-14; spec §2.2)*.
- **Memory: one class at a time.** Each round visits a period's classes one at a time, parallel only
  inside a class (its distinct rows and loci): **no parallel loop over classes or periods**, so at
  most one class's distinct-row table exists at any thread count. The two parallel loops over classes
  today become sequential: step 1 (`penalised_fit.rs:588`) and the error blocks (`:852`). Why: spec
  §3; the climb's lesson, `ssr_fit.rs:2736-2740`.
- **Deterministic:** the same bits on every run and at every thread count. Every sum over loci,
  samples or classes is taken in index order, in every module the fit calls, including
  `SamplePairs::log_likelihood` (`composite_fit.rs:1583-1590`) and `fit_spread_by_pairs`'s objective
  (`:1217-1222`) (spec §3, first trap).
- A class with reads → `Fitted`, provenance `PenalisedCurve`, `samples_fitted_on: None`; with none →
  `Refused { NoSpanningReads }`. No refusal floor (spec §2.2).
- `concentration` is `None` when no locus of the period has two samples with reads, or when the
  estimate is within 1% of an end of its search (`concentration_is_placed`'s rule,
  `composite_fit.rs:83-98`; spec §2.5). `working_concentration` is never written into a `StratumFit`.
- **Asserts one slippage group** in the evidence: the fit pools libraries
  (`one_stutter_for_every_library`, `composite_fit.rs:990`), and `StratumFit::slippage` has one entry.
- Infallible: no input it accepts makes it fail; malformed evidence is a caller bug and panics as
  `gather_strata`'s outputs already do.

### 2.2 Where it is called

`fit_a_cohort` (`run/census_fit.rs:220`) replaces
`ssr_fit::fit_strata_on_sample_subsets(&evidence, &homozygote_excess, &order, tracts)` (`:262`) with
`ssr_fit::penalised_fit::fit_repeat_tracts(&evidence, &homozygote_excess, tracts)`, its parameter
`tracts: &SsrFitConfig` (`:226`) becoming `tracts: &PenalisedFitConfig`. `sample_order` (`:259-260`)
is no longer needed. The CLI's `repeat_tract_config` (`cli/estimate_parameters.rs:362-366`) returns
`PenalisedFitConfig::production()`. The `CohortFit` returned is unchanged.

### 2.3 The file's edge

```rust
/// F = 1 / (1 + concentration); `None` for a concentration that is not finite and positive, so
/// small that F rounds to 1, or so large (about 1.8×10³⁰⁸ and above) that F's own concentration
/// overflows to infinity (owner, 2026-10-14). Only here does the formula appear (companion §2.1).
pub fn between_locus_share_of(concentration: f64) -> Option<BetweenLocusShare>;
/// concentration = (1 − F) / F; finite by the newtype's invariant.
pub fn concentration_of(share: BetweenLocusShare) -> f64;
/// se(F) = se(c) / (1 + c)², the delta method (companion §2.1; lead, 2026-10-14).
pub fn between_locus_share_error_of(concentration: f64, concentration_error: f64) -> f64;

/// **Rewrite a format-1 or format-2 file's parsed TOML into format 3** (companion §2.2, every rule
/// there); an old concentration `between_locus_share_of` refuses is refused with its key path.
fn migrate_to_format_3(file: &mut toml::Value) -> Result<(), ParametersFileError>;
```

`from_toml` (`from_toml.rs:90`) parses to `toml::Value` only to read `format_version`. **Version 3:
the original text is deserialised, errors keep their line** (`:90-101`); **versions 1 and 2 are
migrated and deserialised from the `Value`, errors name the key path.** All three conversions live in
`calling/parameters_file/`; nothing else in the crate computes F.

### 2.4 The CLI

`--str-param-estimates-at-once` (`cli/estimate_parameters.rs:145`) becomes a hidden argument that,
if given, exits with an error saying the repeat-tract fit no longer holds per-stratum tables (owner,
spec §7). Its doc and the `strata_at_once` plumbing (`:362-366`) go. Detecting "was given" changes
its type, which touches the `EstimateParametersArgs` literals in `cli/cross_platform_digests.rs:371`
and `cli/estimate_parameters/tests.rs:98`.

### 2.4a One fixed order for the cohort — owner, 2026-10-15

Spec §3 (the second trap); plan step C10. `OpenPspCohort::open` (`run/psp_caller.rs`) keeps its psps
"in the order the paths were given", which becomes the run's sample order, and numbers every read
group run-wide in that order; the fit's sums over samples and read groups follow that numbering. So
`open_the_censuses` (`cli/estimate_parameters.rs`) puts the cohort into one fixed order **before
anything is summed**: samples sorted by the sample name in each psp's header, and read groups within
a sample by their declared ID (`@RG ID`). The sort lives in `open_the_censuses`, not in the opener,
because the opener is shared with `call-from-psps`, whose run order is the order its paths were
given. The written file then lists samples in name order, which `parameters_file.md` §6 allows:
calling joins the file's rows by sample name and `@RG ID`, not by position. The defect is on main
too; this change reaches main with the branch.

---

## 3. Design decisions — decided

- **One implementation, no trait.** Nothing competes with it in production (spec §1).
- **`StratumOutcome` is reused, not replaced**, because `StratumFits::over` already reads it
  (spec §2.1).
- **No `Derived` from the new fit** (spec §2.2; lead). This closes the former open item: every point
  is `PenalisedCurve`, since a `Curve` source with no curve panics in `the_recorded_curve`
  (`from_run_parameters.rs:854-868`).
- **The climb's record is wrapped, then deleted**, so both fits share `StratumFit` until milestone E.
- **New provenance variants, not reuse of `ThisStratum`** (companion §2.3).
- **`BetweenLocusShare` newtype** *(lead, 2026-10-14)*: F = 0.2 read as a concentration would pass
  every check of a bare `f64`.
- **Migration on the parsed TOML, not serde aliases**: the concentration's value changes, not just
  its key (companion §2.2).
- **The bundle seed uses the edge rule run by run** (`run_lengths_at_edge`, `seed_ssr.rs`; companion
  §3.3).

## 4. Reconciliation with existing code

| new / changed | existing code | how |
|---|---|---|
| `fit_repeat_tracts` | `fit_period_by_penalised_counting`, `penalised_fit.rs:328` | wraps it per period; adds the `StratumOutcome` conversion |
| production settings | `PenalisedFitConfig::default`, `penalised_fit.rs:130-155`; bench overrides `fit_bench.rs:497-515`; `ALLELE_SPAN`, `ssr_fit.rs:950` | `production()` holds the bench's values, `ALLELE_PRIOR_LOCI` and the window |
| step 1 and pairs | `expected_counts`, `ClassCounts`, `SamplePairs` in `composite_fit.rs` | moved to `expected_counts.rs` unchanged |
| F's placement | `concentration_is_placed`, `composite_fit.rs:83-98` | reused |
| outcome | `StratumOutcome`, `ssr_fit.rs:854`; `StratumFit`, `:537` | reused; `concentration: Option`, `climb: Option<ClimbRecord>` |
| provenance | `LevelProvenance`, `ssr_fit.rs:733`; `LevelSource`, `slippage_curve.rs:519`; `ShareSource`, `share_curve.rs:580`; `CurveReach`, `slippage_curve.rs:654` | one variant each; `curve_weight()` optional |
| error sites | `of_the_numbers_emitted`, `stratum_fits.rs:107`; `an_own_fit_error`, `validate.rs:1184`; `warrant_of`, `repeat_tract_parameters.rs:1154` | exhaustive matches accepting `PenalisedCurve` |
| errors | `StratumErrors`, `ssr_fit.rs:611`; `StratumError`, `:632` | one variant; mapping §1.4 |
| calling's tables | `FittedAlleleLengthVariation`, `stratum_fits.rs:379`; rung `:359`; `FittedFrom` `:190`; `allele_length_variation_at` `:781`; `checked_allele_length_variation` `:1022`; `median_concentration` `:1100` | optional concentration and its error; new rung; median in F |
| file | `FORMAT_VERSION`, `mod.rs:192`; rows `:1354`, `:1374`; fallback `:1065`; `LevelSmoothing` `:1199`; `ShareSmoothing` `:1226`; `from_toml`, `from_toml.rs:90` | format 3 + migration + `BetweenLocusShare` |
| provenance in the file | `level_origin_of`, `from_run_parameters.rs:799`; `share_smoothing_of`, `:830`; `a_level_provenance`, `to_run_parameters.rs:609`; `a_share_provenance`, `:644`; `a_level_smoothing`, `to_toml.rs:1281`; `validate.rs:1364-1396` | the `PenalisedCurve` arm in each |
| caller | `fit_a_cohort`, `run/census_fit.rs:220`, the call at `:262` | one call and one parameter type replaced |
| CLI | `str_param_estimates_at_once`, `cli/estimate_parameters.rs:145` | removed with an error |

### 4.1 What goes with the climb (plan step E1) — lead, 2026-10-14

| item | fate |
|---|---|
| the climb (`fit_strata*`, `fit_stratum`, the quadrature, the walks, `SampleSubsets`, `draw_and_apply_curves`'s blend, `ClimbRecord`) | deleted (owner, spec §7) |
| `fit_period_allele_length_variations`, `ssr_fit.rs:3392` (only producer of the pooled rung; test callers only) | deleted; its tests build period tables by hand; the rung stays readable (companion §3.4) |
| `ssr_fit/library_stutter.rs` (uses `dirichlet_points`, `TractLikelihoods`, `:71`, `:107`) | deleted; the library-multiplier work rebuilds on the new fit (spec §1.4) |
| `ssr_fit/genotype_rows.rs:55`, `ssr_fit/neighbour_start.rs:55, 91` (measurement helpers on `TractLikelihoods`, `climb_from_parameters`, `Scorer`) and their examples | deleted with the climb they call |
| bench methods `A` (`fit_bench.rs:284`) and `B_curves` (`:660`) | deleted; their results stay in the research reports |
| `SsrFitConfig`, `ssr_fit.rs:1046` | deleted; `allele_span` lives in `PenalisedFitConfig` (§1.1). Examples that only name it for the window are ported |
| `StratumRefusal::BelowTheFloor`, `ssr_fit.rs:849` (only the climb produces it) | deleted; tests that build one use `NoSpanningReads` |
| `SlippageCurve`, `ShareCurve`, `a_slippage_curve` (`to_run_parameters.rs:612-631`), `LevelSource::Curve`/`Blend`, `curve_weight()` | kept: older files' `ThisPeriodsCurve`/`Blend` still read |

## 5. Open items

None of design. Impl-time confirmation: how a row's `BetweenLocusShare` is deserialised (§1.6).

## 6. Test & bench shape

- Unit tests next to the code, carried from the prototype (spec §4.1), plus: determinism (1 and 18
  threads); the `StratumOutcome` conversion (with and without reads, F at a search end); the error
  mapping of §1.4; migration (format 1, both format-2 spellings); the new rung and its rank; the
  stated F's warrant at one sample; the stated F's median skipping classes under 8 loci with reads,
  and falling back to every class with an F when none has 8 (plan C9); a fitted file round-tripping
  with its errors present.
- Fixtures and matches that change are listed in the plan's steps (A2, A4, B2, B3).
- The bench (`bench_fixtures`, `examples/ng_str_fit_bench.rs`) is the accuracy oracle: production
  settings reproduce `D_tract_prior` on the cells of `compare_D_prior_2026-10-14.md`, against a
  fixture of the prototype's numbers recorded before the code moves.
- The oracle script and `cli::cross_platform_digests` are re-recorded once the move is explained.
