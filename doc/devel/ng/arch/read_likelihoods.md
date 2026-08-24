# ng — the read likelihood models: types & interfaces

*Status: architecture draft (2026-08-21), companion to the spec
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) (the design and its rationale) and to
the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) and
[`module_layout.md`](module_layout.md). One of three coordinated calling arch docs — the siblings
are [`calling_priors.md`](calling_priors.md) and [`calling_em_loop.md`](calling_em_loop.md); §0
says which doc owns each shared type. `ng_step_interfaces.md` §3 sketches step 7 as
`read_log_lik(&self, read, allele, params)`; that sketch is superseded here (§1.2). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md); **STR** in prose ↔
`ssr` in code. Signatures are illustrative; the **contract** is the deliverable. The spec carries
every why.*

## 0. Who owns what — the three calling docs

| shared thing | owner |
|---|---|
| the per-sample evidence views (`GenericSampleEvidence`, `SsrSampleEvidence`), the `Lg` **row** contract, the stutter distribution and the scoring contexts | **this doc** (§2, §3, §4) |
| `Concentration`, the seed builders, the per-sample concentration | [`calling_priors.md`](calling_priors.md) §2–§3 |
| `CandidateAlleles`, `GenotypeTable` + `AlleleId`, `ExpectedAlleleCopies`, `CallingScratch` (which *holds* the `Lg` table), `Phred` | [`calling_em_loop.md`](calling_em_loop.md) §2 |
| `LogProb`, `ErrorRate`, `Ploidy`, `Motif`, `SsrPeriod` | exist in [`src/ng/types.rs`](../../../../src/ng/types.rs) — reconciled, not re-minted |

Split worth stating twice: the likelihood owns the **row** — how one sample's `Lg` values are
computed — and the loop owns the **table** that stores rows for all samples. One type each, no
overlap.

## Module home

`src/ng/calling/likelihood/` — step 7, inside the `calling/` folder that holds steps 6–9
([`module_layout.md`](module_layout.md)). The only swappable seam is the STR emission (the SNP/indel
closed form has no competing implementation, so it is a file with no trait ceremony —
module_layout principle 1a):

```
src/ng/calling/likelihood/
├── mod.rs          – the evidence views, the row contract, the parameter tiers, shared floors
├── generic.rs      – the SNP/indel closed form + contamination mixture + the partial rule
├── ssr.rs          – the STR row: emission cache, three-term mixture, censored tail
└── ssr_emission.rs – SsrEmissionModel trait + StutterSubstitutionEmission (Model A)
                      + ClassicEmissionOracle (Model B, test-only)
```

The stutter distribution itself is **reused from `src/ng/alignment/stutter.rs`**, not duplicated
(§4.2) — the spec fixes ownership of the *definition* here (spec §7) while the built type stays
where both consumers reach it.

## 1. General — what both paths share

### 1.1 The row contract

One call computes **one sample's whole `Lg` row** — one `LogProb` per candidate genotype, parallel
to the loop's `GenotypeTable`. Contract, both paths:

- pure function of (evidence, candidates, frozen parameters); bit-identical at any thread count,
  any observation order sorted as the merge sorts (spec §8);
- fills caller scratch, allocates nothing per sample (spec §8). The scratch shapes
  (`GenericRowScratch`, `SsrRowScratch`) are this doc's; they live as one section inside the
  loop's `CallingScratch` ([`calling_em_loop.md`](calling_em_loop.md) §2);
- an empty evidence row is all-zeros — the prior decides, no branch (spec §3.3);
- a mis-shaped input (row length ≠ genotype count, a candidate with no stratum entry, a non-finite
  parameter) is a caller bug → assertion, structural ones held in release
  ([`per_group_merger.rs:1963`](../../../../src/var_calling/per_group_merger.rs) is the precedent);
- every probability is floored before a logarithm: `MIN_BASE_ERROR = 1e-12`
  ([`contamination_estimation.rs:1449`](../../../../src/var_calling/contamination_estimation.rs)),
  geometric clamps `(0.01, 0.99)`
  ([`alignment/stutter.rs:74`](../../../../src/ng/alignment/stutter.rs)) — imported as named
  constants with their reasons.

### 1.2 Correction to `ng_step_interfaces.md` §3 step 7 — recorded, not applied there

`read_log_lik(&self, read: &LocusRead, allele: &Allele, params) -> LogProb` is superseded twice.
**There are no reads here**: the evidence is a distinct observation with a count and summed
moments (spec §1.4), so the unit argument is an observation, never a `LocusRead` (that type is
itself retired — `ng_step_interfaces.md` §6). **And per-allele scalar returns cannot compose the
formula**: spec §2.1 puts a logarithm around a sum over alleles (the junk/contamination mixture),
so the outer call must return a whole genotype row (why: spec §7's own correction). Consequence
for the recipe: the swappable `Box<dyn ReadLikelihood>` becomes the narrower STR emission seam of
§4.1 — the SNP/indel row has no impl to swap.

### 1.3 The parameter tiers — the middle tier is switchable without touching the model

Spec §6.1's three tiers, as code obligations:

| tier | parameters | code shape |
|---|---|---|
| frozen for the run | error rate + calibration scale, STR substitution rate, the contamination **fraction** | fields of the run-level parameter views (§2.3, §4.1) — nothing downstream may write them |
| **per-locus re-estimable** (off by default) | slippage level, direction split, fall-off | arrive **per call** inside `SsrScoringContext` (§4.1); the emission never asks where they came from. This is the seam that lets the EM loop re-fit them ([`calling_em_loop.md`](calling_em_loop.md) §6.1) with zero changes here — the one constraint spec §6.1 makes binding |
| re-estimated every pass | per-locus allele frequencies | mostly the prior's, and **one term here reads them too**: the contaminating population's frequency for the allele an observation shows (see below) |

**Correction, 2026-08-24, made while A2 was built: the third row used to say "invisible here … no
term of §2.1 reads them", and spec §3.6's correction of the same day makes that false.** The
contamination mixture's second half is the frequency of the observation's own allele at the locus
being called, over the samples in that sample's sequencing batch, recomputed every iteration — so
one term of the SNP/indel row does read a per-locus frequency. The row above it is corrected in
the same breath: it is the contamination **fraction** that is frozen, not contamination as a whole.

**The two halves therefore sit in different tiers, and nothing about that reopens the ruling that
contamination never enters the loop.** That ruling is about the fraction, which is a property of a
library and of nothing else. The second half is a property of the locus, and a per-locus quantity
is what a per-locus loop is for. What it costs is a lookup rather than a fit, because it is the
same number the genotype prior already reads.

### 1.4 Provenance

The model never branches on provenance; it **propagates** it: each scoring context carries the
weakest `Provenance` ([`parameter_estimation/mod.rs:60`](../../../../src/ng/parameter_estimation/mod.rs))
among the parameters that entered it, and the loop copies the per-locus weakest onto
`LocusInference` (spec §8; [`calling_em_loop.md`](calling_em_loop.md) §2).

## 2. Types — the evidence views

### 2.1 Generic path

```rust
/// One sample's evidence at one cohort locus, as the SNP/indel row consumes it.
/// A view over the merge's SampleSupport — see the reconciliation note below.
pub struct GenericSampleEvidence<'a> {
    /// One entry per (allele, read group) this sample's complete reads showed.
    pub supported: &'a [GenericObservation],
    /// q_sum pooled over reads matching NO candidate allele — the support of the table
    /// alleles that candidate selection dropped, folded here by that selection. Cancels
    /// in genotyping, kept for the data likelihood (spec §3.3's q_sum_other).
    pub unmatched_q_sum: f64,
    /// The partial observations, bases + witnessed positions intact — §3's compatibility
    /// rule (spec §5.3) needs them, so they are NOT folded onto alleles. **A set of runs
    /// with holes in it, not one run** (spec §5.3, corrected 2026-08-24), which makes the
    /// restricted projection a gather rather than a subslice.
    pub partials: &'a [PartialObservation<'a>],
}

/// The merge's fold of every read that showed one allele from one read group.
pub struct GenericObservation {
    pub allele: AlleleId,          // a(o) — the CANDIDATE index, not the merge's own
    pub read_group: ReadGroupId,   // part of the identity — spec §2.3's aggregation contract
    pub num_reads: u32,            // n_o
    pub q_sum: f64,                // Σ ln P(error) over those reads
}
```

**Reconciliation note — `partials` now exists on the merge's own row, owned rather than
borrowed** (`calling_prerequisites.md` C1, 2026-08-23). `SampleSupport` carries
`partials: Vec<PartialObservation>` — the witnessed positions **in the cohort locus's
coordinates**, the bases as the mint recorded them, the read group, a read count and a quality
sum. The sketch above types it `PartialObservation<'a>`; it cannot borrow, because the
observations the bases come from live in the `ObservationCache` and the organiser drops that
ground once a locus is released (`arch/cohort_merge.md` §4) while the built locus is still in the
caller's hands. So the view's `&'a [PartialObservation]` borrows the merge's owned rows and there
is one type rather than two. **The field is empty until C2 routes partials into it.**

**Reconciliation note — the read-group requirement this placed on the merge has been met.** The built
`SampleSupport`/`AlleleSupport`
([`cohort_merge/build.rs`](../../../../src/ng/run/cohort_merge/build.rs)) pooled read groups into
one row per allele, and its own doc booked the split as owed to whoever brought the STR path
through. Spec §2.3 made it a **generic-path** requirement too — summing must stop at the
read-group boundary — and it landed there first:
[`calling_prerequisites.md`](../impl_plan/calling_prerequisites.md) B1, 2026-08-23.
`SupportedAllele` now carries `read_group`, so `SampleSupport::supported` is one row per
`(allele, read group)` in ascending pair order, folding to today's shape where a sample has one
group. **What a consumer must not do is add the rows back**: `SampleSupport::pooled_support_for`
exists for the questions that really are about the sample, and its name is the warning.

**Correction, 2026-08-24, made while A1 was built: `allele` above cannot be filled from the merge's
row, and the deleted comment said it could.** `SupportedAllele::allele` indexes the merge's
unification table — every distinct sequence the whole cohort showed, uncapped — and `AlleleId`
indexes the **candidate** table, which selection produces by keeping some of those alleles and
dropping the rest. Dropping allele *k* renumbers every allele above it, as `CandidateAlleles`'
own doc states, so the two numberings agree only until the first prune, and a view that assumed
the identity would score reads against the wrong sequence with nothing saying so.

Two consequences, both now in the built code. The narrowing is **not** this module's:
`GenericObservation::of_supported_allele` takes the candidate id as an argument, and
`fill_from_supported_alleles` takes selection's whole mapping — `&[Option<AlleleId>]` indexed by
merge allele — so the assumption is a parameter rather than a silence. And a row the mapping
drops is not merely skipped: its `q_sum` comes back from the fill, because those reads are part of
`unmatched_q_sum` (spec §3.3) and a function that discards evidence has to say what it discarded.
**Selection still owes the rest of the pool**; the fill returns only the part it can see.

### 2.2 STR path

The evidence is the locus generator's own type, unchanged: `SequenceObservation`
([`locus_generation/mod.rs:295`](../../../../src/ng/locus_generation/mod.rs)) already keys on
`(bases, witness, read_group)` and carries `num_obs` — the aggregation contract holds by
construction. `SsrSampleEvidence<'a>` is a slice of them plus the locus's `SsrDetail`
([`mod.rs:438`](../../../../src/ng/locus_generation/mod.rs)); complete and partial observations are
told apart by `ReadWitness`.

**Correction, 2026-08-24, made while A1 was built. This sentence used to say the generator's
`complete_observations()` ([`mod.rs:134`](../../../../src/ng/locus_generation/mod.rs)) "stays the
only unguarded access", and neither half of that is true.** It is not a guard — the field it reads,
`SampleLocusObservations::observations`, is `pub`, so the iterator is a helpful name and not an
enforcement — and it is no longer the only one: `SsrSampleEvidence` holds a bare slice, so A1
spells the same split again as `complete_observations()` and `partial_observations()` on the view.
What the pair actually buys is that scoring a partial as complete has to be written rather than
fallen into, and that the split has one spelling per type instead of one per caller.

Two shapes those two methods carry that this sketch did not ask for, both earned in review. They
yield **`(position in the slice, observation)`**, because the STR row's emission cache is keyed by
that position (§4.1) and an iterator yielding only the observation would leave the row re-walking
the unguarded field to recover it. And the partial side is an **exhaustive match** on `ReadWitness`
rather than `!= Complete`, so a third variant is a compile error at the one place that decides what
reaches the censored term.

### 2.3 Frozen per-read-group parameters (both paths)

```rust
/// §3.2's calibration: one multiplier per read group so the mean minted per-read error
/// equals the pre-pass's fitted rate. scale == 1.0 (and Provenance::Defaulted) where the
/// pre-pass emitted no rate — visible in the run's output, never silent (spec §3.2).
pub struct ReadGroupCalibration { pub scale: f64, pub provenance: Provenance }

/// §3.6's mixture inputs — per read group, frozen. `None` = no estimate exists (one
/// sample) — absent, not a fitted zero (spec §3.6).
pub struct ContaminationView {
    pub fraction: f64,                          // c; 0.0 where the fit emitted none
    pub markers_with_reads: u64,                // the evidence counts that tell
    pub reads_on_markers: u64,                  //   "measured clean" from "unmeasurable" —
                                                //   both names as ContaminationEstimate has them
}
```

**The three allele-class frequencies this used to carry are deleted (owner, 2026-08-24).** The
mixture's second half is `q(o)`, the contaminating population's frequency of the allele the
observation shows *at this locus*, and that is the loop's own estimate over the samples in this
sample's sequencing batch — not a triple frozen before calling (spec §3.6). So it does not travel on
this view, which holds only frozen per-library facts: it is read where the frequency is, each
iteration. **The parameter pre-pass owes nothing for it**, and the side-pass that was to emit it is
deleted from `impl_plan/calling_prerequisites.md` rather than left blocked.

## 3. The SNP/indel path

One concrete row function — no trait (no bake-off):

```rust
/// Spec §3.3's closed form; §3.6's mixture when any fraction is above zero — the two are
/// the same algebra at c = 0 to a few ulp, so there is no c == 0 branch (spec §3.6).
pub fn genotype_log_likelihood_row(
    evidence: &GenericSampleEvidence<'_>,
    genotypes: &GenotypeTableView<'_>,          // flat views from the loop's table
    calibration: &[ReadGroupCalibration],       // indexed by ReadGroupId
    contamination: &[ContaminationView],        // indexed by ReadGroupId; empty = §3.3 exactly
    error_spread_divisors: &[f64],              // m(a, g), precomputed per (allele, genotype)
    out: &mut [LogProb],
    scratch: &mut GenericRowScratch,
);

/// m(a, g): 3.0 where the observation differs from every allele the genotype carries by
/// a substitution at exactly one position, 1.0 otherwise. A property of the allele pair,
/// computed once per locus over the projected sequences (spec §3.5).
pub fn error_spread_divisors(alleles: &CandidateAlleles, genotypes: &GenotypeTableView<'_>, out: &mut [f64]);
```

**Contract.** No multinomial coefficient (spec §3.4 — a genotype-changing decision, measured by
the change measurement below, not asserted). A read the genotype explains is charged `log(k_a/P)`
only; a read it cannot explain is charged `q_sum + n·(log scale − log m)`; the unmatched pool is
added as a genotype-independent constant and kept for emission (spec §3.3). **Partials** enter by
spec §5.3's compatibility rule: an allele is compatible when its projection restricted to the witnessed
run equals the partial's bases; a partial compatible with none is charged as an error with
`m = 1` (spec §5.3). The aggregation identity — reads looped individually versus the merge's
aggregate, bit for bit — is spec §12 test 9 and is the reason the formula has this shape.

**What asks something of the pre-pass:** the calibration's numerator/denominator accumulator (two
scalars per read group, per surviving route) — recorded there, consumed here as
`ReadGroupErrorRateFit`
([`generic/read_group_error_rate.rs:45`](../../../../src/ng/parameter_estimation/generic/read_group_error_rate.rs));
the minted-quantity function must be the same function the locus generator mints with (spec §3.2;
test 10).

## 4. The STR path

### 4.1 The emission seam — the one swappable surface

```rust
/// Per-call context — the tier-two seam (§1.3): every number arrives here, none is read
/// from global state. Mirrors production's ReadScoringContext
/// (ssr/cohort/read_model/mod.rs:45).
pub struct SsrScoringContext<'a> {
    pub motif: &'a Motif,
    /// Built per (read group, CANDIDATE stratum) — never hoisted out of the candidate
    /// loop; the lookup is hoisted, the values are not (spec §4.4).
    pub stutter: &'a StutterModel,
    pub substitution_rate: ErrorRate,     // per read group × stratum — never the SNP ε (spec §4.3)
    /// Mass the two slip cutoffs discarded for THIS candidate — computed and reported,
    /// never assumed negligible (spec §4.2; test 5).
    pub truncated_mass_lost: f64,
    pub weakest_provenance: Provenance,
}

/// One candidate as the emission sees it: the bases, and the repeat count that keys the
/// stratum lookup (spec §4.4). Built once per locus from CandidateAlleles.
pub struct SsrCandidate<'a> { pub bases: &'a [u8], pub repeat_count: u32 }

/// Lr(observation | one candidate allele) — the only part that differs between models.
pub trait SsrEmissionModel {
    type Scratch: Default;
    /// Probability space, floored; the row takes one log per observation per genotype.
    fn emission(&self, observation: &[u8], candidate: &SsrCandidate<'_>,
                context: &SsrScoringContext<'_>, scratch: &mut Self::Scratch) -> f64;
    /// P(the read saw ≥ this much | candidate) × the letter match — the censored term
    /// (spec §5.2): factorised form on pure candidates, exact sum on interrupted ones.
    fn censored_emission(&self, witnessed_prefix: &[u8], candidate: &SsrCandidate<'_>,
                         context: &SsrScoringContext<'_>, scratch: &mut Self::Scratch) -> f64;
}

pub struct StutterSubstitutionEmission;   // Model A — the default (spec §4.1)
pub struct ClassicEmissionOracle;         // Model B — #[cfg(test)] only, the independent oracle
```

The row wraps the seam, so §4.5.1's on/off measurement and the model bake-off are one build:

```rust
pub fn genotype_log_likelihood_row<Model: SsrEmissionModel>(
    evidence: &SsrSampleEvidence<'_>, candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>, model: &Model,
    contexts: &[SsrScoringContext<'_>],           // one per (read group, candidate)
    outlier_weight: f64, reachable_length_count: u32,
    contamination: Option<SsrContamination<'_>>,  // None ⇒ the two-term mixture
    out: &mut [LogProb], scratch: &mut SsrRowScratch<Model::Scratch>,
);
```

**Contract.** Emissions are cached per `(observation, candidate)` and reused across genotypes —
the cost is `observations × candidates`, not `× genotypes`, and that is the design, not an
optimisation (spec §8). The three-term mixture is §4.5.1's:
`(1 − λ − c)·copy-mixture + λ·uniform + c·seed(o)`; with `contamination: None` it is the two-term
form. `λ` is `DEFAULT_OUTLIER_WEIGHT = 0.01` — **inherited from production and declared
inherited**, not fitted ([`em.rs:138`](../../../../src/ssr/cohort/em.rs) sets it) — and it is
spread over `reachable_length_count`, **a property of the candidate set and the cutoffs, computed
per locus with no cohort in it** (spec §4.5's decided repair of production's cohort-wide `D`,
[`em.rs:393`](../../../../src/ssr/cohort/em.rs)).

```rust
/// §4.5.1's second half: the fraction (frozen, per read group) and the prior's seed
/// shape normalised to sum to one — built once per locus by the loop from
/// genotype_prior::fill_seed_share_per_candidate (calling_priors.md §5), frozen
/// thereafter.
///
/// OPEN: that builder fills one entry per CANDIDATE, and `c · seed(o)` asks for a
/// probability per observed LENGTH. The two differ wherever a locus has more candidates
/// than lengths (two spellings of one length each take the rung's full share), wherever a
/// read lands at a length the prune dropped, and for a censored read, which has no length
/// at all. This field has to say which support it means before the mixture is written;
/// calling_priors.md §5 carries the sizes.
pub struct SsrContamination<'a> { pub fraction: f64, pub length_distribution: &'a [f64] }
```

### 4.2 The stutter distribution — reused, and two cutoffs replace one

The distribution of spec §4.2 is **already built**:
[`alignment/stutter.rs:147`](../../../../src/ng/alignment/stutter.rs)'s `StutterModel`, constructed
from `StutterRates` ([`:82`](../../../../src/ng/alignment/stutter.rs)), evaluated by
`probability(bp_diff, period)` ([`:300`](../../../../src/ng/alignment/stutter.rs)), with the
one-step-share trap already typed against (its constructor doc). Three changes, all in that file
(ng code, not frozen production):

- **Rename the fields to spec §1.3's vocabulary** — `whole_repeat_longer_share`,
  `part_repeat_one_step_share`, … — keeping HipSTR's names in the doc comments for whoever ports
  alongside. *In frame / out of frame* is banned vocabulary in this doc set (spec §1.3), and the
  fields currently carry it (`in_up`, `out_geom`).
- **Two named cutoffs replace the single `MAX_SLIP = 10`**
  ([`alignment/stutter.rs:63`](../../../../src/ng/alignment/stutter.rs), whose own comment records
  the follow-up): `MAX_WHOLE_REPEAT_SLIP = 10` (repeats) and `MAX_PART_REPEAT_SLIP = 10` (base
  pairs) — both inherited from production's provisional 10
  ([`param_estimation.rs:21`](../../../../src/ssr/cohort/param_estimation.rs)) and declared
  inherited (spec §4.2).
- **The truncated mass is computed and reported per candidate** (feeds
  `SsrScoringContext.truncated_mass_lost`; spec §12 test 5).

Construction from the fit's numbers, with the placeholders named as placeholders:

```rust
/// Seven shares from the fit's three numbers (level, shorter share, fall-off) — the
/// part-repeat side is a placeholder pending an estimator with an owner (spec §10).
pub const PART_REPEAT_SHARE_OF_WHOLE: f64 = 0.05;  // placeholder; production's OUT_FRAME_REL
                                                   // (ssr/cohort/read_model/hipstr.rs:44)
// the two one-step shares are tied to one value — placeholder; HipSTR keeps them free
pub fn stutter_rates_for(slippage: &Slippage) -> StutterRates;
```

The numbers come from the pre-pass at the `(read group, stratum)` grain — `Slippage`
([`joint/ssr_fit.rs:83`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)) inside
`StratumFit` ([`:281`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)), the **level**
read off the fitted curve blend
([`joint/slippage_curve.rs:574`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs)
`blend_level`, provenance in `LevelSource`
[`:517`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs)). The lookup is by the
**candidate's** period and repeat count (spec §4.4).

### 4.3 The exact-reference shortcut and the ruler

Not this module's: measuring a read's tract is the alignment module's
([`../spec/alignment.md`](../spec/alignment.md) §4.2), the reference-identical skip lives in the
STR locus generator ([`ssr.rs:858`](../../../../src/ng/locus_generation/ssr.rs)), and the
substitution comparison this model composes is `FlatEmission`
([`alignment/emission.rs:250`](../../../../src/ng/alignment/emission.rs)) under the fitted rate —
per §4.3, never per-read qualities.

## 5. Reconciliation with existing code

Every row read on 2026-08-21.

| ng name | existing code | action |
|---|---|---|
| the shared copy-weighted mixture with outlier floor | `read_given_genotype`, [`likelihood.rs:75`](../../../../src/ssr/cohort/likelihood.rs) | shape ported for both paths; outlier spread made per-locus (spec §4.5) |
| SNP/indel closed form | `standard_log_likelihood`, [`per_group_merger.rs:1948`](../../../../src/var_calling/per_group_merger.rs) | **shape ported, two terms changed**: coefficient dropped (spec §3.4), `÷3` divisor added (spec §3.5) |
| contamination mixture | [`posterior_engine.rs:1475`](../../../../src/var_calling/posterior_engine.rs) (`compute_mixture_log_likelihoods`), c = 0 fallback [`:1509`](../../../../src/var_calling/posterior_engine.rs) | ported without the fallback branch — ng's two forms agree to ulp (spec §3.6) |
| minted per-read error (worse of BQ/MAPQ; min-BQ over window) | [`open_record.rs:792`](../../../../src/pileup/walker/open_record.rs), [`:944`](../../../../src/pileup/walker/open_record.rs) | upstream mint, consumed as `q_sum`; the calibration scale on top is **new** |
| `ReadGroupCalibration`'s fit input | [`generic/read_group_error_rate.rs:45`](../../../../src/ng/parameter_estimation/generic/read_group_error_rate.rs) `ReadGroupErrorRateFit` | consume; the two-scalar accumulator is asked of the pre-pass (spec §3.2) |
| `StutterModel` / `StutterRates` | [`alignment/stutter.rs:147`](../../../../src/ng/alignment/stutter.rs), [`:82`](../../../../src/ng/alignment/stutter.rs); production [`hipstr.rs:53`](../../../../src/ssr/cohort/read_model/hipstr.rs) | **reuse ng's**, with §4.2's three changes; do not port from the GPL HipSTR tree (spec §4.2's licence rule) |
| stutter parameters at the (read group, stratum) grain | `Slippage` [`joint/ssr_fit.rs:83`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs), `StratumFit` [`:281`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs), `blend_level` [`joint/slippage_curve.rs:574`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs) | consume, provenance included |
| placement enumeration (interrupted candidates) | [`ssr/cohort/stutter.rs`](../../../../src/ssr/cohort/stutter.rs) (whole-repeat only; part-repeat resized at tract end) | port with production's split, stated (spec §4.2) |
| substitution comparison | `FlatEmission`, [`alignment/emission.rs:250`](../../../../src/ng/alignment/emission.rs); production `pair_hmm.rs` | **compose**, never re-implement (spec §4.3, §7) |
| the model seam | `ReadLikelihoodModel` + `ReadScoringContext`, [`read_model/mod.rs:63`](../../../../src/ssr/cohort/read_model/mod.rs), [`:45`](../../../../src/ssr/cohort/read_model/mod.rs) | shape ported as `SsrEmissionModel` + `SsrScoringContext`, grown by the censored method and the truncation report |
| Model B comparator | [`classic.rs`](../../../../src/ssr/cohort/read_model/classic.rs) | port **test-only**, exactly as production keeps it (spec §9) |
| generic evidence | `CohortObservation` [`cohort_merge/build.rs:815`](../../../../src/ng/run/cohort_merge/build.rs), `SampleSupport` [`:858`](../../../../src/ng/run/cohort_merge/build.rs), `AlleleSupport` [`:973`](../../../../src/ng/run/cohort_merge/build.rs) | view over them; the `(allele × read group)` split **landed** in [`calling_prerequisites.md`](../impl_plan/calling_prerequisites.md) B1 — `SupportedAllele` carries `read_group` and the rows are one per pair, ascending |
| STR evidence | `SequenceObservation`, [`locus_generation/mod.rs:295`](../../../../src/ng/locus_generation/mod.rs) | reuse as-is; the read-group identity the contract needs is already there ([`:316`](../../../../src/ng/locus_generation/mod.rs)) |
| contamination inputs | `ContaminationEstimate` [`joint/contamination.rs:430`](../../../../src/ng/parameter_estimation/joint/contamination.rs), per-read-group grain [`:238`](../../../../src/ng/parameter_estimation/joint/contamination.rs) | consume as `ContaminationView`; the three allele-class frequencies are asked of the pre-pass's side-pass (spec §3.6) |
| numeric floors | `MIN_BASE_ERROR` [`contamination_estimation.rs:1449`](../../../../src/var_calling/contamination_estimation.rs); geometric clamps [`alignment/stutter.rs:74`](../../../../src/ng/alignment/stutter.rs) | import as named constants with reasons (spec §8) |
| censored term (both paths) | — | **new**; production discards partials ([`locus_tally.rs:91`](../../../../src/ssr/pileup/locus_tally.rs)) |

## 6. Design decisions — decided

- **The only swappable seam is the STR emission — decided.** The SNP/indel closed form has no
  competing impl, so it is a function, not a trait (module_layout 1a); the recipe's step-7 slot is
  `SsrEmissionModel` plus the tier configuration. Why: spec §2.4, §4.1.
- **Contamination on by default on both paths, no `c = 0` branch — decided.** The mixture *is* the
  plain formula at zero to a few ulp; the tolerance is a test constant (spec §12 test 11). Why:
  spec §3.6, §4.5.1.
- **The outlier spread is per-locus (`reachable_length_count`) — decided**, replacing production's
  cohort-wide denominator; goal 4's per-sample property is what it buys. Why: spec §4.5.
- **The stutter distribution is `alignment/stutter.rs`, renamed to spec vocabulary, with two
  cutoffs — decided.** One built implementation with two consumers, not two spellings that drift;
  the file already records both follow-ups its port deferred. Why: spec §4.2, §7.
- **The STR substitution term uses the fitted per-stratum rate, never `q_sum` — decided.** Unit
  mismatch: per-read versus per-base. Why: spec §4.3 (its Q6, closed).
- **Partials scored, never rescuing the merge's keep rule on the generic path — decided.** The
  one-line repeat-path amendment (test a partial against the reference over its witnessed run)
  belongs to the merge, recorded there. Why: spec §5.4.
- **`GenotypeTableView` flat views, not the owned `GenotypeTable` — decided**, for the same reason
  the prior takes flat slices: nothing allocates per sample per row, and the view borrows what
  `calling/genotype_table.rs` owns ([`calling_priors.md`](calling_priors.md) §7).

## 7. Open items

- `OPEN:` spec Q1 — whether the calibration scale needs a floor/ceiling on refused-ladder samples;
  the `ReadGroupCalibration` type stays a free ratio until that comparison runs.
- `OPEN:` spec Q7 — what would switch STR contamination off; the seam keeps both builds one binary.
- **Resolved 2026-08-21, and it is a blocker rather than a confirmation: the built merge discards
  partial observations.** Collation skips every observation whose witness is not `Complete`
  ([`cohort_merge/build.rs:1351`](../../../../src/ng/run/cohort_merge/build.rs)) and projection
  panics rather than pad one ([`:323`](../../../../src/ng/run/cohort_merge/build.rs)), so
  `SampleEvidence.partials` (§2.1) has no source today. **The merge must keep them, keyed and
  projected over the witnessed stretch** — now a named requirement in the spec's ownership table
  (spec §5.4 and §7, corrected 2026-08-21), and a dependency the censored term cannot be built
  without on either path.
- Impl-time confirmation: the concrete indexing of `contexts` per (read group, candidate).

## Test & bench shape

Unit tests beside each file, pinning spec §12's thirteen properties — the distribution-sums-to-one
tripwire (test 4, catches an inverted one-step share silently), the reported truncation loss (5),
bitwise aggregation (9), the same-function calibration check (10), the ulp-tolerance zero-
contamination identity (11), the censored complement against the truncated total, not 1 (12).
`ClassicEmissionOracle` is the independent implementation behind tests, never a shipping arm. The
change measurements (dropped coefficient; STR contamination on/off on the simulator; the
STR-vs-generic path at zero slippage) run through `bench/` on GIAB HG002 and the tomato interval —
spec §12 items 13–17; the two end-to-end regressions are shared with the siblings' benches so the
prior's effect and the likelihood's are attributable apart.
