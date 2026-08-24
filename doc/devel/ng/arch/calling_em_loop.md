# ng — the calling loop: types & interfaces

*Status: architecture draft (2026-08-21), companion to the spec
[`../spec/calling_em_loop.md`](../spec/calling_em_loop.md) (the design and its rationale) and to
the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) and
[`module_layout.md`](module_layout.md). One of three coordinated calling arch docs — the siblings
are [`calling_priors.md`](calling_priors.md) and [`read_likelihoods.md`](read_likelihoods.md); §0
says which doc owns each shared type. `ng_step_interfaces.md` §3 sketches step 9 as
`infer(&self, lik, prior, f)`; that sketch is superseded here (§3.3). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md); **STR** in prose ↔
`ssr` in code. Signatures are illustrative; the **contract** is the deliverable. The spec carries
every why.*

## 0. Who owns what — the three calling docs

| shared thing | owner |
|---|---|
| `CandidateAlleles` (the allele table), `GenotypeTable` + `AlleleId` + `GenotypeIdx` + `Genotype` (the candidate genotypes and their indexing), `ExpectedAlleleCopies`, `CallingScratch` (which holds the `Lg` **table**), `Phred`, `LocusInference`, and the two seam arguments `LocusEvidence` + `FrozenParameters` | **this doc** (§2) |
| the `Lg` **row** contract, the evidence views, the stutter distribution and scoring contexts | [`read_likelihoods.md`](read_likelihoods.md) §1–§2, §4 |
| `Concentration`, the seeds, the per-sample (leave-one-out) concentration | [`calling_priors.md`](calling_priors.md) §2–§3 |
| `LogProb`, `InbreedingF`, `Ploidy` | exist in [`src/ng/types.rs`](../../../../src/ng/types.rs) |

**The shared types live in `calling/mod.rs`, not in `inference/`** — one level up from all four
sub-modules, so each imports downward and no no-import rule is needed between them
([`module_layout.md`](module_layout.md) principle 1b). The prior still takes flat copy-count slices
and the likelihood flat `GenotypeTableView`s, for the no-allocation contract rather than for
dependency reasons.

## Module home

`src/ng/calling/inference/` — step 9, inside the `calling/` folder that holds steps 6–9
([`module_layout.md`](module_layout.md)). A folder of its own because the step has a real bake-off
(spec Q1's arms). **This doc also fixes what sits at the `calling/` level**, because the types it
owns are the ones the other three sub-modules consume:

```
src/ng/calling/
├── mod.rs                     – the vocabulary all four sub-modules share (§2):
│                                CandidateAlleles, ExpectedAlleleCopies, CallingScratch,
│                                LocusEvidence, FrozenParameters, LocusInference
├── genotype_table.rs          – GenotypeTable + GenotypeTableView (the port of production's
│                                GenotypeShape). Beside mod.rs, not inside inference/: the
│                                prior and the likelihood both consume its flat views
└── inference/                 – step 9, this doc's own module
    ├── mod.rs                 – LocusGenotyper seam, CallingLoopConfig
    ├── summarise_condition.rs – arm A: the three nested loops (frequency / slippage / discovery)
    ├── discovery.rs           – the stutter-hidden allele search (mechanism; config-gated)
    └── assignment.rs          – arms B/C: JointAssignmentPrior + the exhaustive/local scorers
```

Two homes, and the split is deliberate. The **calling-only** vocabulary of §2 lives in
`calling/mod.rs`. The **scalars other steps also name** seed `src/ng/types.rs`: `AlleleId(pub u16)`
(index into one locus's `CandidateAlleles`),
`Phred(f32)` (validated ≥ 0, conversions named — the `ng_step_interfaces.md` §1 sketch, now
real), and the opaque `Genotype(Box<[AlleleId]>)` multiset (output vocabulary; the loop's working
currency is a row index into the table).

## 1. General — the loop in one paragraph

Three loops, one inside the next, and the outer two ship switched off — so the default run is one
pass through their bodies and only the frequency loop repeats (spec §2). The expensive object is
the `Lg` table (`samples × genotypes`), built once per set of slippage numbers and reused by every
pass (spec goal 3); the per-pass work is `n` leave-one-out priors plus `n × genotypes` row scores.
The loop is linear in cohort size with no shared state; the M-step's sum runs in fixed sample
order, which is the whole of the determinism contract (spec §8).

## 2. Types

```rust
/// The alleles one locus is called over. REFERENCE is index 0, always present. Owned —
/// a discovery round appends to it and the final prune shrinks it (spec §4.1).
pub struct CandidateAlleles {
    pub alleles: Vec<Box<[u8]>>,   // sequence-resolved; AlleleId indexes this
    pub kind: LocusKind,           // routes the row builder and the prior's seed
}

/// Per-(ploidy, allele-count) genotype indexing, built once per shape and cached —
/// the port of production's GenotypeShape. The flat views are what the prior and the
/// likelihood consume; homozygous_allele_for is the ONE homozygous test
/// (calling_priors.md §3.2).
pub struct GenotypeTable {
    n_genotypes: usize,
    genotype_allele_counts: Vec<u32>,          // n_genotypes × n_alleles, row-major
    log_multinomial_coeffs: Vec<f64>,          // ln C(ploidy; counts) per genotype
    homozygous_allele_for: Vec<Option<AlleleId>>,
}
impl GenotypeTable {
    pub fn build(ploidy: Ploidy, allele_count: usize) -> Arc<Self>;  // cached per shape
    pub fn view(&self) -> GenotypeTableView<'_>;                     // the flat borrow
}
/// A row of one locus's GenotypeTable — the loop's working currency; the owned
/// `Genotype` multiset is minted from the row only at the final pass.
pub struct GenotypeIdx(pub u32);

/// The loop's feedback quantity: posterior-expected copies of each allele, summed over
/// the cohort in fixed sample order. Fractional, never a call (spec §1.3). Handed on in
/// the output — recomputing it from called genotypes gives a different number (spec §9).
pub struct ExpectedAlleleCopies(Vec<f64>);     // parallel to CandidateAlleles

/// Every buffer the per-locus work fills — allocated once per worker, reused per locus
/// (spec §8; the 16%-of-cycles allocator profile is the measured reason).
pub struct CallingScratch {
    pub lg_table: Vec<LogProb>,            // samples × genotypes — the rows the sibling fills
    pub posterior_row: Vec<f64>,           // genotypes
    pub concentration: Vec<f64>,           // alleles — handed to the prior's builder
    pub expected_copies: Vec<f64>,         // alleles, current
    pub expected_copies_prev: Vec<f64>,    // alleles — the convergence comparison
    pub per_sample_copies: Vec<f64>,       // samples × alleles
    pub row_scratch: RowScratch,           // the likelihood's own sections (its §1.1)
}

/// One locus's outcome. Evidence of HOW it was reached travels with it: nothing
/// downstream can otherwise tell a settled loop from a capped one (spec §6), a fitted
/// parameter from a defaulted one, or a seed that met its target from one that could not
/// (calling_priors.md §5).
pub struct LocusInference {
    pub region: GenomeRegion,
    pub alleles: CandidateAlleles,             // final — post-discovery, post-prune
    pub per_sample: Vec<SampleGenotypeCall>,   // run's sample order
    pub cohort_expected_copies: ExpectedAlleleCopies,
    pub converged: bool,                       // false = hit the pass cap; EMITTED, never dropped
    pub passes: u32,                           // instrument for Q4
    pub weakest_provenance: Provenance,
    pub seed_diversity_unreachable: bool,      // calling_priors.md §5's marker
}
pub struct SampleGenotypeCall { pub genotype: Genotype, pub genotype_quality: Phred }
```

**The two arguments the seam takes and this doc has not shaped yet, sketched here because a coder
cannot build against `call_locus` (§3.1) without them.** Both are assemblies of what other documents
own, so their *contents* are settled and only their wrapper is not:

```rust
/// One locus's evidence, per sample, in whichever shape the path uses. The enum is the
/// only place the two paths meet in this module — everything below it is path-pure, and
/// the row builder is chosen by the same discriminant that chose the candidates.
pub enum LocusEvidence<'a> {
    Generic { region: GenomeRegion, per_sample: &'a [GenericSampleEvidence<'a>] },
    Ssr     { region: GenomeRegion, per_sample: &'a [SsrSampleEvidence<'a>], detail: &'a SsrDetail },
}

/// Everything the pre-pass froze, borrowed for the run — the `ModelParams` successor of
/// `ng_step_interfaces.md` §2. Assembled once per run, never written during calling
/// (spec §5). Each field is owned by a parameter-prepass arch doc; this type only gathers
/// them so one borrow crosses the seam.
pub struct FrozenParameters<'a> {
    pub calibration: &'a [ReadGroupCalibration],     // read_likelihoods.md §2.3
    pub contamination: &'a [ContaminationView],      // read_likelihoods.md §2.3
    pub inbreeding: &'a [InbreedingF],               // per sample, run order
    pub seed: &'a SpectrumSeed,                      // calling_priors.md §2.3 — SNP/indel
    /// The (read group, stratum) slippage lookup — `joint/stratum_fits.rs`, built
    /// 2026-08-24 and pinned by `parameter_prepass_joint_fit.md` §1.7. `at` takes the
    /// read group and the **candidate's** period and repeat count, not the reference
    /// tract's (read_likelihoods.md §4.4), and returns the numbers the fit emitted with
    /// their provenance; it does not re-derive the level.
    pub ssr_strata: &'a StratumFits,
    pub ploidy: Ploidy,
}
```

**Contract on both:** `per_sample` is in the run's sample order and every per-sample slice in
`FrozenParameters` is indexed by that same order — one order for the whole run, which is what makes
the M-step's fixed-order sum meaningful (spec §8). `LocusEvidence`'s discriminant and
`CandidateAlleles::kind` must agree; disagreeing is a caller bug and asserts.

`GenotypeCall`'s quality here is the posterior-derived GQ of the final pass (production's
`final_calls`); step 13's `QualityModel` refines quality, it does not replace this field.

### 2.1 Configuration — every switch of spec Q1–Q3 is a value here, not a code path

```rust
pub struct CallingLoopConfig {
    /// On expected copies over chromosomes — the division is load-bearing across the
    /// cohort range (spec §6). Inherited, soft.
    pub convergence_threshold: f64,            // DEFAULT_CONVERGENCE_THRESHOLD = 1e-3 (unitless frequency scale)
    pub max_passes: u32,                       // DEFAULT_MAX_PASSES = 50 (production; observed need 3–5)
    pub max_candidate_alleles: u16,            // DEFAULT_MAX_CANDIDATE_ALLELES = 6 (production Stage-5 default)
    pub slippage_refit: SlippageRefitConfig,   // spec Q2 — §6.1
    pub discovery: DiscoveryConfig,            // spec Q3 — §6.2
}
```

## 3. The loop seam — spec Q1's two axes as two seams, one configuration each

### 3.1 The step-9 seam: how the cohort is handled

```rust
/// One locus, whole cohort → calls. The recipe selects the impl; everything the impls
/// share (evidence, parameters, candidates, scratch) crosses this one boundary.
pub trait LocusGenotyper {
    fn call_locus(
        &self,
        evidence: &LocusEvidence<'_>,          // per-sample views, per LocusKind (sibling §2)
        parameters: &FrozenParameters<'_>,     // pre-pass outputs, read-only
        candidates: CandidateAlleles,
        config: &CallingLoopConfig,
        scratch: &mut CallingScratch,
    ) -> LocusInference;
}

/// Arm A — the default: summarise the others, condition each sample on the summary,
/// leave-one-out subtraction in the prior. §2 of the spec, §4 and §6 here.
pub struct SummariseConditionLoop;

/// Arms B and C: score whole genotype-per-sample assignments under a joint prior.
/// The prior is the SECOND axis (§3.2); the enumeration is a config, so exhaustive
/// (the oracle) and search (the realistic cost) are the same scorer.
pub struct AssignmentGenotyper {
    pub joint_prior: Box<dyn JointAssignmentPrior>,
    pub enumeration: AssignmentEnumeration,
}
pub enum AssignmentEnumeration {
    /// Every assignment, marginalise exactly — small cohorts only; the piece that
    /// separates "the model disagrees" from "the search is too narrow" (spec Q1).
    Exhaustive { max_assignments: u64 },
    /// Start at the likelihood-preferred assignment, climb by bounded moves, sum the
    /// neighbourhood. Store MOVES, not assignment copies — freebayes' copy is what makes
    /// its path quadratic in cohort size, and that is its implementation, not the method
    /// (spec Q1's prediction).
    LocalSearch { max_rounds: u32 },
}
```

### 3.2 The second axis: which joint prior

```rust
/// log-prior of one whole-cohort assignment, from the per-allele chromosome counts it
/// implies plus the per-sample genotypes (the arrangement term needs the multiplicities).
pub trait JointAssignmentPrior {
    fn assignment_log_prior(&self, allele_chromosome_counts: &[u32],
                            genotypes_by_sample: &[GenotypeIdx],
                            table: &GenotypeTableView<'_>) -> LogProb;
}
/// Arm B: the same Dirichlet-multinomial model written jointly — the ported primitive
/// (genetics.rs:127) on the cohort's own counts; no leave-one-out, nothing conditioned.
pub struct CohortDirichletMultinomialPrior;
/// Arm C: freebayes' arrangement term × Ewens' sampling formula. Its total over all
/// assignments is (2+θ)/(1+θ), not 1 — REPORTED (a pinned test), never silently
/// normalised; a normalised variant would be a fourth arm (spec Q1).
pub struct EwensArrangementPrior;
```

An arm is one `(LocusGenotyper impl, JointAssignmentPrior impl, enumeration)` triple in the
recipe: A = `SummariseConditionLoop`; B = `AssignmentGenotyper` + `CohortDirichletMultinomialPrior`;
C = the same genotyper + `EwensArrangementPrior`. The fourth cell (summarise-and-condition ×
Ewens) has no closed form and is deliberately unrepresentable (spec Q1).

**Settled in the spec, 2026-08-21: every arm runs at `F = 0`.** The joint priors are written over
allele counts and have no per-sample slot; the exact composition is a mixture over which homozygous
samples are identical by descent, `2^k` terms at `k` homozygotes, so there is no cheap joint
counterpart to arm A's mixture (spec Q1). Running all three arms at `F = 0` is what makes the
comparison attributable. **So `JointAssignmentPrior` implementations take no inbreeding argument**,
and `LocusGenotyper` implementations for arms B and C reject a non-zero `InbreedingF` rather than
ignoring it — a silently dropped `F` on a selfing panel is the failure this is guarding.
`OPEN:` whether a whole-cohort scorer can carry inbreeding at all, which has to be answered before
either joint arm could replace arm A (spec Q1's follow-up).

### 3.3 Correction to `ng_step_interfaces.md` §3 step 9 — recorded, not applied there

`infer(&self, lik: &[GenotypeLikelihoods], prior: &dyn GenotypePrior, f: &[InbreedingF])` is
superseded. The loop is now **three nested levels** — discovery → slippage re-fit → frequencies —
with the outer two off by default (why: spec §2), so a likelihood table cannot arrive prebuilt as
an argument: the loop itself rebuilds it per slippage round and appends columns per discovery
round. The genotyper therefore takes evidence and frozen parameters and *drives* the two sibling
functions, rather than receiving their outputs; and the prior argument is gone from the signature
because each sample's prior is built per pass from the loop's own expected copies
([`calling_priors.md`](calling_priors.md) §3.1) — there is no locus-constant prior object to pass.

## 4. General — the passes, starting, stopping

- **First pass of every outer round is prior-free** (reads only) — it exists to mint the expected
  copies the prior needs; later rounds restart it because their parameters changed (spec §3).
  Production's `EmStepPhase` is the ported shape. At one sample the flat pass is wasted, not
  wrong; skipping it is spec Q6, and no cohort-size branch is written until Q6 says so (spec §7).
- **Stopping**: largest change in expected copies **divided by cohort chromosomes**, against
  `convergence_threshold`; the division is what keeps the constant meaningful from 1 to several
  thousand samples (spec §6). The slippage rounds stop on their own numbers' movement or the round
  cap — production's rule, not HipSTR's likelihood test, because a likelihood test needs the table
  rebuilt to be read (spec §6).
- **A capped locus is emitted with `converged = false`** — never dropped, never an error; `passes`
  is always emitted (Q4's instrument). Errors proper are caller bugs → assertions (spec §6, §8).

## 5. The SNP/indel path

Exactly one quantity moves: the allele frequencies. The two outer rounds are structurally inert —
the slippage round runs only for `LocusKind::Ssr` loci, and discovery's retrace is defined on
stutter attribution, so both configs are ignored on `Generic` loci rather than half-honoured (spec
§5.1's closing paragraph gives the reasons: no slippage numbers; error rate frozen because a refit
would measure the merge's selection; contamination frozen by grain). Candidate selection is flat —
a cap and a support bar over the merge's already-unified table (spec §4); the selection step's own
design has no spec yet, and this doc only fixes its output type (`CandidateAlleles`).

## 6. The STR path

### 6.1 The slippage rounds — spec Q2's three pull-back strengths, one code path

```rust
pub struct SlippageRefitConfig {
    /// 0 = frozen — THE DEFAULT, and the whole frozen arm: no rounds, no rebuild.
    pub max_rounds: u32,                          // DEFAULT_SLIPPAGE_REFIT_ROUNDS = 0; production caps at 3
    /// Pseudo-count weight pulling the re-fit direction split + fall-off back toward the
    /// per-stratum values. Production's 50; 0.0 is the free (HipSTR-setting) arm.
    pub shape_pull_back_pseudocounts: f64,        // DEFAULT_SHAPE_PULL_BACK = 50.0 (inherited, unmeasured)
    /// Slipped-read weight pulling the level back. Production's 20; 0.0 = free.
    pub level_pull_back_slipped_reads: f64,       // DEFAULT_LEVEL_PULL_BACK = 20.0 (inherited, unmeasured)
    pub round_convergence_threshold: f64,         // largest re-fitted-number change; 1e-3, inherited
}
```

Contract (spec §5.1, decided there): the nested shape — an outer round re-fits, rebuilds the `Lg`
table, and reruns the frequency loop — so frozen / pulled-back / free are three values of this
config, not three implementations. The re-fit reads the **genotype posteriors**, never called
genotypes (HipSTR's choice, production's deferred refinement); it moves **three numbers**, holding
the part-repeat placeholders fixed. The level's pull-back target is **the fitted curve's value at
the cell**, not a cell's own estimate — the curve redefinition postdates the spec's Q2 wording and
is what "far from its stratum" is now measured against (spec §1, Q2). The emission-call count per
locus is `candidates × Σ_s observations × builds` with `builds ≤ rounds + 1` — the property spec
§13 test 5 pins.

### 6.2 Discovery — spec Q3's three settings, one mechanism

```rust
pub enum DiscoveryMode {
    Off,                        // THE DEFAULT
    /// Converge once, freeze the frequencies, discover against them — each round is a
    /// scoring pass; converge once more at the end (the middle setting, spec §4.1).
    AgainstFrozenFrequencies,
    /// Full convergence every round (the outermost repeat of spec §2's pseudocode).
    AgainstFullConvergence,
}
pub struct DiscoveryConfig { pub mode: DiscoveryMode, pub bar: DiscoveryBar, pub max_rounds: u32 }
/// Both halves must clear — a single stray read cannot mint an allele. Values inherited
/// from HipSTR's high-depth setting and SOFT: below ~13 reads only the count binds
/// (spec §4.1, Q3 sweeps them).
pub struct DiscoveryBar { pub min_reads: u32 /* 2 */, pub min_spanning_read_share: f64 /* 0.15 */ }
```

Contract (spec §4, §4.1): discovery runs only **between whole runs of the frequency loop** — the
prior and the STR outlier spread do not survive a candidate being added, the emission columns do
(spec §4's table) — and a round **appends** the new alleles' emission columns, never rebuilds
(spec §13 test 12). After the last round, alleles no sample's best genotype used are pruned and
the frequency loop reruns on what is left. With `Off`, emission counts and genotypes are
bit-for-bit those of the plain loop (spec §13 test 11). Extent never grows — alleles within the
locus, only (spec §1.2).

Candidate selection on this path is rung-structured (production's ladder), and the rung/prior
coupling is the sibling's ([`calling_priors.md`](calling_priors.md) §5); this doc only requires
that discovered sequences enter `CandidateAlleles` through the same admission as selected ones.

## 7. Where it runs

Inside the merge's builder, on the region the builder owns — decided in the spec for the memory
reason; it commutes either way (spec §9). The organiser's `CohortObservation` stream
([`cohort_merge/organise.rs:271`](../../../../src/ng/run/cohort_merge/organise.rs)) is the input
edge; `LocusInference` is what buffers instead.

## 8. Reconciliation with existing code

Every row read on 2026-08-21.

| ng name | existing code | action |
|---|---|---|
| the loop internals (E-step / M-step / delta behind one seam; flat first pass) | `GenotypeEmModel` [`posterior_engine.rs:2635`](../../../../src/var_calling/posterior_engine.rs), `run_em_loop` [`:2733`](../../../../src/var_calling/posterior_engine.rs), `EmStepPhase` [`:2586`](../../../../src/var_calling/posterior_engine.rs) | **shape ported** into `SummariseConditionLoop`; ng needs no inner trait — the paths differ only in the sibling row builders |
| convergence on expected copies over chromosomes | [`posterior_engine.rs:2704–2726`](../../../../src/var_calling/posterior_engine.rs) (the comment carries both halves of the reason) | port arithmetic + reasoning |
| `DEFAULT_CONVERGENCE_THRESHOLD`, `DEFAULT_MAX_PASSES` | [`posterior_engine.rs:86`](../../../../src/var_calling/posterior_engine.rs), [`:96`](../../../../src/var_calling/posterior_engine.rs) | inherit, marked soft (spec Q4) |
| emit-with-flag on non-convergence | [`posterior_engine.rs:26`](../../../../src/var_calling/posterior_engine.rs) (module doc: `DidNotConverge` retired for `converged = false` → `FILTER=EMNoConv`) | port rule + reasoning |
| `GenotypeTable` | `GenotypeShape` + `shape_for` cache, [`posterior_engine/shape.rs:42`](../../../../src/var_calling/posterior_engine/shape.rs), [`:76`](../../../../src/var_calling/posterior_engine/shape.rs) | port; `nonzero_pairs` comes along if profiling asks |
| `DEFAULT_MAX_CANDIDATE_ALLELES = 6` | Stage-5 `max_alleles` default, recorded at [`shape.rs:35`](../../../../src/var_calling/posterior_engine/shape.rs) | inherit, soft |
| `CallingScratch` (lift out of the iteration) | [`posterior_engine.rs:1874`](../../../../src/var_calling/posterior_engine.rs) (allocator ~16% of cycles before the lift) | port with its measured reason |
| the nested slippage rounds | [`em.rs:571–604`](../../../../src/ssr/cohort/em.rs) (re-fit → rebuild → rerun; break-before-rebuild order), config + defaults [`em.rs:60–75`](../../../../src/ssr/cohort/em.rs), [`:136–145`](../../../../src/ssr/cohort/em.rs) (`refit_max_rounds 3`, `theta_shrink 50`, `level_shrink 20`, `lambda 0.01`) | **shape ported, default not** — rounds ship at 0 (spec §5.1) |
| the pull-back | `refine_theta_locus`, [`ssr/cohort/stutter.rs:181`](../../../../src/ssr/cohort/stutter.rs) | port; input changes to genotype posteriors |
| read attribution | `attribute_locus`, [`em.rs:1189`](../../../../src/ssr/cohort/em.rs) (called-genotype input; the soft split its comment defers) | **shape ported, input not** — posteriors (spec §5.1) |
| the frequency M-step + final calls | `run_pi_em` [`em.rs:816`](../../../../src/ssr/cohort/em.rs), `final_calls` [`:857`](../../../../src/ssr/cohort/em.rs) | one ng loop for both paths |
| STR loop differential (not oracle) | π-convergence at `1e-6`, [`em.rs:137`](../../../../src/ssr/cohort/em.rs) | run ng under production's rule, require matching genotypes, then restore and report (spec §10) |
| `CandidateAlleles` (generic input) | `CohortObservation.alleles`, [`cohort_merge/build.rs:815`](../../../../src/ng/run/cohort_merge/build.rs) | select from it (cap + bar); the selection step has no spec — this doc fixes only the output type |
| `CandidateAlleles` (STR input) | `assemble_candidates` + `occupied`, [`candidate_set.rs:193`](../../../../src/ssr/cohort/candidate_set.rs), [`:221`](../../../../src/ssr/cohort/candidate_set.rs) | port the rung ladder when the STR path comes through the merge |
| arm B's prior | [`genetics.rs:127`](../../../../src/genetics.rs) on cohort counts | reuse the sibling's port — one primitive, two callers |
| arm C's prior | [`freebayes/src/Ewens.cpp`](../../../../freebayes/src/Ewens.cpp) (51 lines) | reimplement natively; store moves, not assignment copies (spec Q1) |
| slippage-curve pull-back target | `SlippageCurve` / `blend_level`, [`joint/slippage_curve.rs:141`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs), [`:574`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs) | consume as the level's shrink target (spec Q2 preamble) |
| discovery mechanism | HipSTR's, described in full in spec §4.1 | build **from the spec's description**, not the GPL tree ([`read_likelihoods.md`](../spec/read_likelihoods.md) §4.2's licence rule) |

## 9. Design decisions — decided

- **Arm A is the default; arms B/C and the exhaustive oracle share one scorer — decided.** The two
  axes are two seams (`LocusGenotyper`, `JointAssignmentPrior`); exhaustive-vs-search is an
  enumeration config. Why: spec Q1 (all arms built and compared; B−A isolates the shape, C−B the
  prior).
- **The nested shape for slippage re-fit; frozen is rounds = 0 — decided.** One code path, three
  configurations; the flat shape would rebuild the table per pass. Why: spec §5.1.
- **Re-fit reads posteriors, moves three numbers, shrinks toward the curve — decided.** Why: spec
  §5.1, Q2 preamble.
- **Discovery against converged posteriors, off by default, three modes in one enum — decided.**
  Why: spec §4, §4.1, Q3.
- **No cohort-size branch anywhere in the loop — decided.** One sample terminates by arithmetic
  (leave-one-out is exactly zero); the wasted pass is Q6's measurement, not a rule. Why: spec §7.
- **`ExpectedAlleleCopies` is part of the output — decided.** Site filtering reads it; recomputing
  from calls discards the carried uncertainty. Why: spec §9.
- **Ewens' non-normalisation is pinned, not repaired — decided.** A normalised variant is a fourth
  arm. Why: spec Q1, §13 test 10.

## 10. Open items

- `OPEN:` whether a whole-cohort scorer can carry per-sample inbreeding at all (§3.2). The
  comparison itself is settled — all arms at `F = 0` — but a joint arm cannot replace arm A until
  this is answered, because a caller that cannot express `F` is unusable on a selfing crop.
- `OPEN:` spec Q4–Q7 — thresholds/cap from the pass-count distribution; loop-vs-likelihood cost
  crossover; the one-sample skip; the flat first pass's measured effect. All instrumented by
  `passes` and the emission-call counter; none blocks building.
- Impl-time confirmations: both `FrozenParameters` and `LocusEvidence` are sketched in §2 with their
  ordering contract; what stays open is only the concrete type of each borrowed field, which the
  parameter-prepass arch docs pin. **`StratumFits` was built on 2026-08-24** and is pinned by
  [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §1.7, so this item is closed;
  what its own §1.7 leaves open is a candidate repeat count no kept reference tract occupies, and
  the STR substitution rate, which sits at the same grain and is not in the gather.

## Test & bench shape

Unit tests beside each file, pinning spec §13's twelve properties — the one-sample fixed point
asserted on pass-1-equals-pass-2, bit for bit (test 1); the bitwise M-step mutation check on
summed copies, not on the argmax (test 2); the flat-first-pass trap (test 3); the instrumented
emission-call count `candidates × Σ_s observations × builds`, run at both re-fit settings (test 5);
no-op re-fit collapse (6); zero allocations per pass (7); the production parity oracle and the
STR differential (8, §10); exhaustive-vs-search agreement on enumerable cohorts (9); the Ewens
total (10); discovery's plant/terminate/free-when-off triple (11) and append-only columns (12).
The bench arms live in `bench/`: Q1's three arms × the cost/disagreement report at 1, 63 and 1,000
samples; Q2's three pull-backs and Q3's three modes on the HG002 tandem-repeat bundle and the
tomato silver standard. The end-to-end anchors are the siblings' two regressions — this doc owns
neither number, only that loop changes move them for recorded reasons (spec §13).
