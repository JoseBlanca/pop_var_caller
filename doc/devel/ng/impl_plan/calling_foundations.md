# ng calling foundations — implementation plan

**Status:** draft, 2026-08-21. The build order for the **first `src/ng/calling/` code**: the
scalars calling adds to `types.rs` (`AlleleId`, `Phred`, `Genotype`), the shared vocabulary in
`calling/mod.rs`, and `calling/genotype_table.rs` — the port of production's `GenotypeShape`.
Design is settled in the three calling arch docs —
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §2 owns every type this plan builds —
under [`../arch/module_layout.md`](../arch/module_layout.md) (the one `calling/` folder for steps
6–9) and [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) §1 (newtype
conventions). This plan turns that design into build order; it is **not** a place for new design.

**Small on purpose: it unblocks the fan-out.** Six plans build calling:
[`calling_prerequisites`](calling_prerequisites.md) ∥ `calling_foundations` →
[`calling_prior`](calling_prior.md) ∥ [`calling_read_likelihoods`](calling_read_likelihoods.md) →
[`calling_loop`](calling_loop.md) → [`calling_bakeoffs`](calling_bakeoffs.md). This plan runs **in
parallel with the prerequisites plan** and neither needs the other. Both fan-out plans consume
what it builds — the prior and the likelihood both take the genotype table's flat views — so the
day it merges, the prior plan can branch; the read-likelihoods plan additionally waits on
prerequisites items 1–5.

---

## Scope

**In:** `src/ng/calling/` scaffold; `types.rs` gains `AlleleId`, `Phred`, `Genotype`;
`calling/mod.rs` gains the vocabulary that needs nothing from the sibling plans —
`CandidateAlleles`, `ExpectedAlleleCopies`, `LocusInference` + `SampleGenotypeCall`;
`calling/genotype_table.rs` — `GenotypeTable`, `GenotypeTableView`, `GenotypeIdx`, the per-shape
cache.

**Out (later plans):**

- **`LocusEvidence`, `FrozenParameters`, `CallingScratch`** — they borrow the sibling modules'
  types (`GenericSampleEvidence`, `ReadGroupCalibration`, `SpectrumSeed`, the row scratch), which
  do not exist until the fan-out plans build them; [`calling_loop.md`](calling_loop.md) owns all
  three ([`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §2).
- **All four step sub-modules** — `genotype_prior/` and `likelihood/` are the fan-out plans;
  `inference/` is the loop plan; `allele_candidates/` has **no spec** (both calling specs record
  the gap) and no plan exists for it.
- **`InbreedingF`'s `[0, 1)` tightening** — [`calling_prerequisites.md`](calling_prerequisites.md)
  Milestone A, in parallel.

## Principles (how the order was chosen)

- **Types first, then implementation** (project rule) — the scalars before the table that indexes
  with them.
- **Verify against ground truth.** The genotype table is a port, and production's
  `GenotypeShape` + `shape_for` cache
  ([`posterior_engine/shape.rs:42`](../../../../src/var_calling/posterior_engine/shape.rs),
  [`:76`](../../../../src/var_calling/posterior_engine/shape.rs)) is callable in-crate: the
  north-star test is **value parity with production's shape** across `(ploidy, allele count)`
  grids, not self-consistency.
- **Foundations set the conventions.** `AlleleId` and `Phred` follow the newtype rules the tree
  already uses (unconstrained → `pub` field; constrained → private field + `try_new`); the region
  discipline against the parallel prerequisites branch is stated below and kept.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place)

- `src/ng/types.rs` exists with `LogProb`
  ([`types.rs:250`](../../../../src/ng/types.rs)), `Ploidy`
  ([`:419`](../../../../src/ng/types.rs)), `GenomeRegion` ([`:79`](../../../../src/ng/types.rs));
  `AlleleId`, `Phred`, `Genotype` do **not** — this plan adds them.
- `LocusKind` exists
  ([`locus_generation/mod.rs:422`](../../../../src/ng/locus_generation/mod.rs)) — reused as
  `CandidateAlleles.kind`, not re-minted.
- `Provenance` exists
  ([`parameter_estimation/mod.rs:60`](../../../../src/ng/parameter_estimation/mod.rs)) — carried
  on `LocusInference`.
- Production's `GenotypeShape` — the reuse target **and** the parity oracle.

## Worktree, branch, merge

- **Worktree** `../pop_var_caller-calling-foundations`, **branch** `ng-calling-foundations`, from
  `main`, plain `git worktree add`.
- **Runs in parallel with** `ng-calling-prerequisites`. Shared file: `src/ng/types.rs`. Region
  discipline that avoids the conflict rather than resolving it: this branch **appends** its three
  scalars at the end of the sections they belong to and appends its `DomainError` variant
  (`Phred`'s) **at the end of the enum**; prerequisites edits only the existing `InbreedingF`
  block and inserts its variant beside the existing `InbreedingF` variant. Disjoint regions.
- **Merge order back:** whichever branch finishes first merges first; the second merges `main` in
  and re-runs. `src/ng/calling/` is created only here, so no other phase-1 branch can conflict on
  it.
- **What must merge before the fan-out:** the prior plan branches from `main` once **this** branch
  has merged; the read-likelihoods plan branches once **both** phase-1 branches have merged.

---

## The steps

### Milestone A — the `types.rs` scalars

**A1. `AlleleId` and `Phred`.**  ☐
`AlleleId(pub u16)` — index into one locus's candidate-allele table, unconstrained newtype with
the ergonomic derives. `Phred(f32)` — constrained: validated `≥ 0` and finite via `try_new` + a
new `DomainError` variant, conversions as **named functions** (`Phred::from_log_prob`,
never `as`), per the interfaces doc's sketch, now real. Unit tests: `Phred` boundary both
directions; `from_log_prob` on a hand-computed pair. *Source:* calling_em_loop arch §Module home;
ng_step_interfaces §1.

**A2. `Genotype`.**  ☐
`Genotype(Box<[AlleleId]>)` — the opaque output multiset, alleles stored sorted so equal
genotypes compare equal; a constructor that sorts, `.alleles()` accessor. The loop's *working*
currency is `GenotypeIdx` (Milestone C) — this type is minted only at the final pass, which is why
it is small and owns no arithmetic. Test: construction order does not change equality. *Depends:*
A1. *Source:* calling_em_loop arch §2, §Module home.

> **Checkpoint A:** the three scalars compile with tests; `types.rs` conventions kept (constrained
> vs unconstrained). Pause for review.

### Milestone B — the `calling/` scaffold and its vocabulary

**B1. Scaffold.**  ☐
`src/ng/calling/mod.rs` (declares `genotype_table`; the step sub-modules arrive with their plans)
wired into `ng/mod.rs`. One folder for steps 6–9, per module_layout's dependency argument (keeping
them apart forced a no-import rule). *Source:* module_layout §The tree; calling_em_loop arch
§Module home.

**B2. `CandidateAlleles` and `ExpectedAlleleCopies`.**  ☐
`CandidateAlleles { alleles: Vec<Box<[u8]>>, kind: LocusKind }` — REFERENCE at index 0, always
present; owned, because a discovery round appends and the final prune shrinks (a later plan's
behaviour, this plan's shape). `ExpectedAlleleCopies(Vec<f64>)` — parallel to the allele table;
fractional, never a call. Doc comments carry both contracts. Tests: reference-at-zero invariant.
*Depends:* B1, A1. *Source:* calling_em_loop arch §2; spec §1.3.

**B3. `LocusInference` and `SampleGenotypeCall`.**  ☐
The outcome type as the arch writes it: region, final alleles, per-sample calls in run order,
`cohort_expected_copies`, `converged` (false = capped, **emitted, never dropped**), `passes`,
`weakest_provenance`, `seed_diversity_unreachable`. Plain data, no logic; the loop plan fills it.
*Depends:* B2, A2. *Source:* calling_em_loop arch §2; spec §6, §9.

> **Checkpoint B:** the shared vocabulary compiles and is importable downward by the future
> sub-modules. Pause for review.

### Milestone C — `genotype_table.rs`, the `GenotypeShape` port

**C1. `GenotypeTable` + `GenotypeIdx` + the flat views.**  ☐
The port: per-`(ploidy, allele count)` genotype indexing built once and cached
(`build(ploidy, allele_count) -> Arc<Self>`), holding `genotype_allele_counts`
(`n_genotypes × n_alleles`, row-major), `log_multinomial_coeffs` (`ln C(ploidy; counts)` per
genotype), and `homozygous_allele_for` (`Vec<Option<AlleleId>>` — **the one homozygous test** the
prior consumes; nothing else may decide homozygosity, so the above-diploidy spec has one place to
change). `view()` returns `GenotypeTableView<'_>` — the flat borrow both siblings take.
`nonzero_pairs` comes along only if profiling asks. *Depends:* A1, B1. *Source:* calling_em_loop
arch §2; calling_priors arch §3.2, spec §3.3.

**C2. Parity with production.**  ☐
Across ploidy 2 and 4 and allele counts 1–6: genotype count, every row's allele counts in the
same enumeration order as production, every `ln`-coefficient to floating-point equality, and the
homozygous lookup, all against `GenotypeShape`
([`shape.rs:42`](../../../../src/var_calling/posterior_engine/shape.rs)) called directly. Plus a
cache test: two `build` calls for one shape return one `Arc`. **Own commit, do not bundle** — a
wrong coefficient or a swapped enumeration order is a quietly-wrong prior and posterior for every
locus, not a crash; the oracle is production's shape, green before and after. *Depends:* C1.
*Source:* calling_em_loop arch §8 (reconciliation row).

> **Checkpoint C:** the table matches production value-for-value; the flat views are consumable.
> **The fan-out is unblocked** — the prior plan can branch once this merges. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | newtype boundary tests; named-conversion test on a hand-computed pair |
| B | invariant tests (reference at index 0); compiles as the import root for future sub-modules |
| C | **value parity with production's `GenotypeShape`** over a (ploidy × allele-count) grid — counts, coefficients, homozygous lookup, enumeration order — plus the cache identity |

## Out of scope (next plans)

- **`GenotypePriorModel` and everything in `genotype_prior/`** —
  [`calling_prior.md`](calling_prior.md).
- **The evidence views, the `Lg` row, `likelihood/`** —
  [`calling_read_likelihoods.md`](calling_read_likelihoods.md).
- **`CallingScratch`, `LocusEvidence`, `FrozenParameters`, `inference/`, `CallingLoopConfig`** —
  [`calling_loop.md`](calling_loop.md).
- **`allele_candidates/` (step 6)** — no spec exists; both calling specs record the gap, and
  [`calling_loop.md`](calling_loop.md) carries the consequence (fixture-supplied candidates).
- **The `types.rs` split into concept modules** — module_layout principle 3 says split when the
  file tells you to; adding three scalars is not that moment.
