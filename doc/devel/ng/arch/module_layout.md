# ng — module layout

*Status: architecture draft (2026-07-10), companion to
[`ng_step_interfaces.md`](ng_step_interfaces.md). That doc defines the shared **types
and step traits**; this one is their **physical home** — how the ng caller is laid out
as Rust modules, and the rules that keep the layout serving the lab (§3 of the
[spec](../spec/ng_proposal.md)) rather than fighting it.*

## Where ng lives

ng is a **new module inside `pop_var_caller`**, not a separate crate: `src/ng/`. That
keeps it able to **reuse the existing code directly** (the filters, `ReadLikelihoodModel`,
`GenotypeEmModel`, …; see the reconciliation table in `ng_step_interfaces.md` §6) and
means the eventual port-back of winning steps into the production engine is an
in-crate move, not a dependency dance.

## The tree

```
src/ng/
├── mod.rs            – module declarations + re-exports only (kept minimal)
├── types.rs          – the shared vocabulary, one file to start: the domain newtypes
│                       (Bp, SsrPeriod, LogProb, ids, ReadWeight…) plus the wire structs
│                       (GenomeRegion, LocusKind, Genotype, ModelParams, …). A deliberate
│                       *temporary* catch-all: splits into concept modules (units, locus,
│                       genotype, params) as clusters grow — see principle 3.
├── ref_seq.rs        – reference-sequence access: the RefSeq trait + impls (resident,
│                       streaming, in-memory). Foundational infra (NOT a step), shared by
│                       read filtering (#8), pileup/, BAQ, DUST. Reuses src/fasta. Splits
│                       into ref_seq/ when the impls grow. Spec: ../spec/ref_seq.md.
│
│   # one module per pipeline step — each owns its trait + its swappable impls + tests
├── read/             – steps 1+2, one read-handling module (see principle 1, note):
│                        · aligned_read.rs – AlignedRead: ng's decoded read, and the decode
│                          that builds it. Beside its producer rather than in types.rs,
│                          the way production keeps MappedRead beside its own reader and
│                          ng keeps SampleLocusObservations in locus_generation/: it is
│                          filtering's output, not free-floating vocabulary. Carries the
│                          read group (../spec/read_groups.md §8)
│                        · filtering.rs – step 1: the fixed filtering prelude (a single
│                          file — no bake-off; wraps the bam/alignment_input filters)
│                        · mod.rs + left_align_baq.rs – step 2 (GENERIC PATH ONLY): the
│                          ReadPreparer trait + its impls, side by side. Read preparation is
│                          a generic-only step (../spec/read_preparation.md §1); the STR path
│                          has no preparer — its per-read alignment produces an observation and
│                          lives in locus_generation/ssr.rs, calling alignment/. The alignment
│                          algorithms a preparer calls live in alignment/, below — it composes
│                          one, it does not contain one.
├── alignment/        – NOT a pipeline step: the shared alignment algorithms both step 2 and
│                       step 7 call. Knows nothing about caller steps. mod.rs (the traits) +
│                       one file per algorithm. Spec ../spec/alignment.md
├── region_typing/    – step 3  the typed-region generator: walks the reference and cuts it into
│                       TypedRegion { region, kind } (SsrSegment/SsrBundle/Generic/Satellite). A
│                       concrete iterator, no trait (no bake-off). Spec ../spec/typed_regions.md;
│                       arch typed_regions.md
├── locus_generation/ – the locus generators: turn one typed region into a sample's loci.
│                       LocusGenerator<S> trait + one impl per segment kind, side by side —
│                       a step folder, because alternatives per kind are expected.
│                        · mod.rs   – SampleLocusObservations, SequenceObservation, LocusGenerator,
│                          the dispatch match, NoLoci (the count-only generator)
│                        · ssr.rs   – the STR generator (../spec/locus_generation_ssr.md)
│                        · pileup/  – the generic generator: walks a non-STR stretch, splits
│                          it into loci from the data. Reuse target: production pileup/walker/.
│                       Spec ../spec/locus_generation.md. See *The locus stream*.
├── psp/              – the per-sample store: what one sample's reads showed at every position
│                       a run analysed, written once by the locus generation stage and read back
│                       by the cohort gather. NOT a pipeline step — it is the seam between the
│                       two stages, made durable. Eleven files: the writer, the reader and its
│                       walks, and one per part of the file (header, block, block index, trailer,
│                       footer, chain ids, record). Spec ../spec/psp_file_format.md; the store is
│                       built and measured, and nothing in pipeline.rs writes or reads one yet —
│                       see *The artifact between the two stages*.
├── parameter_estimation/ – step 4  the parameters the caller runs on — noise rates and
│                       rates of variation — estimated from the data before calling, and
│                       emitted as ModelParams. Named for what it owns, not for when it
│                       runs: "pre-pass" says only that something else comes after, and
│                       "prior" would both collide with calling/genotype_prior/ below and overstate
│                       the output, half of which is noise-model terms rather than priors.
│                       SampleSummarizer + CohortEstimator. Specs ../spec/parameter_prepass*.md;
│                       arch parameter_prepass_generic.md (which owns the shared fitting/
│                       machinery) and parameter_prepass_ssr.md. Sub-units: fitting/,
│                       generic/, ssr/, joint/ — STR-ness is a property of an implementation
│                       (principle 2), not a top-level split.
│                        · joint/ – step 4's second route to the same parameters: instead of
│                          folding each sample into histograms, it keeps raw evidence at a
│                          bounded set of positions that is the same in every sample — the
│                          census — and fits the whole cohort at once, which is what buys a
│                          locus's own allele frequency. loci.rs (which positions), census.rs
│                          (what is recorded there), fit.rs and ssr_fit.rs (the two halves of
│                          the estimator), contamination.rs. Arch parameter_prepass_joint_*.md.
├── calling/          – steps 6–9, one folder because they are one question: which alleles,
│                       how probable are the reads under each genotype, how likely is each
│                       genotype before the reads, and what comes out. Principle 1's rule (b)
│                       — the loop *drives* the other three and they share the allele table,
│                       the genotype indexing and the per-locus scratch, so siblings at the
│                       top of ng/ would have forced a no-import rule between them. Steps 10
│                       and 11 are NOT here: phasing and locus filtering read calling's
│                       output rather than taking part in it. Specs ../spec/calling_*.md and
│                       ../spec/read_likelihoods.md; arch calling_priors.md,
│                       read_likelihoods.md, calling_em_loop.md.
│                        · mod.rs           – the vocabulary the four share: CandidateAlleles,
│                          GenotypeTable + AlleleId + GenotypeIdx + Genotype, ExpectedAlleleCopies,
│                          CallingScratch, LocusEvidence, FrozenParameters, LocusInference
│                        · allele_candidates/ – step 6: selection from the merge's table
│                          (flat cap + bar [generic], rung ladder [STR]). No spec yet — both
│                          calling specs record that gap
│                        · likelihood/      – step 7: the Lg row. One seam only, the STR
│                          emission (SsrEmissionModel); the SNP/indel row is a function
│                        · genotype_prior/  – step 8: GenotypePriorModel + the marginalized
│                          default and the plug-in comparator
│                        · inference/       – step 9: LocusGenotyper (summarise-and-condition
│                          [default] / whole-cohort assignment scoring) × JointAssignmentPrior
├── phasing/          – step 10 (physical / SNP-based)
├── locus_filter/     – step 11  the locus-level filters (add more as needed):
│                        · hidden_dup.rs  – ArtifactFilter (11a: paralog / hidden-duplication)
│                        · emission/      – EmissionModel (11b: heuristic, bic, freebayes)
├── allele_representation/ – step 12 AlleleRepresentation
├── quality/          – step 13 QualityModel
│
├── pipeline.rs       – the CallerRecipe + the driver that runs it end-to-end (per-sample
│                       stage, then the cohort gather; the artifact between them is still passed
│                       in memory — psp/ exists but is not wired in)
└── bench/            – the standards harness: gold / silver / synthetic scoring
```

(Steps 5 — STR read-class / spanning — and the read-class machinery live inside
`calling/allele_candidates/` or `locus_generation/ssr.rs`; see *Open items*.)

## Organizing principles

**1. One module per pipeline step — trait, impls, and tests together.** This is the
load-bearing rule. Each step folder holds its trait *and* every competing implementation
*and* their tests, so a step's alternatives sit **side by side**. That is exactly what
the bake-off needs — "swap `calling::allele_candidates::assembly` for `::rung_ladder`,
hold the rest, re-measure." It also satisfies the naming rule ([naming.md](../../../../ai/skills/rust-code-review/code_review/naming.md)):
modules are named for the **concept** they own (`calling`, `allele_candidates`, `genotype_prior`),
never for a layer or pattern (`models`, `services`, `common`, `utils`). *One step may host
sub-modules when the spec frames it as one step with several questions* — `locus_filter/`
(step 11) holds `hidden_dup` (artifact) and `emission` (emit decision) side by side, each
with its own trait and impls.

*Two rules of thumb the read-filtering spec pinned down (`../spec/read_filtering.md`).*
**(a) A step with no bake-off is a file, not a folder.** Read filtering is a fixed
prelude with no competing implementations, so it is a single `read/filtering.rs`, not a
folder — the folder shape earns its keep only when a step has alternatives to sit side
by side. **(b) Tightly-coupled steps may share one folder.** Steps 1 (filtering) and 2
(read preparation) both turn a `MappedRead` into locus evidence and share the same input
type and reference accessor, so they live together in one `read/` module rather than in
two sibling folders. This bends "one folder per step" while keeping its intent: step 2's
`ReadPreparer` implementations still sit side by side within `read/`.

**Steps 6–9 are rule (b)'s second and larger case, and the test that decided it was a
dependency, not a feeling.** As four sibling folders at the top of `src/ng/`, the three
calling arch docs each had to state that the prior and the likelihood must never import
from `inference/`, and the shared types had to cross as flat slices partly to enforce it —
a constraint the tree imposed rather than the design. Under one `calling/` folder the
shared vocabulary sits in `calling/mod.rs` and every sub-module imports downward, so the
rule disappears. **The generalisation worth carrying: when keeping steps apart forces you
to write a no-import rule between them, they are one folder.** (The flat slices stay, for
their own reasons — the no-allocation contract, and `genetics.rs:127`'s deliberate
avoidance of a back-reference into its caller.)

**2. STR-ness is not a separate subtree.** An STR candidate generator is just
`calling/allele_candidates/rung_ladder.rs` sitting next to the generic `assembly.rs`, just
as `locus_generation/ssr.rs` sits next to `locus_generation/pileup/`; the region's kind decides which runs. We do **not**
split the pipeline into
`ssr/` vs `generic/` — that would scatter each step's variants across two trees and make
the per-step comparison awkward. STR-ness is a property of certain *implementations*, not
a top-level division. (STR domain *types* — `SsrMotif`, `SsrPeriod`, `SsrLocus` — still
carry the `Ssr` prefix and live in `units`/`locus` beside their generic peers.)

**3. Shared vocabulary starts in one `types.rs`, splits by concept as it grows.** The
cross-step types (the scalar newtypes + the wire structs) begin together in `types.rs`
rather than pre-split into many small files — the idiomatic Rust rhythm is to let a file
grow and extract a coherent chunk *when it tells you to*. `types.rs` is a common Rust
convention and an honest name for a mixed starting file. naming.md discourages a
*permanent* generic catch-all; we satisfy that by **splitting into concept modules** —
`units` (scalars), `locus`, `genotype`, `params` — as each cluster grows. Start minimal,
end concept-named.

**4. `bench/` is first-class.** The lab's whole purpose is *measuring* (spec §2), so the
standards harness — the representation-normalising truth comparison, the gold/silver/
synthetic scorers — is a real module the pipeline is built around, not a `tests/`
afterthought. A step's winner is decided here.

**5. Reuse over rewrite.** New modules, standing on existing code: `read/filtering.rs` wraps
the filters in [bam/alignment_input.rs](../../../../src/bam/alignment_input.rs),
`likelihood` builds on [ssr/cohort/read_model/](../../../../src/ssr/cohort/read_model/),
`inference` on [var_calling/posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs).
The §6 reconciliation table is the map of what to reuse.

## Anatomy of a step module

Every step folder follows the same shape — the trait in `mod.rs`, one file per
implementation, tests beside the code:

```
calling/allele_candidates/
├── mod.rs          – the CandidateGenerator trait + re-exports of the impls
├── rung_ladder.rs  – the STR repeat-length implementation  (+ #[cfg(test)] tests)
├── assembly.rs     – the generic local-assembly implementation (+ tests)
└── …               – further impls as the bake-off grows
```

New implementation of a step = a new file in that step's folder implementing the trait.
Nothing else in the tree changes — that locality is the point.

## How it assembles

- **`CallerRecipe`** (in `pipeline.rs`; defined in `ng_step_interfaces.md` §4) names one
  impl per step — it *is* one experiment ("freebayes candidates + HipSTR stutter
  likelihood + our cohort prior").
- **`pipeline.rs`** drives a recipe over its inputs: a per-sample stage producing each
  sample's loci, then the cohort gather that merges them — the artifact in between held in
  memory, not written (see *Crate boundary*).
- **`bench/`** scores the pipeline's output against the standards and reports the frontier.

So: the step folders provide the *parts*, the recipe *selects* a set, the pipeline *runs*
it, and bench *judges* it. Swapping one part and re-running is the unit of work.

### The locus stream — where `SampleLocusObservations` is born

`pipeline.rs` is orchestration only; the per-locus units it drives come from a **locus
stream** (`ng_proposal.md` §1, *The locus stream*). `locus_generation/` mints every locus, whatever
its kind, which is what keeps SNP / indel / STR at one level:

```
region_typing/  types the reference into TypedRegions (reference-based; spec ../spec/typed_regions.md)
   └─▶ locus_generation/ dispatches each region on its kind to the generator that handles it:
        ├─ SsrSegment → locus_generation/ssr.rs     → 1 locus, defined from the reference
        ├─ Generic    → locus_generation/pileup/    → many loci, split from the data
        └─ Satellite / SsrBundle → NoLoci → 0 loci, counted with a reason
             └─▶ one stream of SampleLocusObservations
                  └─▶ pipeline.rs feeds it to the per-locus core (steps 6–9)
```

**`locus_generation/` is a step folder** — it owns the `LocusGenerator<S>` trait and every
implementation of it, side by side, which is what principle 1 asks for. ng is an
experimental caller and more than one generator per kind is expected; the segment type
is a parameter on the trait so two generators for the *same* kind stay interchangeable
(`../spec/locus_generation.md` §4).

*This reverses this doc's earlier call* that `pileup/` was infrastructure with "no
swappable-trait bake-off surface of its own", sitting at the tree's top level beside
`pipeline.rs` and `bench/`. It has one: the pileup is a `LocusGenerator` like any other,
and it lives inside `locus_generation/` with its siblings. Recorded rather than quietly changed
because the old placement is what the tree above used to show.

The related open question is **resolved**: `pileup/` is **built from** step 2's
`ReadPreparer`, not a subsumer of it (`../spec/read_preparation.md` §6 — *compose, not
subsume*). The generic path does not opt out of the per-step bake-off.

## Crate boundary and the port-back

ng stays a module inside `pop_var_caller` (spec §3): one thread, reuse freely. The module
tree here is the *research* home; the production modules remain the *scaling* home.

**On the per-sample/cohort split — revised.** This doc originally said "no `.psp` split,
single-phase". ng does adopt production's two-level *shape* — per-sample stage → artifact →
cohort gather — because a `SampleLocusObservations` is one sample's and cohort loci are built by merging
many (`../spec/locus_generation.md` §3).

### The artifact between the two stages

**ng has its own on-disk store for that artifact, and it is written: `src/ng/psp/`.** This doc
said until 2026-08-30 that the artifact "starts in memory and gains a serialization when memory
forces it". Memory forced it: everything the calling stage holds is multiplied by the cohort size,
and the committed range runs to three thousand samples (`../spec/run_streaming.md` §7.1), so the
artifact has to be something each sample's reader streams from rather than something the run
holds. Production's `.psp` is the shape not to copy — its block index alone is 3.8 MB per open
file, 11.5 GB across three thousand samples (§7.2 there). ng's store is a different file
from production's `src/psp/` and shares no code with it; the two use the words *trailer* and
*block* for different things, which `src/ng/psp/mod.rs` spells out for anyone moving between
them.

What exists: a store can be created, pushed to, finished, opened, walked from the start or from a
coordinate, walked with a per-record predicate deciding which bodies get built, re-trailered and
appended to. It is measured — **480 kB an open sample on a human reference against a 500 kB
budget** — and checked field by field against a codec that is not ours over 7.7 million records.
Spec: `../spec/psp_file_format.md`.

**What does not exist yet is the wiring.** `pipeline.rs` still passes the artifact in memory;
nothing in the run writes a store or reads one, and only the measuring examples use the module.
Which run object writes and which opens is `../spec/run_streaming.md`'s to say — it owns the
three run objects — and that is the next step.

## Naming to confirm

- `read/` (steps 1+2) — the merged read-handling module (see principle 1). Step 2 is
  **generic-path only** (`../spec/read_preparation.md` §1); its files are named for the **transform
  they perform** (`left_align.rs`, …), not for the taxonomy pole they sit on — "trust the mapper"
  is an *axis* (`ng_proposal.md` §2) and names a family, not one implementation. The `ReadPreparer`
  trait lives in `read/mod.rs`; its v1 impl is `LeftAlignPreparer` (pass-through + left-align; BAQ is
  deferred sine die and the re-align mode is gated — `../spec/read_preparation.md` §10). (An
  `ssr_delimit.rs` / `SsrDelimitPreparer` and a `reassembly.rs` / `ReassemblyPreparer` were listed here
  until 2026-07-25 / 2026-07-23: the STR "preparer" is not a preparer at all — its per-read alignment
  produces a locus observation, so it lives in `locus_generation/ssr.rs` calling `alignment/`; and local
  reassembly is out of scope. `pair_hmm.rs` was listed too; the pair-HMM lives in `alignment/`.) The
  "prep" abbreviation is gone: a `ReadPreparer` does `prepare_read` and yields a `PreparedRead` (verb,
  agent noun, product).
- `types.rs` (the one shared-types file) — a common Rust convention, honest for a mixed
  starting file; naming.md leans against a *permanent* generic module, so the plan is to
  split it into concept modules (`units`/`locus`/`genotype`/`params`) as it grows
  (principle 3). If you'd rather avoid `types.rs` entirely, start with `units` + concept
  files instead.

## Open items

- **Where step 5 (STR read-class / spanning) lives** — it is STR-only and feeds candidate
  generation; likely a submodule of `calling/allele_candidates/` or of `locus_generation/ssr.rs`, not its own
  top-level step folder. Decide when the STR path is built.
- **How much of the production `pileup/walker/` lifts** into an in-memory context, versus a
  lean rewrite that calls its decompose/active-set core. Decide when `locus_generation/pileup/` is
  built. (*The subsume-or-compose half of this item is closed:* the pileup is **built from**
  step 2's `ReadPreparer`, per `../spec/read_preparation.md` §6.)
- **Feature-gating.** If ng grows heavy, gate it behind a `cargo` feature so the production
  build need not compile the lab. Decide once there is code to gate.
- **`bench/` vs the existing `benchmarks/` tree.** `benchmarks/` holds data + scripts; the
  ng `bench/` module is in-crate scoring code. Keep the boundary explicit so they don't
  blur.
