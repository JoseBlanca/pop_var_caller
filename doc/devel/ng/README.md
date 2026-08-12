# `ng` — next-generation caller docs

All documentation for the **ng** effort (the step-decomposed, benchmark-driven
algorithm lab — see [`spec/ng_proposal.md`](spec/ng_proposal.md)) lives here, split
by document kind:

- **`spec/`** — proposals and design specs (what to build and why).
  - [`ng_proposal.md`](spec/ng_proposal.md) — the plan: the step-decomposed caller
    taxonomy, the benchmark strategy, and the single-phase lab.
  - [`read_filtering.md`](spec/read_filtering.md) — step 1 (the whole-read keep/drop
    prelude) + the ng foundations it settles (skeleton, `types.rs` seed, conventions).
  - [`read_filtering_stages.md`](spec/read_filtering_stages.md) — dividing step 1 into two
    filters and a conversion, so `filtering.rs` is left holding only the keep-or-drop rules
    and the cursor owns the loop. Adds no types: it renames the raw read to say what it is,
    deletes the source trait and the three-way stop, and fixes the boundary — the filters and
    the conversion live *below* what the cursor keeps, or reads get converted a dozen times
    each. Which filters run, and their thresholds, stay in `read_filtering.md`.
    **Settled, no code yet.**
  - [`read_preparation.md`](spec/read_preparation.md) — step 2, the per-read, **locus-independent**
    transform (pass-through / canonicalize / re-align → `PreparedRead`). It is a **generic-path-only**
    step: the STR path has no read preparation — it goes filtering → observation generation, aligning
    each read against its tract in [`locus_generation_ssr.md`](spec/locus_generation_ssr.md). (The old
    `read_preparation_generic.md` and `read_preparation_ssr.md` are retired redirects.)
  - [`alignment.md`](spec/alignment.md) — the **alignment algorithms** both step 2 and step 7 call:
    best-path (one line-up) vs marginal (summed over all line-ups), affine vs repeat-aware, plus
    alignment normalization. Not a pipeline step — it knows no caller. Lists the seven algorithms
    to build and compare.
  - [`typed_regions.md`](spec/typed_regions.md) — step 3 (the typed-region generator): walks
    the reference and cuts it into `TypedRegion`s (SsrSegment / SsrBundle / Generic / Satellite).
    First integrating spec — stands on `RefSeq`, the tandem-repeat scanner, and the STR catalog.
  - [`typed_regions_cli.md`](spec/typed_regions_cli.md) — the **`pop_var_caller_exp`** binary
    (ng's command surface, kept out of the production CLI) and its first subcommand,
    `type-regions`: step 3's walk driven from the command line, writing the genome's
    partition to a text file (kind + motif + repeat count per region).
  - [`ref_seq.md`](spec/ref_seq.md) — the `RefSeq` reference-sequence accessor
    (foundational infra: resident + streaming + in-memory impls). Read filtering #8 and
    the pileup depend on it.
  - [`reference_info.md`](spec/reference_info.md) — reading a reference's info (contig table +
    content digest + reconstructed index) from a FASTA (optionally checked against a `.fai`) or
    a `.fai` alone; with the MD5s, the fasta↔fai check, a caller-held cache, and a `.fai` writer.
    Builds the `ContigList` every `RefSeq` impl takes but none can build.
  - [`alignment_file.md`](spec/alignment_file.md) — **one** alignment file: open-and-validate
    (`SO`, `@SQ`↔reference, index, `@RG SM`) and serve a region as an ordered, filtered read
    stream. Closes the `@SQ` permutation hole and the sort-order check; step 1's input edge.
  - [`sample_reads.md`](spec/sample_reads.md) — **the sample**: k files (usually several
    experiments) merged into one coordinate-ordered stream, with the cross-file checks. Stands
    on `alignment_file.md`.
  - [`read_groups.md`](spec/read_groups.md) — **the read group** (`@RG`) as a first-class object:
    parsed once, identified run-wide, stamped on every read. Makes the library visible to the error
    model, which fits per-base error and stutter per chemistry — a property of the library prep, not
    of the individual. Changes decisions in `alignment_file.md` and `sample_reads.md`.
  - [`locus_generation.md`](spec/locus_generation.md) — the shared shape of locus generation:
    typed regions → **a sample's loci** (`SampleLocusObservations`), the `LocusGenerator` contract,
    the dispatcher. Joins the typed-region walk to read ingestion; ships only the `NoLoci` generator.
  - [`locus_generation_ssr.md`](spec/locus_generation_ssr.md) — the first generator (STR): one
    tract segment → one locus, adapting production `src/ssr/pileup/` and carrying partial (censored)
    observations the old path dropped.
  - **Step 4, the parameter pre-pass, in eight documents** — the noise rates and rates of variation
    a caller runs on, estimated **without first calling genotypes**. Five describe the **per-sample
    route** (walk each sample, fold its loci into histograms, fit from those); the last three
    describe a **second, complete route to the same parameters** (keep raw evidence at the same loci
    in every sample, fit everything once when they are all in). Both are built, and their estimates
    compared:
    [`parameter_prepass.md`](spec/parameter_prepass.md) (read first: what the step produces, why
    production's numbers are biased, the estimator both paths share);
    [`parameter_prepass_generic.md`](spec/parameter_prepass_generic.md) (the SNP/indel path: two
    histograms, the per-base error rate, the sample's rates, the inbreeding coefficient);
    [`parameter_prepass_ssr.md`](spec/parameter_prepass_ssr.md) (the STR path: what stutter
    actually looks like, and the per-locus table its four numbers are fitted from);
    [`parameter_prepass_census_sites.md`](spec/parameter_prepass_census_sites.md) (the two
    censuses — the same loci in every sample, so answers can be compared); and
    [`parameter_prepass_cohort.md`](spec/parameter_prepass_cohort.md) (the gather: diversity, the
    frequency spectrum, contamination, relatedness).
    Then the second route, in three documents that map to three modules:
    [`parameter_prepass_joint_fit.md`](spec/parameter_prepass_joint_fit.md) (**read first** — what the
    route is, what it produces, and the estimator: each locus weighted by its own allele frequency in
    the cohort, so the fit runs once over every sample rather than per sample);
    [`parameter_prepass_joint_loci.md`](spec/parameter_prepass_joint_loci.md) (which loci every sample
    keeps evidence at — uniform for ordinary sites, **equal per repeat-count stratum** for STR loci);
    and [`parameter_prepass_joint_records.md`](spec/parameter_prepass_joint_records.md) (what is
    recorded at each, and how it is encoded).
  - [`repeat_catalog.md`](spec/repeat_catalog.md) — **the reference's tandem-repeat catalog**, built
    inside the pass that already streams the whole FASTA (`reference_info`), by a command whose only
    job is that, so a genome is scanned once instead of once per sample. It exists because the
    parameter pre-pass's **random sample of STR loci per repeat-count stratum** needs the genome
    enumerated first. It records **repeats, not loci** — both spans, period, score, motif and purity —
    so a reader derives the segmentation, and the STR loci under any copy floor it chooses, without
    opening the FASTA. Independent work; step 3 is designed for as a future reader, not wired.
    Companions: [`arch/repeat_catalog.md`](arch/repeat_catalog.md) (types & interfaces) and
    [`impl_plan/repeat_catalog.md`](impl_plan/repeat_catalog.md) (the build order).
  - [`synthetic_validation.md`](spec/synthetic_validation.md) — the generated data the calling
    steps are graded against.
- **`research/`** — measurements, kept apart from the designs they settled so a spec can point at
  a number rather than repeat it.
  - [`parameter_estimator_experiments_2026-08-06.md`](research/parameter_estimator_experiments_2026-08-06.md)
    — what step 4's estimators actually do, from three harnesses in `examples/`. Bias computed
    **exactly** rather than simulated, so "unbiased" is decided rather than estimated. It carries
    the numbers behind the multi-library cell key, the inbreeding coefficient, depth binning, the
    heterozygote `½`, and the STR stutter accumulator — and a list of the findings it overturned,
    including several of its own.
- **`arch/`** — architecture (the shared types and the interfaces implementations
  plug into).
  - [`ng_step_interfaces.md`](arch/ng_step_interfaces.md) — the common domain
    newtypes and one swappable trait per pipeline step.
  - [`module_layout.md`](arch/module_layout.md) — the `src/ng/` module tree: one
    folder per step (trait + impls + tests together), shared vocabulary, `bench/`.
  - [`read_filtering.md`](arch/read_filtering.md) — step 1's types & interfaces,
    distilled (the code-facing companion to the spec).
  - [`read_filtering_stages.md`](arch/read_filtering_stages.md) — the renames, the two verdict
    functions, and the loop as the cursor owns it; companion to
    `spec/read_filtering_stages.md`.
  - [`alignment.md`](arch/alignment.md) — the alignment module's types & interfaces
    (`BestPathAligner`, `MarginalAligner`, `AlignmentNormalizer`, `RepeatSpan`, `StutterModel`;
    seeds `LogProb`); companion to `spec/alignment.md`. Called by two steps, not a step itself.
  - [`typed_regions.md`](arch/typed_regions.md) — step 3's types & interfaces (the
    typed-region generator); companion to `spec/typed_regions.md`.
  - [`typed_regions_cli.md`](arch/typed_regions_cli.md) — the `pop_var_caller_exp` /
    `type-regions` types & interfaces (the CLI that drives step 3); companion to
    `spec/typed_regions_cli.md`.
  - [`reference_info.md`](arch/reference_info.md) — the reference-info reader's types &
    interfaces (`ReferenceInfo`/`ContigInfo`, the cache, the writer, the background verify);
    companion to `spec/reference_info.md`. Foundational infra, not a step.
  - [`alignment_file.md`](arch/alignment_file.md) — `AlignmentFile` (the validated handle) and
    the validate-on-open gate; companion to `spec/alignment_file.md`. **Its region-query half —
    the per-region `RecordSource` impls, the order guard and the reader pool — is superseded by
    `alignment_cursor.md` and was deleted**; that part of both documents is a design record.
  - [`alignment_cursor.md`](arch/alignment_cursor.md) — the long-lived reader that stays
    positioned in one chromosome of one file and keeps the reads it has already decoded and
    filtered: `AlignmentCursor`, `SampleCursor`, the per-format `AlignedReadsReader`s, and the forget
    rule. **The only way to read a BAM or a CRAM in ng.** Companion to
    `spec/alignment_cursor.md`.
  - [`sample_reads.md`](arch/sample_reads.md) — `SampleReads`, the argmin merge and its
    per-read budget; seeds the shared `GenomePosition`. Companion to `spec/sample_reads.md`.
    **The merge it specifies now lives in `sample_cursor.rs` over cursors** rather than over
    per-region streams; the k-way rules are unchanged and the entry point is not.
  - [`read_groups.md`](arch/read_groups.md) — `ReadGroup`/`ReadGroups`, the run-wide `ReadGroupId`,
    the per-open `ReadGroupResolution`, and ng's own `AlignedRead`; companion to
    `spec/read_groups.md`.
  - [`locus_generation.md`](arch/locus_generation.md) — the shared locus-generation types &
    interfaces (`SampleLocusObservations`, `ObservedSequence`, `LocusGenerator<S>`, the dispatcher,
    `NoLoci`); companion to `spec/locus_generation.md`.
  - [`locus_generation_ssr.md`](arch/locus_generation_ssr.md) — the STR generator's types &
    interfaces (`SsrLocus`, `SsrGenerator`, the reservoir cap and its traps); companion to
    `spec/locus_generation_ssr.md`.
  - **The per-site route to step 4's parameters, in three companions to the three specs** —
    [`parameter_prepass_joint_loci.md`](arch/parameter_prepass_joint_loci.md) (`KeptLoci`,
    `SelectionIdentity`, the kept-loci digest; the STR half is two calls on the repeat catalog),
    [`parameter_prepass_joint_records.md`](arch/parameter_prepass_joint_records.md)
    (`SampleRecords`, the five-bit depth code, the STR difference list, `RecordWriter`) and
    [`parameter_prepass_joint_fit.md`](arch/parameter_prepass_joint_fit.md) (`JointFit`, the
    three site classes, `HomozygoteExcess` beside `InbreedingF`, contamination as a value or a
    stated reason there is none). All three land in `src/ng/parameter_estimation/joint/`.
- **`impl_plan/`** — step-by-step implementation plans (build order, not new design).
  - [`foundations.md`](impl_plan/foundations.md) — the first ng code: skeleton,
    `types.rs` seed, and the `RefSeq` accessor (three impls).
  - [`read_filtering.md`](impl_plan/read_filtering.md) — step 1: the `read/` module,
    the cascade, the `RecordSource`/`RawAlignedRead` seam, the `ReadFilter` iterator.
  - [`read_filtering_stages.md`](impl_plan/read_filtering_stages.md) — dividing step 1 into two
    filters and a conversion: the renames, the contig check as a table comparison, the loop
    moving into the cursor, and the two tests output identity cannot see.
  - [`read_input.md`](impl_plan/read_input.md) — step 1's input edge (`read/input/`): the
    validate-on-open gate, the BAM/CRAM region queries, the order guard, and the k-file
    merge. Covers both `alignment_file` and `sample_reads`.
  - [`read_groups.md`](impl_plan/read_groups.md) — the read-group table, opening from it, the
    ng-owned `AlignedRead` that carries the identifier, and per-read-group counts. Modifies
    `read_input.md`'s `AlignmentFile`/`SampleReads`.
  - [`typed_regions.md`](impl_plan/typed_regions.md) — step 3: the catalog rebase/knobs,
    the windowed substrate, and the `region_typing.rs` walk (resident → windowed).
  - [`typed_regions_cli.md`](impl_plan/typed_regions_cli.md) — the `pop_var_caller_exp` /
    `type-regions` build order: binary skeleton, `--min-copies` parser, output writer, and
    the `run_typed_regions` driver.
  - [`reference_info.md`](impl_plan/reference_info.md) — the reference-info reader: types →
    `.fai` reader → FASTA pass (the heart) → writer → cache → the two entry points.
  - [`locus_generation.md`](impl_plan/locus_generation.md) — the shared locus-generation shape:
    the `locus_generation/` module, the locus types, `LocusGenerator<S>` + `NoLoci`, the
    dispatcher, and the iterator (proven with the count-only generator). Ships no real generator.
  - **The alignment module, in three plans, in this order:**
    [`alignment_best_path.md`](impl_plan/alignment_best_path.md) (the module skeleton, the aligner
    types, the repeat delimiter that **unblocks the STR generator**, its banding, the two-penalty
    comparison, and a gated affine aligner);
    [`alignment_marginal.md`](impl_plan/alignment_marginal.md) (`LogProb`, the marginal interface,
    the sequence marginal and the whole-read forward); and
    [`alignment_normalization.md`](impl_plan/alignment_normalization.md) (the normalizer interface,
    three left-aligners, and the property test that grades them against a definition).
  - [`locus_generation_ssr.md`](impl_plan/locus_generation_ssr.md) — the STR generator: the
    prerequisite `flank_bp`→`bundle_threshold` rename, `SsrLocus` + margin fetch, the ported
    reservoir cap, the fetch→align→tally transform, and byte parity vs production. Gated on the
    ng STR aligner.
  - [`parameter_prepass_generic.md`](impl_plan/parameter_prepass_generic.md) — step 4's
    SNP/indel half: the `parameter_estimation/` module, the depth ladder and cell table, the
    two keyed accumulators, the fitting machinery both paths share, and the four numbers a
    sample emits. **Its oracle is not production** — nothing downstream can check these
    parameters and the production code nearest to them is what the step replaces — so every
    milestone is proven against the two research harnesses in `examples/` or against an
    identity. The STR half's plan follows its architecture doc settling; the two censuses and
    the cohort gather still need one.

This mirrors the repo-wide `doc/devel/{specs,architecture,implementation_plans}`
convention but scoped to ng, so the growing set of ng docs stays together.
