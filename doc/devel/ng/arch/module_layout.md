# Module layout

*Status: describes the tree as it stands on 2026-11-27. First drafted on 2026-07-10, when the
caller was being built as `src/ng/` beside an older one; it was revised to the current layout
after the older caller was deleted (2026-09-13) and these modules moved up to the crate root.
Companion to [`ng_step_interfaces.md`](ng_step_interfaces.md), which defines the shared types and
step traits; this one says where they live, and the rules that decided it.*

## Where the caller lives

The caller is the library crate `pop_var_caller`, with its modules at the crate root
(`src/lib.rs` lists them), and the binary `pvc` (`src/main.rs`), whose subcommands live in
`src/cli/`. The design documents keep the name the caller was written under, "ng"
(`doc/devel/ng/`), and a comment that says "ng" means this caller. A comment that says
"production" means the deleted older caller, whose source is whole at commit `74e424ea`.

## The tree

The step numbers are those of the step map in [`ng_proposal.md`](../spec/ng_proposal.md)
("The step map").

```
src/
├── lib.rs              – the module list and the crate doc
├── main.rs             – the `pvc` binary
├── cli/                – one file (and folder, where it grew) per subcommand, plus the
│                         flags and the run set-up the calling commands share (calling_run/)
│
│   # the inputs
├── bam/                – BAM and CRAM: the read vocabulary (CIGAR, flags, filter thresholds)
│                         and the index check; the reader itself is read/input/
├── fasta/              – reference FASTA readers and the contig list
├── regions.rs          – region lists (BED) and the region set a run walks
├── reference_info.rs   – the reference's contig table and digest, as recorded in outputs
├── ref_seq.rs          – reference-sequence access: the RefSeq trait and its implementations.
│                         Shared infrastructure, not a step. Spec ../spec/ref_seq.md
├── raw_chrom_reader.rs – a contig's bases read straight from the FASTA
│
│   # step 3: what each stretch of the reference is
├── tandem_repeat.rs    – the tandem-repeat scanner
├── repeat_catalog/     – every tandem repeat in the reference, found once and written beside
│                         the FASTA; a run reads its typed regions from this file
├── region_typing/      – the typed-region vocabulary (TypedRegion, RegionKind) and the policy
│                         that decides between them. Spec ../spec/typed_regions.md
├── segmentation_inputs.rs – the settings that decide the typing
│
│   # steps 1 and 2: reads
├── read/               – the alignment-file reader (input/), a decoded read (AlignedRead), the
│                         filtering prelude (step 1, filtering.rs) and read preparation for the
│                         generic path (step 2: the ReadPreparer trait, LeftAlignPreparer)
├── alignment/          – NOT a step: the alignment algorithms step 2 and the STR path call.
│                         Knows nothing of the steps that call it. Spec ../spec/alignment.md
│
│   # one sample's loci
├── locus_generation/   – the LocusGenerator trait and one generator per region kind: ssr/ (STR
│                         segments and bundles, one locus each, defined from the reference),
│                         pileup/ (generic sequence, loci split from the data), NoLoci
│                         (satellites, counted).
│                         Spec ../spec/locus_generation.md
├── psp/                – the per-sample store: one file per sample holding what its reads
│                         showed, written by `generate-psps` and read back by the cohort merge.
│                         Not a step — the seam between the per-sample and cohort stages, made
│                         durable. Spec ../spec/psp_file_format.md
├── window_coverage/    – each sample's read depth and GC fraction by window
│
│   # step 4: the parameters the caller runs on
├── parameter_estimation/ – error rates, heterozygosity, inbreeding and stutter, fitted from the
│                         cohort's own loci before calling; joint/ holds the census and the
│                         whole-cohort fits. Specs ../spec/parameter_prepass*.md
├── misread_context.rs  – which misread-after-a-run table entry a base belongs to: the one rule
│                         the parameter fit counts by and calling charges by
│
│   # the cohort
├── cohort_observation.rs – what the cohort merge hands calling: one cohort locus and each
│                         sample's evidence over it
├── calling/            – steps 6–9 and 13, one folder because they are one question (principle
│                         1b): allele_candidates/ (step 6), likelihood/ (step 7, the SsrEmissionModel
│                         seam), genotype_prior/ (step 8, GenotypePriorModel), inference/ (step 9,
│                         LocusGenotyper), quality/ (step 13); parameters_file/ is the file the
│                         fitted parameters are written to and read from
├── genetics.rs         – the population-genetics arithmetic the genotype prior, site quality
│                         and duplication filter share
├── paralog/            – step 11a: the hidden-duplication filter's statistics
├── neighbour_evidence/ – refinement: re-calling a sample with the evidence of the samples it
│                         shares stretches of chromosome with, and the phase written from it
├── read_phase/         – step 10: phase from reads
├── run/                – what a run is made of: the walk over each sample's alignments, the psp
│                         files it can write, the parameter fit, the cohort merge
│                         (cohort_merge/), the duplication filter's wiring (paralog_filter/) and
│                         the calling that follows
├── vcf/                – the one file a run writes: record assembly, field encoding, header,
│                         writer, and the stage that writes phase from reads
├── types.rs            – the shared vocabulary (principle 3)
│
│   # small helpers
├── float.rs, float/    – transcendental functions computed in Rust, the same bits on every platform
├── splitmix.rs         – the hash that thins reads by a shared key
├── memory_limit.rs, memory_limit/ – the most memory the process may hold
├── open_files.rs       – the most files the process may open
├── error_render.rs, progress_text.rs – how an error and a progress line are written
│
│   # test-only
└── golden_catalog.rs, production_post_filter.rs, scanner_parity.rs – recorded output of the
                          deleted caller and of `trf-mod`, kept as fixtures the tests compare
                          against
```

Step 5 (which reads span an STR) has no module of its own: the STR generator classifies each read
against its tract (`locus_generation/ssr/`), and candidate selection decides which partial reads
count (`calling/allele_candidates/deciding_partial_reads.rs`). Step 11b, which loci are written,
is decided in `run/records.rs` from the calls and calling's site quality.

## Organizing principles

**1. One module per pipeline step — trait, implementations and tests together.** A step folder
holds its trait, every competing implementation and their tests, so a step's alternatives sit
side by side and one can be swapped while the rest is held fixed. Modules are named for the
**concept** they own (`calling`, `allele_candidates`, `genotype_prior`), never for a layer or
pattern (`models`, `services`, `common`, `utils`), following the naming rule in
[naming.md](../../../../ai/skills/rust-code-review/code_review/naming.md).

Two rules of thumb the read-filtering spec pinned down (`../spec/read_filtering.md`):

- **(a) A step with no competing implementations is a file, not a folder.** Read filtering is a
  fixed prelude, so it is a single `read/filtering.rs`.
- **(b) Tightly coupled steps may share one folder.** Steps 1 (filtering) and 2 (read
  preparation) both turn a decoded read into locus evidence and share the same input type and
  reference accessor, so they live together in `read/`.

**Steps 6–9 are rule (b)'s second and larger case, and the test that decided it was a
dependency.** As four sibling folders, the three calling design documents each had to state that
the prior and the likelihood must never import from `inference/`, and the shared types had to
cross as flat slices partly to enforce it — a constraint the tree imposed rather than the design.
Under one `calling/` folder the shared vocabulary sits at the top of `calling/` and every
sub-module imports downward, so the rule disappears. **When keeping steps apart forces a
no-import rule between them, they are one folder.**

**2. STR-ness is not a separate subtree.** The STR implementation of a step sits next to the
generic one in that step's folder — `locus_generation/ssr/` beside `locus_generation/pileup/`,
`calling/likelihood/ssr.rs` beside `generic.rs` — and the region's kind decides which runs. The
pipeline is not split into `ssr/` and `generic/` trees, which would scatter each step's variants
across two places. STR domain *types* (`Motif`, `SsrPeriod`) still carry their STR names and live
in `types.rs` beside their generic peers.

**3. Shared vocabulary lives in one `types.rs`.** The cross-step types — the small newtypes that
more than one step passes around — are together in `src/types.rs`, which reads as a dictionary of
short, self-contained types (its module doc lists the five groups). A type moves there when it
gains consumers in more than one step: `Motif` started in `region_typing::segment_criteria` and
moved once `locus_generation` and `alignment` needed it. The file is not split by group, because
every group would then import the shared error type from a sixth file; revisit that when one group
passes about 600 lines of code.

## The locus stream

`locus_generation/` mints every locus, whatever its kind, which keeps SNP, indel and STR at one
level:

```
repeat_catalog/ + region_typing/  type the reference into TypedRegions
   └─▶ locus_generation/ dispatches each region on its kind to the generator that handles it:
        ├─ SsrSegment → locus_generation/ssr/     → 1 locus, defined from the reference
        ├─ SsrBundle  → locus_generation/ssr/     → 1 locus over the bundle's hull
        ├─ Generic    → locus_generation/pileup/  → many loci, split from the data
        └─ Satellite  → NoLoci                    → 0 loci, counted with a reason
             └─▶ each sample's loci, written to its psp file or handed on in memory
                  └─▶ run/cohort_merge/ joins every sample's loci into cohort loci
                       └─▶ calling/ (steps 6–9, 13) └─▶ vcf/
```

The segment type is a parameter on the `LocusGenerator` trait, so two generators for the same
kind stay interchangeable (`../spec/locus_generation.md` §4). The pileup generator is built from
step 2's `ReadPreparer`, not a replacement for it (`../spec/read_preparation.md` §6).

## The per-sample and cohort stages

A run has two stages with an artifact between them: each sample's reads become that sample's
loci, and the cohort merge joins every sample's loci into cohort loci. The artifact can be a file
— the psp store (`src/psp/`), one per sample, written by `pvc generate-psps` and read by
`pvc call-from-psps` — or, in `pvc call-from-alignments`, it is handed from the walk to the merge
in memory. Both routes run the same merge and calling code (`run/`).

The store exists because everything the cohort stage holds is multiplied by the cohort size, and
the committed range runs to thousands of samples (`../spec/run_streaming.md` §7.1), so each
sample's evidence has to be something the merge streams from rather than something the run holds.
