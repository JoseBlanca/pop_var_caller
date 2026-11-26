# Cross-module map — dependencies, duplication, dead code (2026-11-26)

Scope: the whole of `src/`, looking only at what crosses module boundaries. Read-only; nothing
was built or run. The dependency graph was extracted by a script that follows the `mod` tree from
`src/lib.rs`, drops `#[cfg(test)]` items and test-only files, and records every `use crate::…`,
inline `crate::…` path and `super::` path that leaves a top-level module. The script and its
output are in `tmp/cq_cross/` (`deps.py`, `edges.json`, `dead.txt`, `testdups.txt`) if a fixer
wants to rerun it after a change.

## Summary

1. **The layering is mostly sound.** 31 top-level modules, 132 module-to-module edges carrying
   1,852 imports. Removing **8 edges (about 60 imports)** makes the graph acyclic; today they tie
   12 of the 31 modules into two cycles. The worst is `calling` importing the cohort merge's
   output types from `run` (23 imports in 8 files): the calling library depends on the driver.
   Each of the 8 has a cheap fix (move a type down, or move one function up); none changes output.
2. **About 2,600 lines of dead production code in `bam/`**: `segment_reader`, `segment_merge`,
   `bam_input` and `cram_input` have no non-test caller and silence the compiler with a
   module-wide `allow(dead_code)`. Their comments describe a driver deleted with the old caller.
3. **Two module-wide dead-code silencers cover live code**: `calling/allele_candidates/ssr.rs`
   (3,839 lines, the STR candidate selection every run uses) and `read/input/aligned_reads_reader/`
   (the live BAM/CRAM readers). Each hides two or three test-only items and would hide any new
   dead code too.
4. **Duplicated utilities, all small**: slice log-sum-exp (3 identical copies), MD5-to-hex
   (7 copies), progress-line formatting (thousands separators, durations, byte sizes — 9 helpers in
   5 files), contig-name-to-id lookup (3 production copies), `ContigList`→`ContigBounds` (6 copies,
   two policies). Merging the identical ones changes no output. Three near-copies differ in
   floating-point order and **would change output if merged** — listed so nobody merges them.
5. **One wrong message**: `estimate-contamination` reports a contig longer than 4 Gb as "the
   reference has more contigs than a run can address".
6. **Test helpers**: about 640 lines of byte-identical helper bodies repeated across test modules,
   led by a `GenomeRegion` constructor defined in 24 files. The three `test_fixtures.rs` files are
   layered, not duplicated — leave them.

---

## Bugs

### B1 — A too-long contig in `estimate-contamination` is reported as too many contigs

- **Where**: `src/cli/estimate_contamination.rs:232-233` (the variant), `:602-613` (its only use).
- **What**: the conversion of each contig's `u64` length to the `u32` the BED reader takes maps a
  failed `u32::try_from(entry.length)` to `ContigCount`, whose message is "the reference has more
  contigs than a run can address". The real failure is one contig longer than 4,294,967,295 bases.
  The three sibling commands (`typed_regions.rs:512`, `inspect_psp.rs:211`, `run_ground.rs:224`)
  each name the contig and its length.
- **Scenario**: a reference holding one contig over 4 Gb (some conifer and amphibian assemblies) →
  the command stops with a message that sends the user looking at the contig count.
- **Proposal**: a `ContigTooLong { name, length }` variant, or the shared conversion of D4 below.
- **Downside**: none.
- **Effect**: error text only. **Cost**: S. **Confidence**: high (likely — traced end to end).

### B2 — Contig name → id is a linear scan, called once per selected STR locus (suspected, speed only)

- **Where**: the closures at `src/run/gatherer.rs:310`, `src/cli/estimate_parameters.rs:955`,
  `src/cli/estimate_contamination.rs:772`; called per locus at
  `src/parameter_estimation/joint/census.rs:2516` and `src/parameter_estimation/joint/ssr_fit.rs:1149`.
- **What**: each call walks `ContigList::entries` comparing names. `CensusWriter::new` calls it once
  per STR locus the census selected, once per sample. On a reference of a few chromosomes this is
  free; on a scaffold-level assembly with 30,000 contigs each lookup averages 15,000 string
  comparisons, so a census of 100,000 STR loci costs about 1.5 billion comparisons a sample.
  `psp/segmentation_section.rs:432` does the same scan once per stored span.
- **Proposal**: fold into D5 below — one `ContigList::id_of`, backed by a name→index map built once.
- **Downside**: a map per contig list (30,000 entries is about a megabyte).
- **Effect**: none on output; speed gain only on fragmented references. **Cost**: S.
- **Confidence**: medium — the call path is certain, the size of the cost depends on how many STR
  loci a census keeps, which was not measured.

---

## Architecture

### The layering this crate already almost has

Lowest first. A module may import from any layer below its own.

| Layer | Modules | What they are |
|---|---|---|
| 0 | `types`, `float`, `splitmix`, `error_render`, `open_files`, `memory_limit`, `regions`, `fasta`, `bam`, `genetics` | vocabulary, numerics, raw inputs |
| 1 | `raw_chrom_reader`, `ref_seq`, `reference_info`, `tandem_repeat` | the reference |
| 2 | `region_typing`, `repeat_catalog`, `segmentation_inputs` | what each stretch of reference is |
| 3 | `alignment`, `read` | reads turned into evidence |
| 4 | `locus_generation`, then `psp`, `window_coverage` | one sample's loci; storing and measuring them |
| 5 | `parameter_estimation`, `paralog`, then `calling` | fitting, and calling |
| 6 | `read_phase`, `neighbour_evidence`, `vcf` | refinement and the output file |
| 7 | `run` | the drivers |
| 8 | `cli` | the command line |

With the 8 edges below cut, a topological sort of the real graph reproduces this table (it puts
`read_phase` lower, since without its `vcf` import it needs only `float`, `splitmix` and `types`).
The precedent for the fixes is already in the tree: `segmentation_inputs` is a top-level module
"because both sides of the psp seam hold it" (`src/segmentation_inputs.rs:8-12`).

### A1 — `calling` imports the cohort merge's output types from `run`

- **Where**: `src/calling/allele_candidates/mod.rs:44-46`, `allele_candidates/ssr.rs:46-47`,
  `allele_candidates/generic.rs:15`, `allele_candidates/deciding_partial_reads.rs:47`,
  `evidence_shaping.rs:77`, `inference/mod.rs:67`, `likelihood/mod.rs:142`,
  `likelihood/generic.rs:13` — 23 imports in 8 files.
- **What**: `CohortObservation`, `SampleSupport`, `SupportedAllele`, `PartialObservation`,
  `LocusBorder` are defined at `src/run/cohort_merge/build.rs:1289-1688` (about 400 lines with
  their docs); `MinAltObs`, `MinAltReadShare`, `MinAltReads` at `src/run/cohort_merge/mod.rs:574-790`.
  They are calling's *input*: the merge produces them and calling consumes them. Because they live
  in `run`, the calling library cannot be read, tested or reused without the driver, and
  `calling ↔ run` is a cycle.
- **Proposal**: move the struct and enum definitions (and `MinAlt*`, which both the merge's keep
  rule and candidate selection read) to a new top-level `cohort_observation` module at layer 4/5.
  Leave the impl blocks that build them from a `ClosedLocus` in `build.rs` — an inherent impl may
  sit anywhere in the crate. `CohortObservation` carries `Vec<WindowCoverage>`, so the new module
  imports `window_coverage` (layer 4), which is fine once A3 is done.
- **Downside**: a reader of `build.rs` now opens a second file to see the fields it fills. The
  types are long-documented, so the split is along a natural seam.
- **Effect**: none (pure move). **Cost**: M (many `use` lines; mechanical). **Confidence**: high.

### A2 — `run` and `parameter_estimation` import progress formatting from `cli`

- **Where**: `src/run/calling_progress.rs:23` (`PROGRESS`, `elapsed_text`, `with_thousands`);
  `src/parameter_estimation/progress.rs:25` (`PROGRESS`).
- **What**: the progress target constant and two formatters live in `src/cli/run_output.rs:44`,
  `:480`, `:495`. Library code reaching into `cli` is the edge the owner named as a prime target.
- **Proposal**: a small layer-0 module (`progress_text`, say) holding `PROGRESS` and the
  formatters; `cli::run_output` keeps the subscriber and the log file and imports them. This is
  also the natural home for the duplicates in D2.
- **Downside**: one more small module.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### A3 — `window_coverage` imports the merge's cache types from `run`

- **Where**: `src/window_coverage/depth.rs:54` (`Drawn`, `LocusSummary`, `StandIn` from
  `run::cohort_merge::observation_cache`).
- **What**: `for_each_reported_depth` takes a `LocusSummary` (defined
  `src/run/cohort_merge/observation_cache.rs:128`), and `depth.rs:146` implements
  `From<&Drawn> for EvidenceForOneRecord`. Cycle: `run → window_coverage → run`.
- **Proposal**: move `LocusSummary` to `locus_generation` beside `SampleLocusObservations`, which it
  summarises (or into the new `cohort_observation` module of A1); move the `From<&Drawn>` impl into
  `observation_cache.rs`, which already imports `window_coverage`.
- **Downside**: `LocusSummary`'s doc says its 32-byte layout was measured (12% of the calling pass
  when it grew to 40). A move keeps the layout; do not "tidy" its fields into a `GenomeRegion` while
  moving it.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### A4 — `parameter_estimation` and `calling` import each other

- **Where**: `src/parameter_estimation/misread_after_runs.rs:17` and
  `src/parameter_estimation/joint/census_file.rs:57` import
  `calling::likelihood::run_context` (`MisreadTable`, `contexts_along`, `ENTRY_COUNT`, …).
  The other way, `calling` imports `parameter_estimation` 97 times, 65 of them from
  `calling/parameters_file/` and `calling/run_parameters.rs`.
- **What**: `run_context.rs` (547 lines) is the definition of a misread-after-a-run table that the
  pre-pass counts and calling charges (its doc says so, lines 12-14). It is shared vocabulary that
  happens to sit in the consumer.
- **Proposal**: move `src/calling/likelihood/run_context.rs` down into `parameter_estimation`
  (calling already depends on it), or to a top-level `misread_context` module at layer 3.
- **Downside**: calling's likelihood code imports its table from a module named for estimation; a
  top-level module avoids that oddity at the price of one more module.
- **Effect**: none (same crate; inlining unaffected). **Cost**: S. **Confidence**: high.

### A5 — `read` imports the psp header's parameter type

- **Where**: `src/read/filtering.rs:124` (`provenance_parameters`), `:174` (`psp::header::wire_float_of`).
- **What**: `ReadFilterConfig::provenance_parameters` renders the read filters as psp-header
  parameters. That one function closes the cycle `read → psp → locus_generation → read`.
- **Proposal**: move the function to the psp side — a free function `read_filter_parameters(&ReadFilterConfig)`
  in `psp/` or in `run/` where the psp writer is fed. It keeps its exhaustive destructure of the
  struct, so a filter added later still fails to compile until it is recorded.
- **Downside**: the rendering no longer sits beside the struct it renders.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### A6 — `vcf` and the two refinement stages import each other (leave it)

- **Where**: `src/vcf/writer.rs:36-41`, `:707-764` (the writer holds a `NeighbourEvidenceStage` or
  `ReadPhaseStage`); `src/vcf/header.rs:20`, `:772-777` (prints `read_phase` constants);
  the other way `src/neighbour_evidence/line.rs:22-23`, `refine.rs:50`, `stage.rs:42`,
  `src/read_phase/stage.rs:36-37` (`LineFields`, `parse_genotype`, `RecordPlace`, encoders).
- **What**: `vcf` is two things — the line format (`line_fields`, `encode`, `RecordPlace`) that the
  stages rewrite, and a writer that drives the stages. The cycle runs between those two halves.
- **Proposal**: leave it. If `vcf/writer.rs` is reworked anyway: let `VcfWriter` take an
  `Option<Box<dyn LineStage>>` (one stage at most, as the current enum guarantees) built by `run`,
  and hand the header the constants as metadata.
- **Downside of acting**: moves the "which stage runs" decision out of `VcfWriter::create`, where
  `writer_reads_fragments` (`writer.rs:134`) and the writer's tests hold the two answers together.
- **Effect**: none. **Cost**: M. **Confidence**: medium.

### Edges that look wrong and are not

- `alignment → bam`, `locus_generation → bam`: only `CigarOp` and SAM flag constants
  (`bam/alignment_input.rs`). `bam` is layer 0; these point down.
- `calling → read`, `calling → repeat_catalog`, `calling → region_typing` from
  `calling/parameters_file/`: the parameters file records read groups and catalog criteria. Down.
- `paralog → window_coverage` (`paralog/coverage_model.rs:58`): the duplication filter's model is
  built on coverage histograms. Down.
- `calling/parameters_file/` and `calling/run_parameters.rs` (together about 25,000 lines) are the
  parameters-file schema rather than calling. Promoting them to a top-level module would read
  better but cuts no cycle; not worth the churn on its own.

---

## Duplication

### D1 — `ln Σ exp` over a slice, three identical copies

- **Where**: `src/calling/quality/mod.rs:856` (`log_sum_exp_over`),
  `src/parameter_estimation/joint/contamination.rs:1404` (`ln_sum_exp`),
  `src/parameter_estimation/joint/fit.rs:1586` (`ln_sum_exp`).
- **What**: same two passes, same `−∞` guard, same `float::exp`/`float::ln`. 7 lines each.
- **Proposal**: one `float::ln_sum_exp(&[f64])`.
- **Downside**: none.
- **Effect**: none (same operations in the same order). **Cost**: S. **Confidence**: high.
- **Not to merge**: the pairwise `log_sum_exp_2` (`calling/genotype_prior/dirichlet_multinomial.rs:361`,
  already shared crate-wide) and `paralog/locus_score.rs:700` `log_add_exp`, which replaces `log1p`
  with a cubic series below `1.4e-4` for a measured 43% speed-up of the filter's scoring pass. They
  differ in the last bits; merging changes output.

### D2 — Progress and log-line formatting, 9 helpers in 5 files

- **Where**:
  - thousands separators: `src/cli/run_output.rs:480`, `src/neighbour_evidence/peak.rs:408` — identical;
  - durations: `src/cli/run_output.rs:495` (`elapsed_text`), `src/parameter_estimation/progress.rs:155`
    (`duration`, which adds a `350ms` case);
  - byte sizes, five spellings: `run_output.rs:369` (`3.2 GiB`/`850 MiB`), `progress.rs:138` (same
    rule), `progress.rs:150` (from kibibytes), `src/run/paralog_filter/finish.rs:408` (two decimals,
    down to `kiB` and bytes), `peak.rs:422` (**decimal** `GB`);
  - `/proc/self/status` read twice: `run_output.rs:355` (`VmHWM`), `progress.rs:93` (`VmRSS`, `VmHWM`).
- **Proposal**: put them in the layer-0 module of A2. Merge the two identical thousands helpers and
  the duration pair now. Leave the size spellings as they are unless the owner wants one style.
- **Downside**: unifying sizes changes terminal and log text (not output files). `peak.rs` is the
  one place that prints decimal gigabytes while every other message prints binary ones; whether that
  matches how `--max-memory` is read is a question for the owner, not a merge.
- **Effect**: none for the identical ones; log text only for the rest. **Cost**: S. **Confidence**: high.

### D3 — MD5 digest to hex, seven copies

- **Where**: `src/reference_info.rs:79` (`format_md5_hex`), `src/calling/parameters_file/bindings.rs:192`,
  `src/psp/header.rs:1149`, `src/repeat_catalog/parquet_file.rs:444`, `src/repeat_catalog/reader.rs:634`,
  `src/cli/calling_run.rs:1623` (inline), `src/vcf/header.rs:788` (inline).
- **What**: all write 32 lowercase hex characters.
- **Proposal**: one function in `types.rs` (every module can reach it) or `reference_info`. Leave the
  three decoders apart: `psp/header.rs:1153` refuses uppercase by the file format's rule,
  `read/input/open_bam.rs:856` accepts it because a BAM's `M5` may be uppercase, and
  `parquet_file.rs:448` reads back what this crate wrote.
- **Downside**: `bindings.rs:190` documents its copy as "the one place either binding's text is
  produced"; one shared encoder keeps that promise better.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### D4 — `ContigList` → `ContigBounds`, six copies with two policies

- **Where**: refuse a contig over 4 Gb: `src/cli/typed_regions.rs:512`, `src/cli/inspect_psp.rs:211`,
  `src/cli/run_ground.rs:224`, `src/cli/estimate_contamination.rs:602`. Clamp it to `u32::MAX`:
  `src/run/explain.rs:64`, `src/psp/segmentation_section.rs:417`.
- **What**: the same loop, each with its own error type. The clamp/refuse split is deliberate and
  documented at both clamping sites.
- **Proposal**: `ContigBounds::of_list(&ContigList) -> Result<Vec<_>, ContigTooLong>` in
  `regions.rs`; the four refusing callers map the error into their own. Leave the two clamping ones.
  Fixes B1 as a side effect.
- **Downside**: small; the CLI error types stay.
- **Effect**: none (B1's message aside). **Cost**: S. **Confidence**: high.

### D5 — Contig name → `ContigId`, three production copies and five in tests

- **Where**: `src/cli/estimate_parameters.rs:955` (`contig_id_of`), `src/cli/estimate_contamination.rs:772`,
  `src/run/gatherer.rs:310`; test copies in `src/run/census_fit.rs` (lines 1058-1540).
- **Proposal**: `ContigList::id_of(&self, name) -> Option<ContigId>`, with a map if B2 is taken up.
- **Downside**: none. **Effect**: none (speed only with the map). **Cost**: S. **Confidence**: high.

### D6 — Write to a temporary name, then rename: five copies, three behaviours

- **Where**: `src/vcf/writer.rs:494`/`:594` (`.tmp`; fsyncs the file and the directory);
  `src/cli/typed_regions.rs:627`/`:746` (`.tmp`; fsyncs the file); `src/repeat_catalog/parquet_file.rs:184`/`:173`
  (`.tmp`; **no fsync**); `src/cli/generate_psps.rs:782`/`:820` (`.psp.<pid>.partial`; the psp writer
  fsyncs at `src/psp/writer.rs:761`); `src/run/paralog_filter/resume.rs:232` (`.partial`; fsyncs).
- **What**: the same idea with different suffixes and durability. The catalog's missing fsync means a
  crash just after the rename can leave a truncated catalog under its final name; the reader would
  refuse it (no footer), so the consequence is an error, not wrong data.
- **Proposal**: leave the writers — their shapes differ (streaming writers, a closure, a Parquet
  writer). At most share one `temporary_name_beside(path, suffix)`; add an fsync to the catalog if
  the owner cares about that crash window.
- **Downside of merging**: couples five file formats' finish logic for a two-line saving.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### D7 — The global thread pool is set up in six commands

- **Where**: `src/cli/call_from_alignments.rs:676`, `call_from_psps.rs:568`, `estimate_parameters.rs:403`,
  `estimate_contamination.rs:380`, `resume_calling.rs:151`, `repeat_catalog.rs:205`.
- **What**: five copies skip the call when `--threads` is 0; `repeat_catalog.rs` takes 0 as 1 on
  purpose (its header records the pool's size).
- **Proposal**: leave it, or one `cli::set_global_pool(threads)` for the five alike. Low value.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### D8 — Small numeric and geometric helpers

| Copies | Where | Same result? | Proposal |
|---|---|---|---|
| `GenomeRegion` overlap | `src/cli/inspect_psp.rs:206`, `src/run/edge_deletion_gathering.rs:101` | yes | add `GenomeRegion::overlaps` beside `contains` (`types.rs:116`) |
| `logit`, `sigmoid` | `src/paralog/prior.rs:304`, `:315`; `src/paralog/windowed_prior.rs:354`, `:358` | yes, byte-identical | keep one pair in `paralog` |
| Phred → `ln P(error)` | `src/read/prepared_read.rs:150`, `src/locus_generation/pileup/open_record.rs:3645` (table) | yes, same expression | share the table version |
| the same, third copy | `src/locus_generation/ssr.rs:1869` | **no** — multiplies by a pre-divided `LN_10 / 10`, which can round differently in the last bit | **leave**; merging changes output |
| binomial upper tail | `src/calling/quality/artifact_correction.rs:283` (incomplete beta), `src/parameter_estimation/joint/ssr_fit/tract_exclusion.rs:270` (summed series) | **no** — different algorithms | **leave**; merging moves which tracts the fit excludes |
| `logistic` | `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:2939` vs the `sigmoid`s above | **no** — one branch only | **leave** |
| locus-selection hash | `src/repeat_catalog/strata.rs:164` (`hash_locus`), `src/parameter_estimation/joint/loci.rs:276` (`ContigHasher`) | yes — documented as having to agree | move `ContigHasher` into `repeat_catalog` and call it from both, so the agreement is structural |

All "yes" rows: no output change, cost S, confidence high.

---

## Dead code

### X1 — `bam::segment_reader`, `segment_merge`, `bam_input`, `cram_input`: 2,633 lines with no production caller

- **Where**: `src/bam/segment_reader.rs` (1,896 lines; `allow(dead_code)` at `:59`),
  `src/bam/segment_merge.rs` (551; `:28`), `src/bam/bam_input.rs` (84) and `src/bam/cram_input.rs`
  (102), whose only user is `segment_reader.rs:81-82`. Declared at `src/bam/mod.rs:36-57`.
- **What**: outside `bam/`, only doc comments name them. Their comments say they wait for
  `run_pileup` and `AlignmentMergedReader` to be "flipped onto" them; both belonged to the caller
  deleted on 2026-09-13. The live read path is `read/input/` (`open_bam.rs`, `aligned_reads_reader/`).
- **Proposal**: delete the four files, then whatever in `bam/alignment_input.rs` and `bam/errors.rs`
  only they used (for example `ContigInterval::overlaps_record`, and the `AlignmentInputError`
  variants its doc comments at `errors.rs:203-211` tie to `segment_reader`).
- **Downside**: their tests exercise CRAM decoding; `read/input/open_bam.rs` and
  `read/input/test_fixtures.rs` (`indexed_cram`, `multi_container_cram`) already cover CRAM on the
  live path, so check no case is covered only here before deleting.
- **Effect**: none. **Cost**: S–M. **Confidence**: high.

### X2 — `iter_ext` is an empty module

- **Where**: `src/iter_ext/mod.rs` (5 lines of doc naming `BufferedPeekable`, which exists nowhere),
  `src/lib.rs:91`, mentioned at `src/lib.rs:34`.
- **Proposal**: delete. **Effect**: none. **Cost**: S. **Confidence**: high.

### X3 — Module-wide dead-code silencers over live code

- **Where**: `src/calling/allele_candidates/ssr.rs:35` (`#![cfg_attr(not(test), expect(dead_code))]`),
  `src/read/input/aligned_reads_reader/mod.rs:102-106` (`#![allow(dead_code, …)]`).
- **What**: `ssr.rs:30-34` says "a non-test build expects everything here to be dead", but
  `select_ssr` is called from `src/run/callers.rs:1580` on every run. What is still test-only is
  `table_indices_at` (`ssr.rs:170`) and `rescue_occupied_neighbours` (`ssr.rs:686`). In
  `aligned_reads_reader/`, the allow covers the live BAM, CRAM, container and bgzf readers to exempt
  the in-memory reader's constructors.
- **Proposal**: replace each with item-level `#[cfg_attr(not(test), expect(dead_code, reason = …))]`
  on the few test-only items — the idiom the crate already uses at about ten other sites — and correct the
  `ssr.rs` comment.
- **Downside**: none.
- **Effect**: none. **Cost**: S. **Confidence**: high.

### X4 — Unused or test-only items in library code

No non-test user anywhere (confirmed by grep over `src/`, `examples/`, `benches/`, `tests/`):

- `src/tandem_repeat.rs:342-343` `bases_match`, under a bare `#[allow(dead_code)]`;
- `src/bam/alignment_input.rs:431-432` `FilterBucket::BaqRejected` and the `baq_rejected` counter
  (`:174`) — never incremented, and BAQ-processed input is unsupported;
- `src/parameter_estimation/joint/stratum_fits.rs:1073` `each_stutter_multiplier`.

Test-only, superseded:

- `src/parameter_estimation/joint/ssr_fit/read_group_stutter.rs:94` `fit_multiplier` and `:235`
  `scale_to_a_typical_read_group` — the fit calls `fit_multipliers_of_a_typical_one` (`:150`); the
  module doc (`:17`) still names the old one as what fixes the scale.

Test-only `pub` items with no `cfg` (a library crate's `pub` items never trigger the lint):
`calling/likelihood/generic.rs:340` `charges_of`, `calling/likelihood/mod.rs:1152` `with_frequencies`,
`neighbour_evidence/refine.rs:378` `rounds_for`, `neighbour_evidence/stage.rs:332` `VectorStageStop`,
`neighbour_evidence/templates.rs:127` `choose_templates`, `read_phase/fragments.rs:304`/`:313`
`is_compact`/`is_tight`, `tandem_repeat.rs:102` `DEFAULT_PERIODS`, `tandem_repeat.rs:1008` `over_slice`,
`cli/run_output.rs:291` `PARAMETERS_BLOCK_BEGINS`.

- **Proposal**: delete the first three groups; mark the last group `#[cfg(test)]` in one sweep.
- **Downside**: none. **Effect**: none. **Cost**: S. **Confidence**: high.

---

## Test-support duplication

### T1 — About 640 lines of byte-identical test helpers

Measured by normalising whitespace and comparing whole helper bodies in test code across files.

- **`GenomeRegion` constructors**: a `region(…)`/`region_on(…)` helper is defined in **24 files**
  (168 lines); 11 of them are the same `region(start, end)` on contig 0 (for example
  `src/read/input/cursor.rs:1006`, `src/locus_generation/mod.rs:1377`, `src/vcf/tests.rs:11`).
  `run/cohort_merge/mod.rs:61` already has a shared `fixtures` module, yet `build.rs:2556` and
  `close.rs:1106` define `region_on` again.
  **Proposal**: a `const fn GenomeRegion::on(contig, start, end)` in `types.rs`. Library code
  writes the struct literal 219 times outside `tests.rs` files too, so it is not test-only sugar.
- **`run/paralog_filter/` tests**: `scratch`, `quality`, `allele`, `outbred`, `no_call`, `pools`,
  `a_window`, `no_window`, `a_fittable_histogram` repeated across `pass_one/`, `pass_two/`,
  `pass_three/`, `patch/`, `finish/`, `scoring_context/`, `spill_file/` tests (about 100 lines).
  **Proposal**: a `paralog_filter/test_fixtures.rs`, as `run/test_fixtures.rs` did for `run`.
- **Smaller pairs**: `make_read` (38 lines, `locus_generation/pileup/cigar_cursor.rs:963` and
  `decompose.rs:213`), `one_batch` (`calling/mod.rs:4771`, `calling/run_parameters.rs:898`),
  `a_cohorts_read_groups` (`cli/call_from_alignments/tests.rs:475`, `run/report/tests.rs:136`),
  the `ReleaseSpy` reference in `locus_generation/pileup/generator.rs:3000-3026` and
  `locus_generation/ssr.rs:5291-5321`.
- **Downside**: a test that defines its own `region` shows its contig at a glance; a shared helper
  must be named well enough to do the same (`GenomeRegion::on(0, 10, 20)` does).
- **Effect**: none. **Cost**: M for all; S for `GenomeRegion::on` alone. **Confidence**: high.

---

## Readability

### R1 — `lib.rs`'s module list is out of date

- **Where**: `src/lib.rs:7-35`.
- **What**: names `iter_ext` (empty, X2) as a shared helper; omits `float`, `splitmix`,
  `open_files` and `memory_limit`.
- **Proposal**: update when X2 lands. **Effect**: none. **Cost**: S. **Confidence**: high.

---

## Looked at, left alone

- The three `test_fixtures.rs` (`read/input/` 801 lines, `run/` 248, `cli/` 622) build different
  things at different layers — files on disk, catalogs and segmentations, whole cohorts. Not copies.
- `run/cohort_merge/rounds.rs`: `merge_and_call_in_rounds` has no production caller by design — its
  doc keeps it as the stream driver's test oracle; `Frontier` and `CalledRegion` in it are live.
- Two k-way argmin merges (`read/input/sample_cursor.rs:232`, `run/cohort_merge/close.rs:586`) merge
  different things; the third (`bam/segment_merge.rs:190`) goes with X1.
- BED parsing exists once (`regions.rs:113-202`); `region_typing` and `run::explain` wrap it.
- `#[cfg_attr(not(test), expect(dead_code, reason = …))]` on single items (about ten sites) is the crate's
  deliberate idiom for test-read instruments; fine.
- `genetics::lgamma` versus the Lanczos `ln_gamma` in `parameter_estimation/joint/ssr_fit.rs:1079`:
  the second is `#[cfg(test)]` for a simulator.
- `parquet_file.rs:448` `unhex` uses `u8::from_str_radix`, which accepts a leading `+` (`"+f"` reads
  as 15). It reads back digests this crate wrote; not worth a change.
- `paralog::locus_score::log_add_exp` versus `log_sum_exp_2`: deliberately different (D1).
- `psp` importing `locus_generation`, `region_typing` and `segmentation_inputs`: a store of loci
  depending on the loci's types is the right direction.
