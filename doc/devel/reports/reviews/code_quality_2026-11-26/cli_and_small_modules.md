# Code-quality map: `src/cli/`, `src/main.rs` and the smaller library modules (2026-11-26)

Scope: `src/cli/`, `src/main.rs`, `src/neighbour_evidence/`, `src/region_typing/`,
`src/repeat_catalog/`, `src/vcf/`, `src/window_coverage/`, and the single-file modules at the crate
root (`types.rs`, `reference_info.rs`, `ref_seq.rs`, `genetics.rs`, `regions.rs`,
`tandem_repeat.rs`, `production_post_filter.rs`, `scanner_parity.rs`, `golden_catalog.rs`,
`segmentation_inputs.rs`, `open_files.rs`, `memory_limit*`, `float*`, `splitmix.rs`,
`error_render.rs`, `raw_chrom_reader.rs`, `iter_ext/`, `fasta/`). About 80,000 lines, tests
included.

The review only read code; nothing was built or run. Two helper agents read
`neighbour_evidence/` and the `region_typing/` + `repeat_catalog/` + `window_coverage/` group. Their
findings are merged below. I re-read the citations behind every finding marked **verified**.

## Summary, in priority order

1. **No bug found that changes a genotype, QUAL or the record set at default settings.** The VCF
   writer, the record assembly and the neighbour-evidence pass were traced at the usual edges
   (multi-allelic sites, a deletion at a contig's first base, all-no-call samples, zero depth,
   NaN in a likelihood row). They hold.
2. **One suspected output bug at non-default settings.** With `--max-str-len` raised above
   `--max-str-bundle-len`, a run given `--regions` can classify a long satellite as an STR
   segment, where the whole-genome run classifies it as a satellite. The cause is that the
   catalog reader's padding around a requested region ignores the satellite length cap.
3. **The two calling commands declare about 25 identical flags twice** (about 210 lines each in
   `call_from_alignments.rs` and `call_from_psps.rs`). They also copy the last ~150 lines of their
   run function: the scratch directory, passes two and three, the explain table and the report.
   Every new calling flag has to be added in two places, plus about 20 hand-written test argument
   structs.
4. **Library code imports from `cli` in two places, both for log formatting**: the `PROGRESS`
   log target, `with_thousands` and `elapsed_text`. A third module, `neighbour_evidence/peak.rs`,
   copies `with_thousands` byte for byte to avoid that import. Moving these three into a small
   library module removes both wrong-way imports and the copy.
5. **Dead code**: `src/iter_ext/` is an empty module whose doc names a type that no longer exists.
   About 320 lines of `tandem_repeat.rs` (a third region tiler, `RegionScanner` and its helpers)
   have no caller outside their own tests.
6. **`--max-depth` has two default constants of 1,000**, one per command. The psp route matches
   the direct route only while they are equal.
7. **The test-only "production" modules are live tests, not leftovers.** These are
   `golden_catalog.rs`, `production_post_filter.rs`, `scanner_parity.rs` and
   `window_coverage/production_parity.rs`. Keep them. `types.rs` is best left as one file; only its
   module doc should change.

---

## Bugs

### B1. A region read can classify a satellite differently from a whole-genome read when `--max-str-len` exceeds `--max-str-bundle-len` — *suspected*
- **Where:** `src/repeat_catalog/mod.rs:459-479` (`reach_into_window`), called from
  `src/repeat_catalog/reader.rs:213-217` (`rows_for`). Verified.
- **What:**
  - When `--regions` is given, the catalog reader reads only the rows within a padding (`reach`)
    around each requested span. The padding is `max(longest tract on the contig,
    max_str_bundle_len) + bundle_threshold`.
  - A satellite is a merged run of overlapping or abutting tracts longer than `max_str_len`. That
    run can be longer than every single tract in it.
  - The function's own doc (`mod.rs:456-458`) admits that a satellite reaching past the padding
    "comes back cut". When `max_str_len > max_str_bundle_len`, being cut can drop it below the
    satellite cap.
  - Failure scenario:
    - Settings: `--max-str-len 400`, bundle cap 100 (the default), and a contig whose longest
      tract is 120 bp.
    - A run of abutting tracts spans 419 bp.
    - The whole-genome read sees all 419 bp and classifies a satellite.
    - A BED span at the run's last tract pads by 120 + radius. That read sees about 220 bp of the
      run, so the run's last locus is classified as an STR segment and called with the repeat
      model.
  - The defaults set both caps to 100 (`region_typing/mod.rs:425,431`), which is why no test sees
    it.
- **Proposal:** include `max_str_len_bp` in the `max`, and pass it in from `rows_for`. Add a
  differential case that sets the two caps differently.
- **Downside:** larger region reads when the user raises `--max-str-len`, which costs little.
- **Effect:** output changes, at non-default settings only; speed and memory are negligible.
  **Cost** S. **Confidence** medium. The scenario was reasoned through but not run.

### B2. A malformed catalog footer can panic instead of returning `Unreadable` — *likely*
- **Where:** `src/repeat_catalog/parquet_file.rs:448-457` (`unhex`). Verified.
- **What:**
  - `unhex` checks `text.len() == 32`, which counts bytes, then slices `&text[i*2..i*2+2]`.
  - A 32-byte footer string that holds a multi-byte UTF-8 character makes the slice cut through
    that character, and the slice panics.
  - The input is text read from a file. That turns a corrupt-file error into a crash.
- **Proposal:** return `None` unless `text.is_ascii()`. At the same time, replace the two private
  `hex` copies (`parquet_file.rs:444`, `reader.rs:634`) with `reference_info::format_md5_hex`
  (`reference_info.rs:79`).
- **Downside:** none. **Effect:** none on valid files. **Cost** S. **Confidence** high.

### B3. `GQ` and `GP` can disagree on a record that dropped an uncalled alternative — *suspected, low consequence*
- **Where:** `src/vcf/assemble.rs:480-584` (`drop_alternatives_no_sample_calls`) and
  `:586-628` (`rows_over_kept_alleles`).
- **What:**
  - When a record drops an alternative that no sample is called with, each sample's `GP` row is
    cut to the genotypes over the kept alleles and renormalised.
  - The sample's `GQ` keeps the value the loop computed over the full genotype table.
  - So the written `GQ` is lower than `GP` implies; it is the more cautious of the two.
  - The renormalisation was a deliberate choice; leaving `GQ` alone is not documented. On GIAB at
    30×, 9 records in 2,290 dropped an allele (`assemble.rs:448-450`), so few lines are affected.
- **Proposal:** decide which is intended and document it. If `GQ` should match `GP`, recompute it
  from the renormalised posterior.
- **Downside:** recomputing changes `GQ` bytes on those records, so golden files and digests need
  regenerating.
- **Effect:** output changes if fixed. **Cost** S. **Confidence** medium.

### B4. The log's "records re-called" count includes lines where no sample changed, and those lines gain an empty `PS` — *likely, log and cosmetic only*
Found by the neighbour-evidence helper; I did not re-trace it.
- **Where:**
  - `src/neighbour_evidence/line.rs:123-135` (`site_of`).
  - `stage.rs:797-819` and `:1180-1186`.
  - `line.rs:192` (the count) and `line.rs:176`, `:437-440` (`:PS` appended with `.`).
  - The count is logged at `cli/calling_run.rs:634-637`.
- **What:**
  - A record where only one sample has `PL`/`GP`, or where every sample's outcome is "unchanged",
    still goes through `rewrite_line`.
  - The log counts it as re-called.
  - Its line gains `:PS` with `.` in every sample column.
  - Genotypes are unaffected, and `.` is a valid `PS` value.
- **Proposal:** count only lines where some sample's outcome changed. Optionally, pass such lines
  through unchanged.
- **Downside:** the optional part changes output bytes (the `PS` key disappears on those lines).
- **Effect:** the count fix changes the log only. **Cost** S. **Confidence** medium-high.

### B5. `--threads 0` means one thread for `repeat-catalog` and every core for every other command — *likely, minor*
- **Where:** `src/cli/repeat_catalog.rs:52-60` (`default_value_t = 1`, and the doc does not
  mention 0) and `:204-207` (`args.threads.max(1)`). Compare `call_from_alignments.rs:178-180`,
  `resume_calling.rs:51-53` and `estimate_parameters.rs:134-136`: "0 uses every core".
- **What:** a person who types `--threads 0` on every command gets a single-threaded catalog
  build. The output is the same, and the command takes about 3.3 times as long as at six threads on
  human (its own doc gives 40 s against 12 s).
- **Proposal:** either document "0 means 1" on this flag, or make 0 mean every core here too. The
  per-thread memory (1.4 GB) is the reason the default is 1. Keep the default at 1 and change only
  what 0 means.
- **Downside:** with 0 meaning every core, a 64-core machine would hold 64 contigs at once. The
  documentation fix avoids that risk.
- **Effect:** none on output. **Cost** S. **Confidence** high.

### B6. A command-line argument containing a newline breaks the VCF header — *suspected, rare*
- **Where:** `src/cli/provenance.rs:13-18` (`current_command_line` joins the arguments verbatim)
  and `src/vcf/header.rs:683-692`, which writes `##commandline={value}` and
  `##reference={path}` without escaping.
- **What:**
  - A path argument with an embedded `\n` would split the header line in two.
  - The second half becomes a header line that does not start with `##`, and the file is invalid
    VCF.
  - `VcfHeaderMetadata::try_new` checks for empty and duplicate sample and contig names. It does
    not check for control characters.
- **Proposal:** replace `\n`, `\r` and `\t` with spaces (or `\n` escapes) in `current_command_line`
  and the reference path.
- **Downside:** none. **Effect:** none on ordinary input. **Cost** S. **Confidence** medium.

---

## Architecture

### A1. Library code imports log plumbing from `cli` — move the three helpers into the library
- **Where:**
  - `src/parameter_estimation/progress.rs:25` imports `crate::cli::run_output::PROGRESS`.
  - `src/run/calling_progress.rs:23` imports `crate::cli::run_output::{PROGRESS, elapsed_text,
    with_thousands}`.
  - The definitions are at `src/cli/run_output.rs:44`, `:480` and `:495`.
  - `src/neighbour_evidence/peak.rs:408-418` holds a byte-for-byte copy of `with_thousands`, made
    so that it would not import from `cli`.
- **What:**
  - These are the only non-test imports of `crate::cli` from outside `src/cli/`. I checked with
    grep.
  - Every other mention is in a doc comment or a `#[cfg(test)]` block. Tests in
    `calling/parameters_file/from_run_parameters.rs:3281` and `run/census_*.rs` drive the
    subcommands, which is fine for tests.
- **Proposal:**
  - Create a small crate-root module, for example `src/progress_text.rs`, holding `PROGRESS`,
    `with_thousands`, `elapsed_text`, and possibly `size_of_ground` and `peak.rs`'s `gb`.
  - `cli::run_output` re-exports them, so its callers do not change.
  - Delete `peak.rs`'s copy.
- **Downside:** one more tiny module, and the target name would sit apart from the subscriber that
  interprets it. A doc link between the two keeps that readable.
- **Effect:** none. **Cost** S. **Confidence** high.

### A2. The round-size policy lives in `cli` while `run`'s docs describe it — leave it unless the code is touched
- **Where:** `src/cli/calling_run.rs:383-560` (`round_width_for`, `round_shape_for`,
  `narrowest_round_region`, `REGIONS_A_WORKER`). It is referred to from `run/callers.rs:180`,
  `run/cohort_merge/stream.rs:30` and `run/cohort_merge/rounds.rs:27`.
- **What:**
  - This is a memory and parallelism policy for the cohort merge, and it reads `cohort_merge`'s own
    constants (`stream::ROUNDS_HELD`).
  - Only the command layer calls it, and no library code imports it, so there is no wrong-way
    import. Only the doc links point the wrong way.
- **Proposal:** move it to `run::cohort_merge` if that file is being edited for another reason;
  otherwise leave it.
- **Downside:** the module doc of `calling_run.rs` (`:1-23`) argues the opposite position: the
  round width is one decision both calling commands share, and that is what the module is for. The
  argument is reasonable. Leave it.
- **Effect:** none. **Cost** S. **Confidence** medium.

### A3. `window_coverage::depth` and `run::cohort_merge::observation_cache` import each other
Found by the helper agent; I verified the two import lines.
- **Where:**
  - `src/window_coverage/depth.rs:54` imports `run::cohort_merge::observation_cache::{Drawn,
    LocusSummary, StandIn}`.
  - Its only caller is `src/run/cohort_merge/observation_cache.rs:52`.
  - The module doc of `window_coverage/mod.rs` (`:28-32`) acknowledges the cycle.
- **What:**
  - `depth.rs` is 772 lines, about 560 of them tests. It walks the cache's own types.
  - The rest of `window_coverage` depends only on `types` and `float`.
  - It is `pub` for an example (`ng_window_coverage_probe.rs`) that no longer exists.
- **Proposal:** move it under `run/cohort_merge/` (for example `observation_cache/reported_depth.rs`)
  and make it `pub(crate)`.
- **Downside:** the per-record depth rule would no longer sit beside the accumulator it feeds, and
  a reader of spec `window_coverage.md` §3.1 has to look in `run/`.
- **Effect:** none. **Cost** S-M. **Confidence** high.

### A4. `neighbour_evidence` and `vcf` import each other — leave it
- `vcf/writer.rs:36-39` drives the neighbour-evidence stage. `neighbour_evidence/line.rs:21-25` and
  `refine.rs:50` use `vcf::encode` so that the `GQ` the pass tests against is the `GQ` that will be
  written.
- The coupling is real: the stage sits in front of the writer. Rust does not care about cycles
  inside one crate. Breaking this one needs a trait or a file move that makes the code harder to
  follow. Leave it.

### A5. `types.rs` (2,564 lines): leave it as one file and fix its module doc
- **Where:** `src/types.rs:1-4`. Code runs to `:1478`; the rest is tests.
- **What:**
  - The file holds about 30 small newtypes in five groups:
    - coordinates: `ContigId`, `Position`, `GenomePosition`, `GenomeRegion`, `Bp`;
    - qualities: `MapQual`, `BaseQual`, `Phred`, `LogProb`, `SummedLogError`, `MismatchFraction`,
      `ErrorRate`;
    - read groups and batches: `ReadGroupId`, `ChainId`, `BatchId`;
    - population parameters: `InbreedingF`, `ExpectedHeterozygosity`,
      `ExpectedAlternativeFrequency`, `Ploidy`;
    - alleles: `AlleleId`, `Motif`, `SsrPeriod`, `Genotype`, `GenotypeRows`.
  - `DomainError` is shared by all five groups.
  - Each type is short and self-contained, and every one has users elsewhere (checked by grep),
    so there is no dead type.
  - The module doc is stale. It says the file is "seeded here with only what the `RefSeq`
    reference accessor needs" and that it will split by the stale `module_layout.md`.
- **Proposal:** rewrite the module doc to list the five groups in order, and leave the file whole.
  - Splitting it into `types/{coordinates,quality,…}.rs` with re-exports is cheap and
    behaviour-neutral. It would mostly move the 1,080 lines of tests around.
  - The shared `DomainError` would have to live in a sixth file that every other file imports.
  - The file reads easily today as a dictionary of small types.
- **Downside of leaving it:** the file keeps growing. Revisit when one group passes about 600 lines
  of code.
- **Effect:** none. **Cost** S (doc). **Confidence** medium.

### A6. `cli` modules that are test-only or empty in a release build
- `pub mod mode_equivalence;` (`src/cli/mod.rs:19`) contains only a doc and a `#[cfg(test)] mod
  tests` (`mode_equivalence.rs:75`). In a release build it is an empty public module.
- `cross_platform_digests` (`mod.rs:14`) has the same shape.
- **Proposal:** declare both with `#[cfg(test)] mod`. **Downside:** none. **Effect:** none.
  **Cost** S.

---

## Duplication

### D1. The two calling commands declare about 25 identical flags twice
- **Where:**
  - `src/cli/call_from_alignments.rs:155-365` and `src/cli/call_from_psps.rs:135-345`: `ploidy`,
    `min_site_quality`, `low_base_quality_weight`, `threads`, `force`, `log`, `tmp_dir`,
    `keep_tmp_on_error`, `max_generic_alleles`, `max_locus_span`, the three `hidden_dup_*`, the
    seven `refine_*`, and `phase_min_path_agreement`.
  - Also `locus_builder_regions_len` and `explain_loci`, at `:479/:490` and `:362/:384`.
  - A `diff` of the two blocks shows only the `--max-depth` doc paragraph differing, by one
    sentence. Everything else is the same text, defaults and parsers.
  - `resume_calling.rs` repeats 8 of the same flags.
- **What:**
  - About 210 duplicated lines.
  - Each new calling flag is added twice, and then to the hand-written argument structs in tests:
    7 in `call_from_alignments/tests.rs`, 7 in `call_from_psps/tests.rs`, 10 in
    `mode_equivalence.rs`, 3 in `cross_platform_digests.rs`, 5 in `resume_calling/tests.rs`, and
    one each in `run/census_*.rs`.
  - The two `--max-depth` docs disagree on purpose. The two `--refine-*` docs have no reason to.
- **Proposal:**
  - Move the identical flags into one `#[derive(clap::Args)] pub struct CallingFlags` in
    `calling_run.rs`, and `#[command(flatten)]` it into both commands.
  - Keep `--max-depth` and the input flags (`--alignment`/`--psp`, the repeat criteria) in each
    command.
  - For tests, build arguments with `Cli::parse_from([...])`, as `call_from_alignments/tests.rs:24`
    (`args_of`) already does, instead of literal structs. They then follow clap's defaults.
- **Downside:**
  - **The `--help` order may change.** Flattened arguments are listed where the flatten field sits,
    and today the shared flags are interleaved with command-specific ones.
  - Help text is part of the command line that must not change, so the change needs a before/after
    diff of `--help` for both commands. It may need `display_order`, or the command-specific fields
    reordered.
  - The `parameters`/`defaults` `ArgGroup` must stay on the outer struct.
- **Effect:** none on calls. `--help` order is at risk. **Cost** M. **Confidence** high.

### D2. The last ~150 lines of the two calling run functions are copied
- **Where:**
  - `src/cli/call_from_alignments.rs:892-1050`: the scratch directory, the sink, passes two and
    three of the duplicate filter, the explain table, the report and the "done" line.
  - The same sequence is in `src/cli/call_from_psps.rs:786-949`.
  - After renaming the two error enums, `diff` shows about 130 differing lines. Most are layout.
    The real differences are:
    - alignments hands each record over (`call_cohort_handing_each_record_over`), where psps
      prepares spill entries on the pool (`call_cohort_preparing_each_record`);
    - `RunReport::of` against `RunReport::of_a_stored_cohort`;
    - `written` against `(calling, stored)`.
  - The `RefineWindowArgs { … }` literal is also built twice (`call_from_alignments.rs:866`,
    `call_from_psps.rs:760`).
- **What:** `calling_run.rs` exists to stop the two modes drifting (its doc, `:5-10`), but the
  filter-and-report tail still drifts by hand.
- **Proposal:**
  - Add a `calling_run::finish_the_run(...)`, generic over the command's error type through
    `From`.
  - It takes the sink outcome, the coverage histograms, the explanations and the report, and runs
    everything from "take the sink apart" to `the_run_succeeded`.
  - Build `RefineWindowArgs` from the flattened `CallingFlags` (D1).
- **Downside:**
  - The generic error plumbing adds one indirection.
  - The two calling APIs differ in their record callbacks, so the shared part starts only after
    calling returns.
  - Do it after D1, not before.
- **Effect:** none on output. The alignments path's per-record encoding could be checked for
  performance at the same time; the psp path encodes on the pool. **Cost** M.
  **Confidence** medium.

### D3. Two default constants for `--max-depth`, both 1,000
- **Where:** `src/run/depth_ceiling.rs:42` (`DEFAULT_MAX_READ_GROUP_DEPTH`, the calling commands'
  default) and `src/locus_generation/pileup/mod.rs:354` (`DEFAULT_MAX_READS_PER_POSITION`, the
  `generate-psps` default at `cli/generate_psps.rs:154`).
- **What:**
  - The psp route reproduces the direct route only if a psp's stored read cap equals the calling
    run's ceiling. The flag docs say: "At equal values a thinned position is always dropped by
    `--max-depth`".
  - Both doc comments record the same owner decision, dated 2026-09-29.
  - Nothing ties the two numbers together.
- **Proposal:** make one `pub const` the other's value, or add a `const _: () = assert!(…)`.
- **Downside:** none. **Effect:** none. **Cost** S. **Confidence** high.

### D4. `psps_named_by` is written twice, differing only in the error enum
- **Where:** `src/cli/call_from_psps.rs:540-553` and `src/cli/estimate_parameters.rs:1047-1064`.
- **Proposal:** a generic `psp_inputs::psps_for_a_run<E: From<PspArgumentRefusal> +
  From<OpenFilesError>>`, or leave it.
- **Downside:** it is 14 lines, and each copy maps the refusals into its command's own wording.
  **Lean: leave it** unless D1/D2 already touches it. **Cost** S.

### D5. Smaller copies, each worth folding in only when the file is already open
- `FIRST_POSITION_OF_A_CONTIG` is defined twice, at `src/vcf/mod.rs:66` (`pub(crate)`) and
  `src/run/records.rs:65`. **Proposal:** make `run/records.rs` use the `vcf` constant. S, no
  downside.
- The rayon pool setup (`ThreadPoolBuilder::new().num_threads(..).build_global()`) appears 6
  times in `src/cli/`: `call_from_alignments.rs:676`, `call_from_psps.rs:568`,
  `estimate_contamination.rs:380`, `estimate_parameters.rs:403`, `resume_calling.rs:151` and
  `repeat_catalog.rs:205`. A helper would also settle B5. S.
- **Leave it:** `contig_bounds` (`inspect_psp.rs:211`, `typed_regions.rs:512`) and `tmp_path_for`
  (`typed_regions.rs:627`, `vcf/writer.rs:494`). Each is under 20 lines and carries its own
  error type.
- **Leave it:** the repeat-criteria flags (`min_copies` … `min_purity`) in `call_from_alignments`
  and `generate_psps`. Their help text differs on purpose ("called as" against "stored as"), and
  clap cannot vary the doc text of a flattened field.
- Within `region_typing`, `reference_runs.rs:190-198` (`is_a_motif_of_its_own_period`)
  re-implements `segment_criteria.rs:1335-1340` (`is_compound`) plus an A/C/G/T check.
  **Proposal:** `all ACGT && !is_compound(unit)`. S, no downside. (Helper finding.)
- The catalog route restates the rule for which region kinds are clipped at a BED edge
  (`repeat_catalog/segments.rs:101-104`, a `matches!`), against the exhaustive match at
  `region_typing/mod.rs:693-698` (`clips_at_a_bed_edge`). A new `RegionKind` would be checked by
  the compiler only on the test-only scan route. **Proposal:** call `clips_at_a_bed_edge`. S.
  (Helper finding.)
- `neighbour_evidence/model.rs:371-384` and `:469-482` contain the same 14-line "weigh the shares
  by the likelihood" loop. One function would guarantee the bit-equality that the two copies keep
  by hand. It is the pass's hottest loop, so a timing run is needed. S. (Helper finding.)
- Seeded test random generators (`struct …(u64)` with the SplitMix constant) appear about 12
  times in test modules. A `#[cfg(test)]` crate helper would remove them. Low value, so leave it.
  `locus_generation/pileup/read_sampling.rs:91` is a *different* function from `splitmix.rs` (no
  initial add), so leave that one too.

---

## Dead code and migration leftovers

### X1. `src/iter_ext/` is empty
- **Where:** `src/iter_ext/mod.rs:1-5`, which is five lines of doc and no code. It is declared at
  `src/lib.rs:91` and described at `lib.rs:34` as a "small helper the rest shares".
- **What:** its doc says it "hosts `BufferedPeekable`", which no longer exists anywhere in the
  crate. Grep finds no use of `iter_ext` outside `lib.rs`.
- **Proposal:** delete the module and its mention at `lib.rs:34`.
- **Downside:** none. **Effect:** none. **Cost** S. **Confidence** high.

### X2. A third region tiler in `tandem_repeat.rs` has no caller
- **Where:** `src/tandem_repeat.rs:876-1027`: `union_spans`, `build_regions`, `tile` and
  `RegionScanner` at `:1002`. Its types are `SegmentOptions` (`:165`), `RegionSpan` (`:227`),
  `RepeatRegion` (`:239`), `Region` (`:256`) and `ScanError` (`:273`). There are also the
  `DEFAULT_*` constants it alone uses, and `bases_match` (`:340`, `#[allow(dead_code)]`).
- **What:**
  - Grep over `src/`, `examples/`, `tests/` and `benches/` finds no user outside the file's own
    tests.
  - `RegionSpan` appears only in a doc link (`alignment/mod.rs:132,143`).
  - `region_typing/mod.rs:402` still says "`RegionScanner` still does" scan in pieces, which is
    now stale.
  - About 150 lines of code and about 170 of tests.
- **Proposal:** delete it, and fix the two doc references.
- **Downside:** it is `pub` API. No binary, example or test uses it, but the owner may want to
  confirm.
- **Effect:** none. **Cost** S. **Confidence** high.

### X3. `lowDepth` is declared in every header and never written
- **Where:** `src/vcf/header.rs:605` and `src/vcf/mod.rs:723,735`. `run/records.rs:806-810` says
  it is "declared and never written … a file written by an older caller can carry it".
- **What:** the header describes this file, and this caller cannot write the value. The reason
  given (an older caller's files) does not apply to the file being written.
- **Proposal:** owner's decision. Removing the line and the variant changes every header by one
  line.
- **Effect:** **output changes** (one header line). **Cost** S. **Confidence** high.

### X4. The test-only migration modules are live tests — keep them
- `golden_catalog.rs`, `production_post_filter.rs` and `scanner_parity.rs` are all `#[cfg(test)]`
  in `lib.rs:74-81`.
  - `production_post_filter` is used by the port-fidelity differential in
    `region_typing/segment_criteria.rs:2142-2143`.
  - `golden_catalog` is used by `region_typing/mod.rs:3106` and `repeat_catalog/anchor.rs:259`.
- `window_coverage/production_parity.rs` is `#[cfg(test)]` (`window_coverage/mod.rs:43-44`). It
  now compares 200 generated streams against recorded digests, and it is the accumulator's only
  random-stream test.
- None of them ships, and each still catches regressions. Optionally rename
  `production_parity.rs` (for example `recorded_windows.rs`) so the name stops pointing at a caller
  that no longer exists.

---

## Readability

### R1. Stale doc comments that state the wrong behaviour
Each is S, has no output effect, and has high confidence.
- `src/vcf/encode.rs:506-507`, on `genotype_field`, says "Unphased always. ng computes no
  phasing". Phase is written now, by read-backed phasing and by the neighbour-evidence pass
  (`header.rs:614-617`). It is true that the *encoder* writes `/`. Say that a later stage rewrites
  it.
- `src/vcf/assemble.rs:466-472` says that when every sample is homozygous-reference "the record is
  left exactly as it was". The early return at `:503` does not fire there (allele 0 is named), so the code drops every alternative in that case. The case
  never reaches it, because `run/records.rs:79` (`a_written_genotype_carries_an_alternative`)
  removes such loci earlier. Rewrite the paragraph to say the branch is defensive.
- `src/neighbour_evidence/stage.rs:15-17` says the pass is "on by default". The flag defaults to
  `Off` (`cli/call_from_psps.rs:263`). The ≥10-sample, median-inbreeding ≥ 0.5 rule is what
  `--refine-genotypes auto` applies.
- `src/types.rs:1-4`: see A5.
- `window_coverage/mod.rs:17-21` and `region_typing/mod.rs:56-81` still give "not edited" or
  "frozen" as reasons. The freeze ended with production's deletion. Rewrite them as history so a
  reader does not avoid an allowed change. (Helper finding.)

### R2. Five functions of 200 to 400 lines in `neighbour_evidence`
Helper finding.
- The functions:
  - `model.rs:1230-1645` (`walk_with_draw`);
  - `refine.rs:1374-1690` (`solve_round`);
  - `line.rs:155-454` (`rewrite_line`);
  - `stage.rs:981-1234` (`release`);
  - `refine.rs:1695-1903` (`final_solution`).
- **Proposal:** extract the per-site posterior in `walk_with_draw` (`model.rs:1456-1581`) first.
  It is the main obstacle to reading the walk.
- **Downside:** these are hot loops, and the extracted functions take long scratch-buffer argument
  lists. Each extraction needs a timing check; the bit-equality tests guard the output. Leave
  `rewrite_line` alone: it is long but linear.
- **Effect:** speed possible. **Cost** M.

### R3. Allocation in the neighbour-evidence pass — speed only, possible
Helper finding.
- `neighbour_evidence/records.rs:115` allocates a `sample_shares` vector per sample per record:
  about 2,000 small allocations per line at 2,176 samples.
- `refine.rs:838,855` (`template_into`) allocates two vectors per replaced copy, in every round
  from round 2 on.
- **Proposal:** reuse scratch buffers. **Effect:** none on output; faster. **Cost** S-M.

---

## Looked at, left alone

- **VCF encoding of `PL`** when a row is all `-inf` or contains NaN: `f64::max` ignores NaN, and
  `min(i32::MAX)` turns a NaN gap into the cap. It never writes `NaN`.
- **`AF`** when the expected copies sum to 0 writes `0.000000`, not NaN (`encode.rs:340`).
- **Multi-allelic genotype order** for `PL`/`GP`: `GenotypeTable` builds VCF order and two tests
  pin it (`calling/genotype_table.rs:616,640`). After an allele is dropped, the kept rows stay in
  VCF order (`assemble.rs:586-628`).
- **`GT` and `REPCN` alignment**: `Genotype` sorts its alleles at construction, so both fields walk
  the same order.
- **Padding at a contig's first base and after a widened STR**: the right-padding and
  dropped-deletion cases are handled and tested (`assemble.rs:530-660` in tests).
- **Writer order check and the one legal tie** (generic then tract) are correct. A leftover
  `<output>.tmp` after a failed run is overwritten by the next run and never renamed, so it cannot
  be mistaken for a finished file.
- **Header escaping of sample and contig names**: names come from `@RG SM`, psp headers and FASTA
  names. None of those can hold tab or newline, and empty or duplicate names are refused.
- **`regions.rs` BED parsing**: comments, `track`/`browser` lines, empty or inverted spans,
  unknown contigs, overhang (clamped) and overlap merging are all handled.
- **`regions.rs` u32 coordinates**: a contig over 4.29 Gb is refused with an error, not
  truncated. That is a limit, not a bug.
- **`memory_limit.rs`**: cgroup v1/v2 parsing and the "no limit" sentinels look right.
- **`reference_info.rs`**: `get_or_read`, `write_fai`, `read_fai_verify_in_background` and
  `IgnoreBases` are `pub` but used only inside the file. They could be `pub(crate)` or private;
  that is cosmetic.
- **`open_files.rs`**: `suggested_soft_limit` and `DESCRIPTORS_A_PSP_RUN_NEEDS_BESIDES_ITS_PSPS`
  are likewise `pub` and used only in their own file. Cosmetic.
- **`splitmix.rs` and `neighbour_evidence/random.rs`**: they share constants, but `random.rs` must
  match `phaser_reference.py`'s stream bit for bit. Keep the two contracts separate.
- **`TypedRegionCounts` and `CatalogRegionCounts` (and `charge_region` and `tally`)**: these are
  twins on purpose. The scan-against-catalog differential compares them field by field.
- **`cli` used as a library by examples** (`examples/ng_str_*` use
  `cli::estimate_parameters::open_the_censuses` and `cli::run_output`): it is fine. The command
  functions are where that setup lives.
- **`round_shape_for` reads `rayon::current_num_threads()`**, so the round shape depends on the
  thread count. That is documented as changing no answer, and the thread-invariance tests cover
  it.
- **Neighbour-evidence numerics**: emissions stay above zero, `exp(±L/2)` is bounded, a window
  needs at least two samples, and random streams are keyed by window, sample and round, so the
  output does not depend on the thread count. (Helper, traced.)
- **The window accumulator's cursors, `depth_grid.rs` zero cases, and the catalog builder's
  parallel scan order**: correct. (Helper.)
