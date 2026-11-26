# Code-quality map: `src/parameter_estimation/` except the two fits (2026-11-26)

Scope: everything under `src/parameter_estimation/` except `joint/fit.rs`, `joint/fit/`,
`joint/ssr_fit.rs` and `joint/ssr_fit/`. That is `joint/census.rs`, `census_file.rs`,
`census_moments.rs`, `contamination.rs`, `loci.rs`, `sequencing_batches.rs`, `stratum_fits.rs`,
`share_curve.rs`, `slippage_curve.rs`, `fit_trace.rs`, and the top-level `mod.rs`, `calibration.rs`,
`depth_bins.rs`, `misread_after_runs.rs`, `progress.rs`, `repeat_strata.rs`. Read from code only;
nothing was built or run.

## Summary

1. **A census with a few hundred read groups cannot be opened.** The reader takes the first
   mebibyte of a census and decodes its header and section directory from that alone. At
   tomato's 141 repeat strata each read group costs about 5.4 kB of header and directory, so a
   sample with more than about 190 read groups is refused as "not a census", and
   `estimate-parameters` stops on it. The fix is small and changes no output.
2. **The fit computes two population averages on every pass and nothing reads them.**
   `CensusMomentSums` is summed at every census position on every pass of the SNP/indel fit,
   stored on the fit's result, and never read. The report type and the test its documentation
   names do not exist. Either wire it to an output or delete it; deleting it changes no output and
   removes work from the fit's busiest loop.
3. **Two wrong-way imports, both easy to remove.** `misread_after_runs.rs` and `census_file.rs`
   import `calling::likelihood::run_context`, which depends only on `float` and `types`. Moving
   that file out of `calling` removes the only production import from this area into `calling`,
   and with it the two-way dependency. `progress.rs` imports one string constant from
   `cli::run_output`. Moving that constant down, together with three formatting helpers that
   `cli::run_output` duplicates, fixes both problems at once.
4. **Contamination may be slow at thousands of samples.** It builds a sample-by-sample matrix
   (markers × samples² work) and diagonalises it with Jacobi rotation (samples³ work per sweep),
   both on one thread. The code comment says "fifty by fifty at most"; the caller commits to
   several thousand samples. This needs a timing from a large run before anything changes.
5. **Code and documentation left over from two deleted designs.** These are the whole-genome
   histogram route (removed 2026-09-11) and the repeat-tract "climb" fit (removed 2026-10-16).
   What remains: an unreachable refusal in cohort assembly, about 150 test-only lines of
   `depth_bins.rs`, and about a dozen comments describing code that no longer exists.
6. **Duplication:** three copies of `ln_sum_exp`, and the progress-line formatting helpers
   (finding 3).

---

## Bugs

### B1. A sample with many read groups has a census the reader refuses as malformed

- **Where:** `src/parameter_estimation/joint/census_file.rs:454-502` (`open_census_within`),
  `:511` (`HEAD_READ_BYTES = 1 << 20`). Called from `src/run/census_cohort.rs:183`, which every
  psp-mode fit goes through.
- **What:** the reader reads `min(1 MiB, census length)` bytes, then decodes the whole header and
  the whole section directory from that buffer. A census is larger than that buffer, so if the
  header and directory are longer than 1 MiB the cursor runs off the end. `Cursor::take` then
  returns `CensusError::Malformed`, and the user sees "the stream is not a census, or is of a
  version this build does not read". Size per read group, from `encode_header` and `encode_key`:
  - **header:** about 1,100 bytes. Most of it is the misread-after-runs table, 4 + 4 + 19 × 56
    bytes. The rest is the read group's name entry and its claimed-error totals (28 bytes).
  - **directory:** 21 bytes for the read group's ordinary-position section, plus 30 bytes for
    each stratum's section.

  At 141 strata that is about 5.4 kB a read group, so 1 MiB holds about 190 read groups. The
  constant's own comment says "a megabyte covers a directory of about 40,000 sections". The limit
  was known; nothing handles crossing it. The archive survey quoted in `contamination.rs` found a
  sample with 42 libraries; at several lanes a library, its read-group count approaches this.
- **Failure scenario:** a psp for one sample with 200 lane-level read groups. `generate-psps`
  writes its census without complaint. `estimate-parameters` then fails on that sample with
  `CensusNotRead … malformed`, and `regenerate-census` writes the same census again.
- **Proposal:** when decoding fails and the head buffer is shorter than `census.len()`, read
  again with a buffer twice as large (capped at `census.len()`) and retry. A cleaner alternative
  writes the header length up front, but that is a format version bump; avoid it.
- **Downside:** none worth stating. The retry costs one extra read only for censuses that need it.
- **Effect:** none on output, speed or memory for any census that opens today.
- **Cost:** S. **Confidence:** high that the limit is hit at about 190 read groups with 141
  strata (the path is traced end to end). The human stratum count was not checked; with more
  strata the threshold is lower.

### B2. After B1 is fixed, the directory overlap check is quadratic

- **Where:** `census_file.rs:726-735` (`decode_directory`).
- **What:** every directory entry is compared with every later entry. Today the 1 MiB head caps
  this at about 35,000 sections, which is about 6 × 10⁸ comparisons for one sample. With
  hundreds of read groups (once B1 is fixed) it is 10⁹ to 10¹⁰ per sample, repeated for every
  sample of the cohort.
- **Proposal:** sort the extents by offset and compare each one with the next.
- **Downside:** none. **Effect:** none on output. **Cost:** S. **Confidence:** high.

### B3. A census can lose a tract's "reached but crossed no read" count at the start of a walk

- **Where:** `src/parameter_estimation/joint/census.rs:2893-2907` (`CensusWriter::add_ssr`).
- **What:** the reads that reached a tract and produced no observation are added to
  `covering_not_crossing`, charged to "the lowest-numbered group" in `self.ssr`. But `self.ssr`
  gains a read group only when one of that group's observations lands on a kept tract. So
  until one does, the map is empty, and the count and the walked bit are both dropped. The
  generic half avoids this by carrying the declared read groups from construction
  (`read_groups`, line 2426). That count feeds the crossing share in
  `ssr_fit/tract_exclusion.rs:175-184`, and the crossing share decides which tracts are left out
  as duplicated.
- **Failure scenario:** the first kept tracts a sample's walk reaches have reads that stop inside
  the tract and no complete observation. Those reads vanish from their stratum's "reaching"
  total. The size is a handful of reads per sample, so the effect on any stratum is tiny, but it
  depends on which tract a walk happens to meet first.
- **Proposal:** start `ssr` with an empty entry for every declared read group, as `generic` is
  handled. The count then always lands on the lowest declared group.
- **Related documentation error:** the comments claim the walked bit separates "never walked (a
  bug)" from "walked, no read" at a tract (`census.rs:22-25`, `:683-686`, `:849`). It does not.
  `mark_walked` marks ordinary positions only, so a tract no read reached is reported as
  `NeverWalked`. No fit acts on the difference: `ssr_fit.rs:1345` treats everything that is not
  `Crossed` alike. Correct the comments rather than the code.
- **Downside:** changes census bytes for affected samples.
- **Effect:** **output can change**, by a few reads in one stratum's reaching total, and so in
  rare cases which tracts the duplicate filter leaves out. Needs an approved-change note.
- **Cost:** S. **Confidence:** medium (traced; the size is reasoned, not measured).

### B4. Contamination's ancestry axes may take a long time at thousands of samples (performance, suspected)

- **Where:** `contamination.rs:998-1019` (`ancestry_coordinates`), `:1412-1469`
  (`leading_eigenvectors`).
- **What:** two steps grow fast with the number of samples, and both run on one thread:
  - **The similarity matrix** costs markers × samples² / 2 multiply-adds. At 50,000 markers and
    2,000 samples that is 10¹¹.
  - **The eigenvectors** come from cyclic Jacobi rotation, which costs about 6 × samples³ per
    sweep (5 × 10¹⁰ at 2,000 samples), for up to 100 sweeps. The rotation reads the matrix down
    its columns, which is the slow direction for memory.

  The doc comment says "this is fifty by fifty at most". Only 4 eigenvectors are needed.
- **Proposal:** first check the timestamps on the `contamination:` progress line in the log of the
  largest run so far (about 2,169 samples, per `estimation_memory.md`). If the step takes minutes:
  - **Matrix:** build it in parallel over rows. Each cell keeps its marker order, so the result
    is bit-identical.
  - **Eigenvectors:** compute only the 4 leading ones by subspace iteration.
- **Downside:** the eigenvector change alters the coordinates in their last digits, so
  contamination values can change in their last digits.
- **Effect:** the parallel matrix changes no output. **The eigenvector change can change output.**
- **Cost:** M. **Confidence:** low to medium. The complexity is certain; whether it matters in
  practice has not been measured.

---

## Architecture

### A1. Move `run_context` out of `calling`, so this area imports nothing from it

- **Where:** `misread_after_runs.rs:17-20` and `census_file.rs:57` import
  `crate::calling::likelihood::run_context` (547 lines, `src/calling/likelihood/run_context.rs`).
  Its only own imports are `crate::float` and `crate::types::ErrorRate`.
- **What:** this is the only production import from this area into `calling`. Everything else
  in the lead's count of 8 is a doc link: `misread_after_runs.rs:12`, `sequencing_batches.rs:190`,
  `stratum_fits.rs:408,416`, and a test at `census_file.rs:1985`. Because `calling` imports
  `parameter_estimation` heavily, these imports make the two modules depend on each other.
  `stratum_fits.rs:1486` still says "`parameter_estimation` imports nothing from `calling`", which
  shows that was the intended rule.
- **Proposal:** move `run_context.rs` to a leaf module both sides import, for example
  `src/run_context.rs` or a `src/read_error_context/` module. Its subject is one rule both sides
  share: which table entry a base after a homopolymer run belongs to. Nine files change their
  `use` line. The alternative is to put it under `parameter_estimation`, as `depth_bins` and
  `repeat_strata` already are ("shared vocabulary that calling reads too", `mod.rs:32-43`). That
  is cheaper, but it files the charging rule `calling` uses on candidates under the pre-pass.
- **Downside:** one more top-level module, and `run_context`'s documentation is written from the
  calling side and would need one paragraph reworded.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### A2. Move the progress-line target and formatting helpers below `cli`

- **Where:** `progress.rs:25` imports `crate::cli::run_output::PROGRESS`, the event target string
  `"pvc::progress"` defined at `cli/run_output.rs:44`. `run/calling_progress.rs:23` does the same
  from `run`.
- **What:** a library module depends on the command-line layer for one string. The same pair
  of files also holds three duplicated helpers (D1).
- **Proposal:** make a small crate-level leaf module, for example `src/progress.rs`, that holds:
  - the `PROGRESS` target;
  - `StageProgress`, `note` and `detail` (moved from `parameter_estimation/progress.rs`);
  - one copy of each formatter: size, duration, and the memory reading from
    `/proc/self/status`.

  `cli::run_output` (the subscriber that draws the line) imports the target from it, which is
  the right direction. So do `run::calling_progress` and `cli/estimate_parameters.rs`.
- **Downside:** the target name is a contract between the code that emits events and the code
  that draws them, and it would no longer sit beside the drawing code. A doc comment on each side
  pointing at the other covers that.
- **Effect:** none, provided both duration styles are kept (see D1). **Cost:** S.
  **Confidence:** high.

### A3. `crate::run` in this area is documentation links only: leave it

- **Where:** `census_file.rs:92, 945, 1109`.
- **What:** these are rustdoc links to `run::CensusVerdict` and `run::census_freshness`. The two
  real calls into `run` from `parameter_estimation` are in `fit.rs:5881, 8567`, which belong to
  the other agent's area.
- **Proposal:** leave them; they cost no build dependency. **Effect:** none.

### A4. Put the census writer in a child module of `census`

- **Where:** `census.rs:2382-3031` (`CensusWriter` and `thin_to_cap`, about 650 lines).
- **What:** the writer is the only part of `census.rs` that needs `locus_generation`,
  `calibration`, `misread_after_runs` and `ref_seq`. The other 2,350 production lines are data
  types, recording terms, lending of sections, and cohort assembly, which is what the fit,
  the file codec and contamination read. With the writer separate, the census data types no
  longer depend on how a walk produces loci.
- **Proposal:** `census/writer.rs` as `mod writer;` inside `census`. A child module sees its
  parent's private fields, so no visibility changes are needed. Re-export `CensusWriter`.
  **Split nothing else**: the remaining sections of the file are one idea and already separated
  by banners.
- **Downside:** someone reading `add_generic` has to open a second file to see the field layout
  it writes into. The benefit is moderate.
- **Effect:** none. **Cost:** S. **Confidence:** medium.

---

## Duplication

### D1. Progress formatting exists twice, in `parameter_estimation::progress` and `cli::run_output`

- **Where and what:**
  - **size:** `progress.rs:138` (`size`) and `cli/run_output.rs:369` (`size_in_bytes`) are the
    same function: whole MiB below 1 GiB, GiB to one decimal above.
  - **duration:** `progress.rs:155` (`duration`) and `cli/run_output.rs:495` (`elapsed_text`) are
    the same except below one second. The first prints `350ms`, the second `0s`.
  - **memory:** `progress.rs:93-123` (`memory`, `memory_in_status`) and
    `cli/run_output.rs:355` (`peak_resident_bytes`) both parse `VmHWM` out of
    `/proc/self/status`.
- **Proposal:** one copy of each, in A2's module. Keep `duration` and `elapsed_text` as two
  functions, or one with a flag for milliseconds.
- **Downside:** merging the two duration styles into one changes log text for sub-second stages.
  Do that only if the owner accepts it.
- **Effect:** none if both styles are kept. **Cost:** S. **Confidence:** high.

### D2. `ln_sum_exp`: three identical copies

- **Where:** `contamination.rs:1404`, `joint/fit.rs:1586`, `calling/quality/mod.rs:856`
  (`log_sum_exp_over`). All three are the same two-pass arithmetic: find the largest term, then
  `exp` the differences and sum.
- **Proposal:** one `pub fn ln_sum_exp` in `crate::float`. The arithmetic is identical, so
  results are bit-identical.
- **Downside:** touches `fit.rs`, which is the other agent's file; do it as one change across
  both areas.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### D3. Golden-section search exists twice and should stay that way

`contamination.rs:1303` (30 fixed steps, returns the better of the two inner probes) and
`fit.rs:5214` (stops on a width test, returns the bracket's middle) differ in what they return.
Merging them changes outputs. **Leave both.**

---

## Readability (and dead code)

### R1. `CensusMomentSums` is summed on every fit pass and never read

- **Where:** `joint/census_moments.rs` (308 production lines); summed at `fit.rs:4722` and
  `:4994`, merged at `:1974`, stored in `JointFit::census_moments` at `fit.rs:360`.
- **What:** the type has three methods, `over`, `add_position` and `merge`, and no accessor.
  Its fields are private, and no code outside the file reads them. The documentation names
  `CensusMomentSums::finish`, `CensusMomentsReport::of`, `CensusMoments::from_posteriors` and a
  test `the_summed_moments_and_the_stored_array_agree`. None of the four exists. `mod.rs:27` says
  "`joint::census_moments`'s output says so", but there is no output.
  `add_position` loops over every sample at every census position on every pass of the fit.
- **Proposal:** the owner decides between two options:
  - **delete** the accumulation and the module. A future spec-driven report can come back with
    its reader;
  - **finish it:** add `finish` and a report that writes the two averages.
- **Downside of deleting:** the spec (`ordinary_site_prior_moments.md`) asks for these averages.
  Deleting postpones that.
- **Effect:** none on output. Deleting removes an O(samples) loop from every census position on
  every pass of the fit, which may make the fit faster (not measured). Check first that no
  `Debug` dump of `JointFit` reaches a file a user reads.
- **Cost:** S. **Confidence:** high.

### R2. The "two samples claim one read group" refusal cannot happen any more

- **Where:** `census.rs:1987-2016` (`SharedReadGroup`), `:2027` (variant), `:2038`, `:2053-2057`,
  `:2241-2269` (`read_groups_of`).
- **What:** `CohortCensusEvidence::new` renumbers every sample onto fresh run-wide identifiers
  (`onto_run_wide_read_groups`, `:2202`) before `read_groups_of` runs, so two samples can no
  longer hold the same number. The refusal and its 20 lines of explanation describe the state
  before renumbering existed, and no test reaches it.
- **Proposal:** make `read_groups_of` a plain union with a `debug_assert!` that the samples'
  groups do not overlap. Delete the type, the variant and the `From` impl, and update `new`'s
  `# Errors` section.
- **Downside:** none. **Effect:** none. **Cost:** S. **Confidence:** high.

### R3. `depth_bins.rs` is still written for the deleted histogram route

- **Where:** `depth_bins.rs`: module doc `:1-42`; `DepthBinEdges` doc `:166-190` (the shared
  `Arc` that lets `merge` check two histograms use the same ladder); `row_starts` `:201` (marked
  `expect(dead_code)`); `Default` `:204`; `ladder_tops` `:242`; `new` `:275`; `row_start` `:360`;
  `cell_count` `:394`; `max_depth` `:443`.
- **What:** in production only `for_census()` is called (`run/gatherer.rs`,
  `calling/parameters_file/bindings.rs:987`). The twenty-bin ladder, its cell offsets and the
  Arc-identity argument serve only tests. One part still matters: `for_census` takes its ten
  rungs above 124 from the twenty-bin ladder's growth ratio (`widening_ratio(8, 124, 11)`).
- **Proposal:**
  - keep `for_census`, `bin_for`, `depth_range`, `bins` and `bin_count`;
  - keep the three constants as the named inputs to the ratio. **Changing them changes
    `DepthLadderDigest` and makes every census on disk unreadable**, and the existing test that
    pins the census tops must stay;
  - delete the rest, and rewrite the module doc around the census ladder.
- **Downside:** the research justification for the 20-bin ladder moves to the report it came
  from.
- **Effect:** none, provided the census tops are unchanged. **Cost:** M (mostly the
  documentation). **Confidence:** high.

### R4. Comments that describe code that no longer exists

Each takes a line or two to fix; together they mislead a new reader in about a dozen places:

| Where | Says | True now |
|---|---|---|
| `census.rs:32-38` | the byte-level census file "is the next unit" | `census_file.rs` exists |
| `census.rs:124-126` | the ladder is "shared with the histogram route" | route deleted |
| `census.rs:467` | buckets `-4 … +4` | `±8` since 2026-08-20 |
| `census.rs:2388-2390` | the writer is fed "the same locus stream the histogram accumulators are handed" | no histograms |
| `census_file.rs:35` | `version u16 1` | 5 |
| `census_file.rs:48-52` | reading one section "is the next unit of work" | `Sections::fill` does it |
| `contamination.rs:190`, `:1352-1355` | "five-bit code", "above eight reads … a range" | 8-bit code, exact to 124 |
| `contamination.rs:237` | `examples/ng_joint_contamination_control.rs` | no such file (`ng_joint_contamination_harness.rs` exists) |
| `contamination.rs:1413` | "fifty by fifty at most" | cohorts of thousands (see B4) |
| `stratum_fits.rs:1486` | "`parameter_estimation` imports nothing from `calling`" | false until A1 is done |
| `calibration.rs:277` | "the same observations the read-group histogram counts" | no histogram |
| `census_moments.rs:18-24, 53, 73`; `mod.rs:27` | see R1 | |
| `mod.rs:165-167` (test doc) | provenances are "not a scale", none is better | `weaker_of` (`:94`) ranks them |

**Effect:** none. **Cost:** S.

### R5. Contamination carries five research options that no command can select

- **Where:** `contamination.rs:166-295` (`ContaminationConfig`, `ContaminationGrain`,
  `OwnCoordinates`), `:350-402` (`Units`), `:1184-1299` (`fit_alpha`).
- **What:** these settings are set only by `contamination.rs`'s own tests:
  - the `Sample` grain;
  - `OwnCoordinates::UndoneByAlpha` and `MaximisedFreely`;
  - `leave_self_out = true`;
  - `integrate_over_depth_bin = false`;
  - `weight_by_posterior = false`.

  The only production caller, `cli/estimate_contamination.rs:473`, takes the defaults. The
  measurements that rejected each option are already written into the doc comments. Without
  them, `fit_alpha` would be about 25 lines instead of 115, and `Units` could be removed.
- **Proposal:** owner decision. Remove the options and keep their measurements in the module doc
  and the reports, or keep them as a research harness.
- **Downside:** re-measuring a rejected option means writing it again. Tests pinning "the sample
  grain reproduces the earlier numbers" go with it.
- **Effect:** none on the default path. **Cost:** M. **Confidence:** high that nothing reaches
  them.

### R6. `PackedDepthCodes` bit-packs codes that are exactly one byte

- **Where:** `census.rs:185-291`.
- **What:** `set` and `get` handle a code that straddles two bytes, through a 16-bit window,
  a mask and a shift. `DEPTH_CODE_BITS` has been 8 since format version 2, so every code is one
  byte and `bits[index]` is the whole job. The fit calls `get` once per position per read group
  per pass (`fit.rs:1087, 1137, 1221`).
- **Proposal:** store a `Vec<u8>` and index it. Keep the compile-time assertion that the ladder
  fits under the sentinel.
- **Downside:** if the code is ever widened again it has to be generalised again, but that is a
  format version bump anyway.
- **Effect:** none on output, and possibly slightly faster (unmeasured). **Cost:** S.
  **Confidence:** high.

### R7. Support for reading old repeat-tract parameter files is the owner's call

- **Where:** `share_curve.rs` (all 188 production lines), most of `slippage_curve.rs`,
  `stratum_fits.rs`'s pooled rung (`allele_length_variation_by_period`) and the `Curve` and
  `Blend` provenance variants.
- **What:** all of this exists only to read parameters files written by the repeat-tract fit
  deleted on 2026-10-16. The matching code in `calling/parameters_file/` is outside this area.
- **Proposal:** ask the owner whether files from before 2026-10-16 still need to be read. If
  not, remove the code on both sides. **This changes which input files are accepted**, so it
  needs approval.
- **Cost:** M. **Confidence:** high that nothing else produces these values.

### R8. Small items

- `contamination.rs:633`: the bare `100` (fewest markers) should be a named constant, like its
  neighbours `MIN_SAMPLES_WITH_DATA` and `MIN_FREQUENCY` at `:420-423`.
- `census.rs:1332-1338` (`Sections::read_groups_used`) repeats the match inside
  `Sections::keys` (`:1372`); write it as `self.keys().map(SectionKey::read_group).collect()`.

---

## Looked at, left alone

- `census_file.rs` decoding checks the same things `GenericEvidence::from_parts` and
  `SsrEvidence::from_parts` assert. That is deliberate: a file's mistake is an error, a caller's
  mistake is a panic.
- `AlleleLengthVariation::fitted` and `checked_allele_length_variation` both check the class
  count. One runs per lookup and is cheap; the other is the gather-time check. Documented.
- Three spellings of "stratum" (`repeat_strata::Stratum`, `census::Stratum`, the catalog's
  `(u8, u64)`). The trade-off is documented at `repeat_strata.rs:13-19`.
- `repeat_strata.rs` is read only by `calling` and `run`. It is shared vocabulary, by the rule
  `mod.rs:32-43` states.
- `fit_trace.rs` keeps a process-wide sink opened from an environment variable. It is an opt-in
  diagnostic and off by default.
- `sequencing_batches.rs`'s test-only `declared` and `SequencingBatchError` (about 250 lines).
  The feature is specified, unbuilt, and labelled as such in the module doc.
- `loci.rs:576` `SelectionTerms::first_disagreement` writes the field names out again rather than
  using `census::SELECTION_FIELDS`. It is test-only.
- `calibration.rs`'s fixed-point `i128` sum is justified at length and correct.
- `CensusWriter::add_generic` puts an undeclared read group's non-reference reads in `pending`
  (`census.rs:2769-2783`). It cannot happen, because `read/input/read_groups.rs` refuses a record
  whose read group the file does not declare.
- `CensusWriter::mark_walked` is O(kept positions × read groups) with a map lookup each, but runs
  once per walk.
- `contamination.rs`'s two passes over position-length arrays: a measured memory design
  (`estimation_memory.md` §6).
- `stratum_fits.rs`'s nested `checked_spectrum` module exists so that the private fields cannot
  be bypassed; it is intentional.
- `Sections::Backed` opens the file once per lending call rather than holding it open. This is
  deliberate, to stay under the open-file limit at thousands of samples.
