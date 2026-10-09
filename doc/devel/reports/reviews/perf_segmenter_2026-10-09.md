# Performance Review: segmenter
**Date:** 2026-10-09
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the genome segmenter — typing a reference into ordinary ground, STR segments, clusters of repeats and satellites from the stored repeat catalog, and building that catalog
**Verdict:** Apply the listed wins
**Hot-path evidence:** sampling profiles (`perf`, cpu-clock) on GRCh38 and tomato, before and after the fixes; interleaved wall-clock timings; allocation and comparison counts

---

## What this review found

**Typing GRCh38 from its catalog takes 39.7 seconds a pass today; with four measured fixes it takes 3.8 seconds, and the typed regions come out byte for byte the same.** The fixes are on branch `segmenter-perf-combined` (commit `d4027a8e`), ready to be reviewed and merged one at a time.

Three costs made typing slow, and none of them is the statistics of typing:

1. **The catalog file's index is read again for every contig.** A Parquet file ends with an index describing every block of rows in it. The reader re-reads and re-decodes that whole index once per contig. GRCh38 has 2,580 contigs, so the index was decoded 2,580 times a pass; tomato has 13 contigs, which is why tomato never showed it.
2. **Asking for part of the genome cost more than asking for all of it.** When a run is given a BED file, several steps compared every item on a contig (every repeat, every typed region, every locus) with every BED interval on that contig. The 481,622 intervals of the GIAB HG002 confident regions took **476 seconds** to type, against 40 seconds for the whole genome. With the fix they take **3.0 seconds**.
3. **Absorbing STR loci into nearby satellites compared every locus with every satellite, repeatedly.** On GRCh38 that was 9.3 billion comparisons to absorb 15,476 loci.

**How often a run pays this.** Typing happens once per command, never per region or per worker. But three commands type the genome three times: `generate-psps`, `estimate-parameters` and `regenerate-census` each run one typing pass for the run's ground and two more to choose the census positions. And `generate-psps`'s own help text recommends **one invocation per sample** to spread a cohort over cores (`src/cli/generate_psps.rs:126`). So in the psp workflow, typing is paid three times **per sample**. Measured per run of `generate-psps`:

| reference and ground | today | with the fixes |
|---|---|---|
| tomato, whole genome | 2.8 s | 1.8 s |
| GRCh38, whole genome | 123 s | 12 s |
| GRCh38, GIAB confident regions (481,622 intervals) | about 24 minutes (476 s a pass, three passes; extrapolated) | 9 s (measured) |

On a 63-sample human cohort called through psps, that last row is 25 hours of setup before the fixes and under 10 minutes after.

**The "minutes" the owner measured were mostly not typing.** The dump example `ng_typed_region_dump` writes one line at a time to standard output, and inside the dev container every line becomes a separate write to the shared filesystem. GRCh38: 31.6 s with output discarded, **226.7 s** writing to a file in the project, 32.8 s when the same output goes through `| cat` first (identical bytes). Tomato: 1.0 s, 43.6 s and 1.1 s. The shipped `pop_var_caller type-regions` command buffers its output and does not have this problem.

---

## 1. Scope and constraints

- **What was reviewed:** a module group, the genome segmenter. Its entry points, confirmed from the code:
  - **Typing (paid per run):** `RepeatCatalog::genome_segments`, `RepeatCatalog::str_loci` and `RepeatCatalog::sample_loci_per_stratum` ([reader.rs](../../../../src/repeat_catalog/reader.rs#L141)). Per contig they read the contig's rows with `rows_overlapping` ([parquet_file.rs](../../../../src/repeat_catalog/parquet_file.rs#L451)) and type them with `segments_of_contig[_in]` / `loci_of_contig[_in]` ([segments.rs](../../../../src/repeat_catalog/segments.rs#L35)), which call `prefilter`, `admit`, `split_bundles`, `coverage_runs`, `resolve_features` and `fill_generic_gaps` ([region_typing/mod.rs](../../../../src/region_typing/mod.rs#L742), [segment_criteria.rs](../../../../src/region_typing/segment_criteria.rs#L687)).
  - **Who calls typing:** `run_ground::segments_over` / `segments_cut_from` ([run_ground.rs](../../../../src/cli/run_ground.rs#L282)), collected into `Segmentation::build` ([run/segments.rs](../../../../src/run/segments.rs#L76)) — one pass in `call-from-alignments`, `call-from-psps`, `generate-psps`, `estimate-parameters`, `regenerate-census`. `CensusPlan::of_run` → `select_kept_loci` ([loci.rs](../../../../src/parameter_estimation/joint/loci.rs#L827)) adds a second `genome_segments` pass and a `str_loci` pass in `generate-psps`, `estimate-parameters` and `regenerate-census`. `estimate-contamination` runs `select_kept_loci` and a third `genome_segments`.
  - **Catalog build (once per reference, then stored beside the FASTA):** `pop_var_caller repeat-catalog` → `RepeatCatalogBuilder` ([builder.rs](../../../../src/repeat_catalog/builder.rs#L93)), the scanner `find_tandem_repeats` ([tandem_repeat.rs](../../../../src/tandem_repeat.rs#L483)) and the Parquet writer.
- **Reviewed against:** commit `39336d54` on branch `segmenter-perf` (= `3a31a724` of `tract-ceiling` plus the timing example this review added).
- **Targets:** none are written in the specs. The working target used here: typing is fixed setup before the first read is decoded, so it should cost seconds, not minutes, on a 3 Gb reference with thousands of contigs, and it must not grow faster than linearly in the number of BED intervals or contigs. Inputs that must work: references from a few megabases to 3 Gb, from 13 contigs (tomato) to 2,580 (GRCh38 analysis set) and the 10,000+ scaffolds common in plant and fungal assemblies; ground from the whole genome to BEDs of hundreds of thousands of intervals (GIAB confident regions: 481,622; exome kits: about 200,000). Output must stay byte-identical.
- **Hardware:** Apple `container` Linux VM on the Mac, arm64, 8 CPUs (sub-agents capped at 4), release build with fat LTO. No hardware counters are available on this machine.
- **Hot-path evidence available:** `perf record -e cpu-clock --call-graph dwarf` profiles of a whole-genome typing pass on GRCh38 (5,055 samples) and tomato (281 samples) before the fixes, and on GRCh38 after them; wall-clock timings run interleaved; allocation counts (counting allocator) and comparison counts (temporary counters). All raw output is in `tmp/perf_review_2026-10-09_segmenter/`.
- **In-scope files:** `src/repeat_catalog/{reader,segments,parquet_file,mod,builder,criteria,row,strata,tally}.rs`, `src/region_typing/{mod,segment_criteria}.rs`, `src/tandem_repeat.rs`, `src/run/segments.rs`, `src/cli/run_ground.rs`, `src/parameter_estimation/joint/loci.rs` (`select_kept_loci`).
- **Deliberately out of scope:** everything downstream of the typed regions (read walking, calling, cohort merge); `partition_resident`, the live scan kept only as a test oracle; test code; `examples/ng_typed_region_dump.rs`, whose output cost is explained above and noted in section 6.
- **Categories dispatched** (each sub-agent worked in its own git worktree and measured its own fixes):
  - `methodology` — always; there was no segmenter benchmark.
  - `io_and_syscalls` — the profile put 83% of GRCh38 time in reading the catalog file.
  - `hot_loops` — satellite absorption (11% of GRCh38, 20% of tomato) and the per-interval loops.
  - `allocations` — malloc and free were 27% of GRCh38 self time.
  - `data_layout` — a run holds the whole genome's typed regions for its lifetime, a cost the code's own comment calls unmeasured.
  - `concurrency` — typing runs on one thread while the run has a thread budget.
  - Not dispatched: `float_portability`. No proposal changes a floating-point computation; the only float in typing is a purity threshold compared, not computed.

## 2. Verdict

**Apply the listed wins.** Three Hot-path findings are measured, each with identical output on every input tried and with the module tests passing. Together with two smaller Likely findings they take a GRCh38 typing pass from 39.7 s to 3.8 s and the 481,622-interval GIAB pass from 476 s to 3.0 s.

The combined result was re-measured by the orchestrator in one container, baseline and fixed builds alternating, three rounds each (seconds; "census" is the two extra passes `generate-psps` makes):

| input | pass | baseline (3 rounds) | all fixes (3 rounds) |
|---|---|---|---|
| tomato, whole genome | segments | 0.933 / 0.912 / 0.875 | 0.601 / 0.602 / 0.590 |
| tomato, whole genome | census | 1.886 / 1.885 / 1.849 | 1.242 / 1.235 / 1.230 |
| GRCh38, whole genome | segments | 39.41 / 39.68 / 39.69 | 3.93 / 3.83 / 3.80 |
| GRCh38, whole genome | census | 81.89 / 83.79 / 82.99 | 7.99 / 8.13 / 8.07 |
| GRCh38, first 10,000 GIAB intervals | segments | 5.31 / 5.37 / 5.32 | 0.101 / 0.095 / 0.097 |
| GRCh38, all 481,622 GIAB intervals | segments | 476.2 (one earlier run) | 2.97 / 2.98 |
| GRCh38, all 481,622 GIAB intervals | census | not run | 5.96 / 5.97 |

Output identity of the combined build: the typed-region dump for tomato (976,237 lines) and GRCh38 (4,258,536 lines), whole genome, is identical (`cmp`) to the baseline's; `cargo test --release --lib -- repeat_catalog region_typing run::segments` passes 177 tests and `--test ng_repeat_catalog_differential` passes 11. The sub-agents also checked each fix separately on the tomato STR BED and on GRCh38 with 1,000, 10,000 and all 481,622 GIAB intervals.

What is left after the fixes, from a profile of the fixed build on GRCh38 (3.8 s, `perf_combined_*.txt`): 75% is reading the rows — decompressing Parquet pages is 31% on its own — and about 25% is typing. So the next steps are reading less often (finding S1) and reading contigs in parallel (finding L1), not further work on the typing code.

## 3. Measurement plan

In order:

1. **Gate each fix on output identity before timing it.** For each fix, dump the typed regions with `ng_typed_region_dump <fasta> <bed> calling | cat > file` (the `| cat` matters, see above) on tomato whole genome, tomato `benchmarks/ssr_tomato1/ssr_regions.bed`, GRCh38 whole genome and the GRCh38 first 10,000 GIAB intervals, and `cmp` against the baseline. Then `cargo test --release --lib -- repeat_catalog region_typing run::segments` and `--test ng_repeat_catalog_differential`.
2. **Gate on counts where the fix removes countable work** — they are identical run to run, where wall time on this machine varied by 28–51% between runs while other containers were busy (methodology finding):
   - index decodes per pass: 2,580 → 1 on GRCh38 (H1);
   - feature-against-satellite comparisons: 9.3 billion → one sort of about 2.3 million spans on GRCh38 (H3);
   - allocations in a tomato typing pass: 1,933,161 → 984,799 (L4), → 521,511 more with L5.
3. **Wall time, interleaved, in one container:** `tmp/perf_review_2026-10-09_segmenter/combined_timing.sh` is the script that produced the table in section 2; it runs `examples/segmenter_timing.rs` on tomato and GRCh38, whole genome and the 10k GIAB BED, alternating two binaries. A fix is worth merging when its wall-time change is larger than the round-to-round spread on a quiet machine (about 1% in the section 2 run).
4. **Commit a benchmark that does not need the GRCh38 file.** The methodology sub-agent wrote a criterion bench, `benches/segmenter_perf.rs` (`tmp/perf_review_2026-10-09_segmenter/methodology_segmenter_bench.patch`), that builds synthetic catalogs and reproduces both growth patterns: on 8 Mb of sequence, typing took 44 ms with 8 contigs, 1.22 s with 256 and 29.0 s with 2,048 (per-contig cost 5.5 → 14.2 ms); with all intervals on one contig, 23 µs an interval at 1,000 intervals and 82 µs at 64,000. Proposed merge gates: after H1, the per-contig cost at 2,048 contigs within 1.5× of that at 8; after H2, the per-interval cost at 64,000 within 1.5× of that at 1,000.
5. **Re-profile with the allocator the binary ships.** `segmenter_timing` (and `ng_typed_region_dump`) declare no allocator, so this review's profiles ran on glibc malloc; `pop_var_caller` installs mimalloc ([main.rs](../../../../src/main.rs#L19)). Before the fixes, GRCh38 typing took 51–53 s under glibc and 43–45 s under mimalloc on a loaded machine (allocations sub-agent). After H1 almost all of those allocations are gone, so this matters less, but any further allocation work should be measured under mimalloc (the methodology patch adds it to the example).

## 4. Build / toolchain configuration

- **The release profile suits this code** (`lto = "fat"`, `codegen-units = 1`, line tables only, `panic = "abort"`). No build setting addresses a cost that grows with contigs × index size or with items × intervals.
- **`[profile.bench]` sets `debug = true` ([Cargo.toml](../../../../Cargo.toml#L58)), and with fat LTO the bench binary did not finish linking in a 10 GB container after 50 minutes.** With `debug = "line-tables-only"` it linked in 4 min 13 s in 12 GB. Recommendation: set `debug = "line-tables-only"` in `[profile.bench]`; `perf` and `samply` need only line tables to name source lines.

## 5. Code-level findings

### Hot-path

#### H1: [parquet_file.rs:461](../../../../src/repeat_catalog/parquet_file.rs#L461) — The catalog's index is re-read and re-decoded for every contig

- **Confidence:** High
- **Hot-path evidence:** GRCh38 whole-genome typing profile, inclusive time: `83.03% rows_overlapping`, `62.26% ArrowReaderBuilder::try_new_with_options`, `32.58% decode_footer_metadata`, `29.59% parse_index_data`. Reading the rows alone (no typing) took 0.59 s for tomato's 6.4 million rows but 26.3 s for GRCh38's 23.6 million — 92 ns a row against 1.1 µs, the difference being about 9 ms per contig × 2,580 contigs. Allocation count: 512 million of the pass's 521 million allocations (98%) come from these decodes, 198,496 each.
- **Pattern matched:** read once what every read needs (io_and_syscalls).
- **Mechanism:** `rows_overlapping` opens the file and builds a new reader with the page index required, for every contig. Each build decodes the file's whole footer (2.6 MB on GRCh38) and the page index of every row group in the file, then uses one contig's row group. The work per pass is contigs × whole-index size, so it grows with the square of the contig count. The fix loads the metadata once when the catalog is opened (`ArrowReaderMetadata::load`), keeps it in a crate-private `CatalogIndex` on `RepeatCatalog` with a map from contig to its row groups, and builds each contig's reader with `ParquetRecordBatchReaderBuilder::new_with_metadata` (an `Arc` clone).
- **Measured (io_and_syscalls sub-agent, loaded machine, medians of 3):** GRCh38 whole genome, segments pass 49.7 → 12.2 s; rows only 45.4 → 5.8 s; census passes 107.0 → 23.7 s. Tomato and the 10k BED unchanged, as expected (they decode the index at most a few times). Output identical on five inputs.
- **Complexity cost:** one crate-private type of about 90 lines, one field on `RepeatCatalog`, one argument to `rows_overlapping`. The decoded index stays in memory for the catalog's lifetime (2.6 MB of encoded index on GRCh38; decoded size not measured). A malformed catalog is now refused when it is opened rather than at its first contig read; no test depends on the old timing.
- **Fix:** `tmp/perf_review_2026-10-09_segmenter/io_and_syscalls_cached_footer.patch`; in the combined branch as part of `d4027a8e`.

#### H2: [mod.rs:124](../../../../src/repeat_catalog/mod.rs#L124), [reader.rs:358](../../../../src/repeat_catalog/reader.rs#L358), [segments.rs:86](../../../../src/repeat_catalog/segments.rs#L86) and five more sites — Typing cost is items on a contig × BED intervals on that contig

- **Confidence:** High
- **Hot-path evidence:** typing time against the number of GIAB intervals requested, segments pass: 1,000 → 0.48 s; 10,000 → 4.8 s; 50,000 → 34.8 s; 481,622 → 476.2 s, against 30.8 s for the whole genome in the same session. `OwnedScope::spans_on` was 4.5% of the whole-genome profile and 18.7% once H1 was fixed. Tomato with the 14,455-interval STR BED took 2.9 s against 0.7 s for the whole genome.
- **Pattern matched:** recomputing an invariant in a tight loop (hot_loops); a `collect` into a new `Vec` per item (allocations).
- **Mechanism:** eight places each walk the whole interval list for each item:
  - `OwnedScope::spans_on` scans every interval of the run and copies the contig's into a new `Vec`, and `GenomeSegments::tally` calls it once per emitted locus. On the first 50,000 intervals this copied 393.6 GB in total.
  - `segments_of_contig_in` tests every typed region against every interval on the contig.
  - `admit`'s `charged` test, the `repeat_bp` sum and `loci_of_contig_in` walk every interval per row, run or locus.
  - the row filter in `rows_overlapping` walks every interval per row.
  - `contig_walk` calls `ReadScope::touches`, a scan of all intervals, once per contig: 2,580 × 481,622 ≈ 1.2 billion comparisons before typing starts.

  The fix groups the intervals by contig once, when the scope is made, and adds a 20-line `SpanIndex` whose `candidates(from, to)` returns, by two binary searches, the slice of intervals that can overlap a stretch. It only narrows: each caller keeps its own overlap test on the slice, so results cannot change; a list not in the sorted, non-overlapping shape `GenomeRegions` produces gets the old full walk.
- **Measured (hot_loops sub-agent, segments pass, 3 interleaved rounds):** 1,000 intervals 0.64–0.76 → 0.07–0.09 s; 10,000 6.6–7.3 → 0.12–0.13 s; 50,000 45.5–49.1 → 0.46–0.48 s; all 481,622: 476 → 3.93 s; tomato STR BED 3.6–3.9 → 0.64–0.67 s. Output identical (hash of every region, the tally and every locus) at 1,000, 10,000 and 481,622 intervals and on the tomato BED.
- **Complexity cost:** a `ContigSpans` struct and the `SpanIndex` type in `repeat_catalog/mod.rs`; `OwnedScope::Regions` holds a per-contig map instead of a flat list; `admit` takes an `Option<SpanIndex>`. One rule to keep: `candidates` must only narrow, never decide.
- **Fix:** `hot_loops_span_index.patch` (applies after `hot_loops_absorb_sweep.patch`); in the combined branch as part of `74b2c508`. The allocations sub-agent's smaller fix for the tally alone (`allocations_move_kind_and_contig_spans.patch`, 50,000 intervals 51.5 / 47.4 → 36.9 / 35.4 s) is superseded by this one.

#### H3: [region_typing/mod.rs:803](../../../../src/region_typing/mod.rs#L803) and [:899](../../../../src/region_typing/mod.rs#L899) — Satellite absorption compares every locus with every satellite until nothing changes

- **Confidence:** High
- **Hot-path evidence:** `absorb_into` was 11.1% of GRCh38 self time and 20.3% of tomato's; after H1 it was 53% of GRCh38 self time (io_and_syscalls sub-agent's profile of its fixed build). Counted on GRCh38: 9,321,516,976 locus-against-satellite tests to absorb 15,476 loci and bundle members into 13,058 satellites, over 2,791 passes; tomato: 326 million tests for 2,269 absorptions. Each absorption also re-sorts the satellite list.
- **Pattern matched:** an all-pairs test re-run to a fixed point where one ordered pass gives the same answer (hot_loops).
- **Mechanism:** a satellite swallows any locus or bundle member within the bundle radius of it, and the grown satellite can then swallow more; the loop repeats until nothing changes. The final state depends only on which spans are close to which, not on the order they are met, so sorting satellites, loci and members together by start and sweeping once gives the same groups. The sweep keeps a short queue of groups a later locus could still reach, because two satellites join only when they overlap or touch while a locus reaches a satellite up to the bundle radius away; the property test caught that case on its first run.
- **Measured (hot_loops sub-agent):** tomato segments pass, median of 5 interleaved rounds, 0.953 → 0.842 s; GRCh38 whole genome, 3 rounds on a loaded machine, 64.0 / 48.3 / 45.4 → 58.5 / 43.3 / 39.8 s (about 5 s a round, matching its 11% share). Output identical on all six inputs. After all fixes, absorption no longer appears among the top entries of the profile; the sort the sweep needs is 1.6%.
- **Complexity cost:** one private function of about 120 lines with its correctness argument, a small `Group` struct and a queue. The correctness of the sweep is no longer obvious from reading the loop. The old loop is kept, test-only, as the oracle of a property test over 20,000 random layouts (`absorption_by_sweep_is_the_iterations_fixed_point`); deliberately breaking the sweep's touch rule makes it fail after 31 cases.
- **Fix:** `hot_loops_absorb_sweep.patch`; in the combined branch as part of `74b2c508`.

### Likely

#### L1: [reader.rs:397](../../../../src/repeat_catalog/reader.rs#L397) and [:449](../../../../src/repeat_catalog/reader.rs#L449) — Type several contigs at once

- **Confidence:** Medium (measured, but the design question of the thread budget is open)
- **Hot-path evidence:** measured by the concurrency sub-agent with an experimental batch size B (contigs typed in parallel, handed out in contig order), 4 CPUs, loaded machine, medians of 3. GRCh38 whole genome: one at a time 50.9 s, B = 4 17.2 s; with H1's index decoded once, one at a time 11.9 s, B = 4 5.3 s. Tomato 0.98 → 0.43 s. Peak memory at B = 4: 757–776 MB against 269–282 MB.
- **Pattern matched:** independent work items run serially while the process has a thread budget (concurrency).
- **Mechanism:** each contig's typing reads only that contig's rows. The only state that crosses contigs is the running tally (sums) and the order in which loci reach the stratum sampler, and both stay on the consuming thread in contig order. After all fixes, 75% of a pass is reading and decompressing rows, which is exactly what parallel contigs overlap. Output was identical at B ∈ {1, 4, 8, 16, 64} on tomato, GRCh38 whole genome and two BEDs, including the census sample.
- **Complexity cost:** a batch-size parameter on `genome_segments` / `str_loci`, an ordered parallel map, and about 500 MB more peak memory on GRCh38 at 4 threads (about one contig's rows and working data per extra thread; chr2 alone is 1.94 million rows, 124 MB). The batch size must come from the command's thread budget: `call-from-alignments`, `call-from-psps` and `estimate-contamination` size the rayon pool to `--threads` before typing, but `generate-psps`, `estimate-parameters` and `regenerate-census` set no pool, so rayon would take every core — and `generate-psps` is meant to run as many concurrent processes, one a sample. Those three commands should pass 1 until they have a thread flag.
- **Not on the combined branch.** With H1–H3 in, it would take GRCh38 from about 4 s to about 2 s a pass at 4 threads, by the sub-agent's ratio; that is not measured on the combined build.

#### L2: [parquet_file.rs:552](../../../../src/repeat_catalog/parquet_file.rs#L552) — The row filter runs, and decodes three columns twice, when every row passes

- **Confidence:** Medium (the gain is close to the run-to-run spread of the loaded machine it was measured on)
- **Hot-path evidence:** after H1, medians of 4 interleaved runs: reading rows on whole contigs 5.77 → 4.18 s on GRCh38 and 0.85 → 0.69 s on tomato; a segments pass 11.97 → 10.99 s and 1.06 → 0.89 s. BED inputs unchanged.
- **Mechanism:** a run passes whole contigs as intervals, so the filter keeps every row but still decodes the start, end and contig columns once to test them and again to output them. The fix keeps each row group's start/end range from its statistics and skips the filter when one requested window covers it. The skip is exact because integer min/max statistics in Parquet are exact.
- **Complexity cost:** one vector in `CatalogIndex`, two statistics helpers, one branch. In the combined branch (`d4027a8e`).

#### L3: [parquet_file.rs:576](../../../../src/repeat_catalog/parquet_file.rs#L576) — Nine column type checks per row instead of per batch

- **Confidence:** Medium (measured on tomato only)
- **Hot-path evidence:** tomato profile: Arrow's per-row `is_null` / `as_any` / `type_id` / `value` accessors about 8% together. Measured: tomato rows-only median 0.766 → 0.680 s, segments 0.836 → 0.743 s (5 interleaved rounds); lost in the noise on GRCh38 before H1.
- **Mechanism:** `row_from_batch` looks up and type-checks all nine columns for every row. A struct holding the nine typed column references is built once per batch.
- **Complexity cost:** one private struct; the old function stays for a test. In the combined branch (`74b2c508`).

#### L4: [segments.rs:107](../../../../src/repeat_catalog/segments.rs#L107) — Every locus and cluster of repeats is deep-copied, and the original dropped

- **Confidence:** High on the count; no wall-time change could be measured
- **Hot-path evidence:** allocations in a tomato typing pass 1,933,161 → 984,799 with this and the tally part of H2; on GRCh38 the copy alone is 2,114,944 allocations (one per locus and per cluster).
- **Mechanism:** `segments_of_contig_in` consumes its regions but clones each one's kind. A locus, cluster or satellite is emitted at most once, so it can be moved; only ordinary stretches are emitted once per interval, and copying those allocates nothing.
- **Complexity cost:** none — a reordered loop. Not in the combined branch: it overlaps H2's edits to the same function and needs re-applying on top of them. Patch: the `segments.rs` half of `allocations_move_kind_and_contig_spans.patch`.

#### L5: [segment_criteria.rs:145](../../../../src/region_typing/segment_criteria.rs#L145), allocated at [segments.rs:428](../../../../src/repeat_catalog/segments.rs#L428) — Each STR locus owns its own copy of its contig's name

- **Confidence:** High on counts and memory; Medium on wall time
- **Hot-path evidence:** 1,983,306 separate name allocations on GRCh38 (2,580 distinct names), 462,678 on tomato. Sharing one name per contig (`Arc<str>`): allocations in a tomato pass 984,799 → 521,511; memory held by the run's typed regions after they are built, mimalloc, 3 runs each: tomato 241 → 199 MB, GRCh38 395 → 380 MB. Output identical.
- **Mechanism:** `finish_from_row` makes a fresh boxed string for every locus. With `Arc<str>` the name is made once per contig and every locus holds a counted reference; the field stays 16 bytes.
- **Complexity cost:** `SsrSegment::new` takes `Arc<str>` instead of `Box<str>` — a public signature change, though only three call sites needed editing because `"chr1".into()` converts to either. Two sub-agents wrote this patch independently (`allocations_contig_name_arc.patch`, `data_layout_shared_contig_name.patch`); keep one. Not in the combined branch.

#### L6: [run/segments.rs:76](../../../../src/run/segments.rs#L76) — Building the run's list of typed regions briefly needs twice the memory it then holds

- **Confidence:** High on the mechanism; Medium on the size of the peak, which varied between runs
- **Hot-path evidence:** GRCh38: 4,258,530 regions at 80 bytes each, held in a list with room for 8,388,608. Resident memory after the build 395 MB; peak while building it 781 / 989 / 873 MB over three runs, against 391 MB for the same pass when nothing is kept (data_layout sub-agent, mimalloc).
- **Mechanism:** the segment stream gives no length hint, so the list doubles as it grows; the last doubling copies 335 MB into a new 671 MB block while the old one is still live. The unused half is never touched and so never resident; the cost is the peak, which can set a one-sample run's maximum memory. Keeping one list per contig (the iterator already produces them) avoids the copy, but `Segmentation::segments()` would then return a two-level view, and its five consumers index the flat list.
- **Complexity cost:** moderate (an indexing adaptor for the consumers). `shrink_to_fit` after the build does not help: it halves the reserved figure and leaves resident and peak memory unchanged — do not merge it on the strength of an allocation counter.

### Speculative

#### S1: [loci.rs:827](../../../../src/parameter_estimation/joint/loci.rs#L827) — `generate-psps` and `estimate-parameters` read every row of the catalog three times

- **Confidence:** Low as a gain worth its complexity
- **Evidence:** after all fixes, the three passes of a `generate-psps` run on GRCh38 take about 12 s, of which reading rows is 75% (fixed-build profile). Fusing the passes per contig — read a contig once, feed the typing and the locus sampler — would save roughly two reads, about 5 s a run on GRCh38 and well under a second on tomato. The second `genome_segments` pass in `select_kept_loci` keeps only ordinary ground, which the run's own segmentation already holds, so passing that in would also drop one pass.
- **Complexity cost:** an API that hands one contig's rows to several consumers, touching `CensusPlan::of_run`, `select_kept_loci` and the run-ground caller. Worth doing only if the per-sample `generate-psps` setup still matters after H1–H3 and L1.

#### S2: [segment_criteria.rs:145](../../../../src/region_typing/segment_criteria.rs#L145) — A contig id instead of a contig name would make every held region 64 bytes instead of 80

- **Confidence:** Low on whether it is worth it
- **Evidence:** measured sizes: `TypedRegion` 80 bytes, `SsrSegment` 48, of which the name pointer is 16. With a `ContigId`, `SsrSegment` is 32 and `TypedRegion` 64: 341 → 273 MB held on GRCh38, 78 → 63 MB on tomato (arithmetic on measured sizes).
- **Complexity cost:** changes a public type. Two consumers hash the name to seed locus selection (`strata.rs:127`, `locus_generation/ssr.rs:297`) and would need the contig-name table passed in to keep their output identical; about 25 test sites change.

#### S3: [region_typing/mod.rs:1075](../../../../src/region_typing/mod.rs#L1075) — Half the held regions are ordinary stretches the other half imply

- **Confidence:** Low
- **Evidence:** 2,130,555 of 4,258,530 regions on GRCh38 are ordinary stretches (170 MB at 80 bytes). They could be generated from the gaps between the repeat features on iteration.
- **Complexity cost:** high — every consumer of `Segmentation::segments()` would move to an iterator, and the BED-edge clipping rule would be restated in a second place, where a wrong edge silently breaks the partition. Only if a memory budget for one-sample GRCh38 runs asks for it.

### Note

- **The read path's doc comment says the page index skips pages; it does not** ([parquet_file.rs:441](../../../../src/repeat_catalog/parquet_file.rs#L441)). The reader passes only a row filter, which runs after the pages are decoded, and no page selection. For a BED that touches a small part of a contig, a selection built from the offset index would read only the pages the intervals touch. After H2 the BED cases are fast enough that this is a correction of the comment first and an optimisation second.
- **Catalog build** ([builder.rs:199](../../../../src/repeat_catalog/builder.rs#L199)): still 34–41 s for tomato at 1 thread and 16–23 s at 4 threads on the loaded machine (30 s and 13 s on 2026-08-10); the rebuilt file is byte-identical to the stored one at 1 and 4 threads. It waits for the slowest contig of each batch and stops reading the FASTA while a batch is scanned; a bounded pipeline would fix both at the same memory. It runs once per reference, so this is not worth doing now.
- **Small per-contig vectors in typing** (about 81 allocations per contig, 0.04% of GRCh38's total): not worth a scratch-buffer refactor.
- **`StratumCounts`' debug output follows hash-map order** ([strata.rs:15](../../../../src/repeat_catalog/strata.rs#L15)), so two identical runs print it differently. Anything that compares census outputs as text must sort first.

## 6. Out-of-scope observations

- **[ng_typed_region_dump.rs:89](../../../../examples/ng_typed_region_dump.rs#L89) — unbuffered output.** Every `println!` is one write; in the container each write crosses to the Mac's filesystem, which made the GRCh38 dump take 226.7 s instead of 31.6 s. Fix: lock standard output once and wrap it in a `BufWriter`, as `type-regions` does. Same pattern likely in other `examples/*_dump.rs`; a mechanical sweep, separate PR.
- **[combining_stage.rs:347](../../../../src/run/combining_stage.rs#L347) and [:368](../../../../src/run/combining_stage.rs#L368)** build two more per-locus lists from the segmentation for the run's lifetime (about 80 MB on GRCh38, computed from field sizes, not measured) — part of the same fixed per-run memory as L6.
- **Two dev containers were still running when this review ended:** `de77b3bf` (10 GB), a bench link a sub-agent started that was stuck swapping, and `50fe022d` (8 GB), whose origin was not established. Neither was stopped, because stopping a container needs the owner's approval. Both took memory and CPU during the section 2 timings, which were stable regardless (about 1% between rounds).

## 7. What's already good

- **Typing streams one contig at a time** ([reader.rs:397](../../../../src/repeat_catalog/reader.rs#L397)): working memory is one contig's rows (124 MB at GRCh38 chr2), not the genome's 23.6 million, which is also what makes L1 a small change.
- **The pairing of filtered intervals back to their rows is one forward walk, not a lookup per survivor** ([segments.rs:337](../../../../src/repeat_catalog/segments.rs#L340)), and the code says why.
- **The catalog differential and the property tests made it possible to replace a fixed-point loop with a sweep and prove the output identical** — the oracle that the absorption sweep is tested against already existed.

### Author response convention

Address each finding by its identifier with one of: `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`.

---

**Audit trail.** Per-category findings, patches, profiles, timing logs and scripts: `tmp/perf_review_2026-10-09_segmenter/` in the `segmenter-perf` worktree (not committed; `tmp/` is ignored). The combined experiment is committed on branch `segmenter-perf-combined`: `74b2c508` (absorption sweep, interval index, per-batch column checks) and `d4027a8e` (index decoded once, filter skipped when every row passes), both on top of `39336d54`.
