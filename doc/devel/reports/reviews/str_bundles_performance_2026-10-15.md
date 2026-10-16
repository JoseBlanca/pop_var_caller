# Performance Review: the STR bundle caller
**Date:** 2026-10-15
**Reviewer:** rust-performance-review skill (orchestrator), with three read-only category reviews
**Scope:** what calling STR bundles costs, from reading the BAM to writing the record
**Verdict:** Apply the listed wins
**Hot-path evidence:** wall time, CPU time, instructions retired and peak RSS for bundles on and off, and samply profiles of all three commands

---

## Summary

**Bundles cost 6 to 9 in 100 of a run's CPU, and nothing measurable in memory.** Two changes
remove most of that cost. Both keep the output byte for byte the same.

- **Reading the reads: one extra pass over the whole BAM.** The bundle generator opens its own
  read cursor. That cursor decodes every record of every sample, and the segment cursor and the
  pileup walk each decode them as well. This extra decode is 77 in 100 of what bundles add to
  `generate-psps` and 57 in 100 of what they add to HG002's `call-from-alignments`.
- **Calling: work that does not depend on the read is redone for every read.** Three pieces of
  it account for about 70 in 100 of bundle calling:
  - each allele is re-rendered with each slip, for every read;
  - each partial read's layout is rebuilt for every sample;
  - every allele of the observation table is parsed into its runs, and the admitted ones are
    parsed a second time.

**Vocabulary.**
- *Bundle*: two or more tandem repeats within 15 bases of each other, called as one locus over
  their *hull*, the span from the first repeat to the last.
- *ON / OFF*: the default typing, against `--max-str-bundle-len 0`, which types every bundle as a
  satellite that is not called. ON minus OFF is what calling bundles costs.
- *Partial read*: a read that holds one flank and ends inside the hull.
- *Rendering*: an allele's bases rewritten with one slip of one run.
- *Layout*: the pieces and renderings a partial read is compared against, in the order the read
  meets them.

## 1. Scope and constraints

- **Reviewed:** the STR bundle caller at `e4bd1e94` (main), worktree `bundle-perf`, read only.
  Every source file listed below is unchanged.
- **Targets (owner):** wall time and peak memory. A large wall gain bought with a little memory is
  a win. A change that makes the code clearer is accepted even for a small gain. A change that
  complicates the code is rejected unless the gain is large. **So every finding below states its
  effect on clarity.** The caller must scale from 1 to thousands of samples and from 3 to hundreds
  of reads a position.
- **Machine:** Apple M5 Pro, 18 cores, 64 GB. Native release build: fat LTO, one codegen unit,
  `panic = abort`, line tables, mimalloc, rustc 1.98.0. The dev container service was down, so
  everything ran natively.
- **Workloads.** The two simulated workloads are 4 Mb of reference holding 3,000 segments and 500
  bundles, with PCR stutter and true stutter (no parameter fit):
  - `pcr_50s_30x` and `pcr_50s_5x`: 50 samples at 30× and at 5×;
  - `pcr_1s_30x` and `pcr_1s_5x`: one sample, the same ground;
  - GIAB HG002 at 30×, `call-from-alignments` over the tandem-repeat benchmark regions, with
    2,473 bundles.
- **In-scope code:**
  - `src/region_typing/{reference_runs,run_parse}.rs`;
  - `src/locus_generation/ssr_bundle.rs`, and the bundle's path through `src/locus_generation/ssr.rs`;
  - `src/calling/allele_candidates/{ssr.rs (select_ssr_bundle), deciding_partial_reads.rs}`;
  - `src/calling/likelihood/{ssr_emission.rs, ssr.rs}`;
  - `src/calling/inference/repeat_tract_parameters.rs` (`gather_for_bundle`);
  - `src/calling/genotype_prior/seed_ssr.rs`;
  - `src/run/callers.rs` (`call_one_ssr_bundle_locus`).
- **Out of scope**, noted in §6 only: psp reading, the generic path, and the cohort merge's
  scheduling.
- **Categories reviewed:** methodology, allocations, data layout, hot loops, float portability,
  concurrency, and I/O. Each category's full findings are in
  `tmp/perf_review_2026-10-15_str-bundles/`. Owner's deviation from the skill: the category
  reviewers worked read-only in this worktree, with no git worktrees of their own and no builds.
  Every gain marked *inferred* below is therefore a proposal, not a measurement.

## 2. Verdict

**Apply the listed wins.** Every one of findings H1 to H6 rests on a profile that names its site
and has a mechanism that is clear from the code. None of them changes a floating-point result.
**Do H1 first.** It is the largest, and it is the only one that also shortens HG002's wall time
at 8 threads.

## 3. Measurement plan

Two scripts re-run every number below:

- `tmp/perf/bundle_baseline.sh [sim] [hg002] [table]` runs each cell, then prints the table.
  `TAG=<name>` writes a variant build's results beside the baseline's.
- `tmp/perf/profile_bundles.sh` takes the samply profiles and ranks them headlessly.
  - `samply_rank.py` resolves inlined frames through `atos`, ranks self and inclusive time, and
    with `--mark REGEX` gives one family's share and its ranking.
  - `profile_diff.py` subtracts an OFF profile from an ON one.

**Gate every change on instructions retired** (from `/usr/bin/time -l`), which load does not move.
Use wall time as the second check, and require `cmp`-identical psps and VCFs. The host was shared
with another session's 16-thread job during the 5× runs, so **quote the instructions for the 5×
`generate-psps` row, not its wall time.**

### Baseline

Medians of three runs. `generate-psps` was run once per cell. CPU is user plus system time.
"G instr" is billions of instructions retired. Δ is ON minus OFF.

| workload | command | threads | wall ON / OFF (Δ) s | CPU Δ s | G instr ON / OFF (Δ) | peak RSS ON / OFF MB |
|---|---|---|---|---|---|---|
| 50 samples, 30× | generate-psps | 1 | 322 / 291 (31.8) | 32.1 | 9,199 / 8,641 (558, 6.5%) | 189 / 190 |
| 50 samples, 30× | call-from-psps | 1 | 36.7 / 33.7 (3.0) | 3.0 | 669 / 594 (75, 13%) | 262 / 261 |
| 50 samples, 30× | call-from-psps | 8 | 8.6 / 7.2 (1.5) | 3.2 | 673 / 598 (75) | 288 / 284 |
| 50 samples, 5× | generate-psps | 1 | (loaded host) | — | 2,411 / 2,311 (100, 4.3%) | 183 / 183 |
| 50 samples, 5× | call-from-psps | 1 | 17.0 / 16.2 (0.8) | 0.8 | 386 / 362 (24) | 261 / 258 |
| 50 samples, 5× | call-from-psps | 8 | 3.6 / 3.0 (0.6) | 1.0 | 390 / 366 (24) | 280 / 270 |
| 1 sample, 30× | generate-psps | 1 | 6.5 / 5.9 (0.7) | 0.6 | 185 / 174 (11) | 126 / 121 |
| 1 sample, 30× | call-from-psps | 1 | 0.5 / 0.4 (0.1) | 0.1 | 11 / 9 (2) | 24 / 22 |
| 1 sample, 5× | generate-psps | 1 | 1.7 / 1.6 (0.1) | 0.1 | 50 / 48 (2) | 120 / 119 |
| HG002 30× | call-from-alignments | 1 | 33.1 / 30.5 (2.6) | 2.6 | 764 / 710 (54) | 287 / 287 |
| HG002 30× | call-from-alignments | 8 | 30.7 / 28.2 (2.6) | 2.5 | 807 / 754 (53) | 285 / 280 |

**Per-locus cost.**
- **Calling**, at 50 samples and 30×: about 6 ms per bundle locus (3.0 s over 500), or 0.12 ms
  per sample-bundle. At 5× it is about 1.6 ms per bundle locus.
- **Generation**, at 30×: about 1.3 ms per sample-bundle (32 s over 500 × 50). An STR segment
  costs about 0.6 ms per sample in the same profile, in fetching its reads and delimiting them.
  So a bundle costs twice a segment to generate, and the extra BAM decode (H1) is the difference.
- **HG002**: about 1 ms per bundle, generation and calling together.

**Scaling.**
- **With samples, the cost is linear.** Generation costs the same per sample at 1 and at 50
  samples: 11 G instructions per sample at 30×. Calling grows from 2 G at 1 sample to 75 G at
  50, about 1.5 G per sample.
- **With depth, the cost grows less than depth does.** From 5× to 30× (six times the depth),
  generation grows 5.6 times (100 G to 558 G) and calling 3 times (24 G to 75 G). Calling grows
  less because at 5× a larger share of it is per-sample work that depends on no read: laying out
  partial reads (H4) and parsing alleles (H5).
- **Thousands of samples:** nothing found grows faster than the sample count. H4 and H5 are
  per-sample and per-allele costs, so they matter more as cohorts grow, not less.

**Memory: no finding.** Peak RSS ON and OFF differ by at most 11 MB in every cell. That is within
the spread of repeated runs: the 8-thread runs vary by about 15 MB. The psp grows by 173 kB per
sample at 30× (35.58 MB against 35.41 MB, 0.5 in 100). The bundle-scoring buffers are bounded by
the hull length, the read length, and a cap of 16 kept layouts per worker.

### Where the time goes

**Calling** (`call-from-psps`, 50 samples, 30×, 1 thread; samply, 73,350 samples). Bundle work —
any stack holding a frame whose name contains `bundle` — is 8.3 in 100 of the thread, 6,082
samples, about 3.0 s. That agrees with the ON − OFF time. Within the bundle work:

| inclusive share of bundle calling | 30× | 5× |
|---|---|---|
| partial reads: `bundle_censored_emission` | 47% | 42% |
| — the walk, `walk_from_the_held_flank` | 31% | 14% |
| — rebuilding layouts, `CensoredLayout::lay_out` | 14% | 27% |
| complete reads: `bundle_emission` | 26% | 28% |
| rendering slips, `place_slip_in_each_run` (both read kinds) | 21% | 31% |
| — of which `render_tract` | 18.5% | 27% |
| candidate selection, `select_ssr_bundle` (almost all `parse_into_runs`) | 23% | 21% |

The hottest single line is `ssr_emission.rs:1379`, the copy-by-copy append in `render_tract`:
11.7 in 100 of bundle calling. Copying and growing buffers add more:

- `_platform_memmove`: 10.8 in 100, self time;
- `realloc`: 4.3 in 100;
- growing a vector (`reserve` / `grow_amortized`): 9 in 100, inclusive.

**Generation** (`generate-psps`, the first 5 samples at 30×). The main thread takes 65,880 samples
ON and 59,888 OFF, so bundles add 5,992. Of those extra samples:

| where | samples | share |
|---|---|---|
| `fetch_capped_reads` (decoding the BAM) | 4,613 | 77% |
| `delimit` (the flat-gap delimiter) | 1,027 | 17% |
| everything else | 352 | 6% |

Decompression (`inflate`) shows the mechanism directly. With bundles OFF it takes 7,687 samples:
two full decodes, one by the pileup walk and one by the segment cursor. With bundles ON it takes
11,591, one more full decode.

**HG002** (`call-from-alignments`, 1 thread). Bundles add 5,252 samples:

| where | samples |
|---|---|
| generation | 4,333 |
| — `fetch_capped_reads` | 2,985 |
| — `delimit` | 1,304 |
| calling | 774 |

One sample's walk is serial, so at 8 threads the bundle cost is still 2.6 s of wall time.

**8 threads** (`call-from-psps`, 50 samples, 30×). Read off samply's per-thread CPU time:

| | OFF | ON |
|---|---|---|
| time the workers ran | 7.61 s | 9.23 s |
| CPU, all workers | 38.1 s | 42.8 s |
| busy fraction | 0.63 | 0.58 |

At OFF's busy fraction, the extra 4.7 s of CPU would add 0.9 s of wall time. The other 0.7 s is
workers waiting on a round (§5, L4).

## 4. Build / toolchain configuration

- **The build configuration is already what the checklist asks for**: fat LTO,
  `codegen-units = 1`, `panic = abort`, line-table debug info, mimalloc, and the toolchain pinned
  at 1.98.0.
- **L3 [build] — measure noodles-bgzf's `libdeflate` feature.** `Cargo.toml:151`. *Inferred.*
  - **Why:** zlib-rs inflate is the largest self-time row in the generation diff (1,159 samples)
    and 7,687 samples of the OFF thread.
  - **Change:** one feature flag. It helps every decode, bundles or not.
  - **Caveat:** it also switches BGZF compression, so `.vcf.gz` bytes change. Compare the
    decompressed VCF instead.
  - **Clarity:** neutral.

## 5. Code-level findings

Ranked by expected gain. Every finding keeps the psp and VCF bytes identical. Each is *measured*
when a profile names its site and *inferred* when it rests on code reading alone. Each category
file in `tmp/perf_review_2026-10-15_str-bundles/` has the full measurement plan and a suggested
diff.

### Hot-path

**H1. The bundle generator decodes every sample's BAM a third time.**
`src/locus_generation/ssr_bundle.rs:53-58`, `src/run/walker.rs:1708`, `src/read/input/cursor.rs:559-575`.

- **Cost (measured):**
  - `generate-psps`: 77 in 100 of the bundle cost, about 25 s of CPU over 50 samples at 30×, or
    5 in 100 of the whole command;
  - HG002: 57 in 100 of the bundle cost, about 1.5 s of wall time at 1 and at 8 threads.
- **Mechanism:**
  - `SsrBundleGenerator` wraps its own `SsrGenerator`, which opens its own per-contig
    `SampleCursor`.
  - The cursor reuses what it holds when a region starts at or after the last one. It then reads
    forward through every record in between, inflating and decoding each one before dropping it
    (`region_raw_aligned_reads.rs:237`).
  - Bundles are a median 5.8 kb apart on the simulated ground, so the bundle cursor reads the
    whole BAM.
- **Change:** one STR cursor for segments and bundles.
  - Segment and bundle regions come from one stream in coordinate order and use the same flank,
    so a shared cursor never has to jump. A bundle's reads would come from what the cursor
    already holds.
  - Concretely: `measure_reads` takes the delimiter as an argument; the segment generator carries
    the flat-gap delimiter and `find_reference_runs`; `GeneratorSet` routes both region kinds to
    that one generator.
- **Expected gain:** inflate returns to OFF's 7,687 samples. Bundle generation falls from 2 times
  a segment's cost to about 0.4 times.
- **Side effect:** the run report's per-sample "reads kept / dropped" lines change. Today a read
  inside both a segment's span and a bundle's span is counted once per cursor. The psp and VCF
  bytes do not carry these counts.
- **Clarity: more complex, moderately.** `walker.rs` loses a generator and a reference factory,
  but one generator in `locus_generation` would answer two region kinds, a new shape there.
- **The cheaper alternative**, L2, recovers little on the simulated ground.

**H2. `render_tract` builds each rendering one motif copy at a time, in a fresh vector.**
`ssr_emission.rs:1371-1387`, `place_slip_in_each_run` at `:1261`.

- **Cost (measured):** line 1379 alone is 11.7 in 100 of bundle calling, with `memmove` and
  `realloc` on top.
- **Mechanism:** every copy of a 1- or 2-base motif is a separate `extend_from_slice` into a
  vector that starts empty and grows as it goes. Each placement is a new `Vec<u8>` inside a
  `Vec<Vec<u8>>`, which frees the inner vectors when cleared.
- **Change:**
  1. Size the output up front and write each run with one `extend` (a repeated slice).
  2. Hold the renderings in one flat buffer, with a range per rendering, the form
     `CensoredLayout` already uses.
- **Expected gain (inferred):** most of the 18.5 in 100 under `render_tract` at 30×, and 27 in
  100 at 5×. It also helps STR segments, which use the same function.
- **Clarity: simpler.**

**H3. `parse_into_runs` fills its table against its storage order, and allocates a vector per
template state.** `run_parse.rs:247-298`, `build_template` at `:144`.

- **Cost (measured):** 23 in 100 of bundle calling. The hot lines are `:291`, `:259`, `:251` and
  `:266`. The tuple-comparison code (`lt` and `__chaining_lt`) is another 4.7 in 100.
- **Mechanism:**
  - Cells are stored `state × columns + column`, but the fill runs column outer and state inner.
    Each access therefore lands a whole row away from the last.
  - The template allocates one vector per state: about 140 for a 60-base allele with two runs.
  - Both tables are pre-filled and then overwritten.
- **Change:**
  - store the cells column-major;
  - pack the `(edits, distance)` cost into one `u64`;
  - flatten the template's edges into one list;
  - keep only two columns of the cost table;
  - reuse all of it from a scratch across a locus's parses.

  The arithmetic is integer only, so the parse is identical.
- **Clarity: neutral.**

**H4. A partial read's layout is rebuilt for every sample.** `CensoredLayoutInputs::are` at
`ssr_emission.rs:641-655`, `lay_out` at `:2094`.

- **Cost (measured):** `lay_out` is 14 in 100 of bundle calling at 30× and 27 in 100 at 5×.
- **Mechanism:**
  - A kept layout is reused only when the run stutter models match too.
  - The benchmark has 50 libraries, one per sample, each with its own stutter level
    (`libraries.tsv`). So every sample lays every candidate out again, from both flanks.
  - The models enter only through each change's probability. The renderings are the same for
    every sample.
- **Change:** cache the renderings by allele and flank, and keep the probabilities per read
  group. The walk must keep skipping changes whose probability is zero, so the sums stay
  byte-identical.
- **Expected gain (inferred):** most of `lay_out`. It grows with the number of libraries in the
  cohort and matters most at low depth.
- **Clarity: more complex** (two levels of cache). Combined with H5 it can become *the candidate
  laid out once per locus*, which serves complete and partial reads alike, and that combination
  may come out neutral.

**H5. A complete read's emission redoes, for every read, everything that depends only on the
candidate.** `bundle_emission`, `ssr_emission.rs:1458-1581`.

- **Cost (measured):** 26 in 100 of bundle calling.
- **Mechanism:** none of the following depends on the read except through its length change, yet
  all of it is redone for every read and candidate:
  - each run's split into copies and interruptions (`segment_tract_into`, one vector per
    interruption);
  - the slip ranges (`run_slip_ranges`, 1 + R vectors);
  - the slip sets for the length change;
  - each run's renderings, and a `to_vec` of each run that did not slip.

  Only the letter comparison needs the read.
- **Change, in two steps:**
  1. Reuse the buffers, and store an interruption as a range instead of an owned vector.
     Clarity: neutral.
  2. Only if `bundle_emission` is still prominent after step 1: keep a laid-out candidate per
     (candidate, models, length change), and score each candidate's complete reads together, as
     `fill_emissions` already does for partial reads. Clarity: more complex.
- **Expected gain (inferred):** grows with depth. Large at 30×, little at 3×.

**H6. Each admitted allele is parsed twice.** `allele_candidates/ssr.rs:1810-1814` against
`:1773-1777`.

- **Cost (measured as part of H3's 23 in 100):** the split between the two parses is not
  measured.
- **Mechanism:**
  - Admission only picks alleles of the merge's table, whose parses already exist.
  - `key_of` (`:1785`) also allocates a `Vec<u32>` on every call.
  - `vcf/assemble.rs:197` parses the same alleles a third time. That file is out of scope.
- **Change:**
  - take the admitted alleles' parses through the allele remap instead of parsing again;
  - number each allele's rung key once.

  Before implementing, confirm that admission never yields an allele from outside the table.
- **Clarity: simpler.**

### Likely

**L1. Two logarithms per rendering in `substitution_probability`.** `ssr_emission.rs:1087-1114`.
*Inferred.*

- **Mechanism:** it re-derives the per-base scores on every call.
- **Change:** look the scores up once per call, and compare through `MatchedPrefixLetters`. That
  struct's own test pins it bit for bit to `letters_under`.
- The same pattern appears for each rate of each partial read: compare the rate's bits before
  calling `scores_at`.
- **Clarity:** neutral.

**L2. A cursor reads a long forward gap record by record.** `cursor.rs:559`. *Inferred.*

- **Change:** jump through the index when the next region starts more than, say, 64 kb past the
  last.
- **Gain:** little on the simulated ground, where the median gap between bundles is 5.8 kb and
  nine gaps in ten are under 17 kb. It helps HG002 and sparse `--regions` runs.
- **Use:** it is the cheap alternative to H1, or a complement to it on sparse ground.
- **Clarity: slightly more complex.** It adds a speed-only constant to a rule documented as
  having nothing to tune.

**L3.** The `libdeflate` feature; see §4.

**L4. One building region per worker per round, so a region holding a bundle idles the other
workers.** `src/cli/calling_run.rs:456-458`, `src/run/cohort_merge/rounds.rs:118,149`.

- **Cost (measured):** at 8 threads, 0.7 s of the 1.5 s bundle wall cost is idling (§3).
- **Mechanism (inferred):** a round has one region per worker, and a round ends when its slowest
  region ends.
- **Change:** allow 2 to 4 regions per worker over the same ground, so rayon has smaller items to
  balance. Memory stays flat.
- **Scope:** this is the merge's scheduling, not bundle code. It would help any expensive locus.
  The in-scope remedy is to make bundle loci cheaper (H2 to H6).
- **Clarity:** neutral, but one sentence of the `round_shape_for` doc has to change.

**L5. The walk's inner loop.** `walk_from_the_held_flank`, `ssr_emission.rs:2250-2331`; 31 in
100 of bundle calling at 30×. *Inferred.* Three changes:

- keep each change's per-rendering weight in the layout, instead of dividing at every state;
- scan only the offsets actually reached;
- precompute how much each rendering shares with the one before, instead of comparing their
  bytes at `:2309`.

Each keeps the floating-point operations in the same order. Clarity: neutral to slightly more
complex.

**L6. `collect_slip_sets` recomputes its bounds at every recursion node.** `ssr_emission.rs:1799-1841`;
about 6 in 100 of bundle calling. *Inferred.*

- **Change:** hoist the per-depth bounds out of the recursion.
- **Constraint:** keep the left-to-right product, so the 10⁻¹² pruning decision cannot move.
- **Clarity:** neutral.

### Speculative

**S1. Check the slot that matched last before scanning the 16 kept layouts**
(`ssr_emission.rs:1940`). It is small, and H4 changes the cache anyway.

### Note

- **[bench]** No criterion bench covers bundle generation or bundle calling. The two scripts in
  §3, gated on instructions retired, take its place for this round.
- **[profile]** The cursor's counters (`regions_jumping`, `regions_reusing`, `reads_evicted`) are
  printed nowhere, and none counts the records read past without being used. A `records_read`
  counter would make H1 and L2 checkable on a count that is the same on every run.
- **[profile]** `pvc` cannot run under DHAT: `src/main.rs:19-21` installs mimalloc
  unconditionally. About six lines behind the existing `dhat-heap` feature would allow allocation
  counts for H2, H3 and H5.

## 6. Out-of-scope observations

- **`_platform_memmove` is 19 in 100 of `call-from-psps`'s self time** at 50 samples, 30×, one
  thread — more than all bundle calling (8 in 100). Its callers are the psp read path's
  `extend_from_slice` under `resolved_into` (20 in 100 inclusive) and `rebuild` (17 in 100); the
  chain-id replay (`psp::chain_ids::apply_arrivals`) is 15 in 100 inclusive. The attribution of
  memmove between them was not split. A separate review.
- **The BAM path shares no decoded blocks between a sample's cursors**, so with bundles OFF the
  pileup walk and the segment cursor still decode the BAM twice. The CRAM path has a shared
  cache. H1 is the bundle-sized part of this.
- **At 8 threads, HG002's `call-from-alignments` takes almost the same wall time as at 1 thread**
  (30.7 s against 33.1 s). One sample's walk is serial.
- **`rounds.rs:150`'s `map_init`** builds a new calling scratch for every rayon split in every
  round, so worker scratch is not reused across rounds.
- **A stale banner at `callers.rs:7238-7243`** calls the round driver a prototype that is reached
  only under `NG_ROUND_CALL`. It is in fact selected at every thread count above one.

## 7. What's already good

- **Partial-read layouts are cached and compared by their inputs, never by address**
  (`CensoredLayoutInputs`). A layout can only be reused for the allele it was built from. That
  is what makes H4 a change of key rather than a redesign.
- **`fill_emissions` scores candidates outer and reads inner for partial reads**
  (`likelihood/ssr.rs:864-909`), so the per-candidate work can be shared. H5 extends that order
  to complete reads.
- **The 10⁻¹² budget in the forward walk has a stated bound**: the result is at most 10⁻¹² low,
  never high. It is the reason the walk is 31 in 100 of bundle calling and not most of it.

### Author response convention

Answer each finding by its code: `applied in <commit>` / `experiment shows no gain — closing` /
`disputed because …` / `deferred to …` / `won't fix because …`.
