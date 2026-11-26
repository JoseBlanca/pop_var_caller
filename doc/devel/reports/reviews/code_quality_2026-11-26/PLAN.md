# Code-quality program — the work plan (2026-11-26)

This plan turns the eleven area reports beside it into work the next agents can pick up. It does
not restate the reports: each item names the report section that holds the evidence, and a fixer
reads that section before starting. Every bug citation below was re-read against the worktree at
`f84129cf`; one stale range is corrected where it occurs (item 1 of the next section).

## How to use this plan

**The owner's rules, as they apply here.**

- Nothing is an absolute rule. Every item states its downside; "leave it" is a valid outcome,
  and the last section lists what was examined and left.
- The program acts on its own (**auto**) when a change alters no output, slows no step and
  raises no memory. Anything else goes to the owner with its evidence (**owner**, list B).
- **Every refactor commit passes the gate**: `scripts/code_quality_oracle.sh` (being built)
  checks that every catalog, psp, parameters file, VCF and window-coverage file is
  byte-identical to the recorded baseline, and that no step got slower or larger beyond the
  measured noise. Its data are tomato (16, 63 and 1 accessions) and GIAB HG002 at about 290
  reads a position, at `--threads 4`.
- **If the gate flags a package as slower or larger**, the package stops there and the item
  becomes an owner decision, with the gate's figures as its evidence.
- **What the gate cannot see**: log text, error messages, `--help`, and inputs outside its data
  (CRAM, long reads, hundreds of read groups, scaffold-level references, thousands of samples).
  A change touching one of those needs its own test, or a before/after diff in the commit
  message (for `--help`).
- **Bug fixes are separate commits from refactors.** Each bug commit adds a test that fails
  before the fix and passes after. A performance fix with identical output is a refactor: it
  passes the gate instead, and its commit quotes the gate's timing line.

**Report key.** Items cite reports by these short names and the report's own section number.

| Short name | Report | Area |
|---|---|---|
| cand | [calling_likelihood_candidates.md](calling_likelihood_candidates.md) | read likelihood, candidate selection, genotype prior, site quality |
| inf | [calling_inference_core.md](calling_inference_core.md) | `calling/inference/` and the top of `calling/` |
| pfile | [calling_parameters_file.md](calling_parameters_file.md) | `calling/parameters_file/` |
| fits | [parameter_estimation_fits.md](parameter_estimation_fits.md) | the SNP/indel and STR fits |
| est | [parameter_estimation_rest.md](parameter_estimation_rest.md) | census, contamination, the rest of `parameter_estimation/` |
| run | [run_orchestration.md](run_orchestration.md) | `run/` drivers |
| merge | [run_cohort_merge_paralog.md](run_cohort_merge_paralog.md) | cohort merge, duplication filter |
| locus | [locus_generation.md](locus_generation.md) | `locus_generation/` |
| io | [io_psp_read_alignment.md](io_psp_read_alignment.md) | `psp`, `read`, `bam`, `alignment`, `read_phase` |
| cli | [cli_and_small_modules.md](cli_and_small_modules.md) | `cli/` and the small modules |
| cross | [cross_module.md](cross_module.md) | dependency graph, cross-module duplication |

---

## 1. Findings several reports made independently

Two or more reviewers reached these without seeing each other's work, so they carry more
confidence than a single report's finding. Where the reports disagree, the plan's choice is
stated.

| # | Finding | Reports | Where they disagree, and the plan's choice |
|---|---|---|---|
| 1 | **The cohort merge's output types live in `run`, so the calling library imports the driver.** `CohortObservation` and five types it holds (`run/cohort_merge/build.rs:1289-1810`) and the `MinAltObs`/`MinAltReadShare`/`MinAltReads` thresholds (`run/cohort_merge/mod.rs:574-790`) are imported by 9 calling files. | cand A1, inf 2.2, merge A1, cross A1 — **four** | All four: move the definitions to a new top-level `src/cohort_observation.rs`, keep the constructors in `build.rs` as an inherent `impl`, move the thresholds with them. inf also offered `locus_generation/`; cand argues against it. **Plan: top-level.** merge A1 cites the thresholds at `mod.rs:117-383`; the right range is `:574-790`. |
| 2 | **`calling/likelihood/run_context.rs` is the only production edge from `parameter_estimation` into `calling`**, and so the cause of their cycle. It imports only `float` and `types`. | inf 2.1, est A1, cross A4 — three | inf and est: a top-level module. cross: top-level or under `parameter_estimation`. **Plan: top-level, under a name that does not collide with `crate::run`** (for example `misread_context`). |
| 3 | **Library code imports the progress-log target and two formatters from `cli::run_output`** (`run/calling_progress.rs:23`, `parameter_estimation/progress.rs:25`), and `neighbour_evidence/peak.rs:408` copies `with_thousands` to avoid the import. Size, duration and `/proc/self/status` helpers are duplicated beside them. | est A2 + D1, run A2, cli A1, cross A2 + D2 — four | All agree on a small crate-root module that `cli::run_output` re-exports. est and cross both say keep the two duration styles (`350ms` against `0s`) apart. **Plan: as agreed; merging log styles is owner item B23.** |
| 4 | **`ln Σ exp` over a slice is written three times, identically** (`calling/quality/mod.rs:856`, `joint/contamination.rs:1404`, `joint/fit.rs:1586`). | cand A3, fits 12, est D2, cross D1 — four | All: one `float::ln_sum_exp`. cand also moves `log_sum_exp_2` out of the prior's model file; cross says it is already shared and can stay. **Plan: move it too** (site quality then stops importing the prior). All four: **never** merge `paralog/locus_score.rs:700` (`log_add_exp`), which uses a faster series and would change output. |
| 5 | **About 2,600–3,500 lines of `src/bam` are dead**: an older reader stack hidden by module-wide `allow(dead_code)`. | io A1, cross X1 — two | Agree. cross adds: check that no CRAM case is covered only by the dead files' tests before deleting. |
| 6 | **`call-from-alignments` computes a region count its driver never reads**, so the region width depends on `--threads` and the number of covers grows with it. | run A1, merge B1 — two; cli A2 and cli's left-alone list touch it | run says memory would not change with a full-width region; merge says it would rise at high thread counts. cli's left-alone list calls the thread dependence harmless. **Needs a measurement: owner item B14.** |
| 7 | **`window_coverage::depth` and `run::cohort_merge::observation_cache` import each other.** | cli A3, cross A3 — two | cli: move `depth.rs` under `run/`. cross: move `LocusSummary` down and the `From<&Drawn>` impl into `observation_cache.rs`. **Plan: cross's**, which keeps the per-record depth rule beside the accumulator its spec describes. Keep `LocusSummary`'s measured 32-byte layout. |
| 8 | **`read` imports `psp` for one function** (`ReadFilterConfig::provenance_parameters`, `read/filtering.rs:124`). | io A3, cross A5 — two | cross: move it to `psp` or `run`. io: `run`, since both the writer and the reader of those keys live there. **Plan: `run`.** |
| 9 | **`vcf` and the two refinement stages import each other.** | io A4, cross A6, cli A4 — three | cross and cli: leave it. io: move `read_phase/stage.rs` into `vcf/` (cheap; `read_phase` then imports nothing above it). **Plan: do io's move, leave `neighbour_evidence`.** |
| 10 | **Unused or test-only public items in `calling`**: `genotype_log_likelihood_row` (generic), `fill_contaminant_allele_frequencies`, two test-only re-exports, `charges_of`, `with_frequencies`. | cand R2 + T1, inf 4.6, cross X4 — three | Agree: `#[cfg(test)]` and drop the re-exports. |
| 11 | **A module-wide dead-code silencer covers live STR selection** (`allele_candidates/ssr.rs:35`), and its comment calls the module dead. | cross X3, cand R1 (the "public because the measurement calls it" note) | Agree. |
| 12 | **Comments cite the caller deleted on 2026-09-13** (`src/var_calling/…`, `src/ssr/…`, "production's"), deleted examples, and plan step codes. | every report | All say: fix false statements and dead links; do **not** sweep history out wholesale. **History trimming: owner item B21.** |
| 13 | **Empty and caller-less modules**: `src/iter_ext/` (5 lines of doc), the third region tiler in `tandem_repeat.rs:876-1027`. | cli X1 + X2, cross X2 + X4 + R1 | Agree. |
| 14 | **STR fit test-only functions** `fit_multiplier` and `scale_to_a_typical_read_group` (`read_group_stutter.rs:94`, `:235`). | fits 19, cross X4 | Agree. |
| 15 | **MD5-to-hex written seven times.** | cli B2 (proposal part), cross D3 | Agree on one encoder; leave the three decoders apart (cross D3 gives why). |
| 16 | **`ContigList` → `ContigBounds` written six times.** | cross D4 + B1, cli D5 | cli: leave `contig_bounds` (each carries its own error type). cross: one shared function whose error each caller maps. **Plan: cross's**, because it also fixes A38's wrong message. |
| 17 | **Global thread-pool setup in six commands.** | cli D5 + B5, cross D7 | Both lean "leave". **Plan: leave**; the `--threads 0` meaning on `repeat-catalog` is owner item B6. |
| 18 | **Phred → `ln P(error)` written three times.** | locus D4, cross D8 | Both: the `ssr.rs:1869` copy rounds differently and must stay. cross: share the other two. **Plan: share the two.** |
| 19 | **Several SplitMix-like generators.** | locus D4, cli D5 + left-alone | All: leave every one; each feeds a pinned stream. Fix the `ssr.rs:318-320` comment that calls its copy "ported verbatim". |
| 20 | **Two golden-section searches** (`fit.rs:5214`, `contamination.rs:1303`). | fits table + left-alone, est D3 | Both: leave; they return different points. |
| 21 | **Two "round drivers" with similar names.** `rounds::merge_and_call_in_rounds` is a test oracle compiled into release (merge A3, cross left-alone); `callers::call_cohort_from_sources_in_rounds` is psp mode's production driver whose comment calls it a prototype (run D1, R1). | merge, run, cross | Not a disagreement, a naming trap. Read both reports before touching either. |
| 22 | **Test-only cycles**: `fitted_inbreeding_of` called from `parameter_estimation` tests (run A3, est A3); `parameters_file` tests calling `cli` (pfile 4). | run, est, pfile | Agree: cheap to cut. |

---

## 2A. Bugs

Ranked by consequence to calls and runs. **Route** says who decides: *auto* fixes go ahead; *owner*
fixes wait for the decision in list B. "Valid input" means input the caller is meant to accept.

### Calls silently wrong, or data silently lost

| # | Bug | Where | Consequence | Conf. | Route | Test that fails before the fix |
|---|---|---|---|---|---|---|
| A1 | A direct-mode run joins four repeat-tract tables on the parameters file's read-group numbers, not its own (pfile 1). | `calling/parameters_file/to_run_parameters.rs:374-380, 512-529, 541-545, 641-645` | Whenever the BAMs' read-group order differs from name order (one sample, two lanes is enough), one lane's stutter multipliers, substitution rates and misread table score another lane's reads. Wrong STR likelihoods, no message. | high, traced | owner **B1** | Extend `the_same_samples_in_another_order_are_the_same_cohort` (`bindings.rs:1416`) to assert a substitution rate, a stutter multiplier and a misread entry on a moved lane. |
| A2 | A BAM cut at a compressed-block boundary reads as a shorter valid file (io B1). | `read/input/aligned_reads_reader/bgzf_blocks.rs:270-275`, `bam.rs:233-241`; CRAM: `container.rs:552-570` | Every read after the cut is missing; calls there become no-calls or reference calls; the run finishes. | high, traced | auto (damaged input only): an index chunk starting past the file's end, or a block header cut short, is an error; a missing end-of-file marker is a warning. Refusing every file without the marker would be owner. | Truncate a fixture BAM exactly at a block boundary; reading a region past the cut must fail. |
| A3 | Read-group rows written out of order in a parameters file swap calibration, contamination and batch between lanes (pfile 2). | `to_run_parameters.rs:165-178` (the comment at `:176-178` states the opposite) | Wrong base-quality calibration and contamination per lane. Only hand-written or hand-edited files; the spec invites both. | high, traced | auto | Extend `a_file_whose_rows_are_written_in_another_order_is_the_same_file` (`bindings.rs:1648`) to project both files and compare. |
| A4 | A stale index is used without warning, and `.csi` always wins over a newer `.bai` (io B5). | `bam/index_preflight.rs:306-317` | Usually a loud decode failure; possibly silent loss of reads at a contig's tail (suspected). | medium | auto: warn when the index is older than its data; take the newer of `.csi`/`.bai`. | Fixture with an older `.csi` and a newer `.bai`: the `.bai` is chosen. |
| A5 | CRAM positioning can skip an earlier container whose long read reaches the region (io B2). | `read/input/aligned_reads_reader/cram.rs:344-363` | Long or spliced reads missed on a jump; CRAM and BAM of the same data can disagree, and the loss can depend on worker boundaries. | medium, suspected | owner **B9** | Unit test on `first_entry_reaching` with hand-built entries (one at 1,000 spanning 50,000; one at 1,100–1,400; start 20,000): expect index 0. |
| A6 | With `--max-str-len` above `--max-str-bundle-len`, a `--regions` run can call a satellite as an STR segment (cli B1). | `repeat_catalog/mod.rs:459-479`, `reader.rs:213-217` | A region run and a whole-genome run classify the same locus differently. Non-default settings only. | medium, reasoned | owner **B7** | Differential: caps 400/100, a 419-bp run of abutting tracts, region read against whole read. |
| A7 | A census drops a tract's "reached but crossed no read" count until a read group has an observation (est B3). | `parameter_estimation/joint/census.rs:2893-2907` | A handful of reads per sample missing from a stratum's reaching total; in rare cases changes which tracts the duplicate filter excludes. | medium | owner **B4** | Feed a `CensusWriter` a first kept tract with only `reads_without_observation`; assert `covering_not_crossing`. |
| A8 | At about 700× and above, reads touching only a flank use up the STR read cap (locus B7). | `locus_generation/ssr.rs:454-494`, `:985-992` | Fewer informative reads than the cap allows at very deep loci. | medium | owner **B8** | 1,500 reads at one tract, 300 flank-only, cap 1,000: assert 1,000 informative reads kept. |
| A9 | A leading soft clip lets an edge indel survive that the unclipped read drops (io B6). | `alignment/indel_norm.rs:335-350`; `locus_generation/pileup/cigar_cursor.rs:339-340`; `pileup/decompose.rs:85-87` | An indel candidate at an un-normalised anchor, its support split. Frequency unmeasured. | medium mechanism, low size | owner **B10** | `5S 3M 2D 40M` in a homopolymer against `3M 2D 40M`: the same events. The existing test `soft_clip_then_indel_away_from_chrom_start_keeps_indel` pins today's behaviour and would change. |
| A10 | The SNP/indel fit re-fits homozygote excess at quadrature points the pass did not count at (fits 1). | `parameter_estimation/joint/fit.rs:5074-5086` moves the shapes; `:5114-5128` builds the rule from the moved ones | Converged answer is still the maximum; the iteration may take extra passes. Possibly the source of the drift `JUMP_SLACK` tolerates. Output changes in the last digits. | high mechanism | owner **B2** | Unit test: the excess step equals an M-step computed against the rule the counts were gathered at. |
| A11 | STR fit: after a SQUAREM jump the next round reuses pre-jump counts; a fit ending on a jump returns numbers no round checked (fits 2, 3). | `ssr_fit/penalised_fit.rs:1981-1987`, `:1797-1822`, `:1958`, `:2239` | Reported numbers and their curves/standard errors can disagree; small output changes. | medium, suspected | owner **B3** | A fit with `max_rounds` set so the last round jumps: reported numbers equal the state the curves describe. |
| A12 | `GQ` keeps the full-table value after `GP` is renormalised over the kept alleles (cli B3). | `vcf/assemble.rs:480-628` | `GQ` lower than `GP` implies, on about 9 records in 2,290 (GIAB, 30×). | medium | owner **B11** | If recomputed: a record that drops an alternative has `GQ` from its written `GP`. |

### The run stops, or refuses input it should accept

| # | Bug | Where | Consequence | Conf. | Route | Test that fails before the fix |
|---|---|---|---|---|---|---|
| A13 | One read with a deletion longer than `max_record_span` (5,000) stops the run (locus B1). | `pileup/open_record.rs:1788-1795`, `:1706-1713`; message `pileup/errors.rs:71-75` | Long-read data stops mid-genome. | high that it fires | owner **B5** | Walk a read `100M6000D100M`: today `RecordTooWide`; after, the read is dropped and counted. |
| A14 | A read with bases and `QUAL *` stops the run (io B3). | `read/aligned_read.rs:159`; `read/prepared_read.rs:350-356`; `pileup/genome_walk.rs:1024` | BAMs made from FASTA input or some assemblies abort the sample. | medium-high | owner **B5** | A record with `SEQ` and `QUAL *`: the run completes, one read counted as dropped. |
| A15 | A census whose header and directory exceed 1 MiB is refused as "not a census" (est B1). | `parameter_estimation/joint/census_file.rs:454-502`, `:511` | Above about 190 read groups at tomato's 141 strata, `estimate-parameters` stops. | high, traced | auto | Write and reopen a census of 200 read groups × 141 strata. |
| A16 | `sample.bai` (extension replaced) is not found (io B4). | `bam/index_preflight.rs:306-317` | Picard/GATK indexes reported missing, or a redundant `.csi` built (which then feeds A4). | high | auto | Fixture with `x.bai` beside `x.bam`: found. |
| A17 | `--max-locus-span` above 65,535 can panic on a long generic locus (run B2). | `run/records.rs:358-361` | Panic mid-genome at a user-set value. | medium, suspected | owner **B6** | Parse `--max-locus-span 70000`: refused. |
| A18 | `--max-str-len` / `--max-str-bundle-len` above 65,535 saturate `LocusLen` (locus B2). | `locus_generation/witness.rs:324-330`; `ssr.rs:1931-1932` | Depth along a very long tract short at its far end. | medium, suspected | owner **B6** | Parse either option at 70,000: refused. |
| A19 | A malformed catalog footer panics instead of returning `Unreadable` (cli B2). | `repeat_catalog/parquet_file.rs:448-457` | Corrupt file → crash. | high | auto | `unhex` of a 32-byte string holding a multi-byte character returns `None`. |
| A20 | A prefetch-thread panic leaves waiters blocked for ever (io B8). | `read/input/aligned_reads_reader/prefetch.rs:140`, `container.rs:512-523` | A hang instead of an error. | medium | auto | Inject a decoder panic; the waiter receives an error. |
| A21 | An interrupted index build leaves a truncated index at its final path (io B10). | `bam/index_preflight.rs:351-385` | The next run fails until the user deletes it. | high | auto | Fail the build mid-write; no file at the final path. |
| A22 | `read.pos as u32` wraps silently on a CRAM contig over 4.29 Gb (io B11). | `read/prepared_read.rs:68` | Reads placed at wrong positions on giant contigs. | high | auto | A position above `u32::MAX` is an error. |

### Latent: correct today, unguarded

| # | Bug | Where | Consequence | Conf. | Route | Test |
|---|---|---|---|---|---|---|
| A23 | The window search assumes each sample's reaches never decrease; only starts are checked (merge B3). | `run/cohort_merge/observation_cache.rs:2099-2101`, `:2005-2011` | A generator defect would silently evict a record early in direct mode. | medium | auto (one comparison a record) | A source yielding `[100–130], [105–105]` is refused. |
| A24 | A cohort mixing built and kept sources would be misindexed or panic (merge B4). | `observation_cache.rs:1605`, `:1226`, `:1392`, `:1947` | None today: no source mixes. | low reachability | auto (assert, or document the contract) | A source switching shape trips the assertion. |
| A25 | A NaN in the STR fit can neither stop the loop nor fail the jump test (fits 5). | `penalised_fit.rs:1721`, `:1752-1757`, `:1762`, `:1817`, `:2202` | The fit could "converge" on NaN. No input known to reach it. | low reachability | auto (gate must show identical output) | `state()` with a logit at the ±40 clamp stays finite. |
| A26 | The SNP/indel fit reads the census length from the first sample's first section (fits 6). | `fit.rs:1058-1063`; NaN at `:2588-2590` | A first sample with no ordinary-position section makes every rate silently `Defaulted`. Needs a hand-built census. | high mechanism | auto | Cohort whose first sample lends no ordinary section: positions are counted. |
| A27 | The misread-after-runs table skips `validate` (pfile 3). | `validate.rs:128-145`; checks live in `to_run_parameters.rs:580-679` | `estimate-parameters`' self-check does not cover it; NaN standard errors accepted. | medium | auto (refuses only bad hand-edited files) | A file with a NaN `standard_error` in a misread row: `validate` refuses. |
| A28 | An argument with a newline breaks the VCF header (cli B6). | `cli/provenance.rs:13-18`; `vcf/header.rs:683-692` | Invalid VCF for a path containing `\n`. | medium | auto | Header built from an argument with `\n`: every header line starts with `##`. |

### Speed and memory

Output is identical, so each is a refactor gated by `scripts/code_quality_oracle.sh`; the gate's
data do not reach most of these paths, so each commit also quotes a targeted timing.

| # | Problem | Where | Size | Route |
|---|---|---|---|---|
| A29 | `call-from-alignments` narrows regions by thread count and builds them one at a time (run A1, merge B1). | `cli/calling_run.rs:485-519`; `run/callers.rs:867`; compare `run/psp_caller.rs:856-875` | about 8× the cover work at 48 threads against 4 (reasoned) | owner **B14** |
| A30 | Contamination's similarity matrix and Jacobi eigen-solver grow as samples² and samples³ on one thread (est B4). | `joint/contamination.rs:998-1019`, `:1412-1469` | unknown at 2,000 samples | owner **B15** (parallel matrix is auto once measured) |
| A31 | The census directory overlap check is quadratic (est B2). | `census_file.rs:726-735` | 10⁹–10¹⁰ comparisons a sample once A15 is fixed | auto |
| A32 | The STR binomial tail can walk every read count when its first term underflows (fits 7). | `ssr_fit/tract_exclusion.rs:280-292` | ~180,000 `lgamma` calls a tract at 60,000 reads | auto |
| A33 | At the 32,768-read hold limit every arriving read rehashes every held name; deep dropped-read lists are searched linearly (locus B5). | `pileup/active_read_set.rs:526-536`; `pileup/genome_walk.rs:2076-2084` | only above ~1,000 reads a read group a position | auto (+8 bytes a held read; gate) |
| A34 | Each STR locus re-reads ~92 reference bases it just released (locus B6). | `ssr.rs:4218-4222`, `:4406-4417` | one small read a locus a sample | auto |
| A35 | The read cursor holds every read of a region until the next move, deep-copied (io B7). | `read/input/cursor.rs:685`, `:575-592` | ~1 GB a cursor for a 1 Mb segment at 300× | auto, M; first confirm segment lengths |
| A36 | Contig name → id is a linear scan per STR locus (cross B2). | `run/gatherer.rs:310`, `cli/estimate_parameters.rs:955`, `cli/estimate_contamination.rs:772` | ~1.5 × 10⁹ comparisons a sample on a 30,000-contig reference | auto |

### Wrong messages and counters

| # | Bug | Where | Route | Test |
|---|---|---|---|---|
| A37 | A record refused by the output is reported at the locus just called, not the one being written (run B1). | `run/callers.rs:872-880`, `run/psp_caller.rs:843-855` | auto (inside package R2 below) | A `hand_over` failing on a held-back record names that record's locus. |
| A38 | `estimate-contamination` reports a >4 Gb contig as "more contigs than a run can address" (cross B1). | `cli/estimate_contamination.rs:232-233`, `:602-613` | auto | Contig of 5 × 10⁹ bases: message names it. |
| A39 | A warning says "the run will fail" when it will not; two messages name a nonexistent `--max-active-reads` (locus B4). | `pileup/generator.rs:959-962` (passes a real cap, unlike `genome_walk.rs:842-858`); `chain_id_allocator.rs:458-474`; `errors.rs:46-49`, `:63-65` | auto (log/error text) | Error texts name no flag the CLI lacks. |
| A40 | Two walk counters miscount (locus B3). | `pileup/genome_walk.rs:1186-1197`, `:1309-1313`, `:1370-1372` | auto if log-only; `run/callers.rs:2099-2110` documents them, so first check whether they reach the run report (if so, owner) | Two read groups: `positions_short_of_cap` counts per group. |
| A41 | An STR fit that runs out of rounds logs a progress line, not a warning (fits 4). | `penalised_fit/run_log.rs:357-366`; `penalised_fit.rs:1225-1229` | auto (log); a `settled` flag in the file is owner **B3** | Captured log: the line is a warning. |
| A42 | The neighbour-evidence "re-called" count includes unchanged lines (cli B4). | `neighbour_evidence/line.rs:123-135`, `:192` | auto (log); dropping the empty `PS` is owner **B12** | Count over a line where no sample changed is zero. |
| A43 | `--threads 0` means one thread for `repeat-catalog`, every core elsewhere (cli B5). | `cli/repeat_catalog.rs:52-60`, `:204-207` | owner **B6** (help text) | — |
| A44 | `MergedCursors::move_to_region` is documented all-or-nothing and is not (io B9). | `read/input/sample_cursor.rs:193-196` | auto (doc; callers treat the error as fatal) | — |

### Tests that do not test what they claim

| # | Gap | Where | Route |
|---|---|---|---|
| A45 | A comment cites `the_round_driver_writes_what_the_streaming_driver_writes`, which never existed; on a one-core machine nothing compares psp mode's two drivers (run B3). | `run/psp_caller.rs:805-811` | auto: write the test (one psp cohort at one region in flight and at several, VCF bytes compared). |
| A46 | The locus-ownership rule is tested on a test-only copy (`Organiser`); production's `Frontier::claim` displacement branch never runs in a test (merge B2). | `run/cohort_merge/organise.rs:115-492`, `rounds.rs:269-324` | auto (package D4). |
| A47 | The generic read-likelihood tests mostly go through a flankless wrapper production never calls (cand T1). | `calling/likelihood/generic.rs:771` | auto (package D2). |
| A48 | A "parity with production" test compares a literal with a literal (cand T2). | `calling/likelihood/mod.rs:2932` | auto (rename or delete; package X1). |

---

## 2B. Owner decisions

Each item: what is decided, the options, the trade-off, and a recommendation. Evidence is in the
cited report section.

| # | Decision | Options and trade-off | Recommendation |
|---|---|---|---|
| B1 | Fix the parameters-file read-group join in direct mode (A1; pfile 1). | Fix changes output for direct-mode runs whose read-group order differs from name order; psp mode and same-order runs are unchanged. | **Fix.** S. Pair with A3, same code. |
| B2 | Fix the homozygote-excess M-step (A10; fits 1). | Fixed numbers move within the stopping tolerance; fixtures pinning exact values need updating. Possible gain: fewer passes, and `JUMP_SLACK` may become removable. | **Fix**, then rerun the ten-sample cohort in `JUMP_SLACK`'s comment to see if the slack can go. |
| B3 | STR fit: reset the change measures after a jump, restore the fallback when the round limit ends on a jump (A11; fits 2, 3). Also: write a per-motif-length `settled` flag to the parameters file (fits 4), and clamp the level in ordinary rounds as a jump does (fits 21). | The two fixes change output slightly and slow convergence a little. The flag is a format addition. The clamp question changes output wherever the level exceeds 0.95. | **Fix the two defects. Add the flag** when the format next changes. **Measure** how often the level passes 0.95 before deciding the clamp. |
| B4 | Seed the census writer's STR map with every declared read group (A7; est B3). | Census bytes change for affected samples; rarely, which tracts are excluded. | **Fix**; it removes a dependence on which tract a walk meets first. |
| B5 | A read-filter policy for reads that stop the walk: deletions longer than `max_record_span` (A13) and reads with no base qualities (A14). | (a) Drop each before the walk and count it under a new reason — adds tally fields and psp provenance keys, so psp headers change. (b) Drop only the offending event in the walk. (c) Leave: such inputs stop the run. | **(a)** for both. A run that stops on one read is worse than a counted drop; the error message already says an upstream filter should have removed the read. |
| B6 | Bounds and meaning of three options: refuse `--max-locus-span`, `--max-str-len`, `--max-str-bundle-len` above 65,535 (A17, A18); what `--threads 0` means on `repeat-catalog` (A43). | Refusing changes the accepted command line. For `--threads 0`: documenting "0 means 1" keeps behaviour; making it "every core" risks 64 × 1.4 GB on a large machine. | **Refuse above 65,535** (check first whether more code assumes `u16`). **Document** `--threads 0` as one thread. |
| B7 | Include `--max-str-len` in the catalog reader's padding (A6; cli B1). | Output changes at non-default settings only; region reads grow a little. | **Fix.** |
| B8 | Test the STR read's "outside the tract" condition before the read cap (A8; locus B7). | Output changes wherever the cap binds (above ~700×); `reads_discarded_by_cap` changes meaning. | **Fix**, with the counter renamed or a second counter added. |
| B9 | Store each CRAM index entry's running maximum reach (A5; io B2). | Output changes only where reads were wrongly missed; 8 bytes an index entry. | **Fix**; also add a CRAM fixture with one long read. |
| B10 | Treat "first/last operation" as the first/last non-clip operation for edge indels (A9; io B6). | Output changes; a test pins today's behaviour deliberately. | **Measure first**: count HG002 chr20 prepared reads whose first two operations are a clip then an indel. Decide on the number. |
| B11 | `GQ` after dropping an uncalled alternative (A12; cli B3). | Recompute from the renormalised `GP` (changes `GQ` on about 4 records in 1,000), or keep and document that `GQ` is what the site was judged on. | **Keep and document**; the comment at `assemble.rs:478-479` already argues for it. |
| B12 | Pass neighbour-evidence lines that no sample changed through unchanged, without an empty `PS` (cli B4). | Changes output bytes on those lines; today's `.` is valid VCF. | **Leave**; fix only the count (A42). |
| B13 | Remove the `lowDepth` filter line, which this caller never writes (cli X3). | Every VCF header loses one line. | **Remove**, batched with the next approved output change so baselines move once. |
| B14 | `call-from-alignments`' region shape (A29; run A1, merge B1). | (a) Full ground as one region for alignments: width independent of `--threads`; speed and memory unmeasured, and the two reports disagree on memory. (b) Route alignments through the round driver: L, needs the walker to be `Sync`. | **Measure (a)** on the tomato benchmark at 4 and 48 threads (time and peak memory). Adopt it if neither gets worse; fix the doc at `calling_run.rs:466-475` either way. |
| B15 | Contamination at thousands of samples (A30; est B4). | Parallel matrix: identical output. Subspace iteration for 4 eigenvectors: last-digit changes. | **Read the `contamination:` timestamps** in the 2,169-sample run's log first. If minutes: parallel matrix (auto); eigen change only if still slow. |
| B16 | Keep or delete the per-locus slippage re-fit (inf 4.4). | Off by default; measured 7 genotype changes in 6,993 tomato tracts, net +1. ~600 production lines and the nested round loop in `call_locus`. Deleting removes two documented development variables; re-measuring on another cohort would mean rebuilding it. | **Delete.** |
| B17 | Keep or delete contamination's five research options (est R5). | None reachable from a command; `fit_alpha` would shrink from 115 lines to about 25. Re-measuring means rewriting them. | **Delete**, keeping the measurements in the module doc. |
| B18 | `CensusMomentSums`: delete or finish (est R1). | Summed at every census position on every fit pass, never read; the spec `ordinary_site_prior_moments.md` asks for the averages. | **Delete** (fit gets lighter); bring it back with its reader when the report is specified. Check first that no `Debug` dump of `JointFit` reaches a user file. |
| B19 | Keep reading parameters files written before 2026-10-16 (est R7). | Removing `share_curve.rs`, most of `slippage_curve.rs` and two provenance variants changes which inputs are accepted. | **Ask**: if no such files are in use, remove on both sides. |
| B20 | `PspWriter::append` and its tests: delete, or mark as a specified, unwired feature (io A6). | ~440 lines of tests and error variants no run can raise; a future resumable run might want it. | **Delete**; record it as not implemented in spec §6.4. |
| B21 | How much history to take out of comments (cand R5, pfile 14, merge R4, locus R1, io R4). | A sweep risks deleting load-bearing reasons; leaving it keeps comment lines above code lines in many files. | **Auto**: false statements, dead links and plan step codes (`B2`, `C3b`, "Milestone E"). **Only when a file is touched for another reason**: move superseded positions and dated narratives to the cited spec or report. No sweep. |
| B22 | Add an `fsync` to the catalog writer before its rename (cross D6). | A crash just after the rename can leave a truncated catalog, which the reader refuses (error, not wrong data). One fsync per catalog write. | **Add it.** |
| B23 | Unify log text styles: durations below a second, five byte-size spellings, `peak.rs`'s decimal GB (est D1, cross D2). | Log text only; `peak.rs` is the one place printing decimal gigabytes, which may or may not match how `--max-memory` is read. | **Leave the styles**; answer the `peak.rs` question. |

---

## 2C. Autonomous work packages

Each package is one agent session in its own worktree (`../pop_var_caller-cq-<package>`). Each
lists the files it touches; **packages sharing a file are ordered**, never run in parallel.
Effect is "none" for every package unless stated, and every commit passes the gate. A bug fix
inside a package is its own commit with its failing test, before any refactor in that package.

### Phase 0 — autonomous bug fixes, area by area

All seven touch disjoint files and can run in parallel.

| Pkg | Items | Files | Cost |
|---|---|---|---|
| **F1** parameters file | A3, A27; then prepare A1 for the owner (same code) | `calling/parameters_file/{to_run_parameters,validate,bindings}.rs` | S |
| **F2** census file | A15 (retry with a larger head read), A31 (sort extents) | `parameter_estimation/joint/census_file.rs` | S |
| **F3** fits | A25, A26, A32, A41 | `joint/fit.rs`, `joint/ssr_fit/penalised_fit.rs`, `penalised_fit/run_log.rs`, `ssr_fit/tract_exclusion.rs` | S |
| **F4** read input | A2, A4, A16, A20, A21, A22, A44 | `read/input/aligned_reads_reader/{bgzf_blocks,bam,container,prefetch}.rs`, `bam/index_preflight.rs`, `read/prepared_read.rs`, `read/input/sample_cursor.rs` | M |
| **F5** merge | A23, A24, A45 | `run/cohort_merge/observation_cache.rs`; test in `run/psp_caller.rs`'s tests | S |
| **F6** locus generation | A39, A40; A33, A34 as gated refactors | `pileup/{generator,chain_id_allocator,errors,genome_walk,active_read_set}.rs`, `locus_generation/ssr.rs` | S–M |
| **F7** small modules | A19, A28, A38, A42 | `repeat_catalog/parquet_file.rs`, `cli/provenance.rs`, `vcf/header.rs`, `cli/estimate_contamination.rs`, `neighbour_evidence/line.rs` | S |

### Phase 1 — dependency cycles

| Pkg | Changes | Files | Cost | Order |
|---|---|---|---|---|
| **C1** cohort observation (§1 row 1) | Move the six evidence types, `REFERENCE_ALLELE` and the three `MinAlt*` types to `src/cohort_observation.rs`; constructors stay in `build.rs`; `run::cohort_merge` re-exports for one commit, then call sites move. Name the constructors' file in the type's doc. | new `src/cohort_observation.rs`, `lib.rs`, `run/cohort_merge/{build,mod}.rs`, `calling/{evidence_shaping,inference/mod}.rs`, `calling/allele_candidates/{mod,generic,ssr,deciding_partial_reads,ssr_production_differential}.rs`, `calling/likelihood/{mod,generic}.rs`, `tests/ng_calling_loop_calls_genotypes.rs`, `tests/ng_candidate_selection_truth_recall.rs` | M | first; parallel with C3, C4 |
| **C2** three small cuts (§1 rows 2, 3, 22) | `run_context.rs` to a top-level module; the progress target and formatters to a crate-root module (`cli::run_output` re-exports; delete `peak.rs`'s copy; keep both duration styles); `fitted_inbreeding_of` into `parameter_estimation`, re-exported from `run`. One commit each. | `lib.rs`, `calling/likelihood/{run_context,mod,ssr,ssr_emission}.rs`, `calling/{mod,run_parameters}.rs`, `calling/parameters_file/{mod,from_run_parameters,to_run_parameters}.rs`, `parameter_estimation/{misread_after_runs,progress}.rs`, `joint/{census_file,stratum_fits,fit}.rs`, `cli/run_output.rs`, `run/{calling_progress,census_fit}.rs`, `neighbour_evidence/peak.rs` | S | after C1 (`lib.rs`, `likelihood/mod.rs`), C4 (`from_run_parameters.rs`), F1, F2, F3 |
| **C3** window coverage (§1 row 7) | `LocusSummary` to `locus_generation`; the `From<&Drawn>` impl into `observation_cache.rs`; fix `window_coverage/mod.rs:28-32`. | `window_coverage/{depth,mod}.rs`, `run/cohort_merge/observation_cache.rs`, `locus_generation/mod.rs` | S | after F5; parallel with C1 and C4 |
| **C4** read, psp, read phase (§1 rows 8, 9, 22) | `provenance_parameters` and its prefix constants to `run/read_filter_provenance.rs`; `read_phase/stage.rs` to `vcf/`; the three `parameters_file` tests call `validate()` directly (pfile 4). | `read/filtering.rs`, `run/{gatherer,psp_caller,mod}.rs`, new `run/read_filter_provenance.rs`, `psp/header.rs`, `read_phase/{stage,mod}.rs`, `vcf/{mod,writer}.rs`, `calling/parameters_file/from_run_parameters.rs` (tests) | S | after F5 (`psp_caller.rs`); parallel with C1 |

After phase 1, the dependency script in `tmp/cq_cross/` (cross, summary) should report no cycle
except `vcf` ↔ `neighbour_evidence`, which is kept on purpose.

### Phase 2 — dead code

All six can run in parallel once their phase-0 and phase-1 predecessors are merged.

| Pkg | Changes | Files | Cost | After |
|---|---|---|---|---|
| **D1** `src/bam` | Delete `segment_reader`, `segment_merge`, `bam_input`, `cram_input`, `record_specs` and what only they use (io A1, cross X1); `FilterBucket::BaqRejected` (cross X4); the two small `bam` items in io R5; replace the module-wide `allow` in `read/input/aligned_reads_reader/mod.rs:102-106` with item-level ones (cross X3). Check CRAM coverage first (§1 row 5). Retire `MappedRead` (io A2, first step) if time allows. | `src/bam/*`, `read/aligned_read.rs` (parity test), `read/input/aligned_reads_reader/mod.rs`, `doc/devel/ng/arch/alignment_file.md` | M | F4 |
| **D2** calling | Reduce `DiscoveryMode` to `Off`/`BeforeTheLoop` and remove the once-only loop in `call_locus` (inf 4.3); `#[cfg(test)]` the test-only functions and drop their re-exports (§1 row 10; cand R2, T1; inf 4.6 first row); item-level silencers in `allele_candidates/ssr.rs` (cross X3). | `calling/inference/{mod,summarise_condition}.rs`, `calling/mod.rs`, `calling/likelihood/{mod,generic}.rs`, `calling/allele_candidates/ssr.rs` | S | C1, C2 |
| **D3** parameter estimation | Delete `SharedReadGroup` (est R2); the twenty-bin ladder's unused parts, keeping `for_census` and its three constants (est R3 — changing them breaks every census on disk); `TractLikelihoods` and test-only scoring, porting the library-order test (fits 18); dead items and visibility (fits 11, 19; cross X4 `each_stutter_multiplier`). | `joint/{census,stratum_fits}.rs`, `depth_bins.rs`, `joint/ssr_fit.rs`, `ssr_fit/{read_group_stutter,penalised_fit,expected_counts}.rs`, `fit/standard_errors.rs` | M | F3, C2 |
| **D4** merge and run | Port the displacement tests to `Frontier` and delete `Organiser` (merge B2 / A46); `#[cfg(test)]` the test-oracle round driver in `rounds.rs` (merge A3), which also makes merge D2 moot; fold `judge_here` into `judge` (merge D4); drop the test-only fields and round-trip accessors on `AlignedFilesVariantCaller` (run A4). | `run/cohort_merge/{organise,rounds,stream,serial,close}.rs`, `run/callers.rs` (`:224-235`, `:430-521`, tests) | M | F5, C1 |
| **D5** locus generation | Delete `column_census` and `minted_error_census` and their call sites (locus A3; removes an atomic write at nearly every position — the gate should show it faster or equal); delete the two test-only methods (locus D5). | `pileup/{mod,genome_walk,fast_column,minted_error_census,generator,open_record}.rs`, `locus_generation/mod.rs` | S | F6, C3 |
| **D6** small modules | Delete `iter_ext` (cli X1); the third tiler and `bases_match` in `tandem_repeat.rs` (cli X2); `MarginalAligner` and its stand-in tests (io R3); `#[cfg(test)] mod` for `mode_equivalence` and `cross_platform_digests` (cli A6); the test-only `pub` items in cross X4's last group; update `lib.rs`'s module list (cross R1). | `lib.rs`, `iter_ext/`, `tandem_repeat.rs`, `region_typing/mod.rs` (doc at `:402`), `alignment/mod.rs`, `cli/{mod,run_output}.rs`, `neighbour_evidence/{refine,stage,templates}.rs`, `read_phase/fragments.rs` | S | C2, C4, F7 |

### Phase 3 — deduplication

Shared files force four chains: U1 → U2; R1 → R2 → R3; U3 → U7; and U9 after U1, R3, U7 and
U8, because it touches a file of each. Every other pair shares no file and may run in parallel.

| Pkg | Changes | Files | Cost | After |
|---|---|---|---|---|
| **U1** numerics | `float::ln_sum_exp` (§1 row 4); `log_sum_exp_2` moved verbatim to `float`; one Cholesky inverse and one Jacobi rotation in `joint/linear_algebra.rs`, each caller keeping its own thresholds and the STR side's summation order (fits 13, 14). | `float.rs`, `calling/quality/mod.rs`, `calling/genotype_prior/dirichlet_multinomial.rs`, `joint/{fit,contamination}.rs`, `fit/standard_errors.rs`, `ssr_fit/penalised_fit.rs` | S | D3 |
| **U2** fits | One arrow elimination for the standard errors and the Newton step (fits 15; keep loop order); call the read-probability functions in the error-rate M-step (fits 16); the five small copies (fits 17); one name for the concentration (fits 20); named clamp constants (fits 21, constants only). | `joint/fit.rs`, `fit/standard_errors.rs`, `ssr_fit/penalised_fit.rs` | M | U1 |
| **U3** calling | One `cap_and_admit` for both selection paths (cand D1) with `pooled_reads` (cand D4); `partial_read_fits` calls the generic prefix/suffix tests (cand D2); `respelled_over` beside `FlankDeletion` (cand A2); resolve the STR prior once per locus (inf 3.1). cand R3 only if D1's signatures make it free. | `calling/allele_candidates/{mod,generic,ssr,deciding_partial_reads}.rs`, `calling/likelihood/{ssr_emission,generic,ssr}.rs`, `calling/inference/summarise_condition.rs`, `locus_generation/mod.rs` | M | D2, D5 |
| **R1** run drivers, part 1 | One `call_and_assemble` for the per-locus record build shared by both drivers (run D1); one `RunReport` constructor (run D5); assert sample names against walkers (run R2). | `run/callers.rs`, `run/report.rs` | M | D4 |
| **R2** run drivers, part 2 | A37 as a bug commit, then one `RecordRelease` for the hand-over chain (run D2) and one `CallingSettings` struct (run D3). Direct mode keeps a batch of one. | `run/{callers,psp_caller}.rs`, `cli/{call_from_alignments,call_from_psps,calling_run}.rs` | M | R1, C4 |
| **R3** calling commands | Flatten ~25 shared flags into `CallingFlags` (cli D1) — **`--help` of both commands must be byte-identical before and after, else stop and ask**; `finish_the_run` for the shared tail (cli D2); tie the two `--max-depth` defaults with a compile-time assert (cli D3); one `FIRST_POSITION_OF_A_CONTIG` (cli D5). | `cli/{call_from_alignments,call_from_psps,calling_run,resume_calling,mode_equivalence,cross_platform_digests}.rs` and their tests, `run/{records,depth_ceiling}.rs`, `locus_generation/pileup/mod.rs` | M | R2, D5, D6 |
| **U5** merge and filter | One window-lookup helper for the cache and the sealed round (merge D1); one read-ahead loop for passes two and three (merge D3); the histogram codec next to `SampleHistogram` (merge A2); fold `paralog/calibrate.rs` into `calibration.rs` (merge D5); narrow `cohort_merge` visibility (merge A5). | `run/cohort_merge/{observation_cache,observation_cache/round_view,recorded_windows,mod}.rs`, `run/paralog_filter/{pass_two,pass_three,resume}.rs`, `window_coverage/` (codec destination), `paralog/{calibrate,calibration,mod}.rs` | M | D4, C3 |
| **U6** parameters file | One "format-N table in an older file" refusal (pfile 8); one `path()` per keyed row type (pfile 9); the two one-liners (pfile 10). | `calling/parameters_file/{validate,to_run_parameters,mod,to_toml}.rs` | S | F1, C2, and B1 if approved |
| **U7** locus generation | `CursorLedger` shared by the two generators (locus D1); merge the two flank-deletion predicates (locus D2, the S half; the CIGAR iterator is optional); one right-junction predicate (locus D3); one copy of the per-op event code (locus R5); a `Range<u32>` tract in the flank helpers (locus R2). | `locus_generation/{ssr,ssr_bundle,mod}.rs`, `pileup/{generator,open_record,cigar_cursor}.rs` | M | D5, U3 (`locus_generation/mod.rs`) |
| **U8** reads and alignment | Call `indel_norm::left_align_indels` directly and drop the four layers around it (io D1); retire `MappedRead` if D1 did not (io A2); optionally the shared delimiter traceback (io D2). | `read/{left_align,prepared_read,aligned_read}.rs`, `alignment/{indel_norm,left_align_structured,mod,ssr_best_path_flat_gap,ssr_unit_robust}.rs` | M | D1, D6 |
| **U9** shared utilities | One MD5 encoder (§1 row 15); `ContigBounds::of_list` for the four refusing callers (§1 row 16); `ContigList::id_of` with a name map (A36, cross D5); `GenomeRegion::overlaps`, one `logit`/`sigmoid` in `paralog`, `ContigHasher` shared (cross D8 "yes" rows); the Phred table shared by `prepared_read.rs` and `open_record.rs` (§1 row 18); the two `region_typing` re-implementations (cli D5); `PackedDepthCodes` as bytes and est R8. | `reference_info.rs`, `types.rs`, `regions.rs`, `calling/parameters_file/bindings.rs`, `psp/header.rs`, `repeat_catalog/{parquet_file,reader,segments,strata}.rs`, `cli/{calling_run,typed_regions,inspect_psp,run_ground,estimate_contamination,estimate_parameters}.rs`, `vcf/header.rs`, `run/{gatherer,explain,edge_deletion_gathering}.rs`, `psp/segmentation_section.rs`, `paralog/{prior,windowed_prior}.rs`, `region_typing/{reference_runs,mod}.rs`, `repeat_catalog/segments.rs`, `read/prepared_read.rs`, `pileup/open_record.rs`, `joint/{loci,census,contamination}.rs` | M | U1, R3, U7, U8 (last in the phase) |
| **U4** neighbour evidence (gated speed) | One "weigh the shares" loop (cli D5); scratch buffers for the per-sample and per-copy allocations (cli R3). Hot loop: the gate's timing must not move. | `neighbour_evidence/{model,records,refine}.rs` | S–M | D6 |

**Test helpers** (cross T1, inf 3.2, pfile tests): add `GenomeRegion::on` and `Ploidy::DIPLOID`
in U9, then convert helpers only in files a later package opens anyway. No sweep.

### Phase 4 — file splits

Each split is pure moves, one commit per file moved, nothing else in the commit (so `git blame -C`
follows it and a reviewer can see no edit hides in the move).

| Pkg | Changes | After | Parallel with |
|---|---|---|---|
| **S1** `calling/mod.rs` (5,848 lines) | Five files by type family, re-exports in `mod.rs`; rewrite the module doc's false "imports downward" claim (inf 4.1, 2.3). | D2, U3 | S2–S6 |
| **S2** `summarise_condition.rs` (10,664 lines) | Tests into `summarise_condition/tests/{fixtures,…}.rs` first; production into four files only if the test split leaves it hard to read (inf 4.2). | D2, U3 | S1, S3–S6 |
| **S3** `joint/fit.rs` (8,646 lines) | `fit/beta_quadrature.rs`, `fit/maximisation.rs`, `fit/expectation.rs` (fits 8); then `fit_jointly` and `maximise` into phases (fits 9). Cutting `one_position` into functions is gated on speed: drop it if the gate's timing moves. | U1, U2 | S1, S2, S4–S6 |
| **S4** `penalised_fit.rs` | Curve fitting to `penalised_fit/curve.rs`; matrices to the shared `linear_algebra.rs` (fits 22). | U2 | S1–S3, S5, S6 |
| **S5** `locus_generation/ssr.rs` | `ssr/{read_region,classify,tally,read_cap}.rs`; `ssr_bundle.rs` to `ssr/bundle.rs` and the six `pub(super)` plus eight `pub` items narrowed (locus A1, A2). `classify_read`'s extraction (locus R3) only if it passes fewer than about eight values. | U7 | S1–S4, S6 |
| **S6** census writer | `census/writer.rs` as a child module (est A4). | D3, U9 | S1–S5 |

Not split: `ssr_emission.rs` (cand R4, only if someone is about to work there), `types.rs` (cli A5),
the `census_*` files in `run` (run A3 step 2).

### Phase 5 — documentation drift

Scope per B21: false statements, dead links, plan step codes. One package per area, all parallel,
each after its area's last code package.

| Pkg | Items |
|---|---|
| **X1** calling | cand R1, T2 (A48); inf 1.1 (bundles are called), 1.2 ("final prune"), 4.5 (the three `*_parity.rs` docs and the defaults test name; replacing the 732 KB genotype-table fixture with a computed expectation is optional), 4.6 remaining rows. |
| **X2** parameters file | pfile 11, 12, 13; the step codes of pfile 14. |
| **X3** parameter estimation | est R4 table; the walked-bit comments of est B3 (`census.rs:22-25`, `:683-686`, `:849`); fits 19 doc references and the missing-example citations (`ssr_fit.rs:1428-1434`, `run/census_fit.rs:345`, `cli/estimate_parameters.rs:706`, `benches/ng_joint_fit_perf.rs:15`, `Cargo.toml:257`). |
| **X4** run and merge | run R1 table, including renaming `generic_path_generators` to `run_generators`; merge R1, R2, R3; `doc/devel/ng/arch/cohort_merge.md` §4; the `round_shape_for` doc once B14 is decided. |
| **X5** locus generation | locus R1 (false statements and dead `ia/specs` links only), R4; the SplitMix comment (§1 row 19). |
| **X6** io | io R1, R2, R5. |
| **X7** small modules | cli R1; the `types.rs` module doc (cli A5); optionally rename `window_coverage/production_parity.rs` (cli X4). |

---

## 3. Order of execution

```
Phase 0  F1  F2  F3  F4  F5  F6  F7                       all parallel
Phase 1  C1  C3  C4 in parallel, then C2
Phase 2  D1  D2  D3  D4  D5  D6                            all parallel
Phase 3  U1 → U2      R1 → R2 → R3      U3 → U7      U4  U5  U6  U8
         then U9 last
Phase 4  S1  S2  S3  S4  S5  S6                            all parallel
Phase 5  X1 … X7                                           all parallel
```

The order follows the program's rule: cycles first, because they fix the module boundaries the
later packages work inside; dead code next, because it shrinks what every later agent reads;
then deduplication, splits, and documentation. Within each area the auto bug fixes (phase 0)
come before any refactor of that area, so a refactor never has to preserve a known defect.

Owner items run as they are decided, each as its own branch. B1 goes after F1; B2 and B3 after
F3 and before U2/S3; B5 after F4 and F6; B16 (slippage re-fit) before S2, since deleting it
shrinks `summarise_condition.rs`; B17 and B18 before U1, since they shrink `contamination.rs`
and `fit.rs`. Decisions that change output should be merged in one batch where possible, so the
gate's baseline is re-recorded once.

---

## 4. Looked at, left alone

Collected from every report, one line each, so nobody re-opens them. "(report §)" points to the
reasoning.

**Calling**

- `GenotypePriorModel` trait with one implementation: test doubles use the seam (cand).
- The likelihood → genotype-prior import of `COUNT_PATH_DESYNC_THRESHOLD`: one direction, same subject (cand A4).
- `fill_batch_allele_copies`' every-batch check: one scan while every run has one batch (cand).
- `differ_by_one_indel` accepting one contiguous indel of any length: intended (cand).
- `deciding_partial_reads` counting nothing when a flank deletion is present: measured, intended (cand).
- `normalise_batch_rows` not renormalising after its clamp: off by a few × 10⁻¹² (cand).
- Two `Vec`s per pair in the STR misread fold: only for rare flank-deletion alleles (cand).
- `SelectionScratch` destructuring with named `_` fields: makes a new buffer a compile error (cand).
- `ssr_production_differential.rs` and the three `*_parity.rs` files: working regression checks; keeping production-parity tests is the owner's call, not a quality one (cand, inf 4.5).
- Asserts kept in release in selection, likelihood, `CallingScratch`, `score_one_sample`: each guards a silent wrong genotype (cand, inf).
- Leave-one-out subtraction written three times: too short to share (cand D3).
- `LocusGenotyper` trait with one implementation: gives `run` a narrow surface (inf 2.4).
- `RunParameters` and `FrozenParameters`: an owner and a `Copy` view, both used (inf 3.3).
- `non_zero_default` twice: two three-line `const fn`s (inf 3.2).
- Per-locus `GenotypeTable::build`: cached per thread (inf).
- The `1e-12` tolerance in `evidence_shaping.rs:93`, whole-read rounding on the STR path, `num_fwd: 0` for partials: all correct (inf).
- `RunnableCallingLoopConfig` with `Deref`: a validated newtype (inf).
- `repeat_tract_parameters.rs` under `inference/`: either home works (inf).
- `evidence_shaping.rs`'s long module doc: each paragraph records an unenforced join (inf).
- The `SsrEvidenceScratch::rebuild` allocations: measure before acting (inf 4.7).
- `calling/parameters_file/` and `run_parameters.rs` under `calling`: promoting them cuts no cycle (cross).

**Parameters file**

- The module's shape — separate file types, two projections, mirror enums, hand-written writer, `validate` apart from parsing (pfile 5).
- Replacing the row formatters with a custom serde serializer: saves ~300 lines, costs readability (pfile 7).
- The deliberate second copy of variant words, guarded by a test (pfile 6).
- `validate` running twice per run; demotion cloning the file; two identical-bodied functions naming different sources (pfile).
- `every_shape.toml` beside `every_shape_as_written.toml`; `EvidenceCount`'s hand-written `Deserialize`; a derive-based validator crate; `to_string()` against `a_toml_integer` (pfile).

**Fits and parameter estimation**

- Sharing golden-section, SQUAREM or Newton code between the two fits: different stopping rules, output would change (fits table, est D3).
- Thread-count independence of both fits; `newton_range`'s end-point claim; `best.expect`; `Parameters::continuing` not clamping; `count_times_ln`; the unreachable `!position_ln.is_finite()` branch; `step_beta_shapes`' guards; `RuleSlopes::of`'s finite difference; the 1e-300 floor; slippage groups threaded through while one is fitted; one-sample handling (fits).
- `fit_trace.rs`'s process-wide sink: opt-in diagnostic (fits, est).
- No wrong-way dependency in the fits (fits 10); `crate::run` in `census_file.rs` is doc links only (est A3).
- `census_file.rs` re-checking what `from_parts` asserts; class count checked twice; three spellings of "stratum"; `repeat_strata.rs` as shared vocabulary; `sequencing_batches.rs`' test-only, specified-but-unbuilt code; `first_disagreement` in a test; `calibration.rs`'s `i128` sum; `add_generic`'s `pending`; `mark_walked`'s cost; contamination's two passes; `checked_spectrum`'s nested module; `Sections::Backed` reopening the file (est).

**Run and merge**

- The merge continuing after the output refuses a record: wastes time, corrupts nothing (run).
- The round driver's `Drop` flush, `assert_eq!(displaced, 0)`, the psp source's too-large-record guard, `EdgeDeletionGathering`'s `HashMap` and contig-less key, the fragment keeper's thread split, the 11–13-argument waivers, psp mode's three paired lists, a failed-locus progress count (run).
- The walker and the gatherer sharing 15 lines of progress tracking (run D4).
- `stream.rs` shutdown and ordering; `cover` against `cover_in_parallel`; window truncation; `seal_round`; pass two/three shutdown; the windowed prior's serial sum; `ParalogScorePrecompute`'s sizing; infinite-window rejection; the 9-column check; `resume.rs`'s fields; `refuse_malformed_analysed_regions`; reference panics in `LocusReferenceBases::over`; the probe's panic on write error; `paralog`'s imports; `production_parity.rs` (test-only); the quadratic `AlleleGroupTally` scan (merge).
- `crate::cli` in `stream.rs`/`rounds.rs`: doc text only (merge A4).
- The round-size policy living in `cli/calling_run.rs` (cli A2).
- Moving the three `census_*` modules out of `run` (run A3 step 2): churn for findability only.

**Locus generation**

- `begin_region`'s one-region look-back; the satellite right junction; STR loci with zero reads; `q_sum` in read order; `HashMap` before a total sort; one reservoir seed for all samples at a locus; `u32` casts; empty left flank at a contig start; `in_smallest_form`; `slide_reach.rs` indexing; generic and STR read caps different by design; `decompose()` kept as an independent reference; the explicit `Drop`; error variants without `#[from]`; position overflow checks; `AHashMap` never iterated; the three aligner test helpers (locus).
- The two generators' `next_locus` assembly; `read_sampling.rs`'s finaliser; the STR SplitMix increment (locus D4).
- The STR CIGAR walks against `pileup/cigar_cursor.rs`: different questions, should evolve apart (locus D2).

**I/O**

- The psp field table checked against a hand-written codec, pinned by golden bytes; no table-driven codec (io A5).
- Three byte cursors' error mapping, ~30 lines each (io A5).
- psp version checks, `count as u32`, truncation coverage; CRAM MD5; BAM chunk merging; contig checks at open; read-group rules; unconditional secondary/supplementary drops; duplicate detection; MAPQ 255; prefetch determinism; BAM and CRAM reader structs apart; `read_phase` checks; `alignment` numerics; `left_align.rs`; `compute_adaptor_boundary` casts (io).

**CLI and small modules**

- `vcf` ↔ `neighbour_evidence`: the stage sits in front of the writer (cli A4, cross A6).
- `types.rs` as one file; revisit when one group passes ~600 lines of code (cli A5).
- The test-only migration modules (`golden_catalog`, `production_post_filter`, `scanner_parity`, `window_coverage/production_parity.rs`): live tests (cli X4).
- `psps_named_by` twice; `contig_bounds`'s and `tmp_path_for`'s small copies (now partly covered by U9); the repeat-criteria flags whose help text differs on purpose; ~12 seeded test generators (cli D4, D5).
- Thread-pool setup in six commands (cli D5, cross D7).
- Five temp-then-rename writers with different shapes (cross D6).
- VCF `PL` with NaN or all `−∞`, `AF` at zero copies, multi-allelic order, `GT`/`REPCN` order, padding at a contig start, the writer's order check and leftover `.tmp`, header escaping of names, BED parsing, `u32` coordinates, `memory_limit.rs`, `reference_info.rs` and `open_files.rs` visibility, `splitmix.rs` against `random.rs`, the twin region-count types, examples using `cli` as a library, `round_shape_for` reading the thread count (see B14), neighbour-evidence numerics, window accumulator and catalog scan order (cli).
- The three `test_fixtures.rs` files: layered, not copies; two k-way merges of different things; BED parsing once; the item-level `expect(dead_code)` idiom; `genetics::lgamma` against a test-only Lanczos; `unhex` accepting a leading `+`; `psp`'s imports (cross).
- Numeric near-copies that must not merge: `log_add_exp` against `log_sum_exp_2`; the `ssr.rs` Phred conversion; the two binomial tails; `logistic` against `sigmoid` (cross D1, D8).
