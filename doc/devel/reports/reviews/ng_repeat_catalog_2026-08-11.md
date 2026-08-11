# Code Review: ng_repeat_catalog
**Date:** 2026-08-11
**Reviewer:** rust-code-review skill (orchestrator + nine per-category sub-agents, each in its own git worktree)
**Scope:** the reference tandem-repeat catalog — the module, its CLI subcommand, the observer seam in the reference pass, and the differential test
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** a module and its immediate surface — not a diff. The catalog's eight source files as they stand, plus the newest commit's tally.
- **Reviewed against:** commit `8a299f2afdefd6f56399b91265d9e40cb6cffbd9`, branch `ng-cohort-stutter-readgroups`.
- **In-scope files:**
  - [src/ng/repeat_catalog/mod.rs](../../../../src/ng/repeat_catalog/mod.rs)
  - [src/ng/repeat_catalog/builder.rs](../../../../src/ng/repeat_catalog/builder.rs)
  - [src/ng/repeat_catalog/criteria.rs](../../../../src/ng/repeat_catalog/criteria.rs)
  - [src/ng/repeat_catalog/parquet_file.rs](../../../../src/ng/repeat_catalog/parquet_file.rs)
  - [src/ng/repeat_catalog/reader.rs](../../../../src/ng/repeat_catalog/reader.rs)
  - [src/ng/repeat_catalog/row.rs](../../../../src/ng/repeat_catalog/row.rs)
  - [src/ng/repeat_catalog/segments.rs](../../../../src/ng/repeat_catalog/segments.rs)
  - [src/ng/repeat_catalog/strata.rs](../../../../src/ng/repeat_catalog/strata.rs)
  - [src/ng/repeat_catalog/tally.rs](../../../../src/ng/repeat_catalog/tally.rs)
  - [src/pop_var_caller_exp/repeat_catalog.rs](../../../../src/pop_var_caller_exp/repeat_catalog.rs)
  - [src/ng/reference_info.rs](../../../../src/ng/reference_info.rs) — the observer seam only
  - [tests/ng_repeat_catalog_differential.rs](../../../../tests/ng_repeat_catalog_differential.rs)
- **Deliberately out of scope:** `src/ssr/` (frozen production); the internals of `src/ng/region_typing/` (read to judge coupling, not reviewed); `examples/`.
- **Categories dispatched:** reliability (always), errors (always), naming (always), idiomatic (always), refactor_safety (always), smells (always), module_structure (multi-file scope), defaults (the catalog's floors are deliberately not the caller's), unsafe_concurrency (a rayon fan-out in the builder), extras (a parser, a byte-identical output producer, a hot path, and a public API). `tooling` was skipped: the only `Cargo.toml` change in this work is one dependency line, already reviewed when it landed.

## 2. Verdict

**Request-changes.** The classification logic is sound and its differential is the reason to trust the file; what is not sound is the *guarding* around it. Three refusals the design calls fatal do not happen at all, the newest tally is wrong whenever a caller asks for anything smaller than a whole contig, and a satellite crossing a requested edge comes back the wrong size.

## 3. Execution status

Run in the dev container on the main checkout at the reviewed commit:

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | exit 0, clean |
| `cargo clippy --all-targets --all-features` | warnings only, all in `examples/shared/stutter_model.rs`, `examples/shared/stutter_table.rs`, `examples/ng_str_stutter_rate.rs` — pre-existing, out of scope |
| `cargo test --tests` | lib `3304 passed; 0 failed; 9 ignored`; every integration test green; the differential `9 passed` |
| `cargo doc --no-deps` | 12 unresolved intra-doc links, **none in `src/ng/repeat_catalog/`** — pre-existing |
| `cargo audit` | not run |

Sub-agents ran a further **35 mutations** across their worktrees. Findings labeled "Needs verification": **0** — every finding below was reproduced by running code, and each sub-agent reported mutations-run / survived / changed-no-behaviour as required.

## 4. Open questions and assumptions

1. **What does a tally mean when the caller asked for part of a contig?** Today the walk charges the whole contig (its scan span is always the whole contig, whatever the BED said) and the catalog charges only the rows it read. They cannot both be right, and the answer decides B2. Neither is obviously the intended one: the walk's own `spans` counter deliberately counts *requested* regions rather than walked ones, which argues the coverage counter should be scoped the same way. **This one is the owner's** — it changes a number a `--regions` run prints, on both code paths.
2. **Is `RepeatCatalogWriter`'s `.tmp` path meant to be safe against two concurrent builds of one reference?** It is deterministic, so two builds interleave into one file and both rename it into place (Mi-group). Affects nothing today because the command is run by hand.
3. **Is the region read's page-index pruning a requirement or an aspiration?** M17 shows it does not happen; the doc says it does.

## 5. Top 3 priorities

1. **B1 — three fatal refusals never fire.** The tool version is never compared, the scoring weights are never compared, and both are documented as "equal or refuse". A detector change or a re-weighted scan leaves every consumer silently reading a catalog that does not describe what it thinks.
2. **B2 — the tally is short by an order of magnitude for any region read**, and the only test of it asks for whole contigs, the one regime where the two definitions coincide.
3. **B3 — a satellite crossing a requested span's edge comes back clipped**, where the walk emits it whole: a 1,234-base array reported as 100 bases.

## 6. Findings

### Blocker

#### B1: src/ng/repeat_catalog/reader.rs:52 — Two of the three "equal or refuse" checks never run
**Categories:** errors, extras, smells, defaults (convergent — four sub-agents)
**Confidence:** High (proven by probe: a catalog stamped `ancient-0.0.1` opened and served its rows under running version `0.1.0`)

`RepeatCatalogError::ToolVersionDiffers` ([mod.rs:373](../../../../src/ng/repeat_catalog/mod.rs#L373)) is declared and **constructed nowhere in the repository**; `open_checking_against_reference` never compares `header.tool_version` against the running build. `check_scored_with` ([reader.rs:314](../../../../src/ng/repeat_catalog/reader.rs#L314)) is the only site that builds `ScoringWeightsDiffer`, and **nothing calls it** — so `genome_segments`, `str_loci`, `count_loci_per_stratum` and `sample_loci_per_stratum` never refuse a catalog scored with different weights.

Three doc comments assert the opposite, including [parquet_file.rs:47](../../../../src/ng/repeat_catalog/parquet_file.rs#L47) ("the header's `tool_version`, which is checked on read") and [mod.rs:275](../../../../src/ng/repeat_catalog/mod.rs#L275) ("**Equal or refuse**").

**Why it matters:** the spec's whole safety argument is that a wrong file is refused rather than served short. Different scoring weights are a *different set of tracts*, not a subset — no filter reproduces them — so this is the one axis where a silent answer is not merely incomplete but wrong.

**Fix:** compare `tool_version` against `env!("CARGO_PKG_VERSION")` inside `open_checking_against_reference`, and fold `check_scored_with` into `check_serves` so every criteria-taking method runs it. Both need a test; the differential's fixtures currently stamp `"differential"`, so they will need the running version.

#### B2: src/ng/repeat_catalog/tally.rs:145 — The tally is only correct for a whole-contig scope, and the only test of it asks for nothing else
**Categories:** reliability
**Confidence:** High (measured on unmutated code)

`ContigTally::repeat_bp` documents itself as "the whole contig's, not the requested stretches'", and `segments_of_contig_in` ([segments.rs:64](../../../../src/ng/repeat_catalog/segments.rs#L64)) repeats the claim. It is false: `GenomeSegments::next` feeds `segments_of_contig_in` only the rows `rows_for` returned, which for a region scope are the rows inside the *widened windows*, so the coverage and every rejection counter are computed over a fraction of the contig.

Measured on half of the differential's own `chr_random`, with the regions themselves still comparing equal:

| | repeat bases with no locus | copy floor | no clean trim |
|---|---|---|---|
| the walk | 1,822 | 38 | 68 |
| the file | 187 | 0 | 12 |

Asking for the same contig **whole**, both sides report 381 / 38 / 68 — which is why the suite is green.

**Why it matters:** the tally exists so a consumer that stops walking the reference keeps its report. Every sharded or `--regions` run is a region read.

**Fix:** answer open question 1 first. Then either charge the whole contig (read the contig's rows for the tally, or store per-contig totals), or scope the counters to the requested regions on both sides — and add the partial-span case to the differential either way.

#### B3: src/ng/repeat_catalog/segments.rs:85 — A satellite crossing a requested span's edge is clipped by the file and whole from the walk
**Categories:** reliability
**Confidence:** High (measured on unmutated code)

`segments_of_contig_in` emits `SsrSegment` and `SsrBundle` whole and clips everything else. The walk's rule ([region_typing/mod.rs](../../../../src/ng/region_typing/mod.rs), `clips_at_a_bed_edge`) is that **only `Generic` clips** — a satellite's extent *is* its claim, and that function's own doc records this as the case an earlier review got wrong. Asking for a 100 bp window inside a 1.2 kb array:

```
from the file: [ Satellite  600..700  ]
from the walk: [ Satellite  201..1434 ]
```

`a_region_subset_from_the_file_equals_the_walk_over_the_same_spans` asks for two spans of `chr_clean`, a contig with no satellite, so the difference cannot appear.

**Fix:** add `RegionKind::Satellite` to the `whole` arm, and give the region differential a contig that carries an array.

### Major

#### M1: src/ng/repeat_catalog/parquet_file.rs:425 — "This contig has no rows" and "this file is unusable" are the same answer
**Categories:** errors
**Confidence:** High (mutation, behaviour change measured)

`rows_overlapping` returns `Ok(vec![])` when no row group's statistics identify a contig. Disabling statistics in the writer made every read come back empty and `Ok`; `PageIndexPolicy::Required` did not reject the file. **Fix:** return `Unreadable` when a row group has no usable statistics, and keep the empty vector for a contig genuinely absent.

#### M2: src/ng/repeat_catalog/parquet_file.rs:431 — A parquet file with the catalog's header key but a different schema panics
**Categories:** errors
**Confidence:** High (reproduced)

`index out of bounds: the len is 1 but the index is 1`, inside parquet's `ProjectionMask::roots`. `row_from_batch` handles the same downcasts through `Option`; the predicate does not. **Fix:** validate the schema against `catalog_schema()` on open and return `Unreadable`.

#### M3: src/ng/repeat_catalog/mod.rs:305 — `reach_into_window` adds two `u64`s unchecked
**Categories:** errors, extras
**Confidence:** High (reproduced)

`longest + bundle_threshold`. `bundle_threshold` is documented as an axis a reader may set freely; a region read at `u64::MAX` panics in debug and, in release, wraps to a *narrow* window — classification from an incomplete row set, with no error. A header claiming `longest_tract_bp = u64::MAX` yields a reach of 14. Its own callers already use `saturating_*`. **Fix:** `saturating_add`.

#### M4: src/pop_var_caller_exp/repeat_catalog.rs:88 — `--min-flank-bp` below 15 is accepted, and reading the catalog it writes panics
**Categories:** defaults, reliability
**Confidence:** High (probe)

`tally.rs`'s `CatalogRejectionCounts::add` refuses to charge `FlankClamped` and `debug_assert!(false)`s, on the argument that the file's 15 bp flank floor means no stored tract can clamp. Nothing enforces that floor: `--min-flank-bp 0` builds a catalog whose rows *can* clamp, and reading it panics (release: the bases vanish from the breakdown while staying in the total). Separately, `min_flank_bp` passes `serves` and is then **never applied at read time** — a reader asking for 30 bp of flank gets the same loci a reader asking for 15 gets. **Fix:** either enforce a floor on the flag, or make the tally handle the case honestly; and apply the reader's flank at read time or document that it is a build-time axis only.

#### M5: src/pop_var_caller_exp/repeat_catalog.rs:75 — `--max-period 9` is accepted and the header then claims a range the file cannot hold
**Categories:** defaults
**Confidence:** High (probe: header says `1..=9`, stored periods are 1 and 6, `serves(1..=9)` returns `Ok`)

`PeriodRange::new` only rejects `min == 0` and `min > max`; `Motif::new` caps at `MAX_MOTIF_LEN` = 6. A reader asking for period 8 gets an empty answer instead of `CriteriaRefusal::PeriodRange`. **Fix:** reject a ceiling above `MAX_MOTIF_LEN` at the flag.

#### M6: src/ng/repeat_catalog/mod.rs:225 — `trimmed` and `purity` are co-dependent `Option`s, and the decoder can produce the combination that cannot exist
**Categories:** idiomatic
**Confidence:** High

Their own docs say "`None` exactly when `trimmed` is", but `parquet_file.rs` decodes the nullable columns independently, so a file with a cut span and a null purity decodes to `Some`/`None` and `finish_from_row` charges it to `NoCleanTrim` — a locus silently lost, under the wrong reason. `segments.rs` already pays for the shape with a defensive `?` whose comment argues it cannot fire. **Fix:** one `Option<TrimmedTract { span, purity }>`. The sub-agent applied it in its worktree: lib and differential green, on-disk format unchanged.

#### M7: src/ng/repeat_catalog/builder.rs:186 — `--threads N` does not set the number of threads
**Categories:** unsafe_concurrency
**Confidence:** High

It bounds the batch; the scan runs on rayon's process-global pool, which neither `main_exp.rs` nor `run_repeat_catalog` ever sizes — unlike all five other parallel entry points in the crate. `--threads 64` on an 8-core box pays 64 contigs of memory (~16 bytes per base of the largest contig each) for 8 scans, and `RAYON_NUM_THREADS` silently overrides the flag. **Fix:** build a scoped `ThreadPool` of `threads`, or rename the flag to what it bounds.

#### M8: src/ng/repeat_catalog/reader.rs:165 — The window widening is unpinned, and losing it turns a satellite into a locus
**Categories:** reliability, extras
**Confidence:** High (two independent mutations, both survived, both proven to change behaviour)

Deleting the `reach` widening, and separately dropping the `longest_tract_bp` term from `reach_into_window`, each leave the whole suite green. With the first, a window on `chr_bundled` returns `SsrSegment @243..274` where the correct answer is `SsrBundle @200..274`; with the second, a window beside a 1.2 kb array returns `SsrSegment @1405..1434` where the walk returns `Satellite @201..1434`. **`longest_tract_bp` is a whole per-contig header column whose doc calls it "what makes a windowed read exact", and nothing tests it.**

#### M9: src/ng/repeat_catalog/criteria.rs:90 and segments.rs:193 — Two `..Default::default()` tails silently inherit the caller's policy
**Categories:** refactor_safety
**Confidence:** High (mutation: a new field compiles clean at both sites and flags four others)

`StrRepeatCriteria::default` closes with `..SsrSegmentCriteria::default()`, and `repeat_features_of_contig` builds its `TypedRegionConfig` with `..TypedRegionConfig::default()`. The one value that is supposed to be the catalog's deliberately-widened policy would silently take the short-read calling value for any field added later, and `serves` compares only three axes so it would not refuse. The walk's own equivalent carries a "Do not replace with `..`" comment. **The differential cannot catch this** because it builds its walk side the same lax way at four places, so both sides would default together.

#### M10: src/ng/repeat_catalog/parquet_file.rs:311 — `decode_header`'s `>=` guards silently drop extra fields
**Categories:** refactor_safety
**Confidence:** High (mutation: appending a tenth `criteria` field changed the file's bytes and left every test green, including the round-trip test that exists to catch this)

**Fix:** exact-length checks per line, with the count in the error.

#### M11: src/ng/repeat_catalog/segments.rs:237 — `admit`'s score gate is unpinned
**Categories:** extras, reliability
**Confidence:** High (mutation survived, behaviour change proven)

Deleting `row.score < class.min_score` leaves all nine differential tests and every unit test green while turning a below-floor row into a locus.

#### M12: src/ng/repeat_catalog/parquet_file.rs:350 — A non-hex contig digest decodes to "no digest recorded", and that contig's check is then skipped
**Categories:** errors
**Confidence:** High (proven)

A corrupted digest field becomes `md5: None`, and `check_against_reference` skips the comparison for exactly that contig — the `.fai`-only allowance swallowing a malformed file. **Fix:** a malformed digest is `Unreadable`; only an absent one is `None`.

#### M13: src/ng/repeat_catalog/mod.rs:55 — Regions naming a contig the header does not hold are silently ignored
**Categories:** errors
**Confidence:** High (proven: `count_loci_per_stratum` over `ContigId(7)` of a 2-contig catalog returns `Ok` with total 0)

#### M14: src/ng/repeat_catalog/reader.rs:423 — `RegionSegments` is an unreachable second implementation of the region read, and it has no tally
**Categories:** module_structure, smells, idiomatic, errors, naming, defaults (convergent — six sub-agents)
**Confidence:** High (deleted and rebuilt; only unused imports follow, both suites green)

A `pub struct` with a full `Iterator` impl that **nothing ever constructs**, duplicating `GenomeSegments`' region path without the counters. Its private feeder `rows_of_contig_in` is dead with it. Four more items have no caller: `segments_of_rows` (whose `Result` has no `Err` path), `locus_span`, `row::clears_detected_copy_floor`, and `RepeatCatalog::check_scored_with` (B1). Two of them carry doc comments asserting callers that do not exist.

#### M15: src/ng/repeat_catalog/parquet_file.rs:487 — Arrow leaves the module it is promised to stay inside
**Categories:** module_structure, idiomatic, smells
**Confidence:** High

`pub fn row_from_batch(&RecordBatch, usize)` puts `arrow_array::RecordBatch` on the crate's public API, against that file's own module doc. `catalog_schema`, `decode_header` and `rows_overlapping` all compile at reduced visibility too.

#### M16: src/ng/repeat_catalog/mod.rs:18 — Every submodule is `pub mod`, so `dead_code` is blind to the whole folder
**Categories:** module_structure
**Confidence:** High (verified: changing seven of the eight to `pub(crate) mod` needs **zero** import changes anywhere; `cargo check --all-targets` clean)

74 `pub` items are crate-public API, which is why none of M14's five dead items ever warned. Only `criteria` is reached by module path from outside the folder.

#### M17: src/ng/repeat_catalog/reader.rs:150 — The region read decodes the whole contig's coordinate columns; the page index is not doing the work its doc claims
**Categories:** extras
**Confidence:** High (measured on the real tomato catalog)

parquet-rs's `RowFilter` does no page-index pruning of its own predicate columns; `with_row_selection`, the API that does, is never called. A 10 kb, 100 kb and 1 Mb window on tomato chromosome 1 each cost 18 ms, and the floor tracks the contig's length (3.4 ms on a 9.6 Mb contig) rather than the window's. Footer parsing alone is 0.4 ms. **The optimisation is real relative to reading rows, and absent relative to what the doc describes.**

#### M18: src/ng/repeat_catalog/builder.rs:115 — `finish`'s doc promises a check that does not exist
**Categories:** reliability, idiomatic
**Confidence:** High

"Fails … if the pass saw contigs this builder did not." There is no such code; a builder handed a short contig table writes a file whose later row groups no reader can reach.

### Minor

Twenty-nine Minor findings are recorded in full in the per-category files under `tmp/review_2026-08-11_ng-repeat-catalog/`. The ones worth naming here:

- **`ContigTally` is named for where the numbers came from, not what it holds**, and it is a third naming pattern beside `CatalogRegionCounts`/`CatalogRejectionCounts` and `BuildTally`. Worse, the same value is `rejected` on one and `rejected_by_reason` on the other, and both appear in one function. (naming)
- **`Str…` and `Ssr…` are two spellings of one concept inside one module** — `StrLoci` is the iterator whose elements increment `ssr_loci`. (naming)
- **"The seven columns of a catalog row" sits above a schema with nine.** The spec says seven too; it predates storing both spans. (naming, smells)
- **Column positions are bare literals at four sites** that must all track `catalog_schema()`, in the module documented as the one silent-failure site. (naming)
- **Three copy-paste pairs** — two `spans_on`, two `hex`, two `pair_with_rows` — and the copy floor applied twice, once inline and once in a named function whose doc describes a caller that does not exist. (smells, idiomatic)
- **`pair_with_rows`' subsequence precondition is documented and not enforced**; an unmatched survivor is dropped with no assertion, so a future reordering upstream loses loci rather than failing. The precondition does hold today: `prefilter` is a filter-and-collect and `split_bundles` pushes with strictly ascending index. (idiomatic)
- **`repeat_bp_with_no_locus -= bp` is an unchecked subtraction on a public counter** — a debug panic, a wrap in release. (errors, idiomatic, refactor_safety)
- **`CATALOG_MAX_STR_LEN_BP` is the one catalog default no test pins**: 500 → 1000 leaves all tests green while reclassifying a 700 bp array. (defaults)
- **The startup log records three settings and none of the ones that decide the file's contents.** (defaults)
- **The contig-order mismatch message prints the catalog's name twice** where it means the index: "contig chrA of the catalog is `chrA`". (errors, smells)
- **`row_from_batch` re-downcasts nine columns per row** — about 212 million downcasts on GRCh38's 23.6 M rows, hoistable to once per batch. (idiomatic)
- **`RowsByPeriod::for_period(0)` returns the beyond-the-table count** rather than 0. (naming, defaults, errors)

### Nits

Grouped, not enumerated: `IgnoreBases` names a type with a command; `f` carries the header's fields across 65 lines; `class`, `scan`, `pending`, `malformed` and `Admitted` are adjectives or participles where nouns belong; the requested regions are called `wanted`, `spans`, `windows` and `regions` in turn, and `wanted` also means criteria in `serves`; `cap` is undefined in prose and positional beside `seed`; internal labels "(B2)" and "open question 1" appear in doc comments; `impl<'a> ReadScope<'a>` where no method uses `'a`; a deterministic `.tmp` path shared by concurrent builds; `with_threads` coerces 0 to 1 silently; `_ => None` on two internal enums.

## 7. Out of scope observations

- **`segment_criteria.rs`'s gate archaeology is stale.** Its comment says four of classification's five rejection reasons are structurally zero and only the contig-end one is reachable. That was true before the copy floor moved to the cut tract (2026-07-20); measured now on tomato chromosome 1, all four fire — 207,854 bases no-clean-trim, 33,253 purity, 15,743 copy floor, 1,862 compound — and the contig-end one is the zero. Prose fix, in `region_typing`.
- **`prefilter` is fully `pub` while the five helpers beside it are `pub(crate)`.**
- **`TypedRegion` and `SsrSegment` are interchange types living inside a pipeline stage** — used in 7 and 9 files outside `region_typing`. A scheduling item, not a branch edit.
- **`RepeatCatalogWriter` has no `Drop`**, so an abandoned build leaves a `.tmp` behind.

## 8. Missing tests to add now

Grouped by what they hold. The `reliability` sub-agent supplied twelve bodies, ten of which fail on today's code.

**`RepeatCatalog::open_checking_against_reference`**
- `a_catalog_from_another_tool_version_is_refused` — header stamped with a different version; must error naming both versions. Catches B1.
- `a_reader_with_different_scoring_weights_is_refused` — same floors, different match reward; must error. Catches B1.
- `a_parquet_file_with_the_header_key_and_a_foreign_schema_is_unreadable` — catches M2 (currently a panic).
- `a_malformed_contig_digest_is_unreadable_not_absent` — catches M12.

**`GenomeSegments`**
- `the_tally_from_the_file_matches_the_walks_over_a_partial_span` — catches B2.
- `a_satellite_crossing_a_span_edge_comes_out_whole` — catches B3.
- `a_bundle_partner_outside_the_window_still_bundles` and `a_satellite_outside_the_window_still_absorbs_inside_it` — catch M8, and the second is the only test that would hold `longest_tract_bp`.
- `a_region_on_a_contig_the_catalog_does_not_hold_is_refused` — catches M13.

**`admit` / `segments_of_contig`**
- `a_row_below_the_score_floor_is_not_a_locus` — catches M11.

**`StrRepeatCriteria` and the command**
- `the_satellite_cap_default_is_five_hundred_bases` — catches the `CATALOG_MAX_STR_LEN_BP` Minor.
- `a_period_ceiling_above_the_motif_limit_is_refused` — catches M5.
- `a_flank_floor_below_the_catalogs_own_is_refused` — catches M4.

## 9. What's good

- **The differential is the right anchor and it is honest about its one exception** — the flank-floor difference is excluded explicitly and asserted to be edge-confined rather than assumed ([tests/ng_repeat_catalog_differential.rs](../../../../tests/ng_repeat_catalog_differential.rs)).
- **The classification is called, not re-derived.** Each of the six `pub(crate)` helpers widened out of step 3 is used exactly once by the catalog and by nothing else; `segments.rs` restates only the *order* of gates that need bases, which cannot be called. Verified by the module-structure pass.
- **The observer seam kept `reference_info.rs` a leaf** — it imports only `crate::fasta`, std and md5; the trait is infallible and the builder holds its own first error.
- **Byte-identity holds at genome scale, not just on fixtures**: 1-thread and 4-thread tomato builds are both md5 `b91d8037…`, matching the prebuilt file.
- **The coordinate conversion, the row order and `pair_with_rows`' precondition all killed their mutations loudly** — 8, 1 and 6 test failures respectively.

## 10. Commands to re-verify

Run from the project root, in the dev container:

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --tests --all-features -- -D warnings
./scripts/dev.sh cargo test --lib repeat_catalog
./scripts/dev.sh cargo test --test ng_repeat_catalog_differential
```

New once the fixes land: the twelve tests of section 8, and a rebuild of the tomato catalog compared by md5 against `b91d8037…` to confirm the on-disk format did not move.

---

### Author response convention

Address each finding by its identifier (B1, M5, …) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the open questions in section 4 first — question 1 gates B2.

*Per-category files, with every finding in full and the mutation accounting for each: `tmp/review_2026-08-11_ng-repeat-catalog/`.*
