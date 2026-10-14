# Code Review: segmenter_perf
**Date:** 2026-10-09
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the segmenter performance change on branch `segmenter-perf-combined` — `git diff 3a31a724 d4027a8e` (commits `74b2c508`, `d4027a8e`; `39336d54` adds a timing example)
**Status:** Approve-with-changes — the changes are made on this branch, in the commit that adds this report

---

## 1. Scope

**What was reviewed.** A performance change to typing the genome from the repeat catalog (the
parquet file of tandem-repeat tracts, one block of rows per contig). It must change no output. It
removes five costs that grew with the product of two counts:

1. the catalog's footer and page index, which describe every contig, were parsed again for every
   contig read; now they are parsed once when the catalog is opened (`CatalogIndex`);
2. the row filter is skipped when one requested window holds every row of a contig's block;
3. requested BED spans were compared with every item on their contig in eight places; now they
   are grouped by contig once and searched by binary search (`SpanIndex`, `OwnedScope`);
4. absorbing STR loci and bundle members into nearby satellites was a loop repeated until
   nothing changed, each pass comparing every feature with every satellite; now it is one sweep
   over the spans sorted by start (`absorb_by_sweep`), with the old loop kept as the oracle of a
   property test;
5. Arrow column types are checked once per batch instead of once per row.

**Reviewed against:** `d4027a8e` (base `3a31a724`).

**In-scope files:** `src/repeat_catalog/mod.rs`, `src/repeat_catalog/parquet_file.rs`,
`src/repeat_catalog/reader.rs`, `src/repeat_catalog/segments.rs`, `src/region_typing/mod.rs`,
`examples/segmenter_timing.rs` (lightly).

**Deliberately out of scope:** the speed-ups themselves (measured by the performance review,
`doc/devel/reports/reviews/perf_segmenter_2026-10-09.md` on branch `segmenter-perf`); code the diff
does not touch.

**Categories dispatched.** Two isolated sub-agents rather than one per category, because this
machine runs several 16 GB build containers for other sessions and one build container per
category would not fit:

- *reliability* and *refactor_safety* — one agent with its own worktree, building and mutating.
  It was cut off by a session restart before writing its files; what it had run survives as
  probe tests and three proptest counterexamples (kept in
  `tmp/review_2026-10-09_segmenter-perf/reliability_agent_probes/`), and the orchestrator ran
  the reliability work itself (§3, §8): eight mutants, the full suite, the identity oracle.
- *errors, naming, idiomatic, smells, module_structure, defaults, extras, unsafe_concurrency,
  float_portability* — one agent reading only, no build. Its nine files are in
  `tmp/review_2026-10-09_segmenter-perf/`.
- *tooling* — not dispatched: no `Cargo.toml` change.

## 2. Verdict

**Approve-with-changes.** The equivalence the change claims holds on every input the pipeline
produces today: the identity oracle reproduces its recorded output exactly, and the reasoning
behind each binary search and the sweep checks out. Three defects on inputs the pipeline does not
produce today, one error-message regression and the formatting/lint failures are fixed on this
branch (M1–M5 below), each with a test that fails without the fix. What remains is Minor. Merge
after this commit.

## 3. Execution status

All builds in the dev container (`scripts/dev.sh`), release profile.

| Command | Result |
|---|---|
| `cargo fmt --check` on the diff's files, at `d4027a8e` | exit 1 — `parquet_file.rs:530` and `examples/segmenter_timing.rs:53,144` unformatted (M5) |
| `cargo clippy --release --all-targets`, at `d4027a8e` | exit 0, **3 warnings in the library, all new**: `redundant_guards` at `reader.rs:119`, `needless_borrow` at `reader.rs:443` and `:481` (M5). Otherwise only the known warnings in examples `dhat_ng_joint_fit`, `dhat_ng_merge`, `dhat_ng_psp` and the vendored `noodles-cram` |
| `cargo test --release --lib`, at `d4027a8e` | 4,989 passed, 3 failed, 11 ignored. The 3 are the known release-mode should-panic tests (`left_align_repeated`/`left_align_structured` `an_offset_past_the_reference_panics`, `genetics::tests::lgamma_panics_on_non_positive_in_debug`) |
| `cargo test --release --test '*'`, at `d4027a8e` | 41 passed, 0 failed |
| `cargo test --release --lib`, with this review's fixes | **4,995 passed**, the same 3 known failures, 11 ignored (+6: the six tests this review adds) |
| `cargo test --release --test '*'`, with the fixes | 41 passed, 0 failed |
| `cargo clippy --release --all-targets`, with the fixes | exit 0; no library warning; only the known example/vendored warnings |
| `cargo fmt --check` on the diff's files, with the fixes | clean after `rustfmt` (see §10); `src/region_typing/segment_criteria.rs:1617,1627` remain unformatted, pre-existing and outside the diff |
| `scripts/promote_ng_oracle.sh` (tomato, 20 regions, 4 CRAMs, with the fit), with the fixes | **matches the baseline exactly once the reference path is the same** — see below |
| Before/after of each fix | each new test fails without its fix and passes with it (§6) |
| Eight mutants against the `repeat_catalog` and `region_typing` tests | §8 |

**The identity oracle.** The script reported three digests differing from
`scripts/promote_ng_oracle.baseline`: both default-parameter VCFs (`065d4084…` against `1109b1a6…`)
and the VCF called with the fit (`88543d54…` against `155fba62…`). The other five lines —
both parameters files, the window rows, the coverage histograms and the fitted parameters — match.
The three VCFs differ only in their `##reference=` header line, which records the path the
reference was given by: this run passed
`/Users/jose/devel/pop_var_caller-segmenter-perf-combined/tmp/oracle_inputs/S_lycopersicum_chromosomes.4.00.fa`,
the run that recorded the baseline passed `tmp/compound_fix/tomato_ref/S_lycopersicum_chromosomes.4.00.fa`
(its outputs are in `pop_var_caller-tract-ceiling/tmp/compound_fix/oracle/`). With that one line
set to the recording run's path, the three files hash to `1109b1a6…`, `1109b1a6…` and
`155fba62…` — the baseline exactly. Both modes agreed with each other (7,459 records, window
coverage identical over 30,700 rows).

**Not run:** `cargo doc --no-deps`, `cargo audit` (no dependency change), `cargo mutants` (not
installed in the image; eight hand-made mutants instead).

**Findings labelled "Needs verification":** 0.

## 4. Open questions and assumptions

1. **Is `bundle_threshold = 0` meant to be reachable on the catalog path?** `classify` refuses it
   with a release `assert!` (`segment_criteria.rs:908`), but typing from the catalog does not call
   `classify`, and nothing I found on that path checks it; `type-regions --flank-bp 0` sets it.
   M2 is fixed either way; the answer decides whether the catalog path should refuse 0 as the scan
   does. Affects M2.
2. **May a catalog be replaced while a run holds it open?** If a rebuild beside a running cohort
   is a supported workflow, Mi1 should be fixed before merge; otherwise it can wait. Affects Mi1.

## 5. Top 3 priorities

All three are fixed in this commit; they are listed so the author can check the fixes.

1. **M1** — `SpanIndex` searched a list that one inverted span had put out of order, and missed a
   requested span the old scan found.
2. **M2** — the absorption sweep left two abutting satellites where the old loop made one, at a
   bundle radius of 0.
3. **M3** — a catalog from another build was refused as *unreadable* instead of *written by
   another version*, the message that says to rebuild it.

## 6. Findings

### Major

#### M1: src/repeat_catalog/mod.rs#L176 — `SpanIndex` binary-searches starts it never checked are ascending — **fixed**

- **Confidence:** High
- **Categories:** reliability (orchestrator), reliability sub-agent (probe)
- **Assumptions:** inverted spans do not reach this code from a BED file: `GenomeRegions` and
  `RegionSet` reject `end < start` at parse time. They can reach it from a caller that builds
  `GenomeRegion`s itself — the type has public fields and no constructor
  ([types.rs](../../../../src/types.rs#L79-L104)) — and several callers do (`ReadScope::Regions(&spans)`
  in `cli/run_ground.rs:405`, `cli/typed_regions.rs:711`, the examples).
- **Problem:** `SpanIndex::new` called a list searchable when every span is on the contig and each
  span ends before the next one starts. That implies ascending starts only if every span has
  `start <= end`. With `[1,2]`, `[30,5]`, `[6,40]` the pairwise test passes, the starts are 1, 30,
  6, and `partition_point(|s| s.start <= to)` on that non-monotone predicate returns 1 for a
  query ending at 20. `candidates(10, 20)` came back empty; the scan it replaced finds `[6,40]`,
  which overlaps `10..=20`.
- **Why it matters:** a requested region silently drops loci, counts and bases — the wrong
  answer without an error that the change promised not to produce.
- **Fix applied:** the shape check also requires `s.start <= s.end` for every span, so such a list
  falls back to the whole-list walk:
  ```rust
  let ordered = spans
      .iter()
      .all(|s| s.contig == contig && s.start.get() <= s.end.get())
      && spans.windows(2).all(|w| w[0].end.get() < w[1].start.get());
  ```
  Test `candidates_hold_every_overlapping_span_when_a_span_is_inverted` fails without it
  (`[6,40] overlaps 10..=20 but the candidates were []`) and passes with it. A second test,
  `candidates_are_the_overlapping_spans_on_ordered_spans`, pins the inclusive edges on a
  well-formed list.

#### M2: src/region_typing/mod.rs#L1054 — the sweep leaves abutting satellites apart at a bundle radius of 0 — **fixed**

- **Confidence:** High (reproduced; the property test with its radius widened to include 0 fails on it)
- **Categories:** reliability sub-agent (counterexample), orchestrator (cause and fix)
- **Assumptions:** reachable only if the catalog path accepts `bundle_threshold = 0` (open
  question 1).
- **Problem:** the old loop ends every absorption with `merge_runs`, which joins satellites that
  overlap **or abut**. The sweep joins a feature to a group only within the radius (`start <=
  reach + bundle_threshold`). For any radius of at least 1, an abutting span is within it, so the
  two agree. At 0 they do not: satellites `22..=86` and `88..=88` with a member covering `87..=88`.
  The member overlaps only the second satellite; the grown `87..=88` abuts the first. The loop's
  `merge_runs` returns one satellite `22..=88`; the sweep returned two, and its own `debug_assert`
  ("not abutting") would panic in a debug build.
  The property test missed it because its generator draws the radius from `1..30`.
- **Why it matters:** two satellite records where there was one, and a panic in debug builds.
- **Fix applied:** `merge_runs(&mut grown)` after the sweep. Merging two abutting hulls adds no
  base, so no feature the sweep left out comes to overlap the result; for radii of 1 or more it
  changes nothing (no two grown satellites abut). The property test's radius now starts at 0, and
  `absorption_by_sweep_joins_satellites_that_abut_after_growing_at_radius_zero` pins the case.
  Both fail without the fix (the property test shrinks to satellites `414..=529`, `531..=531` and
  member `522..=530`, now saved in `proptest-regressions/region_typing/mod.txt`) and pass with it.

#### M3: src/repeat_catalog/reader.rs#L63 — a catalog from another build was refused as unreadable, not as another version — **fixed**

- **Confidence:** High
- **Categories:** errors, extras
- **Problem:** before the change, `open_checking_against_reference` read the header (no page
  index, no column check), then `check_written_by_this_build`; the column, page-index and
  statistics checks ran at the first contig read. `CatalogIndex::load` now runs first and does
  all three, so a catalog from a build with another layout came back
  `Unreadable { … "not a repeat catalog: its columns are …" }` instead of `ToolVersionDiffers`.
  No test pinned the order: `a_catalog_from_another_tool_version_is_refused` writes this build's
  layout.
- **Why it matters:** *written by another version* tells the user to rebuild the catalog;
  *unreadable or truncated* sends them looking for a corrupt file.
- **Fix applied:** when `load` fails, the header is read without layout checks
  (`header_text_without_layout_checks`) and, if present, decoded and version-checked first; only
  then is the layout error returned. Tests
  `a_catalog_from_another_build_with_another_layout_is_refused_by_version` (fails without the fix:
  it got `Unreadable`) and `a_catalog_of_this_build_with_another_layout_is_unreadable`.

#### M4: src/repeat_catalog/parquet_file.rs#L597 — no test read part of a contig through `repeats()` — **test added**

- **Confidence:** High
- **Categories:** extras, reliability
- **Problem:** every `repeats(ReadScope::Regions(..))` call in the tests asked for a missing
  contig or for `1..=u64::MAX`, which takes the new no-filter shortcut. A shortcut that kept too
  many rows (containment written as overlap, an off-by-one at an edge) would pass every typing
  test, because the typing reads are widened precisely so that extra rows change nothing.
- **Fix applied:** `part_of_a_contig_reads_back_the_rows_overlapping_it` compares eleven window
  sets against a scan of every row: the block's exact extent, a window one base inside each edge,
  windows touching a tract at one base, unsorted and overlapping windows, and `118..=1000`, which
  holds the extent's end but not its start. Mutant m3 (containment turned into overlap) fails it
  (§8). No defect was found; the shortcut is correct (the reasoning is in §9).

#### M5: src/repeat_catalog/reader.rs#L134 — new clippy warnings and unformatted code — **fixed**

- **Confidence:** High
- **Problem:** three library warnings introduced by the diff (at `d4027a8e`: `redundant_guards` at
  `reader.rs:119`, `needless_borrow` at `:443` and `:481`), and `rustfmt` differences in
  `parquet_file.rs:530` and `examples/segmenter_timing.rs`. The skill treats a failing
  verification command as Major.
- **Fix applied:** `Some([]) => {}`, `spans` passed without `&`, `rustfmt` on the diff's files.

### Minor

#### Mi1: src/repeat_catalog/parquet_file.rs#L602 — a catalog replaced after open is read with the old footer

- **Confidence:** Medium. **Categories:** unsafe_concurrency.
- **Problem:** the footer is parsed once at open, but each contig read reopens the file by path
  and applies the old footer's page offsets to whatever file is there now. The writer replaces
  the file by rename, so a rebuild during a run gives most likely a decode error, possibly rows
  decoded from the wrong pages. Before, each read parsed the footer of the file it read.
- **Suggested fix:** keep the `File` opened in `CatalogIndex::load` and hand each read a
  `try_clone()`, so every read sees the inode that was checked; or record length and mtime at open
  and refuse a read when they changed. Open question 2.

#### Mi2: src/region_typing/mod.rs#L986 — five `expect`s in `absorb_by_sweep` without a `PANIC-FREE:` note

- **Confidence:** High. **Categories:** errors.
- **Problem:** lines 986, 1005, 1007, 1045, 1047. All are unreachable; the project marks such
  calls with `// PANIC-FREE:` (70 uses).
- **Suggested fix:** pop with `while let Some(group) = open.front().copied()` and a break,
  index the flag vectors in `retain`, and annotate the two `pop_back`s.

#### Mi3: src/repeat_catalog/parquet_file.rs#L579 — `rows_overlapping(path, index, …)` takes a path and an index that must describe the same file

- **Confidence:** High. **Categories:** idiomatic.
- **Suggested fix:** store the path in `CatalogIndex` and make the read a method on it.

#### Mi4: src/repeat_catalog/mod.rs#L96 — `ContigSpans` built with a placeholder `ordered: true` and corrected in a second loop

- **Confidence:** High. **Categories:** idiomatic.
- **Suggested fix:** group into `BTreeMap<ContigId, Vec<GenomeRegion>>`, then build each
  `ContigSpans` with its real flag (a `ContigSpans::new(spans, contig)`).

#### Mi5: src/repeat_catalog/segments.rs#L87 — the cached shape check is recomputed three times per contig, and `contig_walk` groups the spans a second time

- **Confidence:** High. **Categories:** smells.
- **Problem:** `ContigSpans::ordered` is read only by the tally; `segments_of_contig_in`,
  `loci_of_contig_in` and `repeat_features_of_contig` each rebuild a `SpanIndex`
  (`segments.rs:87,163,207`), and `contig_walk` (`reader.rs:229`) builds a `BTreeSet` of the
  contigs the `OwnedScope` map already keys on. Two derivations of one fact can drift.
- **Suggested fix:** crate-private variants taking a `SpanIndex` obtained from
  `scope.index_on(contig)`; `contig_walk(&OwnedScope)`.

#### Mi6: src/region_typing/mod.rs#L914 — `enum From` shadows the standard `From` trait inside `absorb_by_sweep`

- **Confidence:** High. **Categories:** naming.
- **Suggested fix:** `enum SpanOrigin`.

#### Mi7: src/repeat_catalog/parquet_file.rs#L455 — `groups_by_contig` keyed on a bare `u32`, the block extent an unnamed `(u64, u64)`

- **Confidence:** High. **Categories:** naming.
- **Problem:** the extent `(lo, hi)` is compared against windows `(from, to)` of the same shape
  in the containment test at line 599, the one place a transposition would go unnoticed.
- **Suggested fix:** `HashMap<ContigId, Vec<usize>>`; `struct RowGroupExtent { first_start,
  last_end }`.

#### Mi8: src/region_typing/mod.rs#L902 — `absorb_by_sweep` is about 150 lines in a 3,300-line `mod.rs`

- **Confidence:** High. **Categories:** smells, module_structure.
- **Suggested fix:** move the sweep, its oracle and the property test to
  `src/region_typing/absorption.rs`; give the local `Group` methods for joining and merging.

#### Mi9: src/region_typing/mod.rs#L996 — the sweep's satellite-abuts-satellite rule cannot fire, and its input precondition is unchecked

- **Confidence:** High (mutant m6, §8).
- **Problem:** satellites reach the sweep as merged coverage runs, so no two abut; the
  `span.start <= end + 1` branch for a satellite joining a group by its satellite end never
  decides anything. The sweep's correctness also depends on that input shape: the reliability
  probes found that on unsorted or abutting input satellites the sweep and the loop differ (the
  loop returns the input untouched; the sweep sorts and merges).
- **Suggested fix:** `debug_assert!(satellites.windows(2).all(|w| w[0].end + 1 < w[1].start))` on
  entry, and a doc line saying the input is `merge_runs` output.

#### Mi10: src/repeat_catalog/reader.rs — the layout refusals now made at open are untested one by one

- **Confidence:** High. **Categories:** extras.
- **Problem:** a missing page index, missing statistics and non-`u32` statistics are each refused
  at open, and only the foreign-columns case is tested (by M3's tests).
  `a_parquet_file_that_is_not_a_catalog_is_refused` now trips the column check, not the "no header
  in the footer" branch its doc names.
- **Suggested fix:** one test per refusal, writing the catalog schema with
  `set_statistics_enabled(EnabledStatistics::None)` and so on; correct the old test's doc.

### Nits

- Commit subjects `experiment: hot_loops patches …` and `experiment: I/O patch …` should say what
  the change is; the move of the layout checks to open is a behaviour change worth its own line.
- `parquet_file.rs:706`: `row_from_batch`'s doc still says "Shared with the reader"; it is now a
  test helper.
- `parquet_file.rs:451`: `#[derive(Clone)]` on `CatalogIndex` has no user.
- `parquet_file.rs:687,697`: a one-line note that `as u64` reinterprets an INT64 holding a
  `UInt64` (the casts are correct).
- `examples/segmenter_timing.rs:148`: any second argument other than `calling` silently means the
  catalog criteria; the usage line omits the BED argument.
- Inline `std::collections::…` paths in several places; import once per file.

## 7. Out of scope observations

- **The identity oracle's digests depend on the path the reference is given by**, because the VCF
  header records it (`##reference=`). A run from another worktree, or with an absolute path,
  reports three differing digests for identical calls. `scripts/promote_ng_oracle.sh` could strip
  that line as it strips `##commandline=`, or its header could say to pass the reference as
  `tmp/compound_fix/tomato_ref/…`. Separate change.
- `src/region_typing/segment_criteria.rs:1617,1627` fail `cargo fmt --check` on the base commit.

## 8. Missing tests to add now

Added in this commit (each fails without its fix where there is a fix):

| Test | Catches |
|---|---|
| `repeat_catalog::tests::candidates_hold_every_overlapping_span_when_a_span_is_inverted` | M1 |
| `repeat_catalog::tests::candidates_are_the_overlapping_spans_on_ordered_spans` | off-by-one at either inclusive edge of `candidates` (kills m1, m2) |
| `region_typing::tests::absorption_by_sweep_joins_satellites_that_abut_after_growing_at_radius_zero` | M2 |
| `absorption_by_sweep_is_the_iterations_fixed_point`, radius widened to `0..30` | M2 and any radius-0 divergence |
| `repeat_catalog::reader::tests::part_of_a_contig_reads_back_the_rows_overlapping_it` | M4: over- or under-inclusion in the no-filter shortcut and the window search (kills m3, m4) |
| `repeat_catalog::reader::tests::a_catalog_from_another_build_with_another_layout_is_refused_by_version` | M3 |
| `repeat_catalog::reader::tests::a_catalog_of_this_build_with_another_layout_is_unreadable` | M3's fallback must not swallow a real layout error |

Still to add (Mi10): `a_catalog_without_statistics_is_refused_at_open`,
`a_catalog_without_a_page_index_is_refused_at_open`.

**Mutants.** Eight hand-made mutants, each run against every test whose name contains
`repeat_catalog` or `region_typing`, with this review's tests in place:

| Mutant | Change | Result | Killed by |
|---|---|---|---|
| m1 | `candidates`: `s.end < from` → `s.end <= from` | killed | `candidates_are_the_overlapping_spans_on_ordered_spans` |
| m2 | `candidates`: `s.start <= to` → `s.start < to` | killed | the same |
| m3 | no-filter shortcut: window *contains* the extent → window *overlaps* it | killed | `part_of_a_contig_reads_back_the_rows_overlapping_it` |
| m4 | window search: `to < start` → `to <= start` | killed | the same |
| m5 | sweep closes a group one base early (`<` → `<=`) | killed | the property test and the radius-0 test |
| m6 | sweep: a satellite joins by abutting a satellite (`end + 1` → `end`) | survived — **changed no behaviour** | — (see Mi9: merged input satellites never abut, so the branch cannot decide anything) |
| m7 | sweep: a satellite joins a feature's group only strictly inside the radius | killed | the property test and the radius-0 test |
| m8 | `open`: the version check before the layout error removed | killed | `a_catalog_from_another_build_with_another_layout_is_refused_by_version` |

Eight run, none survived that changed behaviour, one changed no behaviour.

m1–m4 and m8 are each killed by one test, and in every case it is a test this review adds; no
test that came with the diff fails on them. The sweep's mutants (m5, m7) are caught by the
diff's own property test.

## 9. What's good

- The sweep keeps the old loop as its test oracle and checks the two against each other on
  20,000 random layouts per run (`region_typing/mod.rs`, `absorption_by_sweep_is_the_iterations_fixed_point`)
  — the property test is what turned up M2 once its radius was widened.
- `SpanIndex::candidates` is *an answer, never a filter*: it returns a superset and every caller
  keeps its own overlap test, so the search can only cost a missed span, never an extra one, and
  an unsearchable list falls back to the whole list (`repeat_catalog/mod.rs`, the `SpanIndex` doc).
- The window search in `rows_overlapping` (`parquet_file.rs:614-625`) is correct for any windows
  whose starts and ends both ascend, including overlapping and inverted ones: the windows not
  ending before a row starts are a suffix, and its first member has the smallest start. The
  no-filter shortcut (`parquet_file.rs:597`) requires the block to hold one contig (statistics
  minimum equals maximum) and a window to *contain* its whole extent, which implies every row
  overlaps.
- The decoded footer is `Arc`-backed in parquet 59.2 (`ArrowReaderMetadata`), so the per-contig
  clone is three reference-count increments, and `new_with_metadata` is exactly what
  `try_new_with_options` did after loading.
- No hash-map order reaches output: `groups_by_contig` is only looked up, and the new maps and sets
  are ordered (`BTreeMap`, `BTreeSet`).

## 10. Commands to re-verify

```
DEV=/Users/jose/devel/pop_var_caller-segmenter-perf-combined/scripts/dev.sh
$DEV cargo test --release --lib          # 4,995 passed; the 3 known release should-panic failures
$DEV cargo test --release --test '*'     # 41 passed
$DEV cargo clippy --release --all-targets
$DEV cargo fmt --check -- src/repeat_catalog/*.rs src/region_typing/mod.rs examples/segmenter_timing.rs
$DEV cargo test --release --lib -- candidates_ absorption_by_sweep part_of_a_contig another_layout
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks $DEV scripts/promote_ng_oracle.sh \
    <tmp>/S_lycopersicum_chromosomes.4.00.fa <tmp>/S_lycopersicum_chromosomes.4.00.fa.repeats.parquet \
    /Users/jose/devel/pop_var_caller/benchmarks/tomato2/regions_n160_200kb.bed <out> \
    <the four CRAMs named in the baseline>
# then compare with the ##reference= line set to the recording run's path (§3)
```

### Author response convention

Address each finding by its identifier (M1, Mi3, …) with `fixed in <commit>` / `disputed because …`
/ `deferred to <issue>` / `won't fix because …`. Answer the open questions in §4 first.
