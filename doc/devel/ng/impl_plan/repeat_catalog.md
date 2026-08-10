# ng repeat catalog — implementation plan (`src/ng/repeat_catalog/`)

**Status:** draft, 2026-08-10. The build order for the reference's tandem-repeat catalog: the criteria
type, the lifted trim/motif/purity helpers, the observer seam, the Parquet file, the builder, the
reader's five query methods, and the `repeat-catalog` subcommand. The design is settled in
[`../spec/repeat_catalog.md`](../spec/repeat_catalog.md) (the *why*) and
[`../arch/repeat_catalog.md`](../arch/repeat_catalog.md) (the types & interfaces). This turns that
design into build order; it is **not** a place for new design — a step that seems to need a decision
goes back to the spec or the arch doc.

It does **not** touch `src/ssr/` (frozen production) and does not depend on `trf-mod`.

---

## Scope

**In:** `src/ng/repeat_catalog/` — `StrRepeatCriteria` and its refusal, `FoundRepeat`,
`RepeatCatalogHeader`, `RepeatCatalogError`, the Parquet writer/reader, `RepeatCatalogBuilder`, the
eight methods of arch §2.3, and the `pop_var_caller_exp repeat-catalog` subcommand. Two edits land
outside the folder: the `ReferenceBasesObserver` seam in `src/ng/reference_info.rs`, and lifting
`finish_locus`'s trim/motif/purity arithmetic into `pub(crate)` helpers in
`src/ng/region_typing/segment_criteria.rs`.

**Out (later plans / homes):**

- **Step 3 reading the catalog instead of scanning** → a follow-up to
  [`typed_regions.md`](typed_regions.md); it needs `partition_windowed` and `TypedRegionIterator`
  first (spec §8).
- **`type-regions`' partition file becoming a view over this file** →
  [`typed_regions_cli.md`](typed_regions_cli.md) §8 (spec §8).
- **The pre-pass consuming `sample_loci_per_stratum`** →
  [`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md)'s own plan; this plan
  ships the method and its tests, not its caller (spec §8).
- **Comparing the detections against production's `src/ssr/catalog/`** → a report, not a spec
  (spec §8).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The criteria-free arithmetic (trim, motif, purity)
  and the row it produces come first; the file format, the seam and the CLI are plumbing that carries
  them.
- **Round-trip before query.** The writer and reader are proven against each other on synthetic rows
  before any query method exists, so a later failure is a query bug and not a serialisation one.
- **Verify against ground truth.** The north star is **equality with `partition_resident`**
  ([`region_typing/mod.rs:380`](../../../../src/ng/region_typing/mod.rs)) — a live scan over the same
  reference at the same policy — not self-consistency.
- **Isolate the silent-failure steps.** A coordinate conversion (0-based half-open → 1-based
  inclusive), a copy floor measured on the wrong span, and a hash-selection rule are all *quietly
  wrong answers*, not panics. Those steps land as their own commits with their oracle green before
  and after.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds:** all `cargo` via `./scripts/dev.sh` (CLAUDE.md).

## Preconditions (already in place)

- **The scanner:** `find_tandem_repeats` ([`tandem_repeat.rs:483`](../../../../src/ng/tandem_repeat.rs)),
  `ScanParams` (`:123`), `PeriodRange` (`:55`), `RepeatInterval` (`:208`).
- **Classification:** `prefilter` ([`segment_criteria.rs:677`](../../../../src/ng/region_typing/segment_criteria.rs)),
  `classify` (`:849`), `finish_locus` (`:1011`), `minimal_trim` (`:1217`), `MinCopies` (`:355`),
  `SsrSegmentCriteria` (`:478`), `SsrSegment` (`:141`).
- **The reference pass:** `read_reference_info` ([`reference_info.rs:270`](../../../../src/ng/reference_info.rs)),
  `FastaPass` (`:516`), `ContigInfo` (`:51`), `write_fai` (`:816`).
- **The oracle:** `partition_resident` ([`region_typing/mod.rs:380`](../../../../src/ng/region_typing/mod.rs))
  and the fixture reference `tests/data/tandem_repeat/synthetic_ref.fa`.
- **The exp binary:** `PopVarCallerExpCommand` ([`pop_var_caller_exp/cli.rs:22`](../../../../src/pop_var_caller_exp/cli.rs))
  and `typed_regions.rs` as the module shape to mirror.
- **Not yet present:** the `parquet` crate. A1 adds it.

---

## The steps

### Milestone A — the types and the criteria-free arithmetic

**A1. Module scaffold, dependency, types.**  ✅
`pub mod repeat_catalog;` in `ng/mod.rs`; the folder of arch §*Module home*; `FoundRepeat`, `SpanBp`,
`RepeatCatalogHeader`, and the `#[non_exhaustive]` `RepeatCatalogError` with every variant and its doc
comment. Add `parquet` to `Cargo.toml` (no code using it yet). Nouns and errors only, no logic.
*Source:* arch §1.2–§1.4.

**A2. `StrRepeatCriteria` + `serves`.**  ✅
The wrapper over `SsrSegmentCriteria`, the named `pub const` defaults (`[5, 5, 4, 4, 4, 3]`, 15 bp,
500 bp, periods 1..=6), `CriteriaRefusal`, and `serves`. Tests: equal criteria serve; a lower copy
floor at one period refuses **naming that period and both numbers**; a lower flank refuses; a wider
period range refuses; a different purity floor, score floor, satellite cap or bundle radius **serves**
— the mirror case, which is what §4.2 turns on. *Depends:* A1. *Source:* arch §1.1, spec §4.1–§4.3.

**A3. Lift trim, motif and purity out of `finish_locus`.**  ✅  **Own commit, do not bundle.**
`pub(crate)` helpers in `segment_criteria.rs` for the whole-motif trim, the motif slice and the purity
recomputation; `finish_locus` calls them and keeps its behaviour. **Silent** (a changed trim is a
changed locus everywhere, not a panic): the existing `region_typing` test suite must be green before
and after, unchanged. *Depends:* A1. *Source:* arch §5 (row *motif / trim / purity*), spec §3.2, §7.

**A4. `FoundRepeat` from a `RepeatInterval`.**  ✅  **Own commit, do not bundle.**
The conversion at the builder's edge: 0-based half-open → 1-based inclusive, the copy floor **on the
detected span** as `prefilter` measures it, the trim (absent when there is no clean cut), motif,
purity over the trimmed tract, and the 15 bp flank floor against the contig's length. **Silent** (an
off-by-one here is a wrong locus in every consumer). Tests: a hand-checked tract at known coordinates;
a tract with no clean trim keeps its row with `trimmed = None`; a tract 14 bp from a contig end is
dropped and one at 15 bp is kept; a tract whose *trimmed* count falls below the floor but whose
detected count clears it **is kept** (the bundling-preservation rule). *Depends:* A2, A3.
*Source:* arch §2.2, spec §3.1.

> **Checkpoint A:** a repeat interval plus its bases becomes a row, and a policy can be refused with
> the axis named. Pause for review.

### Milestone B — the file

**B1. The Parquet schema and writer.**  ✅
`parquet_file.rs`: the seven columns of arch §3 with their types and encodings, one row group per
contig, the header as footer key-value metadata, atomic write (`.tmp` + rename). Fixed codec, level
and writer-version string. Tests: a file written from synthetic rows opens; its footer metadata
round-trips; two writes of the same rows are **byte-identical**. *Depends:* A1. *Source:* arch §3,
spec §3.5, §6.

**B2. The reader: header and rows.**  ✅
`open_checking_against_reference` (contig table, order, per-contig MD5s, scoring weights, tool
version), `header()`, `contigs()`, and `repeats_in_region` streaming a row group at a time. Tests:
round-trip of B1's rows; **a missing file and a mismatched reference are different errors**, the first
naming the command, the second naming the contig; a reordered contig table with matching lengths is
caught; a file truncated mid-write fails to open. *Depends:* B1. *Source:* arch §2.3, §1.4, spec §4.3,
§10.7.

> **Checkpoint B:** rows survive a write/read round trip, and a file that does not describe this
> reference cannot be opened. Pause for review.

### Milestone C — the builder and the seam

**C1. `ReferenceBasesObserver` in `reference_info.rs`.**  ☐
The trait, `read_reference_info_observing`, and `read_reference_info` reduced to a call with a no-op
observer. `reference_info` gains no import from `repeat_catalog`. Tests: the observer sees
`contig_started → bases* → contig_finished` in file order; the bases it receives, concatenated and
uppercased, equal the contig's; a `.fai`-only read calls nothing; the existing `reference_info` suite
is unchanged and green. *Depends:* A1. *Source:* arch §2.1, spec §2.2.

**C2. `RepeatCatalogBuilder` — sequential.**  ☐
The observer impl: accumulate a contig, `find_tandem_repeats` over it whole, rows through A4, one row
group per contig, `finish` returning the per-period tally. Tests on the fixture reference: a 2 kb
tract comes out as **one row**; repeats at both contig ends behave as A4 specifies; the row order is
(contig, start, period, end). *Depends:* A4, B1, C1. *Source:* arch §2.2, spec §2.3, §10.6.

**C3. `--threads`: contigs in flight, rows in reference order.**  ☐
Up to N contigs scanned at once, completed contigs written in reference order. Tests: the file is
**byte-identical** at 1, 2 and 4 threads, on a fixture whose contigs differ enough in size that they
finish out of order. *Depends:* C2. *Source:* spec §2.4, §10.5.

> **Checkpoint C:** the catalog builds from a reference in one pass, deterministically, at any thread
> count. Pause for review.

### Milestone D — the queries

**D1. `genome_segments` — the segmentation.**  ☐  **Own commit, do not bundle.**
Rows → `prefilter` → bundling → `classify`'s admission, over the stored spans, yielding `TypedRegion`s
that cover the region with no gap. **Silent** (a wrong segmentation is a wrong genotype downstream).
Its oracle is D2, which lands with it. *Depends:* A4, B2. *Source:* arch §2.3, spec §5.1.

**D2. The differential against `partition_resident`.**  ☐  *(lands with D1 — it is D1's guard.)*
Build a catalog over the fixture reference; derive the segmentation; compare against
`partition_resident` at the same policy — regions, kinds, coordinates, motifs, purities identical.
**Run it at several policies**, including ones differing from the build settings on every bounded axis,
and on a fixture carrying overlapping detections (a tract detected at two primitive periods, and two
intersecting tracts). The one stated exception is tracts closer than 15 bp to a contig end, so the
comparison runs at a reader flank of 15 bp or more. *Depends:* D1. *Source:* spec §5.1, §10.1, §9.4.

**D3. `str_loci` and `count_loci_per_stratum`.**  ☐
The loci alone, and the per-stratum tally keyed by (period, **trimmed** span / period). Tests: the
tally equals a count taken by scanning the fixture directly; `str_loci` equals `genome_segments`
filtered to `RegionKind::SsrSegment`; both refuse an under-permissive policy **eagerly**, before any
row is read. *Depends:* D1. *Source:* arch §2.3, spec §5.3, §5.4, §10.3.

**D4. `sample_loci_per_stratum`.**  ☐  **Own commit, do not bundle.**
The `cap` lowest `hash(contig, start, seed)` per stratum, from a bounded heap, returning the counts
from the same pass. **Silent** (a biased sample is a biased parameter estimate, and nothing reports
it). Tests: a stratum with fewer than `cap` loci keeps all of them; the same seed gives the same set
across runs; a different seed gives a different set; the selection is **order-independent** — shuffling
the input order, and merging two halves by taking the lowest `cap` of the union, both give the
identical set. *Depends:* D3. *Source:* arch §2.3, spec §5.3.

> **Checkpoint D:** the file answers every question the pre-pass and step 3 will ask, and its
> segmentation is indistinguishable from a live scan. Pause for review.

### Milestone E — the command

**E1. `pop_var_caller_exp repeat-catalog`.**  ☐
`RepeatCatalogArgs`, `run_repeat_catalog`, `RepeatCatalogCliError`, and the new
`PopVarCallerExpCommand` variant — mirroring `typed_regions.rs`. `--reference` required; `--output`
defaulting to a sibling of the FASTA; `--threads`; the criteria knobs; `--force` guarding an existing
file; the per-period tally printed on completion. Tests: a build writes the file and prints the tally;
a second run without `--force` refuses and leaves the file untouched; a misspelled knob exits non-zero.
*Depends:* C3, D3. *Source:* spec §2.6.

**E2. Measure, and record it.**  ☐
On the tomato reference and one human one: file size, row count by period, wall clock with and without
the scan, and peak RSS. This is what settles spec §9.1, and it is a report under
`ia/reports/implementations/`, not a spec edit. *Depends:* E1. *Source:* spec §10.4, §9.1.

> **Checkpoint E:** the catalog is buildable from the command line and its cost is measured rather
> than assumed. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests: `serves` refuses each bounded axis naming both values and serves each unbounded one; the row conversion on hand-checked coordinates, including no-clean-trim, the 14/15 bp flank pair, and trimmed-below-floor-but-kept; the `region_typing` suite unchanged across the lift |
| B | write→read round trip; byte-identical repeat writes; missing-file vs wrong-reference are distinct named errors; reordered contigs caught; truncated file fails to open |
| C | one row for a 2 kb tract; contig-end cases; **byte-identical output at 1, 2 and 4 threads** on out-of-order contigs; the reference pass's own suite green |
| D | **equality with `partition_resident`** at several policies, with overlapping-detection fixtures; strata tally vs a direct count; eager refusal; sampling order-independence and seed-stability |
| E | CLI writes the file and prints the tally; `--force` guard; measured size, rows by period, wall clock and peak RSS on tomato and human |

## Out of scope (next plans)

- **Step 3 reading the catalog instead of scanning** → a follow-up to [`typed_regions.md`](typed_regions.md),
  once `partition_windowed` and `TypedRegionIterator` exist and D2 is green (spec §8).
- **The `type-regions` partition as a view over this file** → [`typed_regions_cli.md`](typed_regions_cli.md) §8.
- **The pre-pass's use of the sample** → the joint-loci plan (spec §8).
- **A region subset defined by a `--regions` BED at the consumer** → the consumer that first needs it
  (spec §8); the `region` argument exists here, its BED plumbing does not.
