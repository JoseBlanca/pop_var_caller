# Fixes applied: ng_repeat_catalog review of 2026-08-11

**Report reviewed:** [ng_repeat_catalog_2026-08-11.md](ng_repeat_catalog_2026-08-11.md)
**Applied against:** `8a299f2a` on `ng-cohort-stutter-readgroups`
**Validation:** in the dev container — `cargo fmt` clean; `cargo clippy --lib --tests --all-features -- -D warnings` clean; `cargo test --tests` **lib 3,307 passed, 0 failed**, every integration suite green, the differential **10 passed**.

## What changed, in one paragraph

**Three refusals the design calls fatal now happen**: a catalog written by another build of the detector, a region naming a contig the catalog does not hold, and a parquet file that carries the catalog's header key but somebody else's columns. **One wrong answer is fixed**: a satellite crossing a requested edge came back clipped, where a live scan hands it back whole. And **five items nothing reached** are gone, which was only invisible because every submodule was `pub`, so `dead_code` could not see the folder at all.

## Findings table

| ID | Title | Decision | Status |
|---|---|---|---|
| B1 | Tool version and scoring weights never compared | Apply (split) | **Applied with adaptation** |
| B2 | Region-scoped tally is short | **Ask** | **Open — needs a decision** |
| B3 | Satellite clipped at a requested edge | Apply | **Applied** |
| M1 | Missing row-group statistics read as an empty contig | Apply | **Applied** |
| M2 | Foreign schema panics | Apply | **Applied** |
| M3 | `reach_into_window` unchecked addition | Apply | **Applied** |
| M4 | `--min-flank-bp` below the tally's assumption | Apply (split) | **Applied with adaptation** |
| M5 | `--max-period` above the motif limit | Apply | **Applied** |
| M6 | `trimmed`/`purity` co-dependent options | Defer | **Deferred** |
| M7 | `--threads` does not size rayon's pool | Defer | **Deferred** |
| M8 | Window widening unpinned | Apply (partial) | **Applied with adaptation** |
| M9 | Two `..Default::default()` tails | Apply | **Applied** |
| M10 | `decode_header`'s `>=` guards | Apply | **Applied** |
| M11 | `admit`'s score gate unpinned | Defer | **Deferred** |
| M12 | Non-hex digest decodes as absent | Apply | **Applied** |
| M13 | Unknown contig silently ignored | Apply | **Applied** |
| M14 | Five items with no caller | Apply | **Applied** |
| M15 | Arrow leaves `parquet_file.rs` | Apply | **Applied** |
| M16 | Every submodule `pub mod` | Apply | **Applied** |
| M17 | Region read does no page pruning | Defer | **Deferred** |
| M18 | `finish` promises a check it does not perform | Apply | **Applied** |
| Minors, Nits | see the review's section 6 | Mixed | see below |

## Applied, with what actually changed

**B1 — and it is two findings, not one.** The tool-version half was a real defect and is fixed: `open_checking_against_reference` now compares the header against `env!("CARGO_PKG_VERSION")` and refuses, naming both. The builder no longer *takes* a version to stamp — it stamps the running one — because a caller that could choose it could stamp a catalog with a version that never wrote it, and the refusal would then be checking a claim rather than a fact.

The scoring-weights half is **disputed as a defect and fixed as a doc**. `check_scored_with` has no caller because it is for a caller that scans as well as reads, and no such caller exists yet: a reader passing `StrRepeatCriteria` has no weights to offer. What was wrong is `check_serves`'s doc, which claimed to check them. Corrected to say which question it answers and which it does not.

**B3.** `RegionKind::Satellite` joins loci and bundles in the emitted-whole arm, matching the walk's own `clips_at_a_bed_edge`. Pinned by `a_feature_outside_the_window_still_decides_what_is_inside_it`; removing the arm fails it.

**M1, M2, M12.** Three ways a broken file used to read as a valid short one: a row group with no usable statistics now refuses instead of reporting an empty contig, the schema is compared against `catalog_schema()` before any column is addressed by position (it used to panic inside parquet), and a digest field that is present and will not decode is `Unreadable` rather than "no digest recorded" — which used to skip the check for exactly the contig whose record was damaged.

**M3.** `saturating_add`, with the reason in the comment: both terms come from a file this build did not necessarily write.

**M4 — split.** The panic is fixed by making the assumption true rather than by weakening it: the clamp can only fire for a tract abutting a contig's first or last base, so **one** base of flank floor is enough, and the command now refuses `--min-flank-bp 0`. The tally's comment says one base rather than fifteen. The other half — that `min_flank_bp` passes `serves` and is never applied at read time — is **disputed**: applying it at read time would make the catalog drop loci a live scan keeps, and the differential would fail. It is a build-time axis, and `serves` compares it as one.

**M5.** `--max-period` above `MAX_MOTIF_LEN` is refused, naming both numbers.

**M8 — partial, and the shortfall is stated.** The reach widening is now pinned: deleting it fails the new test. The `longest_tract_bp` term inside it is **still unpinned** — a window I placed past a 1.2 kb array came back generic, so the fixture cannot yet reach the case where a row outside the window changes the answer inside it. Recorded as a follow-up rather than claimed.

**M9.** Both `..Default::default()` tails are gone; every field is named, with a comment on the two the catalog has no use for.

**M10.** Every header line's field count is checked exactly, with the count in the message. Appending a field to a line is now an error rather than a silent truncation.

**M13.** A region naming a contig outside the header refuses from all three read surfaces.

**M14, M15, M16.** Seven of the eight submodules are `pub(crate) mod` — which needed no import change anywhere and immediately made `dead_code` fire. `RegionSegments` and its private feeder, `segments_of_rows` and `locus_span` are deleted; `clears_detected_copy_floor` is kept but is now *called* by `row_for_interval`, so the catalog has one copy-floor test rather than two and the named property has a caller its doc can honestly claim. `row_from_batch` is private, `decode_header` and `rows_overlapping` are `pub(crate)`, so no Arrow type is on the public surface.

**M18.** The check its doc promised exists: the builder counts the contigs the pass handed it and refuses to write a header naming more. Counted explicitly, not inferred from `longest_tract`, which a contig holding no repeat never extends.

**Minors applied:** the "seven columns" doc against a nine-column schema (twice); `ReadScope::Regions`' doc, which said an unknown contig is ignored.

## Deferred, with the reason

- **M6 (`trimmed`/`purity` co-dependent options).** The right fix is one `Option<TrimmedTract { span, purity }>`, and the review's sub-agent proved it compiles green. It changes `FoundRepeat`, which is the type every consumer of the catalog will destructure, and the consumer wiring starts next. Better in that commit than in this one.
- **M7 (`--threads` does not size rayon's pool).** A behaviour change to a knob whose memory cost is measured in gigabytes per thread. Needs a decision about whether the flag sizes a pool or is renamed to what it bounds.
- **M11 (`admit`'s score gate unpinned) and the tail of M8.** Both need a fixture that reaches a case the current one cannot. Grouped as one follow-up, because they need the same work.
- **M17 (page pruning).** `RowFilter` does not prune its own predicate columns; the fix is `with_row_selection`, which is a rewrite of the read path plus a bench to hold it. The measured cost today is 18 ms per read on tomato chromosome 1, so nothing is on fire.
- **The Minor and Nit naming items** — `ContigTally`'s name, `Str…` versus `Ssr…`, `class`, `f`, `wanted`/`spans`/`windows`, the bare column indices. Real, and a rename sweep is its own commit with its own diff to read; mixing it into a correctness commit would bury the correctness.

## Follow-ups this leaves

1. **B2 needs the owner's answer** before anything can be applied (see the review's open question 1).
2. A fixture that puts a feature genuinely outside a read window, to finish M8 and M11.
3. The naming sweep, and M6 with the consumer wiring.
4. `segment_criteria.rs`'s stale gate archaeology — out of scope for this module, recorded in the review's section 7.
