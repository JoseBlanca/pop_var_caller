# ng repeat catalog — A1 + A2: the module, the row types, the criteria

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **A1** and **A2**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §1, §3.1, §4.1–§4.3 and
[`arch/repeat_catalog.md`](../../ng/arch/repeat_catalog.md) §1.*

## Plan

A1 (module scaffold, dependency, row/header/error types) and A2 (`StrRepeatCriteria` and its
`serves` comparison) ran as **one loop iteration**, named here and in the commit. A1 alone is pure
scaffold, and `RepeatCatalogHeader` carries a `StrRepeatCriteria` — splitting them would have landed a
header missing the field that makes it useful, then edited it back a commit later. The plan's own
step-writing convention allows a pure-scaffold step to pair with its neighbour when it is said out
loud.

## Assumptions and deviations (all minor, all recorded)

1. **`SpanBp` is named `TractSpan`.** The value is a repeat tract's span, and the arch doc's name said
   only "a span measured in base pairs". No contract changed.
2. **`CriteriaRefusal` lives in `criteria.rs`, not `mod.rs`.** It is the failure of `serves`, so it
   belongs beside it; `mod.rs` re-exports it.
3. **Two facts in the design docs were wrong and are corrected in place** — a citation, not a design
   change: ng's STR locus generator fetches a **15 bp** flank by default
   (`SsrGeneratorConfig::flank_bp` = `DEFAULT_BUNDLE_THRESHOLD` = 15,
   [`locus_generation/ssr.rs:132`](../../../../src/ng/locus_generation/ssr.rs#L132),
   [`segment_criteria.rs:299`](../../../../src/ng/region_typing/segment_criteria.rs#L299)), not 30 as
   spec §4.1 and the §5.2 code block said. The catalog's 15 bp floor therefore **equals** the caller's
   default flank rather than sitting below it: a caller at the default is served exactly, one asking
   for more is served by arithmetic on stored coordinates, and only one asking for less is refused.
   The rule is unchanged; the number it was justified against was.
4. **`parquet` is added with `default-features = false, features = ["arrow", "zstd"]`.** The spec fixes
   the codec, so `snap`/`brotli`/`lz4`/`flate2`/`simdutf8` would be dead weight.

## Changes made

- [`src/ng/repeat_catalog/mod.rs`](../../../../src/ng/repeat_catalog/mod.rs) — module doc (why the
  file exists: the pre-pass's per-stratum sample needs the genome enumerated), `TractSpan` with
  `len_bp` / `repeat_count`, `FoundRepeat` with both spans and `stratum()`, `RepeatCatalogHeader`,
  and the `#[non_exhaustive]` `RepeatCatalogError` with eight variants, each documenting when it
  fires.
- [`src/ng/repeat_catalog/criteria.rs`](../../../../src/ng/repeat_catalog/criteria.rs) —
  `StrRepeatCriteria` (wrapping `SsrSegmentCriteria`, plus `min_flank_bp` and `max_str_len_bp`), the
  catalog's defaults as named `pub const`s, `serves`, and `CriteriaRefusal`.
- [`src/ng/mod.rs`](../../../../src/ng/mod.rs) — `pub mod repeat_catalog;`.
- [`Cargo.toml`](../../../../Cargo.toml) — the `parquet` dependency, with a comment saying what it is
  for and why the default features are off.

**No logic beyond `serves`**: the builder, the reader and the Parquet writer are B and C.

## Tests added

13 unit tests, all green.

- **Spans and strata** (4): an inclusive span's length; whole-copy counting; a repeat with no clean
  trim belongs to no stratum; **the stratum counts the trimmed tract, not the detected one** — the
  rule spec §3.1 added, and the one a future edit could quietly invert.
- **`serves`** (9): exact match and a stricter reader are served; a lower copy floor refuses **naming
  the period and both numbers**; a lower flank refuses; a reader keeping every non-empty flank (step
  3's rule today) refuses; a wider period range refuses at either end; **the four unbounded axes never
  refuse** — the mirror case §4.2 turns on; a floor below the table but outside the built period range
  does not refuse on its own.
- **A guard on the two floor tables**: `the_catalog_floors_are_below_the_calling_floors_everywhere`
  fails if either table moves such that the catalog stops leaving room below the calling floors.

## Validation

Run in the dev container (`./scripts/dev.sh`):

- `cargo fmt` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib repeat_catalog` — **13 passed, 0 failed**.

**One pre-existing failure, unrelated:** `cargo clippy --all-targets` fails on
`examples/ng_inbreeding_harness.rs` (three `dead_code` errors on `RunsFit` fields), committed in
`076cb5e9` and untouched by this work. Reported rather than fixed — it belongs to the parameter
pre-pass work in flight on this branch.

## Tradeoffs and follow-ups

- `TractSpan`'s fields are public and unconstrained; `len_bp` `debug_assert`s `start <= end` rather
  than making the type impossible to build reversed. The check that matters lands in **A4**, where a
  span is built from detector output — the same place `SsrSegment::new` checks its own.
- `StrRepeatCriteria::default()` inherits `SsrSegmentCriteria::default()`'s purity floor, score floor
  and bundle radius. The builder applies none of them, so they are provenance only.
