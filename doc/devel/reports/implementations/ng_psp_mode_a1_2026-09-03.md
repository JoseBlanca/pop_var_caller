# A1 — the psp header carries the run's identity checks: analysed regions and segmentation inputs

**Date:** 2026-09-03. **Plan:** [`run_driver_psp_mode.md`](../../ng/impl_plan/run_driver_psp_mode.md),
Milestone A step A1. **Branch:** `ng-psp-mode`. **Design authority:**
[`run_streaming.md`](../../ng/spec/run_streaming.md) §6.1–§6.2 (what the header records and what a
run refuses on), [`arch/run_streaming.md`](../../ng/arch/run_streaming.md) §4 (the typed field),
[`psp_file_format.md`](../../ng/spec/psp_file_format.md) §3.1 (the header stays plain text).

---

## What landed

**`Header` gains one typed field, `segmentation_inputs: SegmentationInputs`** — the arch §4 shape:
the analysed regions (the cross-cohort check), and the repeat catalog's identity plus the
repeat-tract routing criteria (the file-against-run check). `SegmentationInputs` is the operand
type `first_difference` already compares ([`segments.rs`](../../../src/ng/run/segments.rs)), so
Milestone E's refusals will compare exactly what the file records, with no translation layer to
disagree with it.

**`format_version` stays `(1, 0)`** precisely because no written psp predates these fields: the
store has no callers outside tests and examples, so a required header section today is free, where
after `generate-psps` ships it would be a version negotiation (the plan's "format before wiring"
principle).

### The wire shape (new file, [`segmentation_section.rs`](../../../src/ng/psp/segmentation_section.rs))

One TOML section, three parts, all typed:

- `[[segmentation.analysed-region]]` — one span per row: contig **by name**, `start`/`end`
  1-based inclusive (the pipeline's own convention, `src/regions.rs`). Names anchor back to the
  header's own contig list on decode, so a span cannot quietly index a contig the file does not
  declare.
- `[segmentation.repeat-tract-criteria]` — the routing criteria the walk asked the catalog with:
  period range, per-period copy floors plus the beyond-the-table fallback, purity floor, score
  floor, bundle radius, flank floor, satellite cap.
- `[segmentation.catalog]` — the catalog's header whole: reference digest, tool version,
  `[segmentation.catalog.scan]`, `[segmentation.catalog.built-under]`, and one
  `[[segmentation.catalog.contig]]` row per contig with its longest stored tract **beside it**
  rather than in a parallel list, so the two cannot go different lengths in the file. Recorded
  whole rather than digested because a refusal must name the field that differs (spec §6.1); a
  digest can only say two things disagree.

Decode rebuilds through the same checked constructors the run itself uses — `PeriodRange::new`,
`MinCopies::new`, and a new `RegionSet::from_genomic_order_spans` — so a value the run could not
build cannot arrive from a file either.

### Reused-API adaptations (implementer latitude, recorded)

- **`GenomeRegions` gained a constructor** (`from_normalized_spans`,
  [`region_typing/mod.rs`](../../../src/ng/region_typing/mod.rs)): the only existing constructors
  were `whole_contigs` and BED parsing, and the decode path needs to rebuild a set from recorded
  spans. It **refuses** a list that is not already normalized (out of order, overlapping, or
  touching spans) rather than re-sorting it, because a silently re-normalized set would compare
  unequal to the set that was recorded — in a check whose whole job is that comparison.
  *(As first written this landed in `src/regions.rs`; the review's M5 caught that the file is
  frozen production code (ruling 2026-07-16), so the constructor moved into ng's own wrapper,
  which now owns its span storage and builds through `RegionSet`'s public API —
  `src/regions.rs` is byte-identical to production again.)*
- **The contig-list rules split out of `check_rules` into `check_contigs`**, called once from
  `check_rules` and once by the reader *before* span resolution: resolving spans against a
  duplicated or zero-length contig would otherwise report the span as broken when the contig is.
  One function, called twice, cannot disagree with itself.
- **`BrokenRule`, `digest_of`, `hex_of`, `MAX_TOML_INTEGER` widened to `pub(crate)`** so the
  section module reports rules in the header's own shape.
- **The section serializes before the manifest**, so the manifest's field declarations still close
  the body and the existing cut-at-`[[manifest.field]]` test helper (and anyone reading with
  `head`) keeps every other section intact.
- **A purity floor crosses TOML as its shortest decimal** (`wire_float_of`): the header shows
  `min-purity = 0.93`, not `0.9300000071525574`, and the trip is exact both ways (asserted by
  `a_purity_floor_stays_exact_and_short_through_the_wire_float`).

### What the plan asked that this step does *not* do

Nothing of A2–A4: no read-group table, no `max_record_span`, no read filters in provenance. No
consumer reads the new field yet — E1 is its first reader, exactly as planned.

## Tests

- **Round-trip, field for field**, at the section level
  (`the_section_round_trips_field_for_field`) and through the whole header path — every existing
  header round-trip test now carries the section, because the shared fixtures
  (`a_written_header`, `writer.rs::tests_support::a_header`, `mod.rs`, `record.rs`) embed a
  non-default `SegmentationInputs` (proper sub-spans, purity 0.93, non-default catalog version) so
  a round trip that substituted a default fails rather than passing by coincidence.
- **Byte-order independence** — the existing
  `the_same_header_encodes_to_the_same_bytes_whatever_order_it_was_built_in` now runs over a
  header carrying the section; the section itself holds no map, only ordered lists.
- **Refusals, both sides.** Section-level: an undeclared contig in a span, a copy-floor table of
  the wrong width, an empty period range, a NaN purity, unequal catalog lists, a span past its
  contig, a span off the list (`segmentation_section.rs` tests). Whole-path: writer and reader
  alike refuse a broken section (`a_broken_segmentation_section_is_refused_by_writer_and_reader_alike`).
  Constructor-level: out-of-order / overlapping / touching / out-of-bounds spans
  (`regions.rs` tests, including the `u32::MAX` overflow edge).
- **The readable-body test** now pins the section's spelling: `[[segmentation.analysed-region]]`,
  `min-purity = 0.93`, `[segmentation.catalog.scan]`, `[segmentation.catalog.built-under]`,
  `[[segmentation.catalog.contig]]`.

## Validation

In the dev container (Apple `container`, 2026-09-03):

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --all-targets --all-features` — **all 16 test binaries pass; the lib suite is
  6,050 passed, 0 failed** after the review fixes (the psp module alone: 415, via
  `cargo test --lib 'ng::psp'`).

## A size fact worth carrying to Checkpoint A

The worst-case header test (`a_reference_of_thirty_thousand_scaffolds_still_fits_the_header`) now
rebuilds the section for all 30,000 scaffolds — an analysed span and a catalog row per scaffold —
and prints what it measured: **10,798,518 bytes of the 16,777,187-byte body ceiling** (the test's
own `eprintln`, re-run 2026-09-03 after the review fixes; before this step the fixture only had to
clear 1 MiB). Linearly, digest-carrying assemblies of about 46,000 scaffolds still fit; a
100,000-scaffold draft would be refused loudly at write time by the existing ceiling check. Not a
defect — the refusal is the designed behaviour — but the headroom halved, and the ceiling is one
constant if it ever binds.

## Deviations from the plan

None of substance. The plan's Preconditions cite `src/ng/cohort_merge/observation_cache.rs`; the
real path is `src/ng/run/cohort_merge/observation_cache.rs` (all cited line numbers hold there).
