# ng repeat catalog — the live scan is removed

*Implementation report, 2026-08-11. Completes the migration begun in
[`ng-repeat-catalog-tally_2026-08-11.md`](ng-repeat-catalog-tally_2026-08-11.md). Design:
[`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) and
[`spec/typed_regions.md`](../../ng/spec/typed_regions.md).*

## What changed

ng used to find tandem repeats by scanning the reference during every run. It reads them from a
Parquet catalog built once per reference now, and **the scanning machinery is gone** (owner,
2026-08-11): `TypedRegionIterator`, `partition_windowed`, `BlockWalk`, `SpanWalk`, `WindowTally`,
`ScanSpan` and `TypedRegionError`, about 2,000 lines from
[`region_typing/mod.rs`](../../../../src/ng/region_typing/mod.rs).

That module still owns everything after detection — the region vocabulary, the admission policy, the
satellite cap, the bundling, the generic fill and the tally — all of which the catalog's reader
calls rather than copies. What it no longer owns is any way of *finding* a repeat except from a
contig's bases.

## The last call site, and how it was checked

`parameter_estimation::generic::real_alignments` — four `#[ignore]`d identities that walk a real
alignment file through the real locus generator — was the fourteenth and last consumer. It now opens
the catalog beside the reference.

Verified against `bb27dd42`, the last commit before any consumer moved, on tomato SRR7279481 over
`benchmarks/tomato1/regions.bed` (8.0 Mb in 80 spans, one library). All four identities pass at both
commits and every printed number is the same:

| | value |
|---|---|
| typed regions | 8,845 |
| generic loci | 7,623,391 |
| positions covered | 7,628,517 |
| reads | 78,459,630 |
| occupied cells | 323 |
| fitted error rate | 4.591 × 10⁻³ |
| heterozygosity | 5.59 × 10⁻⁴ |

## Keeping the catalog checkable: `partition_resident_in`

The claim the whole design rests on is *the file gives the same segmentation the bases do*, and a
claim like that needs a second implementation to be checked against. `partition_resident` was that
implementation already, but it answered only for a whole contig and produced no tally — and the
three differential tests that needed a region subset or a tally used the walk instead.

So `partition_resident` gained a sibling. **`partition_resident_in(chrom, contig, bases, config,
wanted)` returns the typed regions inside `wanted` and what they contained**, mirroring the reader's
`segments_of_contig_in`. The two paths differ where it matters — one detects repeats from the bases,
the other reads rows from a file — and share everything after detection, so a disagreement points at
the file rather than at the policy.

### It kills every mutation tried against it

Six deliberate breaks, each run against the region-typing unit tests, the differential, and the
anchor:

| mutation | killed by |
|---|---|
| the requested-span count dropped | region-typing unit tests |
| a rejection charged wherever it is, not only inside the request | the differential |
| the repeat coverage not clipped to the request | the differential |
| the coverage summed per tract instead of merged first | region-typing unit tests |
| a locus clipped at a requested edge like a generic stretch | the differential |
| the whole contig emitted, ignoring the request | region-typing unit tests |

Note which suite catches what: **three of the six are invisible to the unit tests and only the
differential sees them.** A run of the mutations against the unit tests alone reported three
survivors, and they were survivors of a mis-aimed test filter rather than of the code.

## The port anchor moves with its subject

`region_typing/anchor.rs` drove the shipping stack of a scan: a real multi-contig FASTA on disk,
through a file-backed evicting reference, through the walk. It is
[`repeat_catalog/anchor.rs`](../../../../src/ng/repeat_catalog/anchor.rs) now and drives the stack
that actually ships: the same FASTA, through the reference pass and the builder into Parquet, back
out through the reader.

It keeps `.cat` parity against the committed trf-mod-built golden catalog, the partition invariant,
region-subset invariance and the contig-edge cases. Parity is **strict**: all 16 golden loci are
matched, and the file gives 17 loci over the same sequence. An earlier draft excused a golden locus
lost to the file's 15 bp flank floor; measured, that exemption fires 0 times, so it was a branch
nothing took hiding a case that ought to fail loudly.

**What it loses is window-invariance.** The builder scans a contig whole
([`builder.rs`](../../../../src/ng/repeat_catalog/builder.rs)), so there is no window for the answer
to depend on — a satellite of any length comes out as one row, and the three things a windowed scan
needs (a margin carried across each chunk, a rule for which side a straddling detection belongs to,
a cap on the repeat length it can promise to catch whole) do not exist here.

## `type-regions` loses a flag, and one line of output

`--window-bp` sized the walk's memory unit. With no walk it configures nothing, so it is gone, and
with it the `## window_bp:` line from the output header — a header exists to say how a file was
produced, and recording an inert knob says something false.

Measured on tomato SL4.0ch01's first 2 Mb: 5,046 lines before, 5,045 after, and the diff is exactly
that one header line.

The `--max-str-len` / `--flank-bp` cross-flag check stays. It was the walk's error
(`TypedRegionError::MarginNarrowerThanFlank`) and is now the CLI's own
(`TypedRegionsCliError::SatelliteCapBelowBundleRadius`), which is where a flag-pair rule belongs.

## Two things found on the way

**`ng_generic_walk_probe`'s own tests were failing at `2534e6e7`.** The probe reads its regions from
the catalog and refuses to run without one; its synthetic fixture wrote a FASTA and no catalog, so
seven of its nineteen tests aborted at start-up. The fixture builds a catalog beside itself now.

**`RegionScanner` has no caller anywhere**, and did not have one before this change either. It keeps
`scan_windowed`, `WindowCursor`, `ScannedWindow` and `scan_window` alive in
[`tandem_repeat.rs`](../../../../src/ng/tandem_repeat.rs) — roughly 300 lines of windowed-scan
machinery whose only consumer is a region seam nothing consumes. Deleting it is a separate decision
from this one and was left alone.

## Verification

Run in the dev container at `8bb7f07e`:

- `cargo fmt` clean;
- `cargo clippy --lib --tests --all-features -- -D warnings` clean;
- `cargo test --tests`: 3,292 library tests, 9 ignored, and every integration suite passing —
  including the differential's 11 and the new anchor's 4;
- the four `#[ignore]`d tomato identities, before and after, byte-identical (table above).
