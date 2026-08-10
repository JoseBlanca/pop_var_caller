# ng repeat catalog — B1 + B2: the file, and reading it back

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **B1** and **B2**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §3.4, §3.5, §4.3, §6 and
[`arch/repeat_catalog.md`](../../ng/arch/repeat_catalog.md) §2.3, §3.*

## Plan

**B1** writes the file: the schema, one row group per contig, the header in the footer's metadata,
and an atomic rename. **B2** opens one, checks it against the reference this run computed, and
streams rows back. Two commits.

## Assumptions and deviations

1. **The header is tab-separated text, not JSON.** The arch doc said "JSON-encoded under one key".
   Two of the values it carries — the period range and the copy-floor table — have **checked
   constructors**, and any derived deserialisation would rebuild them field by field and skip the
   check. The text encoding decodes through `PeriodRange::new` and `MinCopies::new`, so a header
   claiming a period range of 6..=1 is rejected rather than believed. It also avoids putting
   `Serialize` derives on validated production-adjacent types.
2. **`contig` is stored as the header's index (`UInt32`), not a dictionary of names.** The arch
   doc's decision record already says rows address contigs by index; the names are in the header, and
   inside a row group — one contig — the column is one repeated value that RLE removes. Simpler than
   a string dictionary and closer to the type the reader hands back.
3. **The contig table stores the `.fai` geometry too** (offset, line bases, line width), not just
   name/length/MD5. It costs three integers per contig and lets a decoded `ContigInfo` be the real
   thing rather than one with zeroed fields that no longer describes the FASTA.
4. **`repeats_in_region` takes an `Option<ContigId>` for now**, not an interval set. Interval
   restriction is deferred to the consumer that needs it (spec §8); contig restriction is what the
   builder and the differential need, and it is what one row group per contig buys.
5. **Two knock-on changes outside the module.** `arrow-array` and `arrow-schema` are promoted from
   transitive to explicit, because ng names their types. And one assertion in
   [`open_record.rs`](../../../../src/ng/locus_generation/pileup/open_record.rs) needed its element
   type spelled out: parquet pulls in `serde_json`, whose `impl PartialEq<Value> for u64` makes a
   bare `Vec::new()` ambiguous. Both are recorded in the commits.

## Changes made

- [`src/ng/repeat_catalog/parquet_file.rs`](../../../../src/ng/repeat_catalog/parquet_file.rs) — the
  nine-column schema, `RepeatCatalogWriter` (row group per contig, atomic `.tmp` + rename), the
  header's encode/decode, and `row_from_batch`. **Three writer settings are fixed by us** — codec,
  level, `created_by` — because the default `created_by` carries the parquet crate's version and
  spec §6 asks for byte-identical output.
- [`src/ng/repeat_catalog/reader.rs`](../../../../src/ng/repeat_catalog/reader.rs) — `RepeatCatalog`,
  `open_checking_against_reference`, `header`, `contigs`, `repeats_in_region`.
- [`src/ng/repeat_catalog/mod.rs`](../../../../src/ng/repeat_catalog/mod.rs) — `RowsByPeriod`, the
  build tally the CLI prints.

**Validation is where the FASTA is the source of truth**: the header's contig table is compared
against the digests *this run* computed — names, order, lengths, per-contig MD5s, then the
whole-reference digest. A `.fai`-only read carries no digests, so such a run gets the order and
length checks and nothing pretends otherwise.

## Tests added

16 new tests (39 in the module; 2,914 in the lib, all green).

**B1:** the header round-trips through its text encoding; an unknown line is a malformed file; a
future header format is refused rather than guessed at; rows survive a write and a read; **each
contig is its own row group** (asserted on the file's metadata, not inferred); the header comes back
from the footer; **two writes of the same rows are byte-identical**; and a writer that never
finishes leaves no readable catalog — nothing under the final name, and the leftover `.tmp` cannot be
opened at all.

**B2:** a catalog of this reference opens and its rows come back in file order; one contig can be
read on its own; **a missing catalog names the command that builds one** while **a catalog of another
reference names the contig that differs** — the two failures a caller reacts to differently, kept
apart; reordered contigs of equal length are caught (only the digests and the order can); a matching
contig table with a different whole-reference digest is refused; a truncated catalog will not open;
and a parquet file that is not a catalog is refused rather than read as empty.

## Validation

In the dev container: `cargo fmt` clean; `cargo clippy --lib --tests --all-features -- -D warnings`
clean; `cargo test --lib` **2,914 passed, 0 failed, 5 ignored**.

**Pre-existing and untouched:** `cargo clippy --all-targets` still fails on
`examples/ng_inbreeding_harness.rs` (`076cb5e9`).

## Tradeoffs and follow-ups

- `row_group_contig` reads a row group's contig from the column's **statistics**, which parquet
  writes per row group. If a future writer disabled statistics this would silently fall back to
  scanning every group; the writer enables them by default and the "one contig can be read on its
  own" test would still pass, so this is a performance property without a guard. Worth one when C2's
  builder lands and the file has more than three row groups.
- Interval-level restriction, and the criteria-taking read methods, are milestone D.
