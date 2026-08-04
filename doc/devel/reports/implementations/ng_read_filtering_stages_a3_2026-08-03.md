# ng — read filtering in stages, step A3: the region narrowing, and two names below it

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `db3057a` (A2)
**Plan:** [read_filtering_stages.md](../../ng/impl_plan/read_filtering_stages.md) step **A3**
**Spec:** [read_filtering_stages.md](../../ng/spec/read_filtering_stages.md) §6 ·
**Arch:** [read_filtering_stages.md](../../ng/arch/read_filtering_stages.md) §2

---

## 1. Plan

`RegionRecords` → `RegionRawAlignedReads` (**the file follows the type**, spec §6),
`DecodedContainer::fill_record` → `fill_raw_read`, and the private `RecordIndex` →
`RawReadIndex`. The last of Milestone A, and the same bar as A1 and A2: suite count unmoved,
four dumps byte-identical.

## 2. Assumptions

**One thing had to be got right, and it is the reason this step is not a blind substitution.**

**`region_records` is also a production identifier, and production is frozen.**
`PspReader::region_records` is a public method of the `.psp` reader — 19 occurrences across
`src/psp/reader.rs`, `src/psp/mod.rs`, `src/pop_var_caller/psp_to_pileup.rs`,
`src/var_calling/sample_reader.rs` and `src/regions.rs`, plus two more in
`benches/psp_reader_perf.rs` — and it has nothing to do with alignment records. A repo-wide
substitution would have renamed a production API in a milestone whose whole claim is that it
changes nothing. **The substitution was therefore scoped to a named list of ng files**, and
`src/psp/reader.rs` was checked afterwards to still hold its 14 occurrences. `RecordIndex` and
`fill_record` needed no such care — both are ng-only, and both live in
`aligned_reads_reader/container.rs`.

The other choices are mechanical: `git mv` for the file, so `git log --follow` survives.

## 3. Changes made

- **`src/ng/read/input/region_records.rs` → `region_raw_aligned_reads.rs`** (`git mv`), and the
  module declaration with it.
- **The four identifiers**, across **eight** ng files: `RegionRecords` → `RegionRawAlignedReads`
  (21 sites), the module path (4 sites in ng — production's 19 are occurrences of a *method
  name*, not a module path, and were left alone), `fill_record` → `fill_raw_read` (2),
  `RecordIndex` → `RawReadIndex` (5).
- **The layer diagram in the renamed file re-padded**, and the module's title line now says
  "this region's **raw aligned reads** only", which is what the type is now called.
  `RegionRawAlignedReads` is nine characters longer than `RegionRecords`, so the description
  column had drifted out of line — the same breakage A2's review caught in this very diagram.
  A first attempt bought the width by cutting the gloss to "this region's only"; the review
  caught the dangling possessive, and the column was widened instead so the full phrase fits,
  matching spec §6's own diagram.
- **Seven comment blocks — sixteen lines — re-wrapped** where the longer names pushed them past
  100 columns. `cargo fmt` does not reflow comments, so `--check` stays clean either way and
  nothing would ever have flagged them. Two further formatting edits, neither caused by a name:
  rustfmt reflowed a `use` in `cursor.rs` and added a trailing comma, and the `held` field's doc
  was re-wrapped although it contains no renamed identifier.

## 4. Tests added/updated

**None.** A3 is a substitution; no test's subject changed, and no test expectation contains any
renamed identifier as a string.

## 5. Validation results

Host, debug except the dumps.

| command | result |
|---|---|
| `cargo fmt` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean, 5.57s |
| `cargo test --lib` | **2,839 passed**, 0 failed, 5 ignored — unchanged from A2 |
| `cargo test --lib ng::` | **1,540 passed**, 0 failed, 2 ignored — unchanged from A2 |
| `cargo test --examples` | 52 passed, 0 failed |

**The four acceptance dumps, byte-identical to the `8cf6f03` baseline** by `cmp`: 251,792 /
4,406 / 1,718,914 / 11,945 lines. **The walk probe** prints the anchor exactly —
`loci=236081 observations=251786 reads_admitted=54709` — at `seconds=1.922`.

**Both greps, including the widened one A2's review made this step inherit:**

```
grep -rn  "RegionRecords\|region_records\|fill_record\|RecordIndex" src/ng      → no matches
grep -rniE "region[ -]records|record[ -]index|fill[ -]record|record[ -]readers?" src/ng
                                                                                → one match
```

The one match is `locus_generation/pileup/generator.rs:605`, *"The region records are clamped
to"* — a different sense of the words (the region **that** records are clamped to), not the type
name. Left alone deliberately.

## 6. Tradeoffs and follow-ups

- **The `fill_raw_read` wrapper was collapsed after all, on the review's argument.** It was a
  one-line `pub(crate)` shim over a private `fill`, each with exactly one caller, and the rename
  made the two ends disagree — a reader following `fill_raw_read` → `fill` is told the operation
  is about a raw read and then that it is about nothing in particular. This report first
  deferred it; the reviewer pointed out that **collapsing a 1:1 name-only indirection is a
  rename plus a four-line deletion**, so deferring would push a naming edit into a later
  behaviour diff. Accepted. The surviving method carries `fill`'s doc plus a note on why `out`
  is only the record half of the caller's buffer.
- **Two names the architecture prescribes are disputed by the review**, and were landed as
  specified rather than changed: `RawReadIndex` (arch §2) — "Index" names the container's field,
  not one entry, and `RawRead` abbreviates A1's `RawAlignedRead`, which appears in full 85 times
  in `src/ng`; and `fill_raw_read` (arch §2) — it takes `&mut RecordBuf`, so its caller sets the
  read-group half on the next line. Both are **Checkpoint A questions for the owner**, since
  changing either means amending arch §2's table.
- **Milestone A is complete with this step.** The other things deferred from A1 and A2 — the two
  visibility questions and the one design-doc sweep covering all three renames — are the
  checkpoint's, and are recorded in `PROJECT_STATUS.md`.
