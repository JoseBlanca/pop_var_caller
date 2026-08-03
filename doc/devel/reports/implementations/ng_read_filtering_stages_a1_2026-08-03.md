# ng — read filtering in stages, step A1: the raw read says what it is

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `8cf6f03`
**Plan:** [read_filtering_stages.md](../../ng/impl_plan/read_filtering_stages.md) step **A1**
**Spec:** [read_filtering_stages.md](../../ng/spec/read_filtering_stages.md) §6 ·
**Arch:** [read_filtering_stages.md](../../ng/arch/read_filtering_stages.md) §2, §3.1

---

## 1. Plan

`RawRecord` → `RawAlignedRead` and `NoodlesRawRecord` → `NoodlesRawAlignedRead`, both **moved
out of `read/filtering.rs` into `read/aligned_read.rs`**, beside `AlignedRead` and the
conversion that already lived there. The trait's doc gains the fact that an unmapped read is
one of these.

**No behaviour change at all.** That is the milestone's whole point: renames land alone so
that when B and C change behaviour, the diff *is* the behaviour. The bar is therefore unusual
— the suite count must not move and the four acceptance dumps must be byte-identical.

## 2. Assumptions

Two choices the plan left to the implementer, both recorded rather than silent:

1. **The three `NoodlesRawAlignedRead` tests moved with the type.** `arch` §8 says tests stay
   beside the code and that `filtering.rs`'s tests split by subject; these three are about the
   adapter's flag/MAPQ reads and its refusal to decode an unstamped buffer, not about any
   keep-or-drop rule. They are renamed `noodles_raw_record_*` → `noodles_raw_aligned_read_*`
   and otherwise unchanged, assertion for assertion. **The suite count is unchanged by the
   move** (they moved, none was added or dropped). It is changed by the *review*, separately
   and deliberately — see the fix-application report.
2. **`bam_record` was deleted, not left behind.** It had exactly one call site in
   `filtering.rs` — `noodles_raw_record_reads_flag_mapq_and_decodes`, the first of the three
   moved tests; the other two built their records inline. It went with that test in the shape
   the test actually needs (`record_with_mapq_and_flags` in `aligned_read.rs`), which builds a
   byte-for-byte equivalent record — verified field by field by the review's `refactor_safety`
   agent. Deleting it left five `use` lines (seven noodles items) unused in `filtering.rs`'s
   test module; those went too, which is what `clippy -D warnings` reported. No assertion
   changed.

## 3. Changes made

- **[src/ng/read/aligned_read.rs](../../../../src/ng/read/aligned_read.rs)** — gains
  `RawAlignedRead` (the trait: its name, its opening sentence — now arch §3.1's "one alignment
  record as it comes off the file, undecoded" rather than "a borrowed view of one alignment
  record" — and the new unmapped-read paragraph; every method doc verbatim) and
  `NoodlesRawAlignedRead` with its hand-written `Default` and its `RawAlignedRead` impl, moved
  verbatim apart from the intra-doc link targets that had to follow the move. The module doc
  now opens on **one read in two states** and says plainly that no keep-or-drop rule lives
  here. `decode_record`'s parameter type is spelled `RecordBuf` rather than
  `sam::alignment::RecordBuf` — the same type, one import instead of two.
- **[src/ng/read/filtering.rs](../../../../src/ng/read/filtering.rs)** — loses both types and
  their impls, keeps everything else. Its module doc no longer claims to hold the noodles
  adapter and points at `aligned_read` for the read itself. A `//` note where they stood says
  where they went and why, in the same style as the two notes already there.
- **[src/ng/read/mod.rs](../../../../src/ng/read/mod.rs)** — the re-export moves from the
  `filtering::` group to the `aligned_read::` one, so the public path follows the type.
- **Six call-site files** — `read/input/mod.rs`, `read/input/region_records.rs`, and the four
  `read/input/record_reader/{mod,bam,cram,in_memory}.rs` — take the new names and, where they
  imported the type, the new path.

## 4. Tests added/updated

**None added, none removed.** Three moved and were renamed to match their subject:

| test | what it validates |
|---|---|
| `noodles_raw_aligned_read_reads_flag_mapq_and_decodes` | the two cheap pre-decode reads come off the packed record, and `decode` produces the `AlignedRead` phase two consumes, carrying the stamped read group |
| `noodles_raw_aligned_read_maps_unavailable_mapq_to_zero` | SAM `0xFF` reads as `MapQual(0)`, so any non-zero minimum drops it |
| `noodles_raw_aligned_read_decode_errors_on_a_record_with_no_position` | *as moved:* only that `decode` returns `InvalidData` on a defaulted buffer. **The review proved it never reaches the position path its name claims** — the defaulted buffer's missing read-group stamp fires first. Repurposed during fix application; see the fix report. |

## 5. Validation results

Run on the host, in debug except where noted.

| command | result |
|---|---|
| `cargo fmt` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | `Finished dev profile … in 8.56s`, no warnings |
| `cargo test --lib` | **2,837 passed**, 0 failed, 5 ignored — unchanged from the base |
| `cargo test --lib ng::` | **1,538 passed**, 0 failed, 2 ignored — unchanged from the base |
| `cargo test --examples` | 52 passed, 0 failed |

**The four acceptance dumps, byte-identical to the `8cf6f03` baseline** (`cmp`, not a line
count):

| dump | lines |
|---|---|
| `ng_generic_loci_dump` GRCh38 / HG002 `chr21` | 251,792 |
| `ng_ssr_loci_dump` GRCh38 / HG002 `chr21` | 4,406 |
| `ng_generic_loci_dump` tomato / SRR5079860 `SL4.0ch01` | 1,718,914 |
| `ng_ssr_loci_dump` tomato / SRR5079860 `SL4.0ch01` | 11,945 |

**The walk probe on `chr21`** prints the anchor exactly —
`loci=236081 observations=251786 reads_admitted=54709` — at `seconds=1.880`,
`loci_per_second=125569` (baseline in the same session: `seconds=1.846`, `126k`; the spec's
recorded figure is `1.876`). A rename is not expected to move this and it did not.

**Not the default gate, for reasons already on record:** `cargo test --release` is red on a
clean tree (four tests assert on `debug_assert!` messages release compiles out),
`cargo test --all-targets` aborts on a pre-existing panic in `benches/psp_writer_perf.rs:386`,
and `cargo doc` has 12 pre-existing unresolved links.

## 6. Tradeoffs and follow-ups

- **`filtering.rs` still holds `RecordSource`, `resolve_read_group` and the loop.** A1 moves
  the read, not the seam; `RecordSource` goes at C3 and the loop at C2.
- **`NoodlesRawAlignedRead`'s doc still links `RecordSource::read_next`** — a forward
  reference into `filtering.rs` that C3 deletes. Left pointing at the live item rather than
  softened to prose, so the compiler's intra-doc-link check keeps it honest until then.
- **`read/filtering.rs`'s module doc still describes step 1 as one thing.** Rewriting it is
  C2's, where it stops being true.
