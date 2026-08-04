# ng — read filtering in stages, step A2: the readers say what they read

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `5438927` (A1)
**Plan:** [read_filtering_stages.md](../../ng/impl_plan/read_filtering_stages.md) step **A2**
**Spec:** [read_filtering_stages.md](../../ng/spec/read_filtering_stages.md) §6 ·
**Arch:** [read_filtering_stages.md](../../ng/arch/read_filtering_stages.md) §1, §2

---

## 1. Plan

`RecordReader` → `AlignedReadsReader`, its three arms to `BamAlignedReadsReader` /
`CramAlignedReadsReader` / `InMemoryAlignedReadsReader`, and the module `record_reader/` →
`aligned_reads_reader/`. **The readers deliberately leave "raw" out of their names** (spec §6,
owner's call), so each reader's doc now has to state that what it yields is undecoded — the type
name no longer carries it.

No behaviour change, same bar as A1: suite count unmoved, four dumps byte-identical.

## 2. Assumptions

**None that change direction.** Two mechanical choices worth recording:

1. **The directory moved with `git mv`**, so `git log --follow` stays cheap on all five files.
2. **The substitution ran longest-name-first** — `BamRecordReader`, `CramRecordReader`,
   `InMemoryRecordReader`, then the bare `RecordReader`, then `record_reader`. `RecordReader` is
   a substring of the three arm names, so the other order would have produced
   `BamAlignedReadsReader` from a half-renamed identifier. Verified before and after: the
   family has exactly five members (`grep -o` over `src/`, 11 / 9 / 17 / 30 / 23 occurrences),
   and afterwards the same five counts appear under the new names with **zero** old-name hits.

## 3. Changes made

- **`src/ng/read/input/record_reader/` → `src/ng/read/input/aligned_reads_reader/`** (`git mv`),
  five files: `mod.rs`, `bam.rs`, `cram.rs`, `in_memory.rs`, `container.rs` (the last carried
  no reference and changed not at all).
- **The four type names at all 67 sites, and the module path at 23 more** — 90 substitutions
  across eleven files (the type names span nine of them; `read/input/mod.rs` and
  `test_fixtures.rs` carry only the module path, and `container.rs` carries neither).
  Two of the 90 are inside string literals, and they are the only non-comment bytes in the diff
  that are not identifiers: the `debug_struct("…")` labels in `bam.rs` and `cram.rs`. Nothing
  `{:?}`-formats these readers, so they are unobservable.
- **The undecoded statement, in five places** — the module doc, the enum's own doc, and each of
  the three arms. Each says what the type yields and that the conversion happens above, only for
  the reads that clear the flag/MAPQ filter.
  - The **CRAM** arm's says more, because "undecoded" is ambiguous exactly there: a CRAM
    container *is* decompressed and decoded into `RecordBuf`s before any record can be looked
    at. The doc names the two senses so a reader does not conclude the arm violates the
    contract.
  - The **in-memory** arm's says more for the opposite reason: its records are handed in already
    built, so it would be easy to assume they arrive converted.
- Two doc sentences reworded from "records" to "reads" where the new type name made the old
  wording read oddly (`An [AlignedReadsReader] finds reads and unpacks them`).

## 4. Tests added/updated

**None.** A2 is a substitution; no test's subject changed, and no test expectation contains any
of the renamed identifiers as a string.

## 5. Validation results

Host, debug except the dumps.

| command | result |
|---|---|
| `cargo fmt` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean, 8.11s |
| `cargo test --lib` | **2,839 passed**, 0 failed, 5 ignored — unchanged from A1 |
| `cargo test --lib ng::` | **1,540 passed**, 0 failed, 2 ignored — unchanged from A1 |
| `cargo test --examples` | 52 passed, 0 failed |

**The four acceptance dumps, byte-identical to the `8cf6f03` baseline** by `cmp`: 251,792 /
4,406 / 1,718,914 / 11,945 lines. **The walk probe** prints the anchor exactly —
`loci=236081 observations=251786 reads_admitted=54709` — at `seconds=1.871`.

`grep -rn "RecordReader\|record_reader" src` → **no matches**.

## 6. Tradeoffs and follow-ups

- **`RegionRecords` still has the old vocabulary**; it is A3's, along with
  `DecodedContainer::fill_record` and `RecordIndex`.
- **The contract list in `aligned_reads_reader/mod.rs` still says "records" throughout**, and
  after the review's Mi3 the module doc now *says why*: a record is the encoding of a read, the
  contract is about finding and unpacking those, and that is the only level at which the arms
  differ. The first version of this note defended the choice on grounds ("what a file holds")
  that the enum's own doc simultaneously claimed for the other word — the reviewers caught the
  contradiction, and the fix was to name the distinction rather than pick a winner.
- **The verification grep was too narrow, and the review found what fell through it.**
  `grep -rn "RecordReader\|record_reader" src` matches CamelCase and snake_case only; ten sites
  spelled the type "record reader", including a runtime `unreachable!` message. All ten are
  fixed, and the widened check is `grep -rniE "record[ -]readers?" src/ng`. **A3 must use it
  too** — the same blind spot would certify it the same incomplete way.
