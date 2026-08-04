# ng — read filtering in stages, step B2: `fill_raw_read` fills a whole raw aligned read

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `d5fb526`
**Plan:** [read_filtering_stages.md](../../ng/impl_plan/read_filtering_stages.md) step **B2**
**Arch:** [read_filtering_stages.md](../../ng/arch/read_filtering_stages.md) §2

---

## 1. Plan

`DecodedContainer::fill_raw_read` takes `&mut NoodlesRawAlignedRead` instead of
`&mut RecordBuf`, and sets **both** its fields — the record and the read group — so the CRAM arm
stops stamping the group on the line after the call.

**Added by the owner at Checkpoint A**, from A3's review: the name says it fills a raw aligned
read and it filled half of one, which is the record-versus-read confusion Milestone A removed
everywhere else. A signature change rather than a rename, so it could not travel with A3.

## 2. Assumptions

**None.** The container already held both halves — `PackedReadEntry.owner` is resolved once per
decode — so this moves a line rather than moving information.

## 3. Changes made

- **[container.rs](../../../../src/ng/read/input/aligned_reads_reader/container.rs)** —
  `fill_raw_read(&self, i: usize, raw_read: &mut NoodlesRawAlignedRead)`, setting
  `raw_read.read_group` from the entry and filling `raw_read.record` as before. Its doc now
  says why the read group comes from here: on CRAM it is a container-level number decided at
  decode, not a per-record `RG` tag, which is what makes this arm the documented exception to
  the readers' "records come out raw, read group cleared" contract. Stating it here puts the
  exception in one place instead of two.
- **[cram.rs](../../../../src/ng/read/input/aligned_reads_reader/cram.rs)** — the call site
  loses its second line.
- **`DecodedContainer::read_group(i)` deleted.** It existed for exactly one caller — the line
  just removed — and nothing outside the file needs to ask which group an entry has now.
  (`clippy` would not have caught it: `aligned_reads_reader/mod.rs` carries a module-level
  `#![allow(dead_code)]`, so this was checked by hand.)

## 4. Tests

**This section first said "none added; none needed adding, and that was verified rather than
assumed." The review showed that conclusion was drawn from the one mutation that could not
survive.** One test was added and one strengthened; the full log is below.

`raw_read.read_group = None` trips `RawAlignedRead::decode`'s `Option` guard, so it *cannot*
survive — it proves the stamp is *present*, not that it is *right*. The reviewer ran six more:

| mutation | result |
|---|---|
| `read_group = None` | killed — by both tests §4 originally credited |
| the stamp deleted entirely | killed — by both |
| `Some(ReadGroupId(0))` — a wrong but *valid* group | killed — **by one test only** |
| `Some(self.index[0].owner)` — the wrong entry's group | killed — **by one test only** |
| a stale group from the second container onwards | **SURVIVED all 1,542** |
| `out.data_mut().clear()` deleted | **SURVIVED all 1,542** |
| the deleted `read_group(i)` re-added, uncalled | clippy silent — confirms §3's parenthetical |

**Two things followed.**

**A stale group past a container boundary was invisible, and that is the defect B2 most risks.**
The gap was a hole between two fixtures: the only multi-read-group CRAM is three records, which
noodles writes as one container; the only multi-container CRAM declares one `@RG` and its cursor
tests open it as `ReadGroupResolution::Sole`, an arm that never asks a record which group it is
in. So **no test reached the per-record read-group arm past a container boundary at all**.
Added: `multi_container_cram_two_read_groups` (two libraries of one sample, alternating) and
`a_read_past_the_first_container_carries_its_own_read_group`. Verified: passes on this code,
fails under the stale-group mutation at the second container's first record.

**The value rested on a single test.** Both wrong-value mutants died only to
`a_cursor_keeps_every_read_group_of_its_sample_not_just_one`;
`a_shared_cram_serves_each_open_only_its_own_reads` — the second test §4 credited — collected
`(qname, read_group)` pairs and threw the group away with `.map(|(qname, _)| qname)`. It now
asserts the pairs, so both mutants die twice.

Suite **2,841 → 2,842** (`ng::` 1,542 → 1,543).

## 5. Validation results

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib` | **2,842 passed**, 0 failed, 5 ignored (2,841 + the review's one new test) |
| `cargo test --lib ng::` | **1,543 passed**, 0 failed, 2 ignored |
| `cargo test --examples` | 52 passed, 0 failed |

**The four acceptance dumps are byte-identical** to the `8cf6f03` baseline by `cmp`, and the
walk probe prints the anchor exactly at `seconds=1.845`.

**The two tomato dumps are the ones that matter here.** B2 touches only the CRAM path, and
`ng_generic_loci_dump` / `ng_ssr_loci_dump` on `SL4.0ch01` (1,718,914 and 11,945 lines) run
through it end to end — every read of that chromosome served by `fill_raw_read` and carrying a
read group it set. The two HG002 dumps are BAM and cannot see this change at all.

## 6. Tradeoffs and follow-ups

- **This is the last step of Milestone B.** Checkpoint B is next: the dumps and the anchor
  unchanged (they are), and the walk probe's `seconds` measured before and after (B1's report
  has the A/B; B2 does not move it, being one field assignment relocated).
- **`ReadFilterBuffers` and `with_validated_contigs`** remain, one caller each, for C2.
- **`container.rs` still has no test module**, and the review named three input classes on
  `fill_raw_read` that nothing reaches: a record with **no name** (every fixture record is
  named), an **empty** sequence/quality/CIGAR span, and the clear-and-refill claim its own doc
  makes — no test serves a long read then a short one through one buffer and checks the short
  one carries no tail. A fourth, `out.data_mut().clear()`, has a doc-comment invariant naming
  its own silent failure mode and **deleting it leaves the suite green** (currently unreachable,
  since every caller passes a fresh buffer — latent, not live). Deferred: these are pre-existing
  coverage of pre-existing behaviour, and a new test module is its own piece of work. **Worth
  taking before C2**, which moves the filtering loop and is exactly the change that could hand
  this function a buffer with a history.
