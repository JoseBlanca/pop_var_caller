# ng — read filtering in stages, step B1: the contig check becomes a comparison

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `bfb54dd`
**Plan:** [read_filtering_stages.md](../../ng/impl_plan/read_filtering_stages.md) step **B1**
**Spec:** [read_filtering_stages.md](../../ng/spec/read_filtering_stages.md) §9 Q2 ·
**Arch:** [read_filtering_stages.md](../../ng/arch/read_filtering_stages.md) §6

*Revised after review. The shape the review changed is recorded in §2 and §6; the
[fix report](../reviews/fixes_applied_2026-08-03_v4.md) has the per-finding accounting.*

---

## 1. Plan

`AlignmentFile::cursor` compares the accessor's contig table against the file's with
`ContigList::first_disagreement`, and the per-contig fetch loop stops running. `+ ContigTable`
joins `cursor`'s `R` and propagates to `SampleReads::cursor` and the pileup generator. Its own
commit, with a mutation-verified test that a mismatched accessor is refused.

## 2. What the review changed, before anything else

**The first implementation replaced the loop with one comparison. That was not enough, and the
reviewers proved it two ways.**

**A table comparison does not imply the accessor can serve the bases.** `InMemoryRefSeq` derives
its table from the bytes it holds, so for it equality *does* imply resolvability. But
`ResidentRefSeq::new` and `WindowedRefSeq::new` take the `ContigList` as a **constructor
argument independent of the bytes** — so a `WindowedRefSeq` over a FASTA missing a contig,
behind a table that names it, matched the file perfectly and could not serve a single base.
Measured: the deleted probe rejected that; the first B1 accepted it, and the fault surfaced
mid-stream under a message naming the *BAM*. So the fail-fast the loop provided was **gone, not
moved**, and my "proves strictly more" claim was wrong.

**The fix keeps the cheapness and restores the guarantee: three checks, in order.**

1. the contig is one this file declares (`CursorContigNotInFile`),
2. the accessor's table equals the file's, **order included**
   (`CursorAccessorContigTable`),
3. the accessor can serve **this cursor's own contig** — one zero-length fetch
   (`Reference`).

(3) is what the loop was really for. The loop asked it of *every* contig in the header — ~2,580
on GRCh38, once per cursor — for a property that only matters for the one contig the cursor will
read. Asking it once keeps the fail-fast at **one** `open(2)` instead of 2,580, and revives
`AlignmentFileError::Reference`, which the first implementation had left with no producer.

**The order was also wrong** and is now argument → description → ability; the probe needs the
contig index checked first.

## 3. Assumptions and deviations

**One deviation from the spec, taken on the review's argument and flagged for the owner.**

Spec §9 Q2's snippet reuses `AlignmentFileError::ContigReconcile`. The first implementation did
that, prefixing the detail string to say which check fired. The `errors` reviewer showed the
resulting message is **false at that point in the run**:

> alignment file '…/sample.bam' does not match the reference contig table: the accessor passed
> to cursor() is over a different table: name disagreement at index 0 ('chr1' vs 'not_chr1')

`open` has already proved the file *does* match the reference; what failed is the accessor the
caller wired up. Two clauses that both say "table" contradict each other before the true one
arrives, and the two failures become discriminable only by substring — which the first version
of the new test did, teaching the pattern.

**So B1 adds `AlignmentFileError::CursorAccessorContigTable`.**

**How far this deviates was itself checked, and the first answer here was wrong.** This report
originally said it went "against spec §1's *adds no new error*". B2's review found no such
sentence: spec §1 says "change the meaning of any error", which adding a variant does not do,
and the only "No new error type" statement is **arch §4**, scoped to `ReadFilterError` — a
different enum. What actually exists is §9 Q2's *illustrative* snippet, and arch's own preamble
says "Signatures are illustrative; the **contract** is the deliverable."

So the design authority is **silent** on adding a variant to `AlignmentFileError`, not against
it. Still recorded for the owner at Checkpoint B — but as a gap filled, not a rule broken.

Unchanged from the first implementation: the file is the left operand of (2), so the message
reads *file value vs accessor value*, the direction the open gate prints.

## 4. Changes made

- **[open_bam.rs](../../../../src/ng/read/input/open_bam.rs)** — `cursor` takes
  `R: RawRefSeq + ContigTable` and runs the three checks.
- **[cursor.rs](../../../../src/ng/read/input/cursor.rs)** — `over_records` builds its filter
  through `with_validated_contigs` and is now **infallible**.
- **[mod.rs](../../../../src/ng/read/input/mod.rs)** — the new error variant;
  `AlignmentFileError::Reference`'s doc rewritten (it is produced again, and now scoped to one
  contig); `SampleReads::cursor`'s doc gains the new failure mode.
- **The bound propagated** to `SampleReads::cursor` and the pileup generator's four sites. The
  SSR generator already required `ContigTable`.
- **[filtering.rs](../../../../src/ng/read/filtering.rs)** — `ReadFilter::new` **deleted** (its
  body *is* the loop, and it had no caller left); `RecordSource::header` deleted with it — its
  only caller was the loop. `ReadFilterError::Reference`'s doc rewritten: it explained itself
  through the deleted constructor and read backwards.

**Three consequences worth naming, each a subtraction the change earned:**

- **`filtering.rs` no longer knows what a SAM header is.** `use noodles_sam as sam` left the
  module; what remains needs `RecordBuf` alone.
- **Its tests no longer build one either** — `contig_header` and `one_contig_header` went with
  three imports. Every test in the module had been constructing a header so its `RecordSource`
  could answer `header()` for the probe.
- **The two test doubles lost their `header` field**, and `FakeSource::new` its second
  parameter.

## 5. Tests

**2,839 → 2,841** (`ng::` 1,540 → 1,542). Three deleted, five added.

| test | disposition |
|---|---|
| `read_filter_new_rejects_a_contig_missing_from_the_reference` | **deleted; both its halves have successors** — see below |
| `probe_free_constructor_filters_identically_to_new` | deleted, **no successor possible** — it compared two constructors, and there is one |
| `probe_free_constructor_skips_the_contig_probe_new_would_fail` | deleted, likewise |
| `a_cursor_refuses_an_accessor_over_a_different_contig_table` | added — wrong names, and right names at wrong lengths |
| `a_cursor_refuses_an_accessor_whose_contig_table_is_a_permutation` | added — **the Blocker's fix** |
| `a_cursor_refuses_an_accessor_whose_contig_table_is_shorter_than_the_files` | added — the count-mismatch branch, and the deleted test's actual input class |
| `a_cursor_refuses_an_accessor_that_cannot_serve_its_contig` | added — the deleted test's *resolvability* half |
| `the_fixture_accessors_carry_the_same_contig_table_as_the_fixture_files` | added — a standing guard on the three test call sites that bypass `cursor` |

**The first accounting of the deleted test was half wrong**, and the reviewer caught it: it
asserted `Err(RefSeqError::UnknownContig(ContigId(1)))` — a *resolvability* failure at a *count*
mismatch. Only the table half had moved. Both halves have successors now.

### Mutations, each `grep`-confirmed present before its run and absent after the revert

| mutation | result |
|---|---|
| delete the comparison | **killed** |
| compare **names only** (what the probe effectively did) | **killed**, on the wrong-lengths case |
| **sort both tables before comparing** — same messages, order-sensitivity removed | **killed**, by the permutation test *alone* |
| delete the single-contig probe | **killed**, by the resolvability test *alone* |

**The third is the Blocker the review found.** Before the permutation test existed, an
order-insensitive rewrite — the natural move on a 2,580-entry ordered walk, on a perf branch —
passed **all 1,538 tests**, and a permuted table is precisely what makes filter #8 score every
read against another chromosome's bases with no error at all.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib` | **2,841 passed**, 0 failed, 5 ignored |
| `cargo test --lib ng::` | **1,542 passed**, 0 failed, 2 ignored |
| `cargo test --examples` | 52 passed, 0 failed |
| `cargo doc --no-deps --lib` | 12 unresolved links, all pre-existing, none in a touched file |

**The four acceptance dumps are byte-identical** to the `8cf6f03` baseline by `cmp`, and the
walk probe prints the anchor exactly.

### ⚠ The check found that every cursor test had been passing an accessor over a different table

Applying it turned **23 tests red at once** — a defect in the fixtures, not the check. Every
cursor fixture used `InMemoryRefSeq::from_contigs`, which names contigs `contig0`, `contig1`, …
while the fixture files declare `chr1`, `chr2`. All 23 had been handing `cursor()` an accessor
disagreeing on **every contig name**, and nothing noticed, because the probe of the day fetched
a window per contig and a window resolves whatever the contig is called.

They now go through one shared `fixture_reference_bases()` built from `FIXTURE_CONTIGS`. The
reviewer verified this is a fix and not a weakening three ways: `from_contigs` delegates to
`from_named_contigs` with the *same byte vectors*, so the bases are byte-identical and only the
labels changed; no surviving test had an assertion edited; and reverting the helpers reproduces
exactly 23 failures.

### The measurement Checkpoint B asks for

Clean A/B, same machine, same session, six runs each:

| | runs (`seconds`) | mean |
|---|---|---|
| **before** (`bfb54dd`) | 1.868, 1.858, 1.854, 1.862, 1.862, 1.861 | **1.861** |
| **after** (B1 as shipped) | 1.845, 1.831, 1.833, 1.834, 1.826, 1.836 | **1.834** |

**≈27 ms per cursor, ~1.4 %.** An intermediate measurement of the pre-review implementation
(no probe at all) gave 1.844; the difference between that and the shipped 1.834 is smaller than
the between-session spread, which is what one extra `open(2)` per cursor should look like. The
honest figure is **~20–27 ms**, and it is consistent rather than noise: the slowest *after* run
beats the fastest *before*.

**It is not the ~130 ms the estimate predicted, and the estimate's arithmetic is wrong in a
checkable way.** Spec §9 Q2 takes "roughly 52 µs per open with a shared index" from
[ref_seq.rs:622-655](../../../../src/ng/ref_seq.rs#L622). Read in place, that 52 µs is the cost
of **constructing a `WindowedRefSeq`**, and the same comment breaks it down: 34 µs of it is
cloning the 2,580-entry contig table. The loop constructed no accessors — it called `fetch_into`
on the one it was given — so the table clone was never in its per-contig cost. Multiplying 52 µs
by 2,580 multiplies in a cost the loop did not pay.

*(The same comment attributes the ~18 µs residual to an `open(2)`. The reviewer notes that
cannot be right for `with_shared_index`, which opens nothing at construction — so ~18 µs × 2,580
is not a defensible replacement estimate either. What is defensible is the measurement.)*

**At run scale:** one cursor per file per chromosome, so 50 samples × 25 chromosomes ≈ 1,250
cursors ≈ **25–34 s per run**. Real, worth having, an order of magnitude below the estimate.

Speed was never why this change was made (spec §1: "a constraint here, not a goal"), and with
(3) restored the three checks together do now prove strictly more than the loop did — names,
lengths *and* order, plus the same fetchability for the contig that matters.

## 7. Follow-ups

- **Spec §9 Q2's cost arithmetic should be corrected**, and §9 Q2's snippet updated to the
  variant that shipped. Owner's call — a design document.
- **`RegionRawAlignedReads` no longer has a `header()`**; arch §3.3 lists one among the inherent
  methods C3 keeps. Re-add only if a caller appears; none exists.
- **`ReadFilterBuffers` and `with_validated_contigs`** are alive with one caller each. Spec §10
  assigns their removal to C2.
- **`ResidentRefSeq::new` and `WindowedRefSeq::new` can build a lying accessor** — a
  `ContigList` unrelated to the bytes they will serve. Check (3) contains the damage at the
  cursor; the constructors themselves are an API-design question for another day.
