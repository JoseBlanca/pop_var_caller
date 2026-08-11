# ng repeat catalog — milestone C: the catalog builds from a reference

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **C1**, **C2**, **C3**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §2.1–§2.4, §6 and
[`arch/repeat_catalog.md`](../../ng/arch/repeat_catalog.md) §2.1, §2.2.*

## Plan

**C1** gives the reference pass a seam; **C2** hangs the builder off it; **C3** lets it scan several
contigs at once without that reaching the file. Three commits.

## Assumptions and deviations

1. **The builder holds its first error and reports it at `finish`.** The seam is infallible so that
   `reference_info` stays a leaf that "knows nothing else about ng" — the arch doc's decision — so a
   failed write parks the builder, the rest of the pass is ignored, and `finish` returns the error.
2. **`--threads` is a batch, not a pipeline.** The builder holds finished contigs until it has
   `threads` of them, scans that batch with rayon, then writes the results **in reference order**.
   The arch doc said "hand a completed contig to a worker and let the stream continue"; a batch gives
   the same parallelism with no channel, no completion-order bookkeeping, and a memory bound that is
   obvious by construction (`threads` contigs, never more).
3. **Rows are sorted inside the scan**, not at write time. `find_tandem_repeats` emits period by
   period, so without a sort the file's order would be a property of the detector's loop.
4. **`finish` refuses a reference read without a digest.** A catalog whose header carries no
   reference MD5 could never be validated against anything, so writing one is worse than failing.

## Changes made

- [`src/ng/reference_info.rs`](../../../../src/ng/reference_info.rs) — the `ReferenceBasesObserver`
  trait, `IgnoreBases`, and `read_reference_info_observing`. `read_reference_info` keeps its
  signature and passes `IgnoreBases`, so no existing caller changes. The hooks sit where the facts
  are: `contig_started` where the name is complete, `bases` inside `flush_md5` (**the one place
  uppercased bytes leave the pass**, so the observer sees exactly what the digests see), and
  `contig_finished` where the `ContigInfo` is pushed.
- [`src/ng/repeat_catalog/builder.rs`](../../../../src/ng/repeat_catalog/builder.rs) —
  `RepeatCatalogBuilder`, `BuildTally` (rows kept **and** rejections charged), `with_threads`, and
  the free-standing `scan_contig` that runs on a worker.

## Tests added

13 new tests (52 in the catalog module and 58 in `reference_info`; 2,924 in the lib, all green).

**C1:** the observer sees `start → bases* → finish` per contig in file order; the bytes it receives
are the uppercased bases with terminators already gone; **the observed bases digest back to the
contig MD5 the pass reports** — the property that lets the builder scan what the digest attests to;
and a `.fai`-only read tells the observer nothing at all, because "nothing" and "an empty contig" are
different claims.

**C2:** a planted `(CAG)8` becomes one row at hand-checked 1-based coordinates with the right motif
and stratum; **a 2 kb tract comes out as one row** (spec §2.3's whole point, and the test a windowed
scan would fail); a tract at a contig's very start is not in the file **and is counted**; rows come
out contig-by-contig in reference order with one row group each; rows within a contig are ordered by
start, then period, then end; a digest-less reference cannot be catalogued.

**C3:** `the_thread_count_does_not_change_the_bytes` builds the same three-contig reference at 1, 2
and 4 threads and compares the files byte for byte, plus the tallies. The contigs differ in size on
purpose, so the workers finish out of order and the writer has to put them back.

## Validation

In the dev container: `cargo fmt` clean; `cargo clippy --lib --tests --all-features -- -D warnings`
clean; `cargo test --lib` **2,924 passed, 0 failed, 5 ignored**.

One test assertion was loosened during the run, not the code: the 2 kb fixture asserted an exact
2,000-base tract, but the detector reaches into the flanking filler when its bases happen to continue
the tiling. The assertion is now "one row, at least 2,000 bases", which is what the test is actually
about.

**Pre-existing and untouched:** `cargo clippy --all-targets` still fails on
`examples/ng_inbreeding_harness.rs` (`076cb5e9`).

## Tradeoffs and follow-ups

- A parallel build holds `threads` contigs resident. On tomato that is ~90 MB per worker, on human
  ~250 MB; the CLI in E1 is where a default gets chosen, and E2 measures the real peak.
- `scan_contig` re-allocates a row vector per contig rather than reusing one, because the batch's
  results must all exist before any is written. The allocation is one per contig, against a scan that
  touches every base.
