# ng cohort CRAM ingest — one reference repository per run, not per file

**Date:** 2026-07-29 · **Branch:** `ng-cram-repository-sharing` · **Symptom:** a 51-CRAM
cohort stutter dump is OOM-killed (exit 137) at ~80 s against `scripts/dev.sh`'s 16 GB cap.

## The measurement, first

Peak RSS against file count, `ng_ssr_cohort_stutter` over `tmp/tomato_stutter/ssr_subset.bed`
(556 spans, all 12 tomato chromosomes), tomato reference `S_lycopersicum_chromosomes.4.00.fa`.
Harness: `tmp/mem_curve.py`, reading the child's own `/proc/<pid>/status` `VmHWM`.

| files | peak RSS before | peak RSS after | wall before | wall after |
| ----: | --------------: | -------------: | ----------: | ---------: |
|     1 |       835.1 MiB |      840.1 MiB |      24.0 s |     24.4 s |
|     2 |      1590.3 MiB |      857.4 MiB |      27.5 s |     27.7 s |
|     4 |      3102.0 MiB |      883.9 MiB |      35.2 s |     34.4 s |
|     8 |      6105.8 MiB |      946.9 MiB |      49.4 s |     47.0 s |
|    16 |     12125.8 MiB |     1063.7 MiB |      83.9 s |     79.9 s |
|    32 |      *OOM*      |     1307.0 MiB |       —     |    149.3 s |
|    51 |      *OOM*      |     1590.6 MiB |       —     |    209.6 s |

Before, the curve is a straight line: **752.7 MiB per open file** on an 82 MiB base, fitting every
measured point to better than 0.3%. Sample count is the driver, exactly as reported — and the
per-file constant is not merely "large", it is **the size of the reference**: the FASTA's 12
contigs total 782,520,033 bases = 746.3 MiB, which is 99.2% of the measured slope.

Extrapolated, 51 files ask for 37.6 GiB and cross the 16 GB cap at about the 22nd file — mid-walk,
after rows have been emitted, which is what was observed.

After, the slope is **15.1 MiB per file**, a 50× reduction in the per-file term, and the whole
51-sample cohort peaks at 1.59 GiB — about what *two* files cost before.

## The cause

Confirmed as diagnosed. `AlignmentFile::open` built a `noodles_fasta::Repository` per alignment
file, for CRAM only (`open_bam.rs`; BAM stores its own sequences and needs no reference, which is
why BAM cohorts never hit this). That repository is a whole-contig memoising cache with **no
eviction path** — `noodles-fasta-0.61.0/src/repository.rs` has `get`, `len`, `is_empty` and
`clear`, and `get` inserts the entire contig — so each open CRAM accumulates its own copy of every
contig it decodes against.

The cohort walk is region-outer and sample-inner, and the BED spans all 12 chromosomes, so every
file ends up holding every contig: `files × genome`.

## The fix

`src/ng/read/input/reference.rs` — a new `RunReference`: the reference description a run opens
every file against, **plus the one shared cache of its bases**.

- `AlignmentFile::open` and `SampleReads::open` / `open_only_sample` now take `&RunReference`
  where they took `&ReferenceInfo`. The repository is *asked for*, not built, so there is no
  per-file arm in which a second one could appear.
- The bases open **lazily**, on the first CRAM open, so a BAM-only run against a FASTA with no
  `.fai` still works — behaviour that an eager field would have broken.
- Only successes are memoised, so a transient failure to open the FASTA can retry — the policy
  `ReferenceInfoCache` already states for its own reads.
- `ReferenceBasesError` keeps the two faults apart (`.fai`-only reference vs. unopenable FASTA);
  `open` maps them onto the same `CramNeedsReferenceFasta` / `Open` errors it raised before, so no
  message changed.

`ReferenceInfo` itself is untouched. It is documented as a plain data record whose `fasta_path` is
"provenance passed through for consumers that need to go back to the bases"; hanging a 746 MiB
cache off it would contradict that, break its `Clone`/`PartialEq`/`Eq` derives, and make the
process-lifetime `ReferenceInfoCache` pin a genome. The owner of the bases is the run, so the type
that owns them is the run's.

## Verification

**Output is byte-identical.** Same BED, same sample sets, pre-change release binary vs. post-change:

| samples | rows  | md5 before = after |
| ------: | ----: | :----------------- |
|       8 | 47462 | `62d59ee6…` ✓      |
|      16 | 98509 | ✓                  |

(The precedent for this check is `d422e64`, which verified the CRAM reader-reuse optimisation the
same way on three samples.)

**Tests:** `cargo test` — 2654 lib tests + all integration targets, 0 failures. `cargo clippy
--all-targets` clean apart from pre-existing feature-gated dead-code warnings in
`examples/dhat_ng_merge.rs`.

**No perf regression.** The `DecodedContainer` container cache (`b918fb6`, `region_query.rs`) is a
different cache and is untouched — the only diff in `region_query.rs` is in its tests. Wall time
improves monotonically with file count (−5% at 8, −5% at 16) because each file no longer re-reads
the genome into a cache of its own.

## The two open questions, answered

**Is a bounded cache warranted?** Not for tomato, and not yet. One shared repository still holds
every contig it touches — 746 MiB here, ~3 GiB for a human reference, more for a polyploid crop.
Because the cohort walk is region-outer, every file is on the same contig at the same time, so
production's existing policy — `Repository::clear()` on contig transition, described in
`bam::alignment_input::build_fasta_repository` — would cap the resident bases at the largest single
contig, chr01 at 86.6 MiB: a further 8.6× on this reference. That is the shape `WindowedRefSeq`'s
`evict_before` already has for the walk's own view of the bases.

It is deliberately **not** done here, because eviction needs a signal the input layer does not
have: only the caller knows when the walk has left a contig for good, and clearing early just
re-reads. Adding it means an evict-on-transition entry point on `RunReference` *and* a caller that
calls it — a change to the walk's contract. Worth doing when `genome` is the problem; here
`files × genome` was.

**Does anything else scale per open file?** Yes, but marginally: the residual 15.1 MiB per file is
the reader pool, the parsed index, the `.crai`-by-contig table and the per-worker decoded CRAM
container. At the 64-file cohort now being planned that is ~1.0 GiB of the ~1.8 GiB total, and at
the archive's 2,085 files it would be ~31 GiB — so it is the *next* thing to measure if the file
count grows another order of magnitude, but it is not what was killing this run.

**Thread contention.** Sharing the cache shares its lock, and noodles takes the repository's
*write* lock across the adapter read, so the first fetch of a contig blocks every other reader of
that repository until it is in. That is a convoy of one whole-contig read per contig per run —
against per-file repositories, where each file paid that read in full, total work drops even as the
blocking appears. Every fetch after the first takes the read lock and is concurrent. Today's
callers are single-threaded, so this costs nothing yet; re-time when one is not.

## Consequence for the cohort workflow

`tmp/tomato_stutter/run_batches.sh` and
`benchmarks/ssr_tomato1/scripts/combine_stutter_batches.py` are no longer needed for a cohort of
this size: 51 samples run in one process, so read-group identifiers are minted once and need no
renumbering on `(file, @RG ID)`. The ~64-file cohort fits comfortably (projected ~1.8 GiB).
