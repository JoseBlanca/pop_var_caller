# ng cohort CRAM ingest — one reference repository per run, and one contig of it

**Date:** 2026-07-29 · **Branch:** `ng-cram-repository-sharing` · **Symptom:** a 51-CRAM
cohort stutter dump is OOM-killed (exit 137) at ~80 s against `scripts/dev.sh`'s 16 GB cap.

## The measurement, first

Peak RSS against file count, `ng_ssr_cohort_stutter` over `tmp/tomato_stutter/ssr_subset.bed`
(556 spans, all 12 tomato chromosomes), tomato reference `S_lycopersicum_chromosomes.4.00.fa`.
Harness: `tmp/mem_curve.py`, reading the child's own `/proc/<pid>/status` `VmHWM`.

Two fixes, measured separately: sharing one repository across the run, then bounding that
repository to the contig in hand.

| files | before     | one repository | + one contig |
| ----: | ---------: | -------------: | -----------: |
|     1 |  835.1 MiB |      840.1 MiB |    180.4 MiB |
|     2 | 1590.3 MiB |      857.4 MiB |    186.4 MiB |
|     4 | 3102.0 MiB |      883.9 MiB |    218.3 MiB |
|     8 | 6105.8 MiB |      946.9 MiB |    280.8 MiB |
|    16 |12125.8 MiB |     1063.7 MiB |    395.3 MiB |
|    32 |    *OOM*   |     1307.0 MiB |    638.7 MiB |
|    51 |    *OOM*   |     1590.6 MiB |    908.9 MiB |

Before, the curve is a straight line: **752.7 MiB per open file** on an 82 MiB base, fitting every
measured point to better than 0.3%. Sample count is the driver, exactly as reported — and the
per-file constant is not merely "large", it is **the size of the reference**: the FASTA's 12
contigs total 782,520,033 bases = 746.3 MiB, which is 99.2% of the measured slope.

Extrapolated, 51 files ask for 37.6 GiB and cross the 16 GB cap at about the 22nd file — mid-walk,
after rows have been emitted, which is what was observed.

Sharing the repository turns that per-file 752.7 MiB into **15.1 MiB**, a 50× cut, leaving a fixed
~825 MiB — the genome, resident once. Bounding it to one contig removes that fixed term too:
the base falls to **~166 MiB** (one chromosome, chr01 at 86.6 MiB, plus the process) while the
slope stays where sharing left it, 14.6 MiB per file. That is the signature of having removed
exactly the reference and nothing else, twice over.

51 samples: **37.6 GiB → 909 MiB** (908–914 MiB across four runs), and the run completes.

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

`src/ng/read/input/reference.rs` — a new `OpenReference`: the reference description a run opens
every file against, **plus the one shared cache of its bases**. "Open" is the same "open" as
`AlignmentFile::open` — a `ReferenceInfo` *describes* a reference, an `OpenReference` can be *read
from*. Holding both in one value is not tidiness: they must be the same reference, and two
arguments could be mismatched, decoding every read against the wrong bases with nothing to catch
it.

- `AlignmentFile::open` and `SampleReads::open` / `open_only_sample` now take `&OpenReference`
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

**Output is byte-identical.** Same BED, same sample sets, pre-change release binary vs. post-change
— checked after each half of the fix:

| samples | rows  | md5 before = after |
| ------: | ----: | :----------------- |
|       8 | 47462 | `62d59ee6…` ✓      |
|      16 | 98509 | ✓                  |

(The precedent for this check is `d422e64`, which verified the CRAM reader-reuse optimisation the
same way on three samples.)

**Tests:** `cargo test` — 2656 lib tests + all integration targets, 0 failures, including two new
ones for the contig bound (a transition drops the previous contig; staying on one does not).
`cargo clippy --all-targets` clean apart from pre-existing feature-gated dead-code warnings in
`examples/dhat_ng_merge.rs`.

**No perf regression found.** The `DecodedContainer` container cache (`b918fb6`,
`region_query.rs`) is a different cache and is untouched — the only diff in `region_query.rs` is in
its tests. Sharing the repository *improves* wall time (−5% at 8 and 16 files), since each file no
longer re-reads the genome into a cache of its own.

The contig bound should cost nothing on top, and the structural reason is that the BED has exactly
12 contig runs: 11 transitions, 12 whole-contig reads, whether or not the cache is cleared between
them. The timings neither confirm nor contradict that at the resolution measured. Three clean
51-sample runs of the bounded build gave 225.0, 226.5 and 215.0 s (spread 11.5 s); the one clean
run of the sharing-only build gave 209.6 s. A cost of a few percent cannot be ruled out from that,
and settling it needs repeats of both builds interleaved. **Not established as free — established
as not large.**

## Second half — one contig resident, not the genome

Holding the whole genome was ours too, not noodles'. What noodles needs at any moment is the
contig of the slice it is decoding: `get_slice_reference_sequence`
(`noodles-cram-0.93.0/src/io/reader/container/slice.rs:352`) asks the repository for **one contig by
name** and keeps it as `ReferenceSequence::External { sequence: Arc<Sequence> }` for that slice.
Everything beyond that is cache policy, and an unevicted `Repository`'s policy is "keep every
contig ever touched" — over a genome-wide walk, the genome.

`OpenReference::bases_for_contig(contig)` is now what the query path asks, and it clears the cache
when the run moves to a new contig. `AlignmentFile` holds the `OpenReference` rather than a bare
repository, so the narrowing happens per region query with the contig in hand.

- **No extra reading in an ordered walk.** A contig is cleared only when a query for a *different*
  contig arrives, so chr1 → chr2 → … reads each contig exactly once, as before. A cohort walk is
  region-outer, so the transition happens once per contig for the whole cohort, not once per file.
- **Safe under concurrency and against the container cache.** `clear()` drops map entries only;
  noodles holds each contig by `Arc` for as long as it is decoding against it, and a cached
  `DecodedContainer` pins the same single allocation, so an in-flight or cached decode is never
  broken — only the *next* fetch re-reads.
- **The assumption is contig-ordered access**, which every caller in this tree makes. A
  contig-*parallel* walk would clear a chromosome a neighbour is still reading and re-read it
  repeatedly; such a caller builds its reference with `OpenReference::unbounded` and will want a
  k-slot policy (one repository per worker contig) instead.

## The remaining open question

**Does anything else scale per open file?** Yes, and with the reference gone it is now the whole
curve. Peak memory at 51 files splits into a fixed ~166 MiB — one chromosome plus the process — and
14.6 MiB for every file opened, which at 51 files is 743 MiB, or 82% of the 909 MiB total. That
per-file term is what now sets the ceiling: ~1.1 GiB at the 64-file cohort, ~30 GiB at the
archive's 2,085.

**What that 14.6 MiB is has not been measured.** The per-file state `AlignmentFile` holds is the
reader pool with its filter buffers, the parsed index, the `.crai`-by-contig table, and the decoded
CRAM container each pooled reader keeps. The indexes are ruled out by inspection — these `.crai`
files are ~1.6 KB each — which leaves the decoded container as the likely dominant term, since it
holds a whole container's worth of decoded records. That is a hypothesis, not a result; one heap
profile would settle it, and it is the obvious next measurement.

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
renumbering on `(file, @RG ID)`. The ~64-file cohort projects to ~1.1 GiB.
