# ng repeat catalog — milestone E: the command, and what a real genome costs

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **E1** and **E2**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §2.6, §9.1, §10.4.*

## What the command is

```
pop_var_caller_exp repeat-catalog --reference ref.fa [--output ref.fa.repeats.parquet]
                                  [--threads N] [--force] [scan knobs]
```

It scans a reference once, writes the catalog beside it, and prints what it found. The defaults are
the catalog's — periods 1 to 6, copy floors `[5, 5, 4, 4, 4, 3]`, 15 bp of flank — not step 3's
calling defaults. An existing catalog is left alone unless `--force` says otherwise.

## The measurement (spec §9.1, §10.4)

Two references, built at the defaults. **Wall clock and peak memory are from the Linux dev container**
(`./scripts/dev.sh`, 16 GB / 8 CPU); file sizes are the files on disk.

| | tomato SL4.00 | GRCh38 no-alt + hs38d1 |
|---|---|---|
| FASTA | 795 MB, 13 contigs | 3.15 GB, 2,580 contigs |
| **catalog** | **58 MB** | **213 MB** |
| catalog / FASTA | 7.6% | 7.1% |
| repeats written | 6,375,702 | 23,569,229 |
| wall clock, 1 thread | **30 s** | **110 s** |
| wall clock, 4 threads | **13 s** | not measured |
| peak RSS, 1 thread | **1.49 GB** | **3.89 GB** |
| largest contig | 90 Mb | 248 Mb |

**Repeats by period**, and the shape is the same in both genomes — homopolymers are almost all of it:

| period | tomato | GRCh38 |
|---|---|---|
| 1 | 6,122,751 (96.0%) | 22,324,775 (94.7%) |
| 2 | 107,841 | 556,027 |
| 3 | 95,956 | 258,173 |
| 4 | 18,372 | 224,904 |
| 5 | 5,750 | 98,879 |
| 6 | 25,032 | 106,471 |

**What the floors dropped**, which is the other half of the tally: 196,485,350 detections below the
copy floor on tomato and 720,848,965 on GRCh38 — **31 dropped for every one kept**. The flank floor
dropped 7 tracts on tomato and 857 on GRCh38, all of them within 15 bp of a contig end.

### Open question 1 is answered: the file is small

The spec worried that a permissive scan might make an unaffordable file, and hung the format choice on
it. It did not: **58 MB for tomato and 213 MB for human, about 7% of the FASTA either way, and
9 bytes per stored repeat.** Both are far under the 1 GB trigger §3.5 set, so Parquet stays and no
re-costing is needed.

### `--threads` is free, and byte-identical on a real genome

Tomato at 4 threads is **13 s against 30 s**, and `cmp` says the file is **byte-identical** to the
sequential build. The unit test asserts that on fixtures; this is the same property on 13 real
chromosomes finishing out of order.

## The memory number, and why the first one was wrong

**Measured on macOS, the same tomato build peaks at 10.5 GB. Measured in the Linux container, 1.49 GB.**
Same binary, same input, same output file. The gap is the platform allocator, not live memory: the
scan makes many short-lived allocations per contig, and macOS's allocator holds the pages rather than
returning them, so `/usr/bin/time -l`'s "maximum resident set size" counts them all.

Two things follow. **The Linux figure is the one that means anything** — it is where the archive-scale
runs happen. And **a macOS developer building a catalog should expect the process to look like it is
using ten gigabytes**, which is worth knowing before it alarms someone.

The per-period breakdown, on macOS and therefore inflated, still says where the memory goes: period 1
alone peaks at 3.7 GB, period 2 alone at 8.8 GB, period 6 alone at 7.3 GB. It is **the scan's own
working set over a whole contig**, not the bases and not the rows — a single period costs more than
all six of them cost in rows.

**This falsifies a number in the spec, and the correction is a clean law.** §2.3 says the
whole-contig scan costs "the contig resident while it is scanned — 90 MB for tomato's largest
chromosome, 250 MB for human chromosome 1". The bases are indeed that; the scan around them is not:

| largest contig | peak RSS | per base |
|---|---|---|
| tomato, 90 Mb | 1.49 GB | 16.6 bytes |
| GRCh38, 248 Mb | 3.89 GB | 15.7 bytes |

**About 16 bytes of working memory per base of the largest contig**, not one — the two genomes agree
to within 6%, so this is a rule and not a coincidence. The design decision stands (3.9 GB is
comfortable on the machines this runs on), but the spec's cost statement should say 16 bytes a base
and name the largest contig as what sets the peak. `--threads N` multiplies it by N, which is what
makes the default of 1 the safe one.

## Validation

`cargo fmt` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean apart from two
pre-existing failures in `examples/` (`ng_inbreeding_harness`, `ng_multilib_key_harness`, both
committed before this work); `cargo test --lib` **2,934 passed**; the differential **6 passed**.

Both catalogs were built and are readable; the tomato one opens against its own reference and its rows
carry the motifs the scan found.

## Follow-ups

- **Correct spec §2.3's memory statement** to the measured 1.49 GB (tomato, one thread), and note the
  macOS allocator effect where a reader would otherwise be alarmed.
- **The `--threads` default is 1.** At 4 threads tomato halves its wall clock and multiplies its peak
  by up to 4; whether the command should default higher is a choice about the machines it runs on.
- **Not measured:** the reference pass without the scan, so "what the scan adds" is not isolated. The
  comparison that matters more — one build against a scan per sample — is arithmetic: 30 s once
  against 30 s × N.
