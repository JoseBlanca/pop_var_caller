# ng direct mode, step D3 — genotypes from real reads, and where the time goes

**Date:** 2026-09-01. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step D3, and **Checkpoint D**. **Spec:**
[`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §11 question 7, §12 oracle 3.
**Modules:** `examples/ng_call_cohort_end_to_end.rs`, `src/ng/run/mod.rs`.

---

## ng calls genotypes from CRAM files

Six tomato accessions from `benchmarks/tomato1/crams/` over 400 kb of SL4.0, through the real
repeat catalog, in the development container:

```
# segments: 5469 over 4 analysed region(s)
# samples: 6
# loci called: 8411
# loci the merge declined to assemble for being too wide: 0
# alleles a locus is called over, mean: 1.69 (the reference counts as one)
# assembly check: every one of 6 alignment file(s) matched the reference; 78 of 78 contig
#   checksums could be compared and all agreed
```

**Two things had never happened before.** Every test behind `call_cohort` is a fabricated BAM
of three or four reads over a hundred bases of a reference that is one homopolymer; this is
sequenced DNA, a real catalog, and ground with repeat tracts, gaps and contig changes in it. And
the assembly check had only ever run against fixture headers — here it compared all 78 contig
checksums across six CRAMs and they agree.

**6.0% of the analysed ground is not called, and the run says so rather than calling it
wrongly.** 2,733 of 5,469 typed regions are repeat tracts — half of them — but **24,053 bases
of the 400,000**, which is the number that matters and is why the probe reports the share in
bases. They are charged to *not built yet* and kept apart from *never called (satellite)*,
which is zero here. Repeat-tract candidate selection
is specified and unbuilt, so this is the plan's own scope decision showing up as a number.

## Where a calling run's time goes — spec §11 question 7

Measured over the first two BED intervals, 200 kb at about three reads a position, release
build with `--features merge-timing`, in the development container:

| samples | compressed MB | loci called | `call_cohort` | drawing the readers | assembling | genotyping |
|---|---|---|---|---|---|---|
| 3 | 123.5 | 3,291 | 0.61 s | 589.2 ms (97.3%) | 7.4 ms (1.2%) | 5.8 ms (1.0%) |
| 6 | 216.8 | 4,235 | 1.07 s | 1,040.7 ms (96.8%) | 14.9 ms (1.4%) | 12.9 ms (1.2%) |
| 12 | 383.1 | 5,675 | 1.98 s | 1,895.9 ms (96.0%) | 32.9 ms (1.7%) | 34.1 ms (1.7%) |
| 24 | 672.8 | 8,825 | 3.89 s | 3,671.6 ms (94.3%) | 89.6 ms (2.3%) | 103.1 ms (2.6%) |

**Decoding reads is `call_cohort`.** Assembling and genotyping together are 2.2% of it at three
samples and 4.9% at twenty-four.

**`call_cohort` is not the whole wait, and the first draft of this probe did not say so.**
Reading and checksumming the 795 MB reference, opening the catalog and building the segments
cost **2.73 to 2.81 seconds** across those four runs — constant in the cohort and in the ground
— so at the defaults they are more than half of the 4.8 seconds a person waits. They are now a
row of the output.

**The split is only visible because calling happens inside the builder.** A run's own stopwatch
cannot see inside a merge that returns everything at once; D1's per-locus counter is what
separates the assembling from the genotyping.

## ⚑ What this says about Milestone E

**Milestone E, as the plan describes it, parallelises the 5%.** Both arrangements it weighs —
the merge's own region batching switched on, or the merge left on one thread handing each
finished locus to genotyping workers — speed up *assembling and genotyping* and nothing else.
On these cohorts that is bounded above by 4.9% of `call_cohort`.

**What costs the run is `ObservationCache::cover`, and it is one thread.** The driver a calling
run uses draws every sample forward one after another. `cover_in_parallel`, which sweeps the
cohort's samples concurrently and is documented as reaching the same fixpoint by a different
schedule, **already exists** — and is reached only by the merge's parallel driver, which
`call_cohort` does not use. Measured at three samples: **3.199 s of user CPU against 3.313 s
elapsed**, which is one core of the nine the container had.

**Two rates are stable and one share is not, and the difference is the whole of the finding:**

- **reading costs about 5 ms per compressed megabyte** — 4.77, 4.80, 4.95, 5.46 — so it is
  linear in the bytes opened to within 14% across the range;
- **calling costs about 1 µs per locus per sample** — 1.34, 1.09, 0.98, 0.91 — flat, and
  falling;
- **calling's *share* grows from 2.2% to 4.9%, and the whole of that is the locus count.** More
  accessions segregate more sites, so the count goes 3,291 to 8,825 while the cost of each
  falls.

**That last curve has to flatten** — 200 kb of SL4.0 holds a finite number of segregating sites
— and where it flattens decides whether calling is a tenth or a third of a thousand-sample run.
**This probe cannot see it.** An earlier draft of this report fitted two exponents to the four
cohorts and said "roughly a fifth at a thousand samples"; the review showed the exponents drift
by half within the fitted range itself (calling grows 2.06×, 2.44×, 2.87× per doubling), so
three defensible models give a tenth, a fifth and a third. The honest output is the two stable
rates and the statement that the share turns on locus discovery, which needs a run at a
thousand.

**So the recommendation is not "skip Milestone E".** It is that the milestone's first question
should be whether a calling run may use the cover that already spreads across samples, and its
second — which calling arrangement to build — be settled at the cohort size the caller is meant
to serve rather than at six. **One thing already measured bears on the second**: the two costs
charged per building region rather than per locus, building the per-sample windows and setting
each region's walk up, are 0.0 to 0.3 ms against 7 to 90 ms of assembling, so on this ground
dividing the genome more finely is nearly free.

## What the probe is, and what it is not

**An example rather than a test**, following `examples/ng_open_cohort_descriptors.rs`, which is
Checkpoint A's evidence in the same shape. The reference lives on the container's read-only
`$HOME/genomes` mount, so nothing that needs it can be a test the whole suite runs.

**It defaults to a handful** — six samples over four intervals — because it is meant for the
development loop. `NG_SAMPLES` and `NG_REGIONS` are how the whole cohort is asked for
deliberately.

**It writes no VCF.** The format is settled and the writer is coded; what turns a called locus
into a record needs per-sample evidence that has outlived the locus, and assembling that is not
this step.

**The catalog is opened from an explicit path**, not found beside the reference: `$HOME/genomes`
is mounted read-only, so no catalog can be written next to it.

## Verification

| check | result |
|---|---|
| `cargo test --lib` | 5,813 passed, 13 ignored — unchanged, as the step adds no library code beyond a re-export |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |

**The library change is three names.** `CalledCohort`, `CohortWalkTallies` and
`SampleWalkTallies` join the run module's existing re-exports, so a consumer outside the crate
can name what `call_cohort` returns.

## What the runs and the review changed

**The breakdown printed `0.00` on three of its four rows**, because it was in seconds to two
decimals and the whole run was under a second. A split whose only purpose is to be compared has
to be printed at a scale that distinguishes its parts; it is milliseconds and a share now.

**The review re-measured the table and it reproduced** — every locus count to the digit, every
time within 2%. What it overturned was the *interpretation*, twice:

- **A scaling claim that was an artefact of file sizes.** The draft said reading grows as about
  `n^0.9` in the cohort. The 24 CRAMs range from 8.1 MB to 61.5 MB, and the first three in name
  order are 1.5 times the average, so eight times the samples is only **5.45 times the bytes**.
  Per megabyte the reading is flat. The table now carries the megabytes and the probe prints the
  rate.
- **An extrapolation stated as an answer.** See above: the fitted exponents are not stable
  inside the range they were fitted on, and the superlinearity is locus discovery rather than
  per-locus cost — which the per-locus rate, falling from 1.34 to 0.91 µs, says plainly.

**And about 2.7 seconds of every run was timed nowhere.** The probe began its clock at
`call_cohort`, so it described a minority of what a person waits for while saying "decoding
reads is the run". Both are fixed: the setup is a row, and the claim is scoped to
`call_cohort`.

**Four things a person reads were wrong.** The assembly check printed as a Rust `Debug` struct,
which is the one line answering *were these files aligned to this reference*; three per-sample
numbers were printed as though they partitioned when the third is a subset of the second; the
word *regions* meant three different things in one page of output (analysed intervals, typed
regions, the merge's working windows); and the unbuilt ground was reported as a region count —
half the regions — where the number that matters is the base share, **6.0%**.
