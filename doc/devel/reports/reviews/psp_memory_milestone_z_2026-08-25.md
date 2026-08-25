# Where a cohort run's memory actually goes, and what a store redesign can reach

*Measured 2026-08-25, branch `ng-psp-encoding`. Milestone Z of
[`../../ng/impl_plan/psp_encoding_experiments.md`](../../ng/impl_plan/psp_encoding_experiments.md) —
the step that runs before the encoding experiments to find out how much of the peak they could
possibly move.*

---

## What was asked and what came back

The plan proposes three changes to the per-sample store — larger blocks decompressed
incrementally, the three floating-point quantities stored as counts of a step, and the read names
stored as changes — and the objective is to cut the memory a cohort run needs. Before spending
weeks on encodings, Milestone Z asks a cheaper question: **of the memory a run actually uses, how
much is the store's to give back?**

Three things came back, and the third was not being looked for.

1. **A cohort run's peak is almost entirely a per-sample cost.** Measured over cohorts of 1 to 62
   samples, peak resident memory is `91 MB + 7.53 MB per sample`, a straight line to within half a
   per cent. At a thousand samples the fixed part is one per cent of the total. **Whatever we do
   about memory has to move the per-sample slope; nothing else matters at cohort scale.**
2. **The store's decode buffers are a small part of that slope.** At 50 samples the compressed
   block decode holds 43.6 MB of a 342.6 MB live heap — 12.7%. The larger masses are the per-sample
   columns the merge assembles (29.6%), the merger's own projections (23.2%), and one thing nobody
   had counted.
3. **A quarter of the per-sample cost is the `.psp`'s metadata section, and two thirds of that is a
   duplicate.** The summary it carries is needed and stays; what is held for nothing is the *text it
   was parsed from*, kept alive beside the parsed structure for the rest of the run. Dropping the
   text takes the measured slope from 7.53 to 6.50 MB per sample — **13.7% off the per-sample cost
   of every cohort run** — with identical output.

---

## 1. The shape of the peak

`scripts/cohort_memory_vs_samples.sh` runs the same calling job at six cohort sizes and records
peak resident memory. Tomato at about three reads a position, 62 accessions, four threads, inside
the dev container.

| samples | peak resident |
|---:|---:|
| 1 | 109.3 MB |
| 5 | 136.5 |
| 10 | 160.9 |
| 25 | 260.3 |
| 50 | 458.1 |
| 62 | 573.0 |

A straight line through these is **`91.1 MB + 7.527 MB × samples`**, R² 0.995.

**Read the two numbers separately, because only one of them matters.** The 91 MB is what the run
needs whatever the cohort size — the reference, the low-complexity mask, the fixed working set. The
7.53 MB is paid again for every sample added. At 62 samples the fixed part is a sixth of the total;
at a thousand samples it would be one per cent.

*Extrapolating the line gives 7.4 GB at a thousand samples and 22 GB at three thousand. That is
extrapolation from measurements that stop at 62, on one species at one depth, and the line is very
slightly convex at the top — the last two points sit above it. Treat those as the right order of
magnitude, not as measurements.*

## 2. What holds the memory, at the instant the peak occurs

A heap profile records how many bytes each allocation site holds at the instant the whole
program's live heap is largest. At 50 samples, grouped by which part of the caller allocated them:

| what holds it | at peak | share |
|---|---:|---:|
| per-sample columns assembled for the merge | 101.6 MB | 29.6% |
| the merger's per-group projections | 79.6 | 23.2% |
| **the `.psp` metadata section** | **79.1** | **23.1%** |
| compressed block decode | 43.6 | 12.7% |
| the genotype fit | 27.3 | 8.0% |
| reference spans behind the low-complexity filter | 6.5 | 1.9% |
| everything else | 4.9 | 1.4% |
| **total live heap** | **342.6** | |

Live heap is smaller than resident memory — 342.6 MB against 458.1 MB for the same run — because
the allocator holds pages it is no longer using. **Use the shares to compare parts against each
other and the resident figure for totals; the two are not interchangeable.**

### 2.1 A trap in reading a profile like this, which caught this one

A heap profile's peak is **one instant**, and which instant it is changes with the cohort size. At
10 samples the largest moment of the run falls during the low-complexity mask, where a reference
span copy holds 56% of a much smaller heap; at 50 samples it falls during calling. **Two profiles at
two cohort sizes are therefore not two points on one curve unless both peak in the same phase**, and
comparing them directly was wrong the first time it was tried here.

The one figure that survives the comparison is the one that is identical in both: the metadata.

## 3. The finding: the summary is kept twice — parsed, and as the text it came from

Every `.psp` carries a metadata section — the per-sample summary the genotype prior and the
hidden-paralog filter read, chiefly a coverage-against-GC histogram, stored as text and compressed.
`PspReader` decompresses it when the file is opened and keeps the decompressed bytes in a field for
as long as the reader lives, which is the whole run.

**It costs exactly 1.58 MB per open sample, and the two profiles agree to three digits:**

| cohort | metadata at peak | per sample |
|---:|---:|---:|
| 10 samples | 15.8 MB | 1.58 MB |
| 50 samples | 79.1 MB | 1.58 MB |

It splits in two, and **only one of the two is waste**:

- **0.52 MB per sample — the parsed summary.** A coverage-against-GC histogram and heterozygosity
  counts. **Needed, and needed more later, not less.** The genotype prior and the hidden-paralog
  filter read it today, and ng's own duplication filter will read it when that is built. Nothing
  here proposes touching it.
- **1.05 MB per sample — the text the summary was parsed from.** Read in exactly one place in the
  whole calling path, `pipeline.rs`, at startup, to produce the structure above. Nothing reads it
  again, and the reader holds it until the run ends. **So the sample summary is resident twice for
  the length of the run: once as a parsed histogram and once as the TOML it was decoded from.**

**The waste is the duplication, not the information.** A store that dropped this section would break
the hidden-paralog filter, which hard-fails rather than silently emitting an unfiltered callset when
a sample's summary is missing (`require_paralog_summaries`). That the filter ran, and that the VCF
came back identical, is the evidence that the experiment below keeps everything both consumers
need.

**Against the working budget of 500 kB per open sample this is more than three times over, before a
single record has been read** — and none of it is touched by block size, record layout, read-name
encoding, float precision, or read depth. It is the cohort size and nothing else.

### 3.1 What releasing it is worth, measured rather than argued

Taking the bytes instead of borrowing them — eleven lines, committed on this branch as an
experiment and not proposed for `main` — moves the whole curve:

| samples | main | released | saved |
|---:|---:|---:|---:|
| 1 | 109.3 MB | 106.6 MB | 2.5% |
| 10 | 160.9 | 152.9 | 5.0% |
| 25 | 260.3 | 228.4 | 12.3% |
| 50 | 458.1 | 418.7 | 8.6% |
| 62 | 573.0 | 504.5 | 12.0% |

| | fixed | per sample | R² |
|---|---:|---:|---:|
| main | 91.1 MB | 7.527 MB | 0.995 |
| metadata released | 92.6 MB | 6.495 MB | 0.992 |

**The fixed part does not move and the per-sample part falls by 1.03 MB — 13.7%.** That is the shape
a per-open-sample cost has, and the 1.03 MB agrees with the 1.05 MB the heap profile attributed to
the same site by a completely different method.

The emitted VCF is unchanged: identical record counts at all six cohort sizes. 4,532 library tests
pass. *(Ten `should_panic` tests fail under `--release` on this branch and on `main` alike — they
assert a `debug_assert` firing, which release builds strip. Not related to this change.)*

## 4. What this means for the three encoding experiments

**The ceiling is real but it is not where the plan assumed.** Of the per-sample slope, the part the
*encoding* of records controls — the compressed block decode — is 12.7% of live heap at 50 samples.
The parts a *store redesign in the wider sense* controls are larger: the metadata is 23.1%, and the
per-sample columns the merge assembles are another 29.6%, which is what a reader that hands out one
observation at a time would change.

So, in the order the owner asked for — simple things first, complexity only if measurement forces
it:

1. **Release the metadata.** Measured, 13.7% of the per-sample cost, output unchanged, eleven lines.
   Nothing in the encoding plan is a prerequisite. **Recommended, and it is the owner's call because
   it touches production.**
2. **Stop storing that section as text.** A count matrix written as TOML is several times its own
   size, on disk and again in the parser — 1.05 MB of text for a 0.52 MB histogram. **This gets more
   worth doing rather than less**, because ng will read this summary too: its own duplication filter
   is not built yet but is coming, so the section is permanent and its encoding is worth getting
   right. Unmeasured, and the next cheap thing to measure.
3. **Then the encoding experiments**, whose reachable share is now known rather than assumed.

## 5. What is not yet measured

- **The block-window sweep** — how much of the peak follows the block's span in reference
  coordinates. The script is written (`scripts/psp_block_window_sweep.sh`); rewriting 50 samples at
  three window sizes is the expensive part and it has not been run.
- **Anything above 62 samples.** Every figure here is 1–62 on one species at about three reads a
  position. The committed range goes to thousands of samples and to 300 reads a position, and the
  line's slight convexity at 62 is a reason to measure rather than extrapolate.
- **Whether the per-sample columns behave the same way.** They are the largest single mass and the
  one a streaming reader is supposed to remove; nothing here measured them against cohort size.
- **Wall time.** Two containers from other sessions were running throughout, on a machine with six
  fast cores. The wall-time column in these results is not usable and no conclusion here rests on
  it.

## 6. How to reproduce

```
scripts/cohort_memory_vs_samples.sh \
    --psp-dir  <cohort psp dir> \
    --reference <reference.fa> \
    --out-dir  tmp/milestone_z/nsweep \
    --sizes 1,5,10,25,50,62 --threads 4

./scripts/dev.sh cargo run --release --features dhat-heap --target-dir target-dhat \
    --example dhat_var_calling -- \
    --psp-dir <cohort psp dir> --n-samples 50 --reference <reference.fa> \
    --output tmp/milestone_z/n50.vcf --threads 4 --target-variants-per-chunk 128

uv run tmp/milestone_z/attribute_peak.py tmp/milestone_z/n50_default_dhat.json
```

**`--target-dir target-dhat` is not optional.** Without it the instrumented binary overwrites the
release one at the same path and every later "release" run silently re-executes it five to six times
slower, with nothing in the output saying so. That has already cost one performance review its whole
measurement set.
