# Performance Review: psp_to_vcf
**Date:** 2026-11-05
**Reviewer:** rust-performance-review skill (orchestrator, single reviewer, no fan-out)
**Scope:** `call-from-psps` on a cohort, the hidden-duplication filter and the neighbour-evidence pass on
**Verdict:** Apply the listed wins (applied on branch `perf-psp-vcf`)
**Hot-path evidence:** `samply` profiles with line attribution through inlining (`atos -i`), `/usr/bin/time -l`, resident memory sampled every half second, timers in the round loop

---

## 1. Scope and constraints

**What was reviewed.** `call-from-psps`, the command that reads stored per-sample evidence files
(psps) and writes the VCF, on 59 tomato accessions (tomato2) over 160 regions of 200 kb (32 Mb),
at 5 to 12 reads a position, 6 threads, with fitted parameters. Both features the owner asked for
are on by default and were left at their defaults: the hidden-duplication filter (false-discovery
rate 0.01) and the neighbour-evidence pass (`--refine-genotypes auto`, which runs at this cohort's
median inbreeding of 0.99). Measured on main `e6565a37`, changes on branch `perf-psp-vcf`.
`estimate-parameters`, the fit that runs before this command, is out of scope here: it was
reviewed on 2026-11-04 and left for the owner.

**Targets** (owner, 2026-11-05): wall time and peak memory. A large gain in one at a small cost in
the other is a win; a change that makes the code clearer is kept even for a small gain, and one
that complicates it for a small gain is rejected.

**Why no sub-agent fan-out.** The laptop has six fast cores, other sessions kept its load between
2 and 12, and the owner's rule is that performance measurements share at most six threads. One
reviewer measuring in turn, every comparison alternating the two builds and started only when the
one-minute load was below 3 (`tmp/perf/when_quiet.sh`), gave comparable numbers.

**Where the time and memory went on main** (the run's own timing report, three runs):

| stage | wall | cores busy | what it does |
|---|---|---|---|
| calling | 55 s | about 4 of 6 | read every sample's psp forward, build each locus across the cohort, genotype it, park the record for the filter |
| scoring | 13.6 s | 4.6 | the filter's likelihood ratio, duplicated against single-copy, at every record |
| writing | 53 s | about 5 | the neighbour-evidence pass, then the VCF |
| **run** | **122 to 124 s** | | peak resident memory 1,130 to 1,146 MB, reached during the pass |

Resident memory by stage: about 400 MB through calling, 600 MB while scoring, 1.1 GB in the pass.

## 2. Verdict

**Apply the listed wins** — three commits, each measured alone against the one before it, each
leaving the VCF body bit-identical:

| commit | change | the run, before → after |
|---|---|---|
| `824ec349` | the pass's walk keeps its window-long arrays and reuses them | 123.8 → 115.4, 116.0 s |
| `e933395c` | a round's calls are handed to the writer while the pool works the next round | 116.5, 114.2 → 105.4, 106.0 s |
| `55f01e9b` | the filter's scoring takes one logarithm per group of samples | 107.3, 108.3 → 99.5, 98.6 s |

**Together: the run from about 123 s to about 99 s, a fifth faster, with peak memory level**
(1,123 to 1,138 MB against 1,130 to 1,146 MB). Instructions retired fell from 11,052 G to 9,092 G.

## 3. Measurement plan

What was run, so the numbers can be repeated (`tmp/perf/` in the `perf-psp-vcf` worktree; the
psps, catalog and parameters are in the main checkout's `tmp/parallel/`, regenerated with main on
2026-11-05):

- `ab.sh <binary> <label>` — the command under `/usr/bin/time -l`: wall, peak, instructions, the
  calling and after-calling times from the log, and a digest of the VCF body for bit identity.
- `when_quiet.sh` — `ab.sh` once the one-minute load is below 3.
- `rss.sh` — `ab.sh` with resident memory sampled every half second, so a peak can be put in a stage.
- `prof.py` — self and inclusive CPU by function from a `samply --unstable-presymbolicate`
  profile, over a time window, and `--callers` for who calls a function.
- `lines.py` — self CPU by source line, through inlining, with `atos -i`.
- `alloc_lines.py` — allocator CPU by the line of ours that allocated.

## 4. Build / toolchain configuration

No change. One experiment: mimalloc told never to return memory to the system
(`MIMALLOC_PURGE_DELAY=-1`) gave 121.1 s against 121.6 s and the same instructions, so the
allocator's cost in the pass was its allocations, not its page returns.

## 5. Code-level findings

### Hot-path

**H1: [model.rs](../../../../src/neighbour_evidence/model.rs) — the walk allocated about 25 vectors a site** — applied in `824ec349`.

The pass walks each sample along each window five times; each walk computed every site's
emissions and forward probabilities into fresh vectors (289 pair states a site), and the backward
step and posterior allocated about ten more a site. The profile put 43 CPU-seconds of the pass's
365 in the allocator, 37 of them called from the walk, and the single largest allocation was each
site's 289 pair emissions (9.9 CPU-seconds). Now the emissions and forward probabilities live in
flat arrays a window long (`Emissions`, `Forward`), and the backward step writes over buffers it
keeps for the walk.

**The arrays are kept per thread and reused walk after walk**, and that is what keeps the memory
level. Allocated afresh each walk, the same arrays held 1,287 to 1,365 MB at the pass's peak
against 1,130 to 1,146 MB on main: a window's arrays are megabytes, and the allocator does not hand
a freed large block straight to the next walk the way it reused the small per-site ones.

| | wall | writing stage | instructions | peak |
|---|---|---|---|---|
| main | 123.8 s | 53.7 s | 11,051 G | 1,131 MB |
| flat arrays, reused | 115.4, 116.0 s | 45.0, 45.6 s | 10,517 G | 1,099, 1,121 MB |

Bit-identical: the posterior and phase hash tests of `refine/tests.rs` are unchanged.

**H2: [rounds.rs](../../../../src/run/cohort_merge/rounds.rs) — the pool sat idle while the calling thread handed over each round** — applied in `e933395c`.

The calling stage runs in rounds: for 59 samples a round covers about 8,500 bases in 6 building
regions, and the 32 Mb take about 3,800 rounds of 15 ms. Timers in the loop split the stage's 55 s
into reading every sample's psp forward over the round (25.8 s), building and genotyping the
regions (19.1 s), and handing the calls over in genome order (9.9 s). The hand-over runs on the
calling thread alone — it applies the frontier rule and gives each record to the sink, which on
this path formats it for the filter's spill file — and the pool waited through it.

Now each round's eviction, reading and building is spawned on the pool while the calling thread
hands over the round before. Rounds are still handed over one after another and each in region
order, so the records reach the sink in genome order, and a round that fails is reported after the
one before it was handed over, as before. With one thread the old order stands.

| | calling | the run |
|---|---|---|
| before | 55.6, 54.3 s | 116.5, 114.2 s |
| after | 45.6, 45.7 s | 105.4, 106.0 s |

Peak 1,079 to 1,151 MB against 1,088 to 1,120 MB. **One cost:** the calling thread now works beside
the pool's `--threads` workers, so the run can hold one core more than `--threads` for the length of
a hand-over.

**H3: [locus_score.rs](../../../../src/paralog/locus_score.rs) — the filter's scoring took an `exp` and a `log1p` per sample per grid point** — applied in `55f01e9b`.

At each record the filter compares two explanations of every sample's depth and allele counts: a
real single-copy variant, averaged over 200 allele frequencies and three genotypes, and a hidden
duplication, averaged over 40 carrier frequencies and 7 copy configurations. Each was a sum of
logarithms a sample, each logarithm a log-sum-exp: about 680 `exp` calls a sample a record, 12
billion over this run, which was 60 of the stage's 64 CPU-seconds.

Now each sample's terms are scaled by their largest, multiplied as probabilities over a group of
sixteen samples, and one logarithm taken per group; the priors' probabilities are precomputed once
a run. **A group whose product is not a normal number is summed a sample at a time in log space,
as before.** That happens only where a prior is floored to 1e-300 — an inbreeding coefficient of
exactly 1, or a carrier frequency of zero — because the scaling keeps every term between the
probability of the genotype or branch the reads favour and one, whatever the depth.

| | scoring | the run | instructions |
|---|---|---|---|
| before | 13.7, 13.9 s | 107.3, 108.3 s | 10,510 G |
| after | 3.0, 3.0 s | 99.5, 98.6 s | 9,092 G |

The VCF body is bit-identical and the same 17,504 records are dropped. The new test
`the_products_over_groups_give_the_log_space_sums` holds both explanations to the old log-space
sums within 1e-9 nats on an ordinary locus and on two built to reach each fallback (checked, with a
temporary panic, to reach them; no existing test did). The production differential's widest gap
went from 3.6e-12 to 7.3e-12 nats, under its 1e-10 pin. The module header records this as the
file's fourth departure from production's arithmetic, as its rules ask.

### Likely

None open.

### Note

- **The live-read sets rebuilt per locus**
  ([held_live_sets.rs](../../../../src/run/psp_source/held_live_sets.rs)): 12.6 CPU-seconds of the
  calling stage copy the stored set at the nearest restart point. Counted: about 34 million
  rebuilds, each of about 8 entries and 7 read identifiers, so the cost is about half a
  microsecond a rebuild, most likely cache misses on the stored sets. A per-thread, per-sample
  cache of the last rebuilt set could save at most about 4 s of wall time and would add state to
  code tuned against the 2,169-sample memory failure. Not pursued.
- **The drawn haplotypes** ([model.rs](../../../../src/neighbour_evidence/model.rs)
  `haplotypes_along`): about 10 CPU-seconds of allocation, a joint and two vectors a site a path,
  ten paths a walk. Removing it changes the walk's output type, which `phase.rs`, `refine.rs` and
  about fifteen tests read, for about 2 s of wall time. Not pursued.
- **More building regions per round than threads** (three times as many, the same ground):
  93.3 and 99.3 s against 96.6 s on the final commit. No gain shown — rejected.
- **The calling stage's reading** is now its largest part: psp decompression, record decoding and
  the coverage accumulator, spread over many functions, with the fixpoint over samples adding a
  barrier per sweep.

## 6. Out-of-scope observations

- **The fit** (`estimate-parameters`, 6 min 25 s on these 59 samples) remains three times this
  command; see `perf_pipeline_2026-11-04.md` §3.

## 7. What's already good

- **The walk's arithmetic was already written in a fixed order over slices**
  ([model.rs](../../../../src/neighbour_evidence/model.rs)), so moving its storage changed no bit,
  and the hash tests confirmed it on the first run.
- **The run times its own stages** (the hidden-duplication filter's report: calling, scoring,
  writing), which made every comparison here a per-stage one without instrumentation.

### Author response convention

Address each finding by its identifier with one of: `applied in <commit>` / `experiment shows no
gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`.
