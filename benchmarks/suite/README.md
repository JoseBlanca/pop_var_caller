# The benchmark suite

`benchmarks/suite/run_suite.sh` reruns every benchmark whose figures are quoted about pvc —
accuracy, speed and memory — and regenerates the results table and its summary. Every figure in the
README and the user guide comes from its output. The plan behind it is
[`doc/devel/implementation_plans/benchmark_suite.md`](../../doc/devel/implementation_plans/benchmark_suite.md).

## What it produces

For each pvc commit, under `benchmark_results/pvc/<commit>/` at the repository root (ignored by git):

- **`results.tsv`** — one row per figure. Columns: benchmark, dataset, samples, depth, region size
  (bp), variant class, caller, caller version, pvc commit, stage, metric, value, date, machine, the
  host's mean load while the figure's step ran, and the step that produced it (`source`).
- **`summary.md`** — the same figures as tables, one section per benchmark, each saying which data
  it was measured on.
- **The step directories** — every step's commands (`logs/<name>.cmd`), their output and error
  streams, and a `.done` file holding the step's measurement and rows. The `source` column of a row
  names the directory, so any figure leads back to the exact command that made it.

`benchmark_results/latest` points at the directory of the commit run last.

## Running it

```
benchmarks/suite/run_suite.sh                         # every benchmark, skipping what is done
benchmarks/suite/run_suite.sh --only giab_panel       # one benchmark (comma-separated for several)
benchmarks/suite/run_suite.sh --only hg002_tr --depth 30x --caller pvc
benchmarks/suite/run_suite.sh --status                # steps done, and the next one, per benchmark
benchmarks/suite/run_suite.sh --summary-only          # rebuild results.tsv and summary.md
benchmarks/suite/run_suite.sh --slice --only speed    # the small test configuration
```

It runs on the host and starts the dev container itself, with 6 cores (the laptop's fast ones;
`--threads N` to change) and 32 GB (`BENCH_MEM=48g` to change). It refuses a checkout with
uncommitted changes unless given `--allow-dirty`, which records the commit as `<sha>-dirty`.

From a worktree, the data are read from the main checkout (found automatically, or
`--data-root DIR`), which is mounted read-only; results are written in the worktree.

### Running it over several sessions on a machine used during the day

```
benchmarks/suite/run_suite.sh --stop-at 08:00     # start in the evening; nothing new starts after 08:00
benchmarks/suite/run_suite.sh --stop              # from another terminal: finish this step, then stop
benchmarks/suite/run_suite.sh --redo-busy         # re-time the steps that ran on a busy machine
```

- **It waits for a quiet machine.** Before every caller step whose time is recorded, the suite
  checks that processes other than itself used under two cores in every 5-s sample of the last
  2 minutes, and otherwise waits, checking every 30 s (`--quiet-cores`, `--quiet-minutes` change
  the test; `--no-wait` starts regardless, for accuracy-only sessions — such times are marked ⚠).
  Inputs and scoring do not wait. Two cores, not one: an idle Mac's background services (the
  window server, the editor, a PDF preview service) use about one core between them, on its slower
  cores, without touching the six fast ones the suite runs on.
- **It keeps the Mac awake** while it runs (`caffeinate`).
- **It stops when told or when the time comes**: `--stop`, or `--stop-at HH:MM` (local time; the
  next time the clock reads it), both between steps; Ctrl-C at once.
- **Nothing finished is lost.** Rerun the same command next session; it starts where it stopped.
  A step a machine got busy in the middle of is kept for its results and marked; `--redo-busy`
  re-times it later.

### Stopping and resuming

A run is a list of small steps — one caller's one stage on one configuration, or one scoring pass.

- `run_suite.sh --stop`, from another terminal, ends the run **after** the step that is running.
- Ctrl-C ends it **at once**; the step that was running is discarded.
- Running the same command again resumes: every finished step is skipped.

The longest single steps are expected to be about an hour (GATK at 300× over the tandem-repeat
regions; freebayes on the simulated 50-sample runs at 60×).

### What a new pvc commit reruns

| directory | holds | rerun |
|---|---|---|
| `benchmark_results/inputs/` | inputs the suite makes: downsampled alignments, the simulated genomes and their alignments, region lists, the tandem-repeat tract list | never, unless deleted |
| `benchmark_results/competitors/` | GATK's, freebayes' and HipSTR's calls, time, memory and scores | only with `--rerun-competitors` |
| `benchmark_results/pvc/<commit>/` | pvc's calls, time, memory and scores, and the table and summary | for every new commit |

A competitor's stored calls carry a fingerprint of its version, its command and its inputs. If an
input has changed, the suite stops and names it rather than quietly rerunning hours of work; rerun
with `--rerun-competitors` (narrowed with `--only` and `--caller`) to remake it. `--redo-step
<source>` redoes any one step, for instance a speed figure measured while the laptop was busy (the
summary marks those with ⚠).

## The benchmarks

| name | what | callers |
|---|---|---|
| `giab_panel` | GIAB HG002, HG003, HG004, each alone over its own 100 confident regions, 5× to 300× | pvc (defaults, fitted), GATK, freebayes |
| `giab_trio` | the trio jointly over 100 shared regions, 5×, 10×, 30×, 300×; per-member accuracy and Mendelian errors | the same |
| `hg002_5mb` | HG002 over 1,000 confident intervals (5.08 Mb), 30× and full depth; also the one-human-sample speed figures | the same |
| `hg002_tr` | GIAB HG002 tandem-repeat benchmark, 50,000 regions, 5× to 300×; STR genotypes per tract | pvc (defaults, fitted), GATK, HipSTR |
| `stutter_sim` | simulated 4 Mb genome with 3,000 STR segments and 500 bundles, PCR and PCR-free, 1 and 50 samples, 5× to 60× | pvc (defaults, fitted, true stutter), GATK, freebayes, HipSTR (twice) |
| `speed` | tomato cohort, 1, 8 and 62 accessions over 8 Mb (a 63rd, SRR7279543, has a truncated CRAM and is left out): wall time, CPU time, peak memory of each whole pipeline | pvc (both ways), GATK, freebayes |
| read floor (in `speed` and `hg002_5mb`) | `samtools view -c -M -L` over the same files: the cost of only reading and decoding them, set beside pvc's `generate-psps` and `call-from-alignments` with their ratio | samtools, run once like a competitor |
| `tomato_agreement` | pvc's 62-accession genotypes against GATK's (no truth set) | pvc against GATK |
| `external_chr22` | the colleague's chr22 benchmark; not in a full run, needs `--external-config` (see `bench/external_chr22.py`) | pvc (defaults, fitted), GATK, freebayes |

How each is scored:

- **SNPs and indels:** RTG vcfeval (fetched, pinned, by `fetch_tools.sh`), twice per call set — the
  genotype must match, or only the alleles — split into SNPs and indels, plus each caller's best F1
  over a sweep of its QUAL threshold. Records count as each caller writes them, except freebayes,
  kept at QUAL ≥ 20 as its documentation advises.
- **The three error types:** every call set vcfeval scores is scored again by pvc_eval
  (github.com/JoseBlanca/pvc_eval, `lib/pvc_eval.py`), which sorts each error into *wrong GT (truth
  ref)*, *wrong GT (truth var)* or *no call (./.)*, per SNP, indel, STR and complex locus. vcfeval
  cannot count genotypes left missing; vcfeval stays as the check. pvc_eval runs from a checkout:
  a profile's `BENCH_PVC_EVAL`, or else `tmp/pvc_eval` in the data checkout or `../pvc_eval` beside
  it; its commit is recorded, and its interpreter (free-threaded Python 3.15 with numpy) is in the
  dev image. With no checkout the suite scores with vcfeval alone.
- **STRs:** `benchmarks/lib/tract_qual_experiment.py` on HG002 and `scripts/ng_str_stutter_score.py`
  on the simulation, which compare the rebuilt DNA over each repeat.
- **Speed and the machine's state:** every timed step first reads its input files once, untimed,
  so every tool reads from memory. The host records its load every 5 s; each time figure carries
  the load average at its step's start and end and the cores used meanwhile by processes other than
  the suite (`other_cpu_cores`). A figure whose step ran with other processes above two cores is
  marked ⚠ in the summary.
- **Parallelism:** each caller as its documentation tells a user — one process per sample for the
  one-sample stages (pvc `generate-psps`, GATK HaplotypeCaller GVCFs), one threaded process for the
  cohort stages, freebayes as one process per block of regions. GATK's cohort stages are
  single-threaded, so they run as GATK's joint-calling workflow runs them: one process per shard of
  the regions, as many at once as the suite has cores. They join the GVCFs with GenomicsDBImport,
  GATK's recommended route, on x86-64, and with CombineGVCFs on the arm64 laptop, where
  GenomicsDBImport's native library does not exist; the stage rows name the route.
- **GVCFs made once:** where a benchmark calls nested cohorts (the speed benchmark's first 1, 2,
  5, ... accessions), each sample's GATK GVCF is made once, in a run of the largest cohort, and
  every cohort's GATK step joins its share. That cohort's per-sample stage figures are derived from
  each sample's own measured time and memory, laid out as many at once as the suite has cores, and
  are marked † in the summary.

## DeepVariant, and running on an x86-64 Linux server

DeepVariant (1.10.0) and GLnexus (1.4.1), its joint genotyper for the trio, run in their official
images. The dev container cannot start containers, so `run_suite.sh` starts a job runner on the
host (`host_jobs.sh`) that runs only those two images, through the host's podman or Docker, with the
project mounted and the data read-only, measured by the suite's own code inside them.

- **On this Mac** (arm64) the runner says it cannot, and every DeepVariant step, with the scoring of
  its calls, is listed as skipped at the end of the run.
- **On an x86-64 Linux server with podman**, the same command runs everything, DeepVariant
  included. First time there:

  ```
  podman pull docker.io/google/deepvariant:1.10.0 ghcr.io/dnanexus-rnd/glnexus:v1.4.1
  benchmarks/suite/run_suite.sh --slice --data-root <checkout with the data> --threads <cores>
  ```

  The slice run is the first test of the DeepVariant steps: they were built on the Mac, where they
  cannot run.
- **No data needs copying.** Where the data checkout lacks a human benchmark's files, the suite
  fetches them from GIAB and NCBI as input steps (`lib/giab_sources.py`); the project's region lists
  and downsampling fractions travel in `benchmarks/suite/data/`. The tomato CRAMs and reference come
  from the machine's profile.
- **Machine profiles** (`profiles/<hostname>.sh`, or `--profile NAME`) set what differs per machine:
  `profiles/kimura.sh` gives the suite 32 cores and 80 GB, the tomato collection under
  `/home/joxi/tomato25`, and mounts it. On kimura the tomato cohort is tomato1's 62 accessions found
  by sample name in the collection, then further accessions for larger sample counts
  (`BENCH_TOMATO_SIZES` or `--tomato-sizes 1,2,5,10,20,50,100`).
- **On Linux** the quiet check counts every process except this user's podman containers' (the
  dev container and the DeepVariant jobs), and the profile's core and memory limits are applied to
  podman.
- Figures from the server are the server's: every caller, the competitors included, is rerun there,
  so no figure mixes two machines.

## How long, and how much disk

Estimated for the laptop at 6 cores; the plan's §8 says how each figure was reached.

| | time | disk |
|---|---:|---:|
| first run (inputs, competitors and pvc) | about 19 h (14 to 30) | about 50 GB |
| a refresh for a new pvc commit | about 5½ h (3 to 8) | about 5 GB more |

Most of the first run is the simulated grid (about 12 h). A pvc commit's directory can be deleted
once its figures are no longer wanted.

## Code

- `run_suite.sh` — host side: options, the machine's description, the host-load sampler, the
  container.
- `suite.py` — inside the container: builds pvc and the simulator, runs the steps, writes the
  table and summary.
- `bench/` — one module per benchmark, each a list of steps.
- `lib/` — the step runner (`steps.py`), measurement (`measure.py`), how each caller is run
  (`callers.py`, `arms.py`), inputs (`inputs.py`), scorers (`vcfeval.py`, `mendelian.py`,
  `agreement.py`, `stutter_score_callers.py`, `hipstr_qual.py`), and the table (`results.py`).
