# One script that reruns every benchmark — implementation plan

*Status: plan, 2026-10-18, revised the same day after the owner's review (§2). Branch
`benchmark-suite`, from main at `bbcc9f46`, worktree `../pop_var_caller-benchmark-suite`.*

## 1. What this builds

One command, `benchmarks/suite/run_suite.sh`, that reruns every benchmark we quote about pvc and the
competitors, scores every caller the same way, and writes three things:

- **`results.tsv`**: one row per figure. The columns are benchmark, dataset, samples, depth, region
  size, variant class, caller, caller version, pvc commit, metric, value, date and machine.
- **`summary.md`**: generated from that table, for the README and the user guide to quote. Every
  figure in it names the data it was measured on.
- **`logs/`**: the exact command, its standard error, and its time and memory for every step, so any
  row of the table can be traced to the command that produced it.

The script can run one benchmark, one depth or one caller alone. **The competitors (GATK, freebayes,
HipSTR) are run once.** After that, a refresh reruns only pvc, scores it, and puts it beside the
competitors' stored calls. A competitor is rerun only when asked for by name (`--rerun-competitors`).
It is also rerun when its stored calls are missing, or when its input files have changed.

**A run can be stopped at any point and resumed later without redoing finished work** (§6.1).
Every step that finishes leaves a marker, and a resumed run skips every step that has one.

It starts from the inventory of 2026-10-18
(`pop_var_caller-docs/tmp/docs_inventory/benchmark_inventory.md`). That inventory and three surveys
of the scripts and data behind it are the source of every "exists" and "stale" statement below.

## 2. Settled with the owner

The owner's review of the first draft settled four things:

- **Competitors are run once; a refresh reruns only pvc.**
- **Every run uses the laptop's 6 fast cores only**, not all 18. This replaces the first draft's
  proposal to measure speed at 18 threads.
- **A run must stop and resume without redoing finished work**, because a full run takes about a day
  and the laptop is needed for other work.
- **The 2,169-sample cohort on `kimura` is left out** of this plan, to be discussed separately.

The two proposals below received no objection and are built as written.

**D1. Run every caller inside the dev container, with its virtual machine given the 6 cores the suite
uses and 32 GB.**

- **Why not the host.** GATK and the scorer, vcfeval, are Java programs. This Mac has no Java
  runtime, so neither can run on the host without installing one. Every other tool could run on
  either side.
- **What one environment buys.** All speed and memory figures come from the same Linux machine, with
  the same compiler for pvc and the same tool builds.
- **What it costs.** The figures are measured inside Apple's Linux virtual machine, not natively on
  macOS. The summary will say so beside every speed figure.
- **How big the virtual machine is.** The suite sets its size with the existing `DEV_CPUS` and
  `DEV_MEM` settings: 6 cores, and 32 GB. The 32 GB is because the 63-sample freebayes run was killed
  for memory at 16 GB before (`run_freebayes.sh` header, lines 26–40). It leaves 32 GB for the rest
  of the laptop.
- **The other answer.** Installing a Java runtime on the host would let everything run natively. It
  writes outside the project, so I would not do it without your say.

**D2. Score each caller's records as it writes them: records that pass its own filters, with no
extra quality threshold. The one exception is freebayes, kept at QUAL 20 or more.** Each caller also
gets one more row, its best F1 over a sweep of the quality threshold.

- **Why freebayes is the exception.** freebayes writes every candidate it considers. Its
  documentation recommends a quality filter, and the earlier STR comparison used QUAL 20.
- **pvc.** pvc already applies its own threshold (`--min-site-quality`, default 1) and writes no
  record below it. The old benchmark scripts added a second gate at QUAL 30 on top. The suite drops
  that gate, because a user does not apply it.
- **The best-F1 row.** It answers "how good is each caller when its threshold is tuned", which is
  what the colleague's chr22 report quoted. It is computed from vcfeval's own output, so it adds no
  caller run.
- **The other answer.** Gating every caller at one shared threshold, such as QUAL 30, would treat
  the callers more uniformly. But QUAL means a different thing in each caller, so the same number is
  not the same strictness.

**Six cores, and what the suite cannot guarantee about them.** This laptop has 6 fast cores and 12
slower ones. macOS gives no way to pin a process to the fast cores. The suite gives the virtual
machine 6 cores and runs every caller at 6 threads or 6 processes. macOS normally schedules a busy
virtual machine's 6 threads on the fast cores, but this is not guaranteed. To show whether it held,
the suite records the host's load before and after each benchmark. The summary states every speed
figure as "6 threads". `--threads` can change the count for a later run.

## 3. How things run and how they are measured

### 3.1 Where

- **The driver.** It runs on the host and starts the container once per benchmark through
  `scripts/dev.sh`.
- **Data.** Benchmark data stay in the main checkout and are read-only to the suite:
  - the main checkout is passed in with `DEV_EXTRA_MOUNT`;
  - the human and tomato references come from `~/genomes`, which `dev.sh` already mounts read-only.
- **Outputs.** Everything is written under the worktree.
- **The pvc binary.** It is built in the container (`cargo build --release`) from the checked-out
  commit. If the tree has uncommitted changes, the suite refuses to run unless given `--allow-dirty`,
  and then it records the commit as `<sha>-dirty`.
- **Before each benchmark**, the host side records:
  - the machine: model, 18 cores (6 fast, 12 slower) and 64 GB;
  - the container's core and memory allowance;
  - the host load average, taken from the host and not from the virtual machine.

  The load average is taken again after each benchmark. A run that starts with a host load above 3
  is marked in the table, not refused.

### 3.2 Tools and versions

Every tool's version is read from the tool itself at run time and written into the table.

| tool | version | where | notes |
|---|---|---|---|
| pvc | the checked-out commit | built in the container | |
| GATK | 4.6.0.0 | the dev image | |
| freebayes | 1.3.6 (Debian's package) | the dev image | The host has 1.3.10, but it cannot run inside the Linux container. The earlier freebayes calls on disk came from the host's 1.3.10, so the suite reruns freebayes rather than reusing them. |
| HipSTR | 0.6.2 | the dev image | |
| minimap2, samtools, bcftools 1.16 | | the dev image | bcftools carries its `mendelian` plugin. |
| RTG Tools (vcfeval) | 3.13 | a pinned download | Not in the image. A `fetch_tools.sh` step downloads the pinned release into the repository's `tmp/tools/`, checks its checksum, and runs it with `java -jar`. The `rtg` launcher refuses the container's arm64 processor, as `scripts/ng_str_stutter_score.py` already notes. I chose a download over adding it to the image, because rebuilding the shared image takes about ten minutes and affects other sessions. |

### 3.3 Time and memory

Each step runs under a small Python wrapper, `benchmarks/suite/lib/measure.py`. It uses
`python3-psutil`, which is already in the image, because the image has no GNU `time`.

| case | what is measured |
|---|---|
| **One process** | wall time; user plus system CPU time; peak resident memory. These come exactly from the operating system when the process ends (`wait4`). |
| **Several processes at once** (§3.4) | wall time from the first start to the last finish; CPU time summed; and two memory figures: the largest single process, and the **summed resident memory of all processes alive at once**, sampled every 0.2 s. The summed figure is the real cost of the stage. |

Each measurement is written to the step's log as one JSON line. The table takes its time and memory
rows from there.

### 3.4 Parallelism: each caller as its documentation tells a user to run it

This follows the rule settled in June for these comparisons:

- **One-sample tools run one process per sample**, up to 6 at once, each single-threaded:
  - GATK HaplotypeCaller, making one GVCF per sample;
  - pvc `generate-psps`, whose help says "to spread a cohort over cores, run one command per sample".
- **Cohort tools run as one multi-threaded process:**
  - pvc `estimate-parameters`, `call-from-psps` and `call-from-alignments`, at `--threads 6`;
  - GATK GenomicsDBImport, then GenotypeGVCFs.
- **freebayes has no threads**, so it runs as its own `freebayes-parallel` does: one process per
  region, up to 6 at once.
- **A one-sample run** has nothing to spread over samples:
  - GATK HaplotypeCaller runs as one process with `--native-pair-hmm-threads 6`;
  - freebayes still runs by region.

For cohort joint calling I chose GenomicsDBImport over the CombineGVCFs used before. GATK's Best
Practices name it for cohorts, so it is what a GATK user with 63 samples would run.

**A memory risk for the tomato freebayes run.** At 50 samples, four freebayes processes reached
9.6 GB together, so 6 processes at 63 samples may not fit in 32 GB. The tomato slice test measures one region's
process at 63 samples. If 6 at once would not fit, the suite lowers freebayes' process count and
says so in the table.

## 4. The benchmarks

**Times.** Every time below is an estimate, not a measurement. Each is built from the recorded times
the surveys found, or from the tools' size and depth where none was recorded. The slice tests (§7)
replace these estimates with measured rates before the full run.

### 4.1 GIAB per-sample panel: SNPs and indels, one sample at a time

**What it is.** GIAB HG002, HG003 and HG004 on GRCh38. Each sample has its own 100 regions, between
401 kb and 572 kb a sample. There are seven depths: 5, 10, 15, 20, 30 and 50×, which are
downsampled, and the full reads, called "300×". The truth is GIAB v4.2.1 inside each sample's
confident regions.

**What exists and works.**
- All 21 alignment files, the truth files and the region files are on disk.
- No script in the repository can recreate the downsampled files. Their `.coverage.json` files
  record that they were made with `samtools view -s 42.<fraction>` from a source not in the
  repository. The suite uses them as they are.
- `benchmarks/giab/src/run_ng_per_sample.sh` still calls the current pvc correctly.
- `run_freebayes_per_sample.sh` works, but uses the host's freebayes.

**What is stale or missing.**
- GATK has never been run on this panel, and there is no runner for it.
- Scoring (`score_ng_recall.sh`) matches alleles by position and spelling and ignores genotypes. It
  is not vcfeval.
- `benchmarks/giab/README.md` line 29 names `run_ours_per_sample.sh`, the deleted caller's runner.
- The panel's own copy of the reference has no sequence dictionary, which GATK needs. The suite
  instead uses the identical reference in `~/genomes/h_sapiens/gca_grch38/`, which has one. The
  suite first checks that its index file matches the panel's copy.

**What the suite runs, for each sample and depth:**
- pvc with `--defaults`, through `call-from-alignments`;
- pvc fitted: `generate-psps`, then `estimate-parameters` on that one sample's psp, then
  `call-from-psps --parameters`;
- GATK HaplotypeCaller;
- freebayes.

**Scoring.** vcfeval (§5.1), against the truth, inside the sample's confident regions.

**Estimate.** About 1½ hours, and under 3 GB on disk.

### 4.2 GIAB trio: cohort calling and Mendelian consistency

**What it is.** The trio over 100 regions the three share, 452 kb in total.

**What exists.**
- Full-depth alignments for all three, the truth files and the region file are on disk.
- No runner, no scorer and no Mendelian check exist. All three are written new.
- There are no downsampled depths. The suite makes them at 5, 10 and 30× the same way the panel's
  were made (`samtools view -s 42.<fraction>`). It computes the fraction from each file's measured
  mean depth over the regions. It keeps them under the results directory, so later runs reuse them.

**What the suite runs, at 5, 10, 30× and full depth:**
- pvc on all three samples jointly, with `--defaults` and fitted;
- GATK: three GVCFs, then GenomicsDBImport, then GenotypeGVCFs;
- freebayes on all three samples together.

**Scoring.**
- vcfeval for each member against that member's truth.
- The Mendelian error rate of each caller's joint calls, from bcftools' `mendelian` plugin: of the
  sites where all three members have a called genotype and at least one is not reference, how many
  violate inheritance.
- The same rate for the truth itself, as the floor.

**Estimate.** About 45 minutes, and under 2 GB.

### 4.3 HG002 over 5 Mb: one human sample, accuracy and speed

**What it is.** One HG002 CRAM over 1,000 confident intervals, 5.08 Mb. The truth is GIAB v4.2.1
sliced to them.

**What exists.**
- `benchmarks/human_genome_bottle/` holds the data.
- `benchmarks/lib/run_ng.sh`, `run_gatk.sh` and `run_freebayes.sh` all still work.
- GATK and freebayes calls from June are on disk, but they were made with `-stand-call-conf 0` and
  freebayes 1.3.10, so they are not reused (D2, §3.2).

**Depths.** Full depth (about 300×), plus 30× made by the suite as in §4.2. The suite runs the same
four arms as §4.1.

**Speed.** This is also the "one human sample" speed benchmark (§4.6). Each caller's steps are timed
there.

**Estimate.** About 45 minutes, of which GATK at full depth is about 10–20 minutes. Under 2 GB.

### 4.4 HG002 tandem-repeat benchmark: STR accuracy against HipSTR and GATK

**What it is.** GIAB's HG002 tandem-repeat benchmark v1.0.1, 50,000 regions over 6.09 Mb, at seven
depths from 5× to 300×. The alignment files are on disk.

**What exists and works.**
- `benchmarks/ssr_hg002/src/run_ng_coverages.sh` runs the current pvc.
- `benchmarks/lib/tract_qual_experiment.py` scores any VCF that has QUAL and GT, GATK's included.

**What is stale or missing.**
- `benchmarks/ssr_hg002/README.txt` describes the deleted caller's commands (lines 8–9, 24, 39–46,
  83–86, 120–124) and its headline numbers (lines 138–157).
- The HipSTR calls on disk were made in July over 13,272 regions from the deleted caller's
  catalogue. The suite reruns HipSTR over the benchmark's own regions.
- GATK has never been run over these regions.
- No script fits parameters on this sample.

**What the suite runs, at each depth:**
- pvc with `--defaults`, and fitted;
- HipSTR at its defaults, with `--use-unpaired`, as before;
- GATK HaplotypeCaller.

**Scoring.**
- `tract_qual_experiment.py`: right genotypes per tract, split into homopolymers and repeats of
  period 2 or more.
- vcfeval for the non-repeat SNPs and indels in the same regions.

**Estimate.** About 3 hours. GATK at 300× over 6 Mb is most of it. Under 3 GB.

### 4.5 The simulated stuttering benchmark: STR accuracy with exact truth

**What it is.** A drawn 4 Mb genome with 3,000 STR segments and 500 STR bundles. Reads are drawn with
PCR or PCR-free stutter and mapped with minimap2. The grid has 16 runs: PCR and PCR-free, 1 and 50
samples, and 5, 15, 30 and 60×.

**What exists.**
- The driver `scripts/ng_str_stutter_benchmark.sh` does almost everything:
  - simulate the reads, map them, call with pvc (with `--defaults`, with the true stutter, and
    fitted), call with GATK and with HipSTR, and score;
  - it can resume a run that stopped part-way;
  - it records each stage's time.
- `scripts/ng_str_stutter_score.py` scores STR segments and bundles, and runs vcfeval on the rest.

**What is stale or missing.**
- **No simulated data survive.** They lived in worktrees that are gone, so every run is regenerated
  from its seed.
- **freebayes is missing from the driver.** It was run by hand for the October report. The suite
  adds it: all samples together, QUAL 20 or more.
- **The STR accuracy scripts point at a deleted worktree.** `scripts/ng_str_accuracy/copy_run.sh:18`
  points at `pop_var_caller-tract-ceiling`. They are not needed once the driver has freebayes.
- **Fitted stutter is now the headline arm.** The driver already runs `estimate-parameters` on every
  run. The grid script `ng_str_stutter_first_pass.sh` reused the 30× fit at the other depths; the
  suite fits every run on its own data, as a user would.
- **The true-stutter arm is kept** as a reference row, labelled as such.
- **GATK parallelism.** The driver runs GATK as four processes. It changes to the rule in §3.4.

**Scoring.**
- The driver's scorer, for pairs called right and false calls on reference pairs, with STR segments
  and STR bundles counted apart. A pair is one sample at one STR segment.
- HipSTR is scored on segments only, because it is not given bundles.

**Estimate.** About 13 hours, the bulk of the suite:
- each 50-sample run is about 1–2 hours of simulating, mapping, four callers and the fit;
- the 1-sample runs total about 1½ hours.

**Disk.** About 35 GB if the read alignments are kept (0.9 GB at 50 samples and 5× up to 8.9 GB at
60×). About 5 GB if they are deleted after scoring. The suite keeps them, because a later pvc
refresh calls the same alignments the competitors called. Regenerating them would cost hours, and the
competitors' fingerprints would then have to be checked against the regenerated files.

### 4.6 Speed and peak memory, whole pipelines

- **Tomato cohort.** `benchmarks/tomato1/`: 63 accessions at about 3 reads a position, over 80
  regions of 100 kb (8 Mb of SL4.0). It runs at 1, 8 and 63 samples, so the table shows how each
  caller scales.
- **One human sample.** The two HG002 runs of §4.3, at 30× and full depth.

**Each pipeline is timed step by step and in total:**

| caller | steps timed |
|---|---|
| pvc, in steps | `repeat-catalog` once per reference, `generate-psps`, `estimate-parameters`, `call-from-psps` |
| pvc, in one step | `call-from-alignments` |
| GATK | HaplotypeCaller GVCFs, GenomicsDBImport, GenotypeGVCFs |
| freebayes | one step |

**What exists.**
- `benchmarks/tomato1/scripts/run_perf_all.sh` and its `perf_*.py` measured freebayes and GATK in
  June at 4 threads.
- The pvc half was the deleted caller's (`perf_ours_*`), since removed.
- The suite takes the measuring idea, but not the scripts, because their sizes stop at 50 samples
  and they call deleted commands (`perf_dashboard.py`, `build_cohort_intermediates.sh`).

**Tomato's stale files.**
- The header of `benchmarks/tomato1/bench.config.sh` describes 26 samples over 2 Mb.
- The GATK joint VCF on disk is from a 26-sample cohort that is no longer on disk.
- Only 50 of the 63 GVCFs exist.

All are regenerated.

**Estimate.** About 2½ hours for tomato at 1, 8 and 63 samples. freebayes at 63 samples is most of
it. Under 10 GB.

### 4.7 Agreement with GATK on the tomato cohort

There is no truth set, so this benchmark counts agreement, not accuracy. It uses the 63-sample calls
of §4.6, from pvc fitted, pvc with `--defaults`, and GATK, so it runs no caller of its own.

**How the calls are compared.** Both callers' VCFs are normalised the same way: left-aligned,
multi-allelic records split, and complex records decomposed (`bcftools norm -m -any -a`). The unit
counted is one sample's genotype at one site. Sites are counted separately for SNPs and indels, and
for inside and outside the STR segments pvc types.

**What is reported per sample-site:**
- both callers give the same genotype;
- both call a variant, with different genotypes;
- only pvc calls a variant;
- only GATK calls a variant.

**What exists.** Only a two-version comparison script in another worktree
(`pop_var_caller-ref-reads-in-repeats/tmp/tomato/gatk_cmp.py`), which counts moves toward or away
from GATK, and the allele-level marimo dashboard. The suite's scorer is new
(`benchmarks/suite/lib/agreement.py`).

**Estimate.** Minutes.

### 4.8 The colleague's chr22 trio: not run here

This benchmark is HG002, HG003 and HG004 over all of chromosome 22, at 5, 10, 30× and full depth. Its
data are not on this laptop. The suite has an entry for it that runs when pointed at the data. The
entry takes a configuration file naming:

- the three alignment files at each depth;
- the truth files;
- the confident regions.

It then runs §4.1's arms on them. The build report will list this benchmark as not run.

The 2,169-sample cohort on `kimura` is left out of this plan (§2).

## 5. Scoring

### 5.1 SNPs and indels: vcfeval, one way for every caller

Each run uses vcfeval with `--ref-overlap`, against the truth, inside the benchmark's confident
regions. I chose vcfeval directly over hap.py, because hap.py is built only for x86 and this machine
cannot emulate x86. hap.py's own vcfeval engine gives the same matching.

vcfeval writes the matched and unmatched records as four VCFs:
- truth variants found;
- calls that match the truth;
- false calls;
- truth variants missed.

The suite counts each by SNP and indel. vcfeval does not split by class itself.

**Metrics per caller, sample, depth and class:**
- **recall**: truth variants found, over all truth variants;
- **precision**: calls that match the truth, over all calls;
- **F1**: the balance of the two;
- the three counts behind them;
- the same at the best quality threshold (D2).

**Two kinds of match:**
- **genotype match**, vcfeval's default: the call's alleles and its genotype must both agree with
  the truth;
- **allele match** (`--squash-ploidy`): the right variant, but not necessarily the right genotype.

The gap between the two is the number of genotype errors.

### 5.2 STRs

STR genotypes are judged by the scorers built for them:
- `tract_qual_experiment.py` on HG002;
- `ng_str_stutter_score.py` on the simulation.

These compare whole STR alleles, which a record-by-record matcher cannot. Each is run unchanged
except for what §4.4 and §4.5 list.

## 6. Results, options and layout

**Where results go.** Everything is under `benchmark_results/` at the repository root, which is
added to `.gitignore`. It holds three directories, which change at different rates:

| directory | what it holds | rerun when |
|---|---|---|
| `inputs/` | inputs the suite makes: downsampled alignments (§4.2, §4.3) and the simulated genomes, reads and alignments (§4.5) | never, unless deleted; they are deterministic |
| `competitors/` | each competitor's calls, its time and memory, and its scored rows, one directory per benchmark, configuration and caller | only on `--rerun-competitors`, a missing directory, or a changed fingerprint (below) |
| `pvc/<commit>/` | pvc's calls, time, memory and scored rows for one pvc commit, then that commit's `results.tsv` and `summary.md` | every new pvc commit |

**Results for one commit.** A commit's `results.tsv` joins its own pvc rows with the competitors'
stored rows. Each competitor row keeps the date it was measured, so the table shows that GATK's
figures date from the first run.

**The competitor fingerprint.** Stored calls are reused only if the fingerprint written beside them
still matches. It covers three things:

- the tool's version string;
- the command line, with output paths left out;
- each input file's path, size and modification time. Small files, such as region lists and
  references' indexes, are checksummed instead.

A mismatch is never silently rerun. The suite stops and names the file or version that changed. It
reruns only when `--rerun-competitors` is given, because a competitor rerun costs hours.

`benchmark_results/latest` points at the newest commit's directory.

**A full run.** About 50 GB with the simulated read alignments kept, or about 20 GB without them
(§4.5). A later pvc refresh adds about 3 GB a commit, which is pvc's calls and psps. Old commits'
directories can be deleted freely.

### 6.1 Stopping and resuming

A run is a list of small steps. One step is one caller's one stage on one configuration, for
example:

- GATK HaplotypeCaller on HG003 at 10×;
- pvc `estimate-parameters` on the 50-sample PCR run at 30×;
- scoring one caller's calls.

**How a step finishes.** Each step writes into a temporary directory and renames it into place only
when it has finished and its outputs pass a check: the files exist, are not empty, and a VCF parses.
It then writes a `.done` marker holding its measurement. The longest single steps are expected to be
GATK HaplotypeCaller at 300× over the 6 Mb tandem-repeat regions, and freebayes on 63 tomato
samples, at about an hour each.

**Two ways to stop:**

- **`run_suite.sh --stop`** stops after the step that is running. It creates
  `benchmark_results/STOP`, which the suite checks between steps, and the next run removes it.
- **Ctrl-C, or closing the terminal,** stops at once. The container is stopped, and the step that
  was running is discarded and starts over next time.

**Resuming.** Run the same command again. Every step with a `.done` marker is skipped, and its
result is read from the marker. Nothing finished is rerun, so the cost of stopping is at most the
step that was cut short.

**What stopping does to the speed figures.** A step timed while the laptop was busy is wrong, and the
suite cannot prevent that. A host-side process records the host's load every 30 seconds into the
results directory for as long as the suite runs. Each step's row carries the mean host load over its
run. The summary marks any speed figure whose step ran with a mean load above 7, which is the 6
cores plus one. Such a step can be rerun on its own with `--redo-step <name>`.

**Commands:**

```
benchmarks/suite/run_suite.sh                        # everything not already done; pvc at this commit
benchmarks/suite/run_suite.sh --stop                 # finish the running step, then stop
benchmarks/suite/run_suite.sh --only giab_panel --depth 30 --caller pvc
benchmarks/suite/run_suite.sh --rerun-competitors --only hg002_tr --caller hipstr
benchmarks/suite/run_suite.sh --only stutter_sim --slice     # the small test of §7
benchmarks/suite/run_suite.sh --status               # steps done, remaining, and estimated time left
benchmarks/suite/run_suite.sh --summary-only         # rebuild summary.md from results.tsv
```

**Options.**

- `--only` takes a comma-separated list of benchmarks; `--depth` and `--caller` likewise.
- `--threads` defaults to 6.
- `--data-root` points at the checkout that holds `benchmarks/`. It defaults to this checkout; from a
  worktree it must be the main checkout.

**Files:**

```
benchmarks/suite/
  run_suite.sh            host side: parse options, build pvc, record the machine, start each benchmark
  fetch_tools.sh          the pinned RTG Tools download
  bench/                  one script per benchmark: giab_panel.sh, giab_trio.sh, hg002_5mb.sh,
                          hg002_tr.sh, stutter_sim.sh, speed.sh, tomato_agreement.sh,
                          external_chr22.sh
  lib/
    callers.sh            how each caller is run, so every benchmark runs them the same way
    steps.sh              the step runner of §6.1: markers, temporary directories, the stop file
    host_load.sh          the host-side load recorder of §6.1
    measure.py            time and memory (§3.3)
    vcfeval_score.sh      §5.1
    mendelian.sh          §4.2
    agreement.py          §4.7
    cache.sh              the reuse fingerprint
    summarise.py          results.tsv → summary.md
  README.md               how to run it, what each benchmark measures, and how long it takes
```

## 7. Build order

There is one commit per step. Each step is tested on a small slice: one region (or one 100 kb
stretch), one depth, and one sample where the benchmark allows it. Every arm and the scorer must run
end to end and put rows in the table. No full run is launched.

1. **The frame:**
   - `run_suite.sh` and its options;
   - the step runner, the stop file and the host load recorder;
   - `measure.py`, `cache.sh` and `summarise.py`;
   - `fetch_tools.sh`;
   - `.gitignore`.

   *Slice:* two steps that time `samtools faidx`, run through the whole frame to a table and
   summary. The run is stopped once with `--stop` and once with Ctrl-C, and resumed after each. The
   check is that finished steps are not rerun.
2. **GIAB per-sample panel**, with the shared caller runners and vcfeval scoring. *Slice:* HG002, 10×,
   one region.
3. **GIAB trio** and the Mendelian check. *Slice:* 10×, one region.
4. **HG002 over 5 Mb.** *Slice:* 30×, ten intervals.
5. **HG002 tandem-repeat benchmark.** *Slice:* 10×, the first 200 regions.
6. **Simulated stuttering benchmark**, with the driver changes of §4.5. *Slice:* a 200 kb genome,
   2 samples, 5×, both designs.
7. **Speed and tomato agreement.** *Slice:* 2 samples, one region. This step also measures one
   freebayes process's memory at 63 samples over one region, to settle the §3.4 risk.
8. **The chr22 entry**: written and checked to refuse cleanly when its data are absent.
9. **Stale files fixed:**
   - the `benchmarks/tomato1/bench.config.sh` header;
   - `benchmarks/ssr_hg002/README.txt`;
   - `benchmarks/giab/README.md`;
   - `benchmarks/lib/README.md` and the comment in `run_ng.sh` that says no command writes a fitted
     file;
   - the `run_perf_all.sh` header.

   Each points to the suite. The old helper scripts that only the deleted caller used are listed in
   the suite's README as superseded, not deleted; deleting them is a separate decision.

## 8. A full run

There are two kinds of run. **The first run** makes the inputs and runs every competitor and pvc.
**A refresh** for a new pvc commit runs pvc alone, on the stored inputs, and scores it.

All times are at 6 cores.

| benchmark | first run | refresh (pvc only) | disk, first run |
|---|---:|---:|---:|
| GIAB per-sample panel | 2 h | ½ h | 3 GB |
| GIAB trio | 1 h | ¼ h | 2 GB |
| HG002, 5 Mb | 1 h | ¼ h | 2 GB |
| HG002 tandem repeats | 4 h | 1 h | 3 GB |
| simulated stuttering grid | 16 h | 2½ h | 35 GB |
| speed, tomato | 4 h | ¾ h | 10 GB |
| tomato agreement | minutes | minutes | — |
| **total** | **about 28 h** | **about 5 h** | **about 55 GB** |

- **Where the first run's time goes.** Most of it is the competitors and the simulation, which are
  made once.
- **Where a refresh's time goes.** Mostly pvc's fit and calls on the 50-sample simulated cohorts.
- **Uncertainty.** These figures are scaled from the 18-core estimate of the first draft, by about
  1⅓, because the steps that run several processes at once get a third as many cores. The first run
  could fall anywhere from 20 to 40 hours, and a refresh from 3 to 8. The slice tests give
  per-region and per-sample rates, and both estimates are redone from them before any full run.
- **Running it in parts.** Either kind can be stopped and resumed as often as needed (§6.1). One
  benchmark at a time also works, through `--only`.
- **Space.** The disk holds 487 GB free.
