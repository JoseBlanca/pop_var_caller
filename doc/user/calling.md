# Calling

pvc usually splits calling into two phases:

1. gathering the sequence evidence into one pileup file (a psp file) per sample;
2. the calling itself, from the psp files.

Alternatively, it can call directly from the alignment files (BAM or CRAM). This mode does less,
because the parameters used as priors can only be estimated from psp files.

For the commands to run in common cases — one sample, a trio, an inbred or an outbred cohort, a
cluster, a laptop — see [Recipes](recipes.md).

| | Direct mode | psp mode |
|---|---|---|
| commands | `call-from-alignments` | `generate-psps`, then `estimate-parameters`, then `call-from-psps` |
| parameters | built-in values, or from a previous run | estimated from the psp files |
| reads each alignment | once per call | once per sample, ever |
| adding a sample later | call the whole cohort again from the alignments | write one psp file, then fit and call again |

![Direct mode: every BAM goes into one calling step, then optional imputation (refinement). psp mode: each BAM becomes its own psp file; the psp files go to optional parameter estimation and to calling, then optional imputation. A sample added later needs only its own new psp file](figures/calling_modes.svg)

Given the same samples and the same parameters, the two modes write the same calls.

## Direct mode

```
pvc call-from-alignments --reference genome.fa \
    --alignment s1.bam --alignment s2.bam --alignment s3.bam \
    --defaults --output calls.vcf.gz
```

- `--defaults` uses the built-in parameter values. The alternative, `--parameters`, takes a
  parameters file.
- The VCF's sample columns follow the order of the `--alignment` options.
- `--regions regions.bed` calls only those regions (the BED format is in
  [Preparing the inputs](inputs.md#regions)).
- **Direct mode keeps few cores busy, whatever `--threads` says.** On the tomato cohort over 91 Mb
  of chromosome 1, with 30 threads, it kept about 1.4 cores busy with one sample and about 7 with
  100; on HG002 at 30×, about 1.3. To use more cores, call regions at the same time (see
  [Splitting a run by region](#splitting-a-run-by-region)), or use psp mode.

## psp mode

### 1. One psp file per sample

```
pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/
```

This writes `psps/<sample>.psp` and its log, `psps/<sample>.psp.log`. A sample whose reads are
split over several alignment files gets one psp from all of them: pass every file with
`--alignment`.

Several samples can go to one command, but they are processed one after another. To use more
cores or machines, run one command per sample at the same time, all with the same `--regions`.
psp files written over different regions are refused when they are combined.

If a psp for a sample is already in the directory, the command refuses before reading anything.
`--force` replaces it. A psp is written under a temporary name, `<sample>.psp.<number>.partial`,
and takes its final name only once it is complete. A run that fails deletes its partial file; a
run that is killed leaves it behind for you to delete.

### 2. Estimate the parameters

```
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
```

`--psp` takes psp files or directories; a directory contributes every `.psp` file in it. This step
reads only a small summary stored at the end of each psp, not the whole file.

Options worth knowing:

- `--skip-contamination` if the samples are known to be clean. It saves time, and the file then
  holds no contamination table.
- `--start-from` an earlier parameters file, when samples have been added (see
  [Adding samples](#adding-samples-to-a-cohort)).

The log is `pvc_call_parameters.log` beside the parameters file (change with `--log`).

### 3. Call

```
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

The regions called, and what counts as a repeat, come from the psp files. `call-from-psps` has
no `--regions`; to call less of the genome, write the psps over less.

`--defaults` in place of `--parameters` calls the psps with the built-in values, as in direct
mode. It is the way to call psps at a ploidy other than 2, which `estimate-parameters` does not
fit.

## Adding samples to a cohort

Write a psp for each new sample, with the same `--regions`, repeat catalogue and repeat settings
as the others; psps that differ in any of these are refused together. Then run
`estimate-parameters` and `call-from-psps` again over the whole directory. The existing samples'
alignments are not read again.

The parameters must be estimated again, because a parameters file is bound to the exact set of
samples it was estimated on. `--start-from` shortens that:

```
pvc estimate-parameters --reference genome.fa --psp psps/ \
    --start-from cohort.parameters.toml --output cohort2.parameters.toml
```

The fit then starts from the earlier numbers, and reaches the same answer as a fit from scratch.
On 400 tomato accessions, adding 20 took 43 minutes instead of 1 h 20 min (see
[Options](options.md#estimating-the-parameters-estimate-parameters)).

A psp written by a pvc whose psp format differs from this one is refused, and the message says
so. Write it again from its alignments.

## Large cohorts

- **Use psp mode.** Writing psps can be spread over as many machines as you like, one sample per
  command, and the calling step keeps one open file per sample instead of 8 to 15.
- **Raise the open-file limit** if pvc asks you to (`ulimit -n`). It checks before starting and
  says which value to set.
- **Memory while calling.** pvc builds loci a stretch of genome at a time, and shortens the
  stretch as samples are added, so that it holds about 500,000 read observations at once. Up to
  about 1,000 samples this keeps that part of the memory flat. Beyond that the stretch is at its
  shortest, 500 bases, and memory grows with the number of samples.
- **Disk:** `call-from-psps` reads up to 8 GiB of psp data ahead into the operating system's
  file cache. On a slow disk, `--psp-prefetch-bytes` raises it.
- **Splitting the calling by chromosome** over several machines is in
  [Recipes § A large cohort on a cluster](recipes.md#a-large-cohort-on-a-cluster).

## Splitting a run by region

A run can be cut into regions, each called by its own command at the same time, and the VCFs
joined. In direct mode this uses more cores: on the tomato cohort, 100 samples over 91 Mb, three
commands of 10 threads each took 1 h 25 min, against 3 h 5 min for one command with 30 threads.
Each command holds its own memory, so the three together used 23 GB against 10 GB.

```
pvc call-from-alignments --reference genome.fa --alignment s1.bam --alignment s2.bam \
    --defaults --regions chr01.bed --output chr01.vcf.gz --log chr01.log --threads 10
pvc call-from-alignments --reference genome.fa --alignment s1.bam --alignment s2.bam \
    --defaults --regions chr02.bed --output chr02.vcf.gz --log chr02.log --threads 10

bcftools concat -O z -o calls.vcf.gz chr01.vcf.gz chr02.vcf.gz
bcftools index calls.vcf.gz
```

- **Give each command its own `--log`**, or its own output directory. The log is
  `pvc_calls.log` beside the VCF, and a run will not start over the log of a run that finished
  or is still running, so a second command writing into the same directory stops at once. The
  scratch directories need no care: each run makes its own inside `--tmp-dir`.
- **Give every command the `--alignment` options in the same order**, so that the VCFs have
  the same sample columns, and give the VCFs to `bcftools concat` in the reference's order.
- **Cut between chromosomes where you can.** A cut inside a chromosome is an edge for both
  commands: a repeat that crosses it is called by both, so it is written twice, and an indel that
  crosses it in non-repeat sequence is called by neither.
- **In psp mode**, `call-from-psps` has no `--regions`, so each region needs its own psps, written
  with `generate-psps --regions`. [Recipes](recipes.md#a-large-cohort-on-a-cluster) shows how.

**Splitting can change some calls.** The genotypes are called locus by locus, so a split does
not change them, but two later steps look beyond the locus, and a split run sees only its own
regions:

- the hidden-duplication filter learns each sample's usual depth, and how common hidden
  duplications are, from all the run's calls;
- imputation, when it is on, reads up to a few megabases on each side of a site, and stops at the
  edges of the run's regions.

How many calls this changes has not been measured.

## Output, logs and scratch space

- `--output calls.vcf.gz` (or `.vcf.bgz`) writes a bgzipped VCF, and any other name writes plain
  text. `--output -` writes uncompressed VCF to standard output.
- The VCF is written under a temporary name, `<output>.tmp`, and renamed when the run ends, so a
  finished file is always complete. A failed run can leave `<output>.tmp` behind.
- No command overwrites an existing output or log without `--force`, and it refuses before
  reading any data. The exception is a log left by a failed run, which is replaced.
- Every run writes a log: `pvc_calls.log` beside the VCF (change with `--log`). The other
  commands keep theirs as well: `<sample>.psp.log` beside each psp, `pvc_call_parameters.log`
  beside the parameters file, and `<catalogue>.log` beside the repeat catalogue. The log records
  the command, every warning, the parameter values the run used and its peak memory. The terminal
  shows only what needs attention; `--quiet` reduces it to warnings and errors.
- `--threads` sets the number of threads (0, the default, uses every core). The calls are the
  same at any number.
- **Scratch space.** By default pvc checks every call for hidden duplications (see
  [The VCF](vcf.md#filter)), and that check needs the whole genome's calls before it can
  write any. Until then the calls wait in a scratch directory, by default in the VCF's directory
  (change with `--tmp-dir`), and take several times the size of the final compressed VCF. For the
  same reason, with `--output -` nothing reaches standard output until the whole genome is called.
  `--hidden-dup-fdr 0` turns the check off, and the scratch space with it.
- **A run that stops after calling can be finished without calling again.** With
  `--keep-tmp-on-error`, a run that stops on an error keeps its scratch directory, and its log
  names it. A run stopped because an imputation window did not fit `--refine-window-bytes` keeps
  it even without the flag. A run killed outside pvc's control, for example for memory, leaves it
  behind too. Once calling has ended, that directory holds everything the hidden-duplication check
  and the writing need, and

  ```
  pvc resume-calling --scratch-dir <the directory> --output calls.vcf.gz
  ```

  finishes the run from it, as the stopped run would have. It keeps the stopped run's settings,
  except for `--threads`, the imputation memory options and `--phase-min-path-agreement`, which it
  may set anew (see [Options](options.md#finishing-a-stopped-run-resume-calling)). Nothing it does
  changes what was called. It removes the directory when it finishes.
