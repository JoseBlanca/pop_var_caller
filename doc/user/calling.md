# Calling

pvc calls a cohort in one of two ways. Both write one VCF holding SNPs, indels and STRs, and both
need the [repeat catalogue](inputs.md#the-repeat-catalogue) first.

| | Direct mode | psp mode |
|---|---|---|
| commands | `call-from-alignments` | `generate-psps`, then `estimate-parameters`, then `call-from-psps` |
| parameters | built-in values | estimated from your reads |
| reads each alignment | once per call | once per sample, ever |
| adding a sample later | call the whole cohort again from the alignments | write one psp file, then fit and call again |
| best for | a quick look, test regions, outbred samples | inbred or structured material, large cohorts, cohorts that grow |

**Choose psp mode for inbred lines, selfing crops or structured populations.** Direct mode uses
built-in values. Among them, every sample's inbreeding coefficient is 0, meaning outbred. On inbred
material that assumption makes pvc call too many heterozygotes. psp mode estimates each sample's
inbreeding from its reads, along with the other parameters listed in
[Estimated parameters](parameters.md).

Given the same samples and the same parameters, the two modes write the same calls.

## Direct mode

```
pvc call-from-alignments --reference genome.fa \
    --alignment s1.bam --alignment s2.bam --alignment s3.bam \
    --defaults --output calls.vcf.gz
```

- `--defaults` uses the built-in parameter values. The alternative, `--parameters`, takes a
  parameters file, but a parameters file is bound to the samples it was estimated on (see
  [Estimated parameters](parameters.md#what-a-parameters-file-is-bound-to)). In practice direct
  mode runs with `--defaults`.
- The VCF's sample columns follow the order of the `--alignment` options.
- `--regions regions.bed` calls only those regions.

## psp mode

### 1. One psp file per sample

```
pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/
```

This writes `psps/<sample>.psp` and its log, `psps/<sample>.psp.log`. A sample whose reads are
split over several files gets one psp from all of them: pass every file with `--alignment`.

Several samples can go to one command, but they are processed one after another, each on about
one core. To use more cores or machines, run one command per sample at the same time, **all with
the same `--regions`**, or none. psp files written over different regions are refused when they
are combined.

If a psp for a sample is already in the directory, the command refuses before reading anything.
`--force` replaces it. A psp appears only once it is complete, so an interrupted run leaves no
half-written file.

### 2. Estimate the parameters

```
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
```

`--psp` takes psp files or directories; a directory contributes every `.psp` file in it. This step
reads only a small summary stored at the end of each psp, not the whole file.

On 59 tomato accessions over 32 Mb at about 5.5 reads a position, this step took 16 minutes and
410 MB on an 18-core laptop. One sample took 24 seconds. *(Measured 2026-10-17; these figures
will be replaced by the benchmark suite's.)*

Options worth knowing:

- `--ploidy N` if your samples are not diploid. The calling step then uses the file's ploidy.
- `--inbreeding F` states one inbreeding coefficient for every sample instead of estimating it,
  for example `--inbreeding 1` for fully inbred lines. Leave it out to have each sample's value
  estimated.
- `--skip-contamination` if the samples are known to be clean. It saves time, and the file then
  holds no contamination table.

### 3. Call

```
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

The regions called, and what counts as a repeat, come from the psp files. `call-from-psps` has
no `--regions`; to call less of the genome, write the psps over less.

## Adding samples to a cohort

Write a psp for each new sample, with the same `--regions` as the others. Then run
`estimate-parameters` and `call-from-psps` again over the whole directory. The existing samples'
alignments are not read again.

The parameters must be estimated again, because a parameters file is bound to the exact set of
samples it was estimated on.

## Large cohorts

- **Use psp mode.** Writing psps can be spread over as many machines as you like, one sample per
  command, and the calling step keeps one open file per sample instead of two or three.
- **Raise the open-file limit** if pvc asks you to (`ulimit -n`). It checks before starting and
  says which value to set.
- **Memory while calling** is held to about 500,000 read observations at a time across the
  cohort, which was about 300 MB on the tomato benchmark. It is designed to stay about flat as
  samples are added. So far that has been measured from 1 to 63 samples.
- **Disk:** `call-from-psps` reads up to 8 GiB of psp data ahead into the operating system's
  file cache.

## Output, logs and scratch space

- `--output calls.vcf.gz` (or `.vcf.bgz`) writes a bgzipped VCF, and any other name writes plain
  text. `--output -` writes uncompressed VCF to standard output.
- The VCF is written under a temporary name and renamed when the run ends, so a finished file is
  always complete. A failed run leaves `<output>.tmp` behind.
- No command overwrites an existing output or log without `--force`, and it refuses before
  reading any data.
- Every run writes a log: `pvc_calls.log` beside the VCF (change with `--log`). The log records
  the command, every warning, the parameter values the run used and its peak memory. The terminal
  shows only what needs attention; `--quiet` reduces it to warnings and errors.
- `--threads` sets the number of threads (0, the default, uses every core). The calls are the
  same at any number.
- **Scratch space.** By default pvc checks every call for hidden duplications (see
  [Reading the VCF](vcf.md#filter)), and that check needs the whole genome's calls before it can
  write any. Until then the calls wait in a scratch directory, by default in the VCF's directory
  (change with `--tmp-dir`), and take several times the size of the final compressed VCF. For the
  same reason, with `--output -` nothing reaches standard output until the whole genome is called.
  `--hidden-dup-fdr 0` turns the check off, and the scratch space with it.
