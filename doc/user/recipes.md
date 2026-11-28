# Recipes

A recommended sequence of commands for each common case. Every recipe assumes the repeat
catalogue has been built once for the reference (`pvc repeat-catalog --reference genome.fa`; see
[Preparing the inputs](inputs.md#the-repeat-catalogue)). The two modes are explained in
[Calling](calling.md), and the options in [Options](options.md).

Imputation — filling missing genotypes and correcting thin ones from the samples that share
stretches of chromosome — is what the technical pages call *refinement*, after its option
`--refine-genotypes`. It is off by default, and the recipes say where to turn it on.

Times and memory are given only where they were measured, with the data they were measured on.

## One sample

```
pvc call-from-alignments --reference genome.fa --alignment sample.bam \
    --defaults --output calls.vcf.gz
```

- **Direct mode with the built-in values is the quicker way.** On HG002 over chromosomes 1 to 22
  at 30×, it took 2 h 56 min and 1.3 GB, against 3 h 49 min for psp mode.
- **What psp mode would add is a base-quality calibration measured on your own libraries.** The
  built-in one doubles every base's error probability, a factor chosen on three human Illumina
  genomes. If your libraries are unlike those, use the psp-mode commands of the next recipe with
  one sample. Either way the sample's inbreeding coefficient is taken as 0 (see
  [Estimated parameters](parameters.md#what-is-estimated)), and there is nothing to impute from.
- **It keeps about one core busy**, whatever `--threads` says: 1.3 on HG002 at 30×. To finish
  sooner, call the chromosomes at the same time
  ([Splitting a run by region](calling.md#splitting-a-run-by-region)).

## A trio or a small family

```
pvc generate-psps --reference genome.fa --alignment child.bam --output-dir psps/
pvc generate-psps --reference genome.fa --alignment mother.bam --output-dir psps/
pvc generate-psps --reference genome.fa --alignment father.bam --output-dir psps/
pvc estimate-parameters --reference genome.fa --psp psps/ --output family.parameters.toml
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters family.parameters.toml --output calls.vcf.gz
```

- **pvc takes no pedigree.** The samples are called together as a small cohort, and how they are
  related is not used.
- The three `generate-psps` can run at the same time.
- **Leave imputation off.** `auto` would not run below 10 samples, and on a simulated outbred
  cohort of 10 samples it removed only 4 in 100 of the wrong genotypes (see
  [Genotypes from neighbouring sites](neighbour_evidence.md#when-it-runs)). The VCF is still
  phased from the reads where they link two heterozygotes.

## An inbred or selfing cohort at about 3×

```
pvc generate-psps --reference genome.fa --alignment S1.bam --output-dir psps/
pvc generate-psps --reference genome.fa --alignment S2.bam --output-dir psps/
...
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --refine-genotypes auto --output calls.vcf.gz
```

- **Use psp mode, and leave `--inbreeding` unset**, so that each sample's inbreeding coefficient
  is estimated. On a simulated tomato population of 50 samples at 3×, from fully outbred to fully
  inbred, each group's mean estimate was within 0.015 of the true value, and with the estimates
  84.6% of the non-reference genotypes were right, against 79.2% with the coefficient taken as 0
  ([Accuracy](accuracy.md#simulated-inbreeding-does-a-wrong-f-cost-genotypes)).
- **`--refine-genotypes auto` turns imputation on** when the cohort has at least 10 samples and
  their median inbreeding coefficient is at least 0.5; the VCF header says whether it ran. In
  the 100 samples at 2× of a simulated selfing population (inbreeding coefficient 0.9), it cut
  the missed SNP genotypes from 1,710 to 467 and the genotypes left uncalled from 32,687 to 734,
  while false variant genotypes rose from 73 to 178
  ([Accuracy](accuracy.md#imputation-genotype-refinement-and-filling)). Check
  `<name>.fill_check.tsv` beside the VCF afterwards
  ([the fill check](neighbour_evidence.md#how-good-the-fills-are-the-fill-check)).
- **Imputation memory.** It may hold up to `--refine-window-bytes`, 32 GB by default, at once; on
  2,176 tomato samples at 48 threads the whole process peaked at about 1.3 times that. On a
  smaller machine, lower it (see [A laptop](#a-laptop)).

## An outbred cohort

The commands of the inbred recipe, without `--refine-genotypes`:

```
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

- **Estimate the parameters here too.** In the simulated population of the inbred recipe, the
  outbred group's estimate at 3× was 0.007, against a true 0.005: the fit does not invent
  inbreeding.
- **Leave imputation off**, the default. On simulated outbred cohorts at 3× it removed 18% to 63%
  of the wrong genotypes at 50 samples, but only 4% at 10 samples. Outbred species have not been
  measured. To impute, run a tool that
  phases, such as Beagle or GLIMPSE, on pvc's VCF.

## A large cohort on a cluster

Each sample is processed on its own up to its psp file, the parameters are estimated once for
the whole cohort, and each chromosome is then called on its own.

**1. One job per sample**, writing that sample's psp over the whole genome:

```
pvc generate-psps --reference genome.fa --alignment S1.cram --output-dir psps/
```

**2. One job estimating the parameters**, over every psp:

```
pvc estimate-parameters --reference genome.fa --psp psps/ \
    --output cohort.parameters.toml --threads 16
```

One estimate for the whole genome gives the rarer STR classes more data than one estimate per
chromosome would.

**3. One BED per chromosome**, from the reference's index. Small contigs can go together in one
BED.

```
mkdir -p beds
awk -v OFS='\t' '{ print $1, 0, $2 > ("beds/" $1 ".bed") }' genome.fa.fai
```

**4. One calling job per chromosome**, each reading the same psps:

```
pvc call-from-psps --reference genome.fa --psp psps/ --parameters cohort.parameters.toml \
    --regions beds/chr01.bed --output chr01.vcf.gz --log chr01.log --threads 48
```

Give each job its own `--log`, or its own output directory, so the jobs do not refuse each
other's logs. A job reads only its own chromosome's part of each psp.

**5. Join the chromosomes**, in the reference's order:

```
bcftools concat -O z -o calls.vcf.gz chr01.vcf.gz chr02.vcf.gz ...
bcftools index calls.vcf.gz
```

Every chromosome's VCF has the same sample columns, because a directory's psps are taken in name
order.

- **Raise the open-file limit** on the calling jobs: `call-from-psps` keeps one file open per
  sample. pvc checks the limit before starting and says which value to set.
- **Splitting by chromosome can change a few calls**, because the hidden-duplication filter and
  imputation then see one chromosome at a time
  ([Splitting a run by region](calling.md#splitting-a-run-by-region)). A sex chromosome or an
  organelle then gets a depth model of its own, and a short contig gives the filter few calls to
  learn from.
- **Measured, on one 96-core server**, for 2,176 tomato genomes of 783 Mb each, without
  imputation and not split by chromosome: a psp took 19 minutes and 0.5 GB for the median sample,
  and at most 8 h 22 min and 8.3 GB; the parameters took 34 h at 16 threads; calling took 104 h at
  48 threads; no step used more than 33 GB. Jobs split by chromosome have not been measured.
- **Disk.** The scratch space of each calling job takes several times its final compressed VCF
  (`--tmp-dir` moves it). psp file sizes have not been measured yet.

## A laptop

```
pvc call-from-alignments --reference genome.fa --alignment s1.bam --alignment s2.bam \
    --defaults --threads 4 --output calls.vcf.gz
```

- **Direct mode needs little memory**: on HG002 over chromosomes 1 to 22, 0.9 GB at 5× to 1.8 GB
  at 50×; on the tomato cohort over 91 Mb of chromosome 1, 0.3 GB for one sample and 10 GB for
  100.
- **`--threads`** defaults to every core; a lower number leaves the machine usable.
- **The repeat catalogue** needs about 1.4 GB on the human genome at its default of one thread.
- **With imputation on, lower `--refine-window-bytes`.** Its default, 32 GB, is more than most
  laptops have. Choose a value that still fits when multiplied by 1.3, the factor measured on
  2,176 tomato samples: for example `--refine-window-bytes 8000000000` on a 16 GB machine.
  Before calling, pvc estimates the peak and warns if it is above 85% of the machine's memory,
  naming a value that fits. If a window does not fit even at its shortest, the run stops after
  calling and keeps its scratch directory; `pvc resume-calling` finishes it with a larger value or a smaller `--refine-core-floor`
  ([Output, logs and scratch space](calling.md#output-logs-and-scratch-space)).

## When something goes wrong

- **A run is slow.** Direct mode keeps one or a few cores busy, whatever `--threads` says: split
  the run by region ([Calling](calling.md#splitting-a-run-by-region)) or use psp mode, with one
  `generate-psps` per sample at the same time. psps on a spinning disk read faster with a larger
  `--psp-prefetch-bytes`.
- **A run runs out of memory.** During imputation, lower `--refine-window-bytes`; a run killed
  after calling can be finished with `pvc resume-calling`, which the warning names. During
  calling, above about 1,000 samples memory grows with the number of samples
  ([Calling § Large cohorts](calling.md#large-cohorts)).
- **The disk fills.** The hidden-duplication filter's scratch space takes several times the
  compressed VCF: move it with `--tmp-dir`, or turn the filter off with `--hidden-dup-fdr 0`.
- **A variant you expected is missing.** Where to look:
  - `--explain-loci regions.bed out.tsv`, on the calling command, writes, for each locus in the BED, the candidate alleles,
    the genotypes, the qualities and what happened to the record;
  - `pvc dev inspect-psp --psp sample.psp --regions regions.bed` shows whether the sample's reads
    reached the psp there;
  - `--min-site-quality 0` writes the weak sites, and `--hidden-dup-tag` keeps the calls the
    hidden-duplication filter removes, marked `hiddenDup`;
  - the variant may be where pvc does not call: inside a satellite (`--call-satellites`), a
    deletion longer than `--max-locus-span`, or a locus deeper than `--max-depth` (see
    [Limits](limits.md#what-pvc-does-not-call)), or in reads pvc ignores
    ([Read filters](inputs.md#read-filters)).
