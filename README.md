# pvc

pvc is a variant caller for SNPs, indels and STRs, built for cohorts of one to thousands of
samples.

- **Speed.** On 50 tomato accessions at about 3 reads a position, over 91 Mb of chromosome 1, pvc
  took 40 minutes to 1 h 25 min, depending on how it was run, against 5 hours for GATK and nearly
  10 for freebayes. It called 2,176 whole tomato genomes, at the same depth, in 9 days on one
  server, using at most 33 GB of memory, with imputation off
  ([performance](doc/user/performance.md)).
- **Accuracy.** On GIAB's HG002 at 30×, pvc makes about 11% fewer errors than GATK and a third
  fewer than freebayes; DeepVariant makes about a third fewer than pvc. On GIAB's tandem-repeat
  benchmark pvc and DeepVariant get about the same share of STR genotypes right, 88% and 87% at
  30×, against 82% for GATK ([accuracy](doc/user/accuracy.md)).
- **Imputation**, optional, fills and corrects genotypes from the haplotypes the samples share
  (the technical pages call it refinement). In a simulated F2 of 100 plants at 1×, it cut the wrong
  SNP genotypes about tenfold, from 15,249 to 1,451
  ([imputation](doc/user/neighbour_evidence.md)).
- **Overlapping variants**, such as a SNP inside a deletion, are written as one record whose
  alleles are the whole haplotypes, so pvc never needs the `*` allele, which causes information
  loss in population-genetics analyses ([overlapping variants](doc/user/overlapping_variants.md)).

> **Status: beta.** pvc gives good results on our benchmarks, but we have only started testing it
> on real projects. There is no release yet: the results quoted here were measured on git
> commit 909ea135, and options and the output format may still change. No paper describes pvc
> yet, so there is nothing to cite. We welcome feedback and bug reports: please open an
> [issue](https://github.com/JoseBlanca/pop_var_caller/issues) with what you ran and what you
> saw.

Supported on Linux and macOS. Licence: MIT.

## What it does

- **Calls** SNPs, indels (up to 100 bp by default), and STRs with repeat units of 1 to 6 bases
  spanning up to 100 bp by default.
- **Lets you add samples to a cohort later** without re-reading the alignments of the samples you
  already have. To do that, it generates intermediate pileup files, which are then reused for
  the cohort calling.
- **Estimates the calling model's priors from the data** before calling: base-quality calibration
  and contamination for each read group, the inbreeding coefficient F of each sample,
  the population's diversity, and stutter for each kind of repeat.
- **Filters out calls caused by hidden duplications**: places where the genome holds two copies
  of a sequence but the reference only one, so the reads of both copies pile onto one place and
  look like heterozygous variants.
- **Imputes genotypes**, when asked to (`--refine-genotypes`; off by default): in a cohort, it re-calls each sample from the haplotypes it shares with other samples,
  which corrects genotypes and fills those with no reads.

## Install

pvc is installed from source for now. You need:

- Rust 1.98 or newer, from [rustup](https://rustup.rs);
- a C compiler (`cc` on Linux, the Xcode command-line tools on macOS).

```
cargo install --locked --git https://github.com/JoseBlanca/pop_var_caller
```

This puts the `pvc` binary in `~/.cargo/bin`. A prebuilt Linux binary for testing is planned.

## Quick start

pvc needs the reference FASTA the reads were aligned to, and coordinate-sorted, indexed BAM or
CRAM files. Each file's read groups (`@RG`) must name their sample with `SM`.

First, once per reference, find the reference's tandem repeats:

```
pvc repeat-catalog --reference genome.fa
```

Then call in one of two ways.

![Direct mode: every BAM goes into one calling step, then optional imputation (refinement). psp mode: each BAM becomes its own psp file; the psp files go to optional parameter estimation and to calling, then optional imputation. A sample added later needs only its own new psp file](doc/user/figures/calling_modes.svg)

**Direct mode: one step, with built-in parameter values.** Quickest; reads every sample's
alignments in one run.

```
pvc call-from-alignments --reference genome.fa \
    --alignment s1.bam --alignment s2.bam \
    --defaults --output calls.vcf.gz
```

**psp mode: three steps, with parameters estimated from your data.** First summarise each
sample's reads into a per-sample pileup file (*psp*), then estimate the parameters from the
cohort, then call. A sample added later needs only its own psp.

```
pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/   # once per sample
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

Given the same samples and the same parameters, the two modes write the same calls.

`pvc <command> --help` describes every option. The [user guide](doc/user/calling.md) compares
the two modes, and [Recipes](doc/user/recipes.md) gives a recommended command for each kind of
study: one sample, a trio, an inbred cohort at low depth, an outbred cohort, a large cohort on a
cluster, a laptop.

## Documentation

- **[contents page](doc/user/README.md)** — how to do each task and which options to change.
- **Methods** — the statistical models and how they compare with other callers':
  [doc/user/methods.md](doc/user/methods.md).
- **Bug reports and questions** — please open an
  [issue](https://github.com/JoseBlanca/pop_var_caller/issues).

## When pvc might not be the right tool

- **A single sample at high coverage, with no stutter.** For instance, on HG002 at 30×,
  DeepVariant makes 28,246 errors, false and missed variants together, where pvc makes 41,653.
- **Structural variants.**
- **Polyploids, if you want the parameters estimated.** pvc calls any ploidy from 1 to 16, but
  `estimate-parameters` fits diploid samples only, so other ploidies are called with the built-in
  values (`--defaults`). Imputation also runs on diploid samples only.
- **Long reads.** pvc has been built and tested on Illumina short reads. It does not refuse
  PacBio or Oxford Nanopore reads, but it has not been tested on them.
- **Analyses that need invariant sites.** pvc writes no invariant sites: a site is written only
  if at least one sample carries an alternative allele.

[Limits and refusals](doc/user/limits.md) lists the rest.

## How it was built

Most of pvc's code was written with an AI coding assistant (Claude Code). Each part of the
caller was first specified in a design document ([doc/devel/ng/spec/](doc/devel/ng/spec/)), then
implemented and reviewed, and checked against freebayes, GATK, DeepVariant and HipSTR on real and
simulated data.

## License

MIT — see [LICENSE](LICENSE).
