# pvc

pvc is a variant caller for SNPs, indels and STRs, built for cohorts of one to thousands of
samples. Samples can be added to a cohort later without reading the existing samples' alignments
again.

- **Speed.** On 50 tomato accessions at about 3 reads a position, over 91 Mb of chromosome 1, pvc
  took 40 minutes to 1 h 25 min, depending on how it was run, against 5 hours for GATK and nearly
  10 for freebayes. It called 2,176 whole tomato genomes, at the same depth, in 9 days on one
  server, using at most 33 GB of memory, with imputation off. On one human genome freebayes is
  5 to 7 times as fast as pvc ([performance](doc/user/performance.md)).
- **Accuracy.** On GIAB's HG002 at 30×, pvc makes about 11% fewer errors than GATK and a third
  fewer than freebayes; DeepVariant makes about a third fewer than pvc. At 5× pvc calls the most
  false SNPs of the four. On GIAB's tandem-repeat benchmark pvc and DeepVariant get about the same
  share of STR genotypes right, 88% and 87% at 30×, against 82% for GATK
  ([accuracy](doc/user/accuracy.md)).
- **Imputation**, optional, fills and corrects genotypes from the haplotypes the samples share
  (the technical pages call it refinement). In a simulated F2 of 100 plants at 1×, it cut the wrong
  SNP genotypes about tenfold, from 15,249 to 1,451
  ([imputation](doc/user/neighbour_evidence.md)).
- **Overlapping variants** are written as one record, without the `*` allele that many
  population-genetics tools misread ([overlapping variants](doc/user/overlapping_variants.md)).

> **Status: beta.** pvc gives good results on our benchmarks, but we have only started testing it
> on real projects. There is no release yet: the results quoted here were measured on the git
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
- **Filters out calls caused by hidden duplications** (places where two copies of a sequence
  exist in the genome but only one in the reference), so both copies' reads pile onto one place
  and look like heterozygous variants.
- **Imputes genotypes from the neighbouring loci**, when asked to (`--refine-genotypes`; off by
  default): in a cohort, it re-calls each sample from the haplotypes it shares with other samples,
  which corrects genotypes and fills those with no reads.
- **Writes overlapping variants as one record**, such as a SNP inside a deletion, whose alleles are
  the whole haplotypes, so it never needs the `*` allele. See
  [overlapping variants](doc/user/overlapping_variants.md).

## When pvc is not the right tool

- **One human genome at high depth.** On HG002 at 30×, DeepVariant makes 28,246 errors, false
  and missed variants together, where pvc makes 41,653; and
  freebayes takes 33 minutes where pvc takes 2 h 56 min.
- **Structural variants**, and deletions longer than 100 bases: pvc does not call them. Nor does
  it call satellites, repeats longer than 100 bases.
- **Polyploids, if you want the parameters estimated.** pvc calls any ploidy from 1 to 16, but
  `estimate-parameters` fits diploid samples only, so other ploidies are called with the built-in
  values (`--defaults`). Imputation also runs on diploid samples only.
- **Long reads.** pvc has been built and tested on Illumina short reads. It does not refuse
  PacBio or Oxford Nanopore reads, but it has not been tested on them.
- **Analyses that need invariant sites.** pvc writes no invariant sites and no gVCF: a site is
  written only if at least one sample carries an alternative allele. Statistics such as
  nucleotide diversity (π) need the callable sites from elsewhere, for instance a depth mask.

[Limits and refusals](doc/user/limits.md) lists the rest.

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

`pvc <command> --help` describes every option. The [user guide](doc/user/calling.md) compares
the two modes, and [Recipes](doc/user/recipes.md) gives a recommended command for each kind of
study: one sample, a trio, an inbred cohort at low depth, an outbred cohort, a large cohort on a
cluster, a laptop.

## Documentation

- **User guide** — how to do each task and which options to change: [doc/user/](doc/user/).
- **Methods** — the statistical models and how they compare with other callers':
  [doc/user/methods.md](doc/user/methods.md).
- **Bug reports and questions** — please open an
  [issue](https://github.com/JoseBlanca/pop_var_caller/issues).

## How it was built

Most of pvc's code was written with an AI coding assistant (Claude Code). Each part of the
caller was first specified in a design document ([doc/devel/ng/spec/](doc/devel/ng/spec/)), then
implemented and reviewed, and checked against freebayes, GATK, DeepVariant and HipSTR on real and simulated
data. Where pvc adapts a method from another caller, it was ported from that caller's source
code, which the repository keeps for reference.

## License

MIT — see [LICENSE](LICENSE).
