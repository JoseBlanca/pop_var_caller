# pvc

pvc is a variant caller for SNPs, indels and STRs, built for cohorts of one to thousands of
samples. In cohorts it is faster than GATK and freebayes: with 50 tomato samples it took 40
minutes, against 5 hours for GATK and nearly 10 for freebayes. It called 2,176 whole tomato
genomes in 9 days on one server, using at most 33 GB of memory, and samples can be added to a
cohort later without reading the existing samples' alignments again
([performance](doc/user/performance.md)). On GIAB's HG002 at 30× it is more accurate than GATK
and freebayes and close to DeepVariant, and it gets the most STR genotypes right of the four
([accuracy](doc/user/accuracy.md)). Overlapping variants are written as one record, without the
`*` allele that many population-genetics tools misread
([overlapping variants](doc/user/overlapping_variants.md)). Thanks to its optional
[imputation](doc/user/neighbour_evidence.md) module it is quite capable even in low-coverage
samples: for instance, in a simulated F2 of 100 plants at 1×, it cut the wrong SNP genotypes
twelvefold.

> **Status: under active development.** pvc gives good results on our benchmarks, but we have
> only started testing it on real projects. Options and the output format may still change. We
> welcome feedback: please open an
> [issue](https://github.com/JoseBlanca/pop_var_caller/issues) with what you ran and what you
> saw.

Supported on Linux and macOS.

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
- **Refines genotypes with the neighbouring loci**: in a cohort, it re-calls each sample from the
  haplotypes it shares with other samples, which corrects genotypes and fills those with no reads.
- **Writes overlapping variants as one record**, such as a SNP inside a deletion, whose alleles are
  the whole haplotypes, so it never needs the `*` allele. See
  [overlapping variants](doc/user/overlapping_variants.md).

It does not call structural variants or satellite repeats.

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
the two modes.

## Documentation

- **User guide** — how to do each task and which options to change: [doc/user/](doc/user/)
  *(in preparation)*.
- **Methods** — the statistical models and how they compare with other callers':
  [doc/methods/](doc/methods/) *(in preparation)*.
- **Contributing** — building, testing and the design documents:
  [CONTRIBUTING.md](CONTRIBUTING.md) *(in preparation)*.

## How it was built

Most of pvc's code was written with an AI coding assistant (Claude Code). Each part of the
caller was first specified in a design document ([doc/devel/ng/spec/](doc/devel/ng/spec/)), then
implemented and reviewed, and checked against freebayes, GATK and HipSTR on real and simulated
data. Where pvc adapts a method from another caller, it was ported from that caller's source
code, which the repository keeps for reference.

## License

MIT — see [LICENSE](LICENSE).
