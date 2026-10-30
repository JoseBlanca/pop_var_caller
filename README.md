# pvc

pvc is a variant caller, focused on speed, memory usage and accuracy.

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
  already have. To do that, you generate intermediate pileup files, which are then reused for
  the cohort calling.
- **Estimates the calling model's priors from the data** before calling: base-quality calibration
  and contamination for each read group, the inbreeding coefficient F of each sample,
  the population's diversity, and stutter for each kind of repeat. They are written to a
  parameters file that you can read, edit and reuse.
- **Estimates contamination**: what share of each read group's reads came from another individual.
  Calling then takes it into account.
- **Flags hidden duplications**: places where two copies of a sequence exist in the genome but
  only one in the reference, so both copies' reads pile onto one place and look like
  heterozygous variants.

It does not call structural variants, satellite repeats longer than 100 bp, or deletions longer
than 100 bp.

## How it performs

| | SNP F1, HG002 30× [^giab] | indel F1, HG002 30× [^giab] | STR genotypes right, HG002 30× [^str] | 100 tomato samples: time / memory [^tomato] | [overlapping variants](doc/user/overlapping_variants.md) |
|---|---|---|---|---|---|
| pvc | 0.9941 | 0.9819 | 88% | 1 h 45 min / 11 GB | one record, never `*` |
| GATK | 0.9928 | 0.9839 | 82% | 9 h 50 min / 75 GB | separate records, with `*` |
| freebayes | 0.9933 | 0.9606 | not run | stalled, stopped after 17 h | one record, which also merges variants within 3 bases |
| DeepVariant | 0.9966 | 0.9934 | 87% | not run | separate records, without `*`: a copy the deletion removes is written as missing in a cohort, as reference in a single sample |
| HipSTR | n/a | n/a | 73% | not run | STRs only: one record per repeat |

At larger scale, pvc called more than 2,000 whole tomato genomes in a week using less than 30 GB.

[Accuracy](doc/benchmarks/accuracy.md) has the full results behind the accuracy columns, at every
read depth, with the errors split by kind.

[^giab]: GIAB HG002, the whole of chromosome 20 (64 Mb), reads downsampled to 30×, scored with
    RTG vcfeval against the GIAB truth set: a call counts only if its genotype matches. pvc ran
    its full pipeline, with parameters estimated from the data. Measured 2026-10-29 on a Linux
    server limited to 30 cores and 80 GB (Intel Xeon E5-2667 v4): pvc c3d1ed69, GATK 4.6.0.0,
    freebayes 1.3.6, DeepVariant 1.10.0.

[^str]: GIAB HG002 tandem-repeat benchmark v1.0.1, reads downsampled to 30×: of the 6,995 repeats
    where the truth carries a variant, the share whose genotype the caller gets right, compared as
    repeat lengths (HipSTR reports lengths only). Same server as the SNP and indel columns,
    2026-10-30; HipSTR 0.6.2. freebayes was not run on this benchmark.

[^tomato]: Wall time and peak memory of the whole pipeline: 100 tomato accessions at about 3 reads
    a position over 91 Mb of chromosome 1, on the same server, 2026-10-26/27, pvc c3d1ed69.
    pvc's figure is its three steps (pileup files, parameter estimation, calling); run as a
    single step it took 2 h 35 min and 21 GB. Memory is the total of the processes running at once. freebayes did not finish: some
    regions stopped making progress, possibly because they had too many alleles, and the run was
    stopped after 17 hours.

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
