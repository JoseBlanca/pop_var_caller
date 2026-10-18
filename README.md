# pvc

pvc is a variant caller, focused on speed, memory usage and accuracy.

> **Status: under active development.** pvc gives good results on our benchmarks, but we have
> only started testing it on real projects. Options and the output format may still change. We
> welcome feedback: please open an
> [issue](https://github.com/JoseBlanca/pop_var_caller/issues) with what you ran and what you
> saw.

Supported on Linux and macOS.

## What it aims for

- **Fast and small in memory.** pvc is written in Rust and uses a streaming architecture
  designed for a small memory footprint. For example, it called more than 2,000 whole tomato
  genomes in a week using less than 30 GB.
- **From one sample to thousands, and from 3 reads a position to hundreds.**
- **Accurate calls in model and non-model species alike.** It can take into account
  contamination between samples, per-sample inbreeding coefficients (F) and different error
  priors.
- **STR-aware models.** It uses specialised models for STR and STR-bundle regions.

## How it performs

| what | data | pvc | GATK | freebayes |
|---|---|---|---|---|
| STR genotypes right, of 56,621 at variant sites / false calls at 80,279 reference sites | simulated: 50 samples, 30×, 2,738 STRs, PCR stutter [^str] | 55,978 / 76 | 53,290 / 9,445 | 53,456 / 9,864 |
| SNP recall / precision | human, GIAB HG002 | *to be measured* | *to be measured* | *to be measured* |
| indel recall / precision | human, GIAB HG002 | *to be measured* | *to be measured* | *to be measured* |
| run time and peak memory, whole pipeline | 59 tomato accessions, about 5.5 reads a position, 32 Mb | *to be measured* [^speed] | *to be measured* | *to be measured* |

[^str]: Measured 2026-10-11. In this run pvc was given the true stutter rates rather than
    fitting them. The full grid, one sample and 50, 5× to 60×, PCR and PCR-free, is in
    [doc/devel/reports/implementations/str_edge_deletion_2026-10-11.md](doc/devel/reports/implementations/str_edge_deletion_2026-10-11.md).
    To be rerun on the current version.

[^speed]: What exists so far is pvc alone, on one 18-core Apple laptop (2026-10-17): estimating
    the parameters of the 59 accessions takes 16 minutes with a peak of 410 MB, then calling them
    takes 33 seconds. The time to read the alignments is not included.

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
