# pvc — a population variant caller

pvc calls SNPs, short indels and short tandem repeats (STRs) from aligned reads, for one sample
or for a cohort, and writes one VCF. It does the job freebayes or GATK HaplotypeCaller do. It is
one caller whose models are STR-aware: inside a repeat it calls whole repeat alleles and expects
reads to be off by a repeat unit or two, and everywhere else it calls ordinary SNPs and indels.

> **Status: under active development.** pvc gives good results on our benchmarks, but we have
> only started testing it on real projects. Options and the output format may still change. We
> welcome feedback: please open an
> [issue](https://github.com/JoseBlanca/pop_var_caller/issues) with what you ran and what you
> saw.

Supported on Linux and macOS.

## What it aims for

- **Accurate calls in model and non-model species alike.** A caller has to know some numbers
  about the data before it can call: how often the sequencer is wrong, how diverse the
  population is, how inbred each sample is. Most callers take these as built-in constants tuned
  on human data. pvc estimates them from your own reads, so an inbred crop, a wild population or
  a human trio each get their own values, with no truth set needed.
- **STRs called as repeats.** When DNA is copied inside a repeat, the copy often gains or loses a
  unit, so a sample with 12 copies shows some reads with 11 or 13. This is called stutter. A
  caller that does not expect it reports those reads as new indel alleles. pvc treats each repeat
  as one locus whose alleles are whole repeat sequences, and it estimates the stutter rates from
  your reads, for each class of repeat. Repeats close enough together that one read spans several are called
  together, as one locus.
- **From one sample to thousands, from 3 reads a position to hundreds.** A single low-coverage
  sample is the hardest case and will be the weakest. Where a method needs a cohort, pvc says
  what it does with one sample instead of returning a confident wrong number.
- **Fast and small in memory.** pvc is written in Rust and streams through the genome. Calling
  memory is designed to stay about flat as samples are added. So far that has been measured from
  1 to 63 samples.

## How it performs

Every figure below names the data it was measured on. A figure from one corner of the range
(one species, one depth, simulated or real) describes that corner only.

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

- **Calls** SNPs, indels up to 100 bp, and STRs with repeat units of 1 to 6 bases spanning up to
  100 bp. Each sample's genotype comes with a genotype quality (GQ), and each site with a QUAL.
- **Estimates the model's parameters from the data** before calling: base-quality calibration
  and contamination for each read group, the inbreeding coefficient F of each sample,
  the population's diversity, and stutter for each kind of repeat. They are written to a
  parameters file that you can read, edit and reuse.
- **Estimates contamination**: what share of each read group's reads came from another individual.
  Calling then takes it into account.
- **Flags hidden duplications**: places where two copies of a sequence exist in the genome but
  only one in the reference, so both copies' reads pile onto one place and look like
  heterozygous variants. Common in plant genomes. On by default.
- **Lets you add samples to a cohort later** without re-reading the alignments of the samples you
  already have.

It does not call structural variants, satellite repeats longer than 100 bp, or deletions longer
than 100 bp, and it does not write genotype likelihoods (PL) yet.

## Install

pvc is installed from source for now. You need:

- Rust 1.98 or newer, from [rustup](https://rustup.rs);
- a C compiler (`cc` on Linux, the Xcode command-line tools on macOS), for two compression
  libraries pvc builds.

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

**psp mode: three steps, with parameters estimated from your data.** First summarise each sample's reads
into one file (a *psp* file), then estimate the parameters from the cohort, then call. A sample
added later needs only its own psp.

```
pvc generate-psps --reference genome.fa --alignment s1.bam --output-dir psps/   # once per sample
pvc estimate-parameters --reference genome.fa --psp psps/ --output cohort.parameters.toml
pvc call-from-psps --reference genome.fa --psp psps/ \
    --parameters cohort.parameters.toml --output calls.vcf.gz
```

`pvc <command> --help` describes every option. The [user guide](doc/user/calling.md) says which
mode to choose: for inbred or structured material, use psp mode, because direct mode assumes
every sample is outbred.

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
