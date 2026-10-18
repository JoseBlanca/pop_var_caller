# Options

Most runs need only the options in [Calling](calling.md). This page covers the ones you may need
for unusual material, and says when to change each. `pvc <command> --help` lists every option
with its default.

The calling options apply to both `call-from-alignments` and `call-from-psps` unless the table
says otherwise.

## Material

| Option | Default | Change it when |
|---|---|---|
| `--ploidy` | 2 | your samples are haploid or polyploid. At most 16. In psp mode set it on `estimate-parameters`; the calling step takes it from the parameters file. |
| `estimate-parameters --inbreeding` | estimated per sample | you know how the material was bred, for example `1` for fully inbred lines. It applies the same value to every sample. |
| `estimate-parameters --skip-contamination` | contamination estimated | the samples are known to be clean. |

## Which sites are written

| Option | Default | Change it when |
|---|---|---|
| `--min-site-quality` | 1 | you want every site, including weak ones: `0`. Or you want a stricter set: raise it. The default removes sites whose penalties for looking like sequencing artifacts outweigh their evidence. |
| `--hidden-dup-fdr` | 0.01 | you want no hidden-duplication filtering: `0`. The value is the share of removed calls you accept as true variants removed by mistake. |
| `--hidden-dup-tag` | off | you want to see what the hidden-duplication filter removes: the calls are kept and marked `hiddenDup` in FILTER. |

## Depth and locus size

| Option | Default | Change it when |
|---|---|---|
| `--max-depth` | 1,000 | your data are deeper than 1,000 reads a position in one read group. A locus where any read group of any sample is deeper is not called at all, because depth that high usually means several copies of a sequence collapsed onto one place in the reference. It also caps the reads kept per read group at a position. In psp mode set it on `generate-psps` too. |
| `--max-locus-span` | 100 bases | you need deletions longer than 100 bases and your reads are long enough to span them. This is the widest locus pvc builds, so the longest deletion it can call. |
| `--max-generic-alleles` | 6 | sites with more than 5 alternative alleles matter to you. It counts the reference allele, and the worst-supported extra alleles are left out. STR loci have their own limit of 32. |

## What counts as a repeat

These decide which stretches of the genome are called as STRs and which as ordinary sequence.
They are set on `call-from-alignments` and `generate-psps`. A psp records them, and
`call-from-psps` takes them from there. **Most users should not change them.**

| Option | Default | Meaning |
|---|---|---|
| `--min-copies` | `8,6,6,6,5,4` | fewest copies of the unit for a stretch to count as a repeat, one value for each unit length from 1 to 6 bases. Below it, a stretch is called as ordinary sequence. |
| `--min-period`, `--max-period` | 1, 6 | shortest and longest repeat unit, in bases. A minimum of 1 includes homopolymers. |
| `--max-str-len` | 100 bases | a longer repeat is a satellite and is not called. With 150-base reads, a read spans a repeat and enough sequence on both sides only up to about 90 bases. |
| `--max-str-bundle-len` | 100 bases | a longer STR bundle is a satellite and is not called. `0` calls no bundle: every STR bundle is treated as a satellite and skipped. |
| `--min-purity` | 0.8 | how much of a repeat must match a perfect run of its unit, from 0 to 1. |

They can only ask for less than the repeat catalogue holds. For example, `--min-copies` cannot go
below the catalogue's own minimum of `5,5,4,4,4,3`. A request outside it is refused.

## Options you should not need

- `--build-index-if-missing`: builds missing indexes. Build them with `samtools index` instead
  (see [Preparing the inputs](inputs.md#the-alignments)).
- The options under *Developer* in `--help`. Most help investigate a call, or trade memory for
  speed without changing any call; `--max-reads-per-position` does change calls.
  `--explain-loci regions.bed out.tsv` is the useful one: for
  every locus in the BED file it writes the candidate alleles, the likelihoods, the genotype and
  qualities, and what happened to the record.
- The environment variables `NG_SLIPPAGE_REFIT_ROUNDS`, `NG_SLIPPAGE_REFIT_FREE` and
  `NG_TRACT_DISCOVERY` change how STRs are called. They are for development. Do not set them.
