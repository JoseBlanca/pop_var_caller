# Options

Most runs need only the options in [Calling](calling.md). This page covers some advanced or
unusual ones. `pvc <command> --help` lists every option with its default.

The calling options apply to both `call-from-alignments` and `call-from-psps` unless the table
says otherwise.

## Material

| Option | Default | Change it when |
|---|---|---|
| `--ploidy` | 2 | your samples are haploid or polyploid. At most 16. `estimate-parameters` fits diploid samples only, so at any other ploidy call with `--defaults` and `--ploidy`, in either mode. |
| `estimate-parameters --inbreeding` | estimated per sample | you know the samples' inbreeding coefficient, for example `0.99` for fully inbred lines. It must be below 1: at 1 no heterozygote could be called, whatever the reads show. It applies the same value to every sample. |
| `estimate-parameters --skip-contamination` | contamination estimated | the samples are known to be clean. |

## Which sites are written

| Option | Default | Change it when |
|---|---|---|
| `--min-site-quality` | 1 | you want every site, including weak ones: `0`. Or you want a stricter set: raise it. The default removes sites whose penalties for looking like sequencing artifacts outweigh their evidence. |
| `--low-base-quality-weight` | 1 | you want the low-base-quality penalty (`BQPEN`) off: `0`. A larger value removes more sites whose variant reads are mostly poor bases, and more true variants with them; at 2, HG002 at 30× lost 29 true calls on chromosome 20 where 1 lost 6. |
| `--hidden-dup-fdr` | 0.01 | you want no hidden-duplication filtering: `0`. The value is the share of removed calls you accept as true variants removed by mistake. |
| `--hidden-dup-tag` | off | you want to see what the hidden-duplication filter removes: the calls are kept and marked `hiddenDup` in FILTER. |
| `--hidden-dup-window` | 50000 | Hidden duplications cluster in some parts of the genome, so one genome-wide rate overstates them in easy ground and removes true variants there. With a width in bases, the filter estimates the rate from the calls within half the window on either side of each call, pulled toward the genome-wide rate where few calls are near. `0` uses the one genome-wide rate. |

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
below the catalogue's own minimum of `5,5,4,4,4,3`. A request outside it is refused. To go lower,
build the catalogue again with a lower `repeat-catalog --min-copies` (next section).

### Satellites, and calling without a catalogue

Two options, on `call-from-alignments` and `generate-psps`, change what happens to repeats. Both
are off by default, and they cannot be given together. A psp records the choice, and
`call-from-psps` takes it from there.

| Option | What it does | Use it when |
|---|---|---|
| `--call-satellites` | calls satellites with the SNP/indel model instead of leaving them uncalled. What it adds is mostly SNPs and short indels inside them: a satellite's own length variation mostly stays uncalled, because reads that disagree about its length make loci wider than `--max-locus-span`. The log says how many satellite bases were handed over and how many loci were too wide. | you need the SNPs inside satellites. |
| `--generic-path-only` | calls every repeat with the SNP/indel model and reads no repeat catalogue. Repeat length variants are then called as indels. It cannot be given with `--catalog`, `--call-satellites` or the options above. | the catalogue cannot be built for your reference, or you are comparing pvc with a caller that has no repeat model. |

## The repeat catalogue (`repeat-catalog`)

| Option | Default | Change it when |
|---|---|---|
| `--threads` | 1 | you have memory to spare: each thread holds a whole chromosome. See [Preparing the inputs](inputs.md#the-repeat-catalogue) for the memory measured. |
| `--output` | `<reference>.repeats.parquet` | the reference's directory is read-only. Then pass `--catalog` to every later command. |
| `--log` | `<catalogue>.log` | you want the log elsewhere. |
| `--min-copies` | `5,5,4,4,4,3` | you want the calling commands to be able to ask for fewer copies than this. |
| `--min-period`, `--max-period` | 1, 6 | you want fewer unit lengths scanned. 6 is the most. |
| `--min-flank-bp` | 15 bases | rarely. A repeat with less sequence than this between it and the chromosome's end is not written down, because no read can anchor it. |
| `--match-reward`, `--mismatch-penalty` | 2, 7 | rarely. They are the scan's score for two bases one unit apart that agree, and its penalty for two that disagree. |

## Estimating the parameters (`estimate-parameters`)

`--inbreeding` and `--skip-contamination` are under [Material](#material).

| Option | Default | Change it when |
|---|---|---|
| `--start-from` | none | samples have been added to a cohort already fitted. Give the earlier parameters file: the SNP/indel fit starts from its numbers instead of searching from scratch, and a new sample or read group starts at the median of the others. Every number is still fitted, and the result is the one a fit from scratch reaches, to within about a tenth of a standard error. On 400 tomato accessions at 48 threads, adding 20 samples took the whole command 43 minutes, against 1 h 20 min from scratch. The file is not used if it holds fewer than half the samples or half the read groups; the log says so. |
| `--threads` | 0, every core | you want to leave cores for other work. The numbers are the same at any count. |
| `--log` | `pvc_call_parameters.log` beside the parameters file | you write several parameters files into one directory. |
| `--ploidy` | 2 | only 2 is estimated; at another ploidy, call with `--defaults`. |

## Reading psp files (`call-from-psps`)

| Option | Default | Change it when |
|---|---|---|
| `--regions` | everything the psps cover | you want to call part of what the psps cover, for example one chromosome per job on a cluster. A region the psps do not cover is refused. See [Splitting a run by region](calling.md#splitting-a-run-by-region). |
| `--psp-prefetch-bytes` | 8 GiB | the psps are on a slow disk. It is how much psp data is read ahead into the operating system's file cache; on a spinning disk with thousands of psps, more is faster. `0` turns it off. It changes no call. |

## Finishing a stopped run (`resume-calling`)

`pvc resume-calling --scratch-dir <directory> --output calls.vcf.gz` takes `--threads`, `--force`,
`--log` and `--keep-tmp-on-error` as a calling run does, and the four imputation memory options
(`--refine-window-bytes`, `--refine-window-floor`, `--refine-core-floor`,
`--refine-margin-density`) and `--phase-min-path-agreement`, in place of the stopped run's. See
[Calling § Output, logs and scratch space](calling.md#output-logs-and-scratch-space).

## Options you should not need

- `--build-index-if-missing`: builds missing indexes. Build them with `samtools index` instead
  (see [Preparing the inputs](inputs.md#the-alignments)).
- The options under *Developer* in `--help`. Most help investigate a call, or trade memory for
  speed without changing any call; `--max-reads-per-position` does change calls.
  `--explain-loci regions.bed out.tsv` is the useful one: for
  every locus in the BED file it writes the candidate alleles, the likelihoods, the genotype and
  qualities, and what happened to the record.
- `pvc dev inspect-psp --psp sample.psp --regions regions.bed` writes, as a table, what one psp
  stored at every locus overlapping the BED: how many reads, and from which read group. It tells
  whether a missing variant's reads reached the psp at all. `--output` writes the table to a file
  instead of standard output.
- `pvc dev type-regions` writes how the reference is divided into repeats and ordinary sequence,
  and `pvc dev estimate-contamination` estimates contamination on its own (see
  [Estimated parameters](parameters.md#contamination-on-its-own)). Both are tools for developing
  pvc.
- The environment variables `NG_SLIPPAGE_REFIT_ROUNDS`, `NG_SLIPPAGE_REFIT_FREE` and
  `NG_TRACT_DISCOVERY` change how STRs are called. They are for development. Do not set them.
