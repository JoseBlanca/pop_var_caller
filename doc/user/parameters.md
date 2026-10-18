# Estimated parameters

Before it can call, a caller has to know some numbers about the data: how often the sequencer is
wrong, how diverse the population is, how inbred each sample is, how often a repeat stutters. pvc
estimates them from your reads with `estimate-parameters` and writes them to a parameters file,
a TOML text file you can read and edit.

## What is estimated

| Quantity | For each | What it changes |
|---|---|---|
| Base-quality calibration: a correction to the error rates the sequencer reports | read group | how much a base that disagrees with the genotype counts against it |
| Contamination: the share of reads that came from another individual | read group (a read group with too little data gets its sample's value) | a few reads of another allele are less likely to be called a heterozygote |
| Inbreeding coefficient F | sample | how likely a homozygote is compared with a heterozygote |
| Population diversity | run | how likely a site is to be variable before the reads are seen |
| Stutter rates | class of repeat (unit length × number of copies) | how many reads a unit short or long a genotype is expected to produce |
| Spread of allele lengths between repeats | class of repeat | how likely a repeat is to carry an allele of a new length |
| Substitution rate inside repeats | read group × class of repeat | how much a base difference inside a repeat counts |

**Inbreeding needs a cohort.** One genome alone cannot show whether it carries too few
heterozygotes, so with one sample F is set to 0 and marked as such in the file. If you know how
your material was bred, state it with `estimate-parameters --inbreeding F`.

**So does the strength of the spread of allele lengths between repeats.** With one sample, pvc
estimates which lengths each class of repeat carries, but not how much that varies from one repeat
to the next; it holds those lengths at a between-repeat share of 1/2, marked `defaulted`. A class
with no estimate at all gets an even spread at 1/11, as with `--defaults`.

## The built-in values (`--defaults`)

| Quantity | Built-in value |
|---|---|
| base-quality calibration | none: the sequencer's qualities are taken as they are |
| contamination | none |
| inbreeding F | 0, for every sample: outbred, Hardy–Weinberg proportions |
| stutter | HipSTR's built-in rates: 1 read in 20 a unit short and 1 in 20 a unit long, for every repeat |
| spread of allele lengths between repeats | a between-repeat share of 1/11 |
| substitution rate inside repeats | 1 base in 1,000 |
| weight given to a read that fits no allele at a repeat | 0.20 |

The inbreeding value is wrong in a known direction for inbred or structured material: it makes pvc
call too many heterozygotes. Use psp mode for such material.

## Reading the file

Five kinds of number carry a `warrant` that says where the value came from: the base-quality
calibration, the inbreeding coefficient, the spread of allele lengths between repeats, the
substitution rate inside repeats, and the weight of a read that fits no allele at a repeat.

| `warrant` | Meaning |
|---|---|
| `fitted_here` | estimated from this cohort's reads |
| `borrowed` | this read group had too little data, so the mean of the same sample's other read groups was taken |
| `supplied` | given by you, for example with `--inbreeding` |
| `defaulted` | a built-in value, because nothing could be estimated |

The other sections say where their numbers came from in their own entries.

Contamination appears only for read groups where it was estimated. A read group with no row is
treated as uncontaminated.

## Changing a value

Edit the file and pass it to the calling step with `--parameters`. To change one built-in value in
direct mode, run once with `--defaults`. The log of that run contains the values it used,
between the lines `# ---- parameters file begins ----` and `# ---- parameters file ends ----`.
That block is itself a valid parameters file: copy it out, edit it and pass it back with
`--parameters`.

## What a parameters file is bound to

A parameters file was estimated on particular data, and a calling run checks that it is being
used on the same data:

| If the calling run differs in… | the run… |
|---|---|
| the reference | refuses |
| the set of samples (order does not matter) | refuses |
| the read groups, or a read group's library name | refuses |
| `--ploidy` | refuses; leave `--ploidy` out and the file's ploidy is used |
| what counts as a repeat | warns and continues |

To call a different set of samples, estimate the parameters again on that set.

## Contamination on its own

`pvc dev estimate-contamination` runs the same contamination estimate without the rest of the
fit, straight from alignments, and writes the result as JSON. It exists to compare pvc's estimate
with other tools such as verifyBamID. See [Estimating contamination](../help/estimate_contamination.md).
