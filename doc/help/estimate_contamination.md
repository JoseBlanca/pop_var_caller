# Estimating cross-sample contamination

`pop_var_caller_exp estimate-contamination` estimates, for each sample in a
panel of Illumina alignments, what share of that sample's sequencing reads
came from a different individual. It is a side tool taken out of a variant
caller under development: the caller fits contamination internally before
calling, and this subcommand exposes the same estimator on its own so it can
be compared against tools such as verifyBamID.

Two commands, run in this order:

```
pop_var_caller_exp repeat-catalog --reference ref.fa        # once per reference
pop_var_caller_exp estimate-contamination ...               # per panel of samples
```

## How it works

The tool picks a fixed set of ordinary positions — positions outside tandem
repeats, the same set in every sample — and walks each alignment once,
recording what that sample showed there. It then fits the whole panel at
once: the sequencing error rate, each position's allele frequency in the
population, and each sample's contamination fraction.

**No genotype is ever called.** At every position the tool sums over both
unknowns — which genotype this sample carries, and which genotype the stray
reads came from — weighting each possibility by how likely the population's
allele frequencies make it and how well it explains the reads actually seen.
A position where the reads are ambiguous therefore contributes its ambiguity
rather than a guess. What contamination looks like under that sum is a small,
one-sided share of reads carrying an allele that no genotype of this sample
explains well, at a position where the population does carry that allele —
where a heterozygote's two alleles would instead be balanced.

Two design choices matter for how you should read the results. First, the
allele frequencies come from the samples in the run itself — no outside
reference panel is used. Second, each sample gets its own expected frequency
at each position: the fit summarises the panel's population structure as a
few axes of variation (principal components, computed from how many copies of
each allele the samples are expected to carry rather than from called
genotypes), and a
sample's frequency is a straight line in its coordinates on those axes. A
genetically diverged sample is therefore not judged against the average of
everyone else.

Those two choices produce the two limits below.

## Two limits to know before comparing methods

### It needs a panel, and about a dozen samples is the floor

Because the frequencies are fitted from the run itself, a sample that
supplies too much of its own fitted frequency would be judged against an
echo of its own reads. The share a sample supplies is roughly
`(axes + 1) / samples`; when it exceeds one half, the tool reports that
sample as *not identified* rather than giving it a number. With the default
of four axes, that means twelve samples or more; below that, expect
refusals. `--components 2` lowers the floor to about seven samples, at the
cost of modelling less population structure.

**A single sample cannot be estimated at all.** verifyBamID can score one
sample against an external panel such as 1000 Genomes; this tool cannot. If
you have one alignment, do not spend time on this method.

### The fraction reads low; the ranking is what is sharp

On simulated data at 30 reads a position — 40 samples in four
subpopulations, one sample contaminated at 3% — the contaminated sample came
back at about 0.026 and every one of the 39 clean samples at 0.0000. So the
separation between contaminated and clean is clean, but the value understates
the truth by something like a tenth to a quarter. When comparing against
another method, compare which samples are flagged and how they rank, not the
absolute fractions.

There is also no fixed threshold to quote. The floor a contaminated sample
must stand above is a property of the panel's own depth, size and diversity,
which is why the output reports the panel's median and highest fraction:
judge each sample against where its own panel sits.

One more property to know when lining results up against another method: the
reported fraction is capped at 0.5. The model reads sequence only, so it
cannot tell a sample 20% contaminated from one that is 80% someone else — a
value above a half would be a mirror image, not a stronger claim.

## Requirements

- A reference FASTA with its `.fai` index beside it (build one with
  `samtools faidx ref.fa`).
- One indexed BAM or CRAM per sample, aligned against that reference. Each
  file must hold exactly one sample; a file whose read groups declare two
  samples is refused, because the walk is per file and two samples would be
  pooled into one set of reads.
- A tandem-repeat catalog, built once from the same reference (step 1).

## Step 1 — build the repeat catalog

```
pop_var_caller_exp repeat-catalog --reference GRCh38.fa
```

This scans the reference once and writes `GRCh38.fa.repeats.parquet` beside
it, where `estimate-contamination` finds it automatically. On the whole of
GRCh38 it takes 103 seconds at the default single thread.

Why does a contamination tool need a repeat catalog? Reads crossing a short
tandem repeat often come back one repeat unit longer or shorter than the DNA
they came from, which imitates contamination; those positions are excluded
from the estimate, and the catalog is what says where they are.

`--threads` speeds the scan up but costs memory in proportion, because each
thread holds the whole chromosome it is scanning: six threads on GRCh38
inside a 16 GB container is killed by the kernel, while the default one
thread builds the same catalog comfortably. The output file is byte-identical
at every thread count.

## Step 2 — run the estimate

```
pop_var_caller_exp estimate-contamination \
    --reference GRCh38.fa \
    --alignment sampleA.cram \
    --alignment sampleB.cram \
    --alignment sampleC.cram \
    --alignment sampleD.cram \
    --threads 8 \
    --output contamination.json
```

One `--alignment` flag per sample, twelve or more of them.

Add `--regions my_regions.bed` to look at a stretch of your own choosing,
and `--components 2` on a panel of fewer than about twelve samples.
`pop_var_caller_exp estimate-contamination --help` lists the remaining
flags.

`--threads` bounds both halves of the run: how many alignments are walked at
once, and how wide the fit runs. It also bounds memory, since each walking
thread holds its own view of the reference.

### Which stretch of genome it looks at

Without `--regions`, the tool draws 100 blocks of 100 kb spread across the
contigs — 10 Mb in total — from a seed that is reported in the output, so
two runs over the same reference look at the same bases. Contigs shorter
than 1 Mb are left out of the draw, because small scaffolds carry mismapped
reads out of proportion to their length. With `--regions`, your BED is used
instead.

Ten megabases is enough. What governs the precision is not bases but the
number of positions where the panel *varies* — positions where the samples
in the run do not all show the same base — and the estimate stops improving
around ten thousand of those. For scale: over 8 Mb of a tomato panel, the
tool found 29,328 varying positions among 16 samples and 36,712 among 63.
The count for your run is in the output as `panel.markers`, and the tool
prints a warning if fewer than 1,000 survive — below that, look at more
genome (`--drawn-bp 30000000`, or a wider `--regions`) rather than trusting
the fractions.

## Runtime

- **Measured**: 16 tomato accessions over 8 Mb, on a laptop — 4 s of setup,
  16 s of walking the alignments (samples walked in parallel), 314 s of
  fitting; about five and a half minutes in total.
- **Estimated**: 20 human samples at 30x over the default 10 Mb, on eight
  cores — roughly three to five minutes of walking, with the fit as the long
  pole; on the order of fifteen minutes in total.

The fit dominates, and its cost grows with the number of samples and the
number of positions recorded in each, not with the size of the genome.

## The output

Progress goes to the terminal; the record is the JSON file named by
`--output`. It has three parts: `run` (what was run, in enough detail to
repeat it), `panel` (the panel as a whole), and `samples` (one entry per
sample, ordered highest fraction first, so the sample that stands out is at
the top).

### `run`

| Field | Meaning |
|---|---|
| `version` | Version of the binary that wrote the file. |
| `reference` | Path of the reference FASTA, as given. |
| `catalog` | Path of the repeat catalog that was used. |
| `regions` | The BED path, or a description of the draw (block count, block length, seed). |
| `analysed_bases` | Total bases in the analysed stretch. |
| `kept_positions` | How many ordinary positions were recorded in every sample. |
| `components` | How many axes of variation the frequencies were fitted with. |
| `seed` | The seed of the block draw and the position choice. Same seed, same reference: same bases. |

### `panel`

| Field | Meaning |
|---|---|
| `samples` | How many samples were given. |
| `estimated` | How many got a number. |
| `not_identified` | How many were refused. |
| `markers` | Positions where the panel varies — what the precision of every fraction rests on. |
| `median` | Median fraction among the estimated samples. |
| `highest` | Highest fraction among the estimated samples. |

`median` and `highest` are there because there is no universal threshold: a
contaminated sample is one that stands well above its own panel's spread.

### `samples`

| Field | Meaning |
|---|---|
| `sample` | The sample name, taken from the alignment file's read groups. |
| `contamination` | The share of this sample's reads that came from another individual (at most 0.5), or `null` when the panel cannot answer for it. Reads low — see the limits above. |
| `not_identified` | When `contamination` is `null`, the reason in plain words; otherwise `null`. |
| `own_frequency_share` | How much of its own fitted allele frequency this sample supplied. A fair share is about `(axes + 1) / samples`; above 0.5 the sample is refused, because its estimate would be a reading of its own noise. `null` for a refused sample. |
| `positions_with_reads` | How many of the kept positions this sample had at least one read at. |

A shortened example — three of a run's sixteen samples:

```json
{
  "run": {
    "version": "0.1.0",
    "reference": "GRCh38.fa",
    "catalog": "GRCh38.fa.repeats.parquet",
    "regions": "100 blocks of 100000 bp drawn at seed 42",
    "analysed_bases": 10000000,
    "kept_positions": 1987432,
    "components": 4,
    "seed": 42
  },
  "panel": {
    "samples": 16,
    "estimated": 15,
    "not_identified": 1,
    "markers": 30514,
    "median": 0.0,
    "highest": 0.0312
  },
  "samples": [
    {
      "sample": "S07",
      "contamination": 0.0312,
      "not_identified": null,
      "own_frequency_share": 0.29,
      "positions_with_reads": 1961210
    },
    {
      "sample": "S02",
      "contamination": 0.0,
      "not_identified": null,
      "own_frequency_share": 0.31,
      "positions_with_reads": 1948773
    },
    {
      "sample": "S11",
      "contamination": null,
      "not_identified": "this sample supplies most of its own fitted allele frequency",
      "own_frequency_share": null,
      "positions_with_reads": 1876031
    }
  ]
}
```

Reading this run: fifteen of sixteen samples were estimated over 30,514
varying positions; the panel's median is 0.0 and one sample, S07, stands at
0.0312 while no other sample reaches 0.01 — that gap, and not the value
itself, is what makes S07 the one to investigate. S11 was refused because it supplied most of its own fitted
frequency; on a sixteen-sample panel that usually means it is genetically
far from the rest of the panel, not that anything went wrong with its file.

## Getting help and reporting trouble

`pop_var_caller_exp estimate-contamination --help` describes every flag.
Errors are written to name the file and the reason — a missing `.fai` names
the FASTA and tells you to run `samtools faidx`, a file with two samples
names the file and the count — so the message itself should say what to fix.
If one does not, that is a bug worth reporting along with the full error
text.
