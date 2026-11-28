# The VCF

pvc writes one VCF (version 4.4) holding SNPs, indels and STRs, in genome order. A site appears
only if at least one sample is called with an alternative allele, and its QUAL reaches
`--min-site-quality` (1 by default). pvc writes no invariant sites and no gVCF: a position
where every sample matches the reference is not in the file, so the file cannot tell such a
position from one with no reads. What the reads showed there is kept in the psp files.

## Header

Besides the usual lines, the header records how the file was made:

| Line | Holds |
|---|---|
| `##source=pvc <version>` | the pvc version |
| `##commandline=` | the command as typed |
| `##reference=` | the path of the reference FASTA the run was called against |
| `##pvcLog=` | the name of the run's log file, which holds every warning and the parameter values used |
| `##pvcParameters=` | `defaults`, or the parameters file's name and MD5 checksum |
| `##contig=` | one line per contig of the reference, with its length and, where pvc has it, its MD5 checksum |
| `##paralogFilter=` | the hidden-duplication filter's settings and results, when it ran (the default) |
| `##neighbourEvidence=` | whether the refinement ran and why, and its settings when it did ([Genotypes from neighbouring sites](neighbour_evidence.md)) |
| `##readBackedPhase=` | at ploidy two where the refinement did not run: the settings the phase from reads was written with ([Phase from reads](#phase-from-reads)) |

## QUAL and GQ

- **QUAL** starts as the Phred-scaled probability that no sample in the cohort carries any
  alternative allele at the site, given all the reads. At that stage a value of 20 means a 1 in
  100 chance that the site is not variable at all. Samples without reads do not change it. pvc
  then subtracts three penalties, each a test of whether the alternative reads look like a
  sequencing or mapping artifact, and writes what is left, at least 0 and at most 9999. **The
  written QUAL is therefore a score, not a probability**: it ranks sites, and it keeps the
  probability reading only where every penalty is 0. When QUAL is not 0, `QUAL + ABPEN + SPPEN +
  BQPEN` is the value before penalties.
- **GQ**: the Phred-scaled probability that this sample's genotype is wrong, at most 99.

### The three penalties

Each penalty compares the reads of the site's main alternative allele, pooled over the samples
with a genotype, with what a real variant would show, and charges the Phred-scaled probability of
a result at least that far from it. Each is 0 when the reads look as expected, and grows with
the number of reads when they do not. STR records carry none.

| Penalty | Applies when | What it compares | Size |
|---|---|---|---|
| `ABPEN`, allele balance | there are fewer alternative reads than the called genotypes predict (each heterozygote half its reads, a homozygote all of them); not where they predict 90% or more | the alternative reads against that share, a two-sided binomial test | the test's Phred, with no limit on the number of reads |
| `SPPEN`, strand and read position | 4 or more alternative reads | the share of alternative reads on the forward strand, and the share placed mostly left of the site, each against the reference reads' share at the same site (beta-binomial tests, the reference share given a prior worth 10 reads) | 4 times the larger of the two Phreds; a quarter of that at 4 alternative reads, rising to the whole at 7 |
| `BQPEN`, low base quality | 4 or more alternative reads, and `--low-base-quality-weight` above 0 | the share of alternative reads with a base quality below 20 against the reference reads' share (a one-sided beta-binomial test, prior worth 30 reads at 2 in 100); only an excess is charged | the test's Phred times `--low-base-quality-weight` (1 by default), with the same ramp from 4 to 7 reads as `SPPEN` |

`SPPEN` is weighted 4 times because QUAL grows by about 30 with each alternative read and the
test's Phred by only a few, so at a weight of one a lopsided artifact would pass at any depth.
[Methods § QUAL and GQ](methods.md#9-qual-and-gq) gives the reasoning.

## INFO

| Field | Meaning |
|---|---|
| `AF` | each alternative allele's frequency as estimated by the model, not counted from the genotypes |
| `AC` | copies of each alternative allele in the called genotypes |
| `AN` | allele copies in the called genotypes; samples not called are left out |
| `DP` | the sum of the samples' `DP` |
| `ABPEN` | Phred subtracted from QUAL because the reads' allele balance looks like an artifact |
| `SPPEN` | Phred subtracted from QUAL because the alternative reads come from one strand or one read position |
| `BQPEN` | Phred subtracted from QUAL because more of the alternative reads than of the reference reads have a base quality below 20 — the sequencer misreading the same base in many reads, typically late in the read beside a run of G or C |
| `MQREF`, `MQALT`, `MQDIFF` | mean mapping quality of the reads supporting the reference and each alternative allele, and their difference. A strongly negative `MQDIFF` means the alternative reads map worse, a sign of reads from elsewhere in the genome |
| `STR` | flag: the record is an STR or an STR bundle |
| `RU`, `PERIOD` | the repeat unit and its length, one value per repeat in a bundle |
| `PARALOG_LR`, `PARALOG_POST` | the hidden-duplication filter's evidence, and its probability that the site is a hidden duplication |

`ABPEN`, `SPPEN` and `BQPEN` appear only when there were reads to test, and `BQPEN` only when
`--low-base-quality-weight` is not 0.

## FORMAT

| Field | Meaning |
|---|---|
| `GT` | genotype; written with `\|` where it is phased, which `PS` says, and with `/` elsewhere |
| `GQ` | genotype quality (above) |
| `DP` | every read of this sample at the locus. `DP` minus the sum of `AD` is the reads that match no written allele: stutter, noise, alleles left out |
| `AD` | reads that match each allele exactly, reference first |
| `PL` | genotype likelihoods: for each genotype, in the VCF's genotype order, the probability of the sample's reads under that genotype, Phred-scaled so that the most likely genotype is 0. No prior is in it, but where contamination was estimated the contamination term is ([Methods § 4](methods.md#4-region-types-and-their-read-models)) |
| `GP` | genotype posteriors, in the same order, to three decimals: the likelihoods times the genotype prior, which is built from the other samples' genotypes and the sample's inbreeding coefficient ([Methods § 8](methods.md#8-calling)). Where the refinement ran, its posterior, and the header's `GP` line says so |
| `REPCN` | STR records: each called allele's length in repeat units, in `GT` order |
| `REPCNB` | STR-bundle records: each called allele's copies in every repeat of the bundle, joined by `_`, for example `20_8_8,21_8_8` |
| `PS` | phase set, at ploidy two: the position of the first record of the block a phased genotype belongs to, `.` where the genotype is not phased |

**A sample whose reads cannot tell the genotypes apart is not called** (`./.`, with `GQ` as `.`),
and this includes samples with no reads. Its `DP` and `AD` are still written if it had reads; its
`PL` is all 0 and its `GP` is the prior alone. A sample is never called homozygous reference just
because it has no evidence. A sample is also not called when an allele its own reads supported
was among those left out at a site with too many alleles (`--max-generic-alleles`); its `PL` and
`GP` are then `.`.

## Phase from reads

At ploidy two, pvc writes the phase the reads show. A read pair comes from one molecule, so where it
covers two heterozygous sites of a sample, the alleles it shows there lie on one chromosome copy.
Two heterozygous genotypes are written phased together, `0|1` or `1|0` with one `PS`, when at least
two read pairs that cover both agree on how, at odds of 199 to 1 or more. The odds allow for each
base's quality, and never take a base's error to be below 1 in 100, whatever its quality says.
Only `PASS` genotypes with `GQ` 30 or more are phased this way, because a genotype that is wrongly heterozygous gives a
confident link in a random direction. **Genotypes of one sample with the same `PS` are phased
together**; nothing is said about two different `PS` values. A genotype linked to nothing is left
unphased.

Where the [refinement](neighbour_evidence.md) runs, the reads are combined with the phase it
estimates from the cohort, and its page says how. Elsewhere — one sample, a cohort the
refinement does not run on, or `--refine-genotypes off`, the default — the phase comes from the
reads alone. Read pairs reach about a kilobase, so the phase sets are short: they join
neighbouring heterozygotes, not whole genes. STR records are never phased from reads.

How much is phased depends on depth and on how close the heterozygotes are. Measured on GIAB's
HG002 in 100 random regions (452 kb), the phase scored against the one its parents' genotypes
give:

| depth | heterozygous SNPs phased | links scored | of them wrong |
|---|---|---|---|
| 30× | 197 of 289 (68%) | 106 | 0 |
| 10× | 106 of 289 (37%) | 52 | 0 |
| 3× | 20 of 289 (7%) | 10 | 0 |

Above 30× little more is phased: each sample keeps 32 reads at a site for this. These counts are
too small to show an error rate; on simulated single samples, which have no mapping errors, none of
1,893 links at 10× and 30× was wrong. **Reads from a hidden duplication all err the same way**,
which neither measurement tests: the hidden-duplication filter is what keeps such sites out of
`PASS`, and with it off (`--hidden-dup-fdr 0`) the log warns that the phase can then be confidently
wrong.

## STR records

- An STR allele is written as the whole repeat sequence, and the reference allele is the whole
  reference repeat. There are no symbolic alleles.
- All the alleles at one repeat are in one record.
- Select STR records with `bcftools view -i 'STR=1'`, or exclude them with `-e 'STR=1'`.
- If an allele would be empty (a repeat deleted entirely), every allele gets the reference base
  to its left, and `POS` moves one base left. At a chromosome's first base, the base to the right
  is added instead.
- A deletion that runs from the repeat into the sequence beside it widens the STR record over
  those bases. A SNP or indel that overlaps the widened record is merged into it: the record
  covers both, and its alleles are the whole sequences each copy carries, repeat and SNP or indel
  together. Such a record can hold more than 32 alleles. pvc never writes a `*` allele;
  [Overlapping variants](overlapping_variants.md) explains why this matters.
- In a merged record, which SNP or indel allele goes with which repeat allele is decided from each
  sample's reads where it can be. Where a sample is heterozygous at both and its reads cannot
  decide, the pairing seen most often in the other samples is used, and `GQ` falls to reflect how
  sure that pairing is.

## FILTER

| Value | Meaning |
|---|---|
| `PASS` | no filter applied |
| `EMNoConv` | the model's allele-frequency estimate did not settle in the allowed rounds; the call is written but less certain |
| `tooManyAlleles` | STR record with more candidate alleles than the limit of 32, reference included; only the best-supported were considered |
| `hiddenDup` | written only with `--hidden-dup-tag`: the site looks like a hidden duplication, explained below |

A repeat is always genotyped, even where some reads differ from the reference by part of a repeat
unit; those reads are scored by the stutter model like any others.

**Hidden duplications.** Where the genome holds two copies of a sequence but the reference holds
only one, both copies' reads align to one place. In the samples that carry both copies, a
difference between them looks like a heterozygous variant. At each site pvc weighs two
explanations against each other: a real variant, or a hidden duplication. A duplication predicts
more reads than one copy gives, in the samples that carry it, and an allele balance that matches
the number of copies, never a homozygous alternative genotype. pvc measures each sample's
expected depth for one copy from its own reads. It removes sites that the duplication explains
better, accepting by default that 1 in 100 removed sites was a true variant (`--hidden-dup-fdr`).
If some samples have no depth model, pvc warns on the terminal.
