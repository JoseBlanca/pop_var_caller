# The VCF

pvc writes one VCF (version 4.4) holding SNPs, indels and STRs, in genome order. A site appears
only if at least one sample is called with an alternative allele, and its QUAL reaches
`--min-site-quality` (1 by default).

## Header

Besides the usual lines, the header records how the file was made:

| Line | Holds |
|---|---|
| `##source=pvc <version>` | the pvc version |
| `##commandline=` | the command as typed |
| `##pvcLog=` | the name of the run's log file, which holds every warning and the parameter values used |
| `##pvcParameters=` | `defaults`, or the parameters file's name and MD5 checksum |
| `##paralogFilter=` | the hidden-duplication filter's settings and results, when it ran (the default) |

## QUAL and GQ

- **QUAL**: the Phred-scaled probability that no sample in the cohort carries any alternative
  allele at the site, after subtracting two penalties for looking like a sequencing artifact
  (`ABPEN` and `SPPEN`, below). A QUAL of 20 means a 1 in 100 chance that the site is not
  variable at all. It is at least 0 and at most 9999. QUAL does not grow with the number of
  samples that match the reference.
- **GQ**: the Phred-scaled probability that this sample's genotype is wrong, at most 99.

## INFO

| Field | Meaning |
|---|---|
| `AF` | each alternative allele's frequency as estimated by the model, not counted from the genotypes |
| `AC` | copies of each alternative allele in the called genotypes |
| `AN` | allele copies in the called genotypes; samples not called are left out |
| `DP` | the sum of the samples' `DP` |
| `ABPEN` | Phred subtracted from QUAL because the reads' allele balance looks like an artifact |
| `SPPEN` | Phred subtracted from QUAL because the alternative reads come from one strand or one read position |
| `MQREF`, `MQALT`, `MQDIFF` | mean mapping quality of the reads supporting the reference and each alternative allele, and their difference. A strongly negative `MQDIFF` means the alternative reads map worse, a sign of reads from elsewhere in the genome |
| `STR` | flag: the record is an STR or an STR bundle |
| `RU`, `PERIOD` | the repeat unit and its length, one value per repeat in a bundle |
| `PARALOG_LR`, `PARALOG_POST` | the hidden-duplication filter's evidence, and its probability that the site is a hidden duplication |

`ABPEN` and `SPPEN` appear only when there were reads to test. When QUAL is not 0,
`QUAL + ABPEN + SPPEN` is the quality before penalties.

## FORMAT

| Field | Meaning |
|---|---|
| `GT` | genotype, unphased |
| `GQ` | genotype quality (above) |
| `DP` | every read of this sample at the locus. `DP` minus the sum of `AD` is the reads that match no written allele: stutter, noise, alleles left out |
| `AD` | reads that match each allele exactly, reference first |
| `REPCN` | STR records: each called allele's length in repeat units, in `GT` order |
| `REPCNB` | STR-bundle records: each called allele's copies in every repeat of the bundle, joined by `_`, for example `20_8_8,21_8_8` |

**A sample whose reads cannot tell the genotypes apart is not called** (`./.`, with `GQ` as `.`),
and this includes samples with no reads. Its `DP` and `AD` are still written if it had reads. A
sample is never called homozygous reference just because it has no evidence. A sample is also
not called when an allele its own reads supported was among those left out at a site with too
many alleles (`--max-generic-alleles`).

## STR records

- An STR allele is written as the whole repeat sequence, and the reference allele is the whole
  reference repeat. There are no symbolic alleles.
- All the alleles at one repeat are in one record.
- Select STR records with `bcftools view -i 'STR=1'`, or exclude them with `-e 'STR=1'`.
- If an allele would be empty (a repeat deleted entirely), every allele gets the reference base
  to its left, and `POS` moves one base left. At a chromosome's first base, the base to the right
  is added instead.
- A deletion that runs from the repeat into the sequence beside it widens the STR record over
  those bases. A SNP or indel record that overlaps it marks the overlapping allele as `*`.

## FILTER

| Value | Meaning |
|---|---|
| `PASS` | no filter applied |
| `EMNoConv` | the model's allele-frequency estimate did not settle in the allowed rounds; the call is written but less certain |
| `tooManyAlleles` | STR record with more candidate alleles than the limit of 32, reference included; only the best-supported were considered |
| `hiddenDup` | written only with `--hidden-dup-tag`: the site looks like a hidden duplication, explained below |

The header also declares `notPeriodic` and `lowDepth`, but pvc never writes them. A repeat whose
allele lengths do not step by its repeat unit is called over the reference allele alone, so it is
not written; the log counts these repeats.

**Hidden duplications.** Where the genome holds two copies of a sequence but the reference holds
only one, both copies' reads align to one place. In the samples that carry both copies, a
difference between them looks like a heterozygous variant. At each site pvc weighs two
explanations against each other: a real variant, or a hidden duplication. A duplication predicts
more reads than one copy gives, in the samples that carry it, and an allele balance that matches
the number of copies, never a homozygous alternative genotype. pvc measures each sample's
expected depth for one copy from its own reads. It removes sites that the duplication explains
better, accepting by default that 1 in 100 removed sites was a true variant (`--hidden-dup-fdr`).
If some samples have no depth model, pvc warns on the terminal.
