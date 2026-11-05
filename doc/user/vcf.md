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
| `##neighbourEvidence=` | whether the neighbour-evidence pass ran and why, and its settings when it did ([Genotypes from neighbouring sites](neighbour_evidence.md)) |
| `##readBackedPhase=` | at ploidy two where the neighbour-evidence pass did not run: the settings the phase from reads was written with ([Phase from reads](#phase-from-reads)) |

## QUAL and GQ

- **QUAL**: the Phred-scaled probability that no sample in the cohort carries any alternative
  allele at the site, after subtracting three penalties for looking like a sequencing artifact
  (`ABPEN`, `SPPEN` and `BQPEN`, below). A QUAL of 20 means a 1 in 100 chance that the site is not
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
| `BQPEN` | Phred subtracted from QUAL because more of the alternative reads than of the reference reads have a base quality below 20 — the sequencer misreading the same base in many reads, typically late in the read beside a run of G or C |
| `MQREF`, `MQALT`, `MQDIFF` | mean mapping quality of the reads supporting the reference and each alternative allele, and their difference. A strongly negative `MQDIFF` means the alternative reads map worse, a sign of reads from elsewhere in the genome |
| `STR` | flag: the record is an STR or an STR bundle |
| `RU`, `PERIOD` | the repeat unit and its length, one value per repeat in a bundle |
| `PARALOG_LR`, `PARALOG_POST` | the hidden-duplication filter's evidence, and its probability that the site is a hidden duplication |

`ABPEN`, `SPPEN` and `BQPEN` appear only when there were reads to test, and `BQPEN` only when
`--low-base-quality-weight` is not 0. When QUAL is not 0, `QUAL + ABPEN + SPPEN + BQPEN` is the
quality before penalties.

## FORMAT

| Field | Meaning |
|---|---|
| `GT` | genotype; written with `\|` where it is phased, which `PS` says, and with `/` elsewhere |
| `GQ` | genotype quality (above) |
| `DP` | every read of this sample at the locus. `DP` minus the sum of `AD` is the reads that match no written allele: stutter, noise, alleles left out |
| `AD` | reads that match each allele exactly, reference first |
| `REPCN` | STR records: each called allele's length in repeat units, in `GT` order |
| `REPCNB` | STR-bundle records: each called allele's copies in every repeat of the bundle, joined by `_`, for example `20_8_8,21_8_8` |
| `PS` | phase set, at ploidy two: the position of the first record of the block a phased genotype belongs to, `.` where the genotype is not phased |

**A sample whose reads cannot tell the genotypes apart is not called** (`./.`, with `GQ` as `.`),
and this includes samples with no reads. Its `DP` and `AD` are still written if it had reads. A
sample is never called homozygous reference just because it has no evidence. A sample is also
not called when an allele its own reads supported was among those left out at a site with too
many alleles (`--max-generic-alleles`).

## Phase from reads

At ploidy two, pvc writes the phase the reads show. A read pair comes from one molecule, so where it
covers two heterozygous sites of a sample, the alleles it shows there lie on one chromosome copy.
Two heterozygous genotypes are written phased together, `0|1` or `1|0` with one `PS`, when at least
two read pairs that cover both agree on how, at odds of 199 to 1 or more. The odds allow for each
base's quality and for a further 1 error in 100 that no quality shows. Only `PASS` genotypes with
`GQ` 30 or more are phased this way, because a genotype that is wrongly heterozygous gives a
confident link in a random direction. **Genotypes of one sample with the same `PS` are phased
together**; nothing is said about two different `PS` values. A genotype linked to nothing is left
unphased.

Where the [neighbour-evidence pass](neighbour_evidence.md) runs, the reads are combined with the
phase it estimates from the cohort, and its page says how. Elsewhere — one sample, a cohort the
pass does not run on, or `--refine-genotypes off` — the phase comes from the reads alone. Read
pairs reach about a kilobase, so the phase sets are short: they join neighbouring heterozygotes,
not whole genes. STR records are never phased from reads.

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

The header also declares `lowDepth`, but pvc never writes it. A repeat is always genotyped, even
where some reads differ from the reference by part of a repeat unit; those reads are scored by the
stutter model like any others.

**Hidden duplications.** Where the genome holds two copies of a sequence but the reference holds
only one, both copies' reads align to one place. In the samples that carry both copies, a
difference between them looks like a heterozygous variant. At each site pvc weighs two
explanations against each other: a real variant, or a hidden duplication. A duplication predicts
more reads than one copy gives, in the samples that carry it, and an allele balance that matches
the number of copies, never a homozygous alternative genotype. pvc measures each sample's
expected depth for one copy from its own reads. It removes sites that the duplication explains
better, accepting by default that 1 in 100 removed sites was a true variant (`--hidden-dup-fdr`).
If some samples have no depth model, pvc warns on the terminal.
