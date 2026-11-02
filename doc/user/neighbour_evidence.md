# Genotypes from neighbouring sites and related samples

In a cohort of inbred lines, each sample's chromosomes are, stretch by stretch, copies of the
haplotypes other samples carry. pvc uses that: after calling, it re-calls every sample's genotype
from its own reads together with what the surrounding sites say about which other samples it
resembles there. This is the **neighbour-evidence pass**. It does two things:

- it corrects genotypes that a few reads got wrong, and
- it fills many `./.` — sites where a sample had no reads — with the genotype the samples it
  resembles carry there, when there is enough evidence to do so.

It needs no extra input: it reads the genotype likelihoods (`PL`) and posteriors (`GP`) of the
calls and each sample's inbreeding coefficient from the parameters file.

## When it runs

`--neighbour-evidence` on `pvc call-from-psps` and `pvc call-from-alignments`:

| value | what it does |
|---|---|
| `auto` (default) | runs on a diploid cohort of at least 10 samples whose median inbreeding coefficient is at least 0.5 |
| `on` | runs on any diploid cohort of two or more samples |
| `off` | does not run |

The VCF header states whether it ran and why, in a `##neighbourEvidence=` line.

`pvc refine-genotypes --vcf calls.vcf.gz --parameters params.toml --output refined.vcf.gz` runs
the pass on a VCF pvc wrote without it, so the same calls can be compared with and without it.

**Forcing it on outside what `auto` accepts is at your own risk.** The pass is built for inbred
lines. It does not phase haplotypes, so in an outbred sample — heterozygous over most of its
genome — the samples it copies from give it only a blurred guide. On simulated outbred cohorts it
removed 13% to 24% of the wrong genotypes at 50 samples and did harm at 10: three quarters more wrong genotypes where samples
had reads, in a cohort of 10 drawn from 64 lineages. Outbred species such as human, cattle or
outbred crops have not been measured and are expected to fare no better. For those, use an
imputation tool that phases, such as Beagle or GLIMPSE.

## What it changes in the VCF

- `GT`, `GQ` and `GP` of every sample at every PASS record, from the new posterior. `PL` and `AD`
  stay as the reads gave them, cut to the alleles kept: a filled genotype keeps a flat `PL`, which
  shows that none of its evidence came from its own reads.
- `AC` and `AN` are recounted. An alternative allele no sample carries any more is removed from the
  record, and a record left with no alternative carried is removed.

**A `./.` stays `./.`**, with the pass's posterior as its `GP`, when:

- the sample has no reads at sites on both sides within about one copied stretch (30 to 300 kb);
- the filled genotype's `GQ` would be below 10; or
- the sample is probably copying no other sample there — its haplotype is not carried by anyone
  else in the cohort — so the genotype would only be the cohort's allele frequency.

## How good the fills are: the fill check

Every run reports, for each sample, how often its fills would be wrong, measured on its own data.
At the sites where the sample's own evidence settles its genotype, pvc asks what the pass would
have written had the sample had no reads there, and compares that with the genotype its reads
gave. It is the usual test of an imputation — hide genotypes you know and see whether they come
back — done at every such site.

The result is written beside the VCF as `<name>.fill_check.tsv`:

| column | meaning |
|---|---|
| `sample` | the sample |
| `filled` | the `./.` the pass filled for this sample |
| `sites_checked` | the sites where its genotype was hidden and predicted |
| `fills_wrong_there` | how many of those predictions were wrong |
| `wrong_share` | `fills_wrong_there / sites_checked`, given only from 30 checked sites |

The log gives the median share and the samples with the highest.

**Read it as a best case.** The check runs at sites where the sample had reads; a site no read
reached can be harder, for example inside a repeat. On simulated cohorts at 3× the cohort's
estimated share was close to the true fill error rate: 1.1% against 1.5%, 0.13% against 0.15%; on
outbred cohorts it overstated it, 4.2% against 2.9%.

**At about one read a site there is no check.** Almost no genotype is settled at that depth, so
`wrong_share` is `.` for most samples and the log says no rate was estimated. Those fills are still
written; at 1× they were wrong 2 to 2.5 times in a hundred on simulated inbred cohorts.

## Warnings about poor fills

At the end of a run pvc warns:

- when the fill check, pooled over all samples, is wrong at more than 2% of the sites checked —
  usually a sign that the cohort is outbred, or that its samples share few haplotypes; and
- for each sample whose own rate is above 2% by more than the chance of its few checked sites
  explains. A sample with one wrong site in 30 is not named.

A warning does not stop the run. Look at the named samples — a mislabelled sample, a contaminated
library, or a sample from a population the rest of the cohort does not represent can each cause
it — and consider `--neighbour-evidence off`, or `pvc refine-genotypes` on a cohort without the
samples in question.
