# Genotypes from neighbouring sites and related samples

In a cohort of inbred lines, each sample's chromosomes are, stretch by stretch, copies of the
haplotypes other samples carry. pvc uses that: after calling, it re-calls every sample's genotype
from its own reads together with what the surrounding sites say about which other samples it
resembles there. This is the **refinement**, what is often called imputation. It does two
things:

- it corrects genotypes that a few reads got wrong, and
- it fills many `./.` — sites where a sample had no reads — with the genotype the samples it
  resembles carry there, when there is enough evidence to do so.

It needs no extra input: it reads the genotype likelihoods (`PL`) and posteriors (`GP`) of the
calls and each sample's inbreeding coefficient from the parameters file.

## When it runs

`--refine-genotypes` on `pvc call-from-psps` and `pvc call-from-alignments`:

| value | what it does |
|---|---|
| `off` (default) | does not run |
| `auto` | runs on a diploid cohort of at least 10 samples whose median inbreeding coefficient is at least 0.5 |
| `on` | runs on any diploid cohort of two or more samples |

The VCF header states whether it ran and why, in a `##neighbourEvidence=` line. **Where it does
not run, on a diploid cohort, the phase is written from the reads alone**: `GT` becomes `a|b` and a
`PS` field is added wherever two good fragments link two heterozygotes of a sample, and the header
gains a `##readBackedPhase=` line stating the thresholds (see [the VCF page](vcf.md#phase-from-reads)).
`off`, the default, therefore turns refinement off, not phase.

The refinement runs only inside a calling run. To compare the same cohort with and without it,
call it twice, with `--refine-genotypes on` and `off`.

**Forcing it on outside what `auto` accepts is at your own risk.** The refinement was built for
inbred lines. In an outbred sample — heterozygous over most of its genome — few other samples
carry the same pair of haplotypes over a long stretch, and the samples it copies from give it a
weaker guide. In outbred cohorts it has been measured only on simulated data, and the results
differ with the cohort:

| simulated outbred cohort | pvc version | refinement off → on |
|---|---|---|
| 400 samples from a coalescent simulation (msprime), 100 each at 30×, 10×, 5× and 2× | 2026-11-08 | SNP genotypes missed at 2×: 12,220 → 1,546; fewer at every depth below 30× ([Accuracy](accuracy.md#imputation-genotype-refinement-and-filling)) |
| 40 samples from the same simulation, all at 2× | 2026-11-07 | genotypes right 29,381 → 31,378, but heterozygotes right 7,521 → 7,042 |
| 50 samples at 3×, mosaics of 16 founder haplotypes | 2026-11-28 | wrong genotypes where samples had reads: 4,959 → 1,814 (63% fewer) |
| 50 samples at 3×, mosaics of 64 founder haplotypes | 2026-11-28 | the same: 5,031 → 4,149 (18% fewer) |
| 10 samples at 3×, mosaics of 64 founder haplotypes | 2026-11-28 | the same: 973 → 935 (4% fewer) |

The cohorts differ in more than their size. The 400-sample cohort holds 100 samples at 30×,
whose haplotypes are good templates; in the others every sample is shallow. And with 64
founders and 10 samples, a sample's haplotype at a given place is often carried by no other
sample, so there is little to copy. The code differs too: the first two rows were measured
before a change on 2026-11-28 that stopped the refinement from making small cohorts at 3× worse
(with 10 samples from 64 founders, it had raised the wrong genotypes from 973 to 1,088), and
they have not been run again. Cohorts of real outbred species, such as human, cattle or outbred
crops, have not been measured at all. For those, an imputation tool built for outbred
populations, such as Beagle or GLIMPSE ([Methods § References](methods.md#references)), is the
safer choice.

## What it changes in the VCF

- `GT`, `GQ` and `GP` of every sample at every PASS record, from the new posterior. `PL` and `AD`
  stay as the reads gave them, cut to the alleles kept: a filled genotype keeps a flat `PL`, which
  shows that none of its evidence came from its own reads.
- A heterozygote whose phase the refinement is confident of is written with `|` and a `PS` field
  (next section); `REPCN` follows the written order of `GT`.
- `AC` and `AN` are recounted. An alternative allele no sample carries any more is removed from the
  record, and a record left with no alternative carried is removed.

**A `./.` stays `./.`**, with the refinement's posterior as its `GP`, when:

- the sample has no reads at sites on both sides within about one copied stretch (30 to 300 kb);
- the filled genotype's `GQ` would be below 9; or
- the chance that the sample is copying no other sample there is above 1 in 10 — its haplotype
  may not be carried by anyone else in the cohort — so the genotype would only be the cohort's
  allele frequency.

## Phase: `|` and `PS`

The refinement estimates each sample's two haplotypes, and it writes the phase where it is
confident. A heterozygous genotype is then written `0|1` or `1|0`, the first allele on one
chromosome copy and the second on the other, with a `PS` (*phase set*) field. **Genotypes of one sample with the same
`PS` are phased together**: their first alleles lie on one chromosome. Nothing is said about two
different `PS` values. `PS` is the position of the first phased genotype of its set, the VCF
specification's convention.

Every other genotype is written as before, `0/1`, `0/0`, `./.`, with `PS` `.`. A heterozygote is
left unphased when the refinement is not sure which copy carries which allele, when it is not
sure how it relates to its neighbours, or when it would be a phase set of one.

Phase sets are short and many by design. A link between two neighbouring heterozygotes is written
only when at least 29 of the refinement's 36 drawn estimates of the haplotypes agree on it. On
simulated cohorts at 3× (some samples at 30×), mosaics of 16 or 64 founder haplotypes
(*lineages*), run with pvc of 2026-11-02, counting SNPs only:

| simulated cohort | true heterozygotes phased | links inside phase sets that are wrong |
|---|---|---|
| inbred or partly inbred, 50 to 200 samples, 16 or 64 lineages | 44% to 77% | 0.05% to 0.5% |
| 50 outbred samples, 16 lineages | 24% | 0.4% |
| 50 outbred samples, 64 lineages | 13% | 4.7% |
| 10 outbred samples, 64 lineages | 1% | 7 of 13 |

**Outbred cohorts with many distinct haplotypes are where the phase is least reliable**, as in the
last two rows. The refinement phases few genotypes there, but those it does phase are wrong far
more often than elsewhere. These cohorts were not run again after the changes of 2026-11-06 and
2026-11-07. In the 400-sample outbred population of [Accuracy](accuracy.md#imputation-genotype-refinement-and-filling),
with 100 samples at 30×, pvc of 2026-11-08 phased 121,492 of 145,951 heterozygotes right, and
7 of its 94,860 links were wrong.

`--phase-min-path-agreement` (default `0.8`) is the share of the drawn estimates that must agree for
a link to be written. Lower values give longer phase sets with more wrong links.

### Reads in the phase

The refinement also uses each sample's reads, which show directly which alleles lie on one molecule
([Phase from reads](vcf.md#phase-from-reads)). Between two neighbouring heterozygotes it phases:

- a link the drawn estimates agree on is kept unless the sample's reads contradict it at odds of 9
  to 1 or more, in which case the phase set ends there;
- where the estimates do not agree, reads that link the two at odds of 199 to 1, both genotypes at
  `GQ` 30 or more, join them;
- a heterozygote between the two that the estimates cannot orient is joined to them by reads in
  the same way.

Reads reach about a kilobase, so they change the links between close heterozygotes only. On a
simulated outbred cohort of 50 samples at 3× from 64 lineages, with a SNP every 100 bases, they cut
the wrong links inside phase sets from 5.2% to 3.6%, and the links they made themselves were wrong
15 times in 3,029. On the simulated cohorts above, links made by reads were wrong 1 time in 176, and the links reads broke
were mostly wrong ones in outbred cohorts: 5 of 7 at 50 samples from 64 lineages, where the wrong
links inside phase sets fell from 56 to 52 of about 1,300. Reads link records only next to ones
the estimates do orient, so in a cohort where the estimates orient nothing — a few outbred samples
switched on with `on` — the reads phase nothing either, while the same run with the refinement
off phases from reads alone.

The log reports a **read check**: of the links the estimates made that reads reach with confidence,
how many the reads contradict. It is a check of the population phase that does not come from the
refinement's own model, for links under a kilobase, and it changes nothing.

Tools that read `PS`, such as WhatsHap, or SHAPEIT with its option to use it, treat each phase set
on its own. A tool that ignores `PS` reads every `|` of a sample as one phase along the whole
chromosome, which these phase sets do not claim.

## Memory: windows, budget and thinned margins

The refinement reads the cohort in windows: a *core*, whose genotypes it writes, and a *margin* on
each side, read only to tell which samples each sample resembles there. A window holds every sample's
numbers at every record it spans, so its memory grows with the samples and with the records a
kilobase. Four options, under *Advanced* in `--help`, set how large windows may be:

- `--refine-window-bytes` is the memory one window may hold, 32 GB by default. Windows start at a
  2-Mb core with 1-Mb margins and shrink until they fit. Each thread also holds a few kilobytes a
  window record that the budget does not count, and the process as a whole has used up to about a
  fifth more than the budget.
- `--refine-window-floor` is the shortest a margin may become, 500 kb by default. A sample copies
  stretches of up to about 300 kb from the samples it resembles, so shorter margins cut off
  evidence. A window over budget even at the floor stops the run, and the error says how many
  bytes it needs.
- `--refine-core-floor` lets a core shrink below the margin floor while the margins stay at it,
  for a window that would otherwise stop the run. It defaults to the margin floor, so cores do not
  shrink below it. On 2,176 tomato samples a 125-kb core wrote as many wrong genotypes as a 500-kb
  one and took 60% longer, but on a simulated F2, whose haplotypes run the length of a chromosome,
  short cores wrote up to a fifth more wrong genotypes; it is a last resort.
- `--refine-margin-density`, 5 by default, walks at most that many records a kilobase in the
  margins farther than 50 kb from the core, evenly spaced; the records nearer the core, and the core
  itself, are walked in full, and stretches no denser are walked in full too. A record not walked
  holds no state in the walks and is held compressed until a window walks it. On 2,176 tomato
  samples, where the densest stretches hold about 90 records a kilobase, a window then held 41%
  less memory and took half the time, with as many genotypes wrong; at 1 a kilobase and below, the
  simulated cohorts began to lose genotypes. `0` walks every record.
- Before calling, and again at each window, pvc estimates the refinement's peak memory — the
  window, each thread's arrays and about 16% more for the memory allocator — and compares it with the
  machine's memory, or the job's limit where that is lower (a systemd unit's `MemoryMax`, a
  container's). Above 85% of it the run warns, says what to lower, and goes on; the warning names
  the `pvc resume-calling` command that finishes the run if the system kills it for memory after
  calling. On 2,176 tomato samples at 48 threads the estimate is about 1.3 times the budget.

## How good the fills are: the fill check

Every run reports, for each sample, how often its fills would be wrong, measured on its own data.
At the sites where the sample's own evidence settles its genotype, pvc asks what the refinement would
have written had the sample had no reads there, and compares that with the genotype its reads
gave. It is the usual test of an imputation — hide genotypes you know and see whether they come
back — done on a sample of such sites: at most 50 a sample in each window, and at most 1 in 20 of
the sample's informative sites.

The result is written beside the VCF as `<name>.fill_check.tsv`:

| column | meaning |
|---|---|
| `sample` | the sample |
| `filled` | the `./.` the refinement filled for this sample |
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
it — and consider `--refine-genotypes off`, or a calling run without the samples in question.
