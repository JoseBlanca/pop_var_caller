# False SNPs from low-quality bases beside G and C runs

**Status: BUILT, 2026-11-01.** Option D (§4) chosen by the owner, measured against GATK's
approach (§8), and on by default at weight 1 (`--low-base-quality-weight`, owner 2026-11-01).

## 1. The problem

**On HG002 chromosome 20 at 30×, pvc's largest group of false SNPs that GATK does not make is
built from bases the sequencer itself marked as poor.** Of the 317 false SNPs more than a repeat
tract away, 143 are pvc's alone, and in 72 of those the bases showing the variant average below
Phred 20 while the reference bases at the same sites average 28. Most of the 72 sit in G- or C-rich
sequence, and in most the variant base would lengthen a run of the same base: `GGGG[T]GGG` read as
`GGGG[G]GGG`.

**What a site looks like.** At chr20:9,839,088 nine reads show G at Phred 5–13 and five show the
reference T. The genotype likelihood gives `0/1` a lead of 11 natural-log units over `0/0`, and the
site is written at QUAL 9.9. At chr20:595,931 thirteen of 21 reads show G at an average Phred 15;
the variant reads outnumber the reference reads, so the allele-balance penalty charges nothing and
the site is written at QUAL 119.

The data: HG002, HiSeq 2500, 148-base reads aligned by novoalign, subsampled to 30×; GIAB v4.2.1
truth; pvc at commit `4d00603c` under `--defaults`.

## 2. What the reads show

Four measurements, each made on the same 30× alignment.

**The variant bases are late in the read.** At the 72 sites the variant bases sit, on average, at
sequencing cycle 70 to 120 of 148; the reference bases at the same sites average about cycle 55.
Quality falls along an Illumina read, so this is the same fact as the low quality, seen from the
instrument's side.

**Reported qualities are honest on average, and not beside a run.** Over 112 million bases on
GIAB-confident reference ground, at least 10 bases from any truth variant (seven 0.5 Mb stretches of
chromosome 20), a base the sequencer called Phred 8–12 was wrong at an observed Phred of 12.0, and
one called 13–17 at 17.5. But where the base sits beside a run of four or more of another base,
Phred 13–17 bases were wrong at Phred 13.5 and Phred 18–22 bases at 18.4: **about 2.5 times the
stated rate**.

**A wrong base beside a run is usually the run's base.** pvc's likelihood divides a wrong read's
probability by three, one share for each base it could have turned into
([`read_likelihoods.md`](read_likelihoods.md) §3.5). Beside a run of four or more, 72 to 83 of every
100 wrong bases at Phred 13–22 were the run's base, not 33.

**The errors cluster: many reads at one site are wrong the same way.** This is what decides the
design. Correcting each read's error for the two biases above does not remove the calls. Scored
offline with pvc's formula on the same reads, with a heterozygosity prior of one in a thousand:

| how each read is charged | of the 72 still favouring a variant genotype |
|---|---|
| pvc as built | 72 |
| bases below Phred 13 dropped | 62 |
| the wrong-base share set from the run context (§2's 72–83 in 100) | 56 |
| the error rate and the wrong-base share both set from the run context | 72 |

Thirteen wrong reads in 21 at one site is not what a 1-in-20 per-read error rate produces, however
it is spread. Each read is being treated as independent evidence and they are not. **So the remedy
has to look at the site's reads together**, which is what the artifact correction
([`calling_quality.md`](calling_quality.md) §6) already does for strand and allele balance.

## 3. How GATK and freebayes avoid these calls

**GATK HaplotypeCaller** has three mechanisms that act here, in likely order of effect:

- in the pair-HMM, every base below Phred 18 is set to Phred 6, so all the variant bases here count
  as one-in-four errors (`PairHMMLikelihoodCalculationEngine.java:312`, `PairHMM.java:31`);
- a site needs QUAL 30 to be called (`GenotypeCalculationArgumentCollection.java:28`); 40 of the 72
  have QUAL below 30 in pvc;
- the assembler cuts each read at any base below Phred 10 before building haplotypes, and prunes a
  path fewer than two reads support, so a variant seen only at Phred 5–9 never becomes a candidate
  (`AbstractReadThreadingGraph.java:1039`).

Its base-quality and read-position rank-sum tests are annotations only; they never filter.

**freebayes** has no default base-quality threshold (`--min-base-quality 0`). It charges a read the
genotype cannot explain its full error, not a third of it, and discounts the error side by its
read-dependence factor of 0.9 (`DataLikelihood.cpp:77`, `:159`). Why it misses these sites was not
traced; its haplotype window widening at runs of five or more is the likeliest reason.

Neither caller has a mechanism pvc could port that handles the clustering in §2: GATK's is a quality
cut plus a high calling threshold, both of which act per read or per site without regard to
whether the variant reads differ from the reference reads.

## 4. The options

Each option is scored at 30× against the 72 false SNPs, the other 245 false SNPs, and every one of
the 69,633 SNPs pvc calls correctly on chromosome 20 (43,989 heterozygous, 25,644 homozygous). A
site is counted as removed when its QUAL after the option falls below 1, pvc's write threshold.
These are offline scores on the reads, not pvc runs.

### Option A — a base-quality floor in the likelihood

Drop every base below Phred 13 from the evidence, as bcftools did by default for years and as GATK
does for assembly at Phred 10.

- **30×:** 10 of the 72 removed (table in §2).
- **3×:** a real heterozygote with one good and one poor variant read loses half its evidence.
- **One sample or thousands:** the same per read.

**Not recommended.** It barely touches the problem and costs evidence everywhere.

### Option B — correct each read's error for the run context

Set a read's error rate and its wrong-base share from whether its base would lengthen a run: §2's
measured 2.5 times and 72–83 in 100.

- **30×:** at most 16 of the 72 removed (table in §2).
- **One sample or thousands, 3× or 30×:** the same per read. But the correction table was measured on
  one HiSeq 2500 run; another instrument needs its own, and pvc's error-rate fit has no quality or
  context axis to learn it from ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2).

**Not recommended for this problem.** The divide-by-three is genuinely wrong beside runs and is worth
its own change later, but it does not remove these calls.

### Option C — flag a site whose variant reads average a low quality

Take the variant reads' mean Phred and the reference reads' mean Phred at the site, each pulled
toward the sample's average (Phred 35 on this alignment) by a prior worth `k` reads, and flag the
site when the variant side falls below a threshold `T`. pvc already carries what this needs: every
stored observation keeps the sum of its reads' log error, so the mean is one division, and the
sample's average is what pvc's error-rate calibration already computes.

| prior `k`, threshold `T` | of the 72 | other false SNPs | true SNPs lost |
|---|---|---|---|
| 0, 20 (a plain average) | 67 | 2 | 22 |
| 3, 22 | 54 | 2 | 12 |
| 5, 22 | 34 | 0 | 1 |

- **3×:** the prior protects a site with two or three poor variant reads: two reads at Phred 10
  with `k = 5` average to 28.6 and are not flagged.
- **Thousands of samples:** the reads pool across the cohort, the prior stops mattering, and the
  flag judges the cohort's variant reads' mean. A real variant's reads have the reference reads'
  quality, so it is not flagged however many reads pool. An artifact that recurs in every sample is
  flagged more surely than in one.
- **Cost:** small. No change to the stored pileup files.
- **Weakness:** it is a cut, not a probability. It does not compare the two sides as a test does,
  and its threshold is a number read off one sample.

### Option D — a third artifact penalty: are the variant reads' qualities a fair sample of the site's?

The strand penalty asks whether the variant reads' forward-strand share matches the reference
reads'. This asks the same of the share of reads below Phred 20, with the same machinery
([`calling_quality.md`](calling_quality.md) §6.5): the variant reads' low-quality count is tested
against a beta-binomial whose parameters are the reference reads' low-quality and good counts plus a
prior worth `K` reads centred on 2 in 100 (2.2 in 100 on this alignment; a constant, because calling
does not have each sample's own share). Only an
**excess** of low-quality variant reads is charged. The tail, as a Phred, is multiplied by a weight
and by the strand test's ramp (nothing at three or fewer variant reads, full charge from seven) and
subtracted from QUAL.

| prior `K`, weight | of the 72 | other false SNPs | true SNPs lost |
|---|---|---|---|
| 10, 1 | 29 | 2 | 1 |
| 30, 1 | 43 | 4 | 4 |
| 30, 2 | 54 | 7 | 17 |
| 30, 4 | 65 | 16 | 70 |

- **3× at one sample:** silent below four variant reads, which is the intended behaviour of a test
  with no power.
- **A cohort at 3×:** the counts pool across samples, so a variant present in many samples reaches
  full power. That is the strand test's open question Q1 ([`calling_quality.md`](calling_quality.md)
  §13), unchanged.
- **Thousands of samples:** the weight makes the sign of `QUAL − weight × tail` depend on how lopsided
  the reads are per read, not on how many there are — the argument §6.5 makes for the strand test.
- **The prior has a cost where the reference reads are few and poor.** Twenty reference reads half
  below Phred 20 set an expected share of about 21 in 100 against the prior's pull, so a real variant
  whose reads are as poor as the reference reads is charged 10 to 20 Phred at weight 1. At 300 reads
  a side the pull is from one half to 46 in 100 and the charge a few Phred, against a site quality in
  the thousands.
- **Where the reference reads are poor too** (chr20:595,931, reference reads averaging Phred 11), the
  expectation follows them and charges little. That is the comparison working as designed, and the
  reason some of the 72 survive.
- **Cost:** each stored observation gains one count, reads below Phred 20, beside its forward-strand
  count. That is a change to the stored pileup format (version 1.0 to 1.1) and to the 28 source files
  that carry the strand count from the read walk to the quality step.

## 5. Recommendation

**Option D, with a prior of 30 reads and the weight set by the measurements in §6, between 1 and 2.**
At 30× that removes 47 to 61 false SNPs for 4 to 17 true ones. It is the only option that compares
the variant reads with the reference reads at the same site, which is what keeps it from charging a
real variant in a poor region; it shares the strand test's arithmetic, its ramp and its argument at
cohort scale, so it adds one question to a mechanism already measured rather than a new kind of
mechanism; and it lowers QUAL rather than cutting, so a strong site survives a mild excess.

**If the stored-format change is unwelcome, Option C is the fallback**: about the same trade at 30×,
built in a few files, at the price of a fixed threshold instead of a test.

Either is built behind a switch, off by default until measured. *(Built as Option D; on by default
at weight 1 since the measurements in §8.)*

## 6. How it will be measured

- HG002 chromosome 20 and the GIAB tandem-repeat benchmark at 5× and 30×, scored by vcfeval:
  false calls removed against true calls lost.
- Tomato, 63 accessions at about 3×, against GATK's joint calls: genotypes removed that GATK does
  and does not call. The BAQ-free extract is needed here, since BAQ lowers qualities beside runs and
  pvc does not support BAQ'd input.
- If the 3× cohort result is ambiguous, a simulation that plants clusters of low-quality run-base
  errors late in reads, built on [`examples/ng_tract_simulator.rs`](../../../../examples/ng_tract_simulator.rs).

## 7. Related, not in scope

- **The read-position half of the strand penalty never fires.** A read counts as "placed left" when
  its alignment starts before the site, which is nearly always true; freebayes compares the bases
  left and right of the site instead. Fixing it would not catch these sites as defined, since a
  late-cycle base is on the right of a forward read and the left of a reverse one, but the two
  changes touch the same counts and could share one format change.
- **The divide-by-three beside runs** (Option B) is a real bias, measured in §2, and a candidate for
  its own change to the likelihood.

## 8. Measured: the penalty against GATK's approach (2026-11-01)

**Built**: each stored observation counts its reads whose base quality inside the locus is Phred 1
to 19 (Phred 0 is the losing mate of an overlapping pair and does not count); the count travels
through the merge to the artifact correction, and `--low-base-quality-weight` sets the weight, zero
leaving the test off and the output unchanged. The stored pileup format gains one trailing field;
files written before it read with every count zero.

**GATK's approach, for comparison**: a prototype switch, since deleted, treated every base Phred 1 to
17 as Phred 6 when a read's error was minted, as GATK's pair-HMM does; GATK's QUAL 30 calling
threshold was applied when scoring. The two together are what GATK does that pvc does not.

**HG002**, pvc under `--defaults`, scored by vcfeval (genotype match) against GIAB v4.2.1. Each cell
is false calls removed / true calls lost, against the run with both switches off:

| | chr20 30× | chr20 5× | tandem repeats 30× | tandem repeats 5× |
|---|---|---|---|---|
| baseline (false / missed) | 552 / 757 | 4,616 / 24,614 | 2,412 / 14,268 | 2,540 / 24,059 |
| penalty, weight 0.5 | 27 / 1 | 9 / 0 | 76 / 11 | 7 / 1 |
| **penalty, weight 1** | **53 / 6** | **21 / 2** | **130 / 19** | **9 / 2** |
| penalty, weight 2 | 69 / 29 | 26 / 8 | 194 / 84 | 16 / 4 |
| GATK's base-quality squash | 28 / 4 | 6 / 216 | 119 / 6 | 63 / 66 |
| QUAL 30 threshold | 161 / 303 | 1,180 / 6,518 | 488 / 533 | 777 / 1,975 |
| squash and QUAL 30 | 175 / 334 | 1,091 / 6,702 | 539 / 582 | 764 / 2,038 |

Of the 72 false SNPs this document is about, at chr20 30×: weight 0.5 removes 24, weight 1 removes
43, weight 2 removes 54; the squash removes 22, the QUAL 30 threshold 40, and the two together 52.

**Tomato**, 62 accessions at about 3× over 80 regions of 100 kb, BAQ-free CRAMs, against GATK's joint
calls (made from the BAQ'd CRAMs) on the 56 accessions both VCFs name. The unit is one sample's
genotype at one SNP where either caller gives the sample the allele. Change against the baseline:

| | genotypes only pvc gives | genotypes both give alike | genotypes only GATK gives |
|---|---|---|---|
| baseline | 145,267 | 931,078 | 156,495 |
| penalty, weight 0.5 | −3,327 | −284 | +320 |
| **penalty, weight 1** | **−6,113** | **−541** | **+623** |
| penalty, weight 2 | −10,884 | −1,066 | +1,263 |
| GATK's base-quality squash | −22,495 | −678 | +2,550 |

The squash also moved 2,134 indel genotypes only pvc gives and 1,875 SNP genotypes the two callers
give in different numbers of copies; the penalty moved under 300 of either at any weight.

**What it shows.** The penalty does what it was built for at 30× and stays out of the way at 5×,
where the ramp leaves it nearly silent. GATK's squash removes about half as many false calls at 30×
and costs 216 true calls on chr20 at 5×, because a low-depth heterozygote with one or two poor
variant reads loses most of its evidence. GATK's QUAL 30 threshold is what removes most of the 72
in GATK's hands, and on pvc's quality scale it costs 303 true calls at 30× and 6,518 at 5×.

**Weight 1 is the default (owner, 2026-11-01)**: at 30× it removes 53 false calls on chr20 and 130 on the
tandem-repeat benchmark for 6 and 19 true ones; weight 2 adds 16 and 64 more removals for 23 and 65
more losses.
