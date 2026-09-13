# Phasing and imputation for ng: the algorithm families, the signal they use, and how the panel-based and panel-free methods connect

*Research report, 2026-09-12. This is the integrating document; six sub-reports in this directory
carry the detail and the citations. **Decision taken 2026-09-13**: phasing and imputation stay outside
the caller, an existing pipeline is tested first, and ng writes `PL`/`GP` for it —
[`07_decision_separate_tool_and_route.md`](07_decision_separate_tool_and_route.md) records the ruling,
its reasons, the algorithm, the tools and the test. Where §6 below recommends building an engine in
ng, that recommendation is superseded by 07. Every code claim below is checked against the cloned sources
under `tmp/phasing_research/repos/` (see §9 for how) and every paper claim against the paper the
sub-report cites. Where a number is this report's own arithmetic rather than a measurement, it says
so.*

---

## 0. What this report answers

The owner asked three things about adding a phasing and imputation module to ng, and a fourth
question is implied by the two population cases they described.

1. **What are the families of phasing and imputation algorithms, how do they work, what statistical
   signal does each exploit, and what are their advantages and disadvantages?** — §1 and §2.
2. **How does a reference panel of well-sequenced individuals enter the algorithms, how can a cohort
   of thousands of well-sequenced individuals be its own panel, how do the methods that work with
   every sample at low coverage do it, and is there a continuum between the two?** — §3.
3. **What do the two population cases need?** A segregant population with very few haplotypes where
   every individual may be at low coverage, and a natural population with many haplotypes where some
   individuals are deep and many shallow. — §4.
4. **What does this mean for ng**, which already holds per-read allele evidence (the chain ids) and
   per-sample genotype likelihoods? — §5, §6 and the questions in §7.

The short answers, each expanded below:

- **There is one model, not several families.** Every phasing and imputation method in production
  use today, from PHASE in 2001 to GLIMPSE2 and QUILT2 in 2023–25, is the same hidden Markov model:
  an individual's chromosome is a mosaic of *template* haplotypes, copied in segments, with
  occasional mismatches. The families differ in **what the templates are** (the two parents of a
  cross, a handful of fitted ancestral haplotypes, the other individuals of the cohort, or an external
  panel of thousands) and in **what evidence feeds the model** (hard genotypes, genotype likelihoods,
  or the reads themselves). The two families that are genuinely different — read-backed phasing and
  pedigree transmission — use signals that are not linkage disequilibrium at all.
- **A reference panel is a template set, not an algorithm.** SHAPEIT5 loads panel haplotypes into the
  same matrix as the cohort's own current haplotype estimates and runs one selection and one HMM over
  the union; Eagle2 sets its template library to the cohort's own haplotypes when no panel is given.
  A cohort of thousands of deep-sequenced individuals is its own panel by construction: each
  individual is phased against the others' current estimates, and the estimates are refined over
  10–20 sweeps. Below about 50 samples this stops working for a natural population, and SHAPEIT5
  refuses to run there.
- **Low coverage changes only the emission term.** The hard genotype is replaced by the three genotype
  likelihoods, collapsed to a two-allele likelihood given the current estimate of the other haplotype.
  Everything else stays. STITCH goes one step further and *learns* a small set of template
  haplotypes from all samples' reads at once, which is how it imputes with no panel at all.
- **Yes, there is a continuum, and it has two axes**: how many templates and where they come from,
  and how raw the evidence is. Down the first axis the only change is the template set and the
  recombination clock; along the second, only the emission. What differs between the ends is not
  the mathematics but how much support the data gives the model.
- **The two population cases are the two ends of the template axis.** A biparental cross is the
  model with two templates and a transition rate set by one or a few meioses; the reads of a
  low-coverage individual pool across a founder block hundreds of sites long, so accuracy above 99%
  from 0.1× is measured, not hoped for. A germplasm collection is the model with the cohort as
  templates; the deep individuals, once phased, are the panel for the shallow ones, and accuracy is
  a function of how many deep individuals there are and how diverse the population is.
- **For ng**, the module is a buffered second pass over the called cohort, not a per-record filter;
  it wants genotype likelihoods, which ng computes but does not write; and read-backed phasing from
  the chain ids is worth having inside the imputation model but not as a phasing product at 3×,
  where one adjacent heterozygous pair in five can be linked by any number of short reads.

---

## 1. The signal: shared haplotype segments

Two individuals of one population that carry the same allele at a site usually carry it because they
inherited it from one common ancestor, and with it they inherited the stretch of chromosome around
that site, back to the nearest recombination on each side. **Linkage disequilibrium** — the
non-independence of genotypes at neighbouring sites — is the population-level shadow of that fact.
So a chromosome, looked at over a window of a few hundred kilobases, is a mosaic of pieces each of
which is also present, nearly unchanged, in some other chromosome of the population. The piece
lengths are set by how many generations separate the two copies and by the local recombination
rate: a pair of chromosomes with a common ancestor *g* generations ago share pieces of about
100/(2*g*) centimorgans.

That is the whole signal the statistical methods use. A method that knows the population's
haplotypes can therefore say, for an individual with sparse or noisy data, "over this stretch you
are copying haplotype 17 on one chromosome and haplotype 302 on the other", and from that read off
both the phase (which alleles are on which chromosome) and the missing genotypes (what haplotypes 17
and 302 carry at the sites you did not observe). Phasing and imputation are the same inference.

**A segregant population is the same picture with different numbers.** In an F2 or a set of
recombinant inbred lines every chromosome is a mosaic of the two founders' haplotypes, the pieces are
tens of centimorgans long because only one to a few meioses have happened, and the founders' allele
at every site is either known (parents sequenced) or recoverable from the progeny. In a natural
population the templates are unknown, each is shared with few others over short segments, and the
posterior over "which template here" is decided by tens of sites per segment rather than hundreds.
Sub-report 01 §1 develops this; sub-report 04 §2 does the arithmetic for crosses.

**Two signals that are not linkage disequilibrium.** A read or a read pair that covers two
heterozygous sites says which two alleles are on one molecule: physical linkage, independent of the
population. And a pedigree says which alleles a child received from which parent: Mendelian
transmission. Both are used by their own families of methods (§2.5, §2.6) and both are what remains
at one sample, where the population signal is absent.

---

## 2. The families

The table gives the families as the literature names them; the text after it says how each works,
what it exploits, and where it breaks. Sub-report 01 §2 has the mathematics and the primary
citations; sub-reports 02–05 have the code.

| family | templates | evidence | signal exploited | representative tools |
|---|---|---|---|---|
| Parsimony and EM over haplotype frequencies (historical) | the haplotypes seen so far | hard genotypes, tens of sites | few haplotypes explain many genotypes | Clark 1990; Excoffier & Slatkin 1995 |
| Haplotype-copying HMM, cohort as its own panel | the other individuals' current haplotypes | hard genotypes (deep sequencing or arrays) | LD as shared segments | PHASE, MaCH, SHAPEIT1–5, Beagle 3–5, Eagle2 |
| Haplotype-copying HMM against an external panel | a phased panel of 10³–10⁵ haplotypes | hard genotypes (arrays) | the same, with more and better templates | IMPUTE2–5, minimac3/4, Beagle 5 imputation |
| Haplotype-copying HMM from genotype likelihoods or reads | an external panel (GLIMPSE, QUILT), or panel + the other targets (GLIMPSE1, Beagle 4.1 `gl=`) | genotype likelihoods, or reads with base qualities | the same, plus within-read phase | GLIMPSE1/2, QUILT1/2, Beagle 4.1 `gl=` |
| Haplotype-cluster models | K fitted "ancestral" haplotypes, learned from the sample | hard genotypes (fastPHASE) or reads (STITCH) | the same, with the templates as parameters | fastPHASE, Beagle 3, STITCH |
| Founder-mosaic HMMs for crosses | the 2 (or F) founders, known or inferred from the progeny | read counts or genotype likelihoods | linkage over tens of centimorgans from a known crossing scheme | LB-Impute, TIGER, GBScleanR, RABBIT/magicImpute, STITCH with K = 2 |
| Pedigree transmission | parents' haplotypes, peeled through the pedigree | hard genotypes | Mendelian inheritance | AlphaImpute2, AlphaPlantImpute2, FImpute |
| Read-backed phasing | none | reads spanning two heterozygous sites | physical linkage in one molecule | WhatsHap, HapCUT2; freebayes and GATK within a read length |

### 2.1 The haplotype-copying model, which is most of the table

**In words.** For one chromosome copy of one individual, walk along the sites and ask at each: which
of the K template haplotypes am I copying here? Two rules answer it. Along the chromosome the copied
template is kept with probability 1 − *t* and swapped for a random one with probability *t*, where
*t* grows with the genetic distance between the two sites (a recombination). At each site the copy
either matches the template's allele or, with a small probability, does not (a mutation or an
error). Run forward and backward over the sites and the result is, at every site, a probability
distribution over the K templates; from it come the phased alleles, the genotype posterior at
unobserved sites, and the dosage.

This is Li & Stephens' 2003 model, and the transition probability 1 − exp(−0.04·*N*ₑ·*d*/K), with *d*
in centimorgans and *N*ₑ the effective population size, is in the code of SHAPEIT5, Beagle 5.5 and
GLIMPSE2 verbatim (sub-report 01 §0, sub-report 02 §2 tabulates each tool's exact constants). A
diploid individual is handled either by a diploid HMM over pairs of templates (K² states) or, in every
modern tool, by the one-haplotype-at-a-time trick: fix or sample one chromosome copy, run the haploid
HMM for the other with the evidence explained by the first, alternate.

**What it costs, and why the tools differ.** The forward pass costs K numbers per site per haplotype.
With K equal to the whole cohort or panel that is unaffordable at biobank scale, so every tool's
engineering is in *choosing* the templates: SHAPEIT5 takes the neighbours of the target in a
positional Burrows–Wheeler transform (a sorted index of haplotypes by their suffix at each site);
Eagle2 the 10,000 least-mismatching haplotypes; Beagle folds candidates into a fixed 280 "composite"
haplotypes so that its per-sample cost is independent of cohort size; minimac4 and IMPUTE5 collapse
the panel to the unique haplotypes per block (sub-report 02 §9 compares them). For cohorts of hundreds
to a few thousand none of this is needed: GLIMPSE2 already runs with 2,000 conditioning haplotypes
per sample by default.

**Advantages.** One model serves phasing, imputation of missing genotypes, and imputation at sites the
individual was never typed at; the outputs carry calibrated uncertainty (genotype posteriors and
dosages); accuracy improves with every extra template.

**Disadvantages.** It needs templates. With the cohort as its own panel there is a floor in cohort
size (SHAPEIT5 refuses below 50 samples, sub-report 02 §3.1; Choi et al. 2018 measured 1.9% switch
error with an 85-sample panel against 0.30% with 22,690, sub-report 06 §2). Rare variants exist in
few templates, so both SHAPEIT5 and Beagle 5.2 phase them in a second stage restricted to carriers,
and singletons are placed by a coin flip or a coalescent rule (35% switch error in SHAPEIT5 at
147,754 UK Biobank genomes, sub-report 06 §2). Every tool in this family is biallelic-SNP-only at the
data-structure level (sub-report 03 §7.4, sub-report 06 §7). And, the point that matters for ng's
range: the hard-genotype tools assume a mismatch rate of 10⁻⁴ (SHAPEIT5) to 10⁻² (minimac4), while a
true heterozygote at three reads a site shows only one allele one time in four (sub-report 02 §10,
own arithmetic), so below roughly 10× the input has to be likelihoods.

### 2.2 The same model from genotype likelihoods or reads

**In words.** Replace "the observed allele matches the template" by "the reads are what one would see
if the template's allele were on this chromosome". GLIMPSE2 collapses the three genotype likelihoods
to a two-allele likelihood given the current allele on the other chromosome copy: if the other copy
carries the reference allele, the two candidates for this copy are scored by the likelihoods of the
0/0 and 0/1 genotypes (`GLIMPSE/phase/src/objects/genotype.cpp:94-115`, sub-report 03 §3.1). STITCH
and QUILT go one level lower: each read is scored against the template's allele at every site it
covers, so a read spanning two sites carries its within-read phase into the model
(`STITCH/STITCH/src/haploid.cpp:178-189`; QUILT's emission is copied from STITCH, sub-report 03 §3.2).
QUILT then samples, for every read, which of the two chromosome copies it came from — read-backed
phasing inside the imputation model.

**Advantages.** Calling and imputation become one model, which is the 1000 Genomes Project's design
and the reason 4× sequencing of 1,092 samples reached the accuracy of 15× single-sample calling for
alleles above 1% frequency (sub-report 06 §6.1). The per-read form uses information that site-level
likelihoods throw away: STITCH's own measurement is r² 0.97 with reads against 0.87 without, on 2,073
mice at 0.15× where SNPs are dense (sub-report 03 §3.3).

**Disadvantages.** GLIMPSE2 and QUILT need an external panel; Beagle 4.1's `gl=` mode, the only
production tool that ran the cohort as its own panel from likelihoods, was removed in Beagle 5
because it was about 1,200 times slower than GLIMPSE (sub-report 03 §3.4). GLIMPSE1 used the other
targets' current haplotypes as extra templates; GLIMPSE2 dropped that, and the code confirms it
(`GLIMPSE/phase/src/containers/haplotype_set.cpp:712-714` sizes the index by reference haplotypes
only; sub-report 03 §4.1). All three cap per-site depth (STITCH 50, QUILT 30, GLIMPSE2 40 for its own
pileup), so deep samples are down-sampled.

### 2.3 Haplotype-cluster models: the templates as parameters

**In words.** Assume the population descends from K ancestral haplotypes a few hundred generations
back, and fit those K haplotypes from all samples at once by expectation–maximisation: guess K
templates, assign every chromosome copy of every sample to a mosaic of them, re-estimate each
template's allele at each site from the reads assigned to it, repeat. fastPHASE did this from hard
genotypes in 2006; STITCH does it from reads and is the reference-free low-coverage imputer.

**What the signal is.** The same LD, but the method needs no pairwise match between individuals:
each template is estimated from the pooled reads of every sample that copies it over that stretch.
That is why STITCH works at 0.15× where no individual has enough reads to be anyone's template.

**Advantages.** No panel; runs at any depth; the per-interval recombination rate is learned. **With K
= 2 it is a biparental imputer that learns the two parents from the progeny**: the four states are
the founder-pair genotypes, the learned switch rate is the recombination fraction, and the state
posterior is the founder-origin map (sub-report 03 §5.8, the sub-report's reading of the code).

**Disadvantages.** It needs a cohort large enough to estimate the templates: about 500 low-coverage
samples for r² 0.95 in the published mouse and owl cohorts, and no gain below 100 (sub-report 06
§4.1). K and the number of generations are user choices. The EM has label-switching repairs that are
tuned for K in the tens. It is biallelic-SNP-only.

### 2.4 Founder-mosaic HMMs for crosses

**In words.** The same copying model with the templates fixed to the founders of the cross, the
switch rate derived from the crossing scheme (one meiosis for an F2, the map-expanded rate of a
recombinant inbred line, the pedigree of a multi-parent cross), and the emission from read counts with
an error rate — LB-Impute's switch probability is literally the Li & Stephens jump probability with
K = 2 (sub-report 04 §10). When the parents were not sequenced, the founders are recovered from the
progeny: by maximum parsimony of recombination (Xie et al. 2010), window clustering (FSFHap), joint
Viterbi over founder patterns (GBScleanR), or the K = 2 EM above.

**Why it is a different regime.** Reads pool across a founder block. At 0.5× with founders differing
at one site per 2 kb, a megabase holds about 250 informative reads; two concordant reads separate the
two homozygous founder states at 10,000:1, ten reads rule out heterozygosity at about 1,000:1, so the
resolution for a heterozygous-versus-homozygous block is about 40 kb at 0.5× and 200 kb at 0.1×
(sub-report 04 §3.5, own arithmetic checked against magicImpute's measured break points at 0.053× for
RILs and 0.11× for F2s). Published accuracy: above 99% from 0.1× with about 90 individuals
(LB-Impute), 90% of crossovers placed within 2 kb at 0.4× (TIGER), r² 0.98 on 2,177 medaka F2 at
0.5× (sub-report 06 §5).

**What binds instead.** Cohort size when founders are absent (at least 110 RILs in Xie et al. 2010; a
cliff below 200 F2 in Pierotti et al. 2024); sites where the founders are identical, which carry no
information at any depth; and marker-specific error — per-marker allelic read bias and mis-mapping —
which GBScleanR models and which is worth more than 25 accuracy points at 20× in its F2 simulation
(sub-report 04 §5.3). Recombinant inbred lines keep 1–6% residual heterozygosity, which a haploid
model silently mis-calls (sub-report 04 §7).

### 2.5 Read-backed phasing

**In words.** A read (or pair) covering two heterozygous sites votes for one of the two phasings.
Collect the votes, find the phasing that contradicts the fewest read bases (minimum error correction,
which is hard in general), and report the connected blocks. WhatsHap solves it exactly by dynamic
programming whose cost is exponential only in the read depth, after down-sampling to about 15–20×;
HapCUT2 uses a likelihood and a greedy max-cut (sub-report 05 §1). **The signal is physical, not
populational**, so it works on one sample and it is the only thing that phases a singleton.

**What short reads give.** The reach is one insert, and that is a geometric ceiling: at one
heterozygote per kb with 350-base inserts, 70% of adjacent heterozygous pairs can never be linked by
any number of short reads. Sub-report 05 §2 (own arithmetic, order-of-magnitude checked against
Delaneau et al. 2013's measurement) gives, for 150-base pairs: at 3×, one adjacent pair in five linked,
mostly by a single fragment, blocks of 1.2 heterozygous sites; at 30×, 29 pairs in 100 linked by 3.9
fragments on average. Published gains from adding reads to population phasing are all at 5× or more
on top of a large panel (22% longer switch distance in Delaneau et al. 2013). The one measurement of
read linkage inside a low-coverage LD caller, HapSeq on the 1000 Genomes pilot, cut genotype
discordance by 9–12% relative at about 4×, and by 30% in simulation (sub-report 06 §6.4).

### 2.6 Pedigree transmission

**In words.** With genotyped parents, the child's alleles are assigned to the parental haplotypes by
Mendelian rules, and haplotypes are "peeled" up and down a multi-generation pedigree; population HMMs
fill in what the pedigree leaves open. AlphaImpute2 and AlphaPlantImpute2 are the current tools
(sub-report 04 §6). Pedigree phasing is chromosome-length and near-perfect (0.009% switch error on
98 cattle trios, sub-report 06 §2), but every published accuracy figure for these tools is from array
genotypes, not low-coverage reads. Sub-report 04 §12 recommends treating a breeding population as a set
of biparental families, each with its parents as founders, rather than building a third machine.

---

## 3. The panel, the cohort as its own panel, and the continuum

### 3.1 How a panel enters

A reference panel is a set of phased haplotypes. It enters the copying model as the set of hidden
states: at every site, "which panel haplotype am I copying". Nothing else about the model changes. In
the code, SHAPEIT5 keeps panel haplotypes and the cohort's current haplotype estimates in one matrix
and runs one template selection and one HMM over the union
(`shapeit5/phase_common/src/containers/haplotype_set.cpp:43`, sub-report 01 §0); GLIMPSE2 stores the
panel as a pre-built, run-length-compressed positional Burrows–Wheeler index so that selecting 2,000
conditioning haplotypes per sample per iteration is cheap (sub-report 03 §4.1). What a panel changes
is *how many* good templates exist for a sample and *whether they are already phased*.

### 3.2 The cohort as its own panel

Every population phaser without a panel does the same thing: take the current haplotype estimates of
all *other* individuals as the templates, re-estimate this individual, move to the next, and sweep the
cohort 15–20 times, rebuilding the template index from the new estimates at each sweep. An
individual's own two haplotypes are never in its template set (sub-report 01 §2e cites the loops in
SHAPEIT5, Beagle 5.5 and Eagle2). Beagle 5.2 and later is not a sampler but a deterministic
forward–backward update that freezes the most confident heterozygotes first. A supplied panel is
simply extra rows in the same matrix that never change (sub-report 02 §3.2, §5.2).

So **a cohort of thousands of deep-sequenced individuals imputes and phases itself with no external
panel, using exactly the panel-based code.** The measured cost of cohort size, all with hard
genotypes and no panel (sub-report 06 §2): switch error 2.1% at 5,000 UK Biobank samples falling to
0.35% at 150,000 with Eagle2; 0.12–0.18% at 400,000 for SHAPEIT4, Beagle 5 and Eagle2; no published
curve below 500 samples. At 50–100 samples the nearest number is Choi et al. 2018's 1.9% switch error
for one European sample phased against an 85-sample panel at array density. SHAPEIT5 refuses below 50
samples in total (`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_scaning.cpp:46-48`);
Beagle 5.5 exits at one sample without a panel (sub-report 01 §2e).

### 3.3 Every sample at low coverage

Three answers exist, and they differ in what the templates are:

1. **Templates learned from the pooled reads (STITCH).** K ancestral haplotypes by EM, §2.3. About 500
   samples for r² 0.95, nothing below 100, in populations with few founders (mice, owls, wheat
   collections); 11,670 humans at 1.7× reached r² 0.92 with K = 40 (sub-report 06 §4.1).
2. **Templates = the other targets' current haplotypes, from likelihoods (Beagle 4.1 `gl=`;
   GLIMPSE1's mixed set).** The self-panel loop of §3.2 with the likelihood emission. It is what the
   1000 Genomes Project ran on 179 and then 1,092 samples at 2–6×, and Li et al. 2011's simulation
   found above 98% genotype concordance at 4× "with as few as 60 sequenced individuals" as their own
   panel (sub-report 06 §4.3). No production tool offers it today.
3. **Templates = a deep-sequenced subset of the same study, phased and frozen (the internal panel).**
   Sequence some individuals deep, phase them as in §3.2, hand them to GLIMPSE2 or QUILT2 as the panel
   for the shallow ones. This is the owner's second case exactly, and the literature's sizing
   evidence is all from populations with long shared haplotypes: 50–150 deep individuals reach r² about
   0.9 at 0.5–1× and 200–500 reach 0.95 (barn owl, strawberry, cattle series; sub-report 06 §4.3).
   For a diverse outbred population the only anchor is human: 1,000 population-specific genomes added
   to 1000 Genomes matched a 1,488-genome study panel (Sardinia).

### 3.4 The continuum, answered

Yes. Two axes, and every tool is one cell (sub-report 01 §4, sub-report 03 §10):

| templates ↓ / evidence → | hard genotypes | genotype likelihoods | reads |
|---|---|---|---|
| **2 founders** (cross) | TIGER-style breakpoint callers | LB-Impute, GBScleanR, magicImpute | STITCH with K = 2 |
| **K fitted clusters** | fastPHASE | — | STITCH |
| **the cohort itself** (re-estimated each sweep) | PHASE, MaCH, SHAPEIT1–5, Beagle 3–5, Eagle2 | Beagle 4.1 `gl=` (removed in 5.0); GLIMPSE1 (targets in the set) | none found |
| **external panel** (fixed, phased) | IMPUTE2–5, minimac3/4, Beagle 5 imputation | GLIMPSE1/2; Beagle 4.1 `gl=` with `ref=` | QUILT1/2 |

- **Down a column** the only change is the template set and its recombination clock: two founders and
  one to a few meioses; K learned rows and a learned rate; the cohort's own 2N − 2 haplotypes and a
  rate of 4*N*ₑ/(2N) per Morgan; a fixed panel with its haplotype count in the denominator.
- **Along a row** the only change is the emission: a mismatch probability; the two-allele likelihood
  given the other copy; the product over a read's bases.
- **What is genuinely different between the ends is the support the data gives the model.** With two
  known founders a handful of individuals at 1–3× is enough and the posterior is decided by hundreds of
  sites per block. With the cohort as templates the posterior is decided by tens of sites per segment,
  and accuracy is a function of N, coverage and panel size.
- **The code is not yet the continuum**, and the gap is instructive. Sub-report 03 §10 lists what would
  make GLIMPSE2 run without a panel: let the conditioning set point at target haplotypes and exclude
  a sample's own pair; rebuild the haplotype index every iteration because the targets move (this
  removes the pre-compression that is GLIMPSE2's main speed trick); seed iteration 0 without a panel;
  estimate allele frequencies from the cohort each iteration. That is GLIMPSE1. And what would make
  STITCH use a panel: freeze its learned rows at panel alleles and grow K to hundreds, at which point
  the K² diploid recursion is unaffordable and a read-label scheme is needed — that is QUILT, whose
  emission code is copied from STITCH. So the same machine sits at three points of the axis in three
  codebases, and none of them spans it.

---

## 4. The two population cases

### 4.1 Few haplotypes: a segregant population, every individual shallow

**What the literature says works.** A founder-mosaic HMM (§2.4) with the crossing scheme declared,
the per-site genotype likelihoods as its emission, and per-marker allelic bias and mis-mapping
estimated from the cohort. Above 99% genotype accuracy from 0.1× with about 90 F2 or backcross
individuals is measured (LB-Impute, simulation with 10,000 markers per 100 Mb); magicImpute holds 0.99
down to 0.05× for a two-founder RIL and 0.11× for an F2; TIGER places 90% of crossovers within 2 kb at
0.4×; 2,177 medaka F2 at 0.5× give r² 0.98 with STITCH at K = 16 (sub-report 06 §5).

**What it needs from ng.** Per-site genotype likelihoods, nothing more: the log-likelihood of a
founder-pair state over a block is the sum over its sites of the sample's likelihood for the genotype
that pair implies, so the chain ids add nothing here (sub-report 04 §0, §11). Plus ng's per-allele
read counts per sample, which are exactly the input of GBScleanR's two error estimators.

**What it must say at the edges.** With absent parents and fewer than about 100 progeny, founder
inference is under-determined (Xie et al. 2010's floor of 110 RILs; the cliff below 200 F2 in Pierotti
et al. 2024): report that, rather than a fitted answer. Keep diploid states and let the selfing
generations shrink the heterozygous prior, so residual heterozygosity in RILs is called and not
suppressed; reserve the haploid reduction for declared doubled haploids. Declare the crossing scheme,
then check it: an allele-frequency spectrum concentrated at ½ (¼ and ¾ for a backcross), per-sample
heterozygosity at the scheme's expectation, and a two-template fit that explains the reads as well as
an eight-template fit — sub-report 04 §8 gives six such statistics, of which the first is computable
from ng's cohort expected allele copies at 1× without any genotype call.

### 4.2 Many haplotypes: a collection with some deep and many shallow individuals

**What the literature says works.** Phase the deep individuals against each other (§3.2), freeze
them as the panel, impute the shallow ones from their genotype likelihoods against it (§3.3, route 3).
Depth buys accuracy fast to about 1× and then panel size takes over: with the 54,000-haplotype
Haplotype Reference Consortium panel GLIMPSE1 reaches r² above 0.9 at 0.3× for alleles above 5%
frequency and 0.8 at 1× for alleles at 0.1%; alleles at 0.01% need the 280,000-haplotype UK Biobank
panel (sub-report 06 §3.1). With an internal panel in a low-diversity population, 50–150 deep
individuals give r² about 0.9 at 0.5–1× (§3.3). Adding the shallow samples' own current haplotypes
to the template set (GLIMPSE1's design) helps when the deep set is small, and is the natural later
extension.

**What it must say at the edges.** At one sample there is no template: emit unphased, unimputed
genotypes and say so; never fit a silent template set from one sample (sub-report 01 §6). Rare alleles
carried only by shallow samples cannot be imputed from a panel that lacks them; report them as
uncertain rather than pull them to reference. The tools' identity-by-descent exclusion rules — SHAPEIT5
bans a template whose heterozygotes agree with the target's on more than 75% of a window, Beagle bans
pairs sharing 2 cM of identical genotypes — strip the template sets in highly inbred material (the
tomato panel's median inbreeding coefficient is 0.78) and in any segregant population, where sharing
both haplotypes is the signal, not contamination (sub-report 02 §10). Inbred cohorts with residual
heterozygosity are unmeasured in this literature: every inbred-line study either ran a haploid model
or discarded heterozygous calls before scoring (sub-report 06 §7).

---

## 5. What ng has, and what a module would need

Facts about ng's current shape, checked in the source on 2026-09-12:

- **The stream hands records over one at a time in genome order**, each record already genotyped,
  through a callback that filters and writes
  (`src/ng/run/psp_caller.rs:558-561`, `call_cohort_handing_each_record_over`). A copying HMM needs a
  window of many loci across all samples, so the module is a **buffered stage over a chromosome or a
  large window, or a second pass over the written cohort**, not a per-record filter. The streaming
  spec already rules that coupling across segments "must be a pass over the emitted records, never a
  coupling between in-flight segments" (`doc/devel/ng/spec/run_streaming.md` §4.3).
- **Genotype likelihoods exist in memory during calling but leave no trace in the record.** The
  per-sample evidence rows (`src/ng/calling/likelihood/mod.rs:1421-1446`) and the genotype table
  (`src/ng/calling/genotype_table.rs:134`) produce them; the locus result carries only the per-sample
  call and the cohort's expected allele copies (`src/ng/calling/mod.rs:3032-3050`); the VCF writes
  `GT:GQ:DP:AD` with `GT` always unphased, and `GP` only on request
  (`doc/devel/ng/spec/vcf_output.md` §7). Every low-coverage imputer wants likelihoods, and the
  hard-call tools are the worst choice at 1× (Beagle 5.4 F1 0.85 against 0.99 for QUILT2 in Vi et
  al. 2025, sub-report 06 §8).
- **Chain ids are consumed by the cohort merge and do not reach calling.** The merge uses them to
  follow a read across a locus's records and then hands calling per-(allele, read group) counts and
  quality sums (`src/ng/run/cohort_merge/build.rs:1473-1495`, `AlleleSupport`). Read-backed phasing
  between neighbouring heterozygous loci would need the ids retained, per sample and per allele,
  over an insert-length window, at the point in the merge where the ids and the elongated alleles
  coexist. The cost is per locus, not per position: about 0.5 MB for 1,000 samples at 30×, 150 MB in
  the 3,000-sample, 300×, dense corner (sub-report 05 §6.4, own arithmetic). The reach is one
  segment, because ids are re-minted per walked region and no locus crosses a segment boundary
  (sub-report 05 §6.3).
- **ng's records are multi-allelic and mix SNPs, indels and repeat tracts.** Every tool surveyed is
  biallelic at the data-structure level; GLIMPSE2 zeroes indel likelihoods by default and imputes
  indels only as passengers on a SNP scaffold (sub-report 03 §7.4). This is the clearest advantage an
  in-caller module has: the emission generalises to a per-template allele-probability vector, and
  QUILT2's multi-symbol index shows selection over multi-allelic templates is not a research problem
  (sub-report 03 §11).
- **ng fits a per-sample inbreeding coefficient** and writes it in the parameters file. The copying
  model can take it as a per-sample parameter (haploid behaviour at F near 1, as STITCH's
  `diploid-inbred` does globally), which no surveyed tool does for a mixed cohort (sub-report 06 §9).

---

## 6. What this report recommends

These are the integrating recommendations; each sub-report's own "Questions to discuss" section has
the finer ones, and §7 collects the decisions the owner should take.

1. **One HMM engine, two template sources.** Implement the haploid Li & Stephens recursion once — a
   template matrix, a per-interval transition vector, an emission that takes genotype likelihoods (and
   per-read likelihoods where available) — and vary only where the templates come from: the declared
   or inferred founders of a cross; the phased deep samples of a collection; an external panel when
   one exists. Sub-reports 01 §6, 03 §12 and 04 §12 reach this independently. The one large
   component only the many-haplotype case needs, template selection by a positional Burrows–Wheeler
   index, can wait until cohorts above a few thousand are a real target (GLIMPSE2 runs 2,000 templates
   per sample without it).
2. **Write the likelihoods first.** Emitting `PL` (and `GP`) is a FORMAT field and makes ng usable
   with GLIMPSE2, QUILT2 and STITCH today, through a VCF for the first and a site list for the other
   two; it is also the interface a second pass inside ng reads (sub-reports 02 §10, 03 §11, 06 §9).
3. **Build the segregant case first.** It needs no panel, its model is a few hundred lines, its
   literature shows above 99% from 0.1× with about 90 individuals, and its inputs (likelihoods,
   per-allele counts) already exist. The many-haplotype case is the deep-as-panel route on the same
   engine, plus the self-panel loop later (sub-reports 04 §12, 06 §9).
4. **Retain the chain ids behind a flag and measure before building on them.** The retention point
   and cost are known (§5). What is not known is what read linkage buys at 3× on a plant cohort:
   nobody has measured it. One afternoon on the tomato benchmark and the GIAB trio gives the fraction
   of adjacent heterozygous pairs with a witness, the witness count, and the disagreement rate at 30×
   (sub-report 05 §7). HapSeq's 9–30% relative reduction in discordance at 4× is the number to try to
   reproduce (sub-report 06 §6.4).
5. **Read-backed phasing as a VCF product only for deep samples.** At 3× a `PS` field would be true
   and nearly empty (one pair in five, blocks of 1.2 sites); at 30× a link has about four witnesses
   and its error rate can be read off them. The only mainstream consumer of `PS` is SHAPEIT4, and
   SHAPEIT5 dropped the feature (sub-report 05 §3.3, sub-report 02 §11).
6. **Say what happens at N = 1 and at 20 progeny.** Emit the genotypes unimputed and unphased with a
   diagnostic; never a fitted zero. This is the project's own rule and every surveyed tool either
   refuses or degrades to allele-frequency guessing at that end.
7. **Report the module's own accuracy the way the literature does.** Aggregate r² by allele-frequency
   bin against held-out deep samples, and switch error against a truth set; the truth-free scores
   (INFO, DR2) are the same formula in GLIMPSE2, STITCH and QUILT and should be emitted for
   comparability, but they miss regional failure and are inflated at rare variants (sub-report 06 §1,
   §8).

---

## 7. Questions to discuss

Each carries a recommendation and its trade-off; the sub-report that argues it is named.

1. **Is the module one engine with pluggable templates, or two modules (crosses; collections)?**
   *Recommendation: one engine, two template providers* (§6.1). Trade-off: the founder-mosaic case
   wants a 4-state diploid HMM with no iteration, the collection case wants the one-haplotype-at-a-time
   trick with 10–20 sweeps; a shared engine carries both shapes. (01 §6, 03 §12, 04 §12)
2. **Diploid HMM or one haplotype at a time?** *Recommendation: the trick everywhere, with K = 2
   diploid as a special case.* It is linear in K, takes likelihoods with a three-line transformation,
   and is what GLIMPSE2, QUILT, Beagle 5 and IMPUTE5 do. Trade-off: iteration, and a small diploid
   pass to fix the phase between the two copies. (01 §6)
3. **Which evidence: hard calls, likelihoods, or reads?** *Recommendation: likelihoods as the
   baseline; reads as an option for the many-haplotype case, not the cross.* STITCH measured 0.97
   against 0.87 with and without reads where SNPs are dense, and no difference on humans at 1.7×;
   for a cross the founder pair fixes the phase and reads add nothing. (03 §12, 04 §12)
4. **Where does the module sit?** *Recommendation: a second pass per chromosome over the written
   cohort, reading a columnar store of likelihoods rather than re-reading alignments.* Trade-off: a
   sidecar file per cohort, or `PL` in the VCF and a VCF re-read. (01 §6, 04 §11)
5. **Deep samples as a frozen panel, or a full self-panel loop over everyone?** *Recommendation:
   frozen panel first.* It runs on today's engine and matches GLIMPSE2's design; the self-panel loop
   over shallow samples (GLIMPSE1, Beagle 4.1) is the later extension for collections whose deep set
   is small. Trade-off: rare alleles only the shallow samples carry stay unimputed. (01 §6, 03 §12)
6. **Declare or detect the crossing scheme?** *Recommendation: declare, then check with the
   statistics of sub-report 04 §8, and refuse to run when they contradict the declaration.* Inference
   cannot recover the generation count and fails silently on segregation distortion. (04 §12)
7. **Founders absent and few progeny.** *Recommendation: run founder inference only above about 100
   progeny; below that, accept a user-supplied founder VCF or emit unimputed genotypes with a
   diagnostic.* Trade-off: a user with 40 RILs and no parents gets nothing, and is also the user for
   whom every published method is unvalidated. (04 §12)
8. **Per-marker allele bias and mis-mapping: always on?** *Recommendation: on by default, two
   estimate–decode cycles, gated on cohort size.* It is the largest accuracy lever in the segregant
   literature and ng has the counts. (04 §12)
9. **Multi-allelic templates from day one?** *Recommendation: yes — represent a template site as an
   allele-probability vector, not a bit, with a sparse encoding for rare carriers.* The biallelic
   restriction of every surveyed tool is an early data-structure choice none of them undid.
   Trade-off: memory per template site. (03 §12)
10. **What to do with the inbreeding coefficient.** *Recommendation: a per-sample parameter of the
    HMM.* Trade-off: unmeasured; it would have to be validated on the tomato cohort against its
    20–40× accessions. (06 §9)
11. **Should ng phase within a read length the way freebayes does, by emitting one allele for the
    whole stretch?** *Recommendation: decide this before any `PS` question, because it changes the
    record, not a tag.* ng's multi-base multi-allelic locus is closer to the freebayes form; GATK's
    tags keep single-site records. Trade-off: haplotype alleles inflate allele counts at low depth
    exactly where the genotyper is weakest. (05 §7)
12. **Which measurements to run before designing further.** *Recommendation: two.* (a) The read-link
    census of §6.4 on tomato and GIAB. (b) Down-sample the deep samples of each benchmark cohort to
    0.5×, 1×, 2× and 4×, impute with STITCH and with GLIMPSE2 using the remaining deep samples as the
    panel, and score r² by frequency bin — the protocol every paper uses and the yardstick for any
    module ng builds. (03 §12, 05 §7)

---

## 8. Where the detail is

| file | what it holds |
|---|---|
| [`01_statistical_models.md`](01_statistical_models.md) | the families' mathematics; the continuum table; output quantities; glossary; primary papers |
| [`02_panel_based_tools.md`](02_panel_based_tools.md) | SHAPEIT5/4, Eagle2, Beagle 5.5, minimac4, PBWT, IMPUTE5 at code level: inputs, self-panel loops, state selection, HMM constants, rare variants, chunking, `PS` handling; what ng could feed them today |
| [`03_low_coverage_imputation.md`](03_low_coverage_imputation.md) | GLIMPSE2, QUILT/QUILT2, STITCH, Beagle 4.1 `gl=`: how likelihoods and reads enter, panel representation and selection, STITCH's EM and its K = 2 behaviour, cost, accuracy, the continuum in code |
| [`04_few_haplotype_populations.md`](04_few_haplotype_populations.md) | crosses, multi-parent and pedigree populations: the founder-mosaic HMM, founders from progeny, error modelling, inbred lines, detection statistics, accuracy, what it implies for ng |
| [`05_read_backed_phasing.md`](05_read_backed_phasing.md) | WhatsHap, HapCUT2, SHAPEIT2 reads, SHAPEIT4 `--use-PS`, freebayes and GATK; the short-read arithmetic; VCF conventions; turning ng's chain ids into phase evidence |
| [`06_accuracy_evidence_and_calling_interplay.md`](06_accuracy_evidence_and_calling_interplay.md) | measured accuracy versus depth, N, panel size and frequency; internal panels; segregant populations; imputation coupled to calling (1000 Genomes, HapSeq, polyRAD); failure modes; design implications |
| [`07_decision_separate_tool_and_route.md`](07_decision_separate_tool_and_route.md) | the owner's decision (2026-09-13): no imputation inside ng, and why; the recommended algorithm with its statistics; the existing tools to use instead; the test that decides whether they serve our case |

The cloned sources are under `tmp/phasing_research/repos/` (a shallow clone of each, 2026-09-12; the
official SHAPEIT5 repository is disabled on GitHub, so the copy is a fork, `JosephLalli/shapeit5_old`;
Beagle 5.5's Java sources are from the author's own zip). A shared brief for the sub-reports is at
`tmp/phasing_research/brief.md`, and the ng facts of §5 at `tmp/phasing_research/ng_facts.md`.

## 9. How the citations were checked, and what could not be

Every `path:line` citation in the six sub-reports (about 590 distinct ones) was resolved against the
cloned sources and the ng tree with `tmp/phasing_research/check_citations.sh`: 546 resolve
mechanically, and the remainder were resolved by hand because the sub-report abbreviates the path
(`HMM.R` for `GBScleanR/R/Methods-GbsrGenotypeData_HMM.R`, `functions.R` for
`STITCH/STITCH/R/functions.R`, `gatk/.../AssemblyBasedCallerUtils.java`) and gives the full path at
its first use. All of those point to an existing file long enough to hold the cited line. Two
exceptions: the AlphaPlantImpute2 citations point into its `tinyhouse` submodule, fetched
separately; and the two TASSEL `FILLINImputationPlugin.java` citations in sub-report 04 were read
on GitHub, since TASSEL is not cloned here. This check confirms that a cited line exists, not what
it says; the sub-report authors report having opened each one.

Paper claims were opened by the sub-report authors at the source; each sub-report ends with a list of
what it could not open. The recurring gaps: values that exist only as plotted points (SHAPEIT4 Fig. 2a
below 400,000 samples, Beagle 5.2's size series, STITCH's Fig. 4a cells) are described in the papers'
words and not quoted as numbers; the methods sections of Browning & Browning 2016 (Beagle 4.1 `gl=`)
and Swarts et al. 2014 (FSFHap/FILLIN) were paywalled or unreachable, so those two models are described
from manuals and from the source code; Patterson et al. 2015 (the weighted minimum-error-correction
algorithm) is cited through Garg et al. 2016; and GLIMPSE1's accuracy figures are from its preprint.
