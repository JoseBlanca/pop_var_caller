# The statistical models behind phasing and imputation, and the panel-based / panel-free continuum

Sub-report 01 of the phasing/imputation research (2026-09-12). Reader: a geneticist who knows the
biology and population genetics but is not a statistician. Every claim about a method carries a
paper citation; every claim about code carries a `path:line` into
`tmp/phasing_research/repos/`. Where a paper's full text could not be opened, the report says so.

Companion sub-reports: 04 covers pedigree methods, 05 covers read-backed phasing; both are only
pointed to here.

---

## 0. What this report concludes

1. **Every phasing and imputation method in use today is one model.** An individual's haplotype
   is treated as a mosaic of *template* haplotypes, copied in segments with occasional errors. The
   families differ only in what the templates are: the two parents of a cross, a handful of fitted
   "ancestral" haplotypes, the other individuals of the cohort, or an external panel of thousands.
   The mathematics (a hidden Markov model whose hidden state is "which template am I copying
   here?") is the same at every point of that range; the Li & Stephens 2003 transition and
   emission formulas appear verbatim in SHAPEIT5, Beagle 5, GLIMPSE2, STITCH and QUILT
   (§2c, §4).
2. **A reference panel is not a different algorithm; it is a different template set.** SHAPEIT5
   loads reference haplotypes into the same matrix as the cohort's own haplotypes and runs the same
   selection and HMM over the union (`shapeit5/phase_common/src/containers/haplotype_set.cpp:43`).
   Eagle2 sets its template library to the cohort's own 2N haplotypes when no panel is given
   (`Eagle/src/Eagle.cpp:1427`). What a panel changes is *how many* good templates exist for a
   sample and *whether they are already phased*.
3. **The self-panel methods have a floor in cohort size**, because their only templates are the
   other samples. SHAPEIT5 refuses to run below 50 samples (cohort plus panel)
   (`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_scaning.cpp:46-48`); Eagle1's
   authors recommend SHAPEIT2 below about 15,000 samples (Loh, Palamara & Price 2016, Discussion);
   Eagle2's cohort-only switch error falls from 2.10% at 5,000 samples to 0.35% at 150,000
   (Loh, Danecek et al. 2016, Fig. 5b). At one sample there is no template at all: the only
   signals left are reads (sub-report 05) and pedigree (sub-report 04).
4. **Genotype likelihoods enter as the emission term, nothing else changes.** That is what makes
   low-coverage imputation and genotype calling one model: GLIMPSE2 collapses the three-genotype
   likelihood to a two-allele likelihood given the current estimate of the other haplotype
   (`GLIMPSE/phase/src/objects/genotype.cpp:94-115`); STITCH and QUILT go one step further and
   use each read's base qualities directly (`STITCH/STITCH/src/haploid.cpp:179-186`).
5. **The two population cases the owner cares about are the two ends of one axis.** A biparental
   cross is the model with two templates and a transition rate set by one or a few meioses; a
   germplasm collection is the model with the cohort as templates and a transition rate set by
   thousands of generations of coalescence. One engine with a pluggable template set and a
   pluggable recombination clock serves both (§4, §6).

---

## 1. The signal: shared haplotype segments

### 1.1 Why genotypes at neighbouring sites are not independent

Two chromosomes sampled from a population share a common ancestor at every position. Between that
ancestor and today, recombination has broken each chromosome into pieces inherited from different
ancestral chromosomes. Over a short enough stretch, no recombination has happened since the common
ancestor, so the two chromosomes carry the *same* sequence there, apart from mutations that
occurred since. Such a stretch is **identical by descent** (IBD). Linkage disequilibrium — the
correlation of alleles at nearby sites — is the population-level shadow of these shared segments:
alleles are correlated because the chromosomes carrying them are, in pieces, copies of the same
ancestral chromosomes (Browning & Browning 2012, *Annu Rev Genet* 46:617–633, reviews IBD
detection and its use for phasing and imputation).

The coalescent view says how long the shared pieces are. Two chromosomes whose common ancestor
lived *t* generations ago have had 2*t* meioses to recombine between them, so the expected length
of a shared segment is about 1/(2*t*) Morgans: 100 generations gives about 0.5 cM, 10,000
generations gives about 0.005 cM. Close relatives share long segments; the population at large
shares short ones. In a large outbred population the distance to the nearest common ancestor of two
random chromosomes is of the order of the effective population size *N*ₑ generations, which is why
the useful window is short — typically well under a centimorgan in humans — but never zero.

### 1.2 Haplotypes as mosaics of a few templates

Put the two facts together: over a window shorter than the typical IBD segment, any haplotype in
the sample is a copy of *some* other haplotype in the sample (or of an ancestral haplotype that
several sample members also copy). Over a longer window it is a mosaic — a copy of one template
here, a different template past the recombination break, and so on. This is the picture Li &
Stephens drew (2003, *Genetics* 165:2213–2233, Fig. 2): "*h₄ as having recent shared ancestry with
the haplotype that it copied in each segment ... the copying process is Markov along the chromosome,
with jumps ... occurring at rate ρ/k per physical distance*". The three numbers that govern it are:

- **how many templates** there are (k) — more templates means a closer match exists for any
  segment, but also more state to carry;
- **how often the copied template changes** — set by recombination rate × time since the common
  ancestor, i.e. ρ = 4*N*ₑ*c* per base in a constant-size population (Li & Stephens 2003, eq. A1);
- **how faithfully it is copied** — mutation since the common ancestor plus genotyping error,
  the "imperfect" in "imperfect mosaic".

### 1.3 A segregant population is the same picture with different numbers

In an F2, RIL, backcross, doubled-haploid or MAGIC population every chromosome is a mosaic of the
founder haplotypes, and the number of breakpoints is set by the number of meioses in the pedigree,
not by population history:

- **Templates**: the founders — two for a biparental cross, a handful for MAGIC.
- **Transition rate**: one meiosis' worth of recombination for an F2 or backcross, a few for
  later generations. For recombinant inbred lines, repeated selfing or sib-mating adds map
  expansion: tightly linked loci recombine about twice as often as in one meiosis under selfing
  and about four times under sib mating (Haldane & Waddington 1931, *Genetics* 16:357–374;
  summarised with those factors in Crow 2007, *Genetics* 176:729–732). LB-Impute, the HMM built
  for this case, uses "a hidden Markov model that incorporates marker read coverage to determine
  variable emission probabilities" and a modified Viterbi algorithm over a sliding window of 7
  loci (Fragoso et al. 2016, *Genetics* 202:487–495, abstract).
- **Copy fidelity**: no mutation to speak of; only genotyping or sequencing error.

The consequence for the caller: in a cross, IBD segments are tens of centimorgans long and every
individual shares them with every other. Three reads a position is enough because the model needs
only to decide *which* of two founders is being copied over a window carrying hundreds of
informative sites. In a natural population the segments are short and shared with few others, so
the same coverage buys far less; that is why panel size and cohort size matter there and barely
matter in a cross (§2e, §4).

---

## 2. The families

For each family: what it does in plain words, what signal it uses, what it assumes, where it
breaks, then the mathematics.

### 2a. Parsimony — Clark's algorithm (1990)

**What it does.** Start with the individuals whose haplotypes are unambiguous (homozygous
everywhere, or heterozygous at one site only). Their haplotypes are "known". For every other
individual, look for a known haplotype that fits their genotype; if one fits, the complement is
the other haplotype, and it joins the known list. Repeat. The rule is that populations carry few
haplotypes, so a parsimonious explanation is likely the true one (Clark 1990, *Mol Biol Evol*
7:111–122; the abstract describes inferring alleles "*by taking the remaining sequence after
'subtracting off' the sequencing ladder of each known site*").

**Signal.** The small number of distinct haplotypes in a window — which is LD, used without a
model of it.

**Where it breaks.** It cannot start if no individual is unambiguous (Stephens, Smith & Donnelly
2001, *Am J Hum Genet* 68:978–989, Introduction). It gives no way to choose between several
compatible resolutions, and the answer depends on the order individuals are processed (Browning &
Browning 2011, *Nat Rev Genet* 12:703–714: "*when polymorphisms are not tightly linked there may be
several reasonable haplotype phase assignments ... and the method does not provide a means of
choosing between such assignments*"). On sequences with 60–100 segregating sites it mis-phased 42%
of individuals against 20% for PHASE (Stephens et al. 2001, Table 1).

### 2b. Expectation–maximisation over haplotype frequencies (1995)

**What it does.** Treat the population as having an unknown frequency for every possible
haplotype in the window. Given current frequencies, compute for each individual the probability of
each way its genotype could be split into two haplotypes (the E step, assuming Hardy–Weinberg
proportions); then re-estimate the frequencies from those probabilities (the M step); repeat until
they stop changing (Excoffier & Slatkin 1995, *Mol Biol Evol* 12:921–927).

**Signal.** Again the small number of common haplotypes, now with an explicit frequency for each.

**Where it breaks.** The list of possible haplotypes has 2ᴸ entries for L biallelic sites, so
"*the computing time grows exponentially with the number of polymorphic loci*" (Excoffier &
Slatkin 1995, abstract); in practice it is limited to about 10 markers (Browning & Browning 2011).
It has no notion of recombination or of two haplotypes being *similar*: a haplotype that differs at
one site from a common one gets no credit for the resemblance. It needs several starting points to
avoid local maxima (Excoffier & Slatkin 1995, abstract).

### 2c. The haplotype-copying model — PHASE, Li & Stephens, and everything since

**What it does.** Ask, for one haplotype at a time: given the other haplotypes already known, which
of them is this one copying at each position? The unknown "which template" is a hidden state that
changes along the chromosome at recombination breaks; the observed alleles are the template's
alleles, copied with occasional errors. This is a **hidden Markov model** (HMM): a chain of hidden
states, each depending only on the previous one, each emitting an observation.

PHASE (Stephens, Smith & Donnelly 2001) introduced the coalescent motivation: "*the next sampled
haplotype, h, being obtained by applying a random number of mutations, s, to a randomly chosen
existing haplotype*", and the update scheme that every self-panel method still uses: "*the algorithm
starts with an initial guess H(0) for H, repeatedly chooses an individual at random, and estimates
that individual's haplotypes under the assumption that all the other haplotypes are correctly
reconstructed*" (Methods). PHASE v2 added recombination through Li & Stephens' model and used it
to impute missing genotypes: genotype-imputation error fell from 0.044–0.062 to 0.027–0.038 when
recombination was modelled (Stephens & Scheet 2005, *Am J Hum Genet* 76:449–462, Table 3).

**The mathematics (Li & Stephens 2003, Appendix A).** With k template haplotypes h₁…h_k, a new
haplotype h_{k+1} is generated by a Markov chain X_j ∈ {1…k} over sites j:

- *Transition* (eq. A1). Between sites j and j+1 at physical distance d_j, the copied template
  stays the same with probability exp(−ρ_j d_j / k) and otherwise jumps to a template chosen
  uniformly among the k (including, possibly, the same one). ρ_j = 4*N*ₑ*c*_j, with *c*_j the
  per-base recombination rate. The division by k is the coalescent's statement that with more
  lineages present, any one of them is less likely to be the nearest.
- *Emission* (eq. A2). The copied allele is reproduced exactly with probability k/(k+θ̃) and
  "mutated" with probability θ̃/(k+θ̃), where θ̃ = (Σ_{m=1}^{n−1} 1/m)⁻¹ (eq. A3) is a
  Watterson-style estimate that makes the prior expected number of mutations per site equal to
  one.
- The probability of the whole haplotype is a sum over all copying paths, computed by the
  **forward algorithm** in time linear in the number of sites and linear in k (Appendix A,
  "Computation").
- The likelihood of the sample is the **product of approximate conditionals** (PAC):
  L_PAC(ρ) = π̂(h₁|ρ) π̂(h₂|h₁; ρ) … π̂(h_n|h₁…h_{n−1}; ρ) (eq. 3). It depends on the order the
  haplotypes are added; the authors average over 10 random orders (Fig. 3 legend).

Every modern implementation writes the transition as 1 − exp(−0.04 · *N*ₑ · d_cM / K), where the
0.04 is 4/(100 cM per Morgan): SHAPEIT5 `shapeit5/phase_common/src/objects/hmm_parameters.cpp:46`
(*N*ₑ default 15,000, `.../phaser/phaser_parameters.cpp:64`); Beagle 5.5 phasing
`beagle_src/src/phase/PhaseData.java:59` and imputation `beagle_src/src/imp/ImpData.java:217-220`
(*N*ₑ default 100,000, `beagle_src/src/main/Par.java:102`); GLIMPSE2
`GLIMPSE/phase/src/containers/conditioning_set.cpp:34,190-194` (*N*ₑ default 100,000). The
"mutation" emission became a mismatch probability: SHAPEIT5 hard-codes 10⁻⁴
(`hmm_parameters.cpp:28-29`); Beagle 5.5 starts from Li & Stephens' θ̃ — its initial value is
θ/(2(θ + K)) with θ = 1/(ln K + 0.5) (`beagle_src/src/main/Par.java:494-496`, the harmonic sum
approximated by a logarithm) — and then re-estimates it from the data during burn-in
(`Par.java:103` sets the default to a sentinel meaning "data-dependent";
`phase/ParamEstimates.java:81-89` returns the running estimate; `phase/HmmParamData.java:74-75`
builds the two emission values).

**Three ways to use the same HMM.** Once transitions and emissions are defined, three standard
computations give three different outputs:

- *Forward–backward* gives, at every site, the posterior probability of each template being
  copied. Summing template allele probabilities gives the posterior probability of each allele —
  the imputed dosage — for a haploid HMM, or of each genotype for a diploid one. This is what
  imputers report as GP and DS (§5).
- *Viterbi* gives the single most probable copying path. Useful for one best haplotype; used by
  SHAPEIT5 to produce its final phased output from the accumulated iterations
  (`shapeit5/phase_common/src/phaser/phaser_finalise.cpp:35` → `genotype::solve`, per the code
  reading: max-product forward and argmax traceback).
- *Stochastic traceback* draws a random path in proportion to its posterior probability. This is
  what a Gibbs sampler needs: SHAPEIT5 samples a diploid path forward or backward at random
  (`shapeit5/phase_common/src/objects/genotype/genotype_sweep.cpp:27-30`); GLIMPSE2 samples each
  haplotype from its posterior (`GLIMPSE/phase/src/caller/caller_algorithm.cpp:65,70`).

**"The reference panel enters as the set of hidden states" — concretely.** In IMPUTE v1 the study
individual's two haplotypes are modelled as a mosaic of the HapMap panel, with hidden state = a
*pair* of panel haplotypes, so the state count is N² for N panel haplotypes; "*the computational
burden ... grows quadratically with the number of haplotypes and linearly with the number of SNPs*"
(Howie, Donnelly & Marchini 2009, *PLoS Genet* 5:e1000529, Methods, describing IMPUTE v1;
Marchini et al. 2007, *Nat Genet* 39:906–913 — only its abstract could be opened). IMPUTE2 cut
this to N states by phasing the study sample first and then imputing each *haplotype* against the
panel — the **haploid HMM** — and further to k ≪ N by keeping only the k panel haplotypes nearest
(by Hamming distance) to the current haplotype guess, k = 40 or 80 in their runs (Howie et al.
2009, Methods). In every later tool, "the panel" is literally the list of haplotypes indexed by the
hidden state: SHAPEIT5 stores reference haplotypes as rows 2N onwards of the same bit-matrix as
the cohort (`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_reading.cpp:138-139`,
`containers/haplotype_set.cpp:43`); GLIMPSE2 dereferences `HvarRef`/`ShapRef` for the alleles of
state k at site l (`GLIMPSE/phase/src/containers/conditioning_set.cpp:157-168`, per the existing
code-facts note).

**Why state-space reduction was needed.** Cost per sample is (iterations) × (sites) × (states).
With N = 10⁵ panel haplotypes and 10⁶ sites, the plain O(N) per site is 10¹¹ operations per
sample per iteration. Panels of that size were the point of the Haplotype Reference Consortium
(64,976 haplotypes) and of biobanks, so every 2016–2023 tool is a different answer to "which few
hundred of the N templates are worth carrying at this position":

| Tool | Selection of conditioning haplotypes | Source |
|---|---|---|
| IMPUTE2 (2009) | k nearest by Hamming distance to the current haplotype guess; k = 40–80 | Howie et al. 2009, Methods |
| Beagle 4.1 (2016) | a target-specific set chosen by identity-by-state (IBS) segment matching; model restricted to genotyped markers, linear interpolation between them | Browning & Browning 2016, *Am J Hum Genet* 98:116–126, abstract (methods text not reachable); confirmed by Browning et al. 2018 |
| minimac3 (2016) | collapse haplotypes identical within a block into "unique haplotypes" (m3vcf); iterate only over distinct ones | Das et al. 2016, *Nat Genet* 48:1284–1287, Online Methods |
| Eagle2 (2016) | HapHedge: haplotype prefix trees rooted at markers spaced along the chromosome, built by PBWT; beam search keeps only the most probable diplotype paths | Loh, Danecek et al. 2016, *Nat Genet* 48:1443–1448, Methods; `Eagle/src/HapHedge.hpp:30-53`, beam widths `Eagle/src/EagleParams.cpp:138-140` (100 / 200 state pairs) |
| Beagle 5.0 (2018) | 1,600 "composite reference haplotypes", each a mosaic of IBS segments from many panel haplotypes, so one state carries a different template in different regions | Browning, Zhou & Browning 2018, *Am J Hum Genet* 103:338–348, Methods; `beagle_src/src/main/Par.java:92` |
| SHAPEIT4 (2019) | PBWT of all current haplotype estimates, queried every 8 variants for the P longest-matching neighbours; matches deduplicated into K distinct states | Delaneau et al. 2019, *Nat Commun* 10:5436, Methods |
| IMPUTE5 (2020) | PBWT of the panel at genotyped markers, re-selecting the best L (4 or 8) states every 0.02 cM; 10,000 → 1,000,000 haplotypes costs less than 2× the time | Rubinacci, Delaneau & Marchini 2020, *PLoS Genet* 16:e1009049, abstract and Methods |
| SHAPEIT5 (2023) | as SHAPEIT4 for common variants; for rare variants, neighbours restricted to carriers of the rare allele | Hofmeister et al. 2023, *Nat Genet* 55:1243–1249, Methods; defaults `shapeit5/phase_common/src/phaser/phaser_parameters.cpp:55-58` |
| GLIMPSE2 (2023) | pre-built, run-length-compressed PBWT of the panel stored on disk; `Kpbwt` = 2,000 states | Rubinacci et al. 2023, *Nat Genet* 55:1088–1090; `GLIMPSE/phase/src/caller/caller_parameters.cpp:67-70` |

The **PBWT** (positional Burrows–Wheeler transform) that most of these use is a way of keeping the
haplotypes sorted, at every site, by their sequence *backwards* from that site ("*order the
sequences ... so that their reversed prefixes are ordered*", Durbin 2014, *Bioinformatics*
30:1266–1272, §2). Haplotypes adjacent in that order share the longest recent stretch, so "the
nearest templates here" is a neighbourhood lookup rather than a search; all set-maximal matches
in a collection are found in O(NM) time (§2.3, Algorithm 4).

Two results say the reduction costs nothing in accuracy: "*All the genotype imputation methods that
we evaluated are based on the Li and Stephens probabilistic model and have essentially the same
imputation accuracy*" (Browning et al. 2018, Results), and IMPUTE5's comparison found the same
(Rubinacci et al. 2020, Results, Fig. 3). What differs is time: Beagle 5.0 was 3× / 12× / 43× /
533× faster than the fastest alternative at 10⁴ / 10⁵ / 10⁶ / 10⁷ reference samples (Browning
et al. 2018, Results).

### 2d. Haplotype-cluster models — fastPHASE, Beagle 3, STITCH

**What they do.** Instead of copying from the other sampled haplotypes, copy from K *fitted*
haplotypes that do not exist in the sample: each cluster is a probability of each allele at each
site, and each sampled haplotype is a mosaic of clusters. The clusters and the switching rates are
estimated from the sample by EM. This is the reference-free end: nothing outside the sample is
needed, and the state count is K (typically single digits to tens) rather than 2N.

*fastPHASE* (Scheet & Stephens 2006, *Am J Hum Genet* 78:629–644) lets "*cluster membership of
observed haplotypes ... change continuously along the genome according to a hidden Markov model*";
K ∈ {4…12} was tested and K = 8 "*seemed to perform reasonably well across a range of scenarios*";
fitting used 20 EM starts of up to 25 iterations each, and predictions averaged across the 20 fits
rather than taking the best one. Cost is O(nMK²) — quadratic in K because the diploid state is a
pair of clusters — and 15,532 SNPs × 60 individuals took about 5 hours against 720 for PHASE
(Table 6). Switch error was 0.055 against PHASE's 0.051 on CEPH HapMap data (Table 5).

*Beagle 3* (Browning & Browning 2007, *Am J Hum Genet* 81:1084–1097; Browning & Browning 2009,
*Am J Hum Genet* 84:210–223) dropped the fixed K: the "localized haplotype-cluster model" is a
directed acyclic graph whose edges are clusters, "*the number of clusters and the relationships
between clusters at different positions*" being "*determined largely by the data*". The cycle is:
build the graph from current phased haplotypes, sample new haplotypes for every individual from the
induced diploid HMM, rebuild; "*10 iterations gives good accuracy*". On 5,000 simulated individuals
the switch error was 0.05% against fastPHASE's 0.49% (2007, Table 3). Imputation with a panel of
60 / 100 / 300 / 600 / 1,200 reference individuals gave allele-frequency correlations of
0.9902 / 0.9944 / 0.9976 / 0.9982 / 0.9986 (2009, Results) — the first quantitative curve of
"more templates, better imputation".

*STITCH* (Davies et al. 2016, *Nat Genet* 48:965–969) is fastPHASE's model driven by reads. The
parameters are, in the code's names, `eHapsCurrent_tc` (K × sites: the probability that ancestral
haplotype k carries the alternate allele), `alphaMatCurrent_tc` (K × intervals: which haplotype a
switch lands on) and `sigmaCurrent_m` (per interval: the probability of *not* switching), all
initialised at random or uniform (`STITCH/STITCH/R/functions.R:1350-1353`) and re-estimated each
EM iteration from the forward–backward posteriors (`functions.R:629-631`, `:659-662`). The
switching rate is bounded by the user's `nGen` — generations since founding — times a per-base
rate: `sigma = exp(−nGen · rate · d)` with rate between `minRate` and `maxRate`, default
`expRate` 0.5 cM/Mb (`functions.R:119`, `:4462-4478`); the documentation suggests
`nGen = 4·Nₑ/K` when unknown (`functions.R:6`). The emission is per read: for a read base with
quality Q, the probability of that base given the template allele is 1−10^(−Q/10) if they match
and 10^(−Q/10)/3 otherwise (`STITCH/STITCH/src/haploid.cpp:179-186`); a read spanning several
sites contributes the product. Diploid mode has K² states (`STITCH/STITCH/src/diploid.cpp:128`,
`:182-189`); pseudo-haploid mode is linear in K by estimating for each read which parental
chromosome it came from (Davies et al. 2016). Default 40 EM iterations (`functions.R:115`).
Reported accuracy: outbred mice at 0.15× and N = 2,073 with K = 4, r² = 0.972 against arrays;
Han Chinese at 1.7× and N = 11,670 with K = 40, r² = 0.920; for the mice "*sample size above 500 has
little impact on performance*" (Davies et al. 2016, Results). An optional reference panel only
seeds `eHaps` (`functions.R:1373`, `get_and_initialize_from_reference`); STITCH does not need one.

**Where cluster models break.** K must be chosen, and it stands for the number of distinct
ancestral haplotypes in a window: K = 4 fits a mouse population 100 generations from two founders,
K = 40 an outbred human population (Davies et al. 2016). In a natural population with thousands of
distinct haplotypes over a region, K fitted haplotypes are averages over many real ones, so rare
haplotypes and rare alleles are imputed poorly — STITCH's authors describe it as best for
"*recently bottlenecked*" populations. In a biparental cross the model is exact with K = 2.

### 2e. Iterative self-panel phasing — the cohort as its own panel

**What it does.** Take the current haplotype estimates of all *other* individuals as the templates;
re-estimate this individual; move to the next; sweep the cohort many times. The templates improve
as the sweeps go on because they are being re-estimated too. This is PHASE's Gibbs sampler with
Li & Stephens' HMM as the conditional. "*No external panel needed*" means exactly this: the state
list of the HMM is built from the cohort's own haplotype matrix, and it is rebuilt every iteration.

- *MaCH* (Li et al. 2010, *Genet Epidemiol* 34:816–834): "*describes sampled chromosomes as mosaics
  of each other*"; "*a new pair of haplotypes is sampled for each individual in turn using a Hidden
  Markov Model that describes the haplotype pair as an imperfect mosaic of the other haplotypes*";
  crossover and error parameters re-estimated from the sampled paths each round; 20–100 rounds.
- *SHAPEIT4/5*: the PBWT is recomputed from scratch on the current estimates at every MCMC
  iteration (`shapeit5/phase_common/src/phaser/phaser_algorithm.cpp:133-142`: `H.select()`, phase,
  `H.updateHaplotypes(G)`, transpose); an individual's own two haplotypes are never selected as
  its templates (`shapeit5/phase_common/src/containers/ibd2_tracks.cpp:94`). Schedule
  `5b,1p,1b,1p,1b,1p,5m` — five burn-in, pruning and burn-in alternating, five main iterations
  whose samples are accumulated, 15 in all (`phaser_parameters.cpp:49`; Delaneau et al. 2019,
  Methods). The HMM is diploid but compact: consecutive heterozygous sites are grouped into
  segments of at most 3 hets, giving 8 haplotype hypotheses per segment carried in one AVX2
  register (`shapeit5/phase_common/src/objects/genotype/genotype_build.cpp:52`,
  `shapeit5/common/src/utils/otools.h:79`).
- *Beagle 5.2+*: not a sampler. Each iteration builds composite reference haplotypes "*from the
  estimated haplotypes in other individuals at the start of the iteration*" and updates the phase
  of each still-uncertain heterozygote by forward–backward probabilities, marking the most
  confident as finished ("progressive phasing"); 3 burn-in and 12 phasing iterations by default
  (Browning, Tian, Zhou & Browning 2021, *Am J Hum Genet* 108:1880–1890, Methods;
  `beagle_src/src/main/Par.java:83-84`; 280 states per haplotype, `Par.java:86`,
  `phase/HmmParamData.java:64-69`). In the 5.5 sources the state list is the other samples'
  current estimates combined with the reference panel only if one is given
  (`beagle_src/src/phase/CodedSteps.java:56-61`); each heterozygote takes the more probable of
  its two orientations (`phase/PhaseBaum2.java:327-340`), and a PBWT sweep picks one
  identical-by-state haplotype per target haplotype per step to feed the composite states
  (`phase/PbwtPhaseIbs.java:170-177`). Rare variants (non-major allele frequency < 0.002,
  `Par.java:88`) are phased in a second stage by imputing them onto the phased common-variant
  scaffold (Browning et al. 2021, Methods) — SHAPEIT5's `phase_rare` does the same, restricting the
  templates to carriers of the rare allele (Hofmeister et al. 2023, Methods).
- *Eagle2* without `--vcfRef`: the haplotype library is the cohort's own 2N haplotypes
  (`Eagle/src/Eagle.cpp:1427`: `Nhaps = 2*(Nref==0 ? N : Nref)`), with 10,000 conditioning
  haplotypes (`EagleParams.cpp:109`), an expected copying length of 2 cM (`:113`), and 2–3 PBWT
  rounds at increasing resolution (`EagleMain.cpp:440-455`).

**What it needs from cohort size, and what happens at N = 1, 10, 100, 1,000.** The method's only
templates are the other 2(N−1) haplotypes, so its quality is the quality of the best match among
them. The measured curve is steep:

| Cohort | Evidence |
|---|---|
| N = 1 | No template exists. SHAPEIT5 exits with "*Less than 50 samples is not enough to get reliable phasing*" (`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_scaning.cpp:46-48`); with a reference panel of ≥ 49 samples it runs and the panel supplies every template. Beagle 5.5 exits with "*ERROR: there is only one sample*" once the state count `min(nHaps − 2, 280)` reaches zero (`beagle_src/src/phase/BasicPhaseStates.java:330-332`), where `nHaps` counts target plus reference haplotypes. Without a panel the only signals are reads and pedigree. |
| N = 10 | 18 templates. In a natural population, two random chromosomes share an IBD segment of useful length only rarely, so most windows have no close template and phase is guessed from allele frequencies. SHAPEIT5 refuses (< 50). In a biparental cross, 18 templates already contain both founders many times over, and the model is fully determined — this is where the two population cases part company (§4). |
| N = 100 | Stephens et al. 2001 phased n = 50 microsatellite samples with PHASE (Fig. 4); fastPHASE and PHASE benchmarks used 60 CEPH individuals (Scheet & Stephens 2006, Table 6). Browning & Browning 2011 (Fig. 2) report that at these sizes MACH was the most accurate and BEAGLE only overtook it above about 1,000 individuals. STITCH's fitted-cluster model, which does not need pairwise matches, saw no gain above 500 mice at 0.15× (Davies et al. 2016). |
| N = 1,000–5,000 | Eagle2 cohort-only, UK Biobank arrays: switch error 2.10% at N = 5,000 (SHAPEIT2 2.27%) (Loh, Danecek et al. 2016, Fig. 5b, Suppl. Table 7). Beagle 3 imputation with 1,200 reference individuals: allele-frequency correlation 0.9986 vs 0.9902 with 60 (Browning & Browning 2009). |
| N = 15,000–150,000 | Eagle2 cohort-only: 1.69% at 15,000, 1.15% at 50,000, 0.35% at 150,000 (same source). Eagle1's long-range phasing needs IBD tracts > 4 cM, which exist in numbers only in very large cohorts; its authors "*recommend SHAPEIT2*" at N ≈ 15,000 (Loh, Palamara & Price 2016, *Nat Genet* 48:811–816, Discussion). At 400,000: SHAPEIT4 0.117%, Beagle5 0.125%, Eagle2 0.178% (Delaneau et al. 2019, Fig. 2a). SHAPEIT5 docs recommend `phase_rare` only above 2,000 samples (`shapeit5/docs/docs/documentation/phase_rare.md:20-21`). |

The reason the curve keeps falling is §1.1: the nearest common ancestor of a haplotype and its
closest cohort neighbour gets more recent as N grows, so shared segments get longer and the mosaic
needs fewer, better-supported pieces. The same reasoning says a reference panel helps *small*
cohorts most: Eagle2 phased small European cohorts at 1.36% switch error with the 64,976-haplotype
HRC panel against 3.52% with 1000 Genomes, and "*For very large cohorts (substantially larger than
the reference size), we expect that reference-based phasing will achieve only marginal gains in
accuracy over cohort-based phasing*" (Loh, Danecek et al. 2016, Results and Discussion).

One more detail of SHAPEIT5 matters for a small cohort: with neither `--pbwt-depth` nor
`--pbwt-modulo` given, both are set from N — `depth = max(min(round(9 − log10 N), 8), 2)`,
`modulo = max(min((ln N − ln 50 + 1)·0.01, 0.15), 0.005)` — so N = 50 gets depth 7 and 0.01 cM
(`shapeit5/phase_common/src/phaser/phaser_initialise.cpp:89-93`). The documented "depth 4,
modulo 0.1" is not the behaviour at small N.

### 2f. Read-backed phasing (pointer to sub-report 05)

Reads are a different signal: a read or read pair that covers two heterozygous sites shows which
alleles sit on one molecule, with no population model at all. Its reach is the insert size. In the
statistical methods above it enters as a *scaffold*: a set of heterozygotes whose relative phase
is fixed before the HMM runs. SHAPEIT5 accepts one through `--scaffold`, flagging those sites and
constraining the segment's hypotheses to one orientation
(`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_reading.cpp:197-201`,
`objects/genotype/genotype_build.cpp:136`); Eagle2 reads the VCF `PS` field for the same purpose
(`Eagle/src/SyncedVcfData.cpp:360`); QUILT and STITCH consume the reads themselves (§3). Delaneau
et al. 2019 report that phase sets from sequencing reads cut SHAPEIT4's error on GIAB data by 49%
when combined with a UK Biobank scaffold (Results). Sub-report 05 covers the read-based algorithms
(WhatsHap, HapCUT2) and what ng's chain ids can supply.

### 2g. Pedigree and Mendelian transmission (pointer to sub-report 04)

A parent–offspring trio fixes the phase of every offspring heterozygote at which at least one
parent is homozygous, by inheritance alone; duos fix fewer. Long-range phasing (Kong et al. 2008,
*Nat Genet* 40:1068–1075) extends the idea to distant relatives: chromosomes that agree over a long
stretch at homozygous sites are IBD there, and the IBD partner's phase resolves the target's. In the
HMM methods, pedigree information enters the same way as a read scaffold: SHAPEIT5's `--pedigree`
fixes trio- and duo-resolved heterozygotes as scaffolded sites
(`shapeit5/phase_common/src/containers/genotype_set.cpp:104-145`,
`objects/genotype/genotype_mendel.cpp:31-118`); Eagle1 makes its initial calls from > 4 cM IBD
tracts and then refines them with two HMM passes (Loh, Palamara & Price 2016, Results). In a
segregant population the pedigree *is* the model: every individual descends from the same founders
through a known number of meioses. Sub-report 04 covers this.

---

## 3. Genotype likelihoods in the emission: calling and imputation as one model

### 3.1 What a genotype likelihood is

For one sample at one site, the genotype likelihood is the probability of the observed reads given
each possible genotype: L(g) = Pr(reads | g). With k reads, of which the first l match the
reference, and per-base error ε_j, independence gives L(g) as a product over reads of the
probability of each base under g (Li 2011, *Bioinformatics* 27:2987–2993, eq. 2). VCF carries it
as `GL` (log10) or `PL` (phred-scaled, normalised to a minimum of 0). It is not a genotype and not
a posterior: three reads all showing the alternate allele give a likelihood ratio of about 1,000:1
for hom-alt over het, but no evidence at all between het and hom-ref beyond "not hom-ref".

### 3.2 Replacing the hard genotype

A hard-genotype phaser's emission says: the template's allele matches the called allele with
probability 1−ε, else ε. A likelihood-driven emission says: the probability of the reads given
that this haplotype carries allele a *and* the other haplotype carries whatever it currently
carries. Nothing else in the HMM changes. That is the whole trick, and it was first done for array
intensities: "*the emission probability of the individual's observed signal intensities S for an
HMM state labeled with genotype g is the genotype likelihood P(S|G = g)*" (Browning & Yu 2009,
*Am J Hum Genet* 85:847–861, Methods). Beagle 4.0/4.1 exposed it for sequencing as the `gl=`
argument ("*specifies a VCF file containing a GL or PL (genotype likelihood) format field*",
Beagle 4.0 manual, 3 March 2015, §2.1, text extracted from the PDF); Beagle 5.0 removed it
("*Input VCF file must contain genotype (GT) fields (GL/PL fields are ignored)*", Beagle 5 release
notes; the Beagle 5.5 sources contain no `GL`/`PL` reader). Pasaniuc et al. 2012 (*Nat Genet*
44:631–635) used that Beagle mode to show that 0.1–0.5× sequencing plus imputation with a
1000 Genomes panel matched or beat arrays for association power at equal cost.

### 3.3 The diploid HMM against the one-haplotype-at-a-time trick

With likelihoods over *genotypes*, the natural HMM is diploid: hidden state = the pair of
templates copied by the two haplotypes, emission = L(g) for the genotype those two templates
imply. With K templates that is K² states per site — fastPHASE's O(K²) and IMPUTE v1's O(N²).

The trick that every current low-coverage imputer uses is to fix one haplotype at its current
estimate and run a *haploid* HMM for the other, then swap, and iterate. GLIMPSE2 does it
literally: `makeHaplotypeLikelihoods` takes the three normalised genotype likelihoods and, given the
other haplotype's allele c at this site, returns the two-allele likelihood proportional to
(L(c+0), L(c+1)) — i.e. (L₀₀, L₀₁) when the other haplotype is reference, (L₀₁, L₁₁) when it is
alternate — floored at a minimum (`GLIMPSE/phase/src/objects/genotype.cpp:94-115`). The emission
is then p₀ ∝ HL₀(1−e) + HL₁e, p₁ ∝ HL₀e + HL₁(1−e) with e the imputation error rate
(`models/imputation_hmm.cpp:61-71`). The per-iteration loop is: build HL for haplotype 0 given
haplotype 1, forward–backward, sample haplotype 0; build HL for haplotype 1 given the new
haplotype 0, forward–backward, sample haplotype 1; re-phase the pair with a SHAPEIT-style diploid
segment HMM (`caller/caller_algorithm.cpp:58-72`). Posteriors are accumulated over the main
iterations (1 init, 5 burn-in, 15 main by default, per the code-facts note). Rubinacci et al. 2021
(*Nat Genet* 53:120–126; the published text is paywalled, so the wording below is from the bioRxiv
version, 10.1101/2020.04.14.040329) describe GLIMPSE1's Gibbs sampler this way — "*imputes the two
target haplotypes in turn using the resulting likelihoods as additional layers of emission
probabilities in a haploid version of the Li and Stephens*" model — and draw conditioning
haplotypes from "*the reference panel and from other target individuals*": the 2P panel
haplotypes plus "*the 2Q − 2 haplotypes previously estimated for the other target individuals*",
reduced by PBWT to K = 1,000 states (Online Methods). GLIMPSE2 dropped the target haplotypes from
the conditioning set (code-facts note §4: PBWT arrays are sized to the reference; a panel is
mandatory, `GLIMPSE/phase/src/caller/caller_parameters.cpp:133-134`).

QUILT (Davies et al. 2021, *Nat Genet* 53:1104–1111) pushes the same trick down to the reads:
each read is given a label, maternal or paternal; given the labels, each haplotype is imputed by a
haploid HMM from its reads alone; given the two haplotypes, the labels are re-sampled by Gibbs
(`QUILT/QUILT/src/gibbs-nipt.cpp:733`, `sample_reads_in_grid`; 7 Gibbs samples × 3 seek
iterations, 600-haplotype subset, `QUILT/QUILT/R/quilt.R:110-113`). The read-level emission keeps
the phase information inside a read, which is why QUILT wins at very low coverage and with long
reads: rare-variant r² 0.416 against GLIMPSE's 0.328 at 0.1× (Davies et al. 2021), and with
Nanopore reads at 1× QUILT2 reaches common-variant r² 0.937 against GLIMPSE2's 0.695 (Li,
Albrechtsen & Davies 2025, *Nat Commun* 17:524, Fig. 2B). With short reads QUILT2 is the more
accurate below 0.5× and GLIMPSE2 above 2× for rare variants (same paper). GLIMPSE2 in turn is far cheaper: £0.08 per genome on the UK
Biobank platform against £1.11 for GLIMPSE1 and £242.80 for QUILT 1.0.4 (Rubinacci et al. 2023).
QUILT requires a reference panel (`QUILT/QUILT/R/quilt-prepare-reference.R:148-167`: the only
accepted inputs are a reference VCF or an IMPUTE-format haplotype/legend pair).

STITCH is the panel-free member of this group: the same read-level emission (§2d) but the
templates are the K fitted ancestral haplotypes rather than a panel.

### 3.4 Why this makes calling and imputation one model

The output of the HMM is, per sample and site, a posterior over genotypes that combines the
sample's own reads (through the likelihood) with what the templates say the haplotype should carry
(through the copying path). That is a genotype call with a population-and-linkage prior instead of
a site-wise allele-frequency prior. At high coverage the likelihood dominates and the result is the
ordinary call; at 1× the prior dominates and the result is the imputation; the model does not have
a switch between the two. For ng this is the argument for feeding likelihoods rather than calls
into any such module (the brief's fact 2).

---

## 4. The continuum, answered

Two axes. On the first, the *number and origin of templates*; on the second, the *evidence fed
to the emission*.

| | Hard genotypes (arrays, deep sequencing calls) | Genotype likelihoods (low coverage) | Reads (with base qualities) |
|---|---|---|---|
| **2 founders** (biparental cross) | TIGER-style breakpoint callers; generic HMM with 2 states | LB-Impute (Fragoso 2016): HMM over parental origin, emission from read counts | STITCH with K = 2 |
| **K fitted clusters** (K ≈ 4–40) | fastPHASE (Scheet & Stephens 2006) | — | STITCH (Davies 2016) |
| **The cohort itself** (2N−2 templates, re-estimated each sweep) | PHASE, MaCH, SHAPEIT1–5, Beagle 3–5, Eagle1/2 | Beagle 4.1 `gl=` without a panel (removed in 5.0); GLIMPSE1 (targets in the conditioning set) | — (none found) |
| **External panel** (10³–10⁵ haplotypes, phased) | IMPUTE2/4/5, minimac3/4, Beagle 4.1–5.5 imputation, Eagle2 `--vcfRef` | GLIMPSE1/2; Beagle 4.1 `gl=` with `ref=` | QUILT, QUILT2 |

Reading the table:

- **Down a column, the only change is the template set and its clock.** The HMM state is "which
  template"; the transition is 1 − exp(−rate × distance); the emission is the template's allele
  against the evidence. Two founders: rate = 1/meiosis × generations (map expansion for RILs),
  templates known and phased. K clusters: rate fitted, templates fitted by EM. Cohort: rate =
  4*N*ₑ/(2N) per Morgan, templates re-estimated each sweep. Panel: same rate with the panel's
  haplotype count in the denominator, templates fixed and phased. The GLIMPSE2 code-facts note
  (§9) lists what would have to change to run it panel-free, and the list is bookkeeping — where
  the haplotype matrix comes from, rebuilding the PBWT each iteration, the denominator of the
  recombination rate — not a different model.
- **Along a row, the only change is the emission.** Hard genotype → mismatch probability ε;
  likelihood → the two-allele likelihood given the other haplotype (§3.3); reads → product over
  bases of the quality-derived probability, which additionally carries within-read phase.
- **What is genuinely different between the ends** is not the mathematics but the *support* the
  data gives the model. At the two-founder end the templates are known, every individual copies
  them over tens of centimorgans, and a handful of individuals at 1–3× is enough: the posterior over
  "which founder here" is decided by hundreds of sites per segment. At the natural-population end
  the templates are unknown, each is shared with few others over short segments, and the posterior
  is decided by tens of sites per segment; that is why accuracy there is a function of N (§2e
  table), of coverage (GLIMPSE1 with the 54,330-haplotype HRC panel, imputing downsampled
  1000 Genomes samples: r² > 0.7 at MAF > 5% at 0.1×, r² > 0.9 at MAF > 5% at 0.3×, r² = 0.8 at
  MAF 0.1% at 1×, r² > 0.95 at MAF 0.1% at 8× — Rubinacci et al. 2021, Results and Fig. 2, quoted
  from the bioRxiv text; the published figure was not reachable) and of panel size (GLIMPSE2 with
  the 280,238-haplotype UK Biobank panel reports r² = 0.892 at 0.1× and 0.927 at 1× for variants
  at MAF 0.01%, Rubinacci et al. 2023, Results). Rare variants are the sharpest case: they exist in few templates, so
  only a panel (or a huge cohort) can supply one, and both SHAPEIT5 and Beagle 5.2 treat them in a
  separate stage for that reason.
- **What is only a change of parameterisation**: whether K is fixed (fastPHASE), data-driven
  (Beagle 3), or equal to 2N (SHAPEIT); whether the recombination clock is a genetic map times
  *N*ₑ or a pedigree depth; whether updates are sampled (SHAPEIT, GLIMPSE) or maximised (Beagle
  5.2, STITCH's EM); whether states are selected by Hamming distance, IBS segments or PBWT.

For the owner's two cases the table says: the few-haplotype case sits in the top-left cells and
is well served by two or K templates with a pedigree clock; the many-haplotype case sits in the
bottom two rows and needs either a cohort large enough to be its own panel or an external panel.
The deep-sequenced individuals of a germplasm collection are, once phased, a panel for the
shallow ones — the "cohort as its own panel" row with the deep samples' haplotypes fixed early.

---

## 5. Output quantities and how they are computed

- **Phased genotype (`GT` with `|`).** From the sampled or Viterbi copying path: the allele on
  haplotype 1 is the allele of the template it copies (or the sampled allele given the posterior).
  SHAPEIT5 `phase_common` writes phased GT and nothing else per genotype
  (`shapeit5/phase_common/src/io/haplotype_writer.cpp:79`); `phase_rare` adds `PP`, a per-variant
  phasing confidence for rare heterozygotes (`shapeit5/phase_rare/src/io/haplotype_writer.cpp:79-80`).
  Beagle 5.2 uses the diplotype with the larger forward–backward probability (Browning et al. 2021).
- **Genotype posterior (`GP`).** The three probabilities P(0/0), P(0/1), P(1/1) from
  forward–backward. In the haploid-per-haplotype scheme they are products of the two haplotypes'
  allele posteriors accumulated over main iterations (GLIMPSE2, code-facts note §5). Beagle 5
  writes GP only with `gp=true` (`beagle_src/src/main/Par.java:97-98`; its imputed FORMAT is
  `GT:DS` by default, `imp/ImputedRecBuilder.java:259-265`); QUILT and STITCH always
  (`QUILT/QUILT/R/writers.R:240`).
- **Allele dosage (`DS`).** The expected count of the alternate allele, DS = GP₁ + 2·GP₂, a value
  in [0, 2]; MaCH introduced it as "*the expected number of copies of the minor allele at each
  position (a fractional value between 0 and 2)*" (Li et al. 2010, Methods). It is the quantity
  association tests use because it carries uncertainty without a hard call.
- **Imputation quality per site.** Three related formulas, all "how much of the variance a
  perfectly imputed site would have does this site show":
  - *IMPUTE `INFO`* = 1 − Σᵢ(fᵢ − eᵢ²) / (2Nθ̂(1−θ̂)), with eᵢ = GP₁ + 2GP₂ the dosage,
    fᵢ = GP₁ + 4GP₂, θ̂ the mean dosage / 2 (Marchini & Howie 2010, *Nat Rev Genet* 11:499–511,
    Supplementary S3 eq. 16). GLIMPSE2 (`GLIMPSE/phase/src/io/genotype_writer.cpp:184`) and
    STITCH (`STITCH/STITCH/R/writers.R:224-227`) implement it verbatim, clipped at 0.
  - *MaCH r̂²* = Var(dosage) / (2p(1−p)) (Li et al. 2010, Appendix), the variance of the dosages
    across samples over the binomial variance at Hardy–Weinberg — it "can go above 1"
    (Marchini & Howie 2010, S3).
  - *Beagle `DR2`* ("dosage r²") = (Σp² − (Σp)²/n) / (Σp − (Σp)²/n) over the n haplotype allele
    probabilities p — the variance of the imputed allele probabilities over the variance they would
    have if all were 0 or 1 (`beagle_src/src/imp/ImputedRecBuilder.java:317-329`; header
    `vcf/VcfWriter.java:42`).
  - *Empirical r²*: the squared Pearson correlation between imputed dosage and the true genotype
    at masked sites (Das et al. 2016, Online Methods). This is what accuracy tables report; the
    three estimates above are what a run reports without truth.
- **Phase sets (`PS`).** A VCF integer naming the block within which phase is consistent; blocks
  from different phase sets have no known relative orientation. Statistical phasers do not need it
  on output because they phase whole chromosomes (SHAPEIT5 and GLIMPSE2 `ligate` write GT only,
  no PS — code-facts note §7); read-based phasers must, because their blocks end where reads run
  out (sub-report 05). On input, Eagle2 reads PS as a scaffold (`Eagle/src/SyncedVcfData.cpp:360`).
- **Switch error rate.** The benchmark for phasing: "*the number of switches required [to obtain
  the true phase] divided by the number of opportunities for switch error, which is the number of
  heterozygote markers in the individual's genotype minus 1*" (Browning & Browning 2011, Box 2);
  Eagle's version divides "*the number of phase mismatches at consecutive trio-phased SNPs by the
  total number of trio-phased heterozygous SNPs minus 1*" (Loh, Palamara & Price 2016, Online
  Methods). It counts *changes* of orientation, so one wrong site in the middle of a block costs
  two switches. Its complement, "switch accuracy", traces to Lin, Cutler, Zwick & Chakravarti 2002
  (*Am J Hum Genet* 71:1129–1137; full text not reachable, attribution via Andrés et al. 2007,
  *Genet Epidemiol* 31:659–671).

---

## 6. Questions to discuss

1. **Is one HMM engine with a pluggable template set enough for both population cases?**
   *Recommendation: yes for the model; no for the surrounding machinery.* The transition and
   emission code is identical across the range (§4). What differs is (a) where the templates come
   from and how often they are rebuilt — fixed founders once, versus the cohort's own haplotypes
   every sweep with a PBWT re-sort; (b) the clock — pedigree generations versus 4*N*ₑ/K per Morgan;
   (c) the state count — 2 versus hundreds after selection. The engine should take a template
   matrix, a per-interval transition vector and an emission callback, and know nothing else. The
   template-provider and the clock are two small pluggable pieces; the state-selection step
   (PBWT) is the one large component that only the cohort/panel case needs. Trade-off: a generic
   engine over K states costs 8×K floats per site per window where a two-state engine costs 16;
   irrelevant at K = 2, so the generic engine loses nothing there.

2. **Diploid HMM or the one-haplotype-at-a-time trick?** *Recommendation: the trick, everywhere.*
   It is what GLIMPSE2, QUILT, Beagle 5 and IMPUTE5 all do, it is linear in K, and it consumes
   genotype likelihoods with a three-line transformation (§3.3). The price is iteration: the pair
   must be re-estimated several times (GLIMPSE2: 1 + 5 + 15), and the phase between the two
   haplotypes needs its own small diploid pass (GLIMPSE2's `rephaseHaplotypes`, SHAPEIT's segment
   HMM). For the two-founder case a full diploid HMM has only 4 states and needs no iteration; it
   could be a special case of the same engine with K = 2 and diploid state = pair.

3. **Which evidence to feed: hard calls, likelihoods, or reads?** *Recommendation: likelihoods
   as the baseline, reads as an option for the few-haplotype case.* ng already has per-sample
   genotype likelihoods in memory during calling (brief, fact 2). Reads add within-read phase,
   which QUILT shows matters below 0.5× and with long reads; ng's chain ids would supply exactly
   that (brief, fact 1) but do not survive into calling today (ng facts note). The cost of the read
   route is a retained per-sample, per-allele chain-id list over a window at het sites.

4. **What does the module do at N = 1?** *Recommendation: emit no imputation and no phase unless
   a panel, a pedigree or reads are supplied; never fit a silent template set from one sample.*
   Every self-panel tool either refuses (SHAPEIT5 below 50) or degrades to allele-frequency
   guessing. This matches the project rule that a method needing a cohort must say what it does at
   one sample.

5. **Should the deep-sequenced individuals of a collection be phased first and frozen as a
   panel for the shallow ones?** *Recommendation: yes, as a two-pass design.* Pass one: SHAPEIT-style
   self-panel phasing of the ≥ 20× samples from hard calls. Pass two: GLIMPSE-style
   likelihood-driven imputation of the shallow samples against those haplotypes, optionally adding
   the shallow samples' own current estimates to the template set as GLIMPSE1 did. Trade-off: the
   panel is only as diverse as the deep subset; rare alleles carried only by shallow samples cannot
   be imputed and should be reported as such rather than pulled to reference.

6. **Which recombination clock in the segregant case, and who supplies it?** *Recommendation:
   take the cross design (F2 / RIL-selfing / RIL-sib / BC / DH, generation count) as input and
   derive the per-interval rate from physical distance × a nominal cM/Mb × the design's map
   expansion, then let EM refine it per interval as STITCH does (`sigmaCurrent_m`, bounded by
   `nGen`).* Trade-off: a wrong design declaration biases breakpoint counts; the EM refinement
   limits but does not remove the damage. A genetic map, if one exists, replaces the nominal
   cM/Mb.

7. **How much of the state-selection machinery is worth building for cohorts of hundreds to a few
   thousand?** *Recommendation: start without PBWT selection.* At N ≤ 1,000 the full 2N template
   set is ≤ 2,000 states, which GLIMPSE2 already runs with (`Kpbwt` = 2,000) at a cost of
   O(iterations × N × K × sites). PBWT selection is what makes 10⁵ haplotypes affordable; it can
   be added when a cohort that size is a real target. Trade-off: without selection, cost is
   quadratic in N; at N = 5,000 (10⁴ states) it is already ten times GLIMPSE2's default work per
   sample.

8. **Where in ng's pipeline does the module sit?** The ng facts note settles that it is a buffered
   stage over a window of loci or a second pass over the written cohort, never a per-record filter.
   The open question is whether the first pass should write likelihoods (`PL`) so that a second
   pass can run without re-reading alignments. *Recommendation: yes; it is the cheapest way to keep
   calling and imputation as one model without coupling their execution.*

---

## Glossary

Genetics terms (brief):

- **Identity by descent (IBD).** Two chromosome segments copied, without recombination between
  them, from the same ancestral chromosome. **Identity by state (IBS)** means only that the
  observed alleles agree, whatever the reason.
- **Linkage disequilibrium (LD).** Correlation between alleles at different sites in a population;
  the population-level consequence of shared IBD segments.
- **Map expansion.** In recombinant inbred lines, repeated generations of recombination make
  two loci appear further apart than in one meiosis: about 2× under selfing, 4× under sib mating
  (Haldane & Waddington 1931; Crow 2007).
- **Scaffold.** A set of heterozygous sites whose relative phase is already known (from reads,
  a pedigree, or a previous run) and is held fixed while the rest is estimated.

Statistics and software terms (fuller):

- **Hidden Markov model (HMM).** A model with a hidden state at each position along a sequence,
  where each state depends only on the previous one (the *transition*) and each state produces the
  observation at that position with some probability (the *emission*). Here the hidden state is
  "which template haplotype is being copied", the transition is recombination, and the emission is
  the template's allele reproduced with error.
- **Forward–backward algorithm.** Computes, for every position, the probability of each hidden
  state given *all* observations, by one pass left-to-right (forward) and one right-to-left
  (backward). Cost is linear in sites and in states for the copying model. Gives posteriors.
- **Viterbi algorithm.** Finds the single most probable sequence of hidden states. Same cost.
  Gives one best path, no uncertainty.
- **Stochastic traceback / sampling a path.** Draws one state sequence at random in proportion to
  its posterior probability. What a Gibbs sampler uses so that repeated draws explore the
  uncertainty rather than always returning the same answer.
- **Gibbs sampler / MCMC.** An iterative scheme that updates one unknown at a time (here: one
  individual's haplotypes) by drawing it from its distribution given the current values of all the
  others, and sweeps repeatedly. Early sweeps ("burn-in") are discarded; later ones are averaged.
- **EM (expectation–maximisation).** An iterative scheme for fitting parameters when part of the
  data is hidden: compute the expected hidden values given current parameters (E), then re-fit the
  parameters as if those expectations were observed (M). Deterministic; converges to a local
  optimum; needs several starts.
- **PAC likelihood.** "Product of approximate conditionals": the probability of a set of
  haplotypes written as a chain — the first unconditional, each next one given the previous ones
  via the copying model. Approximate because the true conditional under the coalescent is
  intractable, and order-dependent for the same reason.
- **PBWT (positional Burrows–Wheeler transform).** A data structure that keeps a set of haplotypes
  sorted, at every site, by their sequence running backwards from that site. Haplotypes adjacent in
  the order share the longest recent match, so "the closest templates here" is a neighbourhood
  lookup.
- **Conditioning haplotypes / states.** The subset of templates the HMM is allowed to copy from
  for a given sample in a given region; "state selection" is choosing that subset.
- **Composite reference haplotype.** Beagle 5's device: one HMM state that copies different panel
  haplotypes in different regions, stitched from IBS segments, so 1,600 states can carry segments
  from a panel of millions.
- **Genotype likelihood (GL / PL).** Pr(reads | genotype) for each genotype, from base qualities;
  `PL` is −10·log10 of it, shifted so the best is 0. Evidence, not a call.
- **Genotype posterior (GP).** Pr(genotype | reads, model): the likelihood combined with the prior
  the model supplies (here, the copying model over templates).
- **Dosage (DS).** Expected alternate-allele count, GP₁ + 2·GP₂.
- **Imputation r² / INFO / DR2.** Per-site estimates of how confidently a site was imputed,
  each a ratio of observed to expected dosage variance (§5). *Empirical* r² needs truth data.
- **Switch error rate.** Number of orientation changes needed to turn the inferred phase into the
  true phase, divided by (heterozygous sites − 1).
- **Phase set (PS).** A VCF label for a block of consistently phased sites; different blocks have
  unknown relative orientation.

---

## References

Papers (all opened on 2026-09-12 unless marked; "full text not reachable" means only the abstract
or a secondary source could be opened):

- Andrés AM et al. 2007. Understanding the accuracy of statistical haplotype inference with sequence data of known phase. *Genet Epidemiol* 31(7):659–671. doi:10.1002/gepi.20185.
- Browning BL, Browning SR. 2009. A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals. *Am J Hum Genet* 84(2):210–223. doi:10.1016/j.ajhg.2009.01.005.
- Browning BL, Browning SR. 2016. Genotype imputation with millions of reference samples. *Am J Hum Genet* 98(1):116–126. doi:10.1016/j.ajhg.2015.11.020. (Methods text not reachable; abstract only.)
- Browning BL, Yu Z. 2009. Simultaneous genotype calling and haplotype phasing improves genotype accuracy and reduces false-positive associations for genome-wide association studies. *Am J Hum Genet* 85(6):847–861. doi:10.1016/j.ajhg.2009.11.004.
- Browning BL, Zhou Y, Browning SR. 2018. A one-penny imputed genome from next-generation reference panels. *Am J Hum Genet* 103(3):338–348. doi:10.1016/j.ajhg.2018.07.015.
- Browning BL, Tian X, Zhou Y, Browning SR. 2021. Fast two-stage phasing of large-scale sequence data. *Am J Hum Genet* 108(10):1880–1890. doi:10.1016/j.ajhg.2021.08.005.
- Browning SR, Browning BL. 2007. Rapid and accurate haplotype phasing and missing-data inference for whole-genome association studies by use of localized haplotype clustering. *Am J Hum Genet* 81(5):1084–1097. doi:10.1086/521987.
- Browning SR, Browning BL. 2011. Haplotype phasing: existing methods and new developments. *Nat Rev Genet* 12(10):703–714. doi:10.1038/nrg3054.
- Browning SR, Browning BL. 2012. Identity by descent between distant relatives: detection and applications. *Annu Rev Genet* 46:617–633. doi:10.1146/annurev-genet-110711-155534.
- Clark AG. 1990. Inference of haplotypes from PCR-amplified samples of diploid populations. *Mol Biol Evol* 7(2):111–122. doi:10.1093/oxfordjournals.molbev.a040591. (Abstract only.)
- Crow JF. 2007. Haldane, Bailey, Taylor and recombinant-inbred lines. *Genetics* 176(2):729–732.
- Das S et al. 2016. Next-generation genotype imputation service and methods. *Nat Genet* 48(10):1284–1287. doi:10.1038/ng.3656.
- Davies RW, Flint J, Myers S, Mott R. 2016. Rapid genotype imputation from sequence without reference panels. *Nat Genet* 48(8):965–969. doi:10.1038/ng.3594.
- Davies RW et al. 2021. Rapid genotype imputation from sequence with reference panels. *Nat Genet* 53(7):1104–1111. doi:10.1038/s41588-021-00877-0.
- Delaneau O, Marchini J, Zagury J-F. 2012. A linear complexity phasing method for thousands of genomes. *Nat Methods* 9(2):179–181. doi:10.1038/nmeth.1785. (SHAPEIT v1; abstract only. SHAPEIT2 is Delaneau, Zagury & Marchini 2013, *Nat Methods* 10(1):5–6, doi:10.1038/nmeth.2307.)
- Delaneau O, Zagury J-F, Robinson MR, Marchini JL, Dermitzakis ET. 2019. Accurate, scalable and integrative haplotype estimation. *Nat Commun* 10:5436. doi:10.1038/s41467-019-13225-y.
- Durbin R. 2014. Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT). *Bioinformatics* 30(9):1266–1272. doi:10.1093/bioinformatics/btu014.
- Excoffier L, Slatkin M. 1995. Maximum-likelihood estimation of molecular haplotype frequencies in a diploid population. *Mol Biol Evol* 12(5):921–927. doi:10.1093/oxfordjournals.molbev.a040269. (Abstract only.)
- Fragoso CA, Heffelfinger C, Zhao H, Dellaporta SL. 2016. Imputing genotypes in biallelic populations from low-coverage sequence data. *Genetics* 202(2):487–495. doi:10.1534/genetics.115.182071. (Abstract only.)
- Haldane JBS, Waddington CH. 1931. Inbreeding and linkage. *Genetics* 16(4):357–374.
- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. 2023. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. *Nat Genet* 55(7):1243–1249. doi:10.1038/s41588-023-01415-w.
- Howie BN, Donnelly P, Marchini J. 2009. A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. *PLoS Genet* 5(6):e1000529. doi:10.1371/journal.pgen.1000529.
- Kong A et al. 2008. Detection of sharing by descent, long-range phasing and haplotype imputation. *Nat Genet* 40:1068–1075. doi:10.1038/ng.216. (Abstract only.)
- Li H. 2011. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. *Bioinformatics* 27(21):2987–2993. doi:10.1093/bioinformatics/btr509.
- Li N, Stephens M. 2003. Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data. *Genetics* 165(4):2213–2233.
- Li Y, Willer CJ, Ding J, Scheet P, Abecasis GR. 2010. MaCH: using sequence and genotype data to estimate haplotypes and unobserved genotypes. *Genet Epidemiol* 34(8):816–834. doi:10.1002/gepi.20533.
- Li Z, Albrechtsen A, Davies RW. 2025. Flexible read-aware genotype imputation from sequence using biobank sized reference panels. *Nat Commun* 17(1):524 (online 13 December 2025). doi:10.1038/s41467-025-67218-1. PMC12804713. Preprint: bioRxiv 10.1101/2024.07.18.604149.
- Lin S, Cutler DJ, Zwick ME, Chakravarti A. 2002. Haplotype inference in random population samples. *Am J Hum Genet* 71(5):1129–1137. doi:10.1086/344347. (Full text not reachable.)
- Loh P-R, Danecek P, Palamara PF et al. 2016. Reference-based phasing using the Haplotype Reference Consortium panel. *Nat Genet* 48(11):1443–1448. doi:10.1038/ng.3679.
- Loh P-R, Palamara PF, Price AL. 2016. Fast and accurate long-range phasing in a UK Biobank cohort. *Nat Genet* 48(7):811–816. doi:10.1038/ng.3571.
- Marchini J, Howie B, Myers S, McVean G, Donnelly P. 2007. A new multipoint method for genome-wide association studies by imputation of genotypes. *Nat Genet* 39(7):906–913. doi:10.1038/ng2088. (Abstract only.)
- Marchini J, Howie B. 2010. Genotype imputation for genome-wide association studies. *Nat Rev Genet* 11(7):499–511. doi:10.1038/nrg2796. (Main text paywalled; Supplementary S3 opened.)
- O'Connell J et al. 2016. Haplotype estimation for biobank-scale data sets. *Nat Genet* 48(7):817–820. doi:10.1038/ng.3583. (SHAPEIT3: recursive k-means clustering to 4,000-haplotype clusters, not PBWT.)
- Pasaniuc B et al. 2012. Extremely low-coverage sequencing and imputation increases power for genome-wide association studies. *Nat Genet* 44(6):631–635. doi:10.1038/ng.2283.
- Rubinacci S, Delaneau O, Marchini J. 2020. Genotype imputation using the Positional Burrows Wheeler Transform. *PLoS Genet* 16(11):e1009049. doi:10.1371/journal.pgen.1009049.
- Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. 2021. Efficient phasing and imputation of low-coverage sequencing data using large reference panels. *Nat Genet* 53(1):120–126. doi:10.1038/s41588-020-00756-0. (Published text paywalled; method wording and numbers quoted from the bioRxiv version, 10.1101/2020.04.14.040329, 14 April 2020.)
- Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. 2023. Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes. *Nat Genet* 55(7):1088–1090. doi:10.1038/s41588-023-01438-3.
- Scheet P, Stephens M. 2006. A fast and flexible statistical model for large-scale population genotype data. *Am J Hum Genet* 78(4):629–644. doi:10.1086/502802.
- Stephens M, Smith NJ, Donnelly P. 2001. A new statistical method for haplotype reconstruction from population data. *Am J Hum Genet* 68(4):978–989. doi:10.1086/319501.
- Stephens M, Scheet P. 2005. Accounting for decay of linkage disequilibrium in haplotype inference and missing-data imputation. *Am J Hum Genet* 76(3):449–462. doi:10.1086/428594.
- Sun Q, Li Y. 2026. Advances in haplotype phasing and genotype imputation. *Nat Rev Genet* 27(2):155–169. doi:10.1038/s41576-025-00895-2. (Paywalled; not used for any claim above.)

Software documentation: Beagle 4.0 manual (3 March 2015) and Beagle 4.1 manual (21 January 2017),
faculty.washington.edu/browning/beagle/; Beagle 5 release notes, same site.

Code: all paths are under `tmp/phasing_research/repos/` (shallow clones of 2026-09-12). The
GLIMPSE2 facts marked "code-facts note" come from
`tmp/phasing_research/low_coverage_imputation/glimpse2_code_facts.md`; the lines cited directly
in this report were re-read before citing.
