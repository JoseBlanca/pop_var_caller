# Genotype inference in populations with few haplotypes

Segregant populations (F2, backcross, recombinant inbred lines, doubled haploids, MAGIC and NAM),
breeding populations with a pedigree, and inbred germplasm — from low-coverage sequencing.

Research note for the ng phasing/imputation module. Written 2026-09-12. Code citations point into
`tmp/phasing_research/repos/` (shallow clones of 2026-09-12); every code line quoted was opened and
checked. Paper citations were read from the open-access full text unless a note says otherwise.

## 0. Findings that shape the rest of the report

1. **Every tool in this family is one hidden Markov model.** The hidden state of an individual at a
   marker is *which founder haplotype(s) it carries there*; the state changes along the chromosome
   only at crossovers; the observation is the reads at the marker. The tools differ in three places
   only: how the transition rate is derived from the crossing scheme, how reads enter the emission,
   and how founder haplotypes are obtained when the parents were not sequenced. Sections 3 to 6 are
   organised around exactly those three differences.
2. **Reads pool across a founder block, so coverage per site is almost irrelevant.** With two
   founders that differ at one site in 2 kb, a 0.5× sample has about 250 informative reads per
   megabase. Two concordant reads separate the two homozygous founder states at 10,000:1; ten reads
   separate heterozygous from homozygous at about 1,000:1 (§3.5). That is why magicImpute reaches
   0.99 accuracy from 0.05× in an inbred biparental population and from 0.11× in an F2 (Zheng
   2018, Table 1 and Fig. 2; §9), and why 2,177 medaka F2 at 0.5× give r² = 0.98 (Pierotti 2024).
3. **Per-site genotype likelihoods are the sufficient statistic for this model.** The log-likelihood
   of a founder-pair state over a block is the sum, over its sites, of the sample's genotype
   likelihood for the genotype that founder pair implies. ng already computes those likelihoods
   (`ng_facts.md`). Read-backed phasing (chain ids) is not needed for the founder-mosaic model: phase
   is known once the founder pair is known.
4. **Marker-specific error is the one thing the simple binomial emission gets wrong**, and it is
   the difference between 70% and 95% correct calls at 20× in GBScleanR's F2 simulation (Furuta
   2023, Fig. 1). A per-marker reference-read bias and a per-marker mismapping rate, re-estimated
   from the cohort a few times, fix it. ng has the allele-level read counts to estimate both.
5. **Founders need not be sequenced.** Maximum-parsimony (Xie 2010), window clustering (FSFHap),
   joint Viterbi over founder patterns (GBScleanR), forward–backward maximisation over founder
   configurations (magicImpute), and EM over K ancestral haplotypes (STITCH) all recover the two
   founder haplotypes from the progeny alone. What is lost is only the *label* (which parent is
   which) and, at sites where the founders are identical, nothing can be recovered by anyone.
6. **Inbred lines collapse the model to a haploid one**: K states instead of K², and STITCH's
   `diploid-inbred` mode sets the heterozygote posterior to exactly zero
   (`STITCH/STITCH/R/functions.R:3633-3643`). RILs at F5–F8 still carry 1–6% residual
   heterozygosity, which a haploid model silently mis-calls; the right treatment is a diploid model
   whose selfing generations shrink the heterozygous state (GBScleanR, RABBIT), not a haploid one.
7. **The biparental HMM is Li & Stephens with two templates and a generation-scaled recombination
   rate**, literally: LB-Impute's switch probability `0.5·(1 − e^(−d/D))` is the Li & Stephens
   jump probability with K = 2 (§10). The identity breaks where the two homologous chromosomes stop
   being independent (selfing creates identity-by-descent that RABBIT and GBScleanR model as a
   coupled state) and where founders are unknown (the template set becomes a parameter).
8. **Cohort size, not coverage, is the binding constraint at the low end.** Xie 2010 needed ≥110
   RILs at 0.055× for 99% accuracy (Monte-Carlo, their SI); Pierotti 2024 saw a sharp fall in r²
   below 200 F2 of a single cross (their Fig. 3d). The module has to say what it does at 20
   progeny.
9. **A cohort can be recognised as few-haplotype from the called data** (§8): an allele-frequency
   spectrum concentrated at 0.5 (or 0.25/0.75 for a backcross), pairwise identity between
   individuals that is trimodal over megabase windows, and — the decisive test — a two-haplotype
   fit that explains the reads in a window as well as an eight-haplotype fit does.

## 1. Glossary

Genetics terms are brief; statistics and software terms are explained.

- **Founder / parent haplotype.** One of the chromosomes that entered the cross. A biparental
  population from inbred parents has two founder haplotypes per chromosome; an 8-way MAGIC has
  eight; a cross between two outbred parents has four.
- **Founder-origin state (RABBIT: "origin state"; GBScleanR: "haplotype pattern").** For a diploid
  individual at a marker, the ordered pair (founder of the maternal chromosome, founder of the
  paternal chromosome). With F founders there are F² ordered pairs, F(F+1)/2 unordered pairs.
- **Block, bin, ancestry block.** A run of consecutive markers where the founder-origin state does
  not change. Its ends are crossovers (in the individual) or, for a "bin", the set of positions
  where *any* individual in the cohort recombines (Huang 2009's "recombination bin").
- **Hidden Markov model (HMM).** A model with a hidden state that changes along the chromosome by
  a Markov rule (the next state depends only on the current one) and an observation at every
  marker whose distribution depends only on the state there. Three ingredients: the initial
  distribution, the transition probabilities between consecutive markers, the emission
  probabilities of the observation given the state.
- **Forward–backward.** The algorithm that gives, for every marker, the posterior probability of
  every state given *all* markers on the chromosome. Cost per individual: markers × states²
  (states × jump-targets when the transition has the special "jump anywhere" form).
- **Viterbi.** The algorithm that gives the single most probable state *path* through the whole
  chromosome. Same cost. It gives one answer and no probability per marker.
- **Emission with reads.** The probability of seeing r₁ reads of allele 1 and r₂ of allele 2 given
  the true genotype: a binomial with success probability 1 − ε for a homozygote and ½ for a
  heterozygote, ε being the per-read error rate. The binomial coefficient is the same for all
  genotypes and cancels.
- **Genotype likelihood (GL).** For one sample at one site, P(reads | genotype) for each possible
  genotype. The binomial above *is* a GL. ng computes GLs per sample per locus.
- **Map expansion.** The expected number of crossover junctions per Morgan accumulated on a
  chromosome after g generations. It is the quantity that scales the transition rate with the
  number of generations (§3.3).
- **Junction.** A point on a chromosome where the founder origin changes; the "junction density"
  is junctions per Morgan (Fisher's junction theory, as used by RABBIT and GBScleanR).
- **Identity by descent (IBD) between the two homologues.** After selfing, an individual's two
  chromosomes at a locus often descend from the same founder chromosome. The inbreeding
  coefficient f is the probability of that.
- **Allele read bias (GBScleanR: `bias`, symbol w_m).** At marker m, the probability that a read
  from a heterozygote shows the reference allele. Ideal value ½; departures come from
  amplification bias and mapping bias against the non-reference allele.
- **Mismapping rate (GBScleanR: `mismap`).** At marker m, the probability that a read from a
  homozygote shows the *other* allele for reasons other than sequencing error — reads from a
  paralogue mapped onto this site.
- **False homozygote.** A heterozygous individual whose few reads at a site all happen to show the
  same allele. With d reads it happens with probability 2·(½)^d = (½)^(d−1): 50% at d = 2, 6% at
  d = 5, 1.6% at d = 7.
- **Parent-independent genotyping.** Inferring the founder haplotypes from the progeny alone
  (Xie 2010's term).
- **Segregation distortion.** Allele frequencies in the progeny that depart from the Mendelian
  expectation (½ in an F2) because of selection.
- **Li & Stephens model.** The HMM in which a haplotype is a mosaic of template haplotypes, the
  hidden state is the template being copied, switches happen at a rate proportional to genetic
  distance and land on a uniformly random template, and the emission allows a mismatch with a
  small probability. The engine of every reference-panel imputer.

## 2. The populations and what "few haplotypes" buys

A natural population of N individuals carries up to 2N distinct haplotypes in any window, and an
imputer has to learn them. A segregant population carries a handful, and the crossing scheme fixes
their expected frequencies, the expected heterozygosity and the expected number of crossovers per
chromosome. That fixes the whole prior of the HMM before any data are seen.

| Population | Founder haplotypes per chromosome | Expected genotype frequencies (2 inbred founders A, B) | Heterozygosity | Crossovers per chromosome and generation of meiosis |
|---|---|---|---|---|
| F1 × inbred backcross (BC1) | 2 | AA : AB = 1 : 1 (or BB : AB) | ½ | one meiosis observed |
| F2 | 2 | AA : AB : BB = 1 : 2 : 1 | ½ | two meioses (both gametes of the F1) |
| RIL by selfing to generation F_g | 2 | AA : BB ≈ 1 : 1 | 1/2^(g−1): F6 3.1%, F8 0.8% | accumulates each generation: map expansion → 2× the F2 map at F∞ (Haldane–Waddington) |
| Doubled haploid (DH) from an F1 | 2 | AA : BB = 1 : 1 | 0 | one meiosis |
| 4- or 8-way MAGIC RIL | 4 or 8 | each founder ≈ 1/F | residual only | funnel crosses + selfing; RABBIT computes it |
| NAM | 2 per family, one shared | per family as RIL | residual only | as RIL |
| Breeding population with pedigree | few tens | pedigree-specific | mixed | per meiosis in the pedigree |
| Inbred germplasm (a panel of lines) | many, but each line is one haplotype | AA or BB only | ≈ 0 | not a cross; treated as haploids |

Sources for the expectations: standard Mendelian arithmetic; the F_g heterozygosity is what
GBScleanR's recursion produces by halving the non-IBD probability `a12` every selfing generation
(`GBScleanR/R/Methods-GbsrGenotypeData_HMM.R:846`); the RIL map expansion is what RABBIT's
`mapExpansion` accumulates (`RABBIT/RABBIT_Packages/MagicOrigin/MagicOrigin.m:266-270`).

Three consequences for inference:

- **Every site polymorphic between the founders is informative in every individual**, because
  every individual carries only founder alleles. There is no "rare variant" problem.
- **A block of founder origin is megabases long.** An F2 chromosome of 100 cM carries about two
  crossovers per homologue; a 100 Mb chromosome is then roughly four blocks per homologue. All
  reads in a block vote for the same state.
- **The transition prior is known**, so the model does not have to learn an effective population
  size or a recombination map from the data; it needs the crossing scheme and a genetic (or, as an
  approximation, physical) distance between markers.

## 3. The model when founders are known

### 3.1 In words

Take one progeny individual and walk along a chromosome. At each marker it carries two founder
haplotypes — say A on one homologue and B on the other. Between consecutive markers the pair
stays the same unless a crossover happened in one of the meioses between the founders and this
individual; a crossover swaps one member of the pair for another founder. The reads at the marker
are a noisy look at the genotype the pair implies (A/B → heterozygous; A/A → homozygous for the
allele founder A carries there). Forward–backward turns this into, for every marker, the posterior
probability of each founder pair; that posterior converts to a genotype posterior (and, since
founders are phased by definition, to phased haplotypes). Every tool in this section is this
model; the table lists where each one puts its three ingredients.

| Tool | Hidden states | Transition | Emission | Decoding |
|---|---|---|---|---|
| LB-Impute (Fragoso 2016) | 3: parent 1 hom, parent 2 hom, het (`LB-Impute/.../Impute.java:41,96`) | physical bp distance, one global scale (`FindPath2.java:95-111`) | binomial in read counts with error `-readerr`, floored/capped by `-genotypeerr` (`ImputeOffspring.java:288-317`) | windowed Viterbi, forward and reverse, conflicts → missing |
| TIGER (Rowan 2015) | 3: CC, LL, CL (`TIGER/hmm_prob.pl:124-126`) | empirical, per marker, trained per sample from a sliding-window pre-call (`hmm_prob.pl:152-169`) | empirical confusion matrix between window label and a 6-symbol base call (`hmm_prob.pl:400-434`) | Viterbi (jar), then breakpoint refinement |
| GBScleanR (Furuta 2023) | founder-pair patterns from a simulated crossing scheme (`HMM.R:582-646`) | continuous-time rate matrix from junction densities of the scheme, exponentiated over distance (`HMM.R:892-990`) | binomial with per-marker bias and mismap (`src/gbsrCalcProb.cpp:27-33, 49-154`) | Viterbi for founders jointly with progeny; forward–backward for progeny posteriors |
| magicImpute / RABBIT (Zheng 2015, 2018) | F² ordered founder pairs, or F under `depModel` (`RABBIT/.../MagicDefinition.m:80-84`, `MagicModel.m:123-151`) | continuous-time chain with rates from pedigree junction densities (`MagicModel.m:166-169`) | binomial with base error 10^(−Phred/10) and allelic error ε (`MagicLikelihood.m:305-320`) | forward–backward; founder configurations by forward–backward maximisation |
| STITCH (Davies 2016) | K² (diploid) or K (`diploid-inbred`) ancestral haplotypes (`STITCH/STITCH/R/functions.R:3744-3753`) | `exp(−nGen · d_Morgans)` (`R/genetic-map.R:28-42`) | per read, per base quality (`src/haploid.cpp:174-190`) | forward–backward + EM over haplotypes |
| AlphaPlantImpute2 | ordered pairs of library haplotypes (wheel `CombinedHMM.py:314-324`) | uniform `recomb/n_loci` per marker (`alphaplantimpute2.py:606-607`) | hard genotype with flat error `-error` (default 0.01) | forward–backward ("marginalize") or sampling |
| FSFHap / FILLIN (Swarts 2014) | 5 (parent A, three intermediate het states, parent B) | fixed constants ≈ 0.999 stay (TASSEL source, §4.3) | fixed constants ≈ 0.998 match | Viterbi |

### 3.2 The state space

For F founders and a diploid individual, RABBIT enumerates `nFgl²` ordered pairs
(`RABBIT/RABBIT_Packages/MagicDefinition/MagicDefinition.m:80-84`) and, when the individual is
assumed fully inbred (`depModel`), only F states (`MagicReconstruct/MagicModel.m:143-149`). When
founders are outbred, each founder contributes two haplotypes and `nFGL = 2·nFounder`
(`MagicModel.m:203-206`), so a cross of two outbred parents has 16 ordered-pair states. GBScleanR
does not enumerate; it simulates the crossing scheme generation by generation (`.mateGametes`,
`HMM.R:557-580`) and keeps the set of haplotype-label pairs that the scheme can produce, plus the
set of founder *genotype* patterns: 2^F − 2 patterns for inbred founders, 2^(2F) − 2 for outbred
ones (`HMM.R:510-537`). LB-Impute hard-codes three states and its authors note that n founder
haplotypes would need n(n−1)/2 + n states (Fragoso 2016, eq. 5).

For the biparental, inbred-founder case that ng cares about most: three states (AA, AB, BB), or
four ordered ones if phase within the heterozygote is tracked.

### 3.3 Transitions: how the crossing scheme and the generations enter

**The simplest form (LB-Impute).** The probability of staying in a state across a gap of d bp is
`0.5·(1 + e^(−d/D))` and of leaving it `0.5·(1 − e^(−d/D))`, with D = `-recombdist`
(`LB-Impute/SourceCode/LB-Impute/src/imputation/FindPath2.java:95,103`; paper eqs. 3–4). No
genetic map and no generation count: D is the physical distance at which the two probabilities
meet, 10 Mb in the paper ("expected distance of 50 cM"), 1 Mb in the code default
(`Impute.java:32`). A homozygous-to-other-homozygous jump is charged as two crossovers — the
switch probability squared (`FindPath2.java:99`) — because in an F2 that needs a crossover in
both gametes. The only concession to population type is the `-dr` flag that removes the squaring
for RILs, and the code applies it only to the first window of each chromosome (`FindPath2.java:153,157`
against the unconditional square at `:99`).

**Generation-scaled form (STITCH).** The probability of *no* jump across an interval of d Morgans
is `exp(−nGen · d)` (`STITCH/STITCH/R/genetic-map.R:28-42`); without a genetic map d is
`expRate/100/1e6 · bp` with `expRate` in cM/Mb. For a diploid the two chromosomes are independent,
so the three transition classes (no jump, one jump, two jumps) are `σ², σ(1−σ), (1−σ)²`
(`R/functions.R:4026-4043`). `nGen` is the number of meioses since founding; for an F2 it is 2
(Pierotti 2024 used `nGen = 2`), and STITCH's help says "the algorithm is relatively robust to
this" (`STITCH/STITCH.R:29-31`).

**Junction-density form (RABBIT, GBScleanR).** Both compute, from the crossing scheme, the
expected density per Morgan of each *kind* of founder-origin change on the pair of homologues:
a change on one homologue only, a simultaneous change on both (which happens when the two
homologues are IBD and the junction was inherited from a common ancestor), and so on. RABBIT
packs them into a six-vector `{f, j1122, j1211, j1213, j1222, j1232}` — f the inbreeding
coefficient, the rest junction densities (`RABBIT/RABBIT_Packages/MagicReconstruct/MagicModel.m:189`)
— and builds the maternal and paternal map-expansion rates from them,
`Rm = 2 j1222 + j1122 + j1232`, `Rp = 2 j1211 + j1122 + j1213` (`MagicModel.m:87-88`). Under
`jointModel` the two homologues share one continuous-time chain whose stationary distribution puts
`f/F` on each IBD pair and `(1−f)/(F(F−1))` on each non-IBD pair; under `indepModel` the chain is
the Kronecker sum of two independent single-homologue chains with rates Rm and Rp; under
`depModel` only one homologue is tracked with rate `(Rm+Rp)/2` (`MagicModel.m:123-151`). The
transition matrix over a marker interval of d Morgans is the matrix exponential
`MatrixExp[transitionRate · d]` (`MagicModel.m:166-169`).

How the number of generations enters: the crossing scheme is a list with one entry per
generation, `{"Pairing","Pairing","Selfing","Selfing"}` for a 4-way RIL with two selfing
generations (RABBIT manual, p. 5). `origSummary` sets `nGeneration = Length[mateScheme] + 1` and
recurses generation by generation over population size, coalescence probability, inbreeding and
map expansion (`RABBIT/RABBIT_Packages/MagicOrigin/MagicOrigin.m:295-313`); map expansion is the
running sum of the non-IBD probability, `founderR + Prepend[Most[Accumulate[alpha12]], 0]`
(`MagicOrigin.m:266-270`), so R grows by about one crossover-equivalent per generation until the
homologues become IBD and further selfing adds nothing. GBScleanR's version is the same recursion
in R: `.calcNextJnum` halves the non-IBD probability every selfing generation
(`next_a12 <- 0.5 * prob_df$a12`) and accumulates junctions (`next_r <- prob_df$r + prob_df$a12`)
(`GBScleanR/R/Methods-GbsrGenotypeData_HMM.R:846-855`); `.getXoFreq` turns the counts into the
four rates the rate matrix needs (`HMM.R:884-890`); the matrix is exponentiated per interval with
`rf = distance_bp · 1e-6 · recomb_rate`, `recomb_rate` default 0.04 crossovers per Mb
(`HMM.R:976-979`; `man/estGeno.Rd:11`). The initial distribution is the stationary vector of the
first interval's matrix (`HMM.R:992-998`), which for an F2 is 1 : 2 : 1.

The practical difference between the three forms: LB-Impute's constant needs hand tuning per
population and cannot represent that an F8 RIL has twice the F2 map and one-hundredth of its
heterozygosity; STITCH's `nGen` gets the map expansion right but not the heterozygosity (its
homologues stay independent); the junction form gets both, at the price of asking for the
crossing scheme.

### 3.4 Emissions: how reads enter

Three tools use the same binomial-without-coefficient likelihood on the two allele counts:

- LB-Impute: `check1 = (1−err)^nonref · err^ref`, `check2 = (1−err)^ref · err^nonref`,
  `check3 = 0.5^ref · 0.5^nonref` (`LB-Impute/.../ImputeOffspring.java:288-290`), `err` =
  `-readerr`, default 0.05 (`Impute.java:31`). Then a rescaling: every state's value is divided
  by the best state's and mapped onto `[genotypeerr, 1 − genotypeerr]`
  (`ImputeOffspring.java:300`; `maxprob = 1 − 2·genotypeerror`, `minprob = genotypeerror`,
  `:60-61`). So `-genotypeerr` (default 0.05) is a floor on how far one marker can push the path,
  not an error model — it caps the per-marker likelihood ratio at 19:1. With n reads of one
  allele the heterozygous state keeps `0.05 + 0.90·(0.5/0.95)^n`: 0.52 at n = 1, 0.11 at n = 5,
  0.057 at n = 10.
- magicImpute: `p11 = (1−ε_b)^n1 · ε_b^n2`, `p12 = p21 = (½)^(n1+n2)`, `p22` symmetric, with
  `ε_b = 10^(−minPhredQualScore/10)` (default Phred 30 → 0.001), then convolved with a
  depth-independent allelic error ε (default 0.005) that models a whole-genotype mistake such as
  a misaligned marker (`RABBIT/RABBIT_Packages/MagicReconstruct/MagicLikelihood.m:305-320`;
  manual p. 8; Zheng 2018, Methods).
- GBScleanR: the same binomial but with the heterozygote's reference-read probability equal to
  the per-marker bias `w1` instead of ½, and every probability squeezed into `[ε_seq, 1 − ε_seq]`
  (`GBScleanR/src/gbsrCalcProb.cpp:27-33`, `ε_seq` default 0.0025, `man/estGeno.Rd:12`); on top
  of that a per-marker mismapping mixture that moves probability from each homozygote to the
  heterozygote (`gbsrCalcProb.cpp:49-154`). Section 5 has the estimation.

STITCH works per read rather than per allele count: for each read, and each ancestral haplotype
k, the emission is `e_k · pA + (1 − e_k) · pR`, where `e_k` is haplotype k's frequency of the
alternate allele at that site and `pA`, `pR` come from the base quality (`ε = 10^(−BQ/10)`,
mismatch spread as ε/3) (`STITCH/STITCH/src/haploid.cpp:174-190`). A read spanning several sites
multiplies its per-site terms (the loop over `j`), which is the only place in this family where
read-backed phase information enters — and STITCH caps it at `Jmax` sites per read.

TIGER does something different and worth knowing: it turns the counts at each marker into one of
six symbols — AA, BB, AB, AU, BU, UU, where U means "uncertain" and the homozygous symbols need
five reads of one allele (Rowan 2015, Methods) — then *trains* the emission matrix per sample as
the empirical confusion between a coarse sliding-window genotype and the symbol
(`TIGER/hmm_prob.pl:400-434`). In the shipped model 97.9% of markers emit UU and the informative
signal is the 20 : 1 asymmetry between CU and LU under the two homozygous states
(`TIGER/sample_hmm_model:19-21`).

AlphaPlantImpute2 and AlphaImpute2 take hard genotypes only (0/1/2/missing) with a flat error:
`point_estimates[i,j,k] = 1 − error` if the library pair's genotype equals the observed one, else
`error` (wheel `alphaplantimpute2/tinyhouse/CombinedHMM.py:314-324`); no read counts or genotype
likelihoods are read anywhere (`alphaplantimpute2.py:430` requests `reads=False`). They are array
tools that happen to be used on sequence-derived calls.

**The ng-relevant observation.** In every read-count emission above, the per-site quantity is a
genotype likelihood: P(reads | AA), P(reads | AB), P(reads | BB). The founder-pair state only
selects which of the three to use at each site. A caller that keeps the three GLs per sample per
site has everything the emission needs; it does not need the reads.

### 3.5 How few reads are needed: the arithmetic

Fix a 1 Mb segment inside one founder block of one individual, with 500 sites where the two
founders differ (one per 2 kb) and mean depth 0.5×. Each site is covered by Poisson(0.5) reads,
so the segment holds about 500 × 0.5 = 250 informative reads. Take a per-read error ε = 0.01
(magicImpute uses 0.001; LB-Impute 0.05).

*Homozygous A against homozygous B.* A read showing A has likelihood 1 − ε under AA and ε under
BB: a ratio of 99. Two concordant reads: about 10⁴. Under BB the 250 reads would contain about
2.5 A reads; under AA about 247. The two homozygous states are separated by the first handful of
reads; the remaining 240 are redundant.

*Heterozygous against homozygous A.* Under AB each read shows A with probability ½ (or w_m); under
AA with probability 1 − ε. n reads that all show A have ratio ((1 − ε)/½)^n = 1.98^n in favour of
AA: 10 : 1 after 3 reads, about 950 : 1 after 10, 10⁶ after 20. Conversely a single B read among
A reads has ratio (½)/ε = 50 in favour of AB. So a heterozygous block declares itself as soon as
both alleles are seen, and a homozygous block needs about ten reads to rule out heterozygosity at
1,000 : 1. At 0.5× and one site per 2 kb, ten reads is twenty sites, 40 kb. At 0.1× it is 200 kb.

*Prior.* In an F2 the prior is 1 : 2 : 1, i.e. even odds on het; these ratios are the posterior
odds up to that factor of two. In an F8 RIL the prior on het is 1/128, which is why RIL callers
can run haploid — and why a real residual heterozygote in a RIL needs about seven more reads of
evidence than in an F2 to overcome its prior.

*Breakpoint resolution.* A crossover between two homozygous blocks lies between the last read
showing A and the first read showing B. Informative reads occur every 2 kb / 0.5 = 4 kb on
average at 0.5×, and every 20 kb at 0.1×; the breakpoint interval has that expected width. TIGER
measured this: 261,795 markers over the 119 Mb Arabidopsis genome is one per 455 bp; at 0.1× that
is one informative read per 4.5 kb, and their median crossover resolution at 0.1× was 1,986 bp
(median of an exponential is 0.69 × its mean) — Rowan 2015, Table 2. A crossover between a
heterozygous and a homozygous block is coarser by the factor above: the homozygous side needs
about ten reads before the model believes it, so about 40 kb at 0.5× and 200 kb at 0.1×.

*Coverage floors reported by the tools match this.* magicImpute's accuracy-versus-depth curves
break at 0.053× for a two-inbred-founder RIL and at 0.11× for an F2 from the same founders (Zheng
2018, Table 1): the F2 needs about twice the reads because it has to resolve heterozygosity. LB-
Impute holds > 99% accuracy down to 0.1× on simulated F2 and BC1 with 10,000 markers per 100 Mb
(Fragoso 2016, Fig. 2). Huang 2009 called 150 rice RILs at 0.02× with 1 SNP per 40 kb and a
15-SNP window. Xie 2010 called 238 rice RILs at 0.055× *without parents*.

*What does not scale.* Sites where the founders are identical carry no information at any depth;
a block that contains no founder-differing site is invisible. And a cohort of 20 progeny at 0.1×
has, over the whole cohort, 20 × 0.1 = 2× coverage of each founder haplotype — enough to call the
founders' genotypes only where the two haplotypes can be told apart, which is the founder-unknown
problem of §4.

### 3.6 Decoding and what each tool reports

- LB-Impute runs a *windowed* Viterbi: it enumerates every path through the next `-window` markers
  (default 7), commits the first step of the best path, and slides (`FindPath2.java:85-134`;
  paper Fig. 1). It does this forward and reverse and sets a marker to missing where the two
  disagree (`ImputeOffspring.java:258-269`), unless `-resolveconflicts` picks the higher-scoring
  direction. Markers with no reads are dropped from the chain and filled afterwards from agreeing
  neighbours (`CollapseValuesOffspring.java:36-51`, `ImputeOffspring.java:200-237`). There is no
  posterior; accuracy fell from > 99% at window 7 to 67–71% at window 2 (Fragoso 2016, Fig. S4),
  and run time rose to 17,879 s at 2.5× (Fig. S5).
- GBScleanR runs Viterbi forward and reverse for the founders jointly with the progeny, splices
  the halves at the chromosome midpoint (`HMM.R:1230-1250`), then forward–backward for the progeny
  conditioned on the founder path (`src/gbsrFB.cpp:237-309`), and outputs genotype posteriors,
  best haplotype pairs and, for a two-inbred-founder cross, dosages (`HMM.R:115-120`). Calls with
  posterior below `call_threshold` (0.9) are set missing (`HMM.R:1332-1352`).
- magicImpute outputs posterior genotype probabilities and phased founder-origin posteriors
  (manual p. 9: `p11|p12|p21|p22`).
- STITCH outputs dosages and genotype probabilities per site.

## 4. When founders are not sequenced

The founder haplotypes are two (or F) unknown binary vectors over the polymorphic sites. Every
method below estimates them from the progeny by exploiting the same fact: **within a block, every
progeny is a copy of one founder (homozygous) or of both (heterozygous)**, so across the cohort
the reads at neighbouring sites fall into two clusters, and the clustering is consistent along the
chromosome except at crossovers.

### 4.1 Maximum parsimony of recombination (MPR; Xie 2010)

Rice RILs (Zhenshan 97 × Minghui 63), 238 lines at 0.055× each, parents not sequenced beyond a
low-coverage draft used only to name them (Xie 2010, PNAS 107:10578, doi 10.1073/pnas.1005931107).
Steps, from the paper: SNPs are found from the pooled RIL reads as sites with two alleles; the
chromosome is cut into windows of about 50 SNPs; within a window two arbitrary parental columns
are seeded and alleles are swapped between the columns one at a time (or in larger steps) while
the total number of recombination events implied across all RILs decreases, until no swap
reduces it; the configuration with the fewest recombinations is the founder pair. The result is
then made robust by resampling 50 of the first 100 SNPs in each RIL 20 times and computing a
three-state posterior (parent 1, parent 2, "inferior SNP") per site, which removed 94.1% of the
bad SNPs while keeping 82–84% of the SNPs (their Results). Accuracy: 11,792 of 13,227 retained
chromosome-5 SNPs were confirmed by later deep sequencing of the parents and *all* inferred
parental calls at those sites matched. Assumptions: two founders, fully homozygous lines (0.34%
heterozygous sites were masked), biallelic sites, random sparse errors, and a cohort large enough
— Monte-Carlo gave ≥ 110 RILs for 99% accuracy at that coverage. Progeny were then genotyped with
a three-state Viterbi HMM (1 cM per 244 kb, 3.18% error) and lumped into 143 bins on chromosome
5 with 208.5 kb mean length. The CRAN package `MPR.genotyping` implements it.

### 4.2 Window clustering (FSFHap; Swarts 2014)

FSFHap "works by first identifying parental haplotypes in the progeny" (Fragoso 2016,
Introduction, describing Swarts 2014). In the TASSEL 5 source (GitHub mirror
`zacharymiller90/tassel-ML`, `NucleotideImputationUtils.java`; read through the mirror on
2026-09-12, line numbers as reported there), progeny haplotypes in a window of `window` sites
(default 50, overlap 25) are clustered with a maximum-difference threshold `maxDiff` (default 0),
a haplotype must be seen `minHap` times (default 5) to count, and the two largest clusters become
the parents; `windowLD` first drops sites whose correlation with their neighbours is below
`minR` (0.2); `phet` (default 0.07) is the expected fraction of heterozygous sites used to reject
sites; `bc` selects a backcross variant that checks segregation at 1 : 1 and takes the
major/minor allele directly (`callParentAllelesByWindowForBackcrosses`). The progeny are then
decoded by a 5-state Viterbi with fixed transition constants (0.999 on the diagonal, 10⁻⁴–5·10⁻⁴
off it) and fixed emissions (0.998 for a matching homozygote; the three middle "heterozygous"
states emit 0.6/0.2/0.2, 0.4/0.2/0.4, 0.2/0.2/0.6) and no use of depth — the same constants are
in `FILLINImputationPlugin.java:111-124`. What "minhap 2" in LB-Impute's comparison meant is
lowering that count for low coverage (Fragoso 2016, SI Note A). I could not read Swarts 2014
itself (the publisher, a mirror and the TASSEL wiki all refused the fetch); the paper's accuracy
figures are therefore not quoted here, only what LB-Impute, NOISYmputer and magicImpute measured
against FSFHap (§9).

### 4.3 Joint Viterbi over founder patterns (GBScleanR)

GBScleanR always treats founder genotypes as hidden: at each marker it scores every founder
genotype pattern by the founders' own reads (`calcPemit`, `src/gbsrCalcProb.cpp:197-262`) plus
the *sum over all progeny* of their likelihood under that pattern
(`src/gbsrViterbi.cpp:322-348`), and runs Viterbi over patterns along the chromosome. Founders
with no reads contribute a flat term and are inferred from the progeny alone. The explicit
parentless mode fabricates two dummy founders with `dummy_reads` (default 5) reads of opposite
alleles at every marker (`HMM.R:438-445`); the manual says this "assumes that the given population
is a bi-parental population", that `dummy_reads = 0` should be used for outbred parents, and that
it "is less accurate and has more chance to get a genotype estimate randomly selected from the
equally likely genotype estimates" (`man/estGeno.Rd:79-96`). That randomness is the label
ambiguity: at a marker where the dummy alleles are the wrong way round relative to the previous
marker, the two founder patterns score equally except for the transition term.

### 4.4 Forward–backward maximisation over founder configurations (magicImpute)

Zheng 2018 imputes founders jointly with all offspring: at each locus the candidate founder
haplotype configurations are enumerated (up to `2^(nFounder·(1+[outbred]))`,
`RABBIT/RABBIT_Packages/MagicImpute/MagicImpute.m:175-199`), a forward pass over markers computes
their posterior using every offspring's HMM under an approximation of offspring independence
given the founders, and a backward pass picks the maximum-likelihood configuration
(`parentForwardCalculation`, `parentBackwardMaximize`, `MagicImpute.m:191-196`; Zheng 2018,
Algorithms A and B). A second pass in the reverse direction fixes the chromosome ends. The paper
states the algorithm "is parent independent so that it applies even if some founders' genotypes
are not available", and its founder-imputation step is its run-time bottleneck (Zheng 2018,
Discussion; Table 2: 784 s for 275 maize RILs × 13,912 SNPs against 178 s for Beagle).

### 4.5 EM over K ancestral haplotypes (STITCH)

STITCH never had founders: it alternates forward–backward over the K² (or K) states with an EM
update of the K ancestral haplotypes' allele frequencies (`eHapsCurrent`), 40 iterations by
default (Davies 2016, Online Methods; `STITCH/STITCH.R:183-186`). A biparental population is the
case K = 2 (Pierotti 2024 used K = 16 for ten F2 crosses sharing eight inbred founders,
`nGen = 2`). Identifiability: the labels of the K haplotypes are arbitrary, and with K > the true
number the extra haplotypes soak up noise; the README's advice is to keep K small enough that
"each ancestral haplotype gets at least … 10×" of pooled coverage (`STITCH/README.md:172-174`).

### 4.6 The others

- **TIGER requires known parents.** Its input is `chr pos parent1allele count parent2allele count`
  (`TIGER/README.md:24`), markers come from a user-supplied Col × Ler list, and there is no
  parent-inference step anywhere in the repository (`bamtotiger.sh:10-20,53-56`).
- **Huang 2009 sequenced both parents** and only used the RILs to call; **mpimpute (Huang 2014)**
  imputes *missing* founder genotypes in a MAGIC by auditing the progeny that inherited from that
  founder with high probability (Huang 2014, Genetics 197:401), which needs the other founders to
  be typed.
- **AlphaPlantImpute (Gonen 2018, bioRxiv 10.1101/330027)** phases the parents from the progeny
  heuristically: it lists markers where parents differ, infers linked allele pairs from
  homozygous offspring counts, then tracks parent-of-origin through each offspring — hard
  genotypes only, 22 s per run, accuracy 0.96 for inbred parents at 25k markers and 50
  low-density markers per chromosome (their Results). Its successor AlphaPlantImpute2 instead
  builds a library from high-density individuals by random initial phasing plus repeated HMM
  resampling (`alphaplantimpute2.py:96-131, 190-228`), i.e. the STITCH idea with a fixed set of
  library haplotypes.

### 4.7 What they all assume, and where identifiability stops

- **Biallelic sites, allele frequency near ½ in the cohort** (¼ / ¾ in a backcross). A site whose
  minor allele is rare in the cohort is either a founder-monomorphic site with errors or a new
  mutation; MPR discards it, GBScleanR's `setCallFilter` and NOISYmputer's pre-filters remove it,
  STITCH's EM assigns it to noise.
- **Two haplotypes explain each window.** Within a block the cohort's read vectors must fall into
  two clusters. Where founders are identical over a stretch there is one cluster and nothing to
  infer; the state simply does not change there, which is harmless for genotyping (the genotype
  is the same under any state) and fatal for breakpoint placement.
- **Which parent is which is not identifiable** from the progeny; every method attaches an
  arbitrary label per chromosome (MPR uses the parental draft to orient; GBScleanR draws at random
  among ties, `src/gbsrutil.cpp:146-179`). Two chromosomes' labels are independent. A single
  low-coverage parent, or even a known allele at one site per chromosome, fixes it.
- **Residual heterozygosity in RILs** cannot be distinguished from a genotyping error at low
  coverage unless the model allows a heterozygous *block*: a run of both alleles over many sites.
  Methods that pre-call genotypes under a homozygosity assumption lose it; Zheng 2018 lists this
  as a limitation of magicImpute's own pipeline ("prior transformation of allelic depths to called
  genotypes with homozygosity assumption loses residual heterozygote information"), and the
  medaka F2 study saw its lowest r² (0.694) on the rarest SNPs of a single cross (Pierotti 2024).
- **Cohort size.** Xie 2010: ≥ 110 RILs at 0.055× for 99%; Pierotti 2024: r² falls sharply below
  200 F2 of one cross (Fig. 3d). Nobody reports the founder-unknown case at 20 progeny; §11 says
  what ng should do there.

## 5. Error handling

### 5.1 The error sources and who models which

| Error | Nature | LB-Impute | TIGER | GBScleanR | magicImpute | STITCH | NOISYmputer |
|---|---|---|---|---|---|---|---|
| sequencing error per read | uniform, ~0.1–1% | `-readerr` 0.05 | 1% in base caller | `error_rate` 0.0025 | 10^(−Phred/10) | per-read base quality | two global rates e_a, e_b |
| allele read bias per marker | w_m ≠ ½ from amplification/mapping | no | absorbed by per-sample training | **`bias` per marker, estimated** | no | learned into ancestral haplotype frequencies | no (reviewers raised it) |
| mismapping per marker | paralogue reads onto the site | `-genotypeerr` cap 0.05 | no | **`mismap` per marker, estimated** | allelic error ε 0.005 | no | pre-filter "alien segments", "incoherent SNPs" |
| missing marker (no reads) | | dropped, refilled | UU symbol | flat emission | flat | flat | window of neighbours |
| false homozygote | het seen as hom at low depth | emission keeps het at 0.52 for one read | AU/BU symbols below 5 reads | binomial | binomial | per read | window likelihood |
| segregation distortion | | — | — | not modelled; measured after | assumed absent (stated) | — | — |

### 5.2 Naive window smoothing against a model

Huang 2009 called each RIL with a 15-SNP sliding window and thresholds on the ratio of parental
alleles (≥ 11 : 4 → homozygous indica, 10 : 5 to 3 : 12 → heterozygous, ≤ 2 : 13 → homozygous
japonica), and placed a breakpoint where the ratio crossed 8 : 7 exactly once (Huang 2009,
Methods). Its "expected accuracy of 99.94%" is a Bayesian computation from the measured parental
error rates (4.12% and 0.71%) and the F11 genotype proportions 49.98 : 0.05 : 49.98, not a
measurement against truth; and at 0.02× with one SNP per 40 kb a 15-SNP window is 600 kb, which
is the resolution.

What the window loses relative to the HMM: (a) the window is the same width everywhere, so it is
too wide in read-rich regions and too narrow in read-poor ones, whereas the HMM's evidence per
state is whatever the reads supply; (b) the thresholds ignore depth — 11 : 4 from 15 reads and
from 150 reads are treated alike; (c) a breakpoint inside the window is invisible until the ratio
crosses the threshold, so its position is biased toward the window centre; (d) two windows'
decisions are not forced to be consistent. TIGER is the intermediate design: it uses a 1,000-marker
window only to *label* markers for training and then lets the HMM decide (`TIGER/hmm_prob.pl:152-169`),
with a beta-mixture fit whose only output is the two cut-offs on the window statistic
(`TIGER/beta_mixture_model.R:201-207`, default −50/+50 with clamps at `:210-225`). NOISYmputer
(Triay 2025) is a window method with a likelihood inside the window (binomial over the m SNPs on
each side, m = 30 for rice) and a Kosambi-based gap-filling step; it reaches 99.9% breakpoint
precision on simulated F2 at 0.5–3× and error ≤ 2%, but degrades above 5% error and only tests
for a single transition per support interval (Triay 2025, Results and Limitations).

### 5.3 What GBScleanR's per-marker parameters buy

GBScleanR's bias `w_m` is estimated from the individuals currently called heterozygous at m as
(reference reads per reference allele copy) / (that + alternate reads per alternate allele copy)
— `bias <- ref_prop / (ref_prop + alt_prop)` (`HMM.R:1362-1372`) — and pooled with a second
estimator from homozygotes when the two correlate above 0.7 (`HMM.R:1402-1419`); the mismapping
rates are the fraction of confidently homozygous individuals that nevertheless show the other
allele (`HMM.R:1421-1449`). Four cycles of estimate-genotypes / re-estimate-parameters are the
default (`iter = 4`; vignette: "usually saturates the improvement"). Both are clamped to
`[error_rate, 1 − error_rate]` (`HMM.R:1579-1583`).

The measured effect (Furuta 2023, Genetics 224:iyad055, Fig. 1 and Supplementary Data 1;
simulation of 620 markers on a 50 Mb, 2 Morgan chromosome, founders at 5×, offspring at 0.1×
to 20×): in a 1,000-individual F2 with error-prone markers, GBScleanR exceeds 95% correct calls
at 20× where LB-Impute and magicImpute plateau near 70% — and the competitors' miscall rate
*rises* with depth, because more reads at a biased marker make the wrong genotype more
confident. In the 8-way RIL simulation the tools did not differ, because homozygous populations
give bias no heterozygotes to corrupt. On the real rice F2 (814 individuals, 5,035 SNPs, 0.85×,
*O. sativa* × *O. longistaminata*) GBScleanR was best on 11 of 12 chromosomes, by 10.5 points
over magicImpute on the most biased chromosome (Fig. 3a), and called 14 double crossovers inside
1 Mb against magicImpute's 235 and LB-Impute's 948 (Table 2) — spurious short blocks are the
signature of unmodelled marker error.

### 5.4 Missing data

Zero reads at a marker is not an error; every likelihood-based tool gives it a flat emission and
lets the neighbours decide. LB-Impute is the exception: it removes the marker from the chain and
fills it afterwards only if both flanking states agree (`ImputeOffspring.java:200-237`), which is
why it "left more missing markers (minimum 2.35, maximum 150.57 per event)" around breakpoints
than FSFHap (Fragoso 2016, Results). The consequence for the module: a marker with no reads in
a sample must still be a position in the chain, with emission 1 for every state.

## 6. Multiparental and pedigree methods

### 6.1 RABBIT / magicImpute (Zheng 2015, 2018)

RABBIT's model is §3.3's junction-density chain with F² states. It was built for SNP arrays on
MAGIC (19 Arabidopsis founders), Collaborative Cross and Diversity Outbred mice (Zheng 2015,
Genetics 200:1073, doi 10.1534/genetics.115.177873); its Fig. references there compare
`jointModel` against the two extremes and report a wrongly-called probability of 0.008–0.105
against HAPPY's 0.027–0.225 on simulated data, 93.5% concordance with HAPPY on 703 MAGIC lines,
and about 360 s for 103 pre-CC lines × 14,076 markers in Mathematica. Zheng 2018 added the read
model and founder imputation and measured: break-point depth (the depth below which accuracy
falls) of 0.053× for a two-founder RIL, 0.11× for an F2, 0.21× for an 8-founder MAGIC and 0.21×
for a cross between two outbred parents, with maximum accuracies 0.985–0.995 (Zheng 2018, Table 1,
Fig. 2, 200 simulated offspring each). Scaling: the state space is F² and the founder
configuration space 2^F per locus, so the founder step is exponential in F for outbred founders;
Table 2 gives 3,170 s for 178 rice MAGIC lines × 37,240 SNPs.

### 6.2 AlphaImpute2 (Whalen & Hickey 2020)

Two components, run as a pipeline rather than alternated (`AlphaImpute2/src/alphaimpute2/alphaimpute2.py:334-412`):

- **Multi-locus iterative peeling** over the pedigree: for each individual three probability
  vectors over the four phased genotypes — anterior (from parents), penetrance (own data),
  posterior (from offspring) — and a segregation forward–backward that tracks which parental
  haplotype was inherited at each locus (`AlphaImpute2/peeling_notes.md:11-42, 108-112`).
  Default 4 cycles with the first at a 0.99 calling cut-off (`Heuristic_Peeling.py:52,76`).
- **Population imputation**: not a forward–backward over explicit haplotypes but a particle
  sampler over a positional Burrows–Wheeler transform of the library, where a state is a pair of
  PBWT intervals rather than a pair of haplotypes (`PhasingObjects.py:339,403`); 40 particles per
  phasing cycle, 100 per imputation, 5 phasing cycles (`alphaimpute2.py:152-172`). The error rate
  is hard-coded at 0.01 with the comment `# FLAG: error_rate is hard coded.`
  (`ParticlePhasing.py:175-187`); the transition rate is `map_length / nLoci` per marker.

Order: initial peeling, phase the high-density individuals to build the library, impute
low-density "pseudo-founders" from the library, mask their parents and peel again
(`alphaimpute2.py:365-408`). It takes hard genotypes only — no `-seqfile`, no likelihoods
(`alphaimpute2.py:34-49,455`). Whalen & Hickey report 0.993 accuracy for low-density (10k marker)
individuals in a 107,000-animal pedigree against Beagle 5.1's 0.942, in 105 min for 5,000 markers
of one chromosome (bioRxiv 10.1101/2020.09.16.299677, Table 1). Its relevance to ng is the
architecture: pedigree peeling for individuals whose parents are typed, a haplotype library for
the rest.

### 6.3 AlphaPlantImpute2

A library of phased haplotypes from high-density individuals (two per individual, one under
`-haploid`; `alphaplantimpute2.py:138-143`), refined by ten rounds of HMM resampling in which
each individual is re-phased against a random subsample of `-n_haplotypes` (100) others
(`alphaplantimpute2.py:190-228`). Imputation is a diploid HMM over ordered pairs from a *targeted*
subsample: the chromosome is cut into `-n_windows` (5) bins and in each bin the haplotypes with
fewest opposite-homozygous mismatches to the target are kept (wheel `HaplotypeLibrary.py:153-176`);
with `-founders` the state space is exactly the named founders' haplotypes
(`alphaplantimpute2.py:255`). Transition: `recomb / n_loci` per marker, i.e. one Morgan per
chromosome spread uniformly (`alphaplantimpute2.py:606-607`). No paper describes it; the one
published evaluation (Niehoff 2022, Plant Genome, doi 10.1002/tpg2.20257; read from the bioRxiv
full text 10.1101/2022.03.29.486246) ran it on 1,357 sugar-beet lines from 36 families with
19,519 high-density and 950–2,719 low-density array markers: accuracy 0.90 against a tuned Beagle
5.1's 0.90 (untuned Beagle 0.22), perfect imputation of F1 from DH parents where Beagle gave
0.92, constant accuracy from 5 to 50 offspring per family where Beagle's fell, and 87 min against
Beagle's 187 min. Beagle needed `ne = 4` and thousands of iterations to work on a breeding
population at all (Niehoff 2022, Results).

### 6.4 FSFHap and FILLIN

FILLIN's library is built from inbred lines in blocks of 64 sites; a target is matched to the
donor with mismatch rate below `maximumInbredError` (0.01) or, failing that, to a pair of donors
by a 5-state Viterbi with mismatch rate below `maxHybridErrorRate` (0.003), from at most 20
donor hypotheses (`FILLINImputationPlugin.java:64-80,111-124`, via the GitHub mirror). It scales
linearly in targets and was designed for maize inbreds and NAM (Swarts 2014, abstract). It is a
haplotype-library method with a hard-genotype emission, the same family as AlphaPlantImpute2.

### 6.5 The Practical Haplotype Graph (Bradbury 2022)

The PHG stores founder assemblies as haplotype nodes per reference range and imputes a
low-coverage sample by mapping its reads to the nodes and running an HMM whose states are the
nodes (one path for inbreds, two for diploids) with a recombination transition between adjacent
ranges (Bradbury 2022, Bioinformatics 38:3698, doi 10.1093/bioinformatics/btac410). Reported:
< 0.01 genotype error on maize NAM RILs, and SNP error "approaching one in a million" in
simulation with 26 assemblies and 5,000 RILs; the paper does not tabulate accuracy against
coverage. It is the founder-known HMM with the founders as full assemblies and the reads as the
emission — RABBIT's `depModel` with a pangenome for its states.

### 6.6 Scaling to thousands of individuals

Cost per chromosome is individuals × markers × (states × jump-targets). For a biparental F2 that
is N × M × 4 × 4; for 2,000 individuals and 200,000 markers it is 6·10⁹ multiply-adds, seconds.
GBScleanR's Table 3: 17 s for 1,000 F2 × 620 markers with four estimation cycles, magicImpute 108 s,
LB-Impute 2,086 s; Furuta 2023 attributes LB-Impute's cost to its path enumeration. The founder
step is what grows: GBScleanR's is O(markers × founder-patterns² × individuals × states)
(`src/gbsrViterbi.cpp:286-345`), fine for 2 founders (2 patterns) and prohibitive for 8 outbred
founders (65,534 patterns, `HMM.R:510-537`). Sthapit 2025 ran STITCH on 9,780 wheatgrass genets
at 0.05× with K = 8–28 (Plant Genome, doi 10.1002/tpg2.70139), so the K-haplotype EM copes with
thousands of samples when K is small. AlphaImpute2's particle-PBWT is the answer at hundreds of
thousands, and is only needed there.

## 7. Inbred lines and doubled haploids

When heterozygosity is expected to be zero the model halves: a state is one founder, not a pair,
and a homozygous individual's reads test only "which founder". STITCH's `diploid-inbred` mode
allocates K states instead of K² (`STITCH/STITCH/R/functions.R:3744-3753`), runs the haploid
forward–backward (`functions.R:3289-3292`), and writes the diploid genotype posterior as
(1 − dosage, 0, dosage) — the heterozygote gets probability exactly 0
(`functions.R:3633-3643`). The transition keeps only "stay" and "jump" (`functions.R:4026-4043`).
RABBIT's `depModel` is the same reduction with F states (`MagicModel.m:143-149`); the manual says
"we may set model to depModel for a homozygous population" (RABBIT manual p. 5). AlphaPlantImpute2's
`-haploid` keeps one haplotype per library individual and sets every heterozygous call in the
input to missing with a warning (`alphaplantimpute2.py:299-309`). FILLIN's inbred mode is a
single-donor match. GBScleanR has no haploid mode: an inbred population is expressed as more
selfing generations in the scheme, which drives the heterozygous states' prior toward zero without
removing them (`HMM.R:846`), and `het_parent = FALSE` collapses founder haplotype labels to
founder identities (`HMM.R:599-615`).

What changes in practice:

- **Evidence per state doubles in effect.** Every read is a vote between founders; there is no
  heterozygous state to rule out, so the ten-read requirement of §3.5 disappears and two
  concordant reads per block suffice. This is why magicImpute's break point is 0.053× for the RIL
  and 0.11× for the F2 (Zheng 2018, Table 1), and why Xie 2010 could work at 0.055× and Huang
  2009 at 0.02×.
- **Residual heterozygosity is the failure mode of the haploid model.** A RIL at F6 is
  heterozygous over about 3% of its genome (1/32), in blocks. Under a haploid model such a block
  reads as a run of rapid founder switches or as errors, and the true heterozygote is never
  output (STITCH: P(het) = 0). Huang 2014 states its method "is not recommended for populations
  with high heterozygosity". The diploid model with a selfing-shrunk heterozygous prior (GBScleanR,
  RABBIT `jointModel`) keeps the heterozygote available at a cost of F² states — trivial for two
  founders.
- **Doubled haploids are the clean haploid case** (0 expected heterozygosity, one meiosis), and a
  DH population is where the haploid model is exact. NOISYmputer's authors list DH and RIL as
  future work (Triay 2025); nothing in the family is DH-specific because nothing needs to be.
- **Inbred germplasm** (a panel of lines, not a cross) is the haploid model with an unknown and
  large number of founders — FILLIN's regime, and the bridge to the many-haplotype report: a line
  is one haplotype, the panel of lines is its own reference, and imputation is Li & Stephens
  with haploid states over the other lines.

## 8. Detecting that a cohort is few-haplotype

A module could compute these from a called cohort (genotype calls, cohort allele frequencies,
per-sample heterozygosity) before choosing a regime. Expectations below are Mendelian arithmetic
or reasoning marked as such; none of the tools above performs such a check — every one takes the
population type as input.

1. **Allele-frequency spectrum of segregating sites.** In an F2, RIL or DH from two inbred
   founders every polymorphic site has cohort allele frequency ½ up to sampling
   (SD = √(0.25/2N): 0.035 at N = 100) and segregation distortion; in a BC1, ¼ or ¾; in an
   F-founder MAGIC, multiples of 1/F. A natural cohort's folded spectrum is dense near 0 and
   thin near ½: under a neutral 1/x spectrum, the fraction of SNPs with minor-allele frequency
   above 0.35 in 100 individuals is about ln(0.5/0.35)/ln(100) ≈ 8% (own reasoning). Statistic:
   the fraction of polymorphic sites with cohort allele frequency in [0.35, 0.65] (or near ¼/¾).
   ng's `cohort_expected_copies` gives this without genotype calls, so it works at 1×.
2. **Per-sample heterozygosity, and its uniformity.** F2 and BC1: ½ of polymorphic sites
   heterozygous in every individual; RIL F_g: 1/2^(g−1); DH: 0; natural inbred panel: low but
   varying by line; outbred natural population: 2pq averaged, typically 0.2–0.4 of segregating
   sites. At 1–3× the per-sample estimate must come from the genotype posteriors, not from hard
   calls (a hard call under-calls heterozygotes by the false-homozygote rate of §1).
3. **Pairwise identity over long windows.** For two individuals of a biparental population,
   over a 1 Mb window (a fraction of a block), their genotypes are identical, opposite, or one is
   heterozygous — so the distribution over pairs of "fraction of sites where the two are
   identical" is trimodal at 1, 0, ½ (with heterozygotes) — whereas in a natural population it is
   unimodal around 1 − 2·mean(2pq). This needs windowed genotypes; at 1× use the per-window
   majority call per individual (pool reads over the window, as Huang 2009 does), which is
   reliable exactly when the cohort *is* few-haplotype.
4. **Distinct haplotypes per window among homozygous individuals.** Cluster the windowed
   homozygous genotype vectors; count clusters that hold ≥ 2 individuals. Biparental: ≤ 2 plus
   recombinants (a recombinant inside the window looks like a third haplotype, so the count is
   2 + number of crossovers in the window, small for a 1 Mb window at 1–2 cM). Natural: grows with
   N. Statistic: the fraction of windows where the two largest clusters cover ≥ 95% of homozygous
   individuals.
5. **LD decay.** Between two sites at recombination fraction c, an F2 has r² = (1 − 2c)² (own
   derivation from the 1 : 2 : 1 gametic table), so r² ≈ 0.96 at 1 Mb ≈ 1 cM, and a RIL at
   F∞ has r² = ((1 − 2c)/(1 + 2c))² ≈ 0.92; natural populations are near 0 at 1 Mb. Statistic:
   median r² between site pairs 500 kb–1 Mb apart. Needs genotypes; at low coverage compute it
   from genotype dosages (posterior means), which are unbiased where hard calls are not.
6. **The decisive test: does a two-haplotype model explain the reads?** Run the K = 2 haplotype
   EM of §4.5 over a few windows and compare the total read log-likelihood with K = 8 and with
   the independent-sites model (each site's reads under the cohort allele frequency alone). In a
   biparental cohort K = 2 gains a large amount per read over independent sites and K = 8 gains
   almost nothing more; in a natural cohort K = 2 gains little and K = 8 much more. This is the
   same computation the module would run anyway, so the check costs a few windows of it.
   Reference point: Pierotti 2024 got r² = 0.996 with K = 16 on ten crosses of eight founders,
   i.e. K equal to the true founder count.

A cohort that passes 1, 2 and 6 is a biparental population; the crossing scheme (F2 versus RIL
versus BC) is then read off 1 and 2. What no statistic can supply is the *number of generations*
beyond what heterozygosity implies (an F6 RIL and an F8 RIL differ by 2 percentage points of
heterozygosity); it should be declared by the user and checked, not inferred.

## 9. Published accuracy and the failure modes

All numbers are the authors' own measurements on the data stated; "accuracy" is the fraction of
imputed calls that match the masked or simulated truth unless stated. Coverage is per sample.

| Method | Data | Coverage | Result | Where |
|---|---|---|---|---|
| Huang 2009 (window) | 150 rice RILs F11, parents sequenced, 1 SNP/40 kb | 0.02× | "expected" 99.94% (computed, not measured); breakpoints to 40 kb; 2,334 bins | Results |
| Xie 2010 (MPR + HMM) | 238 rice RILs, parents unknown | 0.055× | 11,792/11,792 inferred parental calls confirmed on chr 5; ≥ 110 RILs for 99% (Monte-Carlo) | Results; SI |
| TIGER | 384 Arabidopsis F2, 261,795 markers | 0.1×, 1×, 10× (simulated 1,000 each); real median 0.6× | 97.5 / 98.7 / 99.3% of crossovers found; median resolution 1,986 / 938 / 0 bp; 90% of crossovers within 27 / 4 / 1 kb; errors 90% false homozygotes, 2.4% in pericentromeres at 0.1× | Tables 2–3 |
| LB-Impute | simulated F2 and BC1, 200 offspring, 10,000 markers / 100 Mb | 0.1×–2.5× | > 99% at every coverage (window 7); 67–71% at window 2; 77–90% of breakpoints correctly flanked; Beagle 4.1: 37–79%; Mendel Impute: 50–99% | Fig. 2, S3, S4, S6 |
| LB-Impute | real maize B73 × CG F2 GBS (11,219 and 127,144 markers), IBM RIL (14,493 markers, 275 lines) | GBS, sub-1× | concordance with masked 7-read calls 94.6% (HincII), 91.4% (RsaI); RIL 97.0%; FSFHap 94.7 / 89.7 / 93.3% | Results, Fig. 3 |
| magicImpute | simulated 200 offspring: RIL, F2, 8-way MAGIC, outbred CP | 0.05×–4× | break point 0.053× / 0.11× / 0.21× / 0.21×; max 0.995 / 0.990 / 0.992 / 0.985; Beagle's break points 0.42× / 3.4× / 0.85× / 3.4× | Table 1, Fig. 2 |
| magicImpute | real maize RIL, maize F2, rice MAGIC, apple CP | down-sampled by 2^−i | 0.980 vs LB-Impute 0.917 (RIL); 0.987 vs 0.986 (F2); +2.5% over mpimpute (MAGIC); 0.94 (apple, 18% map inconsistencies) | Results |
| GBScleanR | simulated F2, outbred F1, 8-way RIL; 10–1,000 individuals; founders 5× | 0.1×–20× | F2 at 20× with error-prone markers: > 95% vs ~70%; +26 points from 10 to 1,000 samples in outbred F1 at 3×; no difference in 8-way RIL | Fig. 1, Supp. Data 1 |
| GBScleanR | real rice F2, 814 individuals, 5,035 SNPs | 0.85× | best on 11/12 chromosomes; +10.5 points over magicImpute on chr 7; double crossovers < 1 Mb: 14 vs 235 vs 948 | Fig. 3, Table 2 |
| NOISYmputer | simulated F2 (300, 66k–220k SNPs), rice F2 WGS (222), maize F2 GBS (91) | 0.5–3×; 3× vs 20× truth | breakpoint precision 99.9%, recall 99.6% (simulated, error ≤ 2%); rice 99% / 97%, median position error 415 bp; map length 208 cM vs LB-Impute 23,436 cM, FSFHap 337,750 cM on the same rice data | Tables 2–3 |
| STITCH (Pierotti 2024) | 2,177 medaka F2 from 8 inbred founders, 10 crosses, K = 16, nGen = 2 | 0.25×, 0.5×, 1×, 1.4× | r² 0.996 at 1.4×; 0.981 at 0.5×; plateau at 0.5× for 2,177, at 1× for 1,000–1,500, none below 500; single cross of 474: 0.964 (0.694 rarest SNPs) | Fig. 3 |
| STITCH (Sthapit 2025) | 9,780 outcrossing wheatgrass genets, K = 8 | 0.01×–0.10×, pipeline at 0.05× | concordance > 0.97, r² ≈ 0.90 at info > 0.8; most gain from 0.01× to 0.04× | Results |
| AlphaPlantImpute2 (Niehoff 2022) | 1,357 sugar-beet lines, arrays 19,519 HD / 950 LD | arrays | 0.90 (tuned Beagle 0.90, untuned 0.22); F1 from DH parents 1.00 vs Beagle 0.92; flat 0.90 from 5 to 50 offspring per family | Results |
| PHG | maize NAM RILs; simulated 26 assemblies × 5,000 RILs | skim (unstated) | genotype error < 0.01; SNP error approaching 10⁻⁶ in simulation | Results |
| RABBIT (arrays) | simulated CC/MAGIC; 703 real MAGIC lines | arrays | wrongly-called probability 0.008–0.105 vs HAPPY 0.027–0.225; 93.5% concordance with HAPPY | Zheng 2015 |

Failure modes the papers name:

- **Near recombination breakpoints.** Every tool is least certain there: LB-Impute leaves 2–150
  markers missing per event (Fig. S3B); Niehoff 2022 says both tools "struggle imputing regions
  around crossovers"; NOISYmputer leaves the support interval unimputed on purpose; Pierotti
  2024's rarest SNPs are the ones whose founder differs from the rest. The interval width is the
  informative-read spacing of §3.5.
- **Regions identical between founders.** No information; the state carries over from the
  flanks, which is right for genotypes and wrong for breakpoints. Huang 2014: "lower accuracy
  with genetically similar founders".
- **Heterozygous regions of RILs.** Lost by any pipeline that hard-calls under homozygosity
  (Zheng 2018, limitation 5) or runs haploid (STITCH `diploid-inbred`, AlphaPlantImpute2
  `-haploid`).
- **Error-prone markers.** Accuracy that *falls* with depth (Furuta 2023, Fig. 1); spurious
  short double crossovers (Table 2); NOISYmputer's degradation above 5% error.
- **Segregation distortion and map errors.** magicImpute assumes neither (Zheng 2018); apple CP's
  ceiling of 0.88–0.94 came from 18.3% map-inconsistent markers.
- **Pericentromeres.** TIGER's 2.4% false-homozygote rate at 0.1× where marker density drops
  abruptly (Rowan 2015).
- **Small cohorts** when founders are unknown (§4.7).

## 10. Relation to Li & Stephens

Is the biparental HMM Li & Stephens with two templates? Yes, in the founder-known, inbred-founder
case, term for term:

- *Templates.* Li & Stephens copies from a panel of K haplotypes; here K = 2 and the panel is
  the two founders. STITCH is exactly Li & Stephens with K unknown haplotypes learned by EM and
  is used on biparental data at K = 2–16.
- *Transition.* Li & Stephens jumps with probability 1 − e^(−ρd) and lands on a uniformly chosen
  template (including the current one). With K = 2 the probability of actually changing template
  is (1 − e^(−ρd))/2 — which is LB-Impute's eq. 4, `0.5·(1 − e^(−d/D))`
  (`FindPath2.java:103`), with ρ = 1/D. The diploid product `σ², σ(1−σ), (1−σ)²` in STITCH
  (`functions.R:4026-4043`) and the `(1−r)², r(1−r), r²` mixture in AlphaPlantImpute2 (wheel
  `CombinedHMM.py:371-383`) are the two-independent-chromosomes form of the same thing.
- *Recombination scale.* Li & Stephens uses ρ = 4N_e·r because the templates are a sample from a
  population whose common ancestors are ~2N_e generations back; the biparental model uses
  ρ = (number of meioses) × r — 2 for an F2, growing with each selfing generation as the map
  expands. STITCH's `nGen` is literally this substitution (`genetic-map.R:28-42`; help:
  "nGen = 4·Ne/K if unsure").
- *Emission.* Li & Stephens allows a mismatch to the template with a small probability θ
  (mutation since the common ancestor). The biparental emission has no mutation term: a mismatch
  is a read error or a marker error, hence the ε and w_m parameters instead of θ. Numerically they
  play the same role; conceptually θ scales with panel size and ε does not.

Where the identity breaks:

1. **The two homologues are not independent after selfing.** Li & Stephens (and STITCH,
   AlphaPlantImpute2) treat a diploid as two independent copying processes. In a RIL the two
   homologues are identical by descent over most of the genome and their junctions are shared; a
   crossover junction inherited from the F2 ancestor appears on *both* homologues at once. RABBIT's
   `jointModel` and GBScleanR's junction recursion model this coupling with an IBD/non-IBD
   component and simultaneous-jump rates (`MagicModel.m:131-133`; `HMM.R:936`). The consequence
   is a correct heterozygosity prior per generation, which the independent model cannot give.
2. **No template ambiguity, no mutation term, no state-space reduction.** With founders known
   the panel is fixed and tiny; nothing in Li & Stephens' machinery for large panels (PBWT,
   state-space collapsing, panel subsetting — the concerns of the many-haplotype report) is
   needed. With founders unknown the model becomes STITCH: the panel is a parameter, labels are
   arbitrary, EM is needed, and cohort size becomes the limiting resource.
3. **A prior on genotype frequencies.** Li & Stephens has a uniform initial state. The biparental
   chain's stationary distribution is 1 : 2 : 1, 1 : 1, or the inbreeding-weighted distribution
   RABBIT computes (`magicStationaryProbXY`), and it matters at chromosome ends and in
   read-poor regions.
4. **Founders are inferred jointly across the cohort**, which Li & Stephens never does: the
   founder step of GBScleanR and magicImpute is a second HMM over founder patterns whose
   emission is the product of all progeny likelihoods (§4.3–4.4).

## 11. What this implies for ng

ng's position: per-sample genotype likelihoods at every locus, per-allele read counts per
sample, cohort allele frequencies, a fitted per-sample inbreeding coefficient, and (in memory
during calling, not in the record) per-read allele evidence; the cohort may be 100–2,000
progeny at 0.5–3×, with parents present at high depth or not at all (`tmp/phasing_research/ng_facts.md`).

1. **The module is a founder-mosaic HMM over the called cohort, and the GLs are its input.** The
   emission at each site for founder-pair state (a, b) is the sample's GL for the genotype
   founders a and b imply there. Nothing else from the reads is needed; chain ids are not needed
   for this regime. ng should carry PL/GL (or the three log-likelihoods) into whatever
   representation the second pass reads — today they are not part of `LocusInference`
   (`ng_facts.md`).
2. **State space and transitions.** For two inbred founders: four ordered states (AA, AB, BA, BB)
   or three unordered; transitions from the crossing scheme by the junction recursion
   (GBScleanR's `.calcNextJnum`/`.getXoFreq` are 50 lines of R and port directly), exponentiated
   over the genetic distance; without a genetic map, physical distance times a user-supplied
   cM/Mb. This gives F2, BCn, RIL F_g, DH and NAM families from one implementation; MAGIC with
   inbred founders is the same code with F² states and 2^F − 2 founder patterns.
3. **Three founder regimes, one code path.** Founders are always hidden variables scored per site
   by (their own GLs, if sequenced) × (product over progeny of the progeny likelihood under each
   founder pattern), decoded by Viterbi along the chromosome, then progeny by forward–backward
   conditioned on the founders — GBScleanR's design (§4.3). Deep parents: their GLs dominate.
   Shallow parents: they contribute what they have. Absent parents: two flat founders; the
   progeny alone decide, labels are arbitrary per chromosome, and the module must say so in the
   output (a per-chromosome "founder label is arbitrary" flag). Below about 100 progeny with
   absent parents, report that founder inference is under-determined rather than emit a fitted
   answer (Xie 2010's ≥ 110; Pierotti 2024's < 200 cliff).
4. **Per-marker bias and mismapping, estimated from ng's own counts.** GBScleanR's two estimators
   need exactly what ng has per (allele, sample): reference and alternate read counts. Two to
   four estimate/decode cycles; clamp to [ε, 1 − ε]. This is the single largest accuracy lever in
   the literature (§5.3) and it is cheap.
5. **Residual heterozygosity is kept, not suppressed.** Diploid states always; the selfing
   generations shrink the heterozygous prior. A DH population is declared as such and gets the
   haploid reduction as an optimisation, not a different model.
6. **Outputs.** Per sample per site: genotype posterior (GP) and dosage; phased GT with one
   phase set per chromosome (the founder assignment *is* the phase); per sample: founder-origin
   blocks with breakpoint intervals (last informative read of one state to first of the next);
   per cohort: bins in Huang 2009's sense if a user wants them. Sites where the founders are
   inferred identical get GP from the founder genotype with no uncertainty reduction claimed.
7. **Cost.** Progeny decoding is N × M × 16 per chromosome — 2,000 × 200,000 × 16 = 6·10⁹
   operations, seconds; founder Viterbi is M × patterns² × N × states, trivial at two founders.
   Memory is the GLs: 2,000 samples × 200,000 sites × 3 floats = 4.8 GB per chromosome if held
   at once, or a checkpointed forward–backward over windows.
8. **Placement.** A second pass over the written cohort per chromosome (the streaming design
   forbids coupling in-flight segments; `run_streaming.md` §4.3 as quoted in `ng_facts.md`).
   The pass needs sites in genome order across all samples, which the VCF or a columnar sidecar
   of GLs provides.
9. **Regime detection.** Compute §8's statistics 1, 2 and 6 by default when a crossing scheme
   is declared, and warn when they disagree with it (an "F2" whose heterozygosity is 0.05 is a
   RIL; a "biparental" cohort whose K = 2 fit is no better than independent sites is not
   biparental). Do not infer the scheme silently.
10. **Multi-allelic loci and STRs.** A biparental cohort has at most two founder alleles per
    site; a third allele in the cohort is a mutation or an error. The module should run on the
    biallelic-between-founders sites and treat others as missing for the HMM, then impute them
    from the founder assignment like any other site.

## 12. Questions to discuss

1. **Dedicated founder-mosaic HMM, or the general K-haplotype EM at small K?** Recommendation:
   dedicated. The general EM (STITCH-like) is the many-haplotype report's business and can be
   pointed at a biparental cohort, but it cannot use the crossing scheme (wrong heterozygosity
   prior, independent homologues), cannot take deep parents as fixed founders, and needs the
   cohort to be large. The dedicated model is a few hundred lines and is exact where it applies.
   Trade-off: two HMM implementations to maintain; mitigated if both share the forward–backward
   and the GL emission.
2. **What to do with absent parents and few progeny.** Recommendation: run founder inference
   only when N ≥ 100 (configurable); below it, emit the progeny genotypes unimputed with a
   diagnostic, or accept a user-supplied founder VCF. Trade-off: users with 40 RILs and no
   parents get nothing; they are also the users for whom every published method is unvalidated.
3. **Per-marker bias and mismapping: always on?** Recommendation: on by default, two cycles,
   with the fixed-parameter escape GBScleanR offers. Trade-off: two extra decodes per chromosome;
   the parameters are absorbed by noise in cohorts under ~50 individuals (Furuta 2023's 10-sample
   runs), so gate on N.
4. **Genetic map input.** Recommendation: accept an optional map; default to physical distance ×
   a per-genome cM/Mb (GBScleanR's 0.04 crossovers/Mb is a rice-scale default and would be wrong
   for tomato's ~1.5 cM/Mb euchromatin and near-zero pericentromeres). Trade-off: a uniform rate
   over-penalises crossovers in hotspots and under-penalises them in cold regions; it costs
   breakpoint placement, not genotype accuracy.
5. **Declare or detect the crossing scheme?** Recommendation: declare, then check with §8's
   statistics and refuse to run when they contradict the declaration by a wide margin.
   Trade-off: one more required input; the alternative (inference) cannot recover generation
   count and would fail silently on segregation distortion.
6. **Pedigree breeding populations.** Recommendation: out of scope for this module; treat each
   family as a biparental cross with its parents as founders, and let the many-haplotype module
   (haplotype library of all parents, FILLIN/AlphaPlantImpute2-style) cover individuals whose
   parents are untyped. Peeling over multi-generation pedigrees (AlphaImpute2) is a third
   machine and the literature's numbers for it are array-based.
7. **Read-backed evidence inside the founder model.** Recommendation: not needed for the
   founder-pair state, so do not retain chain ids for this regime. The only place it would add
   information is a read spanning two founder-differing sites inside a heterozygous block, which
   changes nothing the block already tells. Keep that decision separate from the many-haplotype
   module, where STITCH's per-read multi-site emission (`haploid.cpp:174-190`) does matter.
8. **Breakpoint output.** Recommendation: emit the interval (last read of the old state, first
   read of the new), not a point, and do not refine with filtered-out markers as TIGER does
   unless the marker filter is the module's own. Trade-off: users of bin maps want a point; give
   the midpoint as a second field.

## 13. References

Code (all under `/Users/jose/devel/pop_var_caller/tmp/phasing_research/repos/`, clones of 2026-09-12):
LB-Impute `SourceCode/LB-Impute/src/imputation/*.java`; GBScleanR 2.7.4 `R/Methods-GbsrGenotypeData_HMM.R`,
`src/gbsrCalcProb.cpp`, `src/gbsrViterbi.cpp`, `src/gbsrFB.cpp`, `man/estGeno.Rd`; RABBIT
`RABBIT_Packages/{MagicDefinition,MagicReconstruct,MagicOrigin,MagicImpute}/*.m`, `RABBIT_Manual/RABBIT Manual.pdf`;
STITCH `STITCH/R/functions.R`, `STITCH/R/genetic-map.R`, `STITCH/src/haploid.cpp`, `README.md`, `STITCH.R`;
AlphaPlantImpute2 `alphaplantimpute2/alphaplantimpute2.py` and the vendored wheel
`alphaplantimpute2-1.5.3-py3-none-any.whl` (`tinyhouse/CombinedHMM.py`, `tinyhouse/HaplotypeLibrary.py`);
AlphaImpute2 `src/alphaimpute2/alphaimpute2.py`, `Imputation/ParticlePhasing.py`, `Imputation/PhasingObjects.py`,
`peeling_notes.md`; TIGER `hmm_prob.pl`, `beta_mixture_model.R`, `refine_recombination_break.pl`, `sample_hmm_model`,
`README.md`. TASSEL FSFHap/FILLIN: GitHub mirror `zacharymiller90/tassel-ML`, `src/net/maizegenetics/analysis/imputation/`
(not cloned; read through the mirror).

Papers:

- Bradbury PJ, Casstevens T, Jensen SE, et al. 2022. The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. *Bioinformatics* 38:3698–3702. doi:10.1093/bioinformatics/btac410.
- Davies RW, Flint J, Myers S, Mott R. 2016. Rapid genotype imputation from sequence without reference panels. *Nature Genetics* 48:965–969. doi:10.1038/ng.3594. (Read via the companion paper Nicod J, Davies RW, et al. 2016, *Nature Genetics* 48:912–918, doi:10.1038/ng.3595, PMC4966644, whose Online Methods describe the model; the STITCH paper itself was not fetched.)
- Fragoso CA, Heffelfinger C, Zhao H, Dellaporta SL. 2016. Imputing genotypes in biallelic populations from low-coverage sequence data. *Genetics* 202:487–495. doi:10.1534/genetics.115.182071. (Full text in `tmp/phasing_research/few_haplotype/lb_impute_2016.txt`.)
- Furuta T, Yamamoto T, Ashikari M. 2023. GBScleanR: robust genotyping error correction using a hidden Markov model with error pattern recognition. *Genetics* 224:iyad055. doi:10.1093/genetics/iyad055. PMC10213493.
- Gonen S, Wimmer V, Gaynor RC, Byrne E, Gorjanc G, Hickey JM. 2018. A heuristic method for fast and accurate phasing and imputation of single nucleotide polymorphism data in bi-parental plant populations. *bioRxiv* doi:10.1101/330027.
- Huang X, Feng Q, Qian Q, et al. 2009. High-throughput genotyping by whole-genome resequencing. *Genome Research* 19:1068–1076. doi:10.1101/gr.089516.108. PMC2694477.
- Huang BE, Raghavan C, Mauleon R, Broman KW, Leung H. 2014. Efficient imputation of missing markers in low-coverage genotyping-by-sequencing data from multiparental crosses. *Genetics* 197:401–404. doi:10.1534/genetics.113.158014.
- Niehoff T, Pook T, Gholami M, Beissinger T. 2022. Imputation of low-density marker chip data in plant breeding: evaluation of methods based on sugar beet. *The Plant Genome* 15:e20257. doi:10.1002/tpg2.20257 (read as bioRxiv 10.1101/2022.03.29.486246).
- Pierotti MER, Welz L, Osuna-López F, Fitzgerald T, Wittbrodt J, Birney E. 2024. Genotype imputation in F2 crosses of inbred lines. *Bioinformatics Advances* 4:vbae107. doi:10.1093/bioadv/vbae107.
- Rowan BA, Patel V, Weigel D, Schneeberger K. 2015. Rapid and inexpensive whole-genome genotyping-by-sequencing for crossover localization and fine-scale genetic mapping. *G3* 5:385–398. doi:10.1534/g3.114.016501. PMC4349092.
- Sthapit SR, Crain J, Larson S, Anderson JA, Bajgain P, DeHaan LR, Poland J. 2025. A low-coverage skim-sequencing and imputation pipeline for genomic selection. *The Plant Genome* doi:10.1002/tpg2.70139. PMC12547641.
- Swarts K, Li H, Romero Navarro JA, et al. 2014. Novel methods to optimize genotypic imputation for low-coverage, next-generation sequence data in crop plants. *The Plant Genome* 7(3). doi:10.3835/plantgenome2014.05.0023. (Abstract and TASSEL source only; full text not reachable on 2026-09-12.)
- Triay C, Boizet A, Fragoso C, et al. 2025. Fast and accurate imputation of genotypes from noisy low-coverage sequencing data in bi-parental populations. *PLoS ONE* 20:e0314759. doi:10.1371/journal.pone.0314759. PMC11781708.
- Vi T, Stuart KC, Tan HZ, Lloret-Villas A, Santure AW. 2025. Assessing genotype imputation methods for low-coverage sequencing data in populations with differing relatedness and inbreeding levels. *Molecular Ecology Resources* 25:e70049. doi:10.1111/1755-0998.70049. (Abstract only.)
- Whalen A, Hickey JM. 2020. AlphaImpute2: fast and accurate pedigree and population based imputation for hundreds of thousands of individuals in livestock populations. *bioRxiv* doi:10.1101/2020.09.16.299677.
- Xie W, Feng Q, Yu H, et al. 2010. Parent-independent genotyping for constructing an ultrahigh-density linkage map based on population sequencing. *PNAS* 107:10578–10583. doi:10.1073/pnas.1005931107. PMC2890813.
- Zheng C, Boer MP, van Eeuwijk FA. 2015. Reconstruction of genome ancestry blocks in multiparental populations. *Genetics* 200:1073–1087. doi:10.1534/genetics.115.177873. PMC4574238.
- Zheng C, Boer MP, van Eeuwijk FA. 2018. Accurate genotype imputation in multiparental populations from low-coverage sequence. *Genetics* 210:71–82. doi:10.1534/genetics.118.300885. PMC6116951.
- de Haas et al. 2017. Low-coverage resequencing detects meiotic recombination pattern and features in tomato RILs. *DNA Research* doi:10.1093/dnares/dsx024. PMC5726486. (60 tomato F6 RILs at ~6×; sliding window; 96.8% KASP concordance — cited for the tomato reference point only.)
