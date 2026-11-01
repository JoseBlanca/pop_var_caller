# Using neighbouring loci to improve pvc's genotype calls

*Research report, 2026-11-01, measured on branch `neighbour-hmm` (not kept). It answers the owner's question of the same
day: can a hidden Markov model built for pvc — over the records pvc writes, fed with pvc's own
genotype likelihoods, multi-allelic and STR records included — use the evidence of surrounding loci
to correct weak genotypes, and give a call where a sample has little or no evidence? It is an
evaluation with a prototype, not a design. The phasing and imputation research it builds on is in
[`phasing_imputation/`](phasing_imputation/00_overview.md); the 2026-09-13 decision it revisits is
[`07_decision_separate_tool_and_route.md`](phasing_imputation/07_decision_separate_tool_and_route.md).*

---

## 0. The answers

**Yes for an inbred cohort, no for a single outbred sample.** Measured on tomato1 by cutting nine deep
accessions to low depth and comparing with their full-depth calls, and on a simulated cohort where the
truth is known:

| signal | one sample | cohort | at 3× | at 30× | tomato | HG002 |
|---|---|---|---|---|---|---|
| shared haplotypes | no templates: nothing | removes 36–56% of genotype errors in tomato | large gain | nothing left to fix | large | not available (one sample, no panel) |
| local autozygosity | removes 5% (tomato, 3×); estimates the *F* pvc defaults to 0 | removes 7–25% in tomato, half the heterozygote errors | moderate | nothing left to fix | moderate | none: an outbred genome has no runs |
| read linkage | not measured | not measured | 1 adjacent pair in 5 linkable | reach, not gain | — | reachable misses are mapping problems |

The figures, each with its condition, are in §4–§6. In brief:

- **Shared haplotypes give a call where the sample has no evidence.** In tomato at 3×, true
  non-reference genotypes at sites where the shallow accession had no read go from 85% wrong to 26%
  (§5.2). The model needs a cohort, and in the simple form tested here — other samples' genotypes as
  templates — it is valid only for inbred material; in a simulated outbred cohort it gains nothing
  (§4).
- **Local autozygosity fixes heterozygotes in inbred samples.** pvc's genome-wide *F* of about 0.95
  makes a true heterozygote in a residual heterozygous block lose to the homozygote 43% of the time
  at 3× even with four or more reads; a two-state model of where the sample's runs are halves that
  (§5.2).
- **HG002 at 30× gains nothing from any of the three.** Its missed SNPs with 1–3 variant reads are
  mostly sites where one haplotype's reads are missing (2 of 19 reads, not 2 of 4), and those with no
  variant read need a panel (§6).
- **Both models change genotypes only at records pvc writes.** A locus pvc never built — 42% of true
  non-reference genotypes for one accession alone at 3× — is beyond any second pass (§5.1, §5.2).
- **Where it sits:** a second pass over pvc's records, reading `PL` and `GP`, which pvc does not
  write yet; STR and multi-allelic records need nothing special (§8).
- **Recommendation:** write `PL` and `GP`; build local autozygosity first, then extend it with
  genotype templates for cohorts; leave read linkage. The 2026-09-13 decision keeps its shape — a
  second pass with the VCF as interface — but the pass becomes pvc's own (§9).

---

## 1. Three signals that link neighbouring sites

A genotype call at one site uses the reads at that site and a prior. Three things make the sites
around it informative too, and each needs different data.

- **Shared haplotypes.** Two samples that carry the same allele at a site usually inherited it, with
  the surrounding stretch of chromosome, from one ancestor. So a sample's chromosome is, stretch by
  stretch, a copy of some other chromosome in the population, and the copied chromosome says what the
  sample carries where its own reads are silent. This is the signal of every phasing and imputation
  tool ([`00`](phasing_imputation/00_overview.md) §1). It needs *templates*: other chromosomes to copy
  from — the cohort itself or an external panel. A single sample without a panel has none.
- **Local autozygosity.** In an inbred sample the two copies of a chromosome are, over long runs,
  one ancestral copy inherited twice — *identical by descent*, so the sample is homozygous along the
  whole run. Between runs the copies are independent and the sample can be heterozygous. Where the
  sample is inside such a run can be read off its own likelihoods at the neighbouring sites: a stretch
  where every confident call is homozygous is probably a run. This needs nothing but the one sample.
  The research set does not treat it; it is the one signal here that is new.
- **Read linkage.** A read or read pair that covers two sites shows which alleles sit on one
  molecule. Its reach is one insert; at 3× one adjacent heterozygous pair in five can be linked at all
  ([`00`](phasing_imputation/00_overview.md) §2.5, [`05`](phasing_imputation/05_read_backed_phasing.md)
  §2).

**Where pvc already uses the second signal, genome-wide.** pvc's genotype prior is a mixture over two
cases: with probability *F* the sample's two copies are one copy counted twice, and the genotype is a
homozygote; with probability 1 − *F* they are independent draws from the cohort's allele frequencies
(`src/calling/genotype_prior/dirichlet_multinomial.rs`, `MarginalizedDirichletPrior`). *F* is the
sample's inbreeding coefficient, fitted once over the whole genome. A local autozygosity model keeps
this prior exactly and replaces the one *F* by the probability that *this site* lies in a run, given
the sample's evidence at the sites around it.

---

## 2. What pvc has today, and what a second pass would read

- **pvc writes neither `PL` nor `GP`.** The `FORMAT` fields are `GT:GQ:DP:AD`, plus `REPCN` and
  `REPCNB` at repeat tracts (`src/vcf/header.rs`, `FORMAT_DECLARATIONS`). 07 §1 records that pvc
  owes both; that is still owed.
- **The likelihoods exist, for every record and every path.** In the final calling pass,
  `summarise_final_pass` (`src/calling/inference/summarise_condition.rs`) scores each sample with its
  likelihood row, its prior row and its posterior row all in hand — SNP/indel and STR records alike,
  in VCF genotype order (0/0, 0/1, 1/1, 0/2, …; `src/calling/genotype_table.rs`). They are discarded
  after the call.
- **Everything the prior needs can travel in the VCF.** The prior is rebuilt from the sample's
  *leave-one-out concentration*: per allele, the run's seed plus the copies the *other* samples carry.
  That is the seed (in the parameters file), the cohort's fitted frequency (`AF`, already written)
  and the sample's own expected copies (from `GP`). So `PL` and `GP` are all a second pass needs.
- **The stage is a second pass, as 07 §1 said.** A site's posterior depends on sites on both sides,
  and pvc hands each record over as soon as it is called (`doc/devel/ng/spec/run_streaming.md` §4.3).
  It can be a separate command over the VCF, or a buffered stage that holds a window of written
  records (a few megabases) before writing them.

**The prototype's change to pvc, not kept.** When the environment variable `PVC_PROTOTYPE_ROWS`
named a file, `assemble_record` (`src/vcf/assemble.rs`) appended each record's per-sample
natural-log likelihood rows and leave-one-out concentrations to it. About 100 lines, saved as
[`neighbour_evidence/prototype_rows.patch`](neighbour_evidence/prototype_rows.patch) against main
8578ce52 so the scorer's input can be rebuilt; `PL` and `GP` in the VCF replace it.

---

## 3. The prototype models

Three ways to call one sample's genotype at a record, all from the same likelihood rows and the same
two branches of pvc's prior ([`neighbour_evidence/hmm.py`](neighbour_evidence/hmm.py)):

1. **pvc today.** The mixture with the sample's genome-wide *F*.
2. **Local autozygosity.** A two-state HMM along the sample's records: *in a run* or *not in a run*.
   At each record the evidence for each state is the sample's likelihood averaged over that branch of
   the prior. Between records a state persists with probability exp(−distance / mean run length);
   the run length is chosen by likelihood from 10 kb to 3 Mb, and the share of the genome in runs
   starts at the fitted *F* and is re-estimated. The call at a record mixes the two branches with
   that record's probability of being in a run. Two parameters a sample; works at one sample.
3. **Shared haplotypes.** Each copy of the sample copies one of *K* templates, switching template
   with probability 1 − exp(−distance / segment length). A template is another sample's
   *allele-probability vector* at that record — its expected share of each allele under pvc's own
   call — so a record with five alleles or an STR with eight lengths is handled with no change. The
   autozygous state is kept: in it both copies copy one template. The *K* = 16 templates are the
   samples whose calls agree best with the target's over the region. This is the Li & Stephens model
   of 07 §2 with two simplifications: templates are genotypes rather than phased haplotypes, and the
   sample is scored as a diploid over template pairs (K² + K states) in one pass, with no sweeps.

The genotype-as-template simplification is sound for inbred material, where each accession is nearly
one haplotype, and is expected to fail for outbred samples, where a template that is heterozygous is a
blur of two haplotypes ([`04`](phasing_imputation/04_few_haplotype_populations.md) §7: "a line is one
haplotype, the panel of lines is its own reference"). §4 confirms it.

---

## 4. Simulation: where each signal can help, with the truth known

**What was simulated** ([`neighbour_evidence/sim.py`](neighbour_evidence/sim.py)). A 4 Mb region with a variable site every kilobase.
Every chromosome copy is a patchwork of founder haplotypes; inside autozygous runs the second copy
equals the first. Reads are drawn at a set depth, each from a random copy with 1% error, and the
likelihood rows are built as pvc builds a SNP's. The cohort's leave-one-out concentrations come from
four rounds of the cohort's expected allele copies, as in pvc's frequency loop; a lone sample gets
only a seed (total 10, alternative share 0.5%). Each sample's genome-wide *F* is fitted from its own
likelihoods, as pvc fits it. Three populations:

| population | founder haplotypes | shared segment | share of genome autozygous | mean run |
|---|---|---|---|---|
| inbred (tomato-like) | 8 | 300 kb | 0.8 | 2 Mb |
| outbred (human-like) | 400 | 30 kb | 0.01 | 300 kb |
| consanguineous (outbred, *F* = 1/16) | 400 | 30 kb | 0.0625 | 3 Mb |

**Result.** Wrong genotypes, summed over the scored samples (eight lone samples at one sample; four
samples of each cohort), out of the true non-reference genotypes in brackets:

| population | samples | depth | pvc today | + local autozygosity | + shared haplotypes |
|---|---|---|---|---|---|
| inbred | 1 | 3× | 2,309 (of 8,941) | 2,173 | — |
| inbred | 1 | 10× | 138 | 71 | — |
| inbred | 10 | 3× | 796 (of 4,901) | 651 | 432 |
| inbred | 63 | 1× | 1,976 (of 4,641) | 1,870 | 960 |
| inbred | 63 | 3× | 721 (of 4,901) | 537 | 198 |
| inbred | 63 | 10× | 85 (of 4,849) | 38 | 4 |
| inbred | any | 30× | 0–1 | 0–1 | 0–1 |
| outbred | 1, 10, 63 | 1×–10× | | within 5% of pvc today | within 8% |
| consanguineous | 1, 10, 63 | 1×–10× | | within 1% | within 6% |

What it shows:

- **In inbred material both signals work, and the haplotypes are worth more.** At 63 samples and 3×
  the shared-haplotype model removes 72% of the errors; autozygosity alone removes 26%. At one sample
  only autozygosity is available, and its gain is in heterozygotes: at 3× it cuts heterozygote errors
  from 573 to 369 of 622, while total errors fall only 6%, because most errors at one sample are
  homozygous alternatives with fewer than two reads, which the seed prior holds at reference and no
  neighbour can lift.
- **In outbred material neither moves anything,** at any depth or cohort size. Autozygous runs cover 1%
  of an outbred genome, and genotype templates of outbred samples are blurs of two haplotypes. Getting
  anything from shared haplotypes there needs phased templates, sweeps and a large cohort or panel —
  the full machinery of 07 §2, measured in the literature at hundreds of samples
  ([`06`](phasing_imputation/06_accuracy_evidence_and_calling_interplay.md) §4).
- **At 30× there is nothing left to fix**, in any population.

---

## 5. Tomato: held-out downsampling

### 5.1 What was done

tomato1 is 63 accessions over 80 regions of 100 kb (8 Mb of SL4.0). **It is not a 3× cohort**: depths,
counted from the read files, run from 2.6× to 30× with a median of about 11×, and nine accessions
are at 18× or more. **Six plants were sequenced twice** (two read files with one sample name); pvc
refuses a cohort that holds an individual twice, so one file of each pair was dropped, leaving 57.
The nine deep accessions are the truth set.

| condition | what was called |
|---|---|
| truth | all 57 at native depth; a deep accession's call is truth where pvc's posterior is at least 0.99 |
| nine at 3× | the nine cut to 3×, the other 48 at native depth |
| nine at 1× | the nine cut to 1×, the others native |
| cohort at 3× | every accession above 3× cut to 3× |
| cohort at 1× | every accession above 1× cut to 1× |
| one sample | each of the nine called alone, at 3× and at 1× |

Parameters were fitted separately in each condition, as a real run would fit them. The nine deep
accessions' fitted *F* is 0.82 to 0.99 (cohort median 0.98). Truth is **pvc's own call at 20–30×**,
not an independent truth set: the measurement is how close each model brings a shallow call to what
pvc says at full depth.

**Sites the scorer cannot match.** A truth locus with no record of the same boundaries in the test
call is counted as called reference by every model alike, and reported on its own line. At the nine
cut to 3×, 13,262 truth loci are missing: 4,483 overlap a record whose boundaries moved, and the other
8,779 were never built, because once the nine were cut no sample showed enough evidence. **A second
pass refines records pvc wrote, so a locus pvc never built is out of its reach.** The tables below are
over the matched sites.

### 5.2 Results

Wrong genotypes over the nine accessions, matched sites only (about 187,000 truth sites an accession
at 3×):

| condition | pvc today | + local autozygosity | + shared haplotypes | truth loci not built (non-ref genotypes) |
|---|---|---|---|---|
| nine at 3×, cohort native | 49,736 | 41,217 (−17%) | 27,082 (−46%) | 19,690 |
| nine at 1×, cohort native | 128,179 | 119,027 (−7%) | 56,723 (−56%) | 28,840 |
| cohort at 3× | 50,427 | 37,586 (−25%) | 26,453 (−48%) | 39,920 |
| cohort at 1× | 96,454 | 87,300 (−9%) | 62,008 (−36%) | 87,695 |
| one sample at 3× | 3,608 | 3,422 (−5%) | — | 124,132 |
| one sample at 1× | 8,625 | 8,633 (0%) | — | 239,932 |

By the evidence the shallow sample had at the site (nine at 3×, cohort native):

| sample's reads at the site | true non-ref genotypes | wrong today | + autozygosity | + haplotypes |
|---|---|---|---|---|
| none | 21,962 | 18,721 (85%) | 18,651 | 5,637 (26%) |
| 1–3 | 149,979 | 17,447 (12%) | 14,088 | 11,429 (8%) |
| 4 or more | 104,282 | 12,184 (12%) | 6,138 | 5,788 (6%) |

| true heterozygotes | count | wrong today | + autozygosity | + haplotypes |
|---|---|---|---|---|
| sample has 1–3 reads | 15,671 | 13,405 (86%) | 8,097 | 8,050 |
| sample has 4 or more reads | 27,750 | 11,957 (43%) | 5,693 | 5,425 |

What they show:

- **Shared haplotypes give a call where the sample has no evidence.** Of true non-reference genotypes at
  sites where the shallow sample has no read, pvc gets 85% wrong — it can only say reference — and the
  copying model gets 26% wrong. That is the owner's "call at a site with little or no evidence", and on
  tomato it works: three in four are recovered at 3× and the same share at 1× inside a native-depth
  cohort (87% wrong today, 27% with templates).
- **Local autozygosity fixes heterozygotes.** With *F* about 0.95, pvc's prior holds a sample's
  heterozygote at a twentieth of its outbred weight everywhere, so a true heterozygote in one of the
  residual heterozygous blocks is called homozygous 43% of the time even with four or more reads at
  3×. The local model sees the block and halves that. The copying model, which contains the same
  autozygous state, does no better on heterozygotes; its extra gain is all at sites with little or
  no evidence.
- **The cost is false non-reference calls**, and it is small against the gain. Truth-reference sites
  called non-reference at the nine-at-3× condition: 1,382 today, 2,339 with autozygosity, 4,227 with
  haplotypes — about 2,800 new errors for about 25,500 removed.
- **At whole-cohort 1× the templates are themselves noisy**, and the copying model is worse than
  autozygosity alone at sites the sample has reads for (22,408 against 20,111 wrong at 1–3 reads),
  while still recovering half the non-reference genotypes at sites with no read (33,715 wrong of
  71,821, against 63,047 today). And 30% of the true non-reference genotypes are at loci pvc never
  built at 1×.
- **At one sample the gain is small, and the loss is elsewhere.** pvc cannot fit *F* from one sample —
  the fit needs the cohort's allele frequencies — so each accession called alone is scored as an
  outcrosser (*F* = 0, marked `defaulted` in its parameters file). The local model estimates 0.61 to
  0.72 from the accession's own likelihoods, and moves its errors from 3,608 to 3,422 at 3× and not at
  all at 1×. What dominates is the loci pvc never built: a lone accession gets a record only where its
  own reads show a variant, and 42% of the true non-reference genotypes at 3× (124,132 of 295,915),
  81% at 1×, have none. No second pass can reach those; the simulation's larger one-sample gain in
  heterozygotes (§4) does not reproduce here.
- **STR records behave like the rest**, through the same allele-probability templates: at the nine at
  3×, 1,965 wrong of 27,266 STR truth sites today, 1,850 with autozygosity, 1,381 with haplotypes.

**A separate defect found on the way.** At the cohort cut to 1×, pvc's hidden-duplication filter fitted
a duplication rate of 0.95 and dropped 113,041 of the 118,405 records it scored, writing 5,366. In the
other conditions it drops 12–16%. The prototype rows are taken before the filter, so the figures above
are unaffected; pvc's actual output at that depth is not. The cause was the filter's depth scatter σ₀,
which left out the noise of the reads themselves; the fix is on branch `dup-filter-depth` (0f712586,
not yet merged) and takes the drop at 1× to 12.5%
([`ng_dup_filter_read_sampling_2026-11-01.md`](../../reports/ng_dup_filter_read_sampling_2026-11-01.md)
on that branch).

---

## 6. HG002: the single human sample at 30×

The motivating numbers came from pvc on HG002 chr20 at 30×
(`../pop_var_caller-multiallelic-lik/tmp/generic_errors/`, run 2026-10-31): 423 missed SNPs, 196 with
no read carrying the variant at any mapping quality, 164 more with 1–3 variant reads at mapping quality
20 or more. What none of the three signals can do for them:

- **The 196 with no variant read** need templates, and a single sample without a panel has none (§1).
  An outbred sample also has almost no autozygous runs to borrow from (§4).
- **The 164 with 1–3 variant reads are mostly not short of reads.** 111 of them show fewer than 15%
  variant reads at a mean depth of 19 reads; only 15 have a variant share near a half, at a mean depth
  of 2.4. A heterozygote seen in 2 of 19 reads is a site where one haplotype's reads are missing or
  mis-mapped — the median depth at all 423 missed SNPs is 13, against 30 for the sample — and evidence
  from neighbouring sites cannot restore reads that did not map.
- **Read linkage could reach many of them, and would mostly confirm the imbalance.** 75 of the 164 have
  a correctly called heterozygous site within 150 bases (one read), and 125 within 500 (one read pair);
  at a random position the chance of a heterozygous site within 150 bases is about one in four, so the
  misses cluster. A read joining such a site to a 2-in-19 variant shows which haplotype the variant is
  on, but the reads of that haplotype that show reference are still there. Measuring what linkage buys
  needs per-read allele data that pvc drops at the cohort merge (07 §1, reason 4); this evaluation did
  not build it.

So at HG002 30×, the gain from neighbour evidence is near zero, for reasons the data shows: the misses
are mapping problems, not thin evidence.

---

## 7. What it costs

- **A second pass over the records**, holding a window of them. The autozygous runs and shared
  segments found here are tens to hundreds of kilobases; tomato1's regions are 100 kb, which capped
  what could be learned (the autozygosity fit chose runs of 10–100 kb on every sample, never longer
  than a window). A real pass would work over windows of a few megabases with an overlap.
- **Time, measured in the Python prototype**: per accession over tomato1's 202,489 records, about
  30 seconds for the autozygosity model (twelve runs of the two-state HMM for the run-length choice)
  and 35 seconds for the copying model with 16 templates (two runs). Compiled code would be far
  faster; not measured. Per sample per record, the autozygosity HMM is two numbers and the copying model
  K² + K = 272.
- **Memory, own arithmetic.** At tomato1's density, 26 records a kilobase, a 5 Mb window holds about
  130,000 records. Autozygosity needs one sample's rows at a time. The copying model needs every
  sample's allele shares in the window to choose templates: for 2,000 samples, 130,000 records and
  about two alleles a record, 520 million numbers — about 0.5 GB at one byte each.
- **Template choice scales with the cohort.** Choosing 16 of 56 by agreement over the window is
  trivial; at 2,000 samples it is the step every imputation tool spends its engineering on (07 §2.6).

---

## 8. Where it would sit, and how STR records enter

- **A second pass over pvc's records, reading `PL` and `GP`.** pvc writes neither today; both are owed
  already (07 §1). With them and `AF`, the pass rebuilds each sample's prior exactly (§2). It can be a
  separate pvc command over the VCF, or a buffered stage inside the run that holds a window of written
  records; either is a pass over emitted records, as `run_streaming.md` §4.3 requires.
- **What it changes is the genotype, and with it `GQ` and `GP`.** The record's alleles and site
  quality stay as pvc called them; a locus pvc did not build stays unbuilt.
- **STR and multi-allelic records need nothing special.** The likelihood rows are in VCF order for any
  number of alleles; a template is an allele-probability vector; the autozygous state's prior is the
  same branch pvc already computes. The prototype ran STR segments and bundles through the same code.
  This is the one place a pvc-specific model beats the external tools, which split or drop such
  records (07 §3).

---

## 9. Recommendation

1. **Write `PL` and `GP` first.** It is owed already, about 60–80 lines, and it is the interface the
   pass reads; without it nothing here can run outside this branch.
2. **Local autozygosity is the first model.** Two states and two numbers a sample; it needs nothing but
   the sample's own likelihoods; it replaces a number pvc already uses (*F*) with a local one, inside
   the prior pvc already has. In a tomato cohort at 3× it removes 17–25% of genotype errors and half
   the heterozygote errors. Alone on one accession it removes only 5% at 3×, but it also supplies the
   inbreeding estimate pvc defaults to zero there. On outbred samples it changes nothing, which is the
   correct answer: an outbred sample has no runs to find.
3. **Shared haplotypes, with genotypes as templates, second.** It is the larger gain on tomato — 36–56%
   of errors, three in four non-reference genotypes recovered where the sample has no read — but it
   needs a cohort, and its genotype-template form is only valid for inbred material. Its autozygous
   state is the first model, so it extends that model rather than replacing it. At one sample, or in an
   outbred cohort, it must be off and say so: there it either has no templates or would need phased
   ones.
4. **Read linkage last, or not at all.** At HG002 30× the reachable misses are mapping problems, and at
   3× one adjacent pair in five can be linked ([`05`](phasing_imputation/05_read_backed_phasing.md) §2).
5. **The 2026-09-13 decision should change in one respect.** Its reasons 1–3 stand: the stage is a
   second pass, and the VCF with `PL` is its interface. What changes is that the stage is pvc's own
   rather than an external tool's: the owner has ruled out GLIMPSE2 and Beagle, and the gain measured
   here comes from exactly what those tools cannot take — multi-allelic and STR records, per-sample
   inbreeding, and inbred lines as one-haplotype templates. The engine of 07 §2 is not needed for
   this: no phasing, no sweeps, no positional Burrows–Wheeler index.

**Before building, two measurements the prototype could not make.** Whole chromosomes rather than
100 kb windows, so run and segment lengths can be learned; and a truth that is not pvc's own deep call
— the two deep-plus-shallow replicate pairs in tomato1 (SRR7279542 with SRR7279488, SRR7279514 with
SRR7279529, the same plants sequenced deep and shallow) give a natural low-depth check that needs no
downsampling.

---

## 10. Limits of this evaluation, and how to reproduce it

- **Truth is pvc at 20–30×**, kept where its posterior is at least 0.99. Heterozygous calls in highly
  inbred accessions may include duplicated-region artefacts that appear at both depths; the
  measurement counts agreement with them.
- **The nine downsampled accessions share a cohort with 48 at native depth** in two of the conditions,
  which makes good templates; the whole-cohort conditions are the harder and more realistic test.
- **The models are simplified** (§3) and their few parameters were chosen by likelihood on small grids.
  Nothing was tuned on the truth.
- **The simulation is at the likelihood level**: reads are counts with one error rate, not alignments.

**To reproduce.** Everything is in [`neighbour_evidence/`](neighbour_evidence/):

- `prototype_rows.patch` — the pvc change that writes the likelihood rows; apply it to main 8578ce52
  and build.
- `hmm.py` — the three models; `sim.py` and `sweep.sh` — the simulator and its sweep (host Python
  with numpy, `uv run --with numpy python sim.py inbred 63 3`).
- `run_conditions.sh` and `psps_parallel.sh` — the tomato conditions, run inside the dev container
  with `DEV_EXTRA_MOUNT` pointing at `benchmarks/tomato1`; they were written for a worktree at
  `../pop_var_caller-neighbour-hmm` and name it in their first lines. Steps: `downsample`, psps with
  `psps_parallel.sh`, `deep_psps`, `call_one <condition>`, `single`.
- `eval_tomato.py <calls dir> <condition> <deep names file> [K]` — the scorer; a condition named
  `single_<x>` scores one call per accession and takes the `SRR SRS` map.

The downsampled read files (0.6 GB), the map from read file to sample name and the per-file read
counts are kept, untracked, in the main checkout's `tmp/kept_from_worktrees/neighbour_evidence_2026-11-01/`.
tomato1's own repeat catalog under `benchmarks/tomato1/crams/` predates the 2026-10-30
reference-digest change; build a fresh one with `pvc repeat-catalog`.

---

## 11. The owner's decisions (2026-11-01)

After reading this report the owner ruled:

- **Build it, as pvc's own pass.** The 2026-09-13 decision keeps its shape — a second pass with the
  VCF as interface — and the pass is pvc's.
- **Where it sits:** a window of a few megabases held behind pvc's output, after the
  hidden-duplication filter when that is on; a spill file is not needed.
- **When it runs:** on by default when the cohort is large enough to supply templates and its median
  fitted inbreeding coefficient is high; off otherwise, saying so. The thresholds are set from
  measurements once it is built.
- **The model:** shared haplotypes from genotype templates, with local autozygosity as a state inside
  it rather than a module of its own. Read linkage is left.
- **Gains at one or a few samples are not expected**, and are not a goal.
- **Order of work:** `PL` and `GP`; a read-level cohort simulator with shared haplotypes and
  autozygous runs; a spec for the pass; then the implementation, with the prototype's Python as the
  reference its posteriors are checked against.
