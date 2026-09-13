# 06 — What the literature measures: phasing and imputation accuracy, and how imputation has been coupled to variant calling

Literature report, 2026-09-12. Companion to the other sub-reports under
`doc/devel/ng/research/phasing_imputation/`. The reader is a geneticist who is not a statistician.

## 0. Scope and how to read the numbers

This report collects **measured** accuracy of phasing and imputation as a function of the four things
the owner can control or is stuck with — sequencing depth, cohort size, reference-panel size and allele
frequency — and then asks how imputation has been wired into variant calling. Every number carries its
paper, the data it was measured on, the panel, the metric, and which direction is better. Own
reasoning is marked as such.

Three conventions:

- **Depth** is written as `0.5×`, `1×`, `4×` (mean reads a position). **N** is the number of samples in
  the cohort being phased or imputed. **Panel** is the number of *haplotypes* (two a diploid sample)
  unless the paper counts individuals, in which case it says so.
- **Direction.** For switch error rate and discordance, *lower* is better. For r², concordance, F1,
  N50 and IQS, *higher* is better. Each table header says which.
- **Verification.** I opened every paper cited (PubMed Central or bioRxiv full text) and took numbers
  from the results text, tables or figure captions. Where a number exists only as a plotted point I say
  "figure only" and do not quote a value. Where a fetch returned an approximate reading of a figure I
  discarded it. The list of what I could not verify is in §10.

## 1. The metrics, defined for a geneticist

**Switch error rate (SER).** Take one sample's heterozygous sites in genome order, and the true
haplotypes (from parents, or a trio child). Walk along the estimated haplotypes: each time the
estimate jumps from tracking the maternal chromosome to tracking the paternal one (or back), that is a
*switch*. SER = switches ÷ (number of consecutive pairs of correctly-called heterozygous sites). A
SER of 1% means one phase break every ~100 heterozygous sites. SHAPEIT5's own checker implements
exactly this: it counts, at each phased-and-validated heterozygote, whether the estimated haplotype
keeps or breaks the true relationship with the previous heterozygote, then reports
`n_phasing_errors × 100 / n_phased_hets`
(`tmp/phasing_research/repos/shapeit5/switch/src/models/haplotype_checker.cpp:41-87`). Two things
about SER that matter for reading the literature:

- *It is computed only at heterozygous sites*, so an inbred sample with few heterozygotes has few
  opportunities for error; the same SER means far fewer errors in an inbred than in an outbred sample.
- A single misplaced site counts as *two* switches (in, then back out). Some benchmarks separate these
  as **flip errors** and report "switch" for the long-range breaks only (the 2025 benchmark of Beagle
  5.4, SHAPEIT and Eagle defines a flip as "a single variant's phase is incorrect relative to its
  neighbours" [Beck, Kang and Zöllner 2025/2026, PMC13131814]). SHAPEIT5's checker also writes both
  counts (`haplotype_checker.cpp:108-124`). When a paper says SER = 0.1% at N = 400,000, that is mostly
  flips at rare variants, not long-range breaks.

**N50 of correctly phased segments (and AN50, QAN50).** Read-based phasers produce *blocks*, stretches
of consecutive heterozygotes joined by reads; between blocks the phase is unknown. N50 is the block
length such that half of all phased heterozygotes fall in blocks at least that long. *AN50* adjusts
each block by the fraction of its sites that are phased correctly; Choi et al. 2018 use *QAN50*,
"largest quality-adjusted haplotype block length such that 50% of all heterozygous sites are contained
in haplotype blocks of quality-adjusted length at least QAN50" (their Table 1) [Choi 2018]. A method
can have a very low SER and a very short N50: Illumina-read phasing of NA12878 in that paper has SER
0.096% but QAN50 1.0 kb, while population phasing with the HRC panel has SER 0.297% and QAN50 1.9 Mb
(Choi 2018, Table 2). SER and block length must be read together.

**Genotype concordance.** Fraction of genotypes (0/0, 0/1, 1/1) that agree with the truth set. Its
weakness for rare variants: at a variant with minor allele frequency 0.5%, calling every sample
0/0 is 99% concordant. Ramnarine et al. 2015 show the point formally — concordance "will always
produce a value greater than or equal to IQS" because it ignores agreement by chance, and for rare
variants "concordance and BEAGLE R2 show inflated assessments while IQS and r² reveal poor accuracy"
(their Fig. 2) [Ramnarine 2015]. Papers that report concordance on inbred lines (soybean 97.8% at
0.3×, wheat 97.5% at 0.07×) are reporting a metric that is easy to score high on.

**Non-reference concordance / non-reference discordance (NRC / NRD).** Concordance restricted to
genotypes where either the truth or the call carries the alternate allele. It removes the 0/0-versus-0/0
agreements that inflate plain concordance. GLIMPSE's concordance tool computes it as
`errors × 100 / (errors + correct non-ref genotypes)`
(`tmp/phasing_research/repos/GLIMPSE/concordance/src/containers/call_set_writing.cpp:46-55`).

**Imputation r² — per-variant and aggregate.** The squared Pearson correlation between imputed dosage
(expected count of the alternate allele, 0–2) and the true genotype. Per-variant r² at a rare variant is
a correlation over a handful of carriers and is unstable — Ramnarine et al. report that "a MAF bin
can have a wide range in accuracy values" and that "widely discrepant values for IQS and squared
correlation are attributable to rare and low frequency SNPs" (their Fig. 4). The standard fix is the
**aggregate r²**: pool all (sample, variant) pairs whose variant falls in one allele-frequency bin, and
compute one correlation over the pool. GLIMPSE's authors state it explicitly: "we pooled all validation
and imputed genotypes belonging to the same frequency bin together in order to compute a single squared
Pearson correlation value per bin", over 17 bins from ]0, 0.0002] to ]0.45, 0.5] (Rubinacci 2021,
Methods). Aggregate r² is the metric behind every "r² by MAF" curve in §3–§4. Its known bias: it is
dominated by whichever variants in the bin have most carriers, and it is measured against a
*truth set* that is itself a filtered high-coverage call set (GLIMPSE requires posterior > 0.9999 and
depth ≥ 8× at the truth genotype), so hard-to-call sites are under-represented.

**Why aggregate r² and not concordance is the field's yardstick.** r² has a direct meaning for
downstream power: a sample imputed at r² is worth r² of a perfectly genotyped sample in an association
test, so *effective sample size* = N × r² (Li 2011 Table 1; Pasaniuc 2012). Concordance has no such
reading.

**INFO score (IMPUTE / GLIMPSE / STITCH) and DR2 (Beagle).** Truth-free estimates of r², computed
from the imputed posteriors alone. INFO answers: how much of the variance a perfect genotype would have
at this allele frequency is actually present in the posteriors? GLIMPSE computes it per site as
`1 − (Σ(GP1 + 4·GP2) − Σ DS²) / (n_haps · p · (1 − p))`
(`tmp/phasing_research/repos/GLIMPSE/phase/src/io/genotype_writer.cpp:173-190`); STITCH's header
says its `INFO_SCORE` is "same as IMPUTE info measure I_A"
(`tmp/phasing_research/repos/STITCH/STITCH/R/writers.R:1091`). Beagle's DR2 is declared as the
"estimated squared correlation between estimated REF dose [P(RA) + 2·P(RR)] and true REF dose"
(`tmp/phasing_research/repos/beagle_src/src/vcf/VcfWriter.java:42-44`). Known biases: (i) both are
computed from the model's *own* posteriors, so a confidently wrong model scores high — König et al.
2024 found the software's rsq "comparable between accessible and inaccessible regions" of the human
genome while the true r² was lower in the inaccessible ones (their Supplementary Fig. 1b) [König 2024];
(ii) Beagle's DR2 tracks concordance, not r², at rare variants and is inflated there (Ramnarine 2015,
Fig. 2); (iii) the correlation between estimated and true r² for variants with MAF ≤ 0.5% was
0.78–0.83 in a study-specific-panel setting (Pistis 2015, reported in the search summary; I did not
verify the figure). Use INFO/DR2 to *rank* sites for filtering, not to report accuracy.

**IQS (imputation quality score).** Cohen's kappa between imputed and true genotypes: agreement
beyond chance. Introduced by Lin et al. 2010 to score rare SNPs where concordance is uninformative
[Lin 2010]. Used by the dog and wheatgrass studies below.

## 2. Phasing accuracy versus cohort size, no reference panel

All population phasers (Beagle, SHAPEIT, Eagle) get their signal from other samples in the same run:
the more samples, the longer the shared haplotype stretches they can copy from. The published curves
start at N = 500 and stop at N ≈ 500,000. Lower is better throughout.

| Paper | Data | Validation | N | Tool | SER | Where |
|---|---|---|---|---|---|---|
| Delaneau 2019 (SHAPEIT4) | UKB SNP array, chr20, 18,477 SNPs | 500 trio children | 500 … 400,000 (10 points) | SHAPEIT4, Beagle5, Eagle2, SHAPEIT3 | at N = 400,000: SHAPEIT4 (P = 4) 0.117%, Beagle5 0.125%, Eagle2 0.178%, SHAPEIT3 0.356%; values at N ≤ 200,000 are figure only | Fig. 2a, results text |
| Browning 2021 (Beagle 5.2) | UKB SNP array, 711,651 markers | 1,064 trio offspring | 5,000; 15,000; 50,000; 150,000; 485,301 | Beagle 5.2, SHAPEIT 4.2.1 | "very similar phase error … for all sample sizes"; values figure only | Fig. 2 |
| Browning 2021 (Beagle 5.2) | TOPMed WGS chr20, 7.2 M markers | 217 + 669 trio offspring | 5,000; 10,000; 20,000; 38,387 | Beagle 5.2, SHAPEIT 4.2.1 | "similar error rates"; figure only | Fig. 3 |
| Loh 2016 (Eagle2) | UKB SNP array, chr 1/5/10/15/20 | 70 trio children | 15,000; 30,000; 50,000; 100,000 (as reference) | Eagle2 vs SHAPEIT2 | Eagle2 5–16% lower than SHAPEIT2; gains larger at smaller N | Fig. 2b, Supp. Table 2 |
| Hofmeister 2023 (SHAPEIT5) | UKB WGS chr20 | 31 trios + 432 duos | 147,754 | SHAPEIT5 vs Beagle 5.4 | MAC 11–20: 4.36% vs 8.76%; singletons 35.1% (random = 50%) | Fig. 2a |
| Hofmeister 2023 (SHAPEIT5) | UKB WES | 719 trios + 3,014 duos | 452,644 | SHAPEIT5 vs Beagle 5.4 | MAC 11–20: 2.93% vs 5.18% | Fig. 2b |
| Hofmeister 2023 (SHAPEIT5) | UKB SNP array | trios | 5,000 … 480,000 (9 points) | SHAPEIT5 vs SHAPEIT4 | "significant differences … observed in datasets comprising at least 50,000 samples"; values figure only | Extended Data Fig. 1 |
| Oget-Ebrad 2022 (cattle) | Holstein WGS pedigree | 98 sequenced trios | 98 | Beagle 5.2 / SHAPEIT 4.1 / Beagle 4.1 / Eagle 2.4 | 0.0093% / 0.0192% / 0.140% / 0.249% (median per animal) | Table (scenario 1) |
| Oget-Ebrad 2022 (cattle) | same | same | 264 | same | 0.0027% / 0.0026% / 0.0183% / 0.0451% | Table (scenario 2) |
| Choi 2018 | NA12878 vs GIAB truth, array-density sites | 1 sample | 85 (one EUR sub-population as SHAPEIT2 panel) | SHAPEIT2 | 1.9% (CEU subgroup); EUR (347) 2.0%; all 1000GP (2,503) 1.03%; HRC (22,690) 0.30% | Table 2 |

Reading the table:

- **Between N = 500 and N = 500,000 the SER curves are monotone and shallow.** SHAPEIT4's text says
  error rates "substantially decrease as sample size increases" from N = 500 to 400,000, but the
  end-point at 400,000 is 0.12–0.18% for the modern tools, and SHAPEIT5 finds that SHAPEIT4 and
  SHAPEIT5 only separate at N ≥ 50,000. The large-N gains are concentrated at *rare* variants:
  SHAPEIT5's WGS SER at MAC 11–20 is 4.4% even with 147,754 samples, and 35% for singletons.
- **N = 50–100 has no published curve.** SHAPEIT4's Fig. 2a starts at 500; Browning 2021 at 5,000.
  The closest measurements are Choi 2018's 85-sample panel — 1.9% SER on a European sample at
  array-density sites, versus 1.0% with 2,503 and 0.30% with 22,690 — and the cattle pedigree, where
  98 animals phase at 0.009% because parents and offspring are in the run (identical-by-descent
  segments are chromosome-length). The cattle number is a pedigree number, not a small-cohort number.
- **Own reasoning on N = 50–100 for ng.** The 85-sample European figure (1.9% at ~31,853 array sites
  on chr1) is the right order of magnitude for an unrelated outbred cohort at array density; at
  sequence density the SER per heterozygote falls (more sites between switches) but the number of
  switches per Mb does not. For a highly inbred cohort (tomato accessions, median F ≈ 0.78, `ng_facts.md`)
  the heterozygous sites are sparse and clustered in residually heterozygous blocks, and no paper in
  this survey measures SER there.
- **No plant or non-pedigree small-cohort phasing benchmark surfaced.** Searches for switch error in
  rice, maize, Arabidopsis or tomato cohorts returned nothing quotable. That is a gap, not evidence of
  good or bad performance.

## 3. Imputation accuracy versus coverage, with a reference panel

The low-coverage imputers (Beagle 4.1 `gl=`, GLIMPSE, QUILT, loimpute) take genotype likelihoods or
reads and copy haplotypes from a phased panel. The signal is the panel; depth only decides how well the
sample's reads pick which panel haplotypes to copy. Higher r² is better throughout.

### 3.1 Human, large panels

| Paper | Target | Panel | Depth | Result (aggregate r² unless stated) | Where |
|---|---|---|---|---|---|
| Rubinacci 2021 (GLIMPSE1) | 503 Europeans, 1000G | HRC, 54,330 haplotypes | 0.1× / 0.3× / 1× / 8× | 0.1×: r² > 0.7 at MAF > 5%; 0.3×: r² > 0.9 at MAF > 5%; 1×: r² = 0.8 at MAF 0.1%; 8×: r² > 0.95 at MAF 0.1% | preprint v1 results text, Fig. 2a |
| Rubinacci 2021 (GLIMPSE1) | same, chr1 | HRC | 1× | at MAF 0.1%, GLIMPSE ≈ 0.2 r² above Beagle 4.1 `gl` (the second-best); common variants: GLIMPSE ≈ Beagle 4.1 > GeneImp > STITCH | Fig. 2c |
| Rubinacci 2021 (GLIMPSE1) | NA12878 vs GIAB | HRC | 0.1–8× | phasing SER of the imputed haplotypes 0.72–1.23% | Fig. 2b |
| Rubinacci 2021 (GLIMPSE1) | 61 African Americans (ASW) | HRC | 0.5–1× | 1× beats Omni2.5 array + imputation at MAF < 1%; 0.5× matches it at common variants; advantage larger in ASW than in Europeans (array ascertainment bias) | Fig. 3a–b |
| Rubinacci 2023 (GLIMPSE2) | UKB British samples | UKB, 150,119 genomes (280,238 haplotypes) | 0.1× / 1× | MAF 0.01% bin: GLIMPSE2 0.892 / 0.927; GLIMPSE1 0.561 / 0.725; QUILT 1.0.4 0.728 / 0.925 | results text (Fig. 1) |
| Rubinacci 2023 (GLIMPSE2) | 10,000 UKB British, chr1 | UKB | 0.5× | > 0.1 r² gain over the Axiom array for MAF < 0.01% | Fig. 1c |
| Rubinacci 2023 (GLIMPSE2) | 276 SGDP samples, 129 populations | UKB vs 1000G | 1× | UKB panel cuts NRD by > 67% for Northern Europeans relative to 1000G; other populations figure only | Extended Data Fig. 2 |
| Davies 2021 (QUILT) | NA12878 | HRC | 0.5× / 1× | rare SNPs at 0.5×: QUILT 0.678 (Illumina), GLIMPSE 0.672; common 0.975 / 0.974; rare at 1×: 0.754 (Illumina), 0.741 (ONT), 0.776 (haplotagging); array ≈ 0.694 | results text |
| Davies 2021 (QUILT) | 1000G CEU vs CHB | HRC | 0.25× / 1× | rare r²: CEU 0.63–0.66, CHB 0.581 (GSA array 0.485); CHB at 1×: 0.768 | results text |
| Davies 2025 (QUILT2) | CEU | 1000G (~5,008 haplotypes), HRC (54,330), UKB 200k (400,022) | 0.1–2× | QUILT2 > GLIMPSE2 at ≤ 0.5× on every panel; GLIMPSE2 > QUILT2 at ≥ 2× on rare variants; "very rare" bin = MAF 0.01–0.02%, "rare" = 0.1–0.2%, "common" = 10–20%; per-cell values in Supp. Tables 1–2 (not seen) | Fig. 2a |
| Davies 2025 (QUILT2) | NA12878 ONT long reads | UKB | 1× | common r² 0.937 (QUILT2) vs 0.695 (GLIMPSE2) | results text |
| Wasik 2021 (loimpute) | 79 Europeans, Cambridge UK | 1000G phase 3; HRC | 1× / 0.8× / 0.6× / 0.4× | 1×: mean r² 0.96 at common variants vs PMRA array + minimac2 0.90; over all frequencies 0.93 vs 0.85; 0.4×: 0.91 at common; HRC vs 1000G "marginal", −0.036 in the lowest MAF bin | Fig. 2, Table 1 |
| Pasaniuc 2012 | 909 exomes, off-target reads | 1000G, 762 European haplotypes | 0.24× | mean r² 0.71 (s.d. 0.15); 0.5× on 84 samples in 10 regions: 0.82 | results text |
| Gilly 2019 (HELIC MANOLIS) | 990 at 1× (+249 at 4×) | 10,244 haplotypes: 1000G phase 1 + 249 MANOLIS 4× + 3,781 UK10K | 1× | minor-allele concordance: 97% at MAF > 5%; 73% at 1–5%; 55% at MAF < 1% | Fig. 2 |
| Li 2011 (Thunder) | 1000G pilot CEU, 43–60 per population | none — the cohort itself, ~60 samples | ~4× | dosage r² 84.9% at MAF 1–2%, 88% at 2–5%, ≈ 95% at MAF > 5%; genotypic concordance > 98% | Fig. 6 |

What the human numbers say, in order of practical weight:

1. **Depth buys accuracy fast between 0.1× and 1×, slowly after.** With a 54,000-haplotype panel,
   common variants (MAF > 5%) pass r² 0.9 at 0.3× and rare ones (MAF 0.1%) reach 0.8 at 1× (GLIMPSE1).
   Between 1× and 8× the rare-variant r² moves from 0.8 to > 0.95.
2. **Panel size is the lever for rare variants.** At MAF 0.01% and 1×, moving from GLIMPSE1 to GLIMPSE2
   *on the same 280,238-haplotype panel* moved r² from 0.725 to 0.927; the UKB panel cut non-reference
   discordance by two-thirds against 1000G for Northern Europeans. Wasik's 79 Europeans saw only a
   marginal gain from HRC over 1000G phase 3 at 1×, because at 1× and MAF > 1% the 5,000-haplotype panel
   already saturates.
3. **Population match matters as much as size.** The same HRC panel gives rare-variant r² 0.63–0.66 for
   CEU and 0.58 for CHB at 0.25× (QUILT).
4. **Arrays are beaten at 0.4–1×** on every comparison above (loimpute at 0.4×, GLIMPSE at 0.5×, QUILT at
   1× for rare SNPs).

### 3.2 Non-human species, small panels

Panels here are tens to a few hundred individuals, built by the study itself. Higher is better.

| Paper | Species, target | Panel | Depth | Tool | Result | Where |
|---|---|---|---|---|---|---|
| Happ 2019 | soybean, 114 inbred lines | 99 lines at 17.1× (from 106), 10.8 M homozygous SNPs; heterozygous calls discarded | 0.1–1.0× | Beagle 4.1 | 97.8% concordance at 0.3× on 5% masked calls | results text |
| Lloret-Villas 2023 | Brown Swiss cattle | within-breed 30 / 75 / 150 animals; multibreed 150 | 0.01–4× | GLIMPSE 1.1.1 | F1 > 0.9 at 0.25× with the 150-animal within-breed panel; within-breed 150 beat every multibreed panel at every depth and MAF; panels without the target breed "substantially reduced"; per-cell r² figure only | Abstract, Fig. 4 |
| Teng 2022 | Holstein cattle | 200 / 400 / 600 / 800 / 1,059 animals | several | Beagle 4.1 `gl`, GeneImp, GLIMPSE, QUILT, Reveel, STITCH | all but Reveel > 0.9 accuracy in most scenarios; GLIMPSE, QUILT, STITCH best; recommends STITCH then Beagle when no panel exists | Abstract (full text paywalled) |
| Koorevaar 2025 | octoploid strawberry, 66 test samples from 3 populations | 70 / 127 / 197 / 587 individuals at ≥ 15× | 1× | SHAPEIT5 + GLIMPSE2 | concordance 0.87–0.97 with 70; 0.94–0.98 with 127 or 587; homozygous 0.99 vs heterozygous 0.88 (median); conclusion "≈ 70 genetically representative samples at ≥ 25× … sufficient" | Fig. 7, Table 3 |
| Watowich 2025 | rhesus macaque; gelada | 741 rhesus; 68 gelada at ~11× | 0.1 / 0.5 / 1 / 3 / 10× | loimpute | at 0.5×, median r² 0.92 (rhesus, 741 panel) vs 0.86 (gelada, 68 panel) for MAF 10–50% | Fig. 2c–d |
| Buckley 2022 | dogs, 1× | 676 dogs, 91 breeds | 1× | loimpute | IQS plateau 0.91 unfiltered / 0.95 filtered at MAF > 0.1; heterozygote concordance 90.8% → 96.4% after filtering; breeds in the panel score significantly higher than breeds absent (Table S8) | results text |
| Gundappa 2025 | Atlantic salmon | 365 wild fish | 1–4× | GLIMPSE 1.1.1 | SNVs at 1×: PPV 0.98, recall 0.62 at GP > 0.9; deletions at 1× (merged strategy): 84% recall, 87% accuracy; out-of-panel farmed fish similar (PPV 0.86–0.87) | results text |
| Topaloudis 2026 (barn owl) | 2,800 owls at ~2× (0.2–4.15×) | 50 … 502 individuals | ~2× | GLIMPSE 1.1.1, QUILT 2.0.4 | GLIMPSE r² 0.978 (s.d. 0.008) with the full panel; "all reference panel sizes … above 0.9, … above 200 samples returning more than 0.95"; NRC > 0.95 at MAF ≥ 1–2%, "sharp decrease" below 0.2% | results text, Fig. 3 |
| Vi 2025 | hihi (NZ bird), 283 targets | 30 individuals at 20–30× | 1× | GLIMPSE2, QUILT2, STITCH, GeneImp, Beagle 5.4 | per-variant F1: QUILT2 0.99, GLIMPSE2 0.97, GeneImp 0.95, STITCH 0.93, Beagle 5.4 0.85 | Table 1 |

Three rules of thumb these give for the owner's "many haplotypes" case:

- **A panel of 70–150 deep individuals gets a low-diversity population (a breed, a pedigreed wild
  population, a crop germplasm collection) to r² ≈ 0.9 at 0.5–1×** (strawberry 70; cattle 150; owl 50
  → r² > 0.9, 200 → > 0.95). Every one of these populations has long shared haplotypes; none is a
  diverse outbred human-like population.
- **The panel must contain the target's population**: cattle within-breed 150 beats multibreed 150 at
  every depth; dogs of breeds absent from the panel impute worse; salmon out-of-panel fish lose ~0.1 in
  PPV.
- **Concordance on inbred lines is not r².** Soybean's 97.8% and wheat's 97.5% are on homozygous-only
  calls; the heterozygous-genotype concordance in strawberry is 0.88 against 0.99 homozygous.

## 4. No external panel: the cohort as its own panel

### 4.1 STITCH — reference-free imputation from reads

STITCH assumes the cohort descends from K ancestral haplotypes a few hundred generations back and fits
those haplotypes from all the low-coverage reads at once, so the "panel" is estimated rather than
given. Results, higher r² better:

| Paper | Population | N | Depth | K | Result | Where |
|---|---|---|---|---|---|---|
| Davies 2016 | CFW outbred mice | 2,073 | 0.15× | 4 | r² 0.972 vs array, 0.948 vs 10× WGS unfiltered; 0.981 / 0.974 after INFO > 0.4 and HWE filtering (75–81% of SNPs kept) | results text |
| Davies 2016 | same, downsampled | 100–2,073 | 0.015–0.15× | 4 | "sample size above 500 has little impact on performance" at 0.15×; at 0.06× with all 2,073 "only marginally poorer than 0.15×"; N matters more as depth falls; per-cell values figure only | Fig. 4a |
| Davies 2016 | CONVERGE Han Chinese, chr20 first 10 Mb | 11,670 | 1.7× (0.3–1.7× tested) | 40 | r² 0.920 vs array, 0.949 vs 10× WGS unfiltered; Beagle 4 with the 1000G ASN panel 0.943 vs array "at 7.3× the run time"; sample size "less influence" at 0.3–1.7×, depth "consistently improved" | Fig. 3c–d, Fig. 4b |
| Topaloudis 2026 | pedigreed wild owls | 50 … 2,800 | ~2× | 30 (nGen 1,350) | r² 0.970 (s.d. 0.012) with all 2,800; "diminishing returns … more than 500 samples, which achieved … approximately 0.95"; at N = 50–100 the per-sample r² tracked depth almost perfectly (r > 0.9), i.e. "imputation was rather ineffective" | Fig. 4a, results text |
| Crain 2026 | 1,709 wheat lines | 1,709 | 0.07× | 15, `diploid-inbred` | 97.5% and 97.7% of 67,000 loci correct in two 30×-validated lines; 121,437 markers retained at INFO > 0.8, MAF > 0.05 out of 14 M | results text |
| Sthapit 2025 | 9,780 intermediate wheatgrass genets | 9,780 | 0.05× (0.01–0.10× tested) | 8 (4–28 tested) | concordance > 0.97, IQS > 0.95, r² ≈ 0.90 against 46 genets at 17× (INFO > 0.8) | Fig. 3, Fig. 5 |
| Pierotti 2024 | F2 of 8 inbred lines, 10 crosses | 2,177 | 0.25–1.4× | 16 (nGen 2) | mean r² 0.996 full data; 0.981 at 0.5×; 0.977 in the MAF 0–0.05 bin; a single cross of 474 fish: > 0.964 in most bins | Fig. 3 |
| Teng 2022 | Holstein | ≥ 400 | ≥ 1× | — | "accuracy plateaued when sample size exceeded 400 and depth surpassed 1×" (search-engine summary of the paywalled text; r² 0.98 ± 0.002) | not verified |

The STITCH answer to "how many samples do I need with no panel" is consistent across species: **about
500 low-coverage samples for r² ≈ 0.95 in a population with few founders (K ≤ 30); below 100 the method
does not add information beyond what the reads already say** (owl, mice). Depth and N trade off: 2,073
mice at 0.06× ≈ 0.15×; 100 mice at 0.15× is in the figure but the value is not in the text.

### 4.2 Beagle in genotype-likelihood mode with no panel

Beagle 4.0/4.1 accepts `gl=` (a VCF with GL or PL; any GT is ignored,
`tmp/phasing_research/low_coverage_imputation/beagle_4.1_09Feb16.txt` lines 98–100) and, without
`ref=`, builds haplotypes from the cohort alone. Beagle 5.0 and later dropped this: "Beagle 5.0 does not
infer genotypes from genotype likelihood data, but versions 4.0 or 4.1 can be used for this purpose"
(Beagle 5.0 release page). Measured use without a panel:

- The 1000 Genomes pilot (179 samples, 2–6×, mean 3.56×) and Phase 1 (1,092 samples, mean 5.1×) called
  genotypes by LD-based refinement of likelihoods across the cohort (Beagle, MaCH/Thunder, IMPUTE2,
  SNPTools) — see §6.
- Gilly 2019 ran Beagle 4 twice on 1,239 MANOLIS samples: "a first round of imputation-based genotype
  refinement" then a "reference-free imputation" round, before panel imputation.
- In the head-to-head at 1× with a 54,330-haplotype panel Beagle 4.1 was second only to GLIMPSE at
  common variants and ≈ 0.2 r² below it at MAF 0.1% (Rubinacci 2021, Fig. 2c). Beagle 5.4 (hard
  genotypes, no likelihoods) was the worst tool in every low-coverage scenario of Vi 2025 (F1 0.85 on
  hihi vs 0.99 for QUILT2, Table 1) — it is not a low-coverage tool.

### 4.3 Internal reference panels: sequence a subset deep, impute the rest

This is the owner's second case exactly. What the literature has measured:

| Paper | Design | Finding | Where |
|---|---|---|---|
| Li 2011 | simulation, equal total effort ≈ 12,000× | 3,000 samples at 4× vs 400 at 30×: same SNP discovery (≈ 100% at MAF > 0.5%) and concordance (> 99.67% at MAF > 1%), but effective sample size N·r² 4.8× (MAF 0.1–0.2%) to 7.2× (MAF 2–5%) larger for the shallow design; power for a MAF 0.5% risk allele 75.6% (3,000 at 4×) vs 6.4% (400 at 30×) vs 41.8% (1,000 at 12× + imputation) | Table 1, Table 2 |
| Li 2011 | simulation | "at depth 4×, > 98% genotypic concordance is achieved across the examined MAF spectrum with as few as 60 sequenced individuals" — the cohort itself is the panel | Fig. 2 |
| Pasaniuc 2012 | cost model, $300,000 | 0.1× on 6,800 samples gives an effective N > 4,600, above what arrays buy | results text |
| Pistis 2015 | 2,120 Sardinians at 4.16× (695 families); panel = 1,488 of them | array + SardSeq panel: mean r² 0.94 at MAF ≥ 1%, 0.91 at 0.5–1%; adding 1,000 Sardinians to 1000G reaches the same accuracy as the full SardSeq panel in every bin; 500 does not | Fig. 2, Supp. Table S11 |
| Gilly 2019 | 249 of 1,239 at 4×, rest at 1×; the 249 join a 10,244-haplotype panel | at 1×: 97% minor-allele concordance at MAF > 5%, 73% at 1–5%, 55% at < 1% | Fig. 2 |
| Koorevaar 2025 | strawberry: 70 / 127 / 197 / 587 deep | 70 gives 0.87–0.97 concordance; 127 gives 0.94–0.98; "≈ 70 genetically representative samples at ≥ 25× … sufficient" | Fig. 7 |
| Watowich 2025 | gelada: 68 at ~11× | median r² 0.86 at 0.5× for MAF 10–50% (vs 0.92 with 741 rhesus) | Fig. 2c–d |
| Barn owl 2026 | 50 … 502 deep | > 0.9 with 50; > 0.95 above 200; 0.978 with 502 | Fig. 3 |
| Lloret-Villas 2023 | cattle: 30 / 75 / 150 within-breed | F1 > 0.9 at 0.25× with 150; 150 within-breed > any 150 multibreed | Fig. 4 |
| Vi 2025 | simulated: 20 deep, 80 at 1× | GLIMPSE2 and QUILT2 "outperformed all the other tools in all the scenarios"; every tool's accuracy "improved proportionally to the increase in genetic relatedness" (except QUILT2, already saturated); in the two high-inbreeding, high-relatedness populations all tools > 0.97 F1 | Fig. 2 |

**How many deep samples for a given r² in the shallow ones — the honest summary.** No paper varies
the deep-subset size in a *diverse* population and reports r² at each size; the size series that exist
(owl 50–502, strawberry 70–587, cattle 30–150, Sardinia 500–1,488) are all in populations with long
shared haplotypes, and in those **50–150 deep individuals reach r² 0.9 and 200–500 reach 0.95 at 0.5–2×
target depth**. For an outbred diverse population the only anchor is human: 1,000 population-specific
genomes added to 1000G matched a 1,488-genome panel (Sardinia), and 249 population samples at 4× were
enough to lift a 10,000-haplotype panel for an isolate (MANOLIS). Own reasoning: the number of deep
samples needed scales with the number of distinct haplotypes segregating at the frequency you care about
— a panel of H haplotypes cannot impute an allele carried by fewer than ~2–3 of them (Pistis: rare bins
are where the study-specific panel wins; GLIMPSE2: MAF 0.01% needs 280,000 haplotypes).

## 5. Segregant populations: F2, RIL, backcross, MAGIC at 0.02–1×

Here the model is a two- or few-founder hidden Markov model: each individual is a mosaic of parental
haplotypes with a handful of breakpoints per chromosome, so a marker's genotype is determined by its
neighbours over megabases. Accuracy is limited by breakpoint placement and by allelic read bias, not by
haplotype diversity. Higher is better.

| Paper | Tool | Population | N | Depth / density | Result | Where |
|---|---|---|---|---|---|---|
| Fragoso 2016 | LB-Impute | simulated F2 and F1BC1 | 20 replicates × 40 datasets | 0.1–4× (1,000–40,000 reads a sample) | LB-Impute > 99% correct genotypes at every depth 0.1–2.5×; FSFHap ≥ 99% only at ≥ 0.4× (F2) or ≥ 0.8× (F1BC1) and could not run below 0.4×; Beagle 37.07%, Mendel Impute 50.69% at 0.1× (F1BC1); breakpoints 76.9–90.1% correct (LB-Impute) vs 12.9–80.7% (FSFHap) | Fig. 2 |
| Fragoso 2016 | LB-Impute | maize B73 × Country Gentleman F2 (89–90); IBM RILs (275) | 11,219 / 127,144 / 14,493 markers | GBS | concordance 94.62% ± 3.02, 91.44% ± 6.73, 96.98% ± 5.14 | Fig. 3 |
| Rowan 2015 | TIGER | Arabidopsis Col-0 × Ws-2 F2 | 110 wild-type + 106 recq4a | median 0.4× (mean 0.6×), 261,795 SNPs | 90% of crossovers placed within 2 kb at every depth tested; median resolution 1.986 kb at 0.1×; background genotype error ~0.1%, 2.4% near centromeres/telomeres at 0.1×; 8 of 11 predicted crossovers confirmed by Sanger | Table 2, Table 3 |
| Zheng 2018 | magicImpute | maize F2 (87), maize AI-RIL (275), rice 8-way MAGIC (178), apple outbred CP (87) | 13,912–127,059 markers | offspring reads thinned by 2^−i | accuracy > 0.98 at full depth for maize and rice; accuracy holds down to "break points" of 0.053–0.21× before collapsing, lower than Beagle 4.1's and LB-Impute's; slightly worse than mpimpute only at high depth | results text |
| Furuta 2023 | GBScleanR | simulated F2 / F1 / 8-way RIL (1,000 each); real rice F2 (O. sativa × O. longistaminata) | 814 real | 0.1–20× simulated; real 0.85× over 5,035 SNPs | real data: GBScleanR 10.5–12.5 percentage points more concordant than magicImpute on chr 7, 8, 10; 1.7 points on chr 5; LB-Impute inferred 948 double crossovers vs 14 for GBScleanR; at 0.1× LB-Impute had the most correct calls but also the most miscalls | Fig. 1, Fig. 3, Table 2 |
| Pierotti 2024 | STITCH | medaka F2 of 8 inbred lines | 2,177 (474 in one cross) | 0.25–1.4×, 3.2 M SNPs | r² 0.996 full; 0.981 at 0.5×; 0.964 in one 474-fish cross | Fig. 3 |
| Sthapit 2025 | STITCH | wheatgrass breeding population (not a cross) | 9,780 | 0.05× | concordance > 0.97, r² ≈ 0.90 | Fig. 5 |

Two things carry over to ng:

- **At 0.1–0.5× a two-founder HMM gets > 99% genotypes right in simulation and 91–97% on real GBS
  maize**; the real-data gap is allelic read bias and mis-mapping, which GBScleanR models explicitly
  (marker-specific allele bias `wm` and mis-mapping rates `eref`, `ealt`) and which cost LB-Impute 948
  spurious double crossovers on one rice F2.
- **Breakpoint resolution is set by marker density and depth, not by the model**: TIGER places 90% of
  crossovers within 2 kb at 0.1–0.6× with 262,000 markers on a 120-Mb genome.

## 6. Imputation coupled to variant calling

### 6.1 The standard pipeline: per-sample likelihoods → LD-based genotype refinement

The pattern set by the 1000 Genomes Project in 2010 and still used is: a caller produces per-sample,
per-genotype likelihoods at candidate sites; an LD model across samples turns them into genotypes and
haplotypes.

- **1000G pilot (2010).** 179 samples at mean 3.56× (Table 1). "Statistically phased genotypes were
  derived by using LD structure in addition to sequence information at each site, in part guided by the
  HapMap 3 phased haplotypes." Overall genotype error 1–3% (Fig. 2c); heterozygote accuracy ≈ 90% for
  the lowest-frequency variants, > 95% at intermediate frequencies, 70–80% at the highest frequencies
  (results text); LD-based calling at 4× reached the accuracy of 15× single-sample exome calling
  (Fig. 2d) and called "nearly 15% more sites … with only a modest increase in error rate".
- **1000G Phase 1 (2012).** 1,092 samples at mean 5.1×. "By integrating LD information, genotypes from
  low-coverage data are as accurate as those from high depth exome data for SNPs with frequency > 1%";
  heterozygous-site accuracy > 99% at common SNPs, 95% at 0.5% frequency (Fig. 1b); no gain from LD at
  frequency ≤ 0.1%.
- **Li, Sidore et al. 2011 (Thunder/MaCH).** The HMM that did the refinement; on the pilot CEU at ~4×
  it delivers dosage r² 0.85 (MAF 1–2%) to 0.95 (MAF > 5%) from 60 samples with no external panel
  (Fig. 6). The paper does not tabulate before-versus-after-LD error; the before number quoted in
  secondary sources ("~20% heterozygous error at 4× without imputation") is not in this paper.
- **SNPTools (Wang 2013).** A complete pipeline — effective base depth, site discovery, a BAM-specific
  binomial mixture for likelihoods, then a "constrained Li-Stephens" imputation — on the 1,092 Phase 1
  samples: genotype discordance vs array 0.55% overall, 28.2% at heterozygotes, non-reference error
  1.24%; "superior switch accuracy … compared to Beagle, particularly at distances > 40 kb" (results
  text).
- **Today.** bcftools/GATK/freebayes emit PL/GL; Beagle 4.1 `gl=` or GLIMPSE consume them. GLIMPSE's
  tutorial keeps only bi-allelic SNPs from the panel and notes that "GLIMPSE can impute any type of
  variants as soon it is bi-allelic and has GLs being properly defined"
  (`tmp/phasing_research/repos/GLIMPSE/docs/docs/tutorials/getting_started.md:83,96`). By default
  GLIMPSE2 *discards* the caller's indel likelihoods: `--use-gl-indels` is an "expert setting" because
  "genotype likelihoods from low-coverage data are often miscalibrated, potentially affecting
  neighbouring variants", and indels are imputed "into the haplotype scaffold (assuming flat
  likelihoods)" (`GLIMPSE/docs/docs/documentation/phase.md:60,89`).
- **Illumina DRAGEN + GLIMPSE (IRPv1, 2,489-sample panel).** Marketed as "boosting variant calling":
  at 1× on three GIAB samples, "~90% reduction in the false call rate" mostly from fewer false negatives;
  no evaluation at 30× — the coupling exists only for low-pass data.

### 6.2 GATK CalculateGenotypePosteriors — a population prior with no LD

GATK's Genotype Refinement workflow multiplies the PLs by a prior from (i) a supporting population
callset's allele frequencies, (ii) the pedigree when both parents are present (a Mendelian-violation
penalty of 10⁻⁶), and/or (iii) the callset's own allele counts "if the input VCF contains at least 10
samples"; it then flags posterior GQ < 20. The documentation contains no mention of linkage
disequilibrium or haplotypes, and no accuracy numbers (GATK docs, Genotype Refinement workflow). This
is the same site-wise population prior ng already fits; it is not imputation.

### 6.3 Callers that carry LD inside

- **STITCH and QUILT read BAMs directly** (`STITCH/README.md:101`, `QUILT/README.md:90`) and emit
  genotypes, dosages and posteriors; they are callers of a kind, restricted to a supplied site list
  (`posfile`, `STITCH/README.md:192`) of bi-allelic SNPs.
- **Reveel** (Huang 2016) genotypes a cohort from low-coverage likelihoods with "a novel technique for
  leveraging LD that deviates from previous Markov-based models"; the abstract claims higher accuracy in
  low-frequency allele discovery than Markov methods but gives no number, and it was the one tool that
  failed to exceed 0.9 accuracy in the Holstein comparison (Teng 2022).
- **polyRAD** (Clark 2019) calls genotypes from allelic read counts in diploids and polyploids with
  priors from population structure or from linkage to neighbouring markers. Error reductions versus
  GATK's naive calls: 14.6% (diploid diversity panel), 23.5% (simulated tetraploid), 31.6% (diploid F1
  with linkage prior), 48.0% (tetraploid potato F1 with linkage prior); the LD prior helped in
  self-pollinating or structured panels (soybean, apple) and not in outcrossing *Miscanthus* (Figs. 2–5,
  results text).
- **GBScleanR, LB-Impute, magicImpute, TIGER** work from read counts or likelihoods in crosses (§5).
- **ANGSD/ngsTools** estimate per-site likelihoods and population parameters and do not use LD (not
  re-verified here; stated from the tools' scope).

### 6.4 Read-pair linkage inside an LD caller — the design closest to ng's chain ids

**HapSeq** (Zhi, Wu, Liu and Zhang, Bioinformatics 2012) extends Thunder's HMM with an extra emission
term for "jumping reads" — reads or pairs spanning two adjacent polymorphic sites — so the model sees
which alleles co-occur on a molecule. On the 1000G pilot chr20 the genotype discordance fell from 2.24%
to 1.97% in 47 CEU (12% relative) and from 2.76% to 2.50% in 52 YRI (9% relative) (Table 4); in
simulation at 4× from 0.86% to 0.60% (30% relative), and 43% relative at sites covered by jumping reads
from both directions (Table 1). This is the only published measurement of what read-backed linkage adds
*on top of* LD-based genotype refinement at low coverage, and the answer is a tenth to a third of the
remaining error.

### 6.5 Does anyone "call, impute, then re-call with an LD-informed prior"?

Yes, in three forms, and the gain is documented in each:

1. **The 1000 Genomes pipelines** are exactly this: per-sample likelihoods, LD refinement across the
   cohort, and the refined genotype *is* the call. Gain: 4× low-coverage reaching 15× exome accuracy
   (pilot Fig. 2d); > 99% heterozygote accuracy at common SNPs from 5× (Phase 1).
2. **Gilly 2019's two-round Beagle** on the cohort alone before panel imputation.
3. **Two-pass STITCH/QUILT**: QUILT2 first imputes on common SNPs to choose panel haplotypes, then
   re-loads all reads at all sites for the final Gibbs sampling ("~3× faster and ~4× less RAM" than one
   pass; the accuracy is the same, the pass structure is for cost).

Nobody re-runs a *variant discovery* step with the imputed haplotypes as prior; the second pass always
re-genotypes at fixed sites. HapSeq is the one design that folds read linkage into that second pass.

## 7. Where imputation hurts

**Rare variants in small panels.** The wall is set by carriers in the panel. GLIMPSE2 at MAF 0.01% needs
280,238 haplotypes to reach r² 0.93 at 1× — with GLIMPSE1 on the same panel it was 0.73. 1000G Phase 1
saw no LD gain at ≤ 0.1%. Barn owls (Topaloudis 2026): NRC "sharp decrease" below MAF 0.2% (GLIMPSE), below 1% (STITCH).
Biagini et al. 2025 (Genome Research; HRC, 27,165 individuals) judged MAF < 1% "too low for
downstream analyses" at 0.1–0.3×. MANOLIS at 1× with a 10,244-haplotype panel: 55% minor-allele
concordance at MAF < 1%.

**Population mismatch.** Same panel, different target: rare-variant r² 0.63–0.66 (CEU) vs 0.58 (CHB) at
0.25× (QUILT); cattle within-breed 150 beats multibreed 150 everywhere, and panels without the breed
lose badly; dogs of breeds outside the 676-dog panel impute significantly worse; STR imputation
concordance 97.0% (European) vs 90.6% (African) with an 83%-European panel (Saini 2018); Northern
Europeans gain a two-thirds NRD cut from a matched 150,119-genome panel over 1000G (GLIMPSE2).

**Variant classes the tools refuse or degrade.**

| Tool | Restriction | Evidence |
|---|---|---|
| Eagle2 | multi-allelic site is a hard error ("Either drop or split (bcftools norm -m)") | `Eagle/src/GenoData.cpp:730-733` |
| SHAPEIT5 phase_common | multi-allelic sites removed from the panel at scan time | `shapeit5/phase_common/src/io/genotype_reader/genotype_reader_scaning.cpp:121` |
| GLIMPSE2 | bi-allelic only; indel likelihoods replaced by flat likelihoods unless `--use-gl-indels` | `GLIMPSE/docs/docs/tutorials/getting_started.md:96`; `GLIMPSE/docs/docs/documentation/phase.md:60,89` |
| STITCH / QUILT | sites come from a `posfile` of ref/alt pairs — bi-allelic SNPs | `STITCH/README.md:101,192` |
| Beagle 5.x | phases and imputes multi-allelic markers but takes no genotype likelihoods (4.1 `gl=` does) | Beagle 5.0 release page; `beagle_4.1_09Feb16.txt:98-100` |

Measured accuracy on non-SNP classes: **STRs** imputed from SNP haplotypes with a 1,916-sample
HipSTR-genotyped panel reach 96.7% concordance overall but length r² 0.906 / allelic r² 0.861, ≈ 70%
concordance at the most polymorphic (CODIS) loci (Saini 2018, Table 1, Fig. 2). **Deletions** from 1×
salmon with GLIMPSE + SVtyper: 84% recall, 87% accuracy, versus SNV PPV 0.98 at recall 0.62 (Gundappa
2025). No paper in this survey measures indel imputation r² separately from SNPs; GLIMPSE's authors
distrust low-coverage indel likelihoods enough to zero them by default.

**Regions of low marker density.** König 2024: 7.6–28% of the human autosomes are "inaccessible"
depending on the mask, only 1.4–20.6% of panel variants fall there, true r² is lower there while the
software's own rsq is not (Supplementary Fig. 1b) — the failure is invisible to INFO/DR2.

**Inbred and structured cohorts.** Evidence is thin and points two ways. Relatedness helps every tool
(Vi 2025: all methods > 0.97 F1 in the high-inbreeding, high-relatedness simulations; owls' high
relatedness credited for rare-allele accuracy). But the inbred-line studies all *discard heterozygous
calls* (soybean: heterozygotes filtered before imputation; wheat: STITCH `diploid-inbred`, "no
heterozygosity data calculated") — STITCH's `diploid-inbred` method "assumes all samples are inbred and
invokes an internal haploid mathematical model but outputs diploid genotypes"
(`STITCH/STITCH/R/functions.R:14`). For a cohort at median F ≈ 0.78 with real residual heterozygosity
(the tomato accessions), no paper measures how a diploid HMM behaves; own reasoning: a diploid model
whose transition prior expects two independent haplotypes will spend its evidence explaining
homozygous stretches as "two identical haplotypes" and will be poorly calibrated at the residual
heterozygous blocks, which are exactly the loci of interest.

**Polyploids.** GLIMPSE, QUILT, STITCH, SHAPEIT and Eagle are diploid (GLIMPSE has haploid mode). Octoploid
strawberry was imputed as a diploid after homoeolog filtering (concordance 0.88 at heterozygotes vs 0.99
homozygous; A subgenomes better, Koorevaar Fig. 8). polyRAD is the one tool calling dosage in
autopolyploids with an LD prior (potato F1: 48% fewer errors than GATK). Hexaploid wheatgrass was run
through STITCH as diploid at 0.05× and still reached r² ≈ 0.90 on filtered sites.

## 8. Design implications for ng (measured facts only)

1. **Per-sample genotype likelihoods are the input every low-coverage imputer wants; hard calls are
   not.** Beagle 5.4 on hard genotypes was the worst tool in every 1× scenario (F1 0.85 vs 0.99, Vi 2025
   Table 1); Beagle dropped `gl=` after 4.1 and points low-coverage users back to 4.1 (Beagle 5.0 page).
   ng holds likelihoods in memory but does not emit PL/GL (`ng_facts.md`).
2. **Bi-allelic SNPs are the interface.** Every imputer/phaser surveyed restricts to bi-allelic sites
   (§7 table); GLIMPSE2 zeroes indel likelihoods by default. ng's multi-allelic SNP/indel/STR records
   would need splitting or projection.
3. **Read linkage adds a tenth to a third on top of LD at 4×** (HapSeq: 2.24% → 1.97% CEU, 0.86% → 0.60%
   simulated, Table 4/Table 1). ng's chain ids are the raw material; they do not survive into calling
   today (`ng_facts.md`).
4. **No panel, few founders: ~500 samples for r² 0.95; below 100 nothing is gained** (STITCH mice
   Fig. 4a text; owls Fig. 4a). Two-founder crosses are different: > 99% from 0.1× with 90 individuals
   (LB-Impute Fig. 2), 90% of crossovers within 2 kb at 0.4× (TIGER Table 2).
5. **Internal panel, low-diversity population: 50–150 deep individuals for r² ≈ 0.9 at 0.5–1×, 200–500
   for 0.95** (owl Fig. 3; strawberry Fig. 7; cattle Fig. 4). The panel must contain the target
   subpopulation (cattle, dogs, salmon).
6. **Depth: 0.3× → r² > 0.9 at MAF > 5%; 1× → 0.8 at MAF 0.1%; rare alleles need panel haplotypes, not
   depth** (GLIMPSE1 Fig. 2a; GLIMPSE2 Fig. 1).
7. **Phasing at N = 50–100 without a panel is unmeasured**; the nearest number is 1.9% SER with an
   85-sample panel at array density (Choi 2018 Table 2). Rare-variant phase stays poor even at N = 150,000
   (MAC 11–20: 4.4%, SHAPEIT5 Fig. 2a).
8. **INFO/DR2 do not detect regional failure** (König 2024 Supp. Fig. 1b) and DR2 is inflated at rare
   variants (Ramnarine 2015 Fig. 2). Report aggregate r² by MAF bin against held-out deep samples.
9. **Concordance on inbred material is not evidence of heterozygote accuracy** (strawberry 0.99 vs 0.88;
   soybean and wheat filtered heterozygotes out).
10. **"Call → LD refine → that is the call" is the 1000 Genomes design and gave 4× the accuracy of 15×
    single-sample calling for MAF > 1%** (pilot Fig. 2d; Phase 1 Fig. 1b). No published pipeline uses the
    imputed haplotypes to re-run discovery.

## 9. Questions to discuss

1. **Emit PL/GL in the VCF, or run the LD stage inside ng?** Recommendation: emit PL (and GP) first;
   it makes ng usable with GLIMPSE2/QUILT2/Beagle 4.1 today and costs a FORMAT field. Trade-off: an
   external stage sees only bi-allelic SNPs after splitting, and loses the chain ids (item 3).
2. **Which population case to build for first?** Recommendation: the segregant case. It is the one
   where the literature shows > 99% from 0.1× with ~90 individuals and a two-founder HMM, needs no
   panel, and its model (few founders, allelic bias, mis-mapping) is small. Trade-off: it serves the
   fewest users; the "many haplotypes" case needs either 500+ samples or an internal panel and a
   Li-Stephens machine that is a project in itself.
3. **Internal panel from the deep samples: how many are enough?** Recommendation: do not promise a
   number; expose held-out aggregate r² by MAF bin (leave-one-out among the deep samples) as the
   module's own report, because the literature's 50–150 figure is from low-diversity populations only.
4. **Inbred cohorts.** Recommendation: treat F as a per-sample parameter of the HMM (haploid model for
   samples with F near 1, as STITCH `diploid-inbred` does globally) rather than running a diploid model
   at F = 0.78. Trade-off: no published accuracy for the mixed case; it would have to be measured on
   the tomato cohort against the 20–40× accessions.
5. **Chain ids into calling.** Recommendation: retain, per sample and per heterozygous locus, the
   chain-id list per allele over a window, and measure on tomato and GIAB how much HapSeq-style linkage
   moves genotype discordance before designing anything larger. The published gain (9–30% relative) is
   worth the memory only if it reproduces at 3× on 63 samples.
6. **Indels and STRs.** Recommendation: impute indel and STR genotypes as GLIMPSE does — into the SNP
   haplotype scaffold, with the caller's own likelihoods kept (not flattened) because ng's STR
   likelihoods are model-based, not mpileup-based. Trade-off: unmeasured; the STR panel result (length r²
   0.91, 70% at hypervariable loci) is the ceiling from SNP haplotypes alone.

## 10. What I could not verify

- SHAPEIT4 Fig. 2a values at N = 500–200,000, Beagle 5.2 Fig. 2/3 values, SHAPEIT5 Extended Data
  Fig. 1 values, STITCH Fig. 4a per-cell r², QUILT2 Supplementary Tables 1–2, Lloret-Villas Fig. 4
  per-cell r²: all figure-only; I quote the papers' text statements instead.
- Teng 2022 (J Dairy Sci) full text is paywalled; only the abstract was read. The "plateau above 400
  samples and 1×" sentence comes from a search-engine summary.
- Pistis 2015's estimated-vs-true r² correlation (0.78–0.83) came from a search summary, not the paper.
- Vi 2025 per-population F1 values in the simulations are figure-only; only the hihi Table 1 values and
  the quoted sentences are verified.
- The 1000G pilot's exact LD-refinement methods (Beagle / MaCH / IMPUTE / QCALL) are in its supplement,
  which I did not open; the main text is quoted.
- ANGSD/ngsTools and TASSEL-GBS were not re-read; the statement that they have no LD stage is from
  their scope, not a citation.

## References

- 1000 Genomes Project Consortium 2010. A map of human genome variation from population-scale
  sequencing. Nature 467:1061. doi:10.1038/nature09534. PMC3042601.
- 1000 Genomes Project Consortium 2012. An integrated map of genetic variation from 1,092 human
  genomes. Nature 491:56. doi:10.1038/nature11632. PMC3498066.
- Browning BL, Tian X, Browning SR, Zhou Y 2021. Fast two-stage phasing of large-scale sequence data.
  Am J Hum Genet 108:1880. doi:10.1016/j.ajhg.2021.08.005. PMC8551421.
- Buckley RM et al. 2022. Best practices for analyzing imputed genotypes from low-pass sequencing in
  dogs. Mamm Genome 33:213. doi:10.1007/s00335-021-09914-z. PMC8913487.
- Choi Y, Chan AP, Kirkness E, Telenti A, Schork NJ 2018. Comparison of phasing strategies for whole
  human genomes. PLoS Genet 14:e1007308. doi:10.1371/journal.pgen.1007308. PMC5903673.
- Clark LV, Lipka AE, Sacks EJ 2019. polyRAD: genotype calling with uncertainty from sequencing data in
  polyploids and diploids. G3 9:663. doi:10.1534/g3.118.200913. PMC6404598.
- Crain J et al. 2026. Skim-sequencing for genomic selection in wheat: a comparison of marker
  platforms. Plant Genome. doi:10.1002/tpg2.70197. PMC12871549.
- Davies RW, Flint J, Myers S, Mott R 2016. Rapid genotype imputation from sequence without reference
  panels. Nat Genet 48:965. doi:10.1038/ng.3594. PMC4966640.
- Davies RW et al. 2021. Rapid genotype imputation from sequence with reference panels. Nat Genet
  53:1104. doi:10.1038/s41588-021-00877-0. PMC7611184.
- Davies RW et al. 2025. Flexible read-aware genotype imputation from sequence using biobank sized
  reference panels (QUILT2). Nat Commun. doi:10.1038/s41467-025-67218-1. PMC12804713; preprint
  doi:10.1101/2024.07.18.604149.
- Delaneau O, Zagury J-F, Robinson MR, Marchini JL, Dermitzakis ET 2019. Accurate, scalable and
  integrative haplotype estimation. Nat Commun 10:5436. doi:10.1038/s41467-019-13225-y. PMC6882857.
- Fragoso CA, Heffelfinger C, Zhao H, Dellaporta SL 2016. Imputing genotypes in biallelic populations
  from low-coverage sequence data. Genetics 202:487. doi:10.1534/genetics.115.182071. PMC4788230.
- Furuta T, Yamamoto T, Ashikari M 2023. GBScleanR: robust genotyping error correction using a hidden
  Markov model with error pattern recognition. Genetics 224:iyad055. doi:10.1093/genetics/iyad055.
  PMC10213493.
- Pierotti S, Welz B, Osuna-López M, Fitzgerald T, Wittbrodt J, Birney E 2024. Genotype imputation in F2
  crosses of inbred lines. Bioinform Adv 4:vbae107. doi:10.1093/bioadv/vbae107. PMC11286293.
- Gilly A et al. 2019. Very low-depth whole-genome sequencing in complex trait association studies.
  Bioinformatics 35:2555. doi:10.1093/bioinformatics/bty1032. PMC6662288.
- Gundappa MK et al. 2025. High performance imputation of structural and single nucleotide variants
  using low-coverage whole genome sequencing. Genet Sel Evol. doi:10.1186/s12711-025-00962-6.
  PMC11951665.
- Happ MM, Wang H, Graef GL, Hyten DL 2019. Generating high density, low cost genotype data in soybean.
  G3 9:2153. doi:10.1534/g3.119.400093. PMC6643887.
- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O 2023. Accurate rare variant phasing of
  whole-genome and whole-exome sequencing data in the UK Biobank. Nat Genet 55:1243.
  doi:10.1038/s41588-023-01415-w. PMC10335929.
- Huang L, Wang B, Chen R, Bercovici S, Batzoglou S 2016. Reveel: large-scale population genotyping
  using low-coverage sequencing data. Bioinformatics 32:1686. doi:10.1093/bioinformatics/btv530.
- Koorevaar T et al. 2025. Genotype imputation from low-coverage WGS using haplotype reference panels
  in cultivated strawberry. BMC Genomics. doi:10.1186/s12864-025-12270-w. PMC12629074.
- König E et al. 2024. Impact of the inaccessible genome on genotype imputation and genome-wide
  association studies. Hum Mol Genet. doi:10.1093/hmg/ddae062. PMC11227617.
- Li Y, Sidore C, Kang HM, Boehnke M, Abecasis GR 2011. Low-coverage sequencing: implications for
  design of complex trait association studies. Genome Res 21:940. doi:10.1101/gr.117259.110.
  PMC3106327.
- Lin P et al. 2010. A new statistic to evaluate imputation reliability. PLoS One 5:e9697.
  doi:10.1371/journal.pone.0009697. PMC2837741.
- Lloret-Villas A, Pausch H, Leonard AS 2023. The size and composition of haplotype reference panels
  impact the accuracy of imputation from low-pass sequencing in cattle. Genet Sel Evol 55:33.
  doi:10.1186/s12711-023-00809-y. PMC10173671.
- Loh P-R et al. 2016. Reference-based phasing using the Haplotype Reference Consortium panel. Nat
  Genet 48:1443. doi:10.1038/ng.3679. PMC5096458.
- Oget-Ebrad C, Kadri NK, Moreira GCM, Karim L, Coppieters W, Georges M, Druet T 2022. Benchmarking
  phasing software with a whole-genome sequenced cattle pedigree. BMC Genomics 23:130.
  doi:10.1186/s12864-022-08354-6. PMC8845340.
- Pasaniuc B et al. 2012. Extremely low-coverage sequencing and imputation increases power for
  genome-wide association studies. Nat Genet 44:631. doi:10.1038/ng.2283. PMC3400344.
- Pistis G et al. 2015. Rare variant genotype imputation with thousands of study-specific whole-genome
  sequences: implications for cost-effective study designs. Eur J Hum Genet 23:975. PMC4463504.
- Ramnarine S et al. 2015. When does choice of accuracy measure alter imputation accuracy assessments?
  PLoS One 10:e0137601. doi:10.1371/journal.pone.0137601.
- Topaloudis A, Cumer T, Lavanchy E, … Delaneau O, Goudet J 2026. Benchmarking imputation accuracy in the
  presence or absence of a reference panel. Mol Biol Evol 43:msag094. doi:10.1093/molbev/msag094.
  PMC13122032; preprint doi:10.1101/2025.10.08.680883.
- Rowan BA, Patel V, Weigel D, Schneeberger K 2015. Rapid and inexpensive whole-genome
  genotyping-by-sequencing for crossover localization and fine-scale genetic mapping. G3 5:385.
  doi:10.1534/g3.114.016501. PMC4349092.
- Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O 2021. Efficient phasing and imputation of
  low-coverage sequencing data using large reference panels. Nat Genet 53:120.
  doi:10.1038/s41588-020-00756-0; preprint doi:10.1101/2020.04.14.040329 (v1 text quoted).
- Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O 2023. Imputation of low-coverage sequencing
  data from 150,119 UK Biobank genomes. Nat Genet 55:1088. doi:10.1038/s41588-023-01438-3.
  PMC10335927.
- Saini S, Mitra I, Mousavi N, Fotsing SF, Gymrek M 2018. A reference haplotype panel for genome-wide
  imputation of short tandem repeats. Nat Commun 9:4397. doi:10.1038/s41467-018-06694-0. PMC6199332.
- Sthapit SR et al. 2025. A low-coverage skim-sequencing and imputation pipeline for genomic selection.
  Plant Genome. doi:10.1002/tpg2.70139. PMC12547641.
- Teng J et al. 2022. Assessment of the performance of different imputation methods for low-coverage
  sequencing in Holstein cattle. J Dairy Sci 105:3355. doi:10.3168/jds.2021-21360.
- Vi T et al. 2025. Assessing genotype imputation methods for low-coverage sequencing data in
  populations with differing relatedness and inbreeding levels. Mol Ecol Resour.
  doi:10.1111/1755-0998.70049. PMC12550480.
- Wang Y, Lu J, Yu J, Gibbs RA, Yu F 2013. An integrative variant analysis pipeline for accurate
  genotype/haplotype inference in population NGS data (SNPTools). Genome Res 23:833.
  doi:10.1101/gr.146084.112. PMC3638139.
- Wasik K, Berisa T, Pickrell JK et al. 2021. Comparing low-pass sequencing and genotyping for trait
  mapping in pharmacogenetics. BMC Genomics 22:197. doi:10.1186/s12864-021-07508-2. PMC7981957.
- Watowich MM et al. 2025. Best practices for genotype imputation from low-coverage sequencing data in
  natural populations. Mol Ecol Resour. doi:10.1111/1755-0998.13854. PMC10879460.
- Beck AT, Kang HM, Zöllner S 2025 (preprint; journal version 2026). A benchmark of modern statistical
  phasing methods. PMC13131814; preprint doi:10.1101/2025.06.24.660794.
- Zhi D, Wu J, Liu N, Zhang K 2012. Genotype calling from next-generation sequencing data using haplotype
  information of reads (HapSeq). Bioinformatics 28:938. doi:10.1093/bioinformatics/bts047. PMC3493122.
- Zheng C, Boer MP, van Eeuwijk FA 2018. Accurate genotype imputation in multiparental populations
  from low-coverage sequence. Genetics 210:71. doi:10.1534/genetics.118.300885. PMC6116951.
- Biagini SA, Becelaere S, Aerden M, … Kivisild T 2025. Genotype imputation from low-coverage data for
  medical and population genetic analyses. Genome Res 35:1929. doi:10.1101/gr.280175.124. PMC12400947.
- GATK. Genotype Refinement workflow for germline short variants.
  github.com/broadinstitute/gatk-docs/…/Genotype_Refinement_workflow.md.
- Beagle 5.0 release page. faculty.washington.edu/browning/beagle/b5_0.html.
- Illumina 2024. Boosting variant calling performance using a high-quality reference panel for
  imputing low-coverage sequencing data (DRAGEN + IRPv1).
