# Imputation from low-coverage sequencing: GLIMPSE2, QUILT/QUILT2, STITCH and Beagle 4.1

Research sub-report for the ng phasing/imputation module. Written 2026-09-12 from the cloned
sources and the papers listed in §0. Every claim about code carries a `path:line` reference into
`tmp/phasing_research/repos/` (paths below are relative to that directory); every claim about a
method carries a paper reference. Where a paper could not be opened, §13 says so.

The report answers, for each tool: how the read evidence enters the model, where the haplotypes it
copies from come from, how the iterations are organised, what comes out, what it costs, and what
accuracy its authors published. §10 asks whether the panel-based and the panel-free tools are the
same machine; §11 says what ng could hand to each tool today and where an in-caller module would do
better.

## 0. Sources

Code (shallow clones, 2026-09-12):

| directory | tool and version | language |
|---|---|---|
| `GLIMPSE/` | GLIMPSE2, fork at commit 8671138 (merge of `tfenne/tf_phasing_bitbatch`, so a few comments are local to the fork) | C++ |
| `QUILT/` | QUILT 2.0.4 (`QUILT/DESCRIPTION:4`), commit a8a7ed3, 2026-02-11; contains the QUILT1 path too | R + C++ |
| `STITCH/` | STITCH 1.8.5 (`STITCH/DESCRIPTION:4`), commit e5dc688, 2026-01-19 | R + C++ |

A previous reading of GLIMPSE2 is summarised in
`tmp/phasing_research/low_coverage_imputation/glimpse2_code_facts.md`; the lines cited from it
below were re-opened before citing (the emission, the likelihood collapse, the transition, the
panel gate, the selection defaults, the iteration counts, the INFO score, the chunk defaults and
the SNP-only flag).

Papers:

- STITCH: Davies RW, Flint J, Myers S, Mott R. *Rapid genotype imputation from sequence without
  reference panels.* Nat Genet 2016;48(8):965–969. doi:10.1038/ng.3594 (PMC4966640).
- QUILT: Davies RW, Kucka M, Su D, Shi S, Flanagan M, Cunniff CM, Chan YF, Myers S. *Rapid genotype
  imputation from sequence with reference panels.* Nat Genet 2021;53:1104–1111.
  doi:10.1038/s41588-021-00877-0 (PMC7611184).
- QUILT2: Li Z, Albrechtsen A, Davies RW. *Flexible read-aware genotype imputation from sequence
  using biobank sized reference panels.* Nat Commun 2025, article 524.
  doi:10.1038/s41467-025-67218-1 (PMC12804713; PubMed 41390671). Preprint: bioRxiv
  2024.07.18.604149, *Rapid and accurate genotype imputation from low coverage short read, long read,
  and cell free DNA sequence*.
- GLIMPSE (v1): Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O. *Efficient phasing and
  imputation of low-coverage sequencing data using large reference panels.* Nat Genet
  2021;53:120–126. doi:10.1038/s41588-020-00756-0. Read from the preprint, bioRxiv
  10.1101/2020.04.14.040329 (the journal page redirects to a login).
- GLIMPSE2: Rubinacci S, Hofmeister RJ, Sousa da Mota B, Delaneau O. *Imputation of low-coverage
  sequencing data from 150,119 UK Biobank genomes.* Nat Genet 2023;55:1088–1090.
  doi:10.1038/s41588-023-01438-3 (PMC10335927).
- Beagle 4.1: Browning BL, Browning SR. *Genotype imputation with millions of reference samples.*
  Am J Hum Genet 2016;98(1):116–126. doi:10.1016/j.ajhg.2015.11.020 (PMC4716681). Full text not
  reachable (§13); facts below come from the Beagle 4.0 and 4.1 manuals
  (`tmp/phasing_research/low_coverage_imputation/beagle_4.0_03Mar15.txt`, `beagle_4.1_09Feb16.txt`)
  and the Beagle 5.4 manual (`beagle_5.4_18Mar22.txt`, downloaded today).
- Beagle genotype-likelihood lineage: Browning BL, Yu Z. *Simultaneous genotype calling and haplotype
  phasing improves genotype accuracy and reduces false-positive associations for genome-wide
  association studies.* Am J Hum Genet 2009;85(6):847–861. doi:10.1016/j.ajhg.2009.11.004. Full
  text not reachable (§13).
- loimpute: Wasik K, Berisa T, Pickrell JK, Li JH, Fraser DJ, King K, Cox C. *Comparing low-pass
  sequencing and genotyping for trait mapping in pharmacogenetics.* BMC Genomics 2021;22:197.
  doi:10.1186/s12864-021-07508-2 (bioRxiv 2019, 10.1101/632141). The model is in its Supplementary
  Note, downloaded to `tmp/phasing_research/low_coverage_imputation/wasik2021_loimpute_supp.txt`.
- GeneImp: Spiliopoulou A, Colombo M, Orchard P, Agakov F, McKeigue P. *GeneImp: Fast imputation to
  large reference panels using genotype likelihoods from ultralow coverage sequencing.* Genetics
  2017;206(1):91–104. doi:10.1534/genetics.117.200063 (PMC5419496).
- Reveel: Huang L, Wang B, Chen R, Bercovici S, Batzoglou S. *Reveel: large-scale population
  genotyping using low-coverage sequencing data.* Bioinformatics 2016;32(11):1686–1696.
  doi:10.1093/bioinformatics/btv530.

## 1. Vocabulary

Genetics terms are kept brief; statistics and software terms are explained.

- **Genotype likelihood (GL).** For one sample at one site, the probability of that sample's reads
  under each possible genotype: P(reads | 0/0), P(reads | 0/1), P(reads | 1/1). VCF stores it as
  `GL` (log10) or `PL` (Phred-scaled, normalised so the best genotype is 0). It says how well each
  genotype explains the reads; it is not a probability that the genotype is true.
- **Haplotype likelihood.** The same idea for one chromosome copy: P(reads | this copy carries allele
  0) and P(reads | allele 1). Needed by the haploid models below.
- **Dosage (DS).** Expected number of alternate alleles, 0–2, computed from the genotype posterior:
  P(0/1) + 2·P(1/1).
- **Genotype posterior (GP).** The three probabilities P(0/0), P(0/1), P(1/1) after imputation.
- **Hidden Markov model (HMM).** A model where an unobserved state walks along the chromosome
  (here: *which haplotype the sample's chromosome is currently copying*), switching states with
  probabilities that depend on genetic distance (*transition*), and where the observed data at each
  site (reads or GLs) are explained by the current state (*emission*).
- **Forward–backward.** The standard HMM algorithm that computes, at every site, the posterior
  probability of each hidden state given all the data on the chromosome. Its cost is (number of
  sites) × (number of states) for a haploid model, or × (states)² when a pair of states walks
  together.
- **Li and Stephens copying model.** The HMM used by every tool here: a new chromosome is modelled
  as a mosaic of copies of the *template* haplotypes, with occasional switches (recombinations)
  and occasional mismatches (mutation or error). Li N, Stephens M. Genetics 2003;165:2213–2233.
- **Template haplotypes / conditioning haplotypes.** The set of haplotypes the HMM is allowed to
  copy from. GLIMPSE calls them *conditioning states*; QUILT calls them the *reference panel* or
  the *small reference panel*; STITCH calls them *ancestral* or *founder* haplotypes.
- **Reference panel.** A file of phased haplotypes from well-genotyped individuals, used as the
  template set.
- **PBWT (positional Burrows–Wheeler transform).** An index over a set of haplotypes that sorts
  them, at every site, by the alleles they carry going leftwards from that site; neighbours in the
  sorted order share the longest match ending at that site. Used to find, for a target haplotype,
  the panel haplotypes that match it longest — cheaply.
- **Gibbs sampling.** An iterative scheme where one unknown at a time is redrawn at random from its
  probability given the current values of all the others. After a *burn-in*, the draws from the
  *main* iterations are averaged to estimate posteriors.
- **EM (expectation–maximisation).** An iterative scheme that alternates between (E) computing
  posteriors of the hidden states under the current model parameters and (M) re-estimating the
  parameters from those posteriors. STITCH uses it to learn its ancestral haplotypes.
- **Read label.** In QUILT: for each read, which of the sample's two chromosomes it came from.
- **INFO score.** A per-site imputation quality measure that needs no truth: 1 minus the average
  posterior variance of the dosage across samples divided by its variance under Hardy–Weinberg at
  the estimated frequency. 1 = every sample's dosage is certain; 0 = no information. It is the
  IMPUTE "I_A" measure; §7 gives the formula each tool uses.
- **Switch error.** In phased output, the fraction of consecutive heterozygous sites whose phase is
  flipped relative to the truth.

## 2. The shared machine

All four tools run the same core: a *haploid* Li and Stephens HMM whose hidden state, at each site,
is *which template haplotype this chromosome copy is copying*. The transition probability between
adjacent sites is set by the genetic distance and a population parameter (GLIMPSE: effective
population size Ne; STITCH and QUILT: generations since founding, `nGen`; Beagle: its haplotype
cluster model, see §3.4). The emission probability is the chance of the observed evidence given the
allele carried by the copied template. A diploid sample is handled in one of two ways:

- **Two haploid passes conditioned on each other** (GLIMPSE, Beagle-4-style sampling, QUILT):
  sample or fix one chromosome copy, then run the haploid HMM for the other with the evidence
  "explained away" by the first. QUILT goes one level deeper and assigns each *read* to a copy.
- **One diploid pass over pairs of templates** (STITCH `diploid`, loimpute): the hidden state is a
  pair (k1, k2), the cost is K² states per site, and the emission is for the genotype the pair
  implies.

What differs is **where the templates come from** and **whether they move**:

| tool | templates | do they change during the run? | how many are used per sample |
|---|---|---|---|
| GLIMPSE2 | reference panel only | no (panel is static, PBWT pre-built) | Kinit 1000 then Kpbwt 2000 chosen per sample and per iteration |
| GLIMPSE1 | reference panel + the other target samples' current haplotypes | yes (targets re-sampled each iteration) | K = 1000 default (preprint) |
| QUILT / QUILT2 | reference panel only | no | Ksubset 600 (code; 400 in the 2021 paper) for the Gibbs step; whole panel in the full step |
| STITCH | K ancestral haplotypes learned from all samples' reads (optionally seeded from a panel) | yes — that is the point | all K, for every sample |
| Beagle 4.1 `gl=` | a haplotype-cluster model built from the current haplotype estimates of all samples, plus a panel if given | yes | model built per window |

§10 returns to this table.

## 3. How genotype likelihoods and reads enter

### 3.1 GLIMPSE2: per-site GLs, collapsed to haploid likelihoods

GLIMPSE2 takes either a VCF/BCF with `PL` (default) or `GL` (flag `--input-field-gl`), or BAM/CRAM
files, which it pileups itself with an internal caller (`GLIMPSE/phase/src/caller/caller_parameters.cpp:124-127`;
PL read at `GLIMPSE/phase/src/io/genotype_reader.cpp:467`, GL at `:430`). Each genotype's likelihood
is stored as one byte, Phred-scaled and clamped to 0–255 (`genotype_reader.cpp:456-461`). A site
with missing GLs is marked *flat* and contributes no evidence (`genotype_reader.cpp:478-480`,
`GLIMPSE/phase/src/objects/genotype.cpp:85-86`).

The HMM is haploid, so the diploid triple must become a pair. Given the *other* chromosome copy's
current allele at the site, the triple collapses to two numbers (`genotype.cpp:94-115`): if the
other copy carries allele 0, the haploid likelihoods are proportional to (GL00, GL01); if it carries
allele 1, to (GL01, GL11). Both are floored at `min_gl` (`:112-113`). In the first iteration and for
haploid samples the marginal is used instead (`genotype.cpp:73-81`).

The emission for a template state carrying allele *a* is then the haploid likelihood of *a* mixed
with a small copying-error rate *e* (`--err-imp`, clamped to [1e-12, 1e-3] at
`GLIMPSE/phase/src/caller/caller_initialise.cpp:152`):

```
GLIMPSE/phase/src/models/imputation_hmm.cpp:61-71
p0 = HL0·(1−e) + HL1·e ;  p1 = HL0·e + HL1·(1−e)
Emissions[site][0] = p0/(p0+p1) ;  Emissions[site][1] = p1/(p0+p1)
```

Transition between consecutive polymorphic sites: probability of switching template is
t = 1 − exp(nrho · d_cM), with nrho = −0.04·Ne / n_ref_haps and Ne = 100,000 by default
(`GLIMPSE/phase/src/containers/conditioning_set.cpp:34`, `:190-194`;
`caller_initialise.cpp:142,158`). Without a genetic map, 1 cM per Mb is assumed
(`GLIMPSE/common/src/containers/variant_map.cpp:154`).

Per sample and per iteration (`GLIMPSE/phase/src/caller/caller_algorithm.cpp:58-72`): build haploid
likelihoods for copy 0 conditional on the current copy 1 → forward–backward → sample copy 0 → build
likelihoods for copy 1 conditional on the new copy 0 → forward–backward → sample copy 1 → run a
separate SHAPEIT-style phasing HMM (`GLIMPSE/phase/src/models/phasing_hmm.cpp`) that re-phases the
heterozygous sites of the sampled pair.

Indels in the input are ignored unless `--use-gl-indels` (`genotype_reader.cpp:426`,
`caller_parameters.cpp:56`), and multi-allelic records are skipped (`genotype_reader.cpp:265, :371,
:409`; `GLIMPSE/common/src/io/ref_genotype_reader.cpp:159, :225`). Only SNPs that are the first
record at their position count as "high quality" for the HMM scaffold
(`ref_genotype_reader.cpp:181-183`); other variants are imputed as passengers (§7.4).

### 3.2 QUILT: per-read evidence, and reads that span several sites

QUILT reads BAM/CRAM directly, using STITCH's read extractor. This is the part closest to ng's
per-read allele evidence, so it is described in full.

**The read record.** `STITCH/STITCH/src/bam_access.cpp:124-400` walks each alignment's CIGAR and,
for every SNP of the site list that the aligned bases cover (`:236-238` for the span, `:294-330` for
the CIGAR walk), keeps the base only if it equals the SNP's REF or ALT and its base quality is at
least `bqFilter` (17 by default) after capping the base quality at the read's MAPQ (`:326-329`).
Each kept base becomes one entry: the SNP index and a **signed base quality** — negative if the
base was REF, positive if ALT (`:338-341`). Mates of a pair and reads sharing a `BX` barcode are
merged into one read record by query name (`STITCH/STITCH/R/functions.R:1923-1945`), so one record
can span the whole insert. The record is a 4-element list: number of SNPs minus one, the *central
SNP* (the SNP nearest the middle of the covered SNPs, `functions.R:2522-2526`,
`STITCH/STITCH/R/sample-reads.R:590-591`), the vector of signed base qualities and the vector of
SNP indices (`functions.R:2565`). Sites are downsampled to at most `downsampleToCov` reads
(30 in QUILT, `QUILT/QUILT/R/quilt.R:147`; 50 in STITCH, `functions.R:101`).

**Per-read emission.** For read *r* and template *k*, the probability of the read given the
template is the product over the SNPs *j* in the read of the probability of the observed base
given the template's allele (`STITCH/STITCH/src/haploid.cpp:129-215`, identical code in
`QUILT/QUILT/src/copied-from-stitch.cpp:167-174`). With ε = 10^(−BQ/10):

```
haploid.cpp:178-189
observed REF:  pR = 1−ε, pA = ε/3          observed ALT:  pR = ε/3, pA = 1−ε
eMatRead[k, r] ∏=  θ[k, j]·pA + (1−θ[k, j])·pR
```

Here θ[k, j] is the probability that template *k* carries the ALT allele at SNP *j*: in QUILT it is
0 or 1 from the panel softened by `ref_error` (0.001, `QUILT/QUILT/R/functions.R:2323`); in STITCH
it is the learned ancestral-haplotype parameter (§5). The vector over *k* is then rescaled so its
maximum is 1 and floored at 1/`maxDifferenceBetweenReads` (`haploid.cpp:202-235`), which caps how
much any one read can pull the state posterior — the parameter doc says this "helps prevent against
influence by false positive SNPs" (`functions.R:26`).

**A read spanning several sites.** The read is one observation, attached to the grid cell of its
central SNP; the template is assumed constant over the read (Davies 2021, Methods: "each read has a
central location, and genotypes in sequencing reads are based on the reference haplotype copied at
that central location"). The evidence from all its SNPs multiplies into `eMatRead[·, r]`, and that
whole vector is multiplied into the emission of the central grid cell. So a read carrying REF at
site 1 and ALT at site 2 down-weights every template that is (ALT, REF) or (ALT, ALT) or
(REF, REF) at those two sites — the read-backed phase is used, but a recombination inside a read
is not modelled. STITCH's paper measured the value of this: on the CFW mice at 0.15×, dropping the
read structure (splitting reads one SNP per read) lowered r² against the array from 0.97 to 0.87
(Davies 2016, Supplementary Table 2); on the CONVERGE humans at 1.7× it made no difference, which
the authors attribute to the lower SNP density. The `readAware = FALSE` switch does exactly this
split (`functions.R:20, :492-494`).

**Read labels (QUILT's Gibbs sampler).** QUILT keeps, for each read, a label H(r) ∈ {1, 2} saying
which chromosome copy it came from, and runs *two* haploid HMMs, one per copy, each seeing only the
reads currently labelled to it. The emission of copy *h* at grid *g* is the product of `eMatRead`
over reads labelled *h* whose central SNP is in *g*. Gibbs sampling then revisits the reads one grid
at a time (`QUILT/QUILT/src/gibbs-nipt.cpp:733-1110`, `sample_reads_in_grid`):

- With forward (α) and backward (β) vectors of both copies available at grid *g*, the probability
  of the whole data under the current labels is proportional to Σ_k α_h(k)·β_h(k) for each copy *h*
  (`:846-857`, stored as `pC`).
- Moving read *r* from its current copy to the other one changes those two sums by dividing out or
  multiplying in `eMatRead[·, r]` state by state (`:909-911`; the loops at `:913-960` do the same
  restricted to the templates where the read is informative).
- The two alternatives (keep, flip) are weighted by the product of the copy likelihoods times a
  prior on the label (`:998-1003`) and one is drawn (`:1032-1046`).
- If the label changes, α, the α·β product and the copy's grid emission are updated in place
  (`:1089-1110`); no full forward–backward is rerun for a single flip.

This is why the paper can claim "linear computational complexity in the number of samples, SNPs and
reference haplotypes" for the Gibbs step (Davies 2021, Methods). *Block* Gibbs moves — flipping the
labels of all reads in a stretch of the chromosome at once, to escape a locally consistent but
globally wrong phase — run at the iterations listed in `small_ref_panel_block_gibbs_iterations`
(default `c(3, 6, 9)`, `quilt.R:167`; driver `gibbs-nipt.cpp:2971-3025`).

**From labelled reads to per-site haploid likelihoods.** After the Gibbs step, QUILT runs each copy
against the *full* panel. For that it converts the reads labelled to copy *h* into a 2 × nSNPs
haploid likelihood table: for every base in every read of that copy, multiply the site's
(P(base | REF), P(base | ALT)) pair from the base quality into the site (`QUILT/QUILT/R/functions.R:2015-2024`,
`QUILT/QUILT/R/reference-single.R:19-43`, floor `minGLValue` = 1e-10). At that point read
membership no longer matters — the paper says this explicitly: once reads are known to come from
the same haplotype "membership of a sequenced base in a read no longer matters" (Davies 2021,
Methods). The full-panel HMM then uses these per-site pairs against each distinct 32-SNP panel
segment: `QUILT/QUILT/src/reference-single.cpp:296-322` computes, for every distinct haplotype
pattern in a 32-SNP grid, ∏ over its SNPs of (dR·(1−ref_error) + dA·ref_error) if the pattern has
REF, or (dR·ref_error + dA·(1−ref_error)) if ALT.

So QUILT has *both* interfaces: reads (Gibbs step) and per-site haploid GLs (full-panel step). ng
could feed either.

### 3.3 STITCH: the same per-read emission, but the templates are parameters

STITCH uses the read record and the per-read emission of §3.2 unchanged — the code is STITCH's.
The difference is what θ[k, j] is: not a panel allele but the *k*-th ancestral haplotype's
probability of carrying ALT at SNP *j*, a parameter updated every iteration (§5). In the default
`diploid` method the hidden state is a pair (k1, k2), and the emission at a grid cell is the product
over its reads of

```
STITCH/STITCH/src/diploid.cpp:138-143
eMatGrid[(k1,k2), g] ∏= 0.5·eMatRead[k1, r] + 0.5·eMatRead[k2, r]
```

— each read came from one of the two copies with equal prior probability; the sum over the two
possibilities is taken per read, which is what makes the diploid model exact without read labels
and what makes it cost K² states. The per-cell vector is rescaled and floored at
1/`maxEmissionMatrixDifference` (`diploid.cpp:155-167`).

### 3.4 Beagle 4.1 `gl=`: per-site GLs into a haplotype-cluster model

Beagle 4.0 and 4.1 accept a VCF with `GL` or `PL` through `gl=` (GT ignored) or `gtgl=` (GT used
where present, GL otherwise) (`beagle_4.1_09Feb16.txt:98-104`). The output genotypes in this mode
are **unphased**; phasing needs a second run with `gt=` on the output (`:78-82`). Likelihoods more
than `maxlr` = 5000-fold below the best genotype are set to zero (`:122-124`). The manual notes that
`modelscale` (default 0.8) "when estimating posterior probabilities from genotype likelihood data
… could improve both accuracy and run-time" if increased (`:169-173`).

The model is not Li and Stephens over a fixed panel. Beagle 4.0 builds a *localized haplotype
cluster model* — a graph summarising the haplotype frequencies in the current estimate of all
samples' haplotypes in a window — and samples `nsamples` = 4 haplotype pairs per individual per
iteration from that model weighted by the genotype likelihoods (`beagle_4.0_03Mar15.txt:205-208`,
window of 1200 markers `:208-209`). Beagle 4.1 keeps that as its 10 burn-in iterations and adds 5
iterations of its 4.1 phasing algorithm (`beagle_4.1_09Feb16.txt:144-147`). The 4.1 imputation
paper (Browning & Browning 2016) describes imputation from *phased* target data into a reference
panel with a Li and Stephens model over aggregate markers; the GL-sampling mechanism descends from
Browning & Yu 2009 (genotype calling from array intensity likelihoods by sampling haplotype pairs
from the cluster model). Neither full text could be opened (§13), so this paragraph rests on the
manuals and on the two papers' titles and abstracts only.

Beagle 5 dropped the mode: "Beagle version 5 does not infer genotypes from genotype likelihood
input data, but Beagle versions 4.0 and 4.1 have this capability" (`beagle_5.4_18Mar22.txt:32-33`).

What matters for ng: Beagle 4.x is the one production tool that imputes a cohort *from GLs without
any panel*, using the cohort's own current haplotypes as the template model — the "cohort as its own
panel" case in the brief. Its cost is the problem (§8, §9).

### 3.5 loimpute, GeneImp, Reveel

**loimpute** (Gencove; Wasik 2021, Supplementary Note). A diploid Li and Stephens HMM: the hidden
state is the pair of copied panel haplotypes, transitions use ρ = 0.0001 per base (fixed) and the
panel size M; the emission is *per read*, not per site GL: with Y the number of ALT alleles in the
copied pair, a read showing allele *r* has probability 1−ε if it matches a homozygous state, ε if it
contradicts it, and 0.5 if Y = 1 (ε = 0.001). Cost is O(N·M²), so before running it picks, per
region, the 100 panel haplotypes that share the most rare (< 5 % frequency) ALT alleles with the
target's reads (statistic: matches / ALT count on the haplotype). No Gibbs, no read labels, no
learning. Accuracy (Wasik 2021, Fig. 2; 79 volunteers, Cambridge UK; 1000 Genomes Phase 3 panel):
r² 0.96 at 1× and 0.91 at 0.4× for MAF > 5 %; 0.93 and 0.88 across all frequencies. A Genome
Research 2025 paper by the same group (doi:10.1101/gr.280175.124) is behind a login and was not
read.

**GeneImp** (Spiliopoulou 2017). No HMM at all: in each sliding window the target's two haplotypes
are assumed to be *present unaltered* among the panel haplotypes; a posterior over panel-haplotype
pairs is computed from the product of the GLs across the window, and windows are combined by a
mean-field approximation. Tested only at about 0.5× (0.45–0.76×, 16 samples) against 1000
Genomes Phase 3 (2,504) and a combined 5,334-individual panel: mean r² rises steeply with MAF and
levels off above 5 % (Fig. 1B); it sits about 0.02 below Beagle 4.0 in every MAF bin (Fig. 2) at 15
to 90 times less run time (Table 5: chromosome 22, 16 samples, 1.5 h single-window vs Beagle 4.0
135.6 h).

**Reveel** (Huang 2016). Reference-free, cohort-only, not an HMM: for every SNP it finds the *k*
SNPs in strongest LD within 0.5–1 Mb, estimated from read counts, and iterates
"summarisation–maximisation" (genotype probabilities given the neighbours' genotypes and the reads;
then the most probable genotypes), 3 rounds of 10 iterations with neighbours re-picked between
rounds. It is linear in the number of samples. Reported on simulated 1000G-like data at 7.4×
(!): genotype accuracy 99.84 % at n = 100 rising to 99.98 % at n = 1000, and 99.71 % at 2× with
n = 1000 (Supplementary Table S6); 76–192× faster than Thunder. It is a genotype caller that uses
LD, not a phaser, and its evaluated coverage is far above ng's low end; it matters here only as a
demonstration that cohort-only LD genotyping exists and scales linearly.

## 4. The panel: representation and selection

### 4.1 GLIMPSE2

`GLIMPSE2_split_reference` converts a panel VCF into a per-chunk binary file: variants with MAF
≥ `--sparse-maf` 0.001 go into a dense bit matrix, rarer ones into per-variant lists of minor-allele
carriers (`GLIMPSE/common/src/io/ref_genotype_reader.cpp:171-174`;
`GLIMPSE/common/src/containers/ref_haplotype_set.h:119-121`). The PBWT over the common variants is
built once, run-length compressed and stored in the same file (`ref_haplotype_set.h:42-54,
:123-124`; `ref_haplotype_set.cpp:47-103`), together with the genetic map, so those flags are
refused at imputation time (`GLIMPSE/phase/src/caller/caller_parameters.cpp:236-249`).

Per sample and per iteration, the conditioning set is (`GLIMPSE/phase/src/caller/caller_parameters.cpp:67-73`
for the defaults): iteration 0 — `--Kinit` 1000 haplotypes chosen by shared rare alleles (sort panel
haplotypes by minor-allele count over the rare alleles the sample's GLs suggest, sample, top up
uniformly: `GLIMPSE/phase/src/containers/haplotype_set.cpp:775-833`); later iterations — up to
`--Kpbwt` 2000 haplotypes found by threading the sample's current haplotype pair through the stored
PBWT and taking `--pbwt-depth` 12 neighbours on each side at positions spaced every
`--pbwt-modulo-cm` 0.1 cM (`haplotype_set.cpp:317-394`, `:445-490`). If Kpbwt ≥ panel size the
whole panel is used (`conditioning_set.cpp:75-182`). A `--state-list` can pin haplotypes (`:115-121`).

**The other target samples are not templates in GLIMPSE2.** The PBWT arrays are sized by the
number of *reference* haplotypes (`haplotype_set.cpp:244-246, :272, :330`), the selection guards
`idx < n_ref_haps` (`:712-714`), the emission dereferences only the reference bit matrices
(`conditioning_set.cpp:130, :157-168`), the targets' own haplotypes are used only as the *query*
strings (`haplotype_set.cpp:452, :501-503`), and a panel is mandatory
(`caller_parameters.cpp:133`). The paper confirms the design change: "Conversely from the GLIMPSE1
method, GLIMPSE2 approach is primarily focused on imputation only from the reference panel"
(Rubinacci 2023, Methods). GLIMPSE1 did use them — its preprint states "At each iteration, a new
pair of haplotypes is sampled per target individual by conditioning on haplotypes from the
reference panel and from other target individuals" (Rubinacci 2020 preprint, Methods) — with K =
1000, `pbwt-modulo` 8 and `pbwt-depth` 2 as defaults.

### 4.2 QUILT and QUILT2

The panel is loaded from IMPUTE `.hap/.legend` or a VCF by `QUILT_prepare_reference`, which skips
anything not a bi-allelic SNP at a unique position (`QUILT/QUILT/R/quilt-prepare-reference.R:246`),
and packs it as `rhb_t`: one 32-bit integer per haplotype per **grid** of 32 consecutive SNPs
(`STITCH/STITCH/R/grid.R:16-35`, `quilt-prepare-reference.R:378`). Grids are the unit of the HMM —
transitions happen only between grids, and a read's central SNP places it in a grid. For the
full-panel step, each grid's distinct 32-bit patterns are tabulated (`distinctHapsB`, at most
`nMaxDH` of them, recommended 255 or 511; `hapMatcher` maps each panel haplotype to its pattern;
`quilt-prepare-reference.R:18, :416-426`), so the emission is computed once per distinct pattern
instead of once per haplotype (§3.2).

Selection is iterative. Per sample (`QUILT/QUILT/R/functions.R:575-960`):

1. First "seek" iteration: `Ksubset` (600, `quilt.R:113`) panel haplotypes are drawn **at random**
   (`functions.R:580`).
2. Gibbs sampling of read labels against those Ksubset haplotypes (§3.2), `small_ref_panel_gibbs_iterations`
   = 20 passes (`quilt.R:168`).
3. The two labelled read sets are each run as a haploid HMM against the **whole** panel
   (`impute_using_everything`, `functions.R:1922`). At every `heuristic_match_thin` = 10 % of grids
   the `K_top_matches` = 5 best-matching panel haplotypes are recorded
   (`functions.R:2207-2259`); their union, thinned to `Knew` (600) new haplotypes, replaces a random
   Ksubset − Knew subset of the current set (`functions.R:746, :951`, `everything_select_good_haps`
   `:2262-2311`). In QUILT2 with `use_mspbwt = TRUE` (`quilt.R:173`, default FALSE in this code)
   the whole-panel HMM pass is replaced by an msPBWT query of the current haplotype estimate over
   the 32-SNP symbols, ranking matches by length and rarity (`QUILT/QUILT/R/mspbwt.R:94-131`).
4. Repeat for `n_seek_its` = 3 iterations (`quilt.R:111`); only the last contributes to the
   output unless `n_burn_in_seek_its` is lowered (`quilt.R:248-250`).
5. The whole thing runs `nGibbsSamples` = 7 times from different random starts (`quilt.R:110`,
   loop `functions.R:400`), and the dosages are averaged.

With a panel smaller than Ksubset, the code sets Ksubset = Knew = panel size and n_seek_its = 1
(`quilt.R:451-465`) — the panel is simply used whole.

Target samples' haplotypes are never templates in QUILT: every template index points into `rhb_t`,
and samples are processed one at a time (`get_and_impute_one_sample`, `functions.R:3`). QUILT2's
two-stage mode (`impute_rare_common`, threshold `rare_af_threshold` 0.001; `quilt.R:180-181`)
first runs the whole procedure on common SNPs only, then rebuilds the selected templates with their
rare alleles from per-haplotype lists (`QUILT/QUILT/R/rare_common.R:1-58`) for one final Gibbs pass
over all reads (`rare_common.R:61-107` seeds the labels).

### 4.3 STITCH with a panel

STITCH can read the same panel format, but uses it **only to initialise** the K ancestral
haplotypes; after the first iteration the panel is discarded and only the samples' reads update
the parameters (`STITCH/README.md:190-192`). If the panel has more than K haplotypes a haplotype
EM is run to compress them to K; if exactly K they are used directly; if fewer, the rest are filled
with noise (`README.md:196-200`; `STITCH/STITCH/R/reference.R:113-135, :157-162, :195-197`). Panel
alleles become θ = ε or 1 − ε with ε from `reference_phred` (20 by default).

## 5. STITCH without a panel

### 5.1 What an "ancestral haplotype" is

STITCH's parameter is a K × nSNPs matrix θ[k, j] = P(ALT | ancestral haplotype k, SNP j), plus for
each interval between grids a probability of *not* switching, σ, and a K-vector of where a switch
lands (`alphaMat`), plus a starting distribution. Nothing forces a row of θ to be a real
chromosome: it is the allele profile of the *k*-th hidden source that the reads of all samples are
best explained by. The paper describes "each chromosome in the population as a mosaic of K unknown
founders or ancestral haplotypes" (Davies 2016); the parameter docs call them "founder / mosaic
haplotypes" (`STITCH/STITCH/R/functions.R:4`). After EM, in a population with a small number of
real founders (the CFW mice, descended from eight inbred strains further reduced by a bottleneck),
the rows converge to those founders' haplotypes; in an outbred human population with K = 40 they
are a compression of the local haplotype diversity, and a sample's state posterior is a soft
assignment among them.

### 5.2 Initialisation

Without a panel, θ is drawn uniformly at random in (0, 1) per entry, `alphaMat` and the prior are
1/K, and σ between grids *g* and *g+1* is exp(−nGen · d_cM / 100) with a genetic map or
exp(−nGen · expRate · d_bp / 10⁸) without (expRate 0.5 cM/Mb)
(`STITCH/STITCH/R/functions.R:1350-1353`; `STITCH/STITCH/R/genetic-map.R:28-43`). `S` independent
random starts can be run and averaged (`S = 1` by default, `functions.R:5, :88`).

### 5.3 The E step

Per sample: per-read emissions (§3.2), grid emissions (§3.3), then forward–backward over the K²
pair states (`STITCH/STITCH/src/diploid.cpp:345-430` forward, `:432-500` backward). The transition
uses three numbers per interval — σ², σ(1−σ), (1−σ)² for "neither copy switches", "one switches",
"both switch" (`functions.R:4026-4033`) — and, when a copy switches, it lands on template *k* with
probability `alphaMat[k, g]` (`diploid.cpp:414-420`). The state posterior γ(k1, k2, g) is what the
M step consumes. (The expected-switch accumulator used for `alphaMat` and σ is
`rcpp_make_diploid_jUpdate`, `diploid.cpp:284-340`.)

### 5.4 The M step

Each read's bases are attributed to templates in proportion to the posterior (`diploid.cpp:532-618`):
for SNP *j* on read *r* and template k1, accumulate over k2 the posterior γ(k1, k2) times the
probability that the read came from k1 rather than k2 (a/(a+b) with a = eMatRead[k1, r],
b = eMatRead[k2, r], `:606`), and times the probability that the base implies ALT given θ
(d3 = pA·θ / (pA·θ + pR·(1−θ)), `:596-598`). Summed over all samples and reads, the ratio of the
"ALT-weighted" to the "total" accumulator is the new θ[k1, j] (`functions.R:4485`), clamped to
[emissionThreshold, 1 − emissionThreshold] = [1e-4, 1 − 1e-4] (`:4499-4500`); a SNP with no reads
on a template gets 0.5 or the pileup allele frequency (`:4487-4497`). `alphaMat` is the normalised
expected number of switch landings per template and interval (`diploid.cpp:284-340` accumulates,
`functions.R:4509-4526` normalises with floor `alphaMatThreshold` 1e-4). The switch rate is
re-estimated per interval as σ = exp(−(expected switches)/(2N)) and **bounded** between the rates
implied by nGen·minRate and nGen·maxRate (0.1 and 100 cM/Mb) (`functions.R:4462-4478`); the code
comment says "we estimate the compound parameter T · sigma_t, so we need to use nGen to bound
appropriately" (`:4473-4475`). So nGen is a prior scale and a bound, not a fixed rate.

### 5.5 Heuristics that make EM work

Random starts plus EM get stuck: two templates can be identical over a stretch, a rarely used
template can be empty, and templates can be "shuffled" (template 1 on the left continues as
template 2 on the right, which the transition model then pays for at every sample). STITCH has
three repairs, all on fixed iterations (defaults `functions.R:116-117, :131`):

- **shuffle** (iterations 4, 8, 12, 16; `functions.R:4242`): smooth the estimated switch rate over
  ±`shuffle_bin_radius` (5 kb) and, where it exceeds what nGen·maxRate allows, test whether
  swapping template labels on one side removes the excess (`STITCH/STITCH/R/heuristics.R:287-300,
  :396-420`);
- **refill** (iterations 6, 10, 14, 18; `functions.R:4289`): every 5 kb, templates whose average
  usage is below a threshold are re-seeded from other templates plus noise (`heuristics.R:1-30`);
- **split reads** (iteration 25): reads whose bases fit two templates better than one are split at
  the likely breakpoint (`functions.R:34, :117`).

### 5.6 The three `method`s

- `diploid` — the exact K² model above. The paper found it "computationally prohibitive beyond
  K = 30" (Davies 2016).
- `pseudoHaploid` — two *haploid* HMMs per sample, each read given a soft weight of belonging to
  copy 1 vs copy 2 (`pRgivenH1`, `pRgivenH2`, initialised at random in [0.2, 0.8],
  `functions.R:3088-3095`). The emission of copy 1 for read *r* becomes
  x·eMatRead + (1 − x)·pRgivenH2 with x = pRgivenH1/(pRgivenH1 + pRgivenH2)
  (`STITCH/STITCH/src/haploid.cpp:191-197`), i.e. a read that probably belongs to the other copy
  contributes a near-constant. The weights are re-estimated after each iteration from the two
  copies' posteriors (`functions.R:5221-5270`, bounded at 0.001). Linear in K; less accurate; the
  paper's human runs used 38 pseudo-haploid iterations then 2 diploid ones
  (`switchModelIteration`). The genotype posterior is formed as if the two copies were independent
  (`functions.R:3644-3654`).
- `diploid-inbred` — "is diploid, but fully inbred, so use statistical haploid model"
  (`functions.R:3289-3290`): one haploid HMM (`haploid.cpp:56-100`) with the 2-row transition
  (`functions.R:4034-4038`), and the output genotype is (1 − d, 0, d) from the haploid ALT dosage
  — heterozygotes cannot be called (`functions.R:3633-3643`). This is the right model for RILs and
  doubled haploids.

### 5.7 Choosing K and nGen

From the README (`STITCH/README.md:170-176`): K is "the number of ancestral haplotypes in the model";
larger K "allows for more accurate imputation for large samples and coverages, but takes longer and
accuracy may suffer with lower coverage"; try a few values and judge by validation or by the
distribution of INFO scores; keep K small enough that "each ancestral haplotype gets at least a
certain average X of coverage, say 10X, given your number of samples and average depth" — i.e.
N × depth / K ≳ 10. nGen: "It is probably fine to set it to 4·Ne/K"; if the population "can
reasonably [be] approximated as having been founded some number of generations ago", use that; the
method "should be fairly robust to misspecifications". The paper's choices: K = 4 for the CFW mice
(Supplementary Table 2), K = 40 for CONVERGE with "accuracy declined for K < 40, and was marginally
improved for K > 40" (Supplementary Table 3).

### 5.8 Small N, few real haplotypes: STITCH with K = 2 as a biparental imputer

This part is my reading of the code, not something the authors measured.

With K = 2 the `diploid` model has four states per grid: (1,1), (1,2), (2,1), (2,2) — exactly the
founder-genotype states of an F2 HMM. θ[1, ·] and θ[2, ·] converge to the two parents' haplotypes
(bounded away from 0/1 by 1e-4), the per-interval σ is the estimated probability that a chromosome
copy keeps its founder between grids — the recombination fraction, learned from the whole cohort
and bounded by nGen·(0.1 … 100 cM/Mb) — and `alphaMat` is redundant (a switch has one place to go).
The state posterior γ is the founder-origin posterior at every grid, which is the breakpoint map a
segregant analysis wants; STITCH will write it as the `HD` field if asked (§7). So yes: STITCH
with K = 2 is a biparental imputer that learns the parents from the offspring, with three
differences from a purpose-built one (LB-Impute, GBScleanR): it does not know that the two
founders differ at every site, it has no per-sample error rate, and its label-switching repairs are
tuned for K in the tens. For RILs or DH lines, `diploid-inbred` with K = 2 is a two-state haploid
HMM, the standard RIL genotype caller. For a backcross, one founder's haplotype is present in every
individual and the other in half of the copies; nothing in the model objects to that. For MAGIC, K
= number of founders.

Sample size: the README's 10× rule gives N ≳ 10·K/depth — at 1× and K = 2, about 20 individuals;
the paper found sample size "above 500 has little impact on performance" for the 0.15× mice with
K = 4 (Fig. 4), which is the same rule (500 × 0.15 / 4 ≈ 19 reads per template per site). With
few samples the risk is not the emission update but the *shuffle* repair: it relies on the
estimated switch rate exceeding nGen·maxRate·d, and with N = 20 that estimate is noisy.

## 6. Iteration structure

| tool | outer structure | inner structure | what is averaged for the output |
|---|---|---|---|
| GLIMPSE2 | 1 init + `--burnin` 5 + `--main` 15 (cap 15) iterations over all samples (`GLIMPSE/phase/src/caller/caller_initialise.cpp:51-53`; `caller_parameters.cpp:139-140`) | per sample: two conditional haploid forward–backwards, sample each copy, re-phase (§3.1) | genotype posteriors from HP0 × HP1 accumulated over main iterations only (`GLIMPSE/phase/src/objects/genotype.cpp:168-209`; sparse: stored only if P(0/0) < 0.99999) |
| QUILT | `nGibbsSamples` 7 independent chains × `n_seek_its` 3 seek iterations | per seek iteration: 20 Gibbs passes over reads with block moves at passes 3, 6, 9; then a whole-panel haploid pass per copy | dosage of the last seek iteration of each chain, averaged over chains (`quilt.R:18, :248-250`) |
| STITCH | `niterations` 40 EM rounds over all samples, with shuffle/refill/split at fixed rounds | per sample per round: diploid (K²) or pseudo-haploid forward–backward | the posterior of the final round only (`Options.md` under `switchModelIteration`) |
| Beagle 4.1 `gl=` | 10 burn-in iterations (4.0 sampling, `nsamples` 4 haplotype pairs per individual) + `niterations` 5 (4.1 phasing) | per window of 50,000 markers, overlap 3,000 | posterior genotype probabilities (`beagle_4.1_09Feb16.txt:144-147, :157-160`) |

GLIMPSE2 samples a haplotype pair per iteration (Gibbs); STITCH computes full posteriors and
re-estimates parameters (EM); QUILT mixes both (Gibbs on read labels, deterministic HMM given the
labels). None of the three carries the *cohort* through the iterations except STITCH — in GLIMPSE2
and QUILT the iterations of one sample never see another sample.

## 7. Outputs

### 7.1 Fields

| | GT | GP | DS | other per-sample | per-site |
|---|---|---|---|---|---|
| GLIMPSE2 | phased (from the last sampled pair) | yes | yes | — | `RAF` (panel ALT frequency), `AF`, `INFO` (`GLIMPSE/phase/src/io/genotype_writer.cpp:71-76`) |
| QUILT | phased when `output_gt_phased_genotypes` (default TRUE), else best guess masked below GP 0.9 (`quilt.R:22, :116`) | yes | yes | `HD`: two haploid dosages; `OHD` when truth labels given (`QUILT/QUILT/R/writers.R:64-69`) | `INFO_SCORE`, `EAF`, `HWE`, `ERC`/`EAC`/`PAF` |
| STITCH | best guess masked below 0.9; phased only with the experimental `do_phasing` (`STITCH/STITCH/R/writers.R:1002-1004`, `functions.R:79, :176-179`) | yes | yes | `HD`: K ancestral-haplotype dosages if `output_haplotype_dosages` (`writers.R:1016`, `functions.R:76`) | `INFO_SCORE`, `EAF`, `HWE`, `ERC`/`EAC`/`PAF`, `REF_PANEL` (`writers.R:1063-1076`) |
| Beagle 4.1 `gl=` | unphased | if `gprobs=true` | yes | — | — |

GLIMPSE2's ligation step keeps only `GT` (plus AC/AN) across chunk overlaps and drops DS/GP/INFO
(`GLIMPSE/ligate/src/ligater/ligater_algorithm.cpp:266-273`); users who want dosages keep the
per-chunk files.

### 7.2 INFO score without truth

All three sequencing tools compute the IMPUTE measure from the per-sample posteriors, no truth
needed. With e_i = GP1 + 2·GP2 the dosage of sample *i*, f_i = GP1 + 4·GP2, N samples and θ̂ =
Σe_i / 2N the estimated ALT frequency:

```
INFO = 1 − Σ_i (f_i − e_i²) / (2N · θ̂ · (1 − θ̂))
```

STITCH: `STITCH/STITCH/R/writers.R:529-533` accumulates e and f − e², `:225-227` divides; the header
says "same as IMPUTE info measure I_A" (`:1091`). QUILT: `QUILT/QUILT/R/writers.R:39-56`, identical,
with INFO set to 1 when θ̂ rounds to 0 or 1. GLIMPSE2: `genotype_writer.cpp:174-187`, the same
sum written as 1 − Σ(GP1 + 4·GP2 − DS²)/(2N·f(1−f)), clipped at 0. The numerator is the summed
posterior variance of the dosage, so INFO is 1 minus (uncertainty left) / (uncertainty a
Hardy–Weinberg site has before any data). STITCH's own pass filter is INFO > 0.4 and HWE p >
1e-6 (`functions.R:719`); the paper's post-QC numbers (§9) use that filter.

### 7.3 Phase

QUILT and GLIMPSE2 phase as a by-product (the read labels; the sampled pair). Published switch
error for QUILT at 1× on NA12878: 0.13 % Illumina, 0.09 % ONT, 0.08 % haplotagged; 0.11 % at
0.25× haplotagged (Davies 2021, Fig. 2, Supplementary Table 1). GLIMPSE1's preprint ligates chunks
by minimising Hamming distance over overlaps using up to 16 stored haplotype-pair samples per
variant; GLIMPSE2's ligate flips a downstream chunk when mismatching heterozygous configurations
outnumber matching ones (`ligater_algorithm.cpp:108-120, :194-195`). STITCH's phasing is marked
experimental and disabled from the CLI.

### 7.4 Indels and multi-allelic sites

- GLIMPSE2: bi-allelic only (multi-allelic records skipped, `ref_genotype_reader.cpp:159, :225`).
  Indels are carried by the panel but excluded from the PBWT, the state selection and the HMM
  recursion (`ref_haplotype_set.cpp:77-81`; `haplotype_set.cpp:365-369`;
  `imputation_hmm.cpp:88, :195`); they are imputed "into a haplotype scaffold obtained from
  high-quality SNPs" (Rubinacci 2023, Methods), i.e. by reading off the copied templates, and their
  own GLs are ignored unless `--use-gl-indels`. Accuracy at indels is in Extended Data Fig. 5a of
  that paper.
- QUILT: bi-allelic SNPs only, everything else dropped when the panel is prepared
  (`quilt-prepare-reference.R:246`); the read extractor keeps only bases equal to the single REF
  or ALT base (`bam_access.cpp:329`).
- STITCH: "STITCH only handles bi-allelic SNPs" (`STITCH/Options.md`, `posfile`); validators refuse
  otherwise (`STITCH/STITCH/R/validators.R:174, :183`).
- Beagle 4.1: any VCF record with a GL, but the phasing model is over alleles as symbols; the
  manual does not restrict to SNPs.

## 8. Cost

None of the three repositories ships a complexity table; what follows combines the code structure
with the papers' measurements.

**GLIMPSE2.** Work per sample per iteration is (selected states K ≈ Kpbwt) × (sites where a
selected state carries the minor allele, L_poly) (`conditioning_set.cpp:137-143`); the forward
table is L_poly × K floats per thread (`imputation_hmm.cpp:55-57`). Total ≈ 21 iterations × N × K
× L_poly, plus PBWT queries that are independent of panel size. The panel is loaded once per chunk
from the binary file and is the memory floor. Measured: whole-genome imputation of one 1× genome on
the UK Biobank RAP cost £0.08 with GLIMPSE2 against £1.11 for GLIMPSE1 and £242.80 for QUILT
v1.0.4 (Rubinacci 2023); the method "scales sublinearly in both the number of samples and markers"
(Abstract) and was run on simulated panels of 2 million haplotypes (Supplementary Note 2,
Supplementary Fig. 7). QUILT2's authors measured GLIMPSE2 at about £0.09 per sample for any
coverage with N = 100 samples per run, and note its per-sample cost only flattens with respect to
panel size "for larger number of input samples (N > 100)" (Li 2025, Fig. 3). Chunks: `GLIMPSE2_chunk`
ends a window when it exceeds 2.5 cM, 2 Mb and 20,000 common variants, with a buffer of 0.5 cM,
0.4 Mb and 2,000 (`GLIMPSE/chunk/src/chunker/chunker_parameters.cpp:43-56`).

**QUILT.** Gibbs step: reads × Ksubset per pass × 20 passes × 3 seek iterations × 7 chains. Full
step (QUILT1): grids × (distinct patterns per grid, ≤ nMaxDH) per copy per seek iteration, which
is what makes QUILT1 "approximately linear computational complexity in reference panel size"
(Davies 2021, Supplementary Fig. 2); QUILT2's msPBWT makes the selection "computational time
independent to the panel size" (Li 2025). Memory: `rhb_t` is one bit per haplotype per SNP (4
bytes per 32 SNPs) plus the distinct-pattern tables; the preprint reports "a maximum of 20 GBs of
RAM for N = 400,000 haplotypes" against about 160 GB for QUILT1, and the paper "less than 10 GBs"
for both QUILT2 and GLIMPSE2 at about 1 million haplotypes on chromosome 20, RAM rising linearly
with panel size in both (Li 2025, Fig. 3a). Cost on the UKB RAP with the 200k panel, N = 32 per
run: £0.237, £0.299, £0.371, £0.475 per sample at 0.1×, 0.5×, 1.0×, 2.0× (Fig. 3b) — run time
grows with coverage because the Gibbs step is per read. Chunking is by `quilt_chunk_map` (default
≥ 3 Mb and ≥ 4 cM per chunk, `QUILT/QUILT/R/functions.R:3294`) with a `buffer`; samples are
independent so N only multiplies time.

**STITCH.** `diploid`: 40 rounds × N × grids × K² for the forward–backward, with α and β tables of
K² × grids doubles per sample in flight; `pseudoHaploid`: × K instead of K². All samples must be
read every round, so the read records are kept on disk in bundles (`inputBundleBlockSize`) or in
RAM (`keepSampleReadsInRAM`); the EM accumulators are K × nSNPs. Gridding (`gridWindowSize`) bins
reads into windows to cut the grid count at very low coverage (`functions.R:71`). Measured: on the
CFW mice STITCH was 38× slower than findhap and about 3× slower than Beagle 4 without a panel, but
"imputation for all samples for any of the methods could be performed in less than 48 hours on a
modest computational server" if parallelised by chromosome; on CONVERGE, Beagle with the 1000G
panel took 7.3× longer than STITCH (Davies 2016, Supplementary Tables 2, 5–7).

**Beagle 4.1 `gl=`.** The haplotype-cluster model is rebuilt from all samples every iteration; the
GLIMPSE1 preprint measured Beagle 4.1 with the HRC panel at about 1200× the run time of GLIMPSE
and quadratic in panel size (Fig. 2D), and GeneImp's authors measured Beagle 4.0 at 135.6 h for
chromosome 22 of 16 samples (Table 5). This is why the mode was dropped and why nobody uses it
at scale.

**Reveel**: linear in N; 81 min for 1000 simulated samples against 6,120 min for Thunder.

## 9. Published accuracy

All numbers are r² between imputed dosage and truth unless stated; "rare" and "common" are the
authors' MAF bins.

### 9.1 With a reference panel

**GLIMPSE1** (Rubinacci 2020 preprint, European samples, HRC / 1000G panels): at 0.3× r² > 0.9 for
MAF > 5 %; at 1× r² = 0.8 at MAF 0.1 %; at 8× r² > 0.95 at MAF 0.1 % (Fig. 2A, 0.1×–8×). On
chromosome 1 at 1× GLIMPSE beat Beagle 4.1, the second-best method, by about 0.2 in r² at MAF
0.1 % (Fig. 2C); ~6.5 CPU hours per 1× genome.

**GLIMPSE2** (Rubinacci 2023): chromosome 20, 100 UK Biobank British samples at 1.0×, UKB panel of
280,238 haplotypes — r² 0.892 at MAF 0.01 % (GLIMPSE1 0.561; QUILT v1.0.4 0.925) and 0.927 at MAF
0.1 % (QUILT 0.927) (Fig. 1a); coverage 0.1×–4× on 10,000 UKB British samples, chromosome 1
(Fig. 1c); non-European samples from 1000 Genomes imputed with the UKB panel "similar" to QUILT
(Supplementary Note 2, Supplementary Fig. 4); indels in Extended Data Fig. 5a.

**QUILT** (Davies 2021): NA12878 at 0.5× — rare (0.1–0.2 %) r² 0.678 Illumina, 0.709 haplotagged,
0.669 ONT, GLIMPSE Illumina 0.672; common (20–50 %) 0.975 / 0.980 / 0.937 (Fig. 3A, Supplementary
Table 2). At 0.1× Illumina rare 0.416; at 1.0× rare 0.754, common 0.989. "For lower-coverage data
(0.1X) QUILT is favoured over GLIMPSE"; "for high-coverage data (1.0, 2.0X) GLIMPSE is favoured
over QUILT for the Illumina and HT data types". 59 British 1000G samples at 0.25× (Fig. 4B): rare
0.593 vs GLIMPSE 0.519, common 0.931 vs 0.917. Five 1000G populations with the HRC panel minus the
test samples (Fig. 5): CHB rare r² 0.581 at 0.25×, 0.768 at 1×, 0.840 at 2×, against 0.485 (GSA
array) and 0.43 (UKBB array).

**QUILT2** (Li 2025 and preprint): CEU 1000G samples, chromosome 20, panels 1000G (5,008
haplotypes), HRC (54,330), UKB 200k (400,022) (Fig. 2a): common variants r² above 0.95 with every
panel; very rare (0.01–0.02 %) r² about 0.50–0.65 at the lowest coverages rising past 0.80 at 2×
with the UKB panel; QUILT2 more accurate than GLIMPSE2 at 0.25–0.5×, GLIMPSE2 more accurate at
rare variants at 2×. ONT 1× common r² 0.937 vs GLIMPSE2 0.695 (Fig. 2c). Ancient DNA (Afanasievo)
QUILT2 ahead with the 1000G and HRC panels, equal with UKB (Fig. 2b).

**loimpute** and **GeneImp**: §3.5.

### 9.2 Without a panel — STITCH, by N and K (Davies 2016)

- CFW outbred mice, 2,073 animals at 0.15×, K = 4, nGen ≈ 100: r² 0.972 against the MegaMUGA array
  (21,576 SNPs) before filtering, 0.981 after INFO > 0.4 and HWE p > 1e-6 with 81 % of 5.72 M SNPs
  retained; against 10× sequencing of 4 mice 0.948 / 0.974 (Fig. 2a–b). Beagle 4 without a panel
  reached 0.219 (array) and 0.080 (10×), findhap 0.55 / 0.58. Sample size "above 500 has little
  impact" and 0.06× with all 2,073 is "only marginally poorer" than 0.15× (Fig. 4).
- CONVERGE Han Chinese, 11,670 people at 1.7×, first 10 Mb of chromosome 20, K = 40: r² 0.920
  against the array before filtering, 0.939 after; 0.949 / 0.960 against 10× sequencing of 9 people
  (Fig. 3A–B). Beagle 4 without a panel: 0.930 (array) and 0.886 (10×). Beagle 4 *with* the 1000G
  panel: 0.943 (array) — better than STITCH's 0.922 pre-filter but 7.3× slower (Fig. 3C–D,
  Supplementary Tables 5–6). Accuracy "consistently improve[s] with increasing sequencing depth"
  over 0.3–1.7× and depends less on N than in the mice (Fig. 4); "accuracy declined for K < 40, and
  was marginally improved for K > 40" (Supplementary Table 3).
- Removing read structure: 0.97 → 0.87 in the mice (Supplementary Table 2).

The two cohorts sit at opposite ends of the haplotype-count axis: the mice are a few founders at
very low depth (STITCH's home ground, where Beagle collapses), the humans are outbred at 1.7×
(where STITCH matches Beagle-without-panel and is 0.02 behind Beagle-with-panel). Nothing is
published for STITCH at N in the tens.

## 10. The continuum: one machine, different template sets?

Yes at the level of the equations, no at the level of the code, and the difference is exactly the
template set.

**What would make GLIMPSE2 run without a panel.** The recursion is a haploid Li and Stephens over
an arbitrary haplotype set; the panel-specific parts are the mandatory checks and the *static*
data structures. Concretely (from the earlier reading, `glimpse2_code_facts.md` §9, spot-checked):
the gate at `caller_parameters.cpp:133-134`; the site list, allele counts, major alleles and
common/rare classes all come from the panel (`genotype_reader.cpp:265, :271, :374-376`); the
conditioning set dereferences `HvarRef`/`ShapRef` only (`conditioning_set.cpp:130, :157-168`) while
the targets' mirror structures `HvarTar`/`ShapTar` already exist and are rebuilt every iteration
(`haplotype_set.cpp:171-198`); the PBWT is pre-built over the reference and run-length compressed
(`ref_haplotype_set.cpp:47-103`), indexed by `n_ref_haps`. A panel-free GLIMPSE2 would have to (a)
let the conditioning set point at target haplotypes and exclude a sample's own pair, (b) rebuild
the PBWT over the targets every iteration because they change — this is the largest change and
removes the pre-compression that is GLIMPSE2's main speed trick, (c) seed iteration 0 without
rare-allele sharing against a panel, (d) estimate allele frequencies and the common/rare split from
the cohort each iteration, and (e) size `nrho` by the total haplotype count. That is GLIMPSE1
(§4.1), whose preprint did all of this with K = 1000 and the targets in the PBWT.

**What would make STITCH use a panel.** It already does, in the weakest way: as an initialiser
(§4.3). To use a panel as GLIMPSE does, θ rows would have to be *fixed* panel haplotypes (θ ∈
{ε, 1−ε}) rather than parameters, the M step skipped for those rows, and K would have to grow to
hundreds — at which point the K² diploid recursion is unaffordable and the pseudo-haploid or a
read-label scheme is needed. That is QUILT: the same read record, the same per-read emission, θ
frozen at panel alleles, K² avoided by read labels, and the template subset chosen by matching
instead of learned. QUILT is STITCH with a panel; the shared `copied-from-stitch` sources say as
much.

So the continuum has three points that all share one emission and one transition:

| templates | who learns them | selection per sample | tool |
|---|---|---|---|
| K learned rows of θ | EM over the cohort | none (all K) | STITCH |
| fixed panel rows | nobody | random start → matching (top-K or msPBWT) | QUILT |
| fixed panel rows + other targets' current pairs | Gibbs over the cohort (targets only) | PBWT neighbours | GLIMPSE1; GLIMPSE2 drops the targets |

A design for ng can keep the HMM core (haploid recursion, per-read or per-site emission, genetic-
distance transition) fixed and vary only the template source: learned rows for the few-haplotype
case, the phased deep samples for the many-haplotype case, and an external panel when one exists.

## 11. What ng would hand to each tool and what it gets back

ng has, per locus: per-sample genotype likelihoods over the locus's alleles, and per read (chain
id) the allele it showed. That maps onto the tools' inputs as follows.

**Today, through files.**

- **GLIMPSE2**: a VCF with `PL` at bi-allelic SNPs is enough (`--input-gl`); indel PLs are used
  only with `--use-gl-indels`, multi-allelic records are dropped. It needs a panel, so it serves
  only the many-haplotype case *with* an external panel — or with a panel ng makes by phasing its
  own deep samples (SHAPEIT5 on the 20–40× individuals) and handing that to `GLIMPSE2_split_reference`.
  That "deep samples as the panel" route is the closest thing to the brief's cohort-as-own-panel
  case that runs today, and it is the natural pipeline: ng calls, phases the deep set, imputes the
  shallow set. Returns GP, DS, phased GT, INFO, per chunk.
- **QUILT / QUILT2**: consumes BAM/CRAM plus a `posfile`; it would re-read the alignments ng already
  read. Panel mandatory (same deep-samples route). Returns GP, DS, phased GT, HD, INFO_SCORE,
  bi-allelic SNPs only. Its Gibbs step could in principle take ng's chain-id allele lists directly —
  the `sampleReads` record is (SNP index, signed base quality) per read, which ng can write — but
  only through R internals, not a CLI.
- **STITCH**: BAM/CRAM plus a `posfile`, no panel, bi-allelic SNPs only, `K` and `nGen` chosen by
  the user. This is the only tool that runs the few-haplotype case and the all-shallow case
  without any panel. Returns GP, DS, `HD` (the founder-origin dosages — the breakpoint map for a
  segregant population), INFO_SCORE; phasing experimental.
- **Beagle 4.1 `gl=`**: a VCF with `GL`/`PL`, any ploidy of allele symbols, with or without a panel;
  returns GP/DS, unphased; too slow to be a production route (§8).

**Where an in-caller module has the advantage.**

1. *Reads never re-read.* QUILT and STITCH spend their input stage re-pileuping BAMs to build
   exactly the per-read allele list ng already holds. An in-caller module starts from that list.
2. *Per-read evidence, not base quality.* All per-read emissions here are base quality → ε, with a
   uniform ε/3 for the wrong base. ng has a per-read allele likelihood (from its own alignment
   model); the emission `θ·pA + (1−θ)·pR` takes any (pR, pA) pair, so ng's numbers slot in.
3. *Multi-allelic sites, indels and STRs.* Every tool is bi-allelic-SNP-only at the model level
   (§7.4), because the panel bit matrices, the PBWT symbols and the 2-vector haploid likelihoods
   all assume two alleles. The emission generalises to Σ_a θ[k, a]·P(read | a) with a per-template
   allele-probability vector, and the transition does not care. The msPBWT in QUILT2 is explicitly
   multi-symbol (Li 2025 notes it "can also work with non-binary symbols, such as multi-allelic"),
   so selection over multi-allelic templates is not a research problem. STR loci with ng's
   length-keyed alleles would fit the same slot.
4. *Cohort-level learning with the reads in hand.* STITCH's M step needs, per template and site,
   the ALT-weighted and total read attributions from every sample; ng's locus-by-locus merge is
   where those sums are formed anyway.
5. *Depth range.* All three tools cap per-site depth (STITCH 50, QUILT 30, GLIMPSE2 `--max-depth`
   40 for its own pileup, `GLIMPSE/phase/src/caller/caller_parameters.cpp:83`) to avoid numerical
   trouble; ng's deep samples would be down-sampled by them. An in-caller module
   can keep the deep samples' GLs exact and use them as near-fixed templates.

What ng gets back from such a module is what the tools return: GP and DS per sample per site, a
phased pair per sample (or, for segregants, the founder-origin posterior along the chromosome),
and an INFO score per site — all three tools compute the same INFO, so adopting it costs nothing.

## 12. Questions to discuss

1. **One HMM core, two template sources?** Recommendation: yes — implement the haploid Li and
   Stephens recursion once, with the template set as an abstraction that is either (a) K learned
   allele-probability rows updated by EM (STITCH mode) or (b) fixed phased haplotypes with per-sample
   selection (QUILT/GLIMPSE mode). §10 shows the equations are shared; the cost is that mode (a)
   needs the K² diploid recursion or a read-label scheme, and mode (b) needs a matching index. Trade-
   off: two code paths to test against two sets of external tools, versus the reuse of everything else.
2. **Per-read or per-site evidence?** Recommendation: per-read (chain ids) with read labels, as
   QUILT does, because ng has it for free and STITCH's own measurement (0.97 → 0.87 without reads in
   the mice) says it matters where SNPs are dense — the crop and segregant cases. Trade-off: the
   Gibbs step over reads costs time proportional to read count, and its value is nil where reads
   rarely span two heterozygous sites (humans at 1.7× made no difference).
3. **Few-haplotype populations first.** Recommendation: build mode (a) with small K (2 for
   biparental, founders for MAGIC) and `diploid-inbred` for RILs/DH, since K² is trivial there and
   the `HD` founder-origin output is the deliverable those users want; compare against LB-Impute
   and GBScleanR on the same data. Trade-off: STITCH's label-switching repairs are the fragile part
   at small N and would need their own validation; ng could instead pin θ rows to the two parents'
   genotypes when those are sequenced, which removes the problem.
4. **Many-haplotype populations without an external panel.** Recommendation: phase the deep samples
   (a job for the phasing sub-report's tools) and use them as a fixed panel for the shallow ones —
   the GLIMPSE2/QUILT route — rather than building a GLIMPSE1-style cohort Gibbs sampler first.
   Trade-off: the shallow samples' own haplotypes then never help each other; GLIMPSE1 and Beagle
   4.1 show that helps when the deep set is small, and it could be added later as mode (a) with
   large K in pseudo-haploid form.
5. **Multi-allelic templates from day one?** Recommendation: design the template representation as
   an allele-probability vector per site (not a bit), even if the first implementation only ever
   sees two alleles; the four tools' bi-allelic restriction is a data-structure choice made early
   that none of them undid. Trade-off: memory per template site goes from 1 bit to at least a byte
   per allele; the deep-sample panel for a 1000-sample cohort at 10 M sites would be 2000 × 10 M ×
   (alleles) bytes, so a sparse encoding (GLIMPSE2's rare-carrier lists) is needed from the start.
6. **Quality score.** Recommendation: emit the IMPUTE INFO score exactly as §7.2, so results are
   directly comparable with all three tools' filters (STITCH's own pass rule is INFO > 0.4).
7. **Validation protocol.** Recommendation: down-sample the deep samples of each benchmark cohort to
   0.5×, 1×, 2× and 4×, impute with the module and with STITCH / GLIMPSE2 (deep samples as panel),
   report r² by MAF bin and switch error — the protocol every paper here uses, which makes the
   numbers in §9 the yardstick.

## 13. What could not be verified

- The full texts of Browning & Browning 2016 (AJHG, Beagle 4.1) and Browning & Yu 2009 (AJHG) are
  behind logins (Cell/ScienceDirect 403, PMC reCAPTCHA, EuropePMC no full text). §3.4's description
  of the Beagle 4 GL sampling rests on the 4.0 and 4.1 manuals; the attribution of the mechanism to
  Browning & Yu 2009 is from that paper's title and abstract only.
- The GLIMPSE1 journal version (Nat Genet 2021) redirects to a login; all GLIMPSE1 facts are from
  the bioRxiv preprint, and the "European samples, HRC/1000G" attribution of the Fig. 2A coverage
  curve is as summarised from the preprint, not re-read in the figure legend.
- The QUILT2 volume number returned by the full-text service ("volume 17, 2025") is inconsistent
  (Nat Commun volume 17 is 2026); the DOI 10.1038/s41467-025-67218-1 and PubMed id 41390671 are the
  reliable identifiers.
- The Genome Research 2025 paper on loimpute (doi:10.1101/gr.280175.124) is behind a login; only
  the 2021 BMC Genomics supplementary note describes the model.
- STITCH's supplement (where K and nGen choice is discussed at length, per the README) was not
  opened; the README text is quoted instead.
- The QUILT defaults differ between the 2021 paper (Ksubset 400) and the 2.0.4 code (600); the
  code value is cited.
