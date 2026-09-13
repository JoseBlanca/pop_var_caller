# Decision: phasing and imputation stay outside the variant caller — and the route we take instead

*Decision record, 2026-09-13. Owner's ruling after the research in this directory
([`00_overview.md`](00_overview.md) and sub-reports 01–06). This document records the decision and its
reasons (§1), the algorithm we would want whether we run it through an existing tool or ever build
our own (§2), the existing tools that already implement it and what each needs from ng (§3), and the
measurement that decides whether they serve our case (§4). Every number carries its source in the
sub-reports; the arithmetic marked as such is this document's own.*

---

## 1. The decision and its reasons

**ng will not carry an imputation or phasing stage.** It will write what such a stage needs, and the
stage will be a separate program reading ng's VCF: an existing tool if the test in §4 passes, our own
only if a measured gap in §3 makes that necessary.

The reasons, in order of weight:

1. **The pass is a separate stage in any design.** Phasing and imputation need a window of many loci
   across every sample, whole chromosomes in practice, while ng builds and genotypes one cohort locus
   at a time and hands each record over in genome order as soon as it is called
   (`src/ng/run/psp_caller.rs:558-561`). The streaming spec already forbids coupling across in-flight
   segments and says any such coupling "must be a pass over the emitted records"
   (`doc/devel/ng/spec/run_streaming.md` §4.3). Inside or outside the binary, it is a second pass over
   the called cohort; keeping it outside costs nothing that integration would save.
2. **Everything the model needs travels in a VCF.** The per-sample genotype likelihoods as `PL`, the
   per-sample read counts per allele as `AD`, the cohort allele frequency as `AF`, and ng's
   multi-allelic loci exactly as ng calls them. Reading a VCF back is cheap next to reading alignments,
   which is what the tools that work from reads (QUILT, STITCH) spend their input stage on.
3. **Interchangeability.** With `PL` in the file, the same output feeds GLIMPSE2 today and any tool of
   ours later; either side can be replaced without touching the other. A caller that grows an
   imputation stage ties the two release cycles together and makes the caller the thing that has to
   change when the imputation model does.
4. **The two things that do not travel in a VCF do not justify integration today.**
   - *Per-read allele lists* (ng's chain ids), the raw material of read-backed phasing. At three reads
     a position with short reads, one adjacent heterozygous pair in five can be linked at all, and 70%
     of pairs are out of reach at any depth (sub-report 05 §2, own arithmetic checked against Delaneau
     et al. 2013). The published gain from read linkage inside a low-coverage LD model is 9 to 12%
     relative reduction of genotype discordance at about 4× on humans (HapSeq; sub-report 06 §6.4);
     nothing is measured on a plant cohort. If a census on our data ever shows it is worth having, ng
     can write a small per-sample sidecar of links between neighbouring heterozygous sites behind a
     flag, and a downstream tool can read it. The chain ids are consumed by the cohort merge today and
     do not reach calling (`src/ng/run/cohort_merge/build.rs:1473-1495`); retaining them is a
     contained change at the merge's per-sample fold (sub-report 05 §6.4).
   - *Re-calling genotypes with an LD-informed prior.* No published pipeline does it, because the
     imputed genotype already is that re-call: the 1000 Genomes design is "per-sample likelihoods →
     LD refinement across the cohort → the refined genotype is the call" (sub-report 06 §6.5).

**What ng owes the downstream stage**, and it is small:

- `PL` per sample per record, and `GP` on by default rather than opt-in. Every low-coverage tool takes
  likelihoods and none of the hard-call tools is usable below about 10× (sub-report 02 §10: a true
  heterozygote at three reads shows one allele one time in four, against an assumed mismatch rate of
  1 in 10,000 in SHAPEIT5).
- Likelihoods that are honest about the sample: damage-aware for ancient DNA, model-based at repeat
  tracts, so the downstream stage is not fed a flattened number.
- Records that a biallelic tool can split without losing the caller's meaning (`bcftools norm -m-`
  splits multi-allelic records; what happens to `PL` under the split has to be checked, §4.3).

---

## 2. The algorithm we want, and its statistics

Whether it runs inside GLIMPSE2 or in a tool of ours, this is the model. It is the haplotype-copying
model of Li and Stephens (2003), with the cohort as the set of template haplotypes, the genotype
likelihoods as the evidence, and two passes over a cohort of mixed depth. Sub-report 01 §2c–§3 has the
derivations and the primary papers; sub-report 03 §3–§6 has the code of each step in GLIMPSE2, QUILT
and STITCH.

### 2.1 The picture

One chromosome copy of one accession is, over any stretch of a few hundred kilobases, identical by
descent to some other chromosome of the cohort. So the whole copy is a patchwork of pieces copied from
other accessions' chromosomes, the boundaries being ancestral recombinations, and the rare differences
inside a piece being mutations since the common ancestor or sequencing errors. The model asks at every
site: *which other chromosome is this one copying here?* The imputed allele at an uncovered site is the
allele of the chromosome being copied; the phase is the assignment of alleles to the two copies; the
confidence of both is the confidence of the answer.

### 2.2 The two rules, with their formulas

**Switching (recombination).** Between neighbouring sites *d* centimorgans apart the copied template is
kept with probability 1 − *t* and replaced by a uniformly chosen other template with probability *t*,

    t = 1 − exp(−0.04 · Ne · d / K)

with *N*ₑ the effective population size and K the number of templates. The 0.04 is 4 per Morgan
divided by 100 cM. GLIMPSE2, SHAPEIT5 and Beagle 5.5 use this line verbatim
(`GLIMPSE/phase/src/containers/conditioning_set.cpp:34,190-194`;
`shapeit5/phase_common/src/objects/hmm_parameters.cpp:46`; `beagle_src/src/phase/PhaseData.java:59`);
GLIMPSE2's default *N*ₑ is 100,000 and SHAPEIT5's 15,000, and Beagle re-fits the rate from the data
in its burn-in rounds (sub-report 02 §2). With more templates the closest relative is closer, the
shared pieces are longer, and switches are rarer: K is in the denominator.

**Mismatching (mutation or error).** At a site the copied template's allele is what this copy shows,
except with probability *e*. GLIMPSE2's default *e* is 10⁻⁴, clamped to [10⁻¹², 10⁻³]
(`GLIMPSE/phase/src/caller/caller_initialise.cpp:152`).

### 2.3 The evidence: genotype likelihoods, collapsed to one chromosome copy

The evidence at a site is the genotype likelihood: the probability of the reads under each of the
three genotypes, which ng computes and will write as `PL`. The HMM scores one chromosome copy, so the
three-genotype likelihood is turned into a two-allele one using the current estimate of the other
copy. If the other copy currently carries the reference allele, then "this copy carries reference"
is the homozygous-reference genotype and "this copy carries alternate" is the heterozygote:

    other copy = 0:  HL(0) ∝ L(0/0),  HL(1) ∝ L(0/1)
    other copy = 1:  HL(0) ∝ L(0/1),  HL(1) ∝ L(1/1)

(`GLIMPSE/phase/src/objects/genotype.cpp:94-115`). The emission for a template carrying allele *a* is
then HL(*a*)·(1 − *e*) + HL(1 − *a*)·*e* (`GLIMPSE/phase/src/models/imputation_hmm.cpp:61-71`). A hard
genotype is the special case with one likelihood equal to 1. This one step is what makes the same
code serve a 40× accession, a 5× one and a 0.5× ancient sample.

### 2.4 The computation for one chromosome copy

Walk the sites left to right carrying one number per template, the probability that it is the one
copied here given everything to the left. At each site multiply each template's number by its emission,
then let a fraction *t* of the total leak uniformly into every template (the switch). Walk right to left
the same way. The product of the two walks at a site, normalised, is the posterior over templates given
the whole chromosome; summing it over the templates that carry allele 1 gives the probability that this
copy carries allele 1 there. Cost: K operations per site per copy. The forward step in GLIMPSE2 is

    α_k(l) = ( α_k(l−1) · (1 − t) / Σα(l−1)  +  t / K ) · emit_k(l)

(`GLIMPSE/phase/src/models/imputation_hmm.cpp:146-166`), with the posterior at
`imputation_hmm.cpp:336-351`. Sites where no selected template carries the minor allele are imputed
analytically without the HMM (`:355-365`).

### 2.5 Two copies, and the sweeps

Score copy one with copy two held fixed and draw copy one's alleles; score copy two with the new copy
one and draw it; a small diploid pass over runs of three heterozygous sites settles the phase between the
two (`GLIMPSE/phase/src/models/phasing_hmm.cpp`, sub-report 03 §3.1). That is one visit to one
accession. Because the templates are the other accessions' chromosomes, which are also being inferred,
the cohort is swept: each accession is re-scored against everyone else's current estimate and its
estimate replaced. GLIMPSE2 runs one seeding sweep, five burn-in and fifteen main sweeps, and averages
the posteriors of the main ones (`GLIMPSE/phase/src/caller/caller_initialise.cpp:51-53`,
`genotype.cpp:168-243`). The deep accessions are nearly right from the first sweep and anchor the rest.

### 2.6 Choosing the templates

Scoring all 2N chromosomes at every site is too slow above a few hundred samples: for 2,500 accessions,
20 million sites and 20 sweeps, about 10¹⁶ operations, roughly 120 core-days (own arithmetic). Every
modern tool keeps K to a few hundred with a *positional Burrows–Wheeler transform*: at each site, all
chromosomes sorted by the sequence they show ending at that site, so that the neighbours of a chromosome
in that order are the ones identical to it over the longest stretch to the left, the ones most likely
identical by descent with it there. GLIMPSE2 takes 12 neighbours on each side every 0.1 cM and pools
them, up to 2,000 templates per sample (`GLIMPSE/phase/src/caller/caller_parameters.cpp:67-73`,
`haplotype_set.cpp:317-394`); SHAPEIT5 and Beagle do the same with their own counts (sub-report 02
§3.3, §5.3). In a self-panel run the index is rebuilt every sweep because the templates move; over a
frozen panel it is built once and stored compressed, which is GLIMPSE2's main speed trick.

### 2.7 Two passes over a cohort of mixed depth

- **Pass A — the deep accessions phase and firm each other up.** Sweeps as in §2.5, templates from the
  cohort itself. With hard calls this is SHAPEIT5 or Beagle 5.5; with likelihoods it was Beagle 4.1's
  `gl=` mode, removed in Beagle 5 for its cost (about 1,200 times GLIMPSE's run time, sub-report 03
  §8), and GLIMPSE1's mixed mode. In inbred lines nearly every site is homozygous and phases itself;
  the work is in residual heterozygous blocks and in accessions at 5 to 10×.
- **Pass B — every shallower sample against the frozen result.** Pass A's haplotypes become a panel,
  the index is built once, and each shallow sample iterates only its own two copies against it
  (GLIMPSE2's per-sample loop). Samples are independent, so pass B parallelises and a new sample never
  reruns pass A. Its limit: a haplotype absent from pass A cannot be imputed, and a sample from such a
  lineage is pulled toward the panel's alleles with full confidence (ancient human genomes from
  populations under-represented in the 2,504-genome panel reached 18% error at heterozygous sites,
  Sousa da Mota et al. 2023, [PMC10282092](https://pmc.ncbi.nlm.nih.gov/articles/PMC10282092/)).
  Letting pass-B samples serve as each other's templates (GLIMPSE1's design) is the remedy, at the
  price of rebuilding the index each sweep.

### 2.8 Outputs and the statistics that score them

Per sample per site: the genotype posterior `GP` (three probabilities), the dosage
`DS = GP₁ + 2·GP₂` (expected alternate-allele count, 0 to 2, carrying the uncertainty), and a phased
`GT`. Per site, a truth-free quality score, the IMPUTE `INFO`:

    INFO = 1 − Σᵢ (GP₁ᵢ + 4·GP₂ᵢ − DSᵢ²) / ( 2N · f · (1 − f) ),   f = Σᵢ DSᵢ / 2N

clipped at 0; GLIMPSE2, STITCH and QUILT compute it identically
(`GLIMPSE/phase/src/io/genotype_writer.cpp:174-187`; sub-report 01 §5). It is the fraction of the
variance a perfectly imputed site would show that the imputed dosages do show. Two known biases: it does
not detect regional failure, and it is inflated at rare variants (sub-report 06 §1, §8), so it is a
filter, not a validation.

The statistics that validate against truth (sub-report 06 §1):

- **Aggregate r² by allele-frequency bin**: the squared correlation between imputed dosage and true
  genotype over all sites in a frequency bin, computed over held-out samples. The standard figure in
  every imputation paper; per-variant r² is unstable at low frequency, which is why the bin.
- **Non-reference concordance**: among sites where the truth is not homozygous reference, the fraction
  called correctly; it is not inflated by the sea of easy homozygous-reference sites.
- **Switch error rate** for phasing: the number of places where the phase flips relative to the truth,
  divided by the number of heterozygous sites minus one; needs a truth phase (a trio, or long reads).

---

## 3. Existing tools that already do this

All well maintained unless noted; every restriction below is checked in the cloned code
(`tmp/phasing_research/repos/`, sub-reports 02 and 03).

| tool | role for us | input it needs from ng | restrictions that bite | maintenance |
|---|---|---|---|---|
| **SHAPEIT5** `phase_common` (+ `phase_rare`) | pass A on the deep accessions from hard calls | biallelic VCF with `GT`, ≥ 50 samples, a genetic map or 1 cM/Mb | biallelic only; mismatch assumed 10⁻⁴, so unsuitable below about 10×; excludes a template whose heterozygotes agree with the target's over > 75% of a window (`phase_common/src/objects/compute_job.cpp:85`), a risk in inbred material | official GitHub repository disabled at the time of writing; bioconda package and binaries available |
| **Beagle 5.5** | pass A alternative | VCF with `GT`; multi-allelic records accepted; no minimum cohort | no likelihoods; excludes pairs sharing ≥ 2 cM of identical genotypes (`phase/Ibs2.java:42`) | actively maintained (Feb 2025 release) |
| **GLIMPSE2** | pass B: every shallow sample against pass A's haplotypes | biallelic VCF with `PL` (`--input-gl`); the phased pass-A VCF as panel through `GLIMPSE2_split_reference` | biallelic SNPs; multi-allelic records dropped; indel likelihoods zeroed unless `--use-gl-indels`, and indels imputed only as passengers on the SNP scaffold; panel mandatory; the other targets are never templates | actively maintained; the ancient-DNA field's standard |
| **GLIMPSE1** | pass B when shallow samples must also be each other's templates | same as GLIMPSE2 | biallelic; slower; superseded by GLIMPSE2 but still runs | superseded, not removed |
| **QUILT2** | pass B from alignments rather than likelihoods | BAM/CRAM plus a site list; panel mandatory | biallelic SNPs; re-reads alignments; per-read Gibbs step scales with read count | actively maintained |
| **STITCH** | the no-panel route, treating every sample alike | BAM/CRAM plus a site list; K and generations chosen by hand | biallelic SNPs; caps depth at 50 reads; needs about 500 samples for r² 0.95 and gains nothing below 100 (sub-report 06 §4.1); R + C, slow | maintained |
| **Beagle 4.1** `gl=` | historical: the one tool that ran the whole cohort as its own panel from likelihoods | VCF with `GL`/`PL` | about 1,200 times slower than GLIMPSE; unphased output | unmaintained since 2016 |

**The route to try first**, because both tools are maintained, used at biobank scale, and need nothing
from ng but `PL`:

    ng call-from-psps                  →  cohort.vcf.gz with PL and GP
    bcftools norm -m- ; keep biallelic →  cohort.biallelic.vcf.gz
    SHAPEIT5 phase_common (or Beagle)  →  deep.phased.vcf.gz          (pass A, the ≥ 20× accessions)
    GLIMPSE2_chunk / split_reference   →  the panel, one binary file per chunk
    GLIMPSE2_phase + ligate            →  shallow.imputed.vcf.gz      (pass B, everyone below 20×)

The 500 accessions at 5 to 10× belong in pass B, not pass A: SHAPEIT5's assumed mismatch rate is far
below what a 5× hard call delivers (a true heterozygote shows one allele about 6 times in 100 at 5×,
own arithmetic), and GLIMPSE2 will use their likelihoods.

**What that route cannot give us**, whatever the test says: ng's multi-allelic records, indels
imputed on their own evidence, and repeat-tract genotypes. Whether that loss matters is a count (§4.3),
not an opinion.

---

## 4. How to test whether the existing route serves our case

The protocol every paper in sub-report 06 uses, applied to our data. Each step names what it decides.

### 4.1 Prerequisite in ng

Write `PL` per sample per record and turn `GP` on by default. Check that `bcftools norm -m-` splits ng's
multi-allelic `PL` vectors correctly (the split re-indexes genotypes; a wrong re-indexing silently
corrupts the likelihoods at exactly the records that matter).

### 4.2 Held-out downsampling: does pass B reach the published accuracy?

1. Choose about 50 deep accessions spanning the collection's structure (domesticated, wild relatives,
   any group the ancient samples might resemble). Keep their full-depth genotypes as truth.
2. Downsample their alignments to 0.5×, 1×, 2× and 5×; run ng on the downsampled files to get `PL`.
3. Run pass A on the remaining deep accessions only, so the held-out ones are not in the panel.
4. Run pass B on the downsampled held-out samples.
5. Score aggregate r² by allele-frequency bin (bins 0.1–1%, 1–5%, 5–50%) and non-reference concordance
   against the full-depth truth, at each downsampled depth.

*What decides.* Published yardsticks for internal panels in low-diversity populations are r² about 0.9
at 0.5 to 1× with 50 to 150 deep individuals and 0.95 with 200 to 500 (barn owl, strawberry, cattle
series; sub-report 06 §4.3); with 2,000 deep accessions we should sit at or above the upper figure for
alleles above 5% frequency. A result far below it at 1× is a failure of the route on our material, and
the first suspect is §4.4.

### 4.3 What the biallelic restriction costs: a count

From ng's VCF, count records dropped by the biallelic filter and the reads they carry, split by SNP,
indel and repeat tract; count how many of the dropped records have a cohort allele frequency above 5%.
*What decides.* If the dropped set is a few percent of records and mostly rare, the restriction is
tolerable; if it is a large share of common variation, it is the one gap that argues for a tool of our
own, and it is a gap in the tools' data structures, not in the model (sub-report 03 §7.4, §11).

### 4.4 Inbreeding: are the template sets being stripped?

SHAPEIT5 and Beagle exclude templates that look like the target on both chromosomes over a long
stretch. With median inbreeding about 0.78 across the collection, many accession pairs will. Two checks:

- Run pass A once with the exclusion at its default and once with it relaxed (SHAPEIT5 has no flag for
  it; Beagle's `ibs2` threshold is a parameter in `phase/Ibs2.java`, and the run's own log reports the
  number of conditioning states per window), and compare pass-B accuracy from §4.2 under the two.
- Score pass A's own phasing on the one pair of psps that are one plant sequenced twice: the two runs
  should agree on phase inside every residual heterozygous block; disagreement measures the phasing
  noise in inbred material where no trio truth exists.

*What decides.* If accuracy is unchanged, the rule is harmless here. If relaxing it helps, we have
found the one place where our material breaks the tools' assumptions, and it is a small, contained
change in Beagle's source or a reason to prefer Beagle over SHAPEIT5 for pass A.

### 4.5 Lineage absence: the ancient-sample risk, measured on modern data

Hold out an entire group of deep accessions (a wild species, or the most divergent cluster), downsample
them, and impute them against a pass-A panel that lacks their group. Compare with §4.2's numbers for
the same samples imputed against a panel that includes their group. *What decides.* The gap is the
error an ancient sample from an unrepresented lineage will suffer, and it says whether GLIMPSE1's
mode (shallow samples as each other's templates) has to be part of the route for the ancient set.

### 4.6 Cost

Wall time and peak memory of pass A at 2,000 accessions and of pass B per sample, on the machine that
will run them. GLIMPSE2's published cost is about £0.08 per 1× genome on the UK Biobank platform and
the per-sample cost flattens with respect to panel size only above about 100 samples per run
(sub-report 03 §8); pass A's cost is the PBWT rebuild per sweep, proportional to haplotypes × sites.

### 4.7 Acceptance

The route is adopted if §4.2 meets the yardstick at 1× for alleles above 5% frequency, §4.3's loss is
mostly rare variation, and §4.4 shows no stripping effect. Any one failing names the part to build, and
only that part: a multi-allelic template representation for §4.3, a per-sample inbreeding-aware
template rule for §4.4, the self-panel loop over shallow samples for §4.5. The engine of §2 is not
rebuilt for any of them.
