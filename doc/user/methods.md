# Methods

## 1. Overview

pvc is a population caller: it calls SNPs, indels and STRs in all the samples of a cohort at
once.

![How evidence flows through pvc: a per-sample step that writes one psp file per sample, cohort steps that estimate the parameters and call jointly, and an optional refinement](figures/methods_signal_flow.svg)

*The lane "Per-sample pileup" is done once per sample and never again. The lane "Cohort calling"
is done over the whole cohort, and again when a sample is added. The lane
"Refinement/imputation" is optional. [Calling](calling.md) compares this with direct mode.*

The steps, in order:

1. **Per sample: gather the evidence.** The reads are filtered and corrected, and what they show
   at every locus is written to a per-sample pileup file (a *psp file*). No genotype is called
   and no likelihood is computed here (sections 2 and 3).
2. **Per cohort: estimate the parameters** used as priors in calling: how often each library's
   bases are wrong, how inbred each sample is, how variable the population is, how often STRs
   stutter. These priors are estimated from the cohort's own reads, without calling a genotype
   first (sections 6 and 7). This step is optional: without it, calling uses the built-in
   values of [Estimated parameters](parameters.md#the-built-in-values---defaults), with F = 0
   for every sample.
3. **Per cohort: call.** The genome is divided into region types, each with its own likelihood
   model of how a read arises from an allele (section 4). At STRs pvc calls one locus whose
   alleles are whole repeat sequences: this is the *STR path*. Everywhere else it calls SNPs and
   indels: the *generic path*. At each locus the candidate alleles are
   chosen from all samples (section 5), and a loop alternates between the allele frequencies and
   the samples' genotypes until they agree (section 8). Qualities and site filters follow
   (sections 9 and 10).
4. **Optionally, when asked for: refine.** This is what is often called imputation; these pages
   call it *refinement*. Each sample's genotypes are called again (section 11)
   and phased again (section 12), using the information from the haplotypes that the other
   samples carry along the chromosome.

pvc stores the pileup information in one file per sample, so it does not have to read the
alignment files again each time the cohort is called. GATK's gVCFs serve the same purpose. But
GATK genotypes each sample from its own reads alone. Its joint step only makes the samples'
alleles compatible, so they can be written in one VCF. In pvc, each sample's genotype also uses
what the other samples carry (section 8).

These are the main methods used by pvc:

- Standard read filters, indel left-alignment and read realignment at STRs; no local
  reassembly (section 2).
- Different likelihood models for generic loci (SNPs and indels), STRs and STR bundles;
  satellites are not called by default (section 4).
- Candidate alleles are chosen from all the samples together (section 5).
- Priors are estimated from the cohort's own reads by marginal likelihood, without calling
  genotypes first (section 6).
- Each sample's inbreeding coefficient F is part of the genotype prior, for inbred and selfing
  species (section 7).
- Bayesian genotype calling with an EM loop between allele frequencies and genotypes
  (section 8).
- Site quality with penalties for sequencing artifacts, and a filter for hidden duplications
  (sections 9 and 10).
- Optional refinement from the haplotypes the samples share (section 11).
- Phase from read pairs, extended by the refinement (section 12).

## 2. Read filters and read fixes

### Filters

pvc ignores duplicates, QC failures, secondary, supplementary and unmapped reads, reads with a
mapping quality below 20, reads shorter than 30 bases, reads with a malformed alignment, and
reads where more than 1 aligned base in 10 disagrees with the reference.
[Preparing the inputs](inputs.md#read-filters) gives the exact rules.

On the STR path, a read that spans a repeat is used as a complete observation only if at least
three quarters of its bases there have quality 15 or more.

### Read fixes

pvc adjusts the reads, their qualities and their alignments to avoid some known problems.

- **Indels are left-aligned.** One deletion inside a run of identical bases can be written at
  several positions, and mappers do not always choose the same one. On the generic path, pvc
  moves every indel of every read to its leftmost
  equivalent position, so all the reads that carry one indel show it in one place.
- **Bases read past the end of the fragment are ignored.** When a DNA fragment is shorter than
  the read length, the sequencer reads through the fragment into the adaptor. pvc detects this
  from the pair: the two mates map to opposite strands of one chromosome, and the fragment
  length the mapper reports is shorter than the read. On the generic path, the bases beyond the
  fragment's end are not taken into account.
- **Overlapping mates count once.** Where the two mates of a pair cover the same position, they
  are two readings of one molecule, not two independent observations. If they show the same
  base, one mate keeps it with the sum of the two qualities (at most 200), and the other mate's
  bases there are given quality 0, so they carry no weight. If they disagree, the mate with the
  higher quality keeps its base at 0.8 of its Phred quality, and the other is given quality 0.
  If either shows an indel there, the mate with the lower quality is dropped at that position.
- **A read is as uncertain as its mapping.** On the generic path, pvc gives each read one
  error probability at each locus, written **ε**: the probability that the read shows the
  wrong sequence there. It comes from the read's qualities. A Phred quality *Q* states an error
  probability of 10^(−*Q*/10): Phred 20 is 1 error in 100, Phred 30 is 1 in 1,000. ε is the
  larger of two such probabilities: the one the read's base quality states, and the one its
  mapping quality states. The base quality used is the lowest among the read's bases in the
  locus; for an indel, it is the lowest over the indel and the read bases on each side of it.
- **Base-quality recalibration.** Base qualities are often optimistic. pvc estimates one
  *calibration factor* per read group: how many times more errors its reads make than their
  qualities claim. Every ε of the read group is then multiplied by that factor, except where
  the read differs from an allele by an indel (section 4). With a factor of 2, a base reported
  at Phred 30 (1 error in 1,000) is scored as Phred 27 (1 error in 500). The factor is
  estimated from the cohort's own reads (section 6); with `--defaults` it is 2. For a very poor
  base, the multiplication could give an error above 1 in 2: more likely wrong than right. pvc
  stops it at 1 in 2, where a base counts as no evidence for or against any allele. GATK's
  base-quality recalibration (BQSR; DePristo et al. 2011) solves the same problem with a table
  by reported quality, machine cycle and sequence context, and needs a file of known variants.
  pvc's correction is one number per read group and needs no known variants.
- **At STRs, each read is aligned again.** A mapper places gaps without a model of repeats, so
  it can misjudge where a repeat ends in a read, or whether the read reaches past it at all.
  pvc takes the part of the read that the mapper placed over the repeat, plus 15 bases on each
  side. It aligns that part again against the reference repeat and its flanks, with an aligner
  that makes gaining or losing a whole repeat unit cheap. Then:
  - a read that reaches the flank on both sides is a complete observation, and pvc keeps its
    sequence over the whole repeat;
  - a read that reaches only one side is a partial observation: it shows that the repeat is at
    least as long as the part it covers;
  - a read that reaches neither side is dropped.

Two more corrections act later and are described where they apply. Reads that fit no allele at
an STR are given a small fixed weight (section 4). Variants supported mostly by low-quality
bases are penalised in QUAL (section 9).

### What pvc does not do

- **No BAQ.** BAQ (base alignment quality, from samtools; Li 2011a) lowers the quality of
  bases near possible indels. pvc does not apply it, and does not support reads that already carry it.
  BAQ lowers the bases beside an indel towards quality 0, so pvc would take the reads that show
  an indel in a short run of bases as almost certainly wrong. pvc warns when an alignment's
  header records a `samtools calmd` step that wrote BAQ over the base qualities.
- **No local reassembly.** GATK's HaplotypeCaller (Poplin et al. 2018) rebuilds the reads of
  each active region into haplotypes with a de Bruijn graph and aligns the reads to them again.
  pvc, like
  freebayes, takes the mapper's alignment on the generic path, corrected only by
  left-alignment. This has two consequences. A variant the mapper aligned badly is not rescued.
  And a locus on the generic path can span at most 100 reference bases (`--max-locus-span`), so
  the longest deletion pvc can call there is 99 bases. Inside repeated sequence it is shorter,
  because the deletion's possible positions widen the locus.

## 3. Per-sample pileup files

A psp file holds what one sample's reads show at every position, after the filters and fixes,
but before any statistical model.

- On the generic path there is a record for every position a read covers, including positions
  where every read matches the reference. A record lists each distinct sequence seen. For each
  sequence and each read group, it stores the number of reads, how many are on the forward
  strand, their summed log error probability, and a summary of their mapping qualities. Two
  more counts feed the artifact penalties (section 9): how many of the reads have a base
  quality below 20 at that position, and how many have at least as much of their length to the
  left of the position as to its right.
- On the STR path there is one record per locus, listing each distinct sequence the reads show
  over the repeat and which part of it they covered.
- Reads keep an identifier that the two mates of a pair share. Phasing uses it (section 12).
- At its end, a psp holds a small summary for parameter estimation: the read counts at about
  two million positions and at up to 5,000 repeats of each class, chosen by a fixed
  pseudo-random rule so that they are the same positions in every sample.

A psp holds no genotype and no genotype likelihood. This is the difference from a GATK gVCF. A
gVCF stores, for each sample, likelihoods over the alleles HaplotypeCaller found in that
sample plus a catch-all `<NON_REF>` allele. A psp stores the evidence those likelihoods would
be computed from. The likelihood is computed at calling time. So it uses the candidate alleles
of the whole cohort, and the error rates, contamination and stutter estimated from the whole
cohort.

## 4. Region types and their read models

pvc divides the reference into four types of region. The division depends only on the
reference and on its repeat catalogue ([Preparing the inputs](inputs.md#the-repeat-catalogue)),
so it is the same for every sample.

| Region type | What it is | How it is called |
|---|---|---|
| Generic | everything that is not a repeat, and repeats too short or too impure to qualify | the generic path: SNPs and indels |
| STR segment | a unit of 1 to 6 bases repeated at least 8, 6, 6, 6, 5 or 4 times (by unit length), at least 80% pure, at most 100 bases long | the STR path: one locus, whose alleles are whole repeat sequences |
| STR bundle | two or more repeats within 15 bases of each other, at most 100 bases in all | the STR path: one locus over all of them |
| Satellite | a repeat or bundle longer than 100 bases | not called by default |

The thresholds are options ([Options](options.md#what-counts-as-a-repeat)).

The two paths differ in what they take as an observation and in their read likelihood models.

### The generic path

pvc models a read's errors with the read's error probability ε, defined in section 2: the
error its base quality or its mapping quality states, whichever is larger, multiplied by its
read group's calibration factor. ε is the read's own qualities, turned into a probability and
corrected.

Two alleles appear in the formula. *o* is the observed allele: the sequence the read
shows over the locus. *b* is the allele on the chromosome copy the read came from, which
is unknown. In a sample with genotype *g*, each read comes from one of the genotype's copies,
all equally likely. So the probability of a read that shows *o* is an average over the copies:

```
P(read shows o | g) = (1 − c) · Σ_b (k_b / m) · e(o | b, g)  +  c · q(o)

e(o | o, g) = 1 − ε
e(o | b, g) = ε / 3    when o differs by one substitution from every allele in g
e(o | b, g) = ε        otherwise
```

*e(o | b, g)* is the probability that a read from a copy carrying *b* shows *o*. When the two
are the same allele, the read is right, with probability 1 − ε. When the read differs by one
substitution, the read is wrong in one particular way: an error at that base could have
produced any of the three other bases, so this one gets a third of ε. Any other difference,
such as an indel, is charged the whole of ε. The third is given only when *o* is one
substitution away from every allele the genotype carries: a read one substitution from one
allele of a heterozygote and one indel from the other is charged the whole ε.

*m* is the ploidy, and *k_b* the number of copies of allele *b* in *g*, so *k_b / m* is the
chance that a read came from a copy carrying *b*. The sample's likelihood is the product over
its reads.

The last term is contamination. A share *c* of the read group's reads came from another
individual. That individual shows allele *o* with probability *q(o)*, the allele's frequency in
the rest of the cohort. *q* is re-estimated on every pass of the calling loop (section 8). *c*
and the calibration factor are estimated in section 6; with `--defaults`, *c* is 0 and the
factor is 2.

When *o* differs from the genotype's alleles by one indel, ε is the read's uncalibrated error,
because the calibration is estimated from substitutions.

A read that covers only part of the locus is a partial observation, as on the STR path: a copy
whose allele agrees with the part the read shows gives it 1 − ε, any other copy ε.

Reads of one read group that show the same allele are scored together at the geometric mean of
their error probabilities, which is exact when they share one quality and an approximation
otherwise.

### The STR path and stutter

When a repeat is copied, in the cell or in the PCR of library preparation, the polymerase can
slip and the copy gains or loses whole units. This is stutter. A sample that carries 12
copies of `AC` gives reads with 12, and also some with 11 or 13. The stutter model accounts for
this kind of noise.

An observation at an STR is the read's spelling of the whole repeat. The probability that an
allele *a* gives that spelling is the product of two parts:

- **The length.** If the read is *n* whole units away from the allele:

  ```
  P(n < 0)  = L · u · (1 − ρ) · ρ^(|n| − 1)
  P(n > 0)  = L · (1 − u) · (1 − ρ) · ρ^(n − 1)
  P(n = 0)  = 1 − L − (part-unit changes)
  ```

  *L* is the stutter level, the share of reads that slip by whole units. *u* is the direction
  share: the share of the slipped reads that are shorter than the allele. *ρ* is the fall-off:
  each further unit is *ρ* times as likely as the one before. A length change that is not a
  whole number of units is treated as a sequencing indel. *L*, *u* and *ρ* are estimated per
  class of repeat (unit length × number of copies). *L* is then multiplied by a factor
  estimated for each read group and unit length (section 6).

  Part-unit changes are not estimated. Where the stutter was fitted, they have the same *u* and
  *ρ*, and their total probability is fixed at one twentieth of *L*, on top of *L*: a stated
  value, not a measurement. Where it was not, with `--defaults` and for any class the fit had no
  numbers for, pvc uses HipSTR's built-in values unchanged (Willems et al. 2017): 1 read in 20 a
  unit shorter and 1 in 20 a unit longer (*L* = 0.1), each further unit 20 times rarer, and 1 in
  100 a part of a unit shorter and 1 in 100 longer, so the part-unit changes are one fifth of
  *L*.
- **The bases.** After the allele is brought to the read's length, each base matches with
  probability 1 − ε_rep and differs with ε_rep/3. ε_rep is not the read's base quality. It is a
  substitution rate measured inside repeats, per class of repeat. Sequencers misread the 20
  bases that follow a run of 8 or more identical bases more often than other bases. For a read
  that spans the repeat, pvc adds an extra error rate there, estimated per read group.

The probability of a read is then:

```
P(read o | g) = (1 − λ − c) · Σ_a (k_a / m) · R(o | a)  +  λ · U  +  c · s(o)
```

*R* is the product of the two parts above. *λ* is the weight given to reads that fit no allele
at all: reads from another place in the genome, or chimeric ones. It is 0.07 unless edited in
the parameters file, a value chosen by comparing genotype accuracy at several values. It is
never estimated, and it is not a measurement of how many such reads there are. Such a read is
equally likely to show any length the candidates can reach (*U*). This term limits how much one
odd read can pull a genotype. *c* is contamination, as on the generic path. Here,
though, the contaminant's allele *s(o)* is drawn from the lengths common in that class of
repeat, and it does not change during the calling loop.

### STR bundles

When repeats lie within 15 bases of each other, none has clean sequence on both sides. pvc calls
the bundle as one locus. An allele is the whole sequence from the first repeat to the last.
Within each allele, pvc finds the runs of repeated units it contains. Each run stutters on its
own, with the rates for its unit length and its number of copies in that allele. The
probability of a read is the sum over all the ways the runs could have slipped to give the
read's length.

### STR edges

- A SNP or indel in the sequence beside a repeat is on generic ground and is called on the
  generic path.
- A deletion that runs from a repeat into the sequence beside it is an allele of the STR
  locus. The STR record is widened over the deleted bases, and any SNP or indel record it then
  overlaps is merged into it ([The VCF § STR records](vcf.md#str-records)).
- An interrupted repeat, one that is at least 80% pure but not perfect, is called like any
  other. It is left out when stutter is estimated.

### Satellites

Reads of ordinary length cannot span a repeat longer than about 100 bases with sequence on both
sides, so a satellite's length cannot be read from them. By default pvc calls nothing there.
`--call-satellites` treats satellites as generic ground: SNPs and short indels inside them are
called on the generic path, and changes in the number of units mostly remain uncalled.

## 5. Candidate alleles

pvc is a Bayesian caller, so it needs a way to decide which alleles to consider in calling.

Candidate alleles are chosen over the cohort. An allele is considered if at least one sample
shows it in at least 2 reads and in at least 10% of its reads. The reference allele is always
a candidate.

On the generic path:

- pvc opens a locus only where at least one sample shows something other than the reference:
  at least 2 of its reads, and at least 2% of its reads at that position. Inside a locus, an
  allele still needs the bar above, 2 reads and 10%, to become a candidate.
- Variants that cover some of the same reference bases, such as a SNP inside a deletion, are
  joined into one locus, even when they come from different samples. For example, where the
  reference reads `ACGT`, a deletion of `CGT` and a SNP of `G` to `T` become three alleles of
  one locus: `ACGT`, `A` and `ACTT`. [Overlapping variants](overlapping_variants.md) shows how
  these are written in the VCF.
- A locus longer than 100 reference bases is not called (`--max-locus-span`).
- By default at most 6 alleles are kept, the reference included (`--max-generic-alleles`). They
  are ranked by their best share within one sample, then by the number of samples that passed
  the bar. If a sample passed the bar for an allele that the limit left out, that sample is
  written as not called, because its genotype cannot be expressed.

On the STR path:

- The candidates are repeat sequences seen in the reads, with the same bar of 2 reads and 10%.
  Lengths are chosen first. Each sample nominates the lengths that pass the bar, best-supported
  first, and at most as many as its ploidy. If it nominated fewer, the lengths one unit away
  that its reads show are added. Then every sequence of a nominated length that passes the bar
  is a candidate. If none passes, the best-supported sequence of that length is used.
- A sequence that passes a stricter bar, at least 2 reads and at least 15% of one sample's
  reads, is a candidate even if its length was not nominated. This admits, for example, a
  third length in a diploid sample.
- Some alleles reach into the sequence beside the repeat through a deletion. Such an allele is
  dropped as a misread if it differs by one base from another candidate and makes up less than
  a third of the reads in every sample that shows it.
- At most 32 alleles are kept, the reference included. A record where the limit left alleles
  out is marked `tooManyAlleles`.

## 6. Parameters estimated from the data

Bayesian calling needs some parameters: error rates, contamination, how inbred each sample is,
how variable the population is, how much repeats stutter. `estimate-parameters` measures them on
the cohort and writes them to a parameters file, which the calling step uses as fixed inputs.
[Estimated parameters](parameters.md) describes the file and its built-in values.

### The method: marginal likelihood

To measure an error rate you need to know which reads are errors, which needs the genotypes,
which need the error rate. To avoid this circularity, pvc computes, for any value of the
parameters, the probability of the observed reads. It sums over everything that is unknown at
each position, such as each sample's genotype. This is the *marginal likelihood*. pvc keeps the
parameter values that make the observed reads most probable.

For the positions of the generic path:

```
ℓ(parameters) = Π over positions of
      Σ over kinds of position
      ∫ over the allele frequency f
      Π over samples s  Σ over genotypes g   P(g | f, F_s) · P(reads of s | g, error rates)
```

- The kinds of position are: invariant, fixed for a non-reference base, or segregating. Each
  can be clean or noisy; noisy stands for positions where reads are mismapped.
  Each read group has one error rate at clean positions and another at noisy ones, so the
  noisy positions do not inflate the clean rate. A further kind stands for a duplicated
  sequence: its carriers show the other allele in about half their reads and are never
  homozygous for it.
- The allele frequency *f* is integrated over a Beta distribution whose two shape parameters
  are themselves estimated. It is not fitted position by position.
- *P(g | f, F_s)* is the genotype prior with the sample's inbreeding coefficient (section 7).

The input data are the summaries written at the end of the psp files: about two million
positions, the same in every sample. The search is an EM algorithm (Dempster et al. 1977), sped
up with Newton steps, run from three starting points that differ in the error rates they assume.
It stops when the remaining step could raise the log-likelihood by less than 0.01.
`estimate-parameters` fits diploid samples only; at another ploidy, calling uses the built-in
values.

### What is estimated, and what is fixed

The first fit, the marginal likelihood above, has eight numbers for the whole cohort, two for
each read group and one for each sample. The other numbers come from separate fits or counts.
The last column says whether calling uses the number.

| Quantity | For each | How it is obtained | Used in calling |
|---|---|---|---|
| Error rate at clean positions | read group | the first fit | yes, as the base-quality calibration factor |
| Error rate at noisy positions | read group | the first fit | no |
| Share of noisy positions | cohort | the first fit | no; the contamination fit leaves out the positions this fit judges more likely mismapped than not |
| Homozygote excess, used as the inbreeding coefficient F | sample | the first fit | yes (section 7) |
| Share of invariant positions; share of positions fixed for a non-reference base; the two shapes of the Beta distribution of the segregating allele frequencies | cohort | the first fit | yes, through the prior's seed (section 8) |
| Share of duplicated positions; the two shapes of the Beta distribution of the share of samples that carry a duplication | cohort | the first fit | no |
| Contamination | read group | a second marginal likelihood, after the first | yes |
| Stutter level, direction share and fall-off | class of repeat, with a multiplier of the level per read group and unit length | a marginal likelihood over the repeats | yes |
| Typical allele lengths and their spread between repeats | class of repeat | the same STR fit | yes, as the STR prior's seed (section 8) |
| Substitution rate inside repeats | class of repeat | counted | yes |
| Extra errors after a run of one base | read group | counted | yes |

Some numbers that change genotypes are never estimated:

| Fixed value | Where it acts |
|---|---|
| Weight of an STR read that fits no allele, λ = 0.07 | STR read likelihood (section 4) |
| Part-unit stutter at one twentieth of the fitted level | STR read likelihood (section 4) |
| The built-in values of `--defaults`: calibration factor 2, no contamination, F = 0, HipSTR's stutter, a seed of 1 and 1/1,000, a between-repeat share of 1/11, a substitution rate inside repeats of 1/1,000 ([Estimated parameters](parameters.md#the-built-in-values---defaults)) | wherever the fit gave no number |
| Error rate of a read group without data, 1/1,000 | base-quality calibration |
| The constants of the three QUAL penalties | section 9 |
| The refinement's copying mismatch rate and the weight of its "none of them" template | section 11 |

Every fitted number in the parameters file carries a `standard_error` where one could be
computed: from how sharply the likelihood falls away from the number, or, for a counted rate,
from the spread of the count. A paper can report it beside the value. Two pairs of numbers trade
against each other along what the reads cannot see: the share of invariant positions against
the first Beta shape, and the share of duplicated positions against the first carrier shape.
Where the data do not separate them, the file marks both `separately_determined = false`, gives
neither a standard error, and also gives the share of positions segregating (or duplicated)
above a stated frequency, which is the number to quote.

**Base-quality calibration.** The parameter fit calculates each read group's true error rate
at clean positions. The calibration factor is that rate divided by the rate the reads claim.
The claimed rate is a geometric mean over every read at every position it covers, and what a
read claims is the worse of its base and mapping qualities (section 2). A read group without
data gets a default error rate of 1 in 1,000.

**Contamination.** Each read comes from the sample with probability 1 − *c* and from another
individual with probability *c*. The fit sums over the genotypes of both, and draws them at
different allele frequencies, following verifyBamID2 (Zhang et al. 2020):

- *The sample's own genotype* is drawn, with its F, at allele frequencies of its own. pvc places
  every sample on the cohort's main axes of variation, up to four (a principal-component analysis of
  the samples' genotypes), and at each position fits the allele frequency as a straight line in
  those coordinates, using every sample. A sample from a diverged group thus expects that
  group's alleles. Scored against the pooled frequencies instead, the structure hides
  contamination: on a simulated panel of four subpopulations at an F_ST of 0.2, a
  sample contaminated
  at 3% was then estimated at 0.5%.
- *The contaminant's genotype* is drawn at the cohort's pooled frequencies, without inbreeding:
  a contaminating library is whichever other sample was handled beside it, not one of similar
  ancestry.

The fit keeps the *c*, between 0 and 0.5, that makes the reads most probable, and reports 0
unless a positive value explains the reads better. It uses only positions with reads in at least
8 samples, a pooled allele frequency between 2% and 98%, and no sign of mismapping, and needs
at least 100 of them; with fewer than 8 samples there is no estimate. A sample that alone
supplies more than half of its own fitted frequency, such as one at the far end of an axis, is
reported as *not identified*, and calling treats its read groups as uncontaminated.

**The estimate is biased low.** On a simulated panel of 40 samples in four subpopulations, one
sample contaminated at 3% was estimated at 2.6% at 30 reads a position and at 1.0% at 3 reads,
while the 39 clean samples were at or near 0. It separates a contaminated sample from clean
ones; it understates the share, most at low depth.

**Stutter.** At each sampled repeat the parameter fit sums over the sample's possible pairs of
allele lengths, with the sample's F, and asks which stutter numbers make the observed spread of
read lengths most probable. The numbers are not free for every class: they are smooth curves
along the number of copies, one set per unit length. Three kinds of repeat are left out of this
fit:

- interrupted repeats;
- repeats where more than 1 in 10 of the reads that differ from an allele differ by part of a
  unit;
- repeats whose depth, or excess of reads of other lengths, suggests a duplicated sequence.

The multiplier of each read group is pulled towards 1 when the read group has few reads: one
fitted from fewer than about 20 slipped reads leans on that pull as much as on its reads. It is 1
when a sample has a single read group.

The substitution rate inside repeats is the number of mismatches divided by the number of
bases compared, per class, pooled over all read groups, in reads that show the reference
length. The extra error after a run of 8 or more identical bases is the difference between
reads sequenced through the run and reads that reached the same base from the other side.

## 7. Inbreeding

The inbreeding coefficient F of a sample is the probability that its two copies at a locus
are identical by descent. A selfed line has F near 1, an outbred individual near 0.

For a diploid sample at a site with two alleles of frequencies *p* and *q*:

```
P(AA) = p² + F·p·q        P(AB) = 2·p·q·(1 − F)        P(BB) = q² + F·p·q
```

In general, for any ploidy and any number of alleles, the prior is a mixture: with probability
1 − F the genotype is drawn under random mating, and with probability F the sample is
homozygous for one allele, drawn at its frequency. Above diploidy this is a simplification:
the inbreeding term goes only to genotypes in which every copy is the same allele, so a
tetraploid AABC gets none of it.

With F = 0, the Hardy–Weinberg prior expects 2*pq* heterozygotes. In a selfing species there
are almost no heterozygotes.

F is estimated as one of the parameters of the fit in section 6: the F of each sample is the
value that best explains how few, or how many, heterozygous positions its reads show, given the
allele frequencies of the cohort. This is Wright's inbreeding coefficient (Wright 1922)
measured as a homozygote excess, and three things follow from how it is measured:

- **It is measured against the pooled frequencies of the whole cohort.** In a cohort with
  population structure it therefore also takes up the deficit of heterozygotes that the
  structure causes (the Wahlund effect). It is closer to Wright's F_IT (Wright 1951) than to an
  individual's autozygosity, and can differ from a pedigree F or one from runs of homozygosity.
  It is the quantity a prior built on the pooled frequencies needs.
- **It is at least 0.** An excess of heterozygotes, as in an F1 hybrid or a clonal
  heterozygous crop, is reported as 0.
- **It is below 1.** At F = 1 the prior would give every heterozygote a probability of zero,
  and no read evidence could produce a heterozygous call. For a fully inbred line, 0.99 puts
  20 Phred against a heterozygote, which clean reads can overcome.

With one sample, F cannot be estimated. One genome cannot show whether it has too few
heterozygotes: there are no allele frequencies to compare with. F is set to 0 and marked
`defaulted` in the parameters file. For an inbred line called alone, the user can supply the
value with `estimate-parameters --inbreeding`. The value must be below 1: pvc refuses
`--inbreeding 1`.

The F of the parameters file is used in three places: the calling prior (section 8), the
hidden-duplication filter (section 10), and the refinement (section 11). The stutter and
contamination fits use the F estimated within the same run, not a supplied one.

## 8. Calling

At each locus, two things are unknown and each depends on the other. A sample's genotype prior
depends on the allele frequencies, and the allele frequencies depend on the samples'
genotypes. pvc solves this with an EM loop that alternates between the two.

For each sample *s* and each candidate genotype *g*:

```
prior:       α_s(a)  =  seed(a)  +  Σ over the other samples t of  E[copies of a in t]

             P(g | α_s, F_s)  =  (1 − F_s) · DM(g | α_s)  +  F_s · [g is homozygous for a] · α_s(a) / Σ α_s

posterior:   P(g | reads_s)   ∝  P(g | α_s, F_s) · P(reads_s | g)

update:      E[copies of a in s]  =  Σ over g of  P(g | reads_s) · (copies of a in g)
```

In words:

- The prior of a sample is built from the other samples. α_s(a) is the number of copies
  of allele *a* the rest of the cohort is expected to carry, plus a seed. The sample's own
  copies are left out, so its reads are not counted twice, once as evidence and once as prior.
- DM is the Dirichlet-multinomial distribution: the probability of drawing the genotype's
  copies from a population whose allele counts are α_s. With large counts it is the multinomial
  of Hardy–Weinberg; with small ones it allows for the frequencies being uncertain. This is the
  genotype-sampling distribution of Balding and Nichols (1995).
- The seed is what is known before any sample is seen: a pseudo-count for the reference
  allele and a much smaller one shared among the alternatives. It comes from the distribution
  of allele frequencies fitted in section 6. The two pseudo-counts are chosen so that the
  Dirichlet they define has the same mean alternative-allele frequency and the same expected
  heterozygosity, E[2*f*(1 − *f*)], as the fitted distribution. On a population at neutral
  equilibrium this gives 1 for the reference and θ, the heterozygosity, for the alternatives.
  The alternatives' pseudo-count is shared out evenly among a site's alternative alleles, so a
  site with three alleles is not taken as more variable than one with two. With built-in
  parameters the seed is 1 for the reference and 0.001 shared among the alternatives.
- At an STR, the seed of each length follows how common that length, counted from the reference
  length, is in repeats of that class. Its total is fitted per class too, as how much the allele
  lengths vary from one repeat to the next. The parameters file writes this as the
  *between-repeat share*, 1/(1 + total): a total of 10, the built-in value, is a share of 1/11.
- F_s enters as in section 7.
- The expected copies are recomputed from the posteriors, and the loop repeats.

The loop starts from one pass without a prior. It stops when no allele's expected count moves
by more than one thousandth of the total allele copies in the samples with reads, or after 50
passes. A locus that reaches 50 is written with the filter `EMNoConv`. A sample without reads
at the locus adds nothing to the counts. When contamination was estimated, its allele
frequencies are also updated on each pass.

At the end, each sample's genotype is the one with the highest posterior. A sample whose reads
give every genotype the same likelihood, which includes a sample with no reads, is written
`./.`. It is never called homozygous for the reference for lack of evidence.

The allele frequency written as `AF` is the loop's estimate. Alternative alleles that no
called genotype carries are first removed from the record, and `AF` is rescaled over the
alleles that remain.

### Compared with GATK and freebayes

GATK's HaplotypeCaller in gVCF mode (Poplin et al. 2018) first calls one sample at a time and
writes its genotype likelihoods; no cohort exists yet. Joint genotyping then estimates the cohort's allele
frequencies with a loop like pvc's. It alternates between genotype posteriors given the
frequencies, and frequencies given the expected allele counts. Its prior is a set of fixed
pseudo-counts, and it stops when no expected count moves by more than 0.1. Four things differ.

- In GATK the loop's result is used only for the site: its QUAL, and which alleles are kept. Each
  sample's written genotype and GQ come from its own likelihoods alone. In pvc the written
  genotype is the posterior, so the cohort's allele frequencies change individual genotypes.
  This matters most where a sample has few reads, and least where its reads alone settle the
  genotype.
- GATK's prior pseudo-counts are, by default, constants derived from human diversity (a
  heterozygosity of 1 in 1,000 for SNPs and 1 in 8,000 for indels). pvc's seed is estimated
  from the cohort.
- GATK's prior is the same for every sample. pvc leaves each sample's own expected copies out
  of its prior.
- GATK's prior assumes Hardy–Weinberg. pvc's has F per sample.

freebayes (Garrison and Marth 2012) does not iterate allele frequencies. Instead, it searches
over combinations of genotypes, one per sample. It starts from each sample's most likely genotype and tries
neighbouring combinations, changing one sample at a time. Each combination is scored by the
samples' likelihoods times a prior on the combination. That prior includes the probability of
its allele counts under a neutral population model (Ewens' sampling formula; Ewens 1972) and a
Hardy–Weinberg term. The number of combinations grows exponentially with the number of
samples, so they are searched, not enumerated. The best combination found decides which
alternative alleles the site carries, and QUAL comes from the posteriors of the combinations
explored. By default, each sample's written genotype is then its own most likely genotype over
those alleles. So, as in GATK, the cohort chooses the alleles but not the individual genotypes.
pvc's loop, in contrast, does one pass over the samples per iteration, and its cost grows in
proportion to their number.

## 9. QUAL and GQ

GQ is the Phred-scaled probability that the sample's written genotype is wrong: one minus
its posterior, at most 99.

QUAL is the Phred-scaled probability that no sample in the cohort carries a non-reference
allele, given all the reads. That probability is computed directly. It is not the product of
each sample's probability of being homozygous for the reference. Each sample's likelihoods are
reduced to the likelihood of carrying 0, 1 or 2 non-reference copies (up to the ploidy). These
are combined over the samples into the likelihood of each total number of non-reference copies
in the cohort, the dynamic programme of Li (2011b). A prior is applied to that total, built
from the same seed as in section 8. This value is the posterior probability of a total of
zero, what GATK computes as the probability that the cohort's allele count is zero.

A sample whose reads give every genotype the same likelihood, such as one with no reads, drops
out of the calculation exactly. The prior sets how much read evidence a site needs before its
quality climbs above 0: the prior odds against the site being variable. At θ = 1/1,000 they are
28 Phred at one diploid sample, 23 at 63 and 20 at 3,000; at θ = 1/100, 18, 13 and 10.

Three penalties are then subtracted, each a test for a known kind of artifact. Each looks at
the reads of the main alternative allele, pooled over the called samples, and charges the
Phred-scaled probability of a result at least as far from what is expected.

| Penalty | The artifact | The test |
|---|---|---|
| `ABPEN` | too few alternative reads for the genotypes called | the alternative reads against the share the called genotypes predict, a two-sided binomial tail; only a shortfall is charged, and the test is skipped where the genotypes predict 90% alternative reads or more |
| `SPPEN` | alternative reads on one strand, or all at one end of the reads | the share of alternative reads on the forward strand, and the share placed mostly left of the site (an idea from freebayes), each against the reference reads' share as a two-sided beta-binomial tail; the larger of the two, times 4 |
| `BQPEN` | the same base misread in many reads, typically late in the read beside a run of G or C | the share of alternative reads with a base quality below 20 against the share among the reference reads, a one-sided beta-binomial tail, times `--low-base-quality-weight` (1 by default); only an excess is charged |

The beta-binomial tests allow for the reference reads' share being itself an estimate: it is
the reference reads' counts plus a prior worth 10 reads (an even strand split, or the site's
pooled share for the position) or 30 reads (at 2 low-quality reads in 100). The strand test is
weighted 4 times because QUAL grows by about 30 Phred with each alternative read and the tail of
a one-sided pile-up by a few. `SPPEN` and `BQPEN` are not charged at 3 alternative reads or
fewer, where the tests have no power, and grow linearly to full strength at 7. STR records carry
no penalties.

The written QUAL is the value before penalties minus the three, at least 0 and at most 9999.
It is no longer a probability: it is a score that ranks sites, which equals the probability
above only where every penalty is 0. [The VCF](vcf.md#qual-and-gq) lists the fields.

A site is written when its QUAL after penalties is at least `--min-site-quality` (1 by
default) and at least one sample is called with an alternative allele.

## 10. Site filters

**Hidden duplications.** Where a genome carries two copies of a sequence and the reference
only one, the reads of both copies map to one place. Every difference between the copies then
looks like a heterozygous variant, in every sample that carries both. No read-level model can
see this, because each read is correct.

For every called site, pvc compares two explanations of the reads:

- *A real variant in single-copy sequence.* Every sample has the depth one copy gives, and its
  reads follow its genotype under the prior with its F.
- *A hidden duplication.* The samples that carry it have more depth, in proportion to the
  number of copies, and show the alternative allele in the share of reads the copies predict.
  They are never homozygous for the alternative allele.

Each explanation is a likelihood summed over its unknowns: the allele frequency for the first,
the number of copies and the share of carriers for the second. The depth one copy gives is
measured per sample, from its own reads in 500-base windows, corrected for GC content. The log
ratio of the two likelihoods is written as `PARALOG_LR`.

A ratio says which explanation fits better, not how probable a duplication is. For that, pvc
also estimates how common hidden duplications are among the called sites: over the run, and
locally in windows of 50 kb by default, because duplications tend to cluster. Ratio and rate
together give each site's probability of being a hidden duplication (`PARALOG_POST`). pvc ranks
the sites by that probability and removes them from the top down. It stops before the expected
share of true variants among the removed sites would exceed 1 in 100 (`--hidden-dup-fdr`), a
false discovery rate computed from the sites' probabilities (Newton et al. 2004). With
`--hidden-dup-tag` the removed sites are kept and marked `hiddenDup` instead.

This filter is the reason the calls of a run wait on disk until the whole genome is called
([Calling § scratch space](calling.md#output-logs-and-scratch-space)).

Other filters:

- A locus where any read group of any sample is deeper than `--max-depth` (1,000) is not
  called. Such depth usually means collapsed copies.
- `EMNoConv` marks a site where the loop of section 8 did not settle.
- `tooManyAlleles` marks an STR where the limit of 32 candidates left alleles out.

## 11. Refinement

Calling treats each locus alone. Yet the samples of a cohort share stretches of chromosome,
and that fact can be used to improve the genotypes. For instance, a sample's genotype at a
site with one read could be settled by taking into account the genotypes of the other samples
that share the same haplotype. The optional **refinement**, what is often called imputation,
uses that kind of information. The haplotypes of the other samples act as a prior, and the
sample's own genotype likelihoods as the evidence. Where a sample has no reads at a site, the
posterior is the prior alone: what the haplotypes it copies predict.

Along a chromosome, each of a sample's two copies is modelled as a mosaic of *templates*:
stretches copied from haplotypes that other samples carry. This is the copying model of Li and
Stephens (2003) that most imputation tools use. Three things are particular to pvc:

- The templates come from the cohort itself. For each sample and each window of the
  genome, up to 16 other samples are chosen, those that agree best with it over stretches of
  the window. There is no reference panel. One more template stands for "none of them": it
  carries the cohort's allele frequencies.
- Inbreeding is a state of the model. A sample is, stretch by stretch, either copying two
  templates, one per chromosome copy, or *autozygous*: both copies are the same template. The
  long-run share of the autozygous state is the sample's F from the parameters file. In an
  inbred line almost the whole genome is in that state, and the model then has one template to
  find, not two.
- Distances are physical. No genetic map is used. The typical length of a copied stretch is
  either 30 kb or 300 kb. It is chosen per sample, together with the typical length of its
  autozygous stretches: pvc keeps the values that better fit the sample's first 5,000
  informative sites.

The model's numbers, between two sites *d* bases apart:

- A copy keeps its template with probability e^(−*d*/*S*), where *S* is the typical length of a
  copied stretch (30 or 300 kb, above). A copy that changes template lands on the "none of them"
  template with probability 1 in 100,000, and on the others in proportion to the stretches of
  the window over which each agreed best with the sample.
- The autozygous state forms runs of typical length *R*, chosen per sample from 10 kb, 30 kb,
  100 kb, 300 kb, 1 Mb and 3 Mb, and the gaps between runs have typical length
  *R*·(1 − F)/F, so that over a long distance the share of the autozygous state is F. F is the
  parameters file's, held between 0.01 and 0.99 and not re-estimated.
- A template's allele differs from the copy's with probability 0.01.

The 1 in 100,000 and the 0.01 are fixed values, not estimates.

The refinement works on the calls, not on the alignments.

| Input | Used for |
|---|---|
| the sample's own genotype likelihoods (`PL`) | the evidence at each site |
| the other samples' genotype posteriors (`GP`), then their estimated haplotypes | the templates |
| the allele frequencies (`AF`) | the "none of them" template |
| the sample's F | the share of the autozygous state |
| read pairs that cover two heterozygous sites | the orientation of the sample's drawn haplotypes, and the phase (section 12) |

Reads do not enter a genotype's evidence a second time: they already made the likelihoods.

The genome is taken in windows of 2 Mb, with 1 Mb of margin on each side. On a large cohort the
windows are shortened, down to 500 kb, to keep the refinement within its memory budget
(`--refine-window-bytes`). In each window the refinement makes five rounds. In the first, the
templates are the other samples' genotype posteriors. In each round, every sample is solved
against the others' current haplotypes. A forward–backward pass over its sites gives the
probability of each state. One path of states is drawn at random, and the sample's two
haplotypes are estimated along that path. From the second round on, the templates are those
drawn haplotypes. Where a template copied this same sample in the round before, and the
sample's own reads leave its genotype unsure, that template is read with the cohort's
frequencies in place of the stretch it copied, so a sample does not count its own reads twice.
Drawing paths makes the method a sampler, like Beagle 4.1's, not an exact calculation: the
genotype written is the average of the posteriors of rounds two to five, not a converged
posterior. The random draws are seeded, so a run gives the same output at any number of
threads.

It applies to SNPs, indels and STRs alike, at sites with any number of alleles, in diploid
samples only.

### Compared with Beagle and STITCH

| | Evidence per sample | Where the haplotypes come from | Method |
|---|---|---|---|
| Beagle 4.1 | genotype likelihoods from a VCF | the cohort, and a reference panel of phased haplotypes if given | hidden Markov model of haplotypes; samples haplotypes for each sample in turn |
| STITCH | the reads, at a given list of biallelic SNPs | *K* ancestral haplotypes, estimated; no panel | EM over the ancestral haplotypes; built for low coverage |
| pvc | genotype likelihoods; read pairs for the phase | up to 16 other samples of the cohort per window; no panel | hidden Markov model with an autozygous state; five rounds of sampled haplotypes |

Beagle 4.1 (Browning and Browning 2016) is compared rather than Beagle 5 because it takes
genotype likelihoods. STITCH is Davies et al. (2016). Like Beagle, pvc starts from genotype
likelihoods and needs a caller before it; in pvc, that caller is pvc itself. Like STITCH, it
needs no panel and is aimed at many samples at low depth. Unlike both, it carries each
sample's inbreeding in the model. Unlike STITCH, it is not limited to biallelic SNPs.

The refinement is off by default and runs only when asked for with `--refine-genotypes`.

[Accuracy](accuracy.md#imputation-genotype-refinement-and-filling) compares pvc with GATK
followed by Beagle, and with STITCH, on simulated inbred and outbred cohorts.

**In outbred cohorts the results are mixed, and measured only in simulations.** In an outbred
population of 400 samples simulated with msprime (100 samples each at 30×, 10×, 5× and 2×), the
refinement cut the variant genotypes missed at 2× from 12,220 to 1,546, and 7 of 94,860 phase
links were wrong ([Accuracy](accuracy.md#imputation-genotype-refinement-and-filling)). In a
cohort of 40 samples from the same kind of simulation, all at 2×, it got more genotypes right
(29,381 to 31,378) but fewer heterozygotes (7,521 to 7,042). Earlier runs, of an older version,
on cohorts built from 64 founder haplotypes at 3×, found it removed 13% of the wrong genotypes
at 50 samples and gave 74% more at 10 samples, and 4.7% of its phase links wrong at 50 samples.
[Genotypes from neighbouring sites](neighbour_evidence.md#when-it-runs) gives the details. This
is why `--refine-genotypes auto` runs it only on inbred cohorts.

## 12. Phase

Two heterozygous genotypes of a sample are *phased* when it is known which alleles lie on the
same chromosome copy. pvc writes phase for diploid samples from two sources.

**From reads.** A read pair comes from one molecule. If it covers two heterozygous sites, the
alleles it shows there lie on one copy. For each pair of nearby heterozygous genotypes, pvc
adds up the read pairs that cover both. Each read pair votes for one of the two arrangements.
Its vote weighs less when its bases are more likely to be wrong, and no base counts as more
reliable than 1 error in 100. The two genotypes are linked when the odds reach 199 to 1, which
takes at least two read pairs. Only PASS genotypes with GQ of 30 or more are linked, and STR
records never are. Read pairs reach about a kilobase, so the blocks are short. This is all the
phase there is when the refinement does not run, including with one sample.
[The VCF](vcf.md#phase-from-reads) gives measurements.

**From the cohort.** The refinement estimates each sample's two haplotypes (section 11). But
one estimate is one random draw, and a link between two sites is worth writing only if it is
stable. So the refinement draws 36 paths per sample over its rounds. It links two neighbouring
heterozygous genotypes when at least 80% of the paths agree on the arrangement
(`--phase-min-path-agreement`).

**How the two combine.** The cohort reaches far, over the tens of kilobases a copied stretch
spans, but it knows nothing where a sample's haplotype is rare. Reads reach a kilobase, and
are direct evidence. When the refinement runs:

- Read pairs take part in each haplotype draw. Where a sample's read pairs connect
  heterozygous sites, the draw orients those sites together, weighted by the reads' odds. This
  also improves the haplotypes that later serve as templates for other samples.
- A link the paths agree on is kept unless the sample's reads contradict it at odds of 9 to 1.
  There the block ends.
- Where the paths do not agree, reads at odds of 199 to 1 make the link.

So the refinement phases further than reads alone. In samples whose haplotypes others share,
blocks extend beyond a kilobase, across stretches with no read pair. And reads improve the
cohort's phase where it is weakest. Some read pairs are held out of the draws. The log reports
how often they contradict the links made from the cohort, which checks the cohort's phase
independently of its own model.
[Genotypes from neighbouring sites](neighbour_evidence.md#phase--and-ps) gives the measured
error rates.

## 13. Across the range of inputs

pvc is meant for one sample to several thousand, and for about 3 reads a position to several
hundred. Not every method works equally over that range. This table says what each one does
at the ends.

| Method | One sample | Thousands of samples | About 3× | Hundreds × |
|---|---|---|---|---|
| Read fixes (section 2) | unchanged: all are per read | unchanged | unchanged | at most `--max-depth` reads per read group are kept at a position; a locus deeper than that is not called |
| Candidate alleles (5) | the bar applies to that sample alone | one clear carrier is enough for the whole cohort; at most 6 alleles, so rare ones can be left out at very variable sites | the 2-read minimum decides: a heterozygote needs 2 reads of the allele in some sample | the 10% share decides: alleles under 10% of every sample's reads, such as low-level mosaics, are not considered |
| Error rate and calibration (6) | estimated | estimated; the fit reads every sample at the sampled positions | estimated, with wider uncertainty per read group | estimated |
| Contamination (6) | not estimated, nor with fewer than 8 samples: there is no cohort to give each sample its expected allele frequencies or to draw the contaminant from | estimated | needs positions covered in at least 8 samples | estimated |
| Inbreeding F (7) | not estimable: set to 0 and marked `defaulted`; supply a value below 1 with `estimate-parameters --inbreeding` | estimated per sample | estimated: many positions make up for few reads at each | estimated |
| Population diversity (6) | estimated | estimated | estimated | estimated |
| Stutter (6) | estimated; the spread of allele lengths between repeats is not, and takes a built-in value | estimated | estimated, since reads are pooled over the repeats of a class | estimated |
| Calling loop (8) | the prior is the seed alone, so the loop settles at once | cost in proportion to the samples; rare alleles get a prior that reflects their rarity | the prior carries much of the weight, so F and the allele frequencies matter most here | the reads decide; the prior changes few genotypes |
| QUAL (9) | the probability that this sample is not homozygous for the reference | samples without reads do not change it; the prior asks for less evidence as the cohort grows, 20 Phred at 3,000 samples against 28 at one (at θ = 1/1,000) | `SPPEN` and `BQPEN` are off at 3 alternative reads or fewer | all penalties at full strength |
| Hidden-duplication filter (10) | runs | runs | runs with less power: at low depth a sample's depth at one site is noisy | runs |
| Refinement (11) | cannot run | 16 templates per sample whatever the cohort size; the window shortens to bound memory | where it helps most, in inbred cohorts; in outbred cohorts results are mixed (section 11) | changes little: the reads already settle the genotypes |
| Phase (12) | from reads only | from the cohort and from reads | few read pairs cover two sites: from reads alone, 20 of 289 heterozygous SNPs (7%) were phased in 100 regions of the human sample HG002 at 3× | from reads alone, 197 of 289 (68%) in the same regions at 30×, and little more beyond, because 32 reads per site are kept for phasing |

A single sample at low depth is the hardest case, and the weakest: there is no cohort to give
allele frequencies, F or contamination, and few reads to outweigh the prior. In that case pvc
prefers to leave a genotype uncalled, or a parameter marked `defaulted`, rather than write a
confident value it has no evidence for. [Accuracy](accuracy.md) reports the measured accuracy
in each of these cases.

## References

- Balding DJ, Nichols RA (1995). A method for quantifying differentiation between populations at
  multi-allelic loci and its implications for investigating identity and paternity. *Genetica*
  96:3–12.
- Browning BL, Browning SR (2016). Genotype imputation with millions of reference samples.
  *American Journal of Human Genetics* 98:116–126.
- Davies RW, Flint J, Myers S, Mott R (2016). Rapid genotype imputation from sequence without
  reference panels. *Nature Genetics* 48:965–969.
- Dempster AP, Laird NM, Rubin DB (1977). Maximum likelihood from incomplete data via the EM
  algorithm. *Journal of the Royal Statistical Society, Series B* 39:1–38.
- DePristo MA, Banks E, Poplin R, et al. (2011). A framework for variation discovery and
  genotyping using next-generation DNA sequencing data. *Nature Genetics* 43:491–498.
- Ewens WJ (1972). The sampling theory of selectively neutral alleles. *Theoretical Population
  Biology* 3:87–112.
- Garrison E, Marth G (2012). Haplotype-based variant detection from short-read sequencing.
  arXiv:1207.3907.
- Li H (2011a). Improving SNP discovery by base alignment quality. *Bioinformatics*
  27:1157–1158.
- Li H (2011b). A statistical framework for SNP calling, mutation discovery, association mapping
  and population genetical parameter estimation from sequencing data. *Bioinformatics*
  27:2987–2993.
- Li N, Stephens M (2003). Modeling linkage disequilibrium and identifying recombination hotspots
  using single-nucleotide polymorphism data. *Genetics* 165:2213–2233.
- Newton MA, Noueiry A, Sarkar D, Ahlquist P (2004). Detecting differential gene expression with
  a semiparametric hierarchical mixture method. *Biostatistics* 5:155–176.
- Poplin R, Ruano-Rubio V, DePristo MA, et al. (2018). Scaling accurate genetic variant discovery
  to tens of thousands of samples. bioRxiv 201178.
- Rubinacci S, Ribeiro DM, Hofmeister RJ, Delaneau O (2021). Efficient phasing and imputation of
  low-coverage sequencing data using large reference panels. *Nature Genetics* 53:120–126.
- Willems T, Zielinski D, Yuan J, Gordon A, Gymrek M, Erlich Y (2017). Genome-wide profiling of
  heritable and de novo STR variations. *Nature Methods* 14:590–592.
- Wright S (1922). Coefficients of inbreeding and relationship. *American Naturalist* 56:330–338.
- Wright S (1951). The genetical structure of populations. *Annals of Eugenics* 15:323–354.
- Zhang F, Flickinger M, Taliun SAG, et al. (2020). Ancestry-agnostic estimation of DNA sample
  contamination from sequence reads. *Genome Research* 30:185–194.
