# Read-backed phasing, and how population phasing consumes it

*Research note, 2026-09-12. Part of the phasing/imputation research set; the shared brief is
`tmp/phasing_research/brief.md`. Code citations are `path:line` into the shallow clones under
`tmp/phasing_research/repos/` (paths below are relative to the repository root) and into the
vendored `freebayes/` and `gatk/` trees. The arithmetic in §2 is this note's own, not a
measurement; everything quoted from a paper says so.*

---

## 0. Vocabulary

Genetics terms are assumed. The statistics and software terms this note leans on:

- **Read-backed phasing** (also *haplotype assembly*, *physical phasing*). Deciding which alleles at
  two or more heterozygous sites sit on the same chromosome copy, using only the reads of that one
  sample. A read (or a read pair) that covers two heterozygous sites was sequenced from one
  molecule, so the two alleles it shows were on the same haplotype.
- **Phase-informative read (PIR).** SHAPEIT2's name for such a read: one that "span[s] two or more
  heterozygote genotypes" of the same individual (Delaneau et al. 2013, cited in §3.1).
- **Fragment.** HapCUT2's word for the unit of read evidence after mates and split alignments are
  joined: the list of (site, allele, quality) a molecule reports. WhatsHap calls the same thing a
  read "with a hole in the middle".
- **Phase block / phase set.** A run of heterozygous sites that reads connect into one chain, so
  their alleles are phased relative to each other. Between blocks nothing is known. The VCF field
  `PS` carries the block identity (§5).
- **Switch error rate (SER).** Walking along the true haplotypes, the fraction of consecutive
  heterozygous-site pairs at which the inferred phase flips from one haplotype to the other.
- **MEC, minimum error correction.** The classical objective for read-backed phasing: find two
  haplotypes such that the smallest number of read bases must be called sequencing errors for every
  read to agree with one of them. *Weighted* MEC (wMEC) charges each corrected base its base
  quality instead of one.
- **Fixed-parameter tractable.** A problem that is hard in general but whose running time is
  exponential only in one small quantity (here, the read depth at a site) and linear in everything
  else. Fixing that quantity makes it fast.
- **Max-cut.** Given a graph whose edges carry weights, split the vertices into two sides so that
  the total weight of edges crossing the split is as large as possible. NP-hard in general;
  HapCUT2 uses a greedy approximation.
- **HMM, forward-backward, Gibbs sampling.** Defined in the companion notes on statistical phasing;
  here only the point where read evidence enters those models is described.
- **Genotype likelihood.** For one sample at one site, the probability of the observed reads given
  each possible genotype — the per-sample evidence an imputer takes when hard calls are unreliable.

---

## 1. The signal, the objective, and the two working solutions

### 1.1 What a read pair says

A diploid sample has two haplotypes. At a heterozygous site the reads split between them, and a
single site says nothing about which allele goes with which allele at the next site. A read that
covers **two** heterozygous sites does: it came from one molecule, so the two alleles it shows are
on the same haplotype, up to a sequencing error at either base and up to the read being mapped to
the wrong place. Chaining such reads links sites into blocks; where no read links two neighbouring
heterozygous sites, the chain breaks and a new block starts.

Homozygous sites carry no phase information and every tool drops them before building the chain:
WhatsHap keeps only reads covering at least two variants after merging mates
(`repos/whatshap/whatshap/cli/phase.py:518-520`), and HapCUT2 skips homozygous sites when it adds
graph edges (`repos/HapCUT2/hapcut2-src/variantgraph.c:210`).

### 1.2 Minimum error correction, and why it is hard

Write the reads as a matrix: one row per fragment, one column per heterozygous site, an entry of
0/1 for the allele the fragment shows and a blank where it does not cover the site. If reads were
error-free, the rows could be split into two groups, each group agreeing column by column, and the
two column-wise consensus strings would be the haplotypes. With errors they cannot, and MEC asks for
the smallest number of entries to flip so that they can. Weighted MEC is the same with a cost per
flipped entry: Garg et al. 2016 state it as "flip entries in F to obtain a feasible matrix, while
minimizing the sum of incurred costs, where flipping entry F(j,k) incurs a cost of W(j,k)", with the
weight the phred-scaled base quality (Garg, Martin, Marschall, *Bioinformatics* 32(12):i234–i242,
2016, doi:10.1093/bioinformatics/btw276, Problem 2).

MEC is NP-hard even when every fragment is gapless — one contiguous run of covered sites — and
APX-hard (no polynomial scheme approximates it arbitrarily well) once a fragment may have one gap
(Cilibrasi, van Iersel, Kelk, Tromp, "The complexity of the single individual SNP haplotyping
problem", *Algorithmica* 49(1):13–36, 2007; the arXiv version q-bio/0508012 states: "We prove that
MEC is APX-hard in the 1-gap case and still NP-hard in the gapless case"). The problem was first
formalised and shown intractable in general by Lippert, Schwartz, Lancia, Istrail, *Briefings in
Bioinformatics* 3(1):23–31, 2002, doi:10.1093/bib/3.1.23. Paired-end reads are exactly the one-gap
case.

Two things make the hardness irrelevant in practice. The matrix is very sparse — a fragment covers a
few sites — and at any one site only the fragments that cover it interact, and that number is the
read depth. The two production tools exploit this in different ways.

### 1.3 WhatsHap: exact dynamic programming, exponential only in depth

**What it does.** Walk the heterozygous sites left to right. At each site, consider every way of
splitting the fragments that cover it into two groups (haplotype 1 / haplotype 2); the cost of a
split is the summed quality of the fragment alleles that disagree with their group's majority. Carry
forward the cheapest cost for each split that is consistent with the split at the previous site, and
trace back at the end. With c fragments covering a site there are 2^c splits, so the time is
O(2^c · M) for M sites (Garg et al. 2016, quoted above: "O(2^c·M) time, where M is the number of
variants to be phased and c is the maximum physical coverage"). The original description is
Patterson et al., "WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads",
*J Comput Biol* 22(6):498–509, 2015, doi:10.1089/cmb.2014.0157 (its full text was not open to this
note; the algorithm statements here are taken from Garg 2016 and from the code).

**In the code.** The column loop enumerates every bipartition of the fragments active at the site,
in Gray-code order so that only one fragment changes side per step
(`repos/whatshap/src/pedigreedptable.cpp:238-251`). A fragment's phred score is added to the cost of
the *opposite* allele of whichever side it is placed on
(`repos/whatshap/src/pedigreecolumncostcomputer.cpp:64,67`), and the column's cost is the best
genotype-consistent allele assignment over the two sides
(`repos/whatshap/src/pedigreecolumncostcomputer.cpp:101-115`). The same table serves the pedigree
case with a recombination cost added per bit that flips between columns
(`repos/whatshap/src/pedigreedptable.cpp:286-290`); a single sample is the degenerate family.

**Read downsampling is what makes 2^c affordable.** The DP runs on at most
`--internal-downsampling` fragments per site, default 15
(`repos/whatshap/whatshap/cli/phase.py:1066-1069`), hard-capped at 23
(`repos/whatshap/whatshap/cli/phase.py:1181-1182`), with a warning above 23 including pedigree
members (`:789-793`). The help text says why: "Higher values increase runtime *exponentially* while
possibly improving phasing quality marginally." Selection is a greedy priority queue that scores a
fragment by (new variants it covers that no chosen fragment covers, minus the variants inside its
unsequenced gap; then total variants; then its minimum quality)
(`repos/whatshap/whatshap/readselect.pyx:55-60`, entry `:240`), with a "bridging" pass that adds
fragments needed to keep blocks connected.

**Qualities.** The weight is the per-allele phred score stored on the fragment
(`repos/whatshap/src/read.h:57-58`). In this version of the code the score is a constant 30 for
the default edit-distance allele detection (`repos/whatshap/whatshap/variants.py:843`; the BAM base
qualities are commented out at `:821-824`), and the project's own notes say so: "does not allow us
to make use of base qualities (in fact, the weighted algorithm degenerates into an unweighted
one)" (`repos/whatshap/doc/notes.rst:27-30`). Mapping quality is a hard filter, default 20
(`repos/whatshap/whatshap/cli/phase.py:1070-1071`, applied at `whatshap/variants.py:398`); it never
enters the cost. Mates and supplementary alignments of one read name are merged into one fragment
before selection; where they overlap and agree the qualities are added, where they disagree the
higher-quality allele is kept (`repos/whatshap/whatshap/variants.py:975-984`).

**Blocks.** A block is a connected component of the graph whose vertices are heterozygous sites and
whose edges are "a read exists that covers both"; the block is named by the position of its
leftmost site (`repos/whatshap/whatshap/cli/phase.py:79-81`; union-find whose representative is the
minimum, `repos/whatshap/whatshap/graph.py:42-47`). It is written as `PS = leftmost position + 1`
(1-based) with `GT` phased (`repos/whatshap/whatshap/vcf.py:1141-1145`). A block breaks only where
no *selected* fragment spans two consecutive heterozygous sites — downsampling can therefore break a
block the full read set would have joined, which is what the bridging pass guards against. There is
no confidence-based splitting; WhatsHap writes no `PQ` ("WhatsHap currently does not add this
tag", `repos/whatshap/doc/guide.rst:343-344`), although a per-site cost difference is computed and
discarded (`repos/whatshap/src/pedigreecolumncostcomputer.cpp:163-165`).

### 1.4 HapCUT2: a likelihood, and a greedy max-cut over it

**What it does.** Instead of counting corrected bases, HapCUT2 scores a haplotype pair by the
probability of every fragment given it, with each allele call's error probability from its base
quality, and then improves the haplotype by repeatedly finding a set of sites whose phase, flipped
together, raises that probability most (Edge, Bafna, Bansal, "HapCUT2: robust and accurate
haplotype assembly for diverse sequencing technologies", *Genome Res* 27(5):801–812, 2017,
doi:10.1101/gr.213462.116; the paper states "Rather than the MEC criterion, HapCUT2 uses a
haplotype likelihood model for sequence reads" and gives the fragment probability as its
Equation 1). The greedy procedure descends from the original HapCUT (Bansal & Bafna,
*Bioinformatics* 24(16):i153–i159, 2008, doi:10.1093/bioinformatics/btn298), which phrased MEC as
finding max-cuts "in certain graphs derived from the sequenced fragments".

**Qualities, in the code.** A fragment's likelihood accumulates, per allele call, log10 of the
error probability if the call disagrees with the haplotype and log10(1 − error) if it agrees
(`repos/HapCUT2/hapcut2-src/frag_likelihood.c:15-24`; the two numbers come from the phred
character at `repos/HapCUT2/hapcut2-src/readinputfiles.c:140,143`), and the fragment's total is the
sum over both haplotype assignments, `addlogs(p0, p1)`
(`repos/HapCUT2/hapcut2-src/frag_likelihood.c:65`). Mapping quality enters at fragment extraction:
`extractHAIRS` discards reads below MAPQ 20 (`repos/HapCUT2/hairs-src/extracthairs.c:23`, filter at
`:212`), drops bases below quality 13 (`:22`), and then **caps each allele's base quality at the
read's mapping quality** (`repos/HapCUT2/hairs-src/parsebamread.c:27-32`) — so a MAPQ-25 read can
never be more than a phred-25 witness. Mates are one fragment if they map to the same contig with
an insert no longer than `--maxIS`, default 1,000 bp (`repos/HapCUT2/hairs-src/extracthairs.c:24`,
decision at `:244`); overlapping mates that disagree at a site contribute nothing at that site, and
agreeing ones keep the higher quality (`repos/HapCUT2/hairs-src/hapfragments.c:164,203-215`).

**The cut.** The graph has one vertex per heterozygous site and an edge for every pair of sites
one fragment covers (`repos/HapCUT2/hapcut2-src/variantgraph.c:198-222`; for long reads only
adjacent sites are joined). An edge's weight is the log-odds that the two alleles are in cis versus
in trans given the two base error probabilities
(`repos/HapCUT2/hapcut2-src/variantgraph.c:29-34`). The greedy cut keeps, per fragment, four partial
likelihoods (fragment given the current haplotype, its complement, and the two after the flip) and
adds sites to the flipped side by the change in total likelihood
(`repos/HapCUT2/hapcut2-src/maxcut_lr.c:5-9`, loop at `find_maxcut.c:331-365`). Each block gets
`N/10` restarts capped by `--maxcut_iter` (default 10,000) plus five seeded from the heaviest edges
(`repos/HapCUT2/hapcut2-src/find_maxcut.c:258-260,270`; defaults at `optionparser.c:16-18`), and a
block stops when five rounds pass without improvement. The paper gives the complexity as
O(T · M · (N log N + N · d · V²)) with N sites, d coverage, V sites per fragment — quadratic in
sites per fragment, which is why the long-read mode joins adjacent sites only.

**Blocks, and what breaks them.** A block is a connected component of the fragment graph, named by
the index of its first site (`repos/HapCUT2/hapcut2-src/variantgraph.c:239`); in the VCF `PS` is
that site's position (`repos/HapCUT2/hapcut2-src/output_phasedvcf.c:80`). Two post-processing steps
then use the likelihood: **pruning** blanks a site whose posterior of being correctly phased is below
0.8 — phred 6.98 — (`repos/HapCUT2/hapcut2-src/optionparser.c:21,108-109,197`), and **splitting**
computes, at each site, the posterior that a switch error starts there and cuts the block where the
no-switch posterior falls below 0.8 (`repos/HapCUT2/hapcut2-src/post_processing.c:139-159`). The
per-site confidence is written as `PQ`, the phred of the posterior that "this variant is incorrectly
phased relative to the haplotype", capped at 100 (`repos/HapCUT2/hapcut2-src/output_phasedvcf.c:7,
81-83`), together with `PD`, the number of informative fragments at the site (`:84`).

**Two code facts that matter to anyone reproducing HapCUT2 numbers.** The block-splitting switch
(`--split_blocks`) is commented out of the option parser (`repos/HapCUT2/hapcut2-src/optionparser.c:
110-111`), so splitting only runs in `--error_analysis_mode`; and `--max_IS` is unreachable because
its alias `--mi` collides with `--max_iter` (`optionparser.c:128` vs `:152`).

### 1.5 What separates the two

| | WhatsHap | HapCUT2 |
|---|---|---|
| objective | weighted MEC, solved exactly on the downsampled reads | fragment likelihood, improved greedily |
| depth handling | downsample to 15 fragments a site, hard cap 23 | none; all fragments, time grows with depth and with sites per fragment |
| base quality | weight, but constant 30 in this code | per-call error probability |
| mapping quality | filter at 20 | filter at 20, then caps each base quality |
| block break | no spanning selected read | no spanning fragment; plus posterior-based pruning (and splitting, when enabled) |
| confidence out | none (`PQ` not written) | `PQ` per site, `PD` depth |

Both name a block by its first site and both need heterozygous *calls* as input; neither calls
variants. Neither weighs the chance that a read is mapped to the wrong place as anything but a
threshold — a mis-mapped read pair that passes MAPQ 20 links two sites with full weight.

---

## 2. What short reads give at low coverage — arithmetic, not measurement

**This section is this note's own calculation.** It was run with
`tmp/phasing_research/read_backed/span_table.py` (`uv run`), and nothing in it was measured on data.
The assumptions:

- Reads are 150-base paired-end. Insert (outer fragment) length is Normal with the stated mean and a
  50-base standard deviation, floored at the read length. Both mates align fully.
- Fragment starts are Poisson along the genome, at a rate that gives the stated coverage.
- Heterozygous sites are Poisson along the genome at the stated density, and **both are already
  known** — called from the cohort, not from this sample's own reads. At 1–3× a single sample
  cannot call its own heterozygotes; the only realistic setting is that the cohort supplies the
  site list.
- A fragment **spans** a pair of sites when both fall in its sequenced bases: both in one mate, or
  one in each. Base errors and mis-mapping are ignored (they only lower the numbers).
- Under these assumptions the number of fragments spanning a pair at distance D is Poisson with
  mean (fragments per base) × (number of start offsets from which both sites are sequenced); the
  table averages over D and over the insert length.

Coverage is *base* coverage (a 3× sample has 3 sequenced bases a position). Three densities: one
heterozygote per kb (an outbred plant), one per 1.5 kb (a human of European ancestry), one per 10 kb
(illustrative for inbred material with residual heterozygosity). "P(pair linked)" is the chance that
an adjacent heterozygous pair is spanned by at least one fragment; "E[spanning]" is the expected
number of fragments spanning it; "hets/block" is the mean number of sites in a block (a geometric
chain, singletons included); "block bp" is the mean block span, singletons counting as zero.

| insert | het/kb | cov | P(pair linked) | E[spanning] | hets/block | block bp |
|---:|---:|---:|---:|---:|---:|---:|
| 350 | 1.0 | 1 | 0.10 | 0.13 | 1.11 | 14 |
| 350 | 1.0 | 3 | 0.20 | 0.39 | 1.24 | 34 |
| 350 | 1.0 | 10 | 0.27 | 1.31 | 1.37 | 57 |
| 350 | 1.0 | 30 | 0.29 | 3.93 | 1.40 | 65 |
| 350 | 0.67 | 1 | 0.07 | 0.09 | 1.07 | 9 |
| 350 | 0.67 | 3 | 0.14 | 0.27 | 1.16 | 22 |
| 350 | 0.67 | 10 | 0.19 | 0.90 | 1.23 | 37 |
| 350 | 0.67 | 30 | 0.20 | 2.73 | 1.25 | 42 |
| 350 | 0.1 | 3 | 0.02 | 0.04 | 1.02 | 3 |
| 350 | 0.1 | 30 | 0.03 | 0.44 | 1.03 | 6 |
| 500 | 1.0 | 1 | 0.10 | 0.12 | 1.11 | 20 |
| 500 | 1.0 | 3 | 0.20 | 0.37 | 1.25 | 51 |
| 500 | 1.0 | 10 | 0.30 | 1.24 | 1.42 | 94 |
| 500 | 1.0 | 30 | 0.33 | 3.74 | 1.50 | 115 |
| 500 | 0.67 | 3 | 0.14 | 0.27 | 1.17 | 35 |
| 500 | 0.67 | 30 | 0.24 | 2.62 | 1.31 | 75 |
| 500 | 0.1 | 30 | 0.04 | 0.44 | 1.04 | 11 |

Three things the table says:

1. **The ceiling is geometric, not a matter of depth.** At one heterozygote per kb, the chance that
   the *next* heterozygote lies within one 350-base insert is about 30%, and 70% of adjacent pairs
   can never be linked by any number of short reads. Going from 10× to 30× raises P(pair linked)
   from 0.27 to 0.29; going from 3× to 10× raises it from 0.20 to 0.27. Depth buys *confidence* (the
   expected number of spanning fragments goes 0.4 → 1.3 → 3.9) more than it buys links.
2. **At 3×, most links rest on a single fragment.** With E[spanning] ≈ 0.4 and P(linked) ≈ 0.2, a
   linked pair has on average two spanning fragments, but the distribution is Poisson-shaped: about
   half the linked pairs have exactly one. One fragment with one base error at either site flips
   the link. This is the regime in which HapCUT2's per-site posterior would prune most links below
   phred 7, and WhatsHap's constant weight 30 would report them as certain.
3. **Blocks are one or two heterozygotes long.** Mean block size is 1.2–1.5 sites at every depth for
   short inserts; the mean linked-pair distance is about 140 bp (350-base inserts) to 230 bp
   (500-base inserts). This is the same picture the literature reports for 45× Illumina: HapCUT
   blocks with a length-weighted N50 (QAN50) of 1,035 bp at 148-base reads, against 107,619 bp for
   50× PacBio (Choi et al. 2018, Table 1, cited in §3.4).

A rough check against a measurement. Delaneau et al. 2013 report that in 379 European 1000 Genomes
samples at about 4×, 33.8% of heterozygous sites "were covered by a sequence read that also covered
another heterozygous site in the same individual", and in simulation with 100-base reads and
300-base inserts about 17% at 1× and about 50% at 20×. Their quantity is per *site* and counts a
link to either neighbour (or any farther site), so it is roughly 1 − (1 − p)² of the per-pair
figure here: at 3× and 0.67 het/kb that gives 1 − 0.86² ≈ 0.26, and at 30× about 0.36. Same order,
lower by a third, and the difference is in the direction the assumptions predict (their samples had
more heterozygotes per kb than 0.67 and a mix of libraries). The table should be read as an order of
magnitude, not a prediction.

---

## 3. How population phasing consumes read phase

### 3.1 SHAPEIT2's phase-informative reads: a likelihood multiplied into the sampler

Delaneau, Howie, Cox, Zagury, Marchini, "Haplotype estimation using sequencing reads", *Am J Hum
Genet* 93(4):687–696, 2013, doi:10.1016/j.ajhg.2013.09.002 (SHAPEIT2 is closed source; this is from
the paper). SHAPEIT2's statistical phaser samples, for each individual, a haplotype pair that is
consistent with its genotypes and probable under the population model (the other individuals'
current haplotypes, through a hidden Markov model). The read extension adds a second factor: for a
candidate haplotype pair X, the probability of the reads given X, where each read's alleles are
compared with the haplotype it is assigned to and each base contributes 1 − Q if it matches and Q
if it does not, Q being the base's error probability. The sampler draws haplotype pairs
proportional to the product P(X | H) × P(X | R) — population model times read likelihood — with a
forward–backward pass for each factor, segment by segment. Mates are joined into one read for this
purpose. The read term is soft: a read that disagrees with the population's favoured phase is charged
its base qualities, not vetoed.

Reported gain (their Figure 2 and text; real Illumina data on a trio phased with 1000 Genomes
Europeans as the cohort): without reads, mean distance between switch errors 270.1 kb and SER
0.747%; with the trio's own 5× reads 295.3 kb and 0.686%; with 20× trio reads plus about 4× reads
for the 1000 Genomes samples 328.6 kb and 0.616% — "22% increase in mean distance between switch
errors". Simulated 100-base paired reads at 300-base inserts show diminishing returns above 5×,
mixed insert sizes (300/500/1,000) beat any single size, and long reads are where the method changes
character: 5 kb reads at 4% error and 5× raised the mean switch distance from about 204 kb to
584 kb, and 20 kb reads at 20× cut the switch error at singletons from 49.3% to 2.9%.

### 3.2 SHAPEIT4 `--use-PS`: phase sets as a soft mask on the HMM transitions

SHAPEIT4 replaces the per-read likelihood with something a caller can hand over in a VCF: the
`PS` field. Delaneau, Zagury, Robinson, Marchini, Dermitzakis, "Accurate, scalable and integrative
haplotype estimation", *Nat Commun* 10:5436, 2019, doi:10.1038/s41467-019-13225-y (Methods):
"SHAPEIT4 accommodates phase sets in a probabilistic manner … it assumes an error model in which
paths in the genotype graphs that are consistent with the phase sets receive more weight than paths
that are not", and "The distribution P(D|R) is now controlled by a parameter that defines the
expected error rate in the phase sets and not directly computed from the sequencing reads. By
default, we assume an error rate of 0.0001." The sampled haplotypes "are not necessarily consistent
with all the phase sets which allows correcting phasing errors in sequencing reads".

In the code the option takes that error rate as its argument
(`repos/shapeit4/src/phaser/phaser_parameters.cpp:39`: `("use-PS", bpo::value<double>(), "Informs
phasing using PS field from read based phasing")`). Reading: the `PS` integer is fetched per record
(`repos/shapeit4/src/io/genotype_reader2.cpp:53`), and a heterozygous genotype joins a phase set
only if its `GT` is phased (`|`) and a `PS` is present; otherwise its code is 0, "no set"
(`repos/shapeit4/src/io/genotype_reader2.cpp:62-68`). What is stored per heterozygous site is the
set code plus the two alleles *in the order the input VCF wrote them* — the `0|1` versus `1|0`
orientation (`repos/shapeit4/src/objects/genotype/genotype_header.h:58-67`, the per-individual
vectors at `:96-97`).

The constraint is built per segment of SHAPEIT4's genotype graph: for every candidate transition
(a choice of haplotype path through the segment), the implied haplotype is compared with each phase
set; the first site of the set fixes which haplotype the set's "first allele" is on, and every later
site of the same set must agree with that orientation or the transition is flagged inconsistent
(`repos/shapeit4/src/objects/genotype/genotype_mask.cpp:91-95`). So a phase set is exactly "these
heterozygotes have this *relative* phase"; flipping a whole set is free. Then, before each sampling
iteration, every transition probability is multiplied by (1 − error) if consistent and by *error* if
not, and renormalised (`repos/shapeit4/src/objects/compute_job.cpp:229-235`; called per individual at
`repos/shapeit4/src/phaser/phaser_algorithm.cpp:78`). With the recommended 0.0001, a path that breaks a
phase set pays a factor of 10⁴ against one that keeps it — a strong prior, not a veto. The
documentation adds two practical facts: `PS` must be an integer, since a string `PS` is silently
ignored, and it recommends WhatsHap to produce it (`repos/shapeit4/docs/index.html:288-291`).

The paper's experiment (Results, Figure 4, Table 1): one GIAB sample on chromosome 20 merged with 502
unrelated Europeans and phased against about 800,000 UK Biobank haplotypes. Without phase sets the
SER was 0.82%. Phase sets from Illumina HiSeq covered 57.9% of heterozygous genotypes and brought it
to 0.42%, "49% decrease"; PacBio phase sets covered 90.1% and gave 0.23%; 10x Genomics covered 97.7%
and gave 0.07%. The phase sets were produced by "methods such as WhatsHap" and, for 10x, Long Ranger.

### 3.3 SHAPEIT5, Beagle, Eagle, GLIMPSE, QUILT, STITCH

- **SHAPEIT5 dropped `--use-PS`.** There is no `use-PS` option and no `PS` read anywhere in the
  SHAPEIT5 tree (grep of `repos/shapeit5/` for `use-PS`, `"PS"` and `bcf_get_format_int32` returns
  nothing; `phase_common`'s full option list is
  `repos/shapeit5/phase_common/src/phaser/phaser_parameters.cpp:33-75`). What survives is the
  `--scaffold` mechanism, a *hard* constraint: scaffold heterozygotes are fixed in the genotype
  graph and never resampled (`repos/shapeit5/phase_common/src/io/genotype_reader/
  genotype_reader_reading.cpp:198,216`; SHAPEIT4's equivalent at `repos/shapeit4/src/objects/
  genotype/genotype_build.cpp:49-52`). A scaffold is meant for trio- or panel-phased haplotypes, not
  for read phase, and `phase_rare` requires one.
- **Eagle2 has a hidden, experimental `--usePS`.** It is declared under "Hidden options" as
  `"use FORMAT:PS phase constraints in target VCF: 1=soft, 2=harder"`
  (`repos/Eagle/src/EagleParams.cpp:123-129`). Its convention is not the VCF one: the `PS` integer
  is the *position of an anchor heterozygote*, and its **sign** carries the relative phase
  (`repos/Eagle/src/SyncedVcfData.cpp:360-373`); the constraints are applied in the PBWT haplotype
  search (`repos/Eagle/src/EaglePBWT.cpp:253-258,336-342`). Eagle never reads the `|` of a target
  genotype (its target parser reduces `GT` to a dosage). Treat it as unsupported.
- **Beagle 5.5 reads no `PS`.** It parses the `|` separator, but the only algorithmic use is an
  all-or-nothing test: if *every* genotype in the window is phased, Beagle copies the input phase to
  the output without phasing (`repos/beagle_src/src/main/Main.java:143`); otherwise every
  heterozygote is treated as unphased (`repos/beagle_src/src/phase/PbwtPhaser.java:224-226` builds
  the list from `a1 != a2` alone).
- **GLIMPSE2, QUILT/QUILT2 and STITCH read no `PS` and no target-genotype phase.** GLIMPSE2's target
  input is `GL`/`PL` only (`repos/GLIMPSE/phase/src/io/genotype_reader.cpp:467`); QUILT and STITCH
  take reads directly (§4.1) and have no target-VCF input. A grep for `"PS"` and `phase_set` across
  `repos/QUILT`, `repos/GLIMPSE` and `repos/beagle_src/src` returns nothing.

So among the maintained statistical phasers, **the only one that consumes a caller's `PS` today is
SHAPEIT4**, and its successor removed the feature. The reference-panel imputers consume reads, not
phase sets.

### 3.4 What the literature reports on adding reads to population phasing

- Delaneau 2013 (§3.1): 22% longer mean distance between switch errors on real data at 5–20× short
  reads; long reads are a different regime.
- Delaneau 2019 (§3.2): SER 0.82% → 0.42% with Illumina phase sets covering 58% of heterozygotes;
  → 0.07% with 10x phase sets covering 98%, with 800k reference haplotypes behind the population
  model.
- Choi, Chan, Kirkness, Telenti, Schork, "Comparison of phasing strategies for whole human
  genomes", *PLoS Genet* 14(4):e1007308, 2018, doi:10.1371/journal.pgen.1007308 (Tables 1 and 2;
  one deeply sequenced individual, 45× Illumina 148-base paired-end and 50× PacBio, population
  phasing with SHAPEIT2 against the Haplotype Reference Consortium panel): read-based alone, HapCUT
  on Illumina phased 83.5% of SNVs with SER 0.096% and QAN50 1,035 bp; on PacBio SER 1.015% and
  QAN50 107,619 bp. Population alone: SER 0.297%, QAN50 1.90 Mb. Population plus Illumina reads:
  SER 0.214%, QAN50 2.38 Mb; plus PacBio reads: 0.139% and 3.52 Mb; plus parental genotypes:
  0.025% and 18.9 Mb. Their summary: "up to 85% improvement in QAN50 and 53% reduction in SER".

Two readings of these numbers. First, short reads on their own are *accurate* (SER 0.1% at 45×) but
*short* (1 kb blocks), and the population model is the reverse (megabase blocks, SER 0.3%); adding
the reads mostly fixes the population model's local mistakes, worth a 30% SER reduction at 45×.
Second, every measurement above was made at 45× or more, on a human, with a large panel. Nothing in
the literature reached here reports the gain at 1–3× — where §2 says a fifth of adjacent pairs are
linked, mostly by one fragment.

---

## 4. The reverse direction: reads helping imputation at low coverage

### 4.1 QUILT: reads, not sites, as the unit of evidence

Davies et al., "Rapid genotype imputation from sequence with reference panels", *Nat Genet*
53:1104–1111, 2021, doi:10.1038/s41588-021-00877-0. QUILT's argument for a per-read model is that
site-level genotype likelihoods throw away the correlation a read induces between sites:
"sequencing reads are typically 100-250 bases long, may be paired, and with long read sequencing,
may be many thousands of bases long. As such, information at nearby SNPs may not be independent, by
coming from the same read(s)." Its model assigns each read to one of the two haplotypes by Gibbs
sampling — "for read indexed by v, h_v is drawn equally from {1,2}" at initialisation — and runs a
haploid copying HMM on each haplotype's reads: "If the maternal or paternal origin of each
sequencing read can be determined, imputation becomes much simpler, because maternal and paternal
reads can be separated and haploid imputation performed". Mates are one unit: "Reads were grouped
into molecules, either using paired read information, or using the BX tag from haplotagging."

The code confirms the unit and the qualities. Reads are pulled from the BAM with each allele call's
base quality capped at the read's mapping quality and signed by allele (negative for REF, positive
for ALT) (`repos/STITCH/STITCH/src/bam_access.cpp:326-341`, shared with STITCH); all records of one
read name are concatenated into one read unit, dropped if their span exceeds `iSizeUpperLimit`
(`repos/STITCH/STITCH/src/functions.cpp:320-325,454`; QUILT's default 10⁶,
`repos/QUILT/QUILT.R:227-231`); and the Gibbs step rewrites the read's haplotype label
(`repos/QUILT/QUILT/src/gibbs-nipt.cpp:1089`). The per-read assignment is read-backed phasing and
imputation in one model: the read's alleles at the sites it covers are what decides which haplotype
it is copied onto.

Reported at 0.5× on NA12878 (Figure 3A): r² at rare SNPs 0.678 for Illumina, 0.669 for ONT, 0.709
for haplotagged reads; at 0.1×, 0.416 (Illumina) versus 0.460 (haplotagged). The gain from longer
molecules at the same depth is a few hundredths of r², which is the phase information those
molecules add over 150-base pairs. QUILT's handling of reads is covered in the low-coverage
imputation note; the point here is that **the modern low-coverage imputer already does read-backed
phasing internally, from the BAM, and does not want a `PS` from the caller.**

### 4.2 freebayes: phase by emitting one allele for the whole stretch

freebayes phases nearby variants without any `PS`: "FreeBayes is capable of calling variant
haplotypes shorter than a read length where multiple polymorphisms segregate on the same read. The
maximum distance between polymorphisms phased in this way is determined by the --max-complex-gap,
which defaults to 3bp." (`freebayes/src/Parameters.cpp:44-48`; `--haplotype-length` is the same
option, `:204-207,542-543`; default at `:431`). Per read, adjacent non-reference alleles separated
by at most that many reference bases are merged into one *complex* allele
(`freebayes/src/AlleleParser.cpp:1013-1058`, `RegisteredAlignment::clumpAlleles`), so the record
written has a multi-base REF and ALT and `INFO/TYPE=complex` or `mnp`
(`freebayes/src/ResultData.cpp:368-380`). The phase is therefore carried *inside the allele*: a
sample heterozygous for `ACG→ATA` has both substitutions on one haplotype, and the reads that showed
only one of them are partial observations of it. freebayes writes no `PS` and no `|` (grep for
`"PS"`/`phased` in `freebayes/src/` finds only the help text). `--no-complex` (`-u`) switches the
merging off (`freebayes/src/Parameters.cpp:768`). The reach is at most a read length — the help
suggests up to "half the read length" — and the price is a multi-allelic record with as many
complex alleles as distinct read-level combinations.

### 4.3 GATK HaplotypeCaller: `PID`/`PGT` from the assembled haplotypes

HaplotypeCaller does not look at read pairs at all for phasing; it phases from the local assembly.
After genotyping an active region it maps each called variant to the set of assembled haplotypes
carrying its alternate allele (`gatk/.../haplotypecaller/AssemblyBasedCallerUtils.java:728-739`,
`phaseCalls`). Two variants join a phase group when their alternate alleles sit on exactly the same
haplotypes (in cis) or on complementary, non-overlapping sets (in trans)
(`AssemblyBasedCallerUtils.java:790-908`, `constructPhaseSetMapping`); a variant that would belong to
two incompatible groups aborts phasing for the region (`:843`, `:881`). Each group is written with
`PID = position_ref_alt` of its leftmost variant, `PGT = 0|1` or `1|0` for the alternate's side,
`GT` marked phased, and `PS = position of the first variant`
(`AssemblyBasedCallerUtils.java:938-942,983-1004`; constants `PGT`/`PID` at
`gatk/.../utils/variant/GATKVCFConstants.java:142-143`; header lines at
`HaplotypeCallerEngine.java:552-554`). Only biallelic variants are eligible (`:753-756`), and the
feature is on by default (`--do-not-run-physical-phasing`, `HaplotypeCallerArgumentCollection.java:
249-253`). Because the input is an assembly of the region's reads, mates contribute only as separate
reads, and the reach is an active region — a few hundred bases.

The two callers thus phase the same neighbourhood — one read length — by two conventions: freebayes
by changing the allele, GATK by tagging genotypes. Neither reaches an insert length.

### 4.4 Linked and long reads, as context

The linkage in §2 is bounded by the molecule, and other molecules move that bound: 10x Genomics
linked reads (barcoded short reads from a ~50 kb molecule; HapCUT2's `--10X` mode links reads by the
`BX` tag within 20 kb, QUILT's `bxTagUpperLimit` within 50 kb), PacBio and ONT reads of 10–20 kb
(HapCUT2's `--pacbio`/`--ont` modes with re-alignment at base quality 4), and Hi-C, which links
sites megabases apart. The SHAPEIT4 numbers in §3.2 (58% of heterozygotes in a phase set with
Illumina, 90% with PacBio, 98% with 10x) and Choi 2018's block lengths (1 kb versus 108 kb) are the
measured size of the difference. Everything after this line is about short paired-end reads, which
is what ng's inputs are.

---

## 5. The VCF conventions

VCF 4.3 §1.6.2, Table 2 (reserved genotype keys; `raw.githubusercontent.com/samtools/hts-specs/
master/VCFv4.3.tex`):

- **`GT`** — alleles separated by `/` (unphased) or `|` (phased). Within one sample, the first allele
  of every phased genotype in the same phase set is on the same haplotype, the second on the other.
- **`PS`** — "Phase set, defined as a set of phased genotypes to which this genotype belongs." A
  non-negative 32-bit integer, per sample. Phased genotypes of one sample on one chromosome with the
  same `PS` are in one set; genotypes without `PS` are assumed to be one set. The spec's recommended
  convention is the position of the first variant in the set, which is what WhatsHap
  (`repos/whatshap/whatshap/vcf.py:1141`), HapCUT2 (`repos/HapCUT2/hapcut2-src/output_phasedvcf.c:80`)
  and GATK (`AssemblyBasedCallerUtils.java:942`) all do — so tools agree on numbering by
  construction, not by negotiation, and a downstream phaser only needs the *equality* of `PS`
  values, never their value (SHAPEIT4 re-codes them to dense integers on input,
  `repos/shapeit4/src/io/genotype_reader1.cpp:64-81`).
- **`PQ`** — "Phasing quality, the phred-scaled probability that alleles are ordered incorrectly in
  a heterozygote (against all other members in the phase set)." HapCUT2 writes it; WhatsHap does not;
  no phaser reads it.
- **`HP`** — not a reserved key. It is the GATK `ReadBackedPhasing` convention (`HP = block-haplotype`
  per allele, `GT` left unphased) that WhatsHap can emit with `--tag HP`
  (`repos/whatshap/whatshap/cli/phase.py:1048-1050`). Nothing in this survey reads it.
- **`PID`/`PGT`** — GATK-private (§4.3); `PS` is written alongside.

What a caller must write for SHAPEIT4 to use its read phase: at every heterozygote it phased, `GT`
with `|` in the phase set's orientation and an integer `PS`; at every other heterozygote `GT` with
`/` and `PS` missing (`.`). Both conditions are checked on input
(`repos/shapeit4/src/io/genotype_reader2.cpp:62-68`): a `|` without `PS`, or a `PS` on a `/`
genotype, joins no set. Homozygous sites carry no `PS`.

---

## 6. What ng has, and what a read-backed phasing step would be

### 6.1 The evidence ng already holds

At every locus, for every sample, ng knows which reads showed which allele: each
`SequenceObservation` — one distinct (bases, witness, read group) at the locus — carries the ids of
the reads folded into it (`src/ng/locus_generation/mod.rs:369`, `chain_ids`). A chain id is one
`u64` per read, or per read *pair* with the mates collapsed onto one id via the walker's pending-mate
map (`src/ng/locus_generation/pileup/chain_id_allocator.rs:4-5`), so the unit is already the fragment
WhatsHap and HapCUT2 construct from the BAM, without insert-size gating: two mates are one id however
far apart they map. The field's own doc names the first of its two purposes as "which haplotype a
read came from, which is what lets a later step chain observations at neighbouring loci into one"
(`src/ng/locus_generation/mod.rs:353-354`). Alongside the ids the observation keeps per-read
*moments* — `num_obs` (`:322`), `q_sum`, the summed per-read log error (`:333`), `mapq_sum` and
`mapq_sum_sq` (`:336-338`) — not a quality per read.

### 6.2 Turning that into phase evidence between two heterozygous loci

Take one sample and two loci A and B within an insert length where the sample is heterozygous.
At A the sample's reads are partitioned by allele: sets S(A, a₁), S(A, a₂); at B, S(B, b₁),
S(B, b₂). A read in S(A, a₁) ∩ S(B, b₁) or in S(A, a₂) ∩ S(B, b₂) says *cis*: a₁ with b₁. A read in
S(A, a₁) ∩ S(B, b₂) or S(A, a₂) ∩ S(B, b₁) says *trans*. Reads in only one locus's sets say nothing.
The counts n_cis and n_trans are the entire read-backed signal for the pair, and the intersection is
a merge of two sorted id lists (ids are allocated monotonically, and the lists are short at any
depth ng serves).

**Weighting.** The two counts are not equal witnesses:

- *Base error at either site* flips a read's vote. HapCUT2's fragment likelihood is the right
  shape: a read in the cis set contributes (1 − e_A)(1 − e_B) to the cis hypothesis and roughly
  e_A + e_B to trans, where e is the error probability of the read's base at that site. ng does not
  keep a per-read error at the observation, only `q_sum` over the reads of that observation
  (`src/ng/locus_generation/mod.rs:333`); the mean, `q_sum / num_obs`, is the per-read value
  available without changing the observation's shape. At 3×, where a link often rests on one read,
  the difference between a phred-20 and a phred-35 witness is the difference between a 1-in-100 and a
  3-in-10,000 chance that the link is wrong — worth carrying.
- *Mis-mapping* is worse than base error because it flips *both* alleles coherently and is not
  independent across the reads of a repeat. HapCUT2's cap of base quality at MAPQ is the cheapest
  defensible treatment; the observation's `mapq_sum` gives the mean, and the per-read MAPQ would
  need to be kept if the cap is to be per read.
- *Both sites must be genuine heterozygotes in this sample.* At 3× a sample's own genotype is
  uncertain; a read pair linking a₁–b₁ at a site where the sample is really homozygous for a₁ is
  not phase evidence. The link's weight has to be multiplied by the posterior that the sample is
  heterozygous at both — which is exactly the number the cohort genotyper produces and the reason
  SHAPEIT4 only joins a phase set at sites already called heterozygous.

Put together, the evidence per adjacent pair is a log-odds cis versus trans, and a block is a chain
of pairs whose log-odds exceeds a threshold. That is HapCUT2's post-processing without the max-cut:
with blocks of 1.2–1.5 sites (§2) there is no combinatorial problem to solve, because a block almost
never has a third site to make the pairwise answers disagree. WhatsHap's dynamic programming and
HapCUT2's cut earn their keep at 10 kb reads; for short reads at any depth the pairwise likelihood
*is* the algorithm.

**How many pairs, from §2.** At one heterozygote per kb and 350-base inserts: at 3×, one adjacent
pair in five is linked, on average by 0.4 fragments — so a 63-sample cohort with 100,000 heterozygous
sites per sample would produce about 20,000 linked pairs per sample, half of them on one read pair.
At 30×, 29 in 100 pairs linked, by 3.9 fragments on average; a link then has enough witnesses for
its own error rate to be estimated from the disagreeing minority. In both cases the *same* 70% of
pairs are out of reach.

### 6.3 Reach: where chain ids are comparable

Two documents bound what a phasing step may compare:

- **Within one sample only.** "a chain id is comparable within one sample and means nothing without
  it" (`src/ng/run/cohort_merge/build.rs:2109`); each sample's ids are its own file's, allocated from
  zero (`src/ng/locus_generation/pileup/chain_id_allocator.rs:6-7`). This is no restriction:
  read-backed phasing is by definition within a sample.
- **Within one walked region, and by guarantee within one locus.** The id names "a read within one
  walk, and a read that straddles the boundary between two walked regions is met twice and named
  twice. Nothing downstream links a read across such a boundary, because a segment is never cut and
  no locus crosses one" (`src/ng/locus_generation/mod.rs:361-365`, quoting
  `doc/devel/ng/spec/run_streaming.md` §4.3). The encoding spec promises even less — "the stored form
  is free to be anything that preserves equality over a window of one locus"
  (`doc/devel/ng/spec/psp_chain_id_encoding.md` §1.1) — although the encoding as built restates the
  full live set at every block and keeps the real ids (`psp_chain_id_encoding.md` §4, "Random access
  is preserved by restating, not by chaining across blocks").

So the reach of read-backed phasing inside ng is **one segment**: two heterozygous loci in the same
generic segment can be linked; two on either side of a repeat-tract boundary cannot, even when a read
pair spans the tract. `run_streaming.md` §4.3 says what the alternative would be: "it must be a pass
over the emitted records, never a coupling between in-flight segments". Since the segments at the
routing floor are "probably kilobases long" (§4.4, unmeasured) and the linkable distance is an insert
length, the loss is the pairs whose two sites straddle a boundary — small, but it should be counted
when the segment lengths are measured (open question 1 of `run_streaming.md` §11). A phasing step
that wants more than a locus must also promote the within-locus equality guarantee of
`psp_chain_id_encoding.md` §1.1 to a within-segment one; the code already satisfies it, the spec does
not yet say so.

### 6.4 What the merge does with the ids today, and what retaining them costs

The cohort merge consumes the ids to elongate a read across a locus's records and to tell "covered
and reference" from "never here", then hands calling only per-(allele, read group) counts:
`AlleleSupport { num_reads, num_fwd, q_sum, mapq_sum, mapq_sum_sq, placed_left }`
(`src/ng/run/cohort_merge/build.rs:1473-1495`). No id survives into calling, so a phasing step cannot
sit after the genotyper as a per-record filter without a change upstream: the ids must be retained,
per sample and per allele, for every locus that might turn out heterozygous, until every locus within
an insert length has been genotyped.

The cost is small at any point in the range, because it is per *locus*, not per position, and only
over a window:

- Per sample per retained locus: one id per read at the locus, 8 bytes each (4 if re-coded to an
  index within the window). At 3× that is about 3 reads → 24 bytes; at 30×, 240 bytes; at 300×,
  2.4 kB.
- Loci in the window: those within the maximum insert length (say 1 kb) — at one candidate
  heterozygote per kb, one or two; at ten per kb (a highly diverse region, or if candidates rather
  than called heterozygotes are retained), ten or twenty.
- Cohort: multiply by samples. 1,000 samples × 30× × 2 loci ≈ 0.5 MB; 3,000 samples × 300× × 20 loci
  ≈ 150 MB, and that corner is the one where per-position depth caps already apply.

Against the 500 kB-per-open-psp budget of `run_streaming.md` §7.2, the 30× cohort case adds about
one part in a thousand; the deep-and-dense corner is the only one that needs a cap on retained loci.
The retention has to happen where the ids still exist: in the merge's per-sample fold, before
`AlleleSupport` is summed, and the retained lists must be keyed by the *allele the read showed at
the cohort locus* after elongation, which is what the merge already computes.

What the step would emit is per sample: for each adjacent heterozygous pair inside a segment, a
log-odds cis/trans and the count of witnesses; from those, `GT` with `|` and an integer `PS` at the
linked sites, `/` and no `PS` elsewhere (§5), and `PQ` from the log-odds if a consumer ever wants it.
That output is exactly what SHAPEIT4 consumes (§3.2) — and, per §3.3, what no other current phaser
does.

---

## 7. Questions to discuss

1. **Is read-backed phasing worth doing at 3× short reads?** *Recommendation: not as a phasing
   product; possibly as an input to the cohort model, if that model can take it.* At 3× and one
   heterozygote per kb, §2 says one adjacent pair in five is linked, mostly by a single read pair,
   into blocks of 1.2 sites; a `PS` output at that depth would be true but nearly empty. The only
   published gain at low depth is Delaneau 2013's 22% longer switch distance, and that was with 5–20×
   reads on top of 1000 Genomes depth for the cohort. The value at 3× is inside an imputation model —
   QUILT's per-read labels are that — not in a VCF field. Trade-off: doing it inside the model means
   building the read model in ng rather than exporting; not doing it forfeits whatever a fifth of
   the pairs is worth, which nobody has measured at this depth on a plant cohort.
2. **A `PS` output for deep samples (≥20×)?** *Recommendation: yes, cheaply, once the ids are
   retained.* At 30× a link has about four witnesses and its error rate can be read off them; the
   pairwise log-odds is a few dozen lines on top of the retained lists, and the output is the standard
   the one consumer that exists (SHAPEIT4) reads. Trade-off: SHAPEIT5 dropped the feature, so the
   consumer is a tool whose successor does not want the field; the durable value is for the owner's
   own downstream, or for haplotype-allele reporting (next item).
3. **Should ng phase the way freebayes does — by emitting the read-level haplotype allele — at least
   within a read length?** *Recommendation: decide this before the `PS` question, because it changes
   the record, not a tag.* freebayes' complex alleles carry the phase in the allele and cost a
   multi-allelic record; GATK's `PID`/`PGT` keep single-site records and tag them. ng's locus shape
   is already multi-allelic and multi-base (`CandidateAlleles`), so the freebayes form is closer to
   what exists, but it also determines how the cohort model sees a nearby pair — as one locus with
   four alleles or as two loci with a phase tag. Trade-off: one allele per haplotype combination
   inflates allele counts at low depth exactly where the genotyper is weakest.
4. **Retain ids in the merge now, or when a consumer exists?** *Recommendation: now, behind a flag,
   because the cost is a few hundred kilobytes and the retention point (the per-sample fold before
   `AlleleSupport`) is the one place the ids and the elongated alleles coexist.* Trade-off: a data
   structure kept for a consumer that may not come; mitigated by making it the same list the merge
   already holds, kept one locus longer.
5. **Which quality to carry per read: none, the mean, or the read's own?** *Recommendation: the mean
   at first (`q_sum/num_obs`, `mapq_sum/num_obs`), with HapCUT2's cap of base quality at MAPQ.* A
   per-read quality would change `SequenceObservation`'s shape for every consumer; the mean loses
   only the within-observation spread, which at ng's depths is a handful of reads. Trade-off: a
   single bad read in an observation of three is hidden by the mean, and at 3× that read may be the
   only witness of a link.
6. **Promote the chain-id equality guarantee from one locus to one segment.** *Recommendation: do
   it in `psp_chain_id_encoding.md` §1.1 when a consumer appears; the code already satisfies it.*
   Trade-off: it forbids a future encoding that renumbers per block, which §4 of that spec already
   argued against.
7. **Measure rather than derive.** §2 is arithmetic. One afternoon on `benchmarks/tomato1` (63
   accessions at about 3×) and on the GIAB trio would give the measured fraction of adjacent
   heterozygous pairs with n_cis + n_trans ≥ 1, the distribution of witnesses, and the disagreement
   rate n_min/(n_cis + n_trans) at 30× — the last being the empirical link error rate the `--use-PS`
   argument wants. Recommendation: run it before deciding items 1 and 2; both hinge on numbers this
   note could only derive.
