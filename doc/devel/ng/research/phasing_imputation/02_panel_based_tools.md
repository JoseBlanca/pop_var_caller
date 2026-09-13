# How the mainstream phasing and imputation tools are built: a code-level survey

Sub-report 2 of the phasing/imputation research. Written 2026-09-12 from the shallow clones under
`tmp/phasing_research/repos/` (all `path:line` citations below are relative to that directory) and
from the papers listed in the references. The tools covered are the ones written for
**well-genotyped data** — array or high-coverage sequencing, hard genotype calls in — and the
reader should keep that in mind throughout: none of them takes a genotype likelihood. The
low-coverage imputers (GLIMPSE, QUILT, STITCH) are in sub-report 3; the few-haplotype
(segregant-population) tools are in sub-report 4.

Tools: SHAPEIT5 (with SHAPEIT4 where it differs), Eagle2, Beagle 5.5, minimac4, PBWT, IMPUTE5
(closed source; paper only).

---

## 1. Vocabulary

Genetics terms are kept brief; statistics and software terms are spelled out.

- **Phasing.** Deciding, at every heterozygous site of a sample, which of the two alleles sits on
  which of the two chromosome copies (haplotypes). A *switch error* is a place where the inferred
  haplotype swaps from the true maternal to the true paternal copy; the *switch error rate* (SER)
  is the number of such swaps divided by the number of consecutive heterozygous pairs.
- **Imputation.** Filling in genotypes at sites a sample was not typed at (or was typed badly at)
  by copying them from haplotypes that match the sample's own haplotypes nearby.
- **Reference panel.** A set of haplotypes phased and genotyped beforehand, used as the thing
  being copied from. *Panel-free* (or *within-cohort*) means the other samples of the same study
  play that role.
- **Conditioning set / states.** The haplotypes a given target haplotype is allowed to copy from
  at a given position. Every tool here shrinks this set from "all haplotypes" to a few hundred or
  thousand; how it does so is the main thing that distinguishes the tools.
- **Li and Stephens model.** The statistical picture all these tools share: a haplotype is a
  mosaic of pieces copied from other haplotypes, the copied-from haplotype changes with a
  probability tied to genetic distance (a *recombination* or *switch* parameter), and the copy is
  imperfect with a small *mismatch* (error/mutation) probability. See §2.
- **Hidden Markov model (HMM).** The machinery for the mosaic picture: the hidden state at each
  site is "which haplotype is being copied"; the *forward* pass accumulates, site by site, the
  probability of the data up to that site for every state; the *backward* pass does the same from
  the right; multiplying the two at a site gives the *posterior* probability of each state there.
  *Viterbi* finds the single most probable path instead. *Sampling* draws a path at random in
  proportion to its probability (used inside iterative schemes so each iteration sees a fresh
  draw).
- **Diploid vs haploid HMM.** A *diploid* HMM tracks a pair of copied-from haplotypes (one per
  chromosome copy), so its state space is the square of the conditioning set; a *haploid* HMM
  phases one haplotype at a time against a fixed estimate of the other.
- **PBWT (positional Burrows–Wheeler transform).** Durbin's data structure: at every site, the
  haplotypes are kept sorted by the sequence of alleles they carry *to the left* of that site, so
  that haplotypes with the longest shared suffix ending there sit next to each other. Updating the
  order from one site to the next costs one pass over the haplotypes. A *divergence array* records,
  for each neighbour pair in that order, where their match begins.
- **IBS / IBD / IBD2.** Identical by state: same alleles. Identical by descent: same alleles
  because inherited from a common ancestor. IBD2: two samples share *both* haplotypes over a
  region (siblings, duplicates, inbred lines). All tools here refuse to copy from an IBD2 partner,
  because that partner's phase is not independent evidence.
- **Genotype graph / segments (SHAPEIT).** SHAPEIT cuts each sample's chromosome into segments
  holding at most three heterozygous sites; within a segment the 2³ = 8 possible haplotypes are
  all kept, and the phasing problem becomes choosing one *pair* of those 8 per segment and linking
  the choices across segment boundaries.
- **Composite haplotype (Beagle).** A mosaic built *before* running the HMM: a single HMM state
  that copies from haplotype A for a while, then from haplotype B, and so on, so that a fixed
  number of states can still cover many different haplotypes along the chromosome.
- **Scaffold.** A set of already-phased heterozygous genotypes that a tool treats as fixed
  constraints (from trios, from a previous run, or from a reference panel).
- **Phase set (PS).** The VCF FORMAT field that groups heterozygous calls whose relative phase was
  determined from reads; alleles within one PS block are written `0|1` / `1|0` consistently.
- **Dosage (DS), genotype probability (GP), haploid dosage (HDS).** Imputation outputs: the
  expected number of alternate alleles (0–2), the three genotype probabilities, and the per-haplotype
  probability of the alternate allele.
- **DR2 / R2 / INFO.** Per-site imputation quality estimates: the estimated squared correlation
  between the imputed dose and the true (unknown) dose, computed from the spread of the doses
  themselves (§5.6, §6.6).

---

## 2. The model they all share

Every tool in this report answers the same question for one target haplotype at a time: *which of
these K other haplotypes is it copying, here?* The answer is a probability distribution over the
K haplotypes at every site. Two rules drive it. Along the chromosome the copied-from haplotype is
kept with probability `1 − t` and replaced by a random one of the K with probability `t`, where `t`
grows with genetic distance. At each site the copy either matches the observed allele (probability
close to 1) or does not (a small mismatch probability). The forward pass turns those two rules
into `K` numbers per site; the cost is `O(sites × K)` per haplotype, which is why every tool works
so hard to keep `K` small.

The exact numbers each tool plugs in:

| tool | switch probability between neighbouring sites | mismatch (error) probability | jump target |
|---|---|---|---|
| SHAPEIT5 | `t = 1 − exp(−0.04·Ne·d/Nhap)`, `d` in cM, `Ne = 15,000`, `Nhap = 2·(targets+panel)` (`shapeit5/phase_common/src/objects/hmm_parameters.cpp:46`, `:80` of `phaser/phaser_initialise.cpp`) | `ed/ee = 0.0001/0.9999` (`hmm_parameters.cpp:28-29`, used at `models/haplotype_segment_single.h:133`) | uniform over the K conditioning haplotypes (`haplotype_segment_single.h:146-154`) |
| Beagle 5.5 | `1 − exp(−ρ·d)`, `ρ = 0.04·ne/nHaps`, `ne = 100,000` default, then **re-estimated** during burn-in (`beagle_src/src/phase/PhaseData.java:59`, `vcf/MarkerMap.java:221-231`, `main/Par.java:102`) | `θ/(2(θ+nHaps))`, `θ = 1/(ln nHaps + 0.5)` (`main/Par.java:494-497`), re-estimated upward only (`phase/PhaseLS.java:110-112`) | uniform over the 280 states (`phase/HmmUpdater.java:61-65`) |
| Eagle2 | coalescent IBD-length rule: the chance the copied segment (mean length `2 cM`) ends between two split sites is `1/(1+d₁/a)² − 1/(1+d₂/a)²` (`Eagle/src/DipTreePBWT.cpp:117-125`; `expectIBDcM` default at `EagleParams.cpp:113`) | per-genotype error `0.003` (`EagleParams.cpp:117`, applied at `DipTreePBWT.cpp:354-356`) | any haplotype in the HapHedge prefix tree (§4.3) |
| minimac4 | `1 − exp(−cM/100)` from the map, multiplied through untyped sites, floor `1e-5` (`Minimac4/src/recombination.hpp:28`, `input_prep.cpp:257-267`, `prog_args.hpp:45`) | `1 − err + prandom` vs `prandom = err·AF(observed) + 1e-5`, `err = 0.01` (`hidden_markov_model.cpp:265-285`, `prog_args.hpp:46`, `hidden_markov_model.hpp:51`) | uniform, weighted by how many panel haplotypes each unique haplotype stands for (`hidden_markov_model.cpp:250-254`) |
| PBWT | no HMM; votes weighted by match length (§7) | none | — |
| IMPUTE5 | `ρ = 4·Ne·r/N` per the paper (Rubinacci 2020, Methods) | "relatively insensitive to the mutation parameter" (Rubinacci 2020) | uniform over the selected states |

`Ne` is the effective population size; all tools take it as a constant and only Beagle re-fits it
from the data (it reports the fitted value back as `ne = 25·ρ·nHaps`,
`beagle_src/src/phase/PhaseData.java:131-133`).

---

## 3. SHAPEIT5 (and SHAPEIT4)

SHAPEIT5 is two programs run in sequence plus a stitcher: `phase_common` phases the common
variants (the authors use MAF ≥ 0.1%), `phase_rare` then places each rare heterozygote onto the
resulting haplotypes, and `ligate` joins chromosome chunks (Hofmeister 2023). The `phase_common`
core is SHAPEIT4's algorithm (Delaneau 2019) with a parallel PBWT.

### 3.1 Inputs

`phase_common` reads the target VCF/BCF through htslib's `bcf_gt_allele` — hard `GT` only
(`shapeit5/phase_common/src/io/genotype_reader/genotype_reader_reading.cpp:87-101`). A missing
allele on either side marks the genotype missing (`:92`, `VAR_SET_MIS` at `:97`), and missing
genotypes are then imputed as part of phasing (§3.4). There is no reader for `PL`/`GL` anywhere
in `phase_common` or `phase_rare` (a `grep` for `"PL"`/`"GL"` over both source trees returns
nothing). Multi-allelic records are dropped at the scan (`genotype_reader_scaning.cpp:65`), and
so are target sites absent from the panel when a panel is given (`:68`). Fewer than 50 samples
(targets plus panel) is a hard error (`:46-49`).

A genetic map is optional; without one the distance is set to 1 cM per Mb
(`shapeit5/phase_common/src/containers/variant_map.cpp:166-167`, chosen at
`phaser/phaser_initialise.cpp:75-79`). A reference panel (`--reference`) must be phased and
without missing genotypes (`genotype_reader_reading.cpp:136-137`). A scaffold (`--scaffold`) is
consumed only at sites where the target is heterozygous and the scaffold is a phased het
(`:195-200`); a pedigree file is turned into a scaffold for children (`phaser_initialise.cpp:66-71`).

`phase_rare` reads two files: the phased scaffold (missing not allowed,
`shapeit5/phase_rare/src/io/genotype_reader/genotype_reader_reading.cpp:59`) and the unphased
file, from which it keeps only the rare sites and stores, per site, only the samples that are not
homozygous for the major allele — heterozygous, homozygous minor, or missing
(`:106-116`; the 32-bit record is `objects/sparse_genotype.h:39-45`).

### 3.2 Panel versus self: how the cohort conditions on itself

The conditioning set is a bit matrix of *all* haplotypes, targets first, panel after
(`shapeit5/phase_common/src/containers/haplotype_set.h:38-42`; `n_hap = 2·(targets + panel)` at
`haplotype_set.cpp:43`; panel rows written at `genotype_reader_reading.cpp:138-139`). Each target
sample's current haplotype estimate lives in that same matrix, so the cohort is its own panel by
construction; a supplied panel is simply extra rows that never change.

The loop (`shapeit5/phase_common/src/phaser/phaser_algorithm.cpp:122-150`) runs the iteration
scheme `5b,1p,1b,1p,1b,1p,5m` (`phaser_parameters.cpp:49`): 15 iterations, of which the first
five are burn-in (`b`), three are pruning (`p`) and the last five are main (`m`). Each iteration:
rebuild the PBWT and choose neighbours (`H.select()`, `:133`), run the HMM for every sample and
sample one haplotype pair per sample (`phaseWindow()`, `:135`; per-sample sampling at `:92-100`),
merge the IBD2 bans (`:137`), write the sampled haplotypes of the *targets* back into the bit
matrix (`H.updateHaplotypes(G)`, `:139`; only the hets and missing cells change,
`haplotype_set.cpp:49-62`) and transpose for the next PBWT (`:142`). In pruning iterations the
genotype graph is simplified (§3.4). In main iterations the per-segment transition probabilities
are accumulated (`objects/genotype/genotype_sweep.cpp:125-141`) and at the end a most-probable
path through the accumulated probabilities is taken as the answer
(`genotype_sweep.cpp:87-123`, called from `phaser/phaser_finalise.cpp:35`).

Before the first iteration the haplotypes are given a starting phase by a heuristic PBWT sweep
(`H.solve(&G)`, `phaser_initialise.cpp:106`): at each site each unresolved het is phased by adding
up the alleles of its two PBWT neighbours on each side, with a threshold that is lowered until
everything is decided (`containers/conditioning_set/conditioning_set_solve.cpp:88-147`). This is
Durbin's `phase` heuristic transplanted (compare `pbwt/pbwtImpute.c:325-348`, §7). The switch that
should disable it is declared as `--mcmc-noinit` (`phaser_parameters.cpp:51`) but tested as
`pbwt-disable-init` (`phaser_initialise.cpp:106`), so in this fork the sweep always runs.

A sample never conditions on itself: `ibd2_tracks::noIBD2` returns *false* when both haplotypes
belong to the same individual (`containers/ibd2_tracks.cpp:91-94`), and the neighbour selection
consults it for every candidate (`conditioning_set_selection.cpp:120,125`).

### 3.3 State-space reduction: PBWT neighbours at storage points

The PBWT over all `n_hap` haplotypes is rebuilt every iteration by the standard Durbin update
(`shapeit5/phase_common/src/containers/conditioning_set/conditioning_set_selection.cpp:77-105`,
the same code as `pbwt/pbwtCore.c:485-508`), restricted to sites with minor allele count ≥ 5 and
missing rate ≤ 10% (`--pbwt-mac`, `--pbwt-mdr`, `phaser_parameters.cpp:57-58`;
`conditioning_set_managment.cpp:90`). At one randomly chosen site per `--pbwt-modulo` cM
(`:91`, random pick at `conditioning_set_selection.cpp:154-162`) the tool stores, for every
target haplotype, its `--pbwt-depth` nearest neighbours above and below in the sorted order,
taking whichever side has the longer match first and skipping same-individual and IBD2 pairs
(`conditioning_set_selection.cpp:108-148`).

Defaults are chosen from the total sample count when the user does not set them
(`phaser_initialise.cpp:89-93`): `depth = clamp(round(9 − log₁₀ N), 2, 8)` and
`modulo = clamp((ln N − ln 50 + 1)·0.01, 0.005, 0.15)` cM. My arithmetic on those formulas: 500
samples give depth 6 and a storage point every 0.033 cM; 50,000 samples give depth 4 every
0.079 cM. So the *number of stored neighbours per cM falls as the cohort grows* — the
paper's "K shrinks as N increases" (Delaneau 2019, Discussion) is this rule.

For one sample and one phasing window the conditioning set `K` is the union of the stored
neighbours of both of its haplotypes over all storage points inside the window, de-duplicated
(`objects/compute_job.cpp:56-72`). Two protections follow: if two states belong to one other
individual whose heterozygous genotypes agree with the target's on more than 75% of the window,
that pair is treated as IBD2, removed, and banned for 4 cM either side from then on
(`compute_job.cpp:76-107`, `MAX_OVERLAP_HETS` at `:27`; expansion in `ibd2_tracks.cpp:51-64`);
and if fewer than two states remain, 100 haplotypes are drawn round-robin from a shuffled list
(`compute_job.cpp:110-121`, `N_RANDOM_HAPS` at `:28`). The size actually used is printed per
iteration as `K=mean±sd` (`phaser_algorithm.cpp:119`); it is not a fixed number.

**Cost.** The PBWT pass is one sweep over all haplotypes per site per iteration: linear in
`n_hap × sites`. The HMM per sample per window is `sites_in_window × K × 8` float operations
(the 8 is the per-segment haplotype count, §3.4). Memory has two big terms: the two bit-matrix
copies, `2 × n_hap × sites` bits (`haplotype_set.cpp:45-46`), and the neighbour index,
`(depth+1) × storage_points × 2·N_targets` 32-bit integers
(`conditioning_set_managment.cpp:121`). For 50,000 samples and 200,000 common variants that is
5 GB for the bit matrices; with a 20 cM chunk (253 storage points at 0.079 cM) and depth 4 the
index is another 0.5 GB. This is why the authors chunk WGS runs at 20 cM
(`shapeit5/resources/chunks/b38/README`, `docs/docs/tutorials/simulated.md:41`) and why the paper
reports 30.6–52.3 GB for 400,000 array samples on chromosome 20 (Delaneau 2019, memory table).

### 3.4 The HMM: eight haplotypes at a time, sampled, over a genotype graph

**Explained.** SHAPEIT does not run a diploid HMM over `K²` pairs. Instead each sample's
chromosome is cut into segments containing at most three heterozygous sites
(`shapeit5/phase_common/src/objects/genotype/genotype_build.cpp:44-68`; a segment closes when a
fourth het would make `predicted_unfold == 4`, `:52`, or when it holds 22 ambiguous sites,
`MAX_AMB`, `common/src/utils/otools.h:80`). Three hets give 2³ = 8 possible haplotypes
(`HAP_NUMBER 8`, `otools.h:79`); the sample's true pair of haplotypes must be two of those eight
that are complementary at every het. A 64-bit mask per segment (`Diplotypes`,
`genotype_build.cpp:133-151`; pair encoded as `DIP_HAP0/DIP_HAP1`,
`objects/genotype/genotype_header.h:31-32`) says which of the 8×8 pairs are still allowed —
that is the *genotype graph*. Scaffold hets and pedigree-phased hets remove pairs from the mask
before any HMM runs (`genotype_build.cpp:107-113,136`).

Within a segment the HMM is haploid but vectorised: one AVX register carries the 8 haplotypes'
probabilities for one conditioning state, so the forward update for all 8 is one fused
multiply-add per state (`models/haplotype_segment_single.h:151-158` for a homozygous site,
`:221-228` for a het where haplotypes with allele 0 and allele 1 get different emissions
`g0/g1`, `:308-319` for a missing site where all 8 get the same). Underflow in single precision
triggers a rerun in double (`phaser/phaser_algorithm.cpp:55-73`).

At a segment boundary the forward and backward quantities are combined into a probability for
every allowed *pair transition* between the 8×8 pairs of the two segments: a haplotype-level
8×8 table first (`TRANS_HAP`, `haplotype_segment_single.h:355-373`), then the product of the two
haplotypes' entries for each allowed pair of pairs (`TRANS_DIP_MULT`, `:376-391`; additive
fallback on underflow `:394-409`). Those per-boundary tables are the output of one HMM run
(`SET_OTHER_TRANS`, `haplotype_segment_single.cpp:215-230`; forward loop `:87-140`, backward
`:142-199`). The sample's new haplotype pair is then *sampled* segment by segment through the
tables, forwards or backwards at random (`genotype_sweep.cpp:27-85`). Missing genotypes get a
per-haplotype allele probability from the same pass (`IMPUTE`,
`haplotype_segment_single.h:412-430`) and are sampled with it.

Pruning iterations merge neighbouring segments when almost all (`--mcmc-prune 0.999`,
`phaser_parameters.cpp:50`) of the transition mass sits on one pair, so later iterations see fewer
and longer segments (`objects/genotype/genotype_prune.cpp:46-` onward; the merged count is
reported as "Trimming" at `phaser_algorithm.cpp:144-147`).

One detail that matters for rare variants: at a site with MAF < 0.1% (`RARE_VARIANT_FREQ`,
`otools.h:78`; flagged at `hmm_parameters.cpp:50-54`) where the sample is homozygous for the
*common* allele, the forward update is skipped and the previous site is kept as the reference
point for the next transition (`haplotype_segment_single.cpp:110,118`;
`haplotype_segment_single.h:142-144`). Rare sites therefore cost nothing for non-carriers, and
carry information only for carriers — the prelude to `phase_rare`'s reasoning.

The phasing window is chosen per sample: the segment list is split recursively at random points
until a piece would have fewer than 4 segments, fewer than 100 variants or less than
`--hmm-window` (4 cM) (`containers/window_set.cpp:48-51`, default at `phaser_parameters.cpp:63`).

### 3.5 Rare variants: `phase_rare`

**Why the common-variant HMM is not enough.** A rare heterozygote has, by definition, almost no
other carriers, so among a sample's K conditioning haplotypes chosen at *common* sites there is
usually nobody who carries the rare allele at all; the emission at that site then says nothing
about which of the two haplotypes should get it. `phase_common` skips such sites for
non-carriers (§3.4) and, for a carrier, the sampled phase at an uninformative site is close to a
coin flip. The paper states the rare-variant conditioning set must satisfy two properties: the
haplotypes "belong to samples being locally identical-by-descent (IBD) with the target sample"
and "they are polymorphic at the rare variant" (Hofmeister 2023, Methods). `phase_common` cannot
supply the second.

**What `phase_rare` does.** It takes the scaffold (phased common variants) as fixed and, for each
target *haplotype*, builds a conditioning set from two sources during two PBWT sweeps over the
scaffold, forward then backward (`shapeit5/phase_rare/src/containers/conditioning_set/conditioning_set_selection.cpp:50-102`):

1. At storage points every `--pbwt-modulo` cM, the `--pbwt-depth-common` (default 2,
   `phaser/phaser_parameters.cpp:50`) nearest PBWT neighbours, exactly as in `phase_common`
   (`storeCommon`, `:172-219`).
2. At every rare site with more than one carrier, the carriers' haplotypes are looked up in the
   *current* PBWT order (their rank `R`), sorted by that rank, and each carrier haplotype receives
   its `2·--pbwt-depth-rare` nearest *carrier* haplotypes in that order (`storeRare`,
   `:134-170`; `--pbwt-depth-rare` default 2, `:51`). Nearest in PBWT rank means longest shared
   suffix at the common variants, i.e. locally IBD — property (1) — and being a carrier is
   property (2).

Same-individual and IBD2 pairs are excluded (`checkIBD2`, `conditioning_set_header.h:112-119`;
IBD2 pairs found by a genotype-level PBWT requiring ≥ 2.5 cM, ≥ 1 Mb and ≥ 100 scaffold sites of
identical genotypes, `conditioning_set_ibd2.cpp:27-89`, thresholds at `:66`). Sets are padded
with random haplotypes up to 50 (`conditioning_set_selection.cpp:112-120`); the mean size is
logged (`:130`).

The HMM is then **haploid and per haplotype**: forward and backward over the *scaffold* sites
only, with the fixed conditioning set as states, transition `t` from the map and emission
`{1, ed/ee}` (`models/hmm_scaffold/hmm_scaffold_main.cpp:54-109`, `:111-209`; `match_prob` at
`:28`). During the backward pass, at each scaffold interval, every unphased rare genotype of that
sample that falls between scaffold sites `vs` and `vs+1` (the map is built by
`genotype_set_managment.cpp:147-156`) is phased by `phaseLiAndStephens`
(`hmm_scaffold_main.cpp:189-196`). That function walks the conditioning haplotypes and the
carrier list together and adds each state's posterior mass (average of the mass at `vs` and at
`vs+1`, `:194`) to "carries the minor allele" if that state's sample is a non-missing carrier and
to "carries the major allele" otherwise (`containers/genotype_set/genotype_set_phasing.cpp:31-53`).
After both haplotypes of the sample have been processed, the two per-haplotype minor-allele
probabilities `p₀`, `p₁` are combined into the four ordered genotype probabilities with the error
model (for a het only `0|1` and `1|0` survive: `objects/sparse_genotype.h:132-136`) and the
maximum is taken (`:144-151`); if the winning probability is below 0.5001 the het is left
unphased for the next step (`hmm_scaffold_main.cpp:194`; `genotype_set_phasing.cpp:68-71`).

**Singletons** (and anything still unphased) get the coalescent rule: the Viterbi path of each
haplotype (`hmm_scaffold_main.cpp:211-269`) gives, at every position, the length in cM of the
segment currently being copied; the minor allele is put on the haplotype whose copied segment is
*shorter* there (`genotype_set_phasing.cpp:76-135`; the comparison `w0 > w1` at `:124`). The
idea is that a recent mutation sits on a haplotype that shares a *shorter* stretch with everyone
else. The confidence written for such calls is 0.5 unless `--score-singletons` is given (`:133-134`).
The paper reports singleton SER of 35.1% (duos) and 36.6% (trios) on the 147k-sample WGS set,
against 4.36% for alleles seen 11–20 times (Hofmeister 2023, Results).

### 3.6 Outputs

`phase_common` writes phased `GT` only, plus `AC`/`AN`
(`shapeit5/phase_common/src/io/haplotype_writer.cpp:68-84`); no `GP`, no `PS`, no per-site
quality. `phase_rare` writes phased `GT` at every site and, at rare sites, a per-sample `PP`
("Phasing confidence", header at `shapeit5/phase_rare/src/io/haplotype_writer.cpp:79-80`) equal to
the winning ordered-genotype probability rounded to three decimals (`:107-124`; missing for
non-carriers, `:112`). Missing genotypes come out as hard calls in both programs.

### 3.7 Windows, chunks and ligation

Within a run, `phase_common` phases each sample in the per-sample windows of §3.4 and the PBWT
in parallel chunks of `--pbwt-window` (4 cM) with a 0.5 cM look-back buffer
(`conditioning_set_managment.cpp:101-117`). Across runs, the chromosome is cut into overlapping
chunks — the shipped coordinates are 20 cM for `phase_common` and 4 cM for `phase_rare`, each
column-3 region carrying buffers that overlap the neighbours
(`shapeit5/resources/chunks/b38/README`); the UK Biobank WGS run used chunks of about 4.5 Mb with
250 kb buffers for `phase_rare` (Hofmeister 2023, Methods).

`ligate` stitches `phase_common` chunks. For each sample it counts, over the overlap of two
consecutive chunks, how many phased hets agree in orientation between the two files and how many
are flipped (`shapeit5/ligate/src/ligater/ligater_algorithm.cpp:105-124`); if flips outnumber
agreements the second chunk's haplotypes are swapped for that sample from then on
(`:185-189`, swap applied at `:87-103`), except for pedigree-scaffolded samples which are never
swapped (`:187`). The output takes the first half of the overlap from the left chunk and the
second half from the right one (`:207`, `:473-475`). A per-sample "phaseQ" from the
agreement fraction is printed, not written to the VCF (`:193-199`). `phase_rare` chunks need no
ligation because rare variants are phased onto the shared scaffold and buffers are discarded
(`docs/docs/tutorials/simulated.md:112`).

Memory is proportional to `n_hap × sites_in_chunk` (bit matrices) plus the neighbour index; time
to `n_hap × sites` per iteration for the PBWT plus `Σ_samples sites × K × 8` for the HMM.

### 3.8 Read-aware input: SHAPEIT4's `--use-PS`, gone in SHAPEIT5

SHAPEIT4 accepts `--use-PS <error rate>` (`shapeit4/src/phaser/phaser_parameters.cpp:39`). The
reader pulls the `PS` FORMAT field and, for each heterozygous genotype that is written phased and
has a PS code, records the two alleles and the code (`shapeit4/src/io/genotype_reader2.cpp:52-68`;
`pushPS` at `objects/genotype/genotype_header.h:123-126`). Before phasing, `genotype::mask()`
rebuilds each PS block as a partial haplotype over the ambiguous sites of the segment and marks
every pair transition whose haplotype contradicts a block (`ProbabilityMask[t] = false`,
`objects/genotype/genotype_mask.cpp:24-99`). During each iteration the HMM's transition
probabilities are multiplied by `1 − error` where the mask is true and by `error` where it is
false, then renormalised (`objects/compute_job.cpp:220-240`, called at
`phaser/phaser_algorithm.cpp:78`; the mask is recomputed after each pruning, `:133`). So a phase
set is a soft prior on the genotype graph, and the population model can overrule it at the cost
`error`. The paper gives 0.0001 as the default PS error and reports, on the 1000 Genomes
chromosome 20 set (503 samples), 0.82% SER with Illumina-read phase sets alone falling to 0.42%
when a UK Biobank scaffold is added, and 0.07% with 10x Genomics phase sets plus a trio scaffold
(Delaneau 2019, Results on integrating phase sets).

SHAPEIT5 has no such option: the `phase_common` option list (`phaser_parameters.cpp:37-75`)
contains no `use-PS`, a `grep` for `use-PS`/`"PS"` over `phase_common/src` and `phase_rare/src`
finds nothing, and the documentation pages for both programs do not mention phase sets. The
only way to hand SHAPEIT5 read-derived phase is as a *scaffold* file of phased hets
(`--scaffold`), which it treats as hard constraints (§3.1).

---

## 4. Eagle2

Eagle2 (Loh 2016) phases one target sample at a time against a fixed set of haplotypes, using a
PBWT-derived prefix-tree structure (the *HapHedge*) and a *beam search* instead of a full HMM: it
keeps only the most probable few haplotype-pair paths as it moves left to right.

### 4.1 Inputs

VCF/BCF with hard `GT` (`Eagle/src/SyncedVcfData.cpp:60-92`; no `PL`/`GL` reader exists), or
PLINK files in non-reference mode (`EagleParams.cpp:60-65`). A genetic map file is **required**
(`EagleParams.cpp:48`). Missing target genotypes are allowed (`SyncedVcfData.cpp:81-92`) and are
imputed unless `--noImpMissing` (`EagleParams.cpp:86`). Missing panel genotypes are set to the
reference allele and unphased panel genotypes are phased at random, with a warning
(`SyncedVcfData.cpp:345,419-423`).

### 4.2 Panel versus self

With a panel (`--vcfRef`/`--vcfTarget`) the number of iterations is chosen from the target/panel
size ratio: 1 if targets < panel/2, 2 if targets < 2·panel, else 3 (`EagleMain.cpp:85-95`). From
the second iteration on, `useTargetHaps` is true (`:138-141`), so the targets' haplotypes from the
previous iteration join the pool (`Eagle::initRefIter` copies them into the shared bit matrix,
`Eagle.cpp:3527-3537`; `getNlib` returns all of `N` rather than `Nref` after iteration 1,
`:3618`). Without a panel (`--vcf` or PLINK), everything is self-conditioning: an optional
Eagle1-style long-range-phasing step gives an initial phase (`EagleMain.cpp:352-370`; skipped and
replaced by random phase when `--pbwtOnly`, `:371-376`, which the option text says is automatic
for sequence data, `EagleParams.cpp:119`), then `2 + pbwtOnly` PBWT iterations (`:441-445`) in
which each sample is re-phased against everyone's current haplotypes. Samples are processed in
10 batches and each batch's new haplotypes are committed before the next batch starts
(`:457-486`), so later batches condition on already-updated earlier ones within the same
iteration. The sample is always excluded from its own conditioning set (`EaglePBWT.cpp:62,81`).

### 4.3 State-space reduction: K best haplotypes, then a prefix tree

For each target, `findMinErrDipHap` counts, for every haplotype in the pool, the number of sites
where its allele contradicts the target's *homozygous* genotypes (bit-parallel, 64 sites per
word, `EaglePBWT.cpp:55-85`) and keeps the `--Kpbwt` best (default 10,000, `EagleParams.cpp:109`;
halved in non-final iterations, `EagleMain.cpp:139,452`). Those K haplotypes are re-encoded at the
target's *split sites* — its heterozygous sites plus one homozygous site every 0.5 cM
(`EaglePBWT.cpp:158-185`) — and turned into a HapHedge: a sequence of prefix trees, one rooted at
each split site, in which a path from the root is a set of haplotypes sharing that allele sequence
(`HapHedge.hpp:76-129`). The beam search then never enumerates haplotypes: a state is a *node* of
a tree (a group of identical haplotypes over the recent history) and the recombination move
re-roots at a later tree (`DipTreePBWT.hpp:32-61`).

Cost per target is `O(M·K)` to build the hedge and `O(M·H·P)` for the search, `H` the history
length and `P` the beam width (Loh 2016, Methods); the K-selection is `O(N·M/64)` per target, i.e.
quadratic in N over the cohort but with a small constant. At 1,000 haplotypes K is all of them;
at 100,000 it is 10,000 (and 2,500–5,000 in the rougher iterations).

### 4.4 The HMM: a diploid beam search with a confidence-fixing second pass

The search is diploid — each beam entry holds a maternal and a paternal path — but never wider
than the beam: 30 paths with history 30 in the first, unconstrained pass (`EaglePBWT.cpp:274-276`)
and 50 paths with history 100 in the second (`:348-350`). Each extension multiplies in the
coalescent switch probability (§2) and `0.003` per genotype contradiction
(`DipTreePBWT.cpp:299-360`). The relative phase of two consecutive hets is read off the beam by
summing path probabilities over the two orientations (`callProbAA`, `:439-466`). After the first
pass the most confident half (and up to 90%, if the call probability is within 0.01 of 0 or 1)
of the het-pair calls are frozen as constraints (`EaglePBWT.cpp:282,306-319`), the fine pass
re-runs under those constraints, and in the last iteration a reverse-direction pass is combined
with the forward one (`:381-424`). The per-sample "phase confidence" printed at the end is the
mean call confidence (`:369-375`, `EagleMain.cpp:196-200`).

Missing genotypes are imputed from haplotypes sampled out of the beam: ten sampled reference
haplotype pairs per inter-split interval (`:193`, `:359-361`), each voting for its allele at the
missing site, with recombination points between forward and reverse samples placed where the
homozygous genotypes are best explained (`:515-672`). A heterozygote at a site that is monomorphic
among the K haplotypes (a *singleton* relative to the pool) puts its rare allele on the haplotype
whose sampled copying segments are shorter on average (`:675-687`) — the same coalescent rule as
SHAPEIT5's, one paper earlier.

### 4.5 Outputs

Phased `GT` only (`Eagle.cpp:3310-3345`); missing genotypes replaced by hard imputed calls; no
probabilities in the VCF. The paper reports, with the HRC panel (32,470 samples), a SER of 1.36%
against 3.52% for SHAPEIT2 on CEU trios, and genome-wide phasing in about 1.5 minutes per sample
(Loh 2016, Results).

### 4.6 Windows and chunking

Eagle phases a whole chromosome (or `--bpStart/--bpEnd` with `--bpFlanking`) in one process;
windows exist only as the beam's history. Memory is dominated by the haplotype bit matrix
(`2N × M` bits, twice for the transpose) plus the per-target HapHedge over `K × splits`.

### 4.7 Read-aware input: `--usePS`

Hidden option `--usePS 1|2` (`EagleParams.cpp:128`). The reader collects, per sample, pairs of
sites whose `PS` codes match together with their relative orientation (`SyncedVcfData.cpp:359-391`,
counts of usable constraints printed). In `runPBWT` those pairs become relative-phase constraints on
the beam: with `1` they constrain only the first (unconstrained) pass (`EaglePBWT.cpp:253-272`),
with `2` also the fine pass (`:336-342`). A constraint is a hard filter on which diploid
extensions are allowed at that split (`DipTreePBWT.cpp:318-348`), not a probability, so a wrong
read-based phase cannot be overruled where it is applied. The fraction of constraints respected
is printed in trio-check mode (`EaglePBWT.cpp:731-741`).

---

## 5. Beagle 5.5

Beagle (Browning 2018 for imputation; Browning 2021 for two-stage phasing) is the one tool here
whose per-sample cost does not depend on cohort size at all: it always runs an HMM with a fixed
number of states — 280 for phasing, 1,600 for imputation — and spends its effort on building
those states well.

### 5.1 Inputs

`gt=` is required and is the only genotype input; the accepted argument list in
`beagle_src/src/main/Par.java` (lines 128-171) is `gt ref out ped map chrom excludesamples
excludemarkers burnin iterations initial-lr phase-states step-scale rare impute imp-states
imp-segment imp-step imp-nsteps cluster ap gp em ne err window window-markers overlap buffer seed
nthreads truth`. There is no `gl=`/`gtgl=`: genotype-likelihood input, which Beagle 4.1 had, is
gone from the 5.x code. (I could not fetch the 5.5 manual PDF to quote its sentence on this; the
argument list is the evidence.) If one allele of a genotype is missing both are treated as
missing (`vcf/VcfRecGTParser.java:30-31,269-271`); missing genotypes are imputed during phasing
(§5.4). Multi-allelic markers are supported throughout (the writer loops over `nAlleles`,
`imp/ImputedRecBuilder.java:149-178`). A genetic map is optional (`Par.java:214`); without it
distances come from base pairs. A reference panel (`ref=`, VCF or bref3) must be phased.

### 5.2 Panel versus self

Each window is phased against `allHaps`, which is the current estimate of every target haplotype
plus the panel haplotypes if any (`phase/CodedSteps.java:56-61`). The initial phase comes from a
PBWT sweep (`phase/PbwtPhaser.java:70-84`, in overlapping 0.5 cM-buffered pieces per thread,
`:177-195`). Then `burnin` (3) plus `iterations` (12) rounds (`Par.java:83-84`;
`main/Main.java:169-182`); burn-in ends early once fewer than 1% of hets change phase in a round
(`Main.java:171,178-180`). Every round rebuilds the PBWT over `allHaps`, so the cohort conditions
on its own latest haplotypes; the sample itself is excluded because `Ibs2.areIbs2` answers *true*
for a sample against itself (`phase/Ibs2.java:210-212`) and the selection rejects IBS2 partners
(§5.3). With a panel, its haplotypes simply sit in `allHaps` and never change.

### 5.3 State-space reduction: PBWT candidates → composite haplotypes

The chromosome is divided into *steps* of `step-scale` (3) times the median inter-marker genetic
distance (`phase/FixedPhaseData.java:131`). A PBWT is run over the markers step by step
(forward in odd rounds, backward in even, `phase/PhaseLS.java:84-88`), and at each step each
target haplotype looks at the block of `T` candidates around it in the sorted order — as many as
still share a match reaching this step — and picks one at random that is not IBS2 with the sample
in that step (`phase/PbwtPhaseIbs.java:188-235`). `T` is 100 during burn-in, then falls linearly
from 90 to 5 across the phasing rounds (`phase/PbwtIbsData.java:31-34,91-102`): early rounds
explore, late rounds exploit.

`BasicPhaseStates` turns that per-step stream of IBS haplotypes into at most `phase-states`
(280, `Par.java:86`) *composite* haplotypes: each slot is one HMM state; a new IBS haplotype that
is not yet a slot either takes an empty slot or replaces the slot whose haplotype was last seen
IBS longest ago, and the replaced haplotype's segment is copied into the slot up to the midpoint
between its last IBS step and the new one (`phase/BasicPhaseStates.java:189-233`; a slot is only
recycled when full or after `max(200 steps, 1 cM)` without IBS, `:80`). The result is a mosaic
per state whose pieces are each locally IBS with the target. Random haplotypes fill an empty
queue (`:204-206`).

IBS2 exclusion: segments where two samples' genotypes are consistent with sharing both haplotypes
over ≥ 2 cM (`phase/Ibs2.java:42`, built from markers with MAF ≥ 0.1, `phase/Ibs2Markers.java:41`)
are recorded per sample pair and consulted at every selection.

**Cost.** The HMM per sample is `markers × 280` (phasing) or `clusters × 1,600` (imputation)
regardless of N; the PBWT is `O(N·markers)` per round. So a cohort of 100,000 haplotypes costs
100 times the PBWT work of 1,000 but the same HMM work per sample — the paper's "computation time
scales linearly with sample size" (Browning 2021).

### 5.4 The HMM: two haploid chains, hets decided one by one, phase frozen progressively

`PhaseBaum2` keeps three forward vectors of 280 floats: `fwd[0]` ignores phase, `fwd[1]` and
`fwd[2]` follow the two haplotypes (`phase/PhaseBaum2.java:89-90`). The update is the plain Li and
Stephens step `fwd[k] ← em·((1−p)·fwd[k]/sum + p/K)` (`phase/HmmUpdater.java:56-68`), on marker
*clusters* (markers within `cluster` 0.005 cM merged, `Par.java:96`) with the mismatch probability
scaled by cluster length and capped at 0.5 (`PhaseBaum2.java:190-195`). A backward pass first
stores the backward vectors at every unphased het (`:151-179`); the forward pass then, at each
unphased het, compares `p11·p22` against `p12·p21` (the two orientations of that het given both
chains) and *swaps* the two chains from there on if the flipped orientation wins (`:324-349`,
`:247-254`). A het is marked *phased* — removed from future rounds — when the winning ratio
exceeds the round's likelihood-ratio threshold, which is infinite during burn-in, decays
geometrically from `initial-lr` (100,000) to 4, and is 1 in the last round
(`phase/PhaseData.java:64-79`). Missing genotypes are imputed per haplotype by the argmax of
state-mass per allele (`PhaseBaum2.java:256-303`).

Parameters are fitted, not fixed: with `em=true` (default, `Par.java:101`) the switch and mismatch
rates are re-estimated in each burn-in round from expected switch counts over up to 500 samples
(`phase/PhaseLS.java:56-64,104-147`; `phase/HmmParamData.java:151-182`;
`phase/ParamEstimates.java:81-113`).

### 5.5 Rare variants: stage 2

If more than 25% of the window's markers are *low-frequency* — at most `max(3, ⌊rare·N⌋)`
carriers with `rare = 0.002` (`Par.java:88`; `phase/FixedPhaseData.java:54,123-128,229`) — stage 1
runs on the common markers only and stage 2 places the rare ones. Stage 2 runs one haploid HMM
per target haplotype over the stage-1 markers with states chosen by a separate PBWT pass that
also accepts *low-frequency-allele sharing* as evidence (`phase/LowFreqPbwtPhaseIbs.java:93-149`,
`STAGE2_CANDIDATES = 10` at `PbwtIbsData.java:34`) and keeps the state posteriors at every stage-1
marker (`phase/Stage2Baum.java:92-97`). For a rare het between stage-1 markers `A` and `B` the
posterior is interpolated linearly between them (`FixedPhaseData.java:316-334`; weight at
`Stage2Baum.java:150-151`) and each state's mass is credited to whichever of the target's two
alleles that state's *sample* carries at the rare site (`:136-170`): the orientation with the
larger product `p(allele₁ on hap 1)·p(allele₂ on hap 2)` wins, ties broken at random
(`:115-126`). Missing rare genotypes are imputed the same way (`:172-`). No singleton rule: a
singleton's carrier list contains only the sample itself, so both products are equal and the
phase is a coin flip (`:120`).

### 5.6 Outputs

Phasing alone writes phased `GT`. When a panel has markers the target lacks, imputation runs
(`Main.java:191-210`): per target haplotype a 1,600-state forward–backward over marker clusters
(`imp/ImpLSBaum.java:80-142`) whose states are composite haplotypes assembled from panel
haplotypes that are IBS with the target over 0.1 cM steps, up to 7 steps
(`imp/ImpIbs.java:57-71`, `imp/ImpStates.java:107-151`; `imp-step 0.1`, `imp-nsteps 7`,
`imp-segment 6.0`, `Par.java:93-95`), and the state posteriors at ungenotyped markers are
interpolated linearly from the flanking clusters (`imp/StateProbsFactory.java:27`,
`imp/ImputedVcfWriter.java:198-200`). The writer emits `GT` (argmax per haplotype), `DS` (sum of
the two haplotypes' alternate-allele probabilities), optional `AP1/AP2` and `GP` (products of
the two haplotypes' probabilities), and per site `DR2` and `AF`
(`imp/ImputedRecBuilder.java:137-179,279-303`). `DR2` is
`(Σ dose² − (Σ dose)²/n) / (Σ dose − (Σ dose)²/n)` over haplotypes (`:317-329`): the variance of
the imputed haploid dose divided by the variance a perfectly-imputed 0/1 allele with the same
mean would have — an estimate of the squared correlation between imputed and true dose that needs
no truth data.

### 5.7 Windows and stitching

Windows of `window` cM (40) with `overlap` cM (2) (`Par.java:104-106`;
`vcf/TargSlidingWindow.java:174-177`). Beagle does not ligate afterwards: the phased haplotypes of
the overlap are carried into the next window as fixed, already-phased genotypes
(`SplicedGT`, `phase/FixedPhaseData.java:117-118`; `Main.java:234-241`), and each window's output
is cut at the midpoint of its overlaps (`vcf/MarkerIndices.java:73,76`). Memory is per window:
`allHaps` (bit-packed, `N × window_markers`) plus per-thread `280 × markers` mismatch tables.

### 5.8 Read-aware input

None. `Par.java` has no phase-set option and the GT parser records only alleles and whether the
genotype is phased (`vcf/VcfRecGTParser.java:153-154,219-224`); an input that is already fully
phased is passed through without phasing (`Main.java:143-146`).

---

## 6. minimac4

minimac4 imputes *already phased* target haplotypes against a reference panel stored in the
`msav`/m3vcf "unique haplotype block" format (Das 2016). It does not phase and has no
panel-free mode.

### 6.1 Inputs

Target: an indexed VCF/BCF/SAV with `GT` (`Minimac4/src/input_prep.cpp:117-167`; the `GT`
vector is read at `:134` and split per alternate allele at `:150-163`). I found no check that the
genotypes are phased (`grep phased` over `input_prep.cpp` and `main.cpp` returns nothing); the
two alleles are used in the order written, as the two haplotypes
(`hidden_markov_model.cpp:109`, `gt[hap_idx]`). Missing alleles are skipped in the emission
(`:109-111`, `observed >= 0`). Reference: must be an indexed MVCF/msav (`input_prep.cpp:82-83,185,195`).
Genetic map optional (`main.cpp:57`); without it, the switch probabilities stored in the reference
file are used (`input_prep.cpp:257`).

### 6.2 State-space reduction: unique haplotypes per block

The reference is cut into blocks in which identical haplotypes are collapsed to one
representative with a *cardinality* (how many panel haplotypes it stands for) and a map from each
panel haplotype to its representative (`unique_haplotype.cpp:10-47`; `unique_map_` and
`cardinalities_` at `unique_haplotype.hpp:20-21`). Block boundaries are placed where adding a
variant would worsen the block's compression ratio, once the block has at least
`min-block-size` variants (`unique_haplotype.cpp:482-513`; CLI bounds 10 and 65,535,
`prog_args.hpp:128-129`; the typed-sites copy uses 16–512, `main.cpp:55`). The HMM state count at
a site is the number of unique haplotypes in that block — small where the panel has few distinct
haplotypes, up to the full panel where every haplotype differs. The paper measured a sevenfold
cost increase for a twentyfold panel increase (1,000 to 20,000 individuals), i.e. sub-linear
(Das 2016, Results); against minimac2 on HRC (32,390 samples) it reports 31 versus 925 CPU-hours
and 0.55 versus 9.31 GB.

### 6.3 The HMM

Haploid, forward then backward, over the *typed* sites only (`hidden_markov_model.cpp:34-119`,
`:121-226`). The forward vector is initialised proportional to cardinality (`:20-32`); the
transition keeps `(1 − recom)` of each state and adds `recom × cardinality / n_templates`
(`transpose`, `:228-263`); the emission is `1 − err + prandom` for a match and `prandom`
otherwise, with `prandom = err × AF(observed allele) + 1e-5` (`condition`, `:265-285`). At a
block boundary the mass of each representative is redistributed to the next block's
representatives through the per-haplotype "junction proportions", which is what makes the
collapsed computation exact for the full panel (`:71-101`, `:153-179`). A "no-recombination"
copy of the vector is carried alongside so that at output time the mass that arrived *without* a
switch can be told apart from the mass that jumped in (`probs_norecom`).

### 6.4 Imputation at untyped sites

At each typed site the posterior over representatives is `constants[i]·left_nr·right_nr +
(left·right − left_nr·right_nr)/n` (`impute_typed_site`, `:287-377`), giving the typed-site dose
and a leave-one-out dose that removes the site's own emission (`:348-374`). For the untyped sites
between two typed sites, only the representatives with posterior above `prob_threshold` (0.01,
`prog_args.hpp:42`) are expanded to their panel haplotypes and re-collapsed in the full-panel
block (`s3_to_s1_probs`, `s1_to_s2_probs`, `:384-448`); the alternate-allele probability is the
mass on carriers plus, for the mass below threshold, the panel allele frequency
(`impute`, `:450-567`, frequency term at `:551`). Doses are binned to 1/1000 (`bin_scalar_`,
`hidden_markov_model.hpp:55`).

### 6.5 Outputs

`HDS` by default; `GT`, `DS`, `GP` on request (`prog_args.hpp:104`; `dosage_writer.cpp:73-79`).
`GP` is the product of the two haploid doses (`dosage_writer.cpp:721-744`). Per site `R2 =
var(HDS) / (AF·(1−AF))` over haplotypes (`:550-558`) — the same construction as Beagle's `DR2` —
and, at typed sites, `ER2`, the squared correlation between leave-one-out imputed dose and
observed allele (`:560-574`); `--min-r2` filters output sites (`:373-379`).

### 6.6 Chunking

`--chunk` 20 Mb regions with `--overlap` 3 Mb on each side used as HMM context and then dropped
(`prog_args.hpp:101,113`; `main.cpp:34-43,293-306`); a chunk with too few typed sites relative to
reference sites (`--min-ratio 1e-4`) is skipped or fails (`prog_args.hpp:116-117`). Memory is
proportional to `typed_sites × unique_haplotypes_per_block` per thread.

---

## 7. PBWT (Durbin's tool)

The `pbwt` program is the reference implementation of the data structure (Durbin 2014,
DOI 10.1093/bioinformatics/btu014) with experimental phasing and imputation commands. Its
importance for this survey is that every tool above embeds its core update, and that its
imputation uses **no HMM at all** — a useful baseline for what "matching" alone buys.

- **Core.** `pbwtCursorForwardsAD` (`pbwt/pbwtCore.c:485-508`) is Algorithm 2 of the paper: one
  pass over the `M` haplotypes reorders them by the allele at site `k` (stable counting sort into
  `a`) and updates the divergence array `d`. The paper's statements: `O(M)` per site to maintain
  the order, `O(NM)` to report all set-maximal matches, and run-length compression "more than a
  factor of a hundred smaller than gzip" on large panels (Durbin 2014, Results).
- **Inputs.** `-readVcfGT` reads biallelic diploid `GT` (`pbwt/pbwtMain.c:213`).
  `-readVcfPL` exists (`:214`) but its body reads the `PL` vectors and only prints the first ten
  (`pbwt/pbwtHtslib.c:174-224`) — a stub; there is no likelihood-aware algorithm. `-readGeneticMap`
  loads a map (`pbwtMain.c:271`).
- **Panel-free phasing (`-phase n`).** A backward sweep then a forward sweep
  (`pbwt/pbwtImpute.c:374-396`); at each site each het is decided by the alleles of its two PBWT
  neighbours on each side, summed over the main PBWT, `n` sparse PBWTs (built from every n-th
  site) and the reverse sweep's result, with a threshold lowered until all hets resolve, then by
  match-length-weighted scores (`phaseSweep`, `:288-372`; `score0` `:260-267`, `score1`
  `:269-276`). Each site's decision is committed before moving on, so the cohort's haplotypes
  condition each other site by site rather than through iterations.
- **Reference phasing (`-referencePhase`).** A dynamic programme keeping one cell per reference
  haplotype per target with several "extension" scoring variants (`referencePhase4`,
  `:905-` onward; cells allocated `:913-925`).
- **Imputation (`-referenceImpute`).** Set-maximal matches between each target haplotype and the
  panel at the shared sites (`matchSequencesSweep`, `:1146`); at each panel site the target's
  allele probability is the match-length-weighted vote of the matching panel haplotypes, with
  weight `(k − start)·(end − k)` per match (`:1205-1216`), falling back to the panel frequency when
  no match covers the site (`:1217-1222`). A per-site `imputeInfo` is a correlation-like statistic
  between the vote probability and the hard call (`:1241-1249`). `-imputeMissing` runs the same
  routine self-against-self on the complete sites (`:1323-1371`).
- **Outputs.** VCF/BCF with `GT`, `DS`, `ADS` and a haploid-dosage `DR2` INFO field
  (`pbwt/pbwtHtslib.c:276-278`).

---

## 8. IMPUTE5 (paper only)

IMPUTE5 (Rubinacci, Delaneau & Marchini 2020, PLoS Genet, DOI 10.1371/journal.pgen.1009049) is
closed source; the following is from the paper. It imputes pre-phased target haplotypes against a
panel in its `imp5` format (rare variants with MAF < 1/256 stored as carrier index lists, common
ones as bit vectors) with a genetic map. For each target haplotype it selects conditioning states
with the PBWT, either the `L` haplotypes with smallest divergence or `L/2` neighbours on each side,
re-selected every 0.02 cM (default), with `L = 4`; the state set therefore changes along the
chromosome and its size varies by target. The HMM is Li and Stephens over the typed markers only,
forward–backward, with panel-marker probabilities interpolated linearly between typed markers
("delayed lazy imputation" skips monomorphic sites). Because the number of selected states falls
as the panel grows, the reported cost is sub-linear: from 10,000 to 1,000,000 panel haplotypes
(100×) the run time grew 2.5×; on HRC (about 65,000 haplotypes) it was 30× faster than minimac4
and 3× faster than Beagle 5.1 with lower memory, at equal accuracy (Rubinacci 2020, Results). The
accuracy metric in the paper is r² between the true masked allele and the posterior allele
probability, binned by MAF; the output includes an INFO score of the same family as `DR2`.

---

## 9. Comparison

| | SHAPEIT5 | Eagle2 | Beagle 5.5 | minimac4 | PBWT | IMPUTE5 |
|---|---|---|---|---|---|---|
| genotype input | hard `GT`, biallelic, ≥ 50 samples | hard `GT` (VCF/PLINK) | hard `GT`, multi-allelic OK | phased `GT` | biallelic `GT`; `PL` reader is a stub | phased haplotypes |
| genotype likelihoods | no | no | no (dropped after 4.1) | no | no | no |
| genetic map | optional (1 cM/Mb default) | required | optional | optional | optional | required |
| missing genotypes | imputed during phasing | imputed (hard) | imputed during phasing | skipped in emission | `-imputeMissing` | — |
| panel-free | yes (cohort = panel) | yes (2–3 PBWT iterations) | yes | no | yes (`-phase`) | no |
| state selection | PBWT neighbours at storage points, `depth` 2–8, per window union | K = 10,000 by Hamming distance, then HapHedge prefix trees | PBWT IBS candidates (100→5 per step) folded into 280 composite haplotypes | unique haplotypes per block | set-maximal matches | PBWT `L = 4` neighbours every 0.02 cM |
| states per site | variable (logged); ~hundreds | ≤ K, but the beam visits ≤ 50 pairs | 280 (phasing), 1,600 (imputation) | unique haplotypes in block | number of overlapping matches | variable, small |
| HMM | haploid ×8 haplotypes per segment, diploid at segment joins; MCMC sampling, 15 iterations | diploid beam search, 2 passes + reverse | two haploid chains, hets decided sequentially, 3 + 12 rounds, phase frozen by LR schedule | haploid forward–backward | none (votes) | haploid forward–backward |
| parameters | fixed (Ne 15,000, err 1e-4) | fixed (IBD 2 cM, err 0.003) | EM for switch and mismatch rates | fixed (err 0.01, map) | — | fixed |
| rare variants | `phase_rare`: carriers' IBD neighbours + scaffold HMM; singletons by shortest match | singleton-in-pool hets to shorter haplotype | stage 2: carrier-aware states, interpolated posteriors; singletons random | n/a | n/a | rare panel sites as index lists |
| complexity | PBWT `O(N·M)`/iteration; HMM `Σ M_w·K·8`; memory `2·N_hap·M` bits + index | selection `O(N·M/64)` per target; `O(M·K)` + `O(M·H·P)` | PBWT `O(N·M)`/round; HMM `M·280` per sample, independent of N | sub-linear in panel (7× for 20×) | `O(N·M)` | sub-linear (2.5× for 100×) |
| outputs | phased `GT`; `PP` at rare sites | phased `GT` | phased `GT`; `DS`, `GP`, `AP`, `DR2`, `AF` when imputing | `HDS`, `DS`, `GP`, `GT`, `R2`, `ER2` | `GT`, `DS`, `ADS`, `DR2` | haplotypes, dosages, INFO |
| chunk stitching | `ligate` by het agreement in overlap | none needed (whole chromosome) | overlap phased haplotypes spliced into next window | overlap dropped | — | buffer dropped |
| read-based phase input | none (SHAPEIT4: `--use-PS` soft prior) | `--usePS` hard constraint (hidden) | none | — | — | — |

---

## 10. What a caller would have to hand these tools, and what it would get back

What ng emits today (from `tmp/phasing_research/ng_facts.md`): one record per locus with
per-sample `GT` written unphased, `GQ`, optional `GP`, no `PL`/`GL`; loci are multi-allelic
`CandidateAlleles` records mixing SNPs, indels and STRs; per-sample genotype likelihoods and
per-read chain ids exist in memory during calling but are not written.

**Feeding each tool today.**

- *SHAPEIT5.* Needs a biallelic VCF (split multi-allelic records with `bcftools norm -m-`),
  ≥ 50 samples, a region per run, ideally a map. It would take ng's unphased `GT` and treat
  missing calls as things to impute. Lost: everything ng knows beyond the hard call — `GQ`, the
  likelihoods, the read-backed phase between neighbouring hets. Back: phased `GT`, plus `PP` at
  rare sites from `phase_rare`. ng's cohort allele frequency could drive the common/rare split
  (`--filter-maf` on `phase_common`, `phaser_parameters.cpp:69`).
- *SHAPEIT4.* Same, but ng could add a `PS` field with `0|1`/`1|0` at hets spanned by a read
  pair and hand the read-phase in as a soft prior (`--use-PS 1e-4`). This is the one route where
  ng's chain-id evidence would be *used* by a mainstream tool. It requires ng to retain per-sample,
  per-allele chain ids across neighbouring loci, which it does not today (the note in
  `ng_facts.md` estimates one `u64` per read per allele per het locus).
- *Eagle2.* Needs a genetic map file (hard requirement), VCF `GT`; accepts the same `PS` input as
  a hard constraint (`--usePS`). Back: phased `GT` only; missing calls replaced by hard imputed
  calls.
- *Beagle 5.5.* The least demanding: multi-allelic records pass, no map needed, no minimum sample
  count beyond "more than one" (`phase/LowFreqPhaseStates.java:253`). Back: phased `GT`, and if a
  panel with extra sites is given, `DS`/`GP`/`DR2`. With no panel and no extra sites there is no
  dosage output — the imputation of *missing* calls inside phasing returns hard alleles.
- *minimac4 / IMPUTE5.* Only after phasing, and only with an external panel; irrelevant to the
  panel-free case.
- *PBWT.* Research code; biallelic `GT`; the `PL` command does nothing.

**What is lost, quantified.** All six take the hard call as truth with a mismatch probability of
0.0001 (SHAPEIT5), 0.003 per genotype (Eagle2), about 1 in 15,000 at 1,000 haplotypes (Beagle's
Li–Stephens default, `Par.java:494-497`, then re-fitted) or 0.01 (minimac4). For a sample at
3 reads a site, the chance that a true heterozygote shows only one allele is 2·(½)³ = 1 in 4 (my
arithmetic), so a quarter of its hets arrive as false homozygotes: an error rate two to three
orders of magnitude above what these models assume. Beagle would at least re-fit its mismatch rate
upward (`PhaseLS.java:110-112`); the others would mis-phase and mis-impute around every such
site. This is the boundary between this sub-report and sub-report 3: below roughly 10× the input
must be likelihoods, and none of the tools here can take them.

**Two ng-specific frictions.** (1) IBD2/IBS2 exclusion. SHAPEIT5 bans a partner whose hets agree
on > 75% of a window (`compute_job.cpp:85`); Beagle bans pairs sharing ≥ 2 cM of consistent
genotypes (`Ibs2.java:42`); `phase_rare` bans ≥ 2.5 cM / ≥ 1 Mb genotype matches
(`conditioning_set_ibd2.cpp:66`). In a highly inbred collection (the tomato panel's median
inbreeding coefficient is 0.78) or in any segregant population, most pairs are IBD2 over long
stretches, so these rules would strip the conditioning sets down to the random fill-ins — the tools
would run, but on the wrong evidence. (2) Minimum cohort. SHAPEIT5 refuses below 50 samples;
Eagle2 and Beagle run at any size but with one sample there is nothing to copy from — a
single-sample run of any of these is meaningless, which is the case ng must handle by emitting
"absent", not by calling them.

---

## 11. Questions to discuss

1. **Should ng target hard-call phasers at all, or only likelihood-based ones?** The tools here
   are the right choice for the *many-haplotypes, deep-coverage* part of a cohort (20–40×) and
   nothing else. Recommendation: treat them as an *export format* problem (biallelic split, `PS`
   for SHAPEIT4/Eagle2), not as the module; put design effort into the likelihood-taking family
   (sub-report 3). Trade-off: the export costs nothing to add and lets a deep-coverage cohort be
   phased by SHAPEIT5 or Beagle today; the price is two code paths for phasing output.

2. **If ng grows its own phaser, which design to copy?** Beagle's: fixed 280 composite states, two
   haploid chains, progressive freezing, EM for the two parameters, multi-allelic native,
   per-sample cost independent of N, ~1,300 lines for the HMM plus state builder
   (`PhaseBaum2.java` 350, `BasicPhaseStates.java` 350, `PbwtPhaseIbs.java` 271, `HmmUpdater`,
   `PbwtDivUpdater`). SHAPEIT5's genotype-graph design is faster on AVX but its 8-haplotype
   segments and 64-bit diplotype masks are a larger port and are biallelic only. Recommendation:
   Beagle's design, with the emission generalised to take a genotype likelihood instead of a
   0/1 mismatch — the one change none of these tools made, and the one that would let the same
   code serve 3× samples. Trade-off: Beagle's per-round PBWT is `O(N·M)`; at thousands of samples
   and millions of sites it is the dominant cost, and 15 rounds of it.

3. **What to do about IBD2/IBS2 rules in inbred and segregant material.** Every tool here treats
   "shares both haplotypes" as contamination to be excluded; in an F2 or RIL population it is the
   signal. Recommendation: for the few-haplotype case do not adapt these tools; use the founder-
   haplotype HMMs of sub-report 4, where the state is a founder, not a cohort haplotype.
   Trade-off: two phasing models to maintain, one per population case.

4. **Which read-phase interface to standardise on.** Only SHAPEIT4 (`--use-PS`, soft) and Eagle2
   (`--usePS`, hard) consume `PS`; SHAPEIT5 and Beagle ignore it. Recommendation: if ng emits
   read-backed phase, emit standard `PS` + `0|1` so both consumers work, and keep the retention of
   chain ids across neighbouring hets (the memory estimate in `ng_facts.md`) as the enabling
   change. Trade-off: `PS` conveys the phase but not its confidence; SHAPEIT4's single global
   error rate is the only knob.

5. **Rare-variant phasing needs the AF split that ng already computes.** Both SHAPEIT5 and Beagle
   phase rare variants in a second stage that conditions on *carriers* found through the PBWT
   over common variants. ng's cohort expected allele copies from the EM loop is exactly the
   statistic that picks the stage; and ng's per-sample likelihoods would improve the carrier
   test (a low-coverage carrier is precisely the genotype most often mis-called).
   Recommendation: if a phaser is built, build the two-stage shape from the start (common
   scaffold, then rare placement by carrier-restricted states) rather than one HMM over all sites.
   Trade-off: singletons stay unphased or coin-flipped in all of these designs (35% SER in SHAPEIT5);
   only read evidence phases a singleton, which again argues for retaining chain ids.

**Not verified.** The Beagle 5.5 manual PDF (the vendor's own sentence on genotype-likelihood
input) could not be fetched; the claim rests on the 5.5 argument list. The Loh 2016 UK Biobank
150k timing came back from the paper summary with an ambiguous unit and is omitted. The SHAPEIT4
paper describes querying the PBWT "every 8 variants" with `P` neighbours; the SHAPEIT5 code
queries every `pbwt-modulo` cM with `depth` neighbours, so the paper-era numbers are not the
code's. Hofmeister 2023 does not state `Ne`; 15,000 is the code default.

---

## References

- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. Accurate rare variant phasing of
  whole-genome and whole-exome sequencing data in the UK Biobank. *Nat Genet* 55, 1243–1249 (2023).
  DOI 10.1038/s41588-023-01415-w (PMC10335929).
- Delaneau O, Zagury JF, Robinson MR, Marchini JL, Dermitzakis ET. Accurate, scalable and
  integrative haplotype estimation. *Nat Commun* 10, 5436 (2019). DOI 10.1038/s41467-019-13225-y
  (PMC6882857).
- Loh PR, Danecek P, Palamara PF, et al. Reference-based phasing using the Haplotype Reference
  Consortium panel. *Nat Genet* 48, 1443–1448 (2016). DOI 10.1038/ng.3679 (PMC5096458).
- Browning BL, Zhou Y, Browning SR. A One-Penny Imputed Genome from Next-Generation Reference
  Panels. *Am J Hum Genet* 103, 338–348 (2018). DOI 10.1016/j.ajhg.2018.07.015 (PMC6128308).
- Browning BL, Tian X, Zhou Y, Browning SR. Fast two-stage phasing of large-scale sequence data.
  *Am J Hum Genet* 108, 1880–1890 (2021). DOI 10.1016/j.ajhg.2021.08.005 (PMC8551421).
- Das S, Forer L, Schönherr S, et al. Next-generation genotype imputation service and methods.
  *Nat Genet* 48, 1284–1287 (2016). DOI 10.1038/ng.3656 (PMC5157836).
- Durbin R. Efficient haplotype matching and storage using the positional Burrows–Wheeler
  transform (PBWT). *Bioinformatics* 30, 1266–1272 (2014). DOI 10.1093/bioinformatics/btu014.
- Rubinacci S, Delaneau O, Marchini J. Genotype imputation using the Positional Burrows Wheeler
  Transform. *PLoS Genet* 16, e1009049 (2020). DOI 10.1371/journal.pgen.1009049.
