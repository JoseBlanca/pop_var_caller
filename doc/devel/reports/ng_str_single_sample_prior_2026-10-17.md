# The one-sample STR prior pushes near ties toward a homozygote

*2026-10-17, branch `str-single-sample-prior` (from main 354754b8). Prompted by a colleague's GIAB
benchmark of pvc 0e905bbb (`tmp/giab_findings/colleague_report_pvc0e905bbb.md` in the main checkout,
§4, groups b and c).*

## 1. Verdict

**The colleague's finding is real, and the cause is the stated concentration of the STR prior, not
an allele frequency fitted on the sample's own reads.** When a run has nothing fitted to seed a
repeat tract's genotype prior from, the prior holds one chromosome's worth of belief spread over the
locus's candidate lengths. That total corresponds to a between-locus share F of 0.5 (F is the chance
that two chromosomes drawn from the population carry the same length, `F = 1/(1 + total)`). A prior
that small makes a homozygote more probable than a heterozygote at every tract:

| candidate lengths at the locus | prior P(each homozygote) : P(each heterozygote) |
|---|---|
| 2 | 1.5 : 1 |
| 3 | 2 : 1 |
| K | (K + 1)/2 : 1 |

At one outbred sample that is wrong in a known direction. At the biallelic STR loci pvc calls in the
local GIAB data, the truth is mostly heterozygous (§3.2), while the prior gives a heterozygote 25%.

**Recommendation: state a concentration of 10 in place of 1 (F = 1/11 ≈ 0.09).** On the local GIAB
data it makes 17 STR genotypes right and 1 wrong across 5x, 10x and 30x in `--defaults` runs, moves
nothing at 300x, and leaves a fitted cohort's calls byte-identical. Its cost is at the other end of
the range: a fitted *one-sample* run at 5x loses 4 STR genotypes net, and a `--defaults` run on an
inbred cohort changes 58 of about 1.3 million sample genotypes. §5 sets out the decision.

**Group b of the colleague's report (the likelihood itself prefers 1/1) is a separate defect.**
The change fixes 3 of the 24 such cases at 5x locally; the rest have likelihood margins of 1–4
natural-log units, more than any prior on this scale can move.

## 2. The mechanism

### 2.1 Where the seed comes from at one sample

A repeat tract's genotype prior is a Dirichlet-multinomial: each candidate length `j` carries a
concentration `α_j`, a count of chromosomes' worth of belief that the length is common
(`calling_priors.md` §3.1). The STR seed sets `α_j = total × share_j`
(`src/calling/genotype_prior/seed_ssr.rs`, `fill_ssr_seed`). Which total and which shares depends on
the rung of the tract ladder the tract lands on (`population_diversity.md` §4.4,
`str_fit_to_calling.md` §3.2):

| run | shares | total |
|---|---|---|
| `--defaults` (nothing fitted) | flat, `1/K` over the K candidates | `STATED_FLAT_CONCENTRATION` = 1 |
| fitted on one sample | the stratum's own fitted shares | the stated value, written in the file as `fallback_between_locus_share = 0.5, defaulted` |
| fitted on a cohort | the stratum's own fitted shares | the stratum's own fitted total, or the run's median where it has none |

The one-sample fit produces shares but no total, because at one sample the between-locus share and
the sample's own homozygote excess cannot be told apart. So both one-sample cases rest on the same
stated number.

### 2.2 What the leave-one-out term contributes: nothing, as the colleague suspected it might not

The prior's per-locus update adds the other samples' allele copies to the seed
(`calling_priors.md` §6). At one sample that term is exactly zero, so the prior is the seed. The `AF`
written in the record is the posterior mean, an output, not an input to the prior
(`src/vcf/encode.rs`).

### 2.3 Confirmed on real loci

A temporary dump (the patch in `tmp/experiment_knobs.patch`, not committed) printed the seed, the
prior over genotypes, the likelihoods and the posterior at every locus of the GIAB runs. Every STR
locus in a `--defaults` run is seeded `[0.5, 0.5]` at two candidates and `[0.333, 0.333, 0.333]` at
three, and the log-prior row at two candidates is `[−0.288, −0.693, −0.288]` for 0/0, 0/1, 1/1:
P = 0.375, 0.25, 0.375.

The three local cases of the colleague's group c (all at 5x) are all decided by that 0.405
log-unit gap:

| sample | site | log-likelihood 0/1 | log-likelihood 1/1 | written | GQ | AF |
|---|---|---|---|---|---|---|
| HG002 | chr1:16362505 `C>CA` | −4.11 | −4.29 | 1/1 | 4 | 0.78 |
| HG002 | chr1:157979968 `CA>C` | −4.15 | −4.35 | 1/1 | 3 | 0.77 |
| HG003 | chr11:123862748 `T>TA` | −6.20 | −6.39 | 1/1 | 3 | 0.78 |

These match the colleague's example (chr22:21887695, GQ 3.6, AF 0.78): a likelihood lead of about
0.2 for 0/1, overturned by the prior's 0.405 lead for 1/1.

**The fitted one-sample case is different in direction.** The fitted shares put most of the mass on
the reference length: at HG002 30x, two-candidate tracts are seeded with α for the reference between
0.60 and 0.93 and α for the other length between 0.02 and 0.13. That prior pulls toward 0/0 and
then 0/1, not toward 1/1. So "α₀ is small" does not hold in either case: under `--defaults` there are
no shares, and under a one-sample fit the reference share is the large one.

## 3. Measurements

### 3.1 The data and the score

- **GIAB per-sample set**: HG002, HG003 and HG004, each called alone over its own 100 random
  regions, at 5x, 10x, 30x and 300x (`benchmarks/giab/per_sample/`). These are 100 random regions,
  not whole chr22, so the counts are a few hundred STR alleles a depth where the colleague had
  thousands; only 3 group-c cases fall in them.
- **Score**: every truth allele inside a locus pvc called on the STR path, matched to pvc's record
  on position and alleles after splitting and left-aligning both. An allele scores right when pvc
  wrote the same number of copies as the truth. Counts are pooled over the three samples.
- **Tomato**: 63 accessions over the first 10 regions of `benchmarks/tomato1/regions.bed` (1 Mb),
  called once with `--defaults` and once with parameters fitted on those 63 psps.

### 3.2 What a calibrated prior would say

At the biallelic STR loci pvc called (two candidate lengths), the truth genotype was:

| depth | 0/0 | 0/1 | 1/1 |
|---|---:|---:|---:|
| 5x | 3 | 46 | 44 |
| 10x | 9 | 67 | 49 |
| 30x | 5 | 74 | 49 |

At 30x that is 4% 0/0, 58% 0/1 and 38% 1/1, against a prior of 37.5%, 25% and 37.5%. The prior is
wrong about both ratios that decide a near tie: 0/1 against 1/1 (truth 1.5 : 1, prior 0.67 : 1),
and 0/1 against 0/0 (truth about 15 : 1, prior 0.67 : 1). The second ratio rules out anchoring the
seed on the reference length, which would make 0/0 the likeliest genotype.

### 3.3 GIAB, `--defaults`, one sample at a time

STR truth alleles with the right genotype, of 146 at 5x, 183 at 10x, 190 at 30x and 192 at 300x:

| seed | 5x | 10x | 30x | 300x |
|---|---:|---:|---:|---:|
| total 1, flat (today) | 84 | 138 | 169 | 178 |
| total 2, flat | 87 | 139 | 171 | |
| total 4, flat | 90 | 141 | 171 | |
| **total 10, flat** | **90** | **145** | **172** | **178** |
| total 30, flat | 90 | 147 | 172 | |
| one chromosome per candidate (uniform genotype prior) | 89 | 139 | 171 | |
| three chromosomes per candidate | 90 | 143 | 171 | |
| reference anchored (α = 1 for the reference, 0.5 shared by the rest) | 88 | 137 | 169 | |

None of the settings changes the number of false alleles (a written allele the truth lacks: 6, 10
and 7) or of missed truth alleles by more than four.

**Paired, total 10 against today's total 1**, allele by allele:

| depth | genotypes that changed | became right | became wrong |
|---|---:|---:|---:|
| 5x | 8 | 7 | 1 |
| 10x | 7 | 7 | 0 |
| 30x | 3 | 3 | 0 |
| 300x | 0 | 0 | 0 |

The one that became wrong is a true 1/1 at HG004 chr20:53620520 at 5x, now called 0/1. All three
group-c cases are among the ones fixed. Truth heterozygotes written as homozygotes fall from 27 to
21 at 5x, from 13 to 9 at 10x and from 3 to 1 at 30x.

**No SNP or indel genotype changes.** Outside STR records, the only differences are the paralog
posterior in its fourth significant digit (the paralog rate is fitted over the run's records, STR
records included) and the QUAL of one STR bundle with no reference run (chr1:90667573, which is
seeded from the same flat rung and carries no `STR` flag); its genotype is unchanged.

### 3.4 GIAB, fitted on one sample

The same three samples, each fitted on its own psp and called with that fit; the stated value
edited in the file from F = 0.5 to F = 1/11:

| stated F | 5x | 10x | 30x |
|---|---:|---:|---:|
| 0.5 (today) | 88 | 135 | 169 |
| 1/11 | 84 | 141 | 173 |
| paired: became right / became wrong | 7 / 11 | 11 / 7 | 4 / 0 |

At 5x the larger total strengthens the reference-heavy shares, so true 1/1 calls go to 0/1: truth
homozygotes written as heterozygotes go from 4 to 15. They are uncertain calls: 12 of the 15 have
GQ 9 or below, and the largest is 29.

### 3.5 Tomato, 63 inbred accessions

| run | records | changed by the new value |
|---|---:|---|
| fitted parameters | 17,220 | **none** — the VCF is byte-identical |
| `--defaults` | 20,742 | 58 sample genotypes of about 1.3 million, all at STR loci; 8 STR records gone, 15 new |

A fitted cohort is untouched because its file carries the cohort's own median
(`fallback_between_locus_share = 0.41, fitted_here`), so the stated constant is never read. Under
`--defaults` the 58 changes are 0/1→1/1 19, 0/0→0/1 15, 1/1→0/1 13 and 11 others; heterozygous
sample calls go from 72,408 to 72,420. Tomato has no genotype truth, so these cannot be scored. On an
inbred panel extra heterozygotes are more likely wrong than right, but a `--defaults` run on such a
panel already assumes every sample outbred (its inbreeding coefficient defaults to zero,
`parameters_file/defaults.rs`), and these 12 are small beside what that assumption does.

## 4. The options weighed

The design principle at issue is that a one-sample run must not come back confident and wrong
(`design_principles.md` §0). Today's prior is not confident — the flips it causes sit at GQ 3–4 — but
it is systematically wrong: it asserts that an outbred genome is homozygous at three tracts in four.

| option | one outbred sample | one inbred sample | cohort | verdict |
|---|---|---|---|---|
| **A larger stated total (F lower)** | fixes the measured flips; 17 right, 1 wrong under `--defaults` | more 0/1 at near ties when nobody declared its F; a declared F restores the homozygote excess through the inbreeding term, where it belongs | fitted: unchanged; `--defaults`: 58 of 1.3 M genotypes move | **recommended** |
| Estimate F from the genome's own heterozygosity | would calibrate to the sample | needs the sample's F, which one sample cannot separate from the between-locus share | — | needs the fit to assume the sample's F, which the owner's 2026-08-31 ruling forbids; does nothing for `--defaults`, where no fit runs |
| Anchor the STR seed on the reference, as the generic seed is | worse: 0/0 is the rarest truth class at these loci (5 of 128 at 30x) | — | — | rejected by §3.2 and §3.3 |
| One chromosome per candidate (uniform genotype prior) | 5 right at 5x, 1 at 10x, 2 at 30x | — | — | better than today, weaker than a total of 10 |

**Why a larger total is the right reading of the model, not a tuning.** The sample's own inbreeding
coefficient is already a separate term of the prior (`calling_priors.md` §3.2). A seed total of 1
puts a second homozygote excess of F = 0.5 under it, in every run that has no fitted total, whatever
the sample. A small F in the seed leaves homozygote excess to the term built for it.

**Why 10 and not larger.** Between 10 and 30 the results are flat (§3.3), and a total much larger
than the number of candidates makes the seed close to "every candidate length equally common",
which is what a flat shape already asserts. At two candidates a total of 10 puts the prior at
0.27 / 0.45 / 0.27, a heterozygote-to-homozygote ratio of 1.67 against 1.5 measured at 30x. It is
still a stated number measured on one human trio, and should be marked soft as the current one is.

## 5. The decision needed

*Taken by the owner on 2026-10-17: yes (§6).*

**Change the stated concentration of the tract ladder from 1 to 10 (F from 0.5 to 1/11)?**

- **If yes**: `STATED_FLAT_CONCENTRATION` becomes 10, the defaulted `fallback_between_locus_share`
  becomes 1/11, and the specs that state 0.5 (`str_fit_to_calling.md` §3.5,
  `population_diversity.md` §4.4, `parameters_file/defaults.rs`) are updated, along with the tests
  that pin the constant. The tomato identity baseline's `--defaults` lines move; its fitted lines do
  not.
- **If no**: one-sample STR calls keep leaning homozygous at near ties, as in the colleague's 167
  records at 5x and 56 at 10x.
- **The trade-off**: a fitted single sample at 5x loses 4 STR genotypes net (§3.4). If that matters
  more than the `--defaults` gain, the second rung (own shares, stated total) can keep 1 while the
  bottom rung takes 10, at the cost of two stated numbers where there is one.

Not done here: group b (§1) needs the stutter model under `--defaults`, which is a separate
investigation.

## 6. Applied

**The owner took the recommendation on 2026-10-17**, and the branch now states a concentration of 10
(`STATED_FLAT_CONCENTRATION`), written in the file as `fallback_between_locus_share = 1/11,
defaulted`. Re-run with the built binary, not the experiment switch:

- **GIAB, `--defaults`, one sample at a time**, 5x, 10x, 30x and 300x: every VCF identical, record
  for record, to the experiment run of §3.3.
- **GIAB, fitted on one sample**: each fit, redone with the new binary, writes the new value, and
  its calls are identical to the edited-file runs of §3.4.
- **Tomato, 63 accessions**: the refitted file differs from the old one only in the comment that
  names the built-in value, and its 17,220 calls are identical. The `--defaults` calls are identical
  to the experiment run of §3.5.
- **The three cross-platform checksum tests** (`src/cli/cross_platform_digests.rs`) were
  re-recorded. Each fixture's fit has a median of its own, so the value is not read; the parameters
  file moves by that comment alone, and the VCF by the parameters file's md5 in its header.
  Putting the old md5 back reproduces each old VCF checksum exactly.
- **The tomato identity oracle** (`scripts/promote_ng_oracle.sh`: four tomato2 accessions, 20
  regions, about 5 reads a window), run on the macOS host after merging main 7a1e2417, whose own
  binary reproduced every recorded checksum first. With the fit, the calls do not move: only the
  comment naming the built-in value, and with it the parameters file's md5 in the VCF header. With
  default parameters, 138 of about 30,000 sample genotypes change, all at STR loci, and heterozygous
  calls go from 3,478 to 3,525. Three STR records are new, and one SNP record is gone: it sits inside
  a poly-A tract where the second sample moves from 0/0 to 0/1 for the deletion. The baseline was
  re-recorded with these explanations (`scripts/promote_ng_oracle.baseline`).
  *Corrected 2026-10-18:* an earlier version of this section said main did not reproduce the
  baseline. That run used tomato1 accessions, not the tomato2 ones the baseline is recorded on.

## 7. Reproducing

All under the worktree's `tmp/`: `experiment_knobs.patch` (the dump and the three environment
switches used for the settings above: `PVC_EXP_STATED_CONC`, `PVC_EXP_PER_ALLELE`,
`PVC_EXP_REF_ANCHOR`), `run_one.sh` and `sweep.sh` (the GIAB runs), `fitted_one.sh` (the one-sample
fits), `score_str.py`, `paired.py`, `truth_classes.py` and `cohort_diff.py` (the scoring).

## 8. The second rung needs its own value (2026-10-18)

**The value of 10 adopted in §6 is right for a flat shape and wrong for a stratum's own fitted
shares**, and until this section both read it. A one-sample fit gives every stratum its shares and
no concentration (§2.1), so every repeat tract of such a run is held at the stated value over shares
that put most of the mass on the reference length. At ten chromosomes those shares become a strong
claim that each copy carries the reference length, and true homozygous-alternative tracts are called
heterozygous.

**Measured on the GIAB HG002 tandem-repeat benchmark** (v1.0.1, 50,000 regions, truth from the
assembly; `benchmarks/ssr_hg002/`), one sample at a time, the sample fitted on its own 30x reads and
called at each depth. Correct genotypes, summing truth heterozygous, homozygous-alternative and
homozygous-reference tracts:

| | 5x | 10x | 30x |
|---|---:|---:|---:|
| `--defaults` (flat shape at 10) | 2,822 | 4,599 | 7,101 |
| own fit, shares held at 10 (main until this change) | 2,487 | 4,476 | 7,046 |
| own fit, shares held at 1 (this change) | 2,758 | 4,518 | 7,056 |

At 5x the larger value cost 471 correct homozygous-alternative calls (1,286 against 815) to gain
205 heterozygous ones.

**What the truth says the value should be.** Over the catalog's tracts in those regions, the
assembly gives each tract HG002's two allele lengths. Per stratum (motif period and reference
repeat count) that yields the true share of chromosomes at each length offset, and how often the
two copies match, `P_same`. For a Dirichlet with shares `m` and concentration `A`,
`P_same = (A Σm² + 1)/(A + 1)`, so the concentration that reproduces the genome is
`A = (1 − P_same)/(P_same − Σm²)` — valid for an outbred genome, as HG002 is:

| stratum (period, reference repeats) | tracts | true share at the reference length | fitted share | concentration from the truth |
|---|---:|---:|---:|---:|
| homopolymer, 10 repeats | 1,489 | 0.911 | 0.905 | 0.97 |
| homopolymer, 12 repeats | 1,202 | 0.797 | 0.793 | 1.20 |
| homopolymer, 15 repeats | 727 | 0.700 | 0.692 | 1.04 |
| homopolymer, 20 repeats | 251 | 0.590 | 0.588 | 1.16 |
| dinucleotide, 8 repeats | 620 | 0.890 | 0.893 | 1.16 |
| dinucleotide, 12 repeats | 330 | 0.691 | 0.644 | 1.28 |
| dinucleotide, 18 repeats | 188 | 0.468 | 0.356 | 1.69 |

Over every homopolymer stratum with 30 or more tracts the concentration runs from 0.9 to 1.6, and
over the dinucleotide ones from 1.0 to 2.5. The fitted shares match the truth's to about 0.02 in most
homopolymer strata; long dinucleotide tracts are under-estimated at the reference length (0.36
against 0.47 at 18 repeats). So the shares are right, and the value that holds them should be
about one.

**The change.** The tract ladder now states two values where no stratum has a fitted concentration:
`STATED_FLAT_CONCENTRATION` = 10 for the flat shape, and `STATED_OWN_SHARES_CONCENTRATION` = 1 for a
stratum's own shares. Where any stratum has a fitted concentration, both hold the run's median, as
the one value did, so a fitted cohort is unchanged. The parameters file carries the second value
under `fallback_between_locus_share_with_own_shares` (format 5); a file written before it holds both
rungs at `fallback_between_locus_share`, as it was written to mean.

**What this does not settle.** Even with the shares held at one, the fitted model does not beat
`--defaults` (table above), and it should: its shares are the genome's own. The fitted run loses
truth heterozygotes to the prior, 76 at 5x and 149 at 10x, where the likelihood preferred the
heterozygote and the prior moved the call. With a prior that matches the genome, that points at the
read likelihood being too ready to explain reads that do not fit a homozygote as stutter or outliers
— the same thing that makes the colleague's group b. That is the next investigation.
