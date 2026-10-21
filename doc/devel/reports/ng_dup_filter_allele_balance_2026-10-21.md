# The duplication filter at 300x: what removes real heterozygotes, and what an allele-balance fix can do

*2026-10-21, branch `dup-filter-allele-balance` (from main bbcc9f46). The owner adopted the
recommendation the same day; §9 records what was applied. Prompted by a colleague's GIAB chromosome-22 benchmark of pvc
bbcc9f46 (`tmp/inform_to_pop_var/README.md` in the main checkout, §5, §6 and §8).*

## 1. Verdict

**The read split removes these heterozygotes, not excess depth.** The colleague reported 98 true
heterozygous SNPs removed by the hidden-duplication filter at 300x outside low-mappability and
segmental-duplication regions. Her three examples have likelihood ratios 19.4, 39.1 and 9.8, and
the real scorer reproduces all three only with a window depth of **1.03 to 1.09 copies**: one copy.
At those depths, moving one site's 245 reads from an even split to 59 ALT reads adds 35 to 44 to
the score, and raising the window from 1.0 to 1.3 copies adds 4 to 15 (§2.2).

**These are not the ordinary scatter of heterozygotes, though, which is what the brief assumed.** On clean ground at 300x, 1,096 real heterozygotes in the local GIAB data scatter 1.5
times as widely as a binomial, and **not one has an ALT share below 0.35** (§3.1). Her three
examples sit at 0.24, 0.18 and 0.27. At all three, GIAB filtered most of its own Illumina
callsets (3 of 4, 4 of 5 and 5 of 6), and the truth rests on PacBio reads, which split close to
even. One site lies in a duplication found in HG002's own 300x alignments, and two in tandem
repeats 201 to 10,000 bases long (§3.2). These are
real variants whose short reads include another copy's reference reads.

**The false calls the filter removes at 300x look the same.** In the local 300x data, all eight
false SNPs the filter removes sit at a window depth of 0.70 to 1.28 copies, with ALT shares of 0.14
to 0.25, and a removed true heterozygote sits among three of them within 100 bases (§5.2). In one
sample's reads nothing separates the two kinds: the filter's gain on duplicated ground and its cost
in true heterozygotes come from the same signal.

**Recommendation: give both hypotheses' allele term the measured read-split scatter, as a fixed
beta-binomial dispersion ρ = 0.0015, and nothing larger.** It corrects a real modelling error. A
strict binomial lets the evidence from one site's split grow without limit with depth, so at
normal depth the highest score a split can reach is 0.7 at 30 reads, 43 at 300 and 209 at 1,000.
The measured dispersion caps that at 23 at 300 reads and 53 at 1,000 (§2.3). It moves no verdict at
30x in the GIAB data. On 59 tomato accessions at about 5x, with the cohort's fitted parameters, it
un-flags 37 of 2,914 flagged records and flags none new, all within 0.75 of the cut. At 300x it
keeps the one true heterozygote main removes locally and still removes 7 of main's 8 false SNPs
(§5).

**A larger dispersion is not recommended.** One that made depth the deciding signal, which is
what production's spec §7 says the filter does, stops every 300x removal on the local data. On
the colleague's chromosome 22 that would most likely give back much of the drop in false
positives, from 1,134 to 252. How much the recommended change recovers there, out of her 98 and
the 368 inside duplicated regions, can only be measured on her files. §7 lists them.

**The spec says something the filter does not do, whichever choice is made.** Production's spec
§7 states that "a lone heterozygote at *normal* coverage is kept", and the single-sample
architecture §3 says the score self-gates because "the reads are few". At high depth neither holds,
and the filter's useful work on duplicated ground at 300x depends on that. A proposed amendment to
`doc/devel/ng/spec/hidden_paralog_filter.md` is in that document's §11.

## 2. The mechanism

### 2.1 What the score compares

At each record and in each sample, the filter weighs two stories against each other
(`src/paralog/locus_score.rs`).

- **One copy, a real variant.** The window's depth is about one copy's, and the reads split as
  the genotype says: almost none ALT for a homozygous reference sample, half for a heterozygote,
  almost all for a homozygous ALT one. The split is scored as a **binomial**: each read is an
  independent coin toss at that share.
- **Collapsed copies.** The sample carries `T` copies of a locus the reference holds once
  (`T` = 3, 4, 6 or 8), and `m` of them carry the ALT base. The depth is `T/2` copies' worth and the
  ALT share is `m/T`. Only `m = 1` and `m ≈ T/2` are scored, so every share this story offers is
  ½ or less: 1/3, 1/4, 1/6, 1/8 or ½.

The score is the log of how much better the second story explains the sample's depth and split
than the first. A run then estimates how common duplications are from all its records' scores,
the rate written as π, and removes records above a cut set by the target false-discovery rate.

### 2.2 Confirmed with the real scorer

A small program, `examples/dup_filter_allele_balance_grid.rs`, calls the scorer on one sample with
the inbreeding coefficient at zero, as `--defaults` sets it. The local 300x GIAB samples fitted a
one-copy scatter σ₀ of 0.084 to 0.092 (σ₀ is how much a one-copy window's depth varies, as a
fraction of its mean).

**The three examples.** The colleague's run's σ₀ and window depths are not known here. With σ₀
from 0.080 to 0.092, the only window depths that give her three ratios are these:

| site (HG002, truth 0/1) | reads REF / ALT | ALT share | her ratio | window depth that reproduces it |
|---|---|---|---|---|
| chr22:19,084,734 G>A | 186 / 59 | 0.24 | 19.4 | 1.05–1.09 copies |
| chr22:19,161,147 C>T | 256 / 58 | 0.18 | 39.1 | 1.03–1.08 copies |
| chr22:23,668,977 A>G | 144 / 53 | 0.27 | 9.8 | 1.05–1.09 copies |

At σ₀ = 0.12 the depths would be lower still, 0.82 to 0.90. A depth above 1.1 copies would need σ₀
below 0.08, finer than any local 300x sample fitted.

**Split against depth, one sample of 245 reads, σ₀ = 0.14:**

| window depth (copies) | 122 ALT (even split) | 59 ALT (her first example) |
|---|---|---|
| 1.0 | −15.1 | +23.2 |
| 1.1 | −12.4 | +25.0 |
| 1.2 | −9.5 | +27.0 |
| 1.3 | −6.3 | +29.4 |

Down a column, depth moves the score by 6 to 9. Across a row, the split moves it by 36 to 38. At
σ₀ = 0.10 the depth effect is 12 to 15 and the split effect 41 to 44; at σ₀ = 0.20 they are 4 and
35. (Scores are natural-log likelihood ratios.)

**A depth error is ruled out locally.** The 300x samples fitted one copy at 302, 307 and 310 reads
a window, against window medians of 305, 309 and 334 measured in the depth-scale report
(`ng_coverage_depth_scale_2026-10-17.md` §3.1). The chromosome-22 fit is not known; the line the
run log prints for it is the first thing asked of the colleague (§7).

### 2.3 How the verdict depends on depth

**For a fixed ALT share, the binomial's evidence grows in proportion to the read count.** An
even split of 30 reads and a 0.31 share are two standard deviations apart; at 276 reads they are
six. Under main's score, these are the most REF-heavy ALT shares a one-sample heterozygote can show
before it reaches a cut of 9.8 (σ₀ = 0.092; "–" means no share between 0.05 and 0.5 reaches it):

| reads | 1.0 copies | 1.1 | 1.2 | 1.3 | 1.5 |
|---|---|---|---|---|---|
| 3 to 30 | – | – | – | – | 0.33 to 0.50 |
| 50 | – | – | – | 0.16 | 0.48 |
| 100 | – | 0.15 | 0.22 | 0.29 | 0.45 |
| 200 | 0.25 | 0.28 | 0.32 | 0.35 | 0.43 |
| 300 | 0.31 | 0.32 | 0.35 | 0.37 | 0.42 |
| 400 | 0.33 | 0.35 | 0.36 | 0.38 | 0.42 |

The cut of 9.8 is the highest the colleague's run can have used, since her third example was
removed at that ratio. At a cut of 2.45 (local HG002 at 300x) the 300-read row is 0.34, 0.36,
0.38, 0.41, 0.50. **ALT-heavy shares are never removed below 1.5 copies.**

**The highest score the split alone can reach, at one copy's depth**, by dispersion ρ. ρ is
defined in §3.1. ρ = 0 is main's binomial.

| ρ | 30 reads | 300 reads | 1,000 reads | at 1.3 copies, 300 reads | at 1.5 copies, 30 reads |
|---|---|---|---|---|---|
| 0 (main) | 0.7 | 43.1 | 208.8 | 60.7 | 17.7 |
| 0.0015 (measured) | 0.7 | 22.6 | 52.9 | 36.2 | 17.5 |
| 0.005 | 0.7 | 7.4 | 14.6 | 21.0 | 16.9 |
| 0.01 | 0.7 | 0.8 | 2.5 | 13.6 | 16.3 |
| 0.02 | 0.7 | 0.7 | 0.7 | 8.4 | 15.6 |
| 0.05 | 0.7 | 0.7 | 0.7 | 4.6 | 14.5 |

A beta-binomial weighs a site's split as if it held at most about 1/ρ reads: 670 at ρ = 0.0015,
50 at ρ = 0.02. The last column is detection by depth, a 30-read site at one and a half copies,
and no dispersion tried moves it by more than 3.2.

### 2.4 The asymmetry

**The score reads the ALT share, every alternative allele summed, and not the minor-allele
share.** So a REF-heavy split counts toward collapsed copies and an ALT-heavy one toward a real
variant: 53 ALT of 197 reads at one copy scores +13.5, and 144 ALT of 197 scores −15.1.

The physical reason holds for the low copy numbers that matter here. The sample has two copies of
the reference locus plus the extra ones, so at most `T − 2` of `T` copies can differ from the
reference: ALT shares of 1/3 at three copies and ½ at four. Only at six and eight copies could a
collapse produce an ALT-heavy split (up to 2/3 and 3/4), and those configurations are not among
the two the score keeps.

**Reference bias pushes real heterozygotes the same way.** At 300x, 760 of the 1,096 clean-ground
heterozygotes lean REF and 314 lean ALT, with a mean ALT share of 0.485. A dispersion fitted with the
mean held at ½ absorbs this (§3.1). For crops it is worth checking at high depth: a divergent
introgressed allele maps worse than a reference one, so its reference bias could be stronger. No
high-depth tomato data is available here.

## 3. What real heterozygotes look like

### 3.1 The scatter on clean ground

GIAB truth heterozygous SNPs in each sample's 100 local benchmark regions, outside the GIAB v3.3
union of low-mappability and segmental-duplication regions, where pvc wrote a matching record.
The read counts are pvc's `AD`.

**ρ, the dispersion, is the correlation between two reads at one site.** The ALT count then has
variance `n·v(1−v)·(1 + (n−1)·ρ)`, so ρ = 0 is a binomial, and at 300 reads ρ = 0.0015 widens the
binomial's variance by a factor of 1.45.

| | 30x | 300x |
|---|---|---|
| sites (three samples) | 1,116 | 1,096 |
| median reads | 28 | 296 |
| scatter around ½, as a multiple of a binomial's | 1.07 | 1.50 |
| ALT share below 0.35 / above 0.65 | 93 / 59 | **0 / 0** |
| ρ by maximum likelihood, mean held at ½ | 0.0012 | 0.0015 |
| ρ with the mean free (mean 0.485) | 0: the binomial suffices | 0.0005 |

Per sample at 300x, with the mean held at ½: HG002 0.0017, HG003 0.0012, HG004 0.0014. At 30x the
binomial already explains the scatter, which is why a dispersion cannot be measured there and does
not matter there.

**The current score never comes near removing these.** At 300x the highest ratio among the 1,096
is −2.75, and none is above zero.

### 3.2 The REF-heavy tail is a different population

The full GIAB v4.2.1 truth for HG002 on chromosome 22, in its confident regions and outside the same
union: 23,919 heterozygous SNPs. GIAB's `ADALL` field sums reads over every dataset (Illumina,
PacBio, 10x).

| ALT share over all datasets | sites | half or more of the Illumina callsets filtered by GIAB |
|---|---|---|
| below 0.36 | 32 | 12 (38%) |
| 0.36 to 0.42 | 347 | 16 (5%) |
| 0.42 to 0.58 | 23,305 | 192 (0.8%) |

Thirty-two in HG002 is the size of the colleague's 98 shared over three samples, though her shares
are from pvc's Illumina reads alone and these include PacBio's. A REF-heavy split is 45 times as
likely as an even one to sit where GIAB's own short-read callers were judged unreliable. The
colleague's three examples are among the 32:

| site | ALT share, all datasets | ALT share, datasets whose callsets called it | Illumina callsets filtered | GIAB's region note |
|---|---|---|---|---|
| chr22:19,084,734 | 38 / 98 = 0.39 | 28 / 51 = 0.55 | 3 of 4 | tandem repeat 201–10,000 bp |
| chr22:19,161,147 | 106 / 313 = 0.34 | 26 / 42 = 0.62 | 4 of 5 | tandem repeat 201–10,000 bp |
| chr22:23,668,977 | 76 / 234 = 0.32 | 27 / 59 = 0.46 | 5 of 6 | `AJtrio-HG002…bilkentuniv…dups`, a duplication called from HG002's own 300x alignments |

## 4. The design

**A beta-binomial allele term, under both stories, with a fixed dispersion ρ = 0.0015.** In
words: a real heterozygote's reads do not split like independent coin tosses at exactly ½. Library
amplification, mapping and the local sequence make some sites lean one way, and the measured lean
is small. The beta-binomial says the site's own share varies a little around ½, by the measured
amount, and the reads are tossed at that share. The same applies at ε (a homozygote's error reads)
and at each collapsed configuration's `m/T`.

The decisions:

- **Fixed, not fitted.** The three samples' values agree within a factor of 1.4 (0.0012 to 0.0017).
  The dispersion has no effect at 30 reads or fewer, since it widens the binomial there by less than
  5%, so a fit at low depth would have nothing to learn and nothing to change. Fitting it per sample
  needs one more read of the spill before scoring, because pass two scores each record as it reads
  it, and its value is unproven. **What would change this:** a second high-depth dataset with
  a different library preparation (a PCR-amplified one, say) whose scatter is more than twice
  this. Then fit per sample, by the method of moments, from records with ALT shares between 0.25
  and 0.75 and at least 50 reads: that window is the bulk of the heterozygotes, which duplications
  barely touch. Per sample rather than per library, because the scorer's counts are per sample.
- **Both stories, not only the one-copy story.** Read sampling does not depend on which copy a
  read came from. Measured, the difference is small: at 300 reads and one copy the crossing moves
  from 0.25 (one-copy side only) to 0.26 (both). Both is the fairer form, and with a binomial left on
  the collapsed side a split close to `m/T` keeps a small advantage that grows with depth.
- **Not a larger, deliberately conservative ρ.** A dispersion of 0.01 or more would make depth the
  deciding signal again, and the spec's sentence true. That is the same move production's spec §7
  declined for depth, when a heavy-tailed depth likelihood was rejected: the tail is the signal.
  §5.2 shows this tail is what the filter catches at 300x.
- **One sample to thousands.** The dispersion acts on each sample's term, so with many samples the
  evidence still adds up across carriers. A collapsed copy shared by a cohort is still found by its
  consistent split. At one sample it is the only allele evidence there is, which is where the cap
  matters.
- **3x to 300x.** At 3 reads the variance factor is 1.003; at 5 to 10 reads a sample, over a
  59-sample cohort, the small per-sample changes add up to 1.3% fewer flags, all at the edge of
  the cut (§5.4). At 300 reads the factor is 1.45, and the highest score a split alone can reach
  halves (§2.3). The change is smooth in depth: no switch, no threshold.
- **Cost.** Six log-gamma calls per hypothesis state and eleven states (three genotypes, the
  non-carrier, seven configurations): 66 per sample-record. On 59 tomato accessions the scoring
  pass took 1.76 s against main's 1.61 s (9%), and scoring is 11% of that run. Half the calls are
  constants per state and can be precomputed if it matters. At thousands of samples it stays
  small beside the collapsed story's 40 × 7 log-sum steps per sample-record.

## 5. Measurements

These were measured with a prototype that read the dispersion under each story from two
environment variables; with both unset the score was main's, and the paralog module's 240 tests
passed, including the differential against production's frozen answers. The variables were removed
when the change was adopted (§9).

### 5.1 GIAB, each sample alone, `--defaults`

`generate-psps` then `call-from-psps --defaults --hidden-dup-tag` over each sample's 100 benchmark
regions (572 kb for HG002). Every sample fitted a coverage model at both depths. Counts are SNPs only,
matched to the truth by position and alleles. A *false SNP* is a non-reference SNP record whose
allele is not in the truth. *Inside* means inside the GIAB low-mappability or segmental-duplication
union.

| depth | dispersion | true hets removed (inside / outside) | false SNPs removed (inside / outside) |
|---|---|---|---|
| 30x | main | 0 / 0 of 97 / 1,116 | 0 / 0 of 6 / 2 |
| 30x | 0.0015, 0.01, 0.02 | 0 / 0 | 0 / 0 |
| 300x | main | **1** / 0 of 96 / 1,096 | **7 / 1** of 10 / 1 |
| 300x | 0.0015 | **0** / 0 | **6 / 1** |
| 300x | 0.01 | 0 / 0 | **0 / 0** |
| 300x | 0.02 | 0 / 0 | **0 / 0** |

The fitted duplication rate and cut, per sample (HG002, HG003, HG004):

| | 300x main | 300x ρ 0.0015 | 300x ρ 0.01 and 0.02 |
|---|---|---|---|
| π | 0.0079, 0, 0.0028 | 0.0070, 0, 0.0014 | 0, 0, 0 |
| cut | 2.45, 30.35, −8.55 | 4.15, 32.25, 2.65 | 28.4 to 32.3 |

At 30x every arm fits π = 0 and cuts at 26 to 32. A negative cut means that no record falls between
the run's ordinary records and its two duplication-like ones; it removes nothing more than a cut
just below those two. These panels are clean ground, chosen from GIAB's confident regions, with
about 1,000 SNPs a sample, so their π is not the chromosome's. Across the 2,334 records per depth,
ρ = 0.0015 moves the 30x scores by 0.3 on average and 1.6 at most.

### 5.2 The nine records main removes at 300x

The window depth is solved from main's ratio and the sample's σ₀, as in §2.2.

| sample | site | truth | REF / ALT | ALT share | window depth | main | ρ 0.0015 | ρ 0.01 | ρ 0.02 |
|---|---|---|---|---|---|---|---|---|---|
| HG002 | chr1:88,342,227 | none | 236 / 76 | 0.24 | 1.07 | **28.0** | **16.7** | 0.4 | −3.6 |
| HG002 | chr1:149,023,063 | none | 398 / 72 | 0.15 | 1.05 | **82.3** | **14.6** | 0.7 | 0.7 |
| HG002 | chr1:149,023,086 | none | 309 / 103 | 0.25 | 1.05 | **37.0** | **19.7** | −0.4 | −4.6 |
| HG002 | chr2:800,840 | **0/1** | 179 / 49 | 0.21 | 0.77 | **11.6** | 4.0 | −4.2 | −2.6 |
| HG002 | chr2:800,897 | none | 233 / 54 | 0.19 | 0.77 | **25.0** | **12.3** | 0.5 | 0.1 |
| HG002 | chr2:800,931 | none | 254 / 55 | 0.18 | 0.81 | **32.0** | **16.5** | 0.7 | 0.5 |
| HG002 | chr2:800,934 | none | 258 / 54 | 0.17 | 0.81 | **33.6** | **14.9** | 0.7 | 0.6 |
| HG004 | chr12:42,532,813 | none | 257 / 43 | 0.14 | 1.28 | **54.2** | **17.1** | 4.2 | 3.3 |
| HG004 | chr13:45,646,237 | none | 167 / 35 | 0.17 | 0.70 | **10.0** | 2.5 | 0.7 | 0.6 |

Bold is removed. All but chr12 lie inside the union. **The true heterozygote is recovered at 4.0
against a cut of 4.15**, which is a narrow margin to build on. The false SNPs keep scores of 12 to
20 under the measured dispersion, and fall to about zero under 0.01.

### 5.3 Tomato, 59 accessions, about 5x

The 59 tomato2 accessions over the first 40 of the benchmark's 160 regions (8 Mb). Depth per
sample runs 4.3 to 13.5 reads, median 5.5. `--defaults`, tag mode, the same psps for both arms.

| | main | ρ 0.0015 |
|---|---|---|
| records scored | 40,078 | 40,078 |
| flagged | 2,565 | 2,557 |
| fitted π | 0.0813 | 0.0812 |
| cut | 4.45 | 4.45 |

Twenty-four records change verdict: 16 flagged only by main, 8 only with the dispersion. All 24
have ratios between 4.11 and 4.78 in both arms, within 0.4 of the cut, and mean depths of 5.5 to
10.2 reads a sample. This is records on the edge of the cut moving by a few tenths, not a change in
what the filter finds.

### 5.4 Tomato, 59 accessions, fitted parameters

The same psps, with parameters fitted on them by `estimate-parameters` (15 minutes; among them each
accession's inbreeding coefficient, which for these selfing lines is what makes a heterozygote
unexpected). This is the run a user would make. The two arms ran one after the other, alone on the
machine.

| | main | ρ 0.0015 on both stories | on the one-copy story only | on the collapsed story only |
|---|---|---|---|---|
| flagged, of 33,424 scored | 2,914 | 2,877 | 2,873 | 2,926 |
| flagged only by main / only by this arm | | 37 / 0 | 41 / 0 | 0 / 12 |
| fitted π | 0.0990 | 0.0982 | | |
| cut | 3.95 | 3.95 | | |
| scoring pass | 1.61 s | 1.76 s | | |

The 37 sat between 3.90 and 4.64 under main and fall to 3.47–3.89, under the cut bin of 3.90 to
4.00. Mean depths there are 2.8 to 11.7 reads a sample. Each sample's split moves by little, but
over 59 samples the small changes add up to between 0.1 and 1.0. The two one-sided columns show
why the dispersion goes on both stories: whichever story alone gets the wider allele model gains
at the edge of the cut. With both, the net effect at this depth is 1.3% fewer flags, all at the
edge.

## 6. The decision needed

**Adopt ρ = 0.0015 as a fixed dispersion on both stories' allele terms: I recommend it.** What it
changes and what it does not:

- **Changes:** at 300x and above, a heterozygote's split counts for at most about 670 reads, so the
  evidence from a REF-heavy split stops growing without limit. Locally it keeps one true
  heterozygote main removes and gives up one of eight false-SNP removals. At 1,000x the highest
  score a split can reach at one copy falls from 209 to 53.
- **Does not change:** any GIAB verdict at 30x, where scores move by 0.3 on average. Detection
  by depth. The introgression argument: an introgression is single copy, so its reads split like
  any heterozygote's, which is now modelled as measured.
- **Changes a little at low depth:** on tomato at about 5x it un-flags 37 of 2,914 records, all
  within 0.75 of the cut.
- **Does not fix:** the colleague's 98. Real heterozygotes whose short reads include another copy's
  reads cannot be told from false calls at one sample by depth and split. Some of them will be
  recovered, those whose scores sit just above her cut. How many, and how many false-positive
  removals are given up in exchange, needs her files.

**Not recommended: a dispersion of 0.01 or more**, which would make depth the deciding signal. On
the local 300x data it ends every removal, and every removal there is a false SNP at about one
copy's depth.

**Either way, the spec should stop saying that a heterozygote at normal depth is kept.** The
proposed amendment is in `doc/devel/ng/spec/hidden_paralog_filter.md` §11.

## 7. What we need from the colleague

In order of use. One sample is enough for the first answer; HG002 is the one with examples.

1. **The 300x psps from bbcc9f46, for HG002, HG003 and HG004.** This branch reads them as they are
   (psp format 1, unchanged since 8b54bf5a). With them, both arms run in tag mode on her exact
   inputs, and every count in §5.1 can be made on chromosome 22.
2. **The `call-from-psps` logs of those 300x runs**, or the lines `sample … fitted one copy at …
   reads a window, scatter …` and `duplication rate … cut at …`. They give σ₀ and the one-copy depth
   (and so rule a depth-scale error in or out), and the cut her three examples were removed at.
3. **The list of the 98 removed heterozygotes** (sample, position), and if available the 368
   added false negatives inside low-mappability or segmental-duplication regions, and the 882
   false positives the filter removed there. Without them, a rerun can only be compared to the
   truth, not to her classification.
4. **Only if the psps cannot be shared:** the 300x chromosome-22 CRAM for HG002. A psp can be
   made from it here.
5. The 30x psps for the same sample, to confirm the 30x verdicts do not move on her data either.

## 8. Reproducing

All scratch is under this worktree's `tmp/`.

- The scoring program: `cargo run --release --example dup_filter_allele_balance_grid -- --sigma0
  0.092 --lr-cut 9.8 [--dispersion R] [--site ALT,TOTAL,COPIES …]` (the dispersion defaults to
  the run's, 0.0015; 0 is main's binomial); `tmp/giab/crossings.py SIGMA0 CUT [--dispersion R]`
  builds the §2.3 tables from it.
- GIAB runs: `tmp/giab/run_panel.sh COVERAGE LABEL`, which ran the prototype's arms through `PVC_ENV`;
  `tmp/giab/analyse.py` joins a VCF to the truth and the union; `tmp/giab/het_spread.py` fits ρ.
  The union is GIAB stratifications v3.3, `GRCh38_alllowmapandsegdupregions.bed.gz`.
- Tomato: `tmp/tomato2/generate.sh` (psps), then `call-from-psps` as in §5.3;
  `tmp/tomato2/flips.py` lists the records whose verdict differs between two runs.

## 9. Applied (2026-10-21)

The owner adopted ρ = 0.0015 on both stories.

- **The code.** `DEFAULT_ALLELE_DISPERSION = 0.0015` and the `allele_dispersion` field in
  `src/paralog/model_params.rs`. In `src/paralog/locus_score.rs` every allele factor is an
  `AlleleFactor` built once per pass, so the three log-gamma terms that do not depend on the read
  counts are computed once rather than per sample-record. The environment variables are gone.
- **The tests.** Production's differential and the hand-computed likelihood test now run at
  dispersion 0, which pins the binomial path bit for bit. New tests check the factor against its
  rising-factorial form, check that dispersion 0 is the binomial bit for bit, and check that the
  default caps the split (a 0.2 share of 1,000 reads at one copy scores 154 under the binomial and
  52 with the default) while 3-read sites and a 30-read site at 1.5 copies stay within 0.1.
- **Re-recorded.** The three cross-platform calls checksums, where only `PARALOG_LR` moved (by
  0.03 to 2.4, every record still passing), and the tomato identity oracle's three VCF digests. In
  the oracle only `PARALOG_LR`/`PARALOG_POST` and π's fourth decimal moved. With default
  parameters one record more is kept (7,508 to 7,509); with the fit, one record trades places in
  each direction. Main reproduced every old value on the same host first.
- **The spec.** `doc/devel/ng/spec/hidden_paralog_filter.md` §1.2, §3.2 and §4 carry the change, and
  §11 records the ruling.
- **Not yet.** The colleague's chromosome-22 measurement (§7), which says how many of the 98 this
  recovers and how many false-positive removals it gives up.
