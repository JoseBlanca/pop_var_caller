# The hidden-duplication filter on HG002 chromosome 20: which true SNPs it removes, and why

*2026-11-01, branch `dup-filter-misses` from main 8578ce52. §1–§6 are the investigation, made
with the filter unchanged; §7 is the change the owner then asked for, built and measured. HG002 alone, GIAB chromosome-20 benchmark at 5x and 30x, `--defaults`, one thread,
`--hidden-dup-tag` so that flagged records stay in the file. Counts are inside GIAB's confident
regions unless stated; the filter itself is `doc/devel/ng/spec/hidden_paralog_filter.md`.*

## 1. Verdict

**At 30x the filter does its job; at 5x it removes more true calls than false ones.** Over every
variant type in the confident regions, scored by vcfeval with the filter applied and not applied:

| depth | true calls removed | false calls removed | records flagged on the whole chromosome |
|---|---|---|---|
| 30x | 17 | 92 | 18,237 of 122,228 scored |
| 5x | 94 | 26 | 13,449 of 89,612 scored |

Almost everything the filter flags lies outside the confident regions, where nothing here can say
whether it was right.

**The useful work is on duplicated ground and the cost is on clean ground.** Split by GIAB's union
of low-mappability and segmental-duplication regions (GIAB stratifications v3.3), the SNP calls the
filter removes are:

| depth | duplicated ground: true / false | clean ground: true / false |
|---|---|---|
| 30x | 3 / 83 | **13 / 1** |
| 5x | 1 / 21 | **79 / 3** |

**Three causes, and they stack** (§3):

1. **The ordinary tail of one-copy depth.** At 30x the 16 true SNPs removed sit in windows read at
   1.25 to 1.59 times one copy's depth, where a three-copy duplication (1.5 times the depth, one
   read in three ALT) fits them better than a heterozygote. The depth model is right about how
   often clean ground reaches that depth; these are its tail, not a mistake in it.
2. **The local duplication rate is pulled towards a run-wide rate that hard ground inflates.** Each
   window's rate is pulled towards 0.145 as if 10 extra records sat in it. 0.145 is high because 15%
   of chromosome 20's records are duplication-like, almost all outside the confident regions. On
   clean ground the pull sets a floor of 0.02 to 0.03 under every window's rate. Without it, 8 of
   the 16 fall below the cut.
3. **At 5x the depth model's tail is too thin.** In all one-copy windows on clean ground, raw
   depth reaches twice the median 18 times in 10,000; the model's Normal allows 3. Of the 79 true
   SNPs removed on clean ground at 5x, 77 sit at 1.8 copies or more at 5x and at 0.89 to 1.52
   copies (median 1.25) at 30x: one copy, read twice as deep by chance.

**Recommendation, for the owner to decide: aim the pull at a typical window's rate instead of
the run-wide rate, after re-measuring it on the sets the windowed rate was adopted on** (§5).
Recomputed from this run's scores, it lowers the true SNPs removed on clean ground from 13 to 5 or
6 at 30x and from 79 to 49 at 5x, It leaves the false SNPs removed on duplicated ground at 83 at
30x and lowers them from 21 to 16 or 17 at 5x. The 5x depth tail needs a change to the depth term, which wants its own
design; §5.3.

## 2. The 16 true SNPs removed at 30x

Every one is a heterozygote in the truth and is called 0/1 by pvc; GATK and freebayes call all 16
(the 2026-11-01 trace of chromosome-20 errors). Window depth is the scorer's relative copy number, the 500-base
window's mean depth over the coverage model's one-copy depth at its GC. σ₀, the spread of a
one-copy window's relative depth in this sample, is 0.143. The rate is the window's fitted
duplication rate as run, then recomputed with no pull towards the run-wide rate.

| position | ground | window depth | ALT / total reads | likelihood ratio | rate as run | rate, no pull | posterior as run | posterior, no pull |
|---|---|---|---|---|---|---|---|---|
| 3,302,248 | clean | 1.43 | 16 / 47 | 4.56 | 0.030 | 0.005 | 0.75 | 0.35 |
| 3,409,098 | clean | 1.54 | 17 / 39 | 4.38 | 0.030 | 0.004 | 0.71 | 0.23 |
| 5,611,257 | clean | 1.50 | 11 / 41 | 7.88 | 0.024 | 0.010 | 0.98 | 0.97 |
| 17,968,809 | clean | 1.45 | 9 / 35 | 6.42 | 0.020 | 0.009 | 0.93 | 0.84 |
| 25,896,374 | segdup/lowmap | 1.25 | 17 / 51 | 1.04 | 0.469 | 0.482 | 0.71 | 0.73 |
| 26,244,580 | segdup/lowmap | 1.47 | 17 / 38 | 2.56 | 0.173 | 0.177 | 0.73 | 0.74 |
| 33,822,290 | segdup/lowmap | 1.58 | 14 / 36 | 6.73 | 0.054 | 0.036 | 0.98 | 0.97 |
| 37,009,683 | clean | 1.49 | 16 / 36 | 2.99 | 0.143 | 0.142 | 0.77 | 0.77 |
| 37,040,632 | clean | 1.59 | 24 / 48 | 4.69 | 0.048 | 0.017 | 0.85 | 0.66 |
| 37,281,228 | clean | 1.48 | 15 / 38 | 3.94 | 0.049 | 0.010 | 0.73 | 0.35 |
| 47,644,375 | clean | 1.52 | 9 / 29 | 6.38 | 0.103 | 0.090 | 0.99 | 0.98 |
| 55,006,232 | clean | 1.30 | 13 / 51 | 4.73 | 0.024 | 0.003 | 0.73 | 0.27 |
| 62,142,625 | clean | 1.47 | 11 / 31 | 4.46 | 0.036 | 0.015 | 0.76 | 0.57 |
| 62,142,644 | clean | 1.47 | 14 / 37 | 4.12 | 0.036 | 0.015 | 0.70 | 0.49 |
| 62,528,316 | clean | 1.55 | 13 / 33 | 5.81 | 0.019 | 0.006 | 0.87 | 0.67 |
| 63,062,254 | clean | 1.36 | 8 / 33 | 4.39 | 0.045 | 0.038 | 0.79 | 0.76 |

The run's cut is a posterior of 0.69. A record's posterior is σ(likelihood ratio + logit rate), so
the two columns that move are the rate and, through it, the posterior.

## 3. The causes

### 3.1 The depth tail, and a split that leans the same way

**Each site is 2 to 4 σ₀ above one copy, and that is how often clean ground gets there.** The
true SNPs called on clean ground at 30x, by window depth, against the Normal(1, σ₀) the score
assumes for one copy:

| window depth | true SNPs called | the Normal expects | removed by the filter |
|---|---|---|---|
| 1.2 – 1.3 | 3,624 | 4,172 | 0 |
| 1.3 – 1.4 | 841 | 1,016 | 2 |
| 1.4 – 1.5 | 123 | 154 | 7 |
| 1.5 – 1.6 | 19 | 15 | 4 |
| 1.6 – 1.8 | 1 | 1 | 0 |

So at 30x the depth model is not the fault: it describes this tail correctly. What makes these
sites look duplicated is depth near 1.5 copies *together with* an ALT share between 0.24 and 0.45,
which is where the three-copy story (1.5 times the depth, one read in three ALT) sits. Of the 983
true SNPs at 1.3 to 1.6 copies, the filter removes 13: the ones whose split also leans. The
likelihood ratios this gives, 1 to 8, are moderate evidence, and they cross the cut only with the
rate's help (§3.2) and the cut's (§3.3).

### 3.2 The pull towards the run-wide rate

**A window's rate is estimated from the records within 25,000 bases of its bin, and then pulled
towards the run-wide rate as if 10 more records sat there at that rate** (spec §3.3). The
run-wide rate is 0.145 at 30x and 0.160 at 5x: one record in seven on chromosome 20 is
duplication-like, nearly all outside the confident regions. A window around these sites holds 40 to 160
records, so 10 records at 0.145 alone give it a rate of about 1.45 / 60 ≈ 0.024, whatever its own
records say. For 10 of the 12 clean sites its own records say 0.003 to 0.017 (the table in §2).

On a likelihood ratio of 4.5 the difference is the verdict: σ(4.5 + logit 0.024) = 0.70 against
σ(4.5 + logit 0.005) = 0.31. The four sites that keep a high rate without the pull are a different
case: 25,896,374 and 26,244,580 lie in duplicated ground beside the centromere, in windows where
the records themselves say 0.48 and 0.18, and 37,009,683 and 47,644,375 sit in windows at 0.14 and
0.09.

### 3.3 The cut is a mean, so it flags individual records at 0.69

**The rule flags the longest run of records, from the most likely duplication down, whose
average chance of being a real variant is at or below the target of 1 in 100.** At 30x the run
holds 18,237 records, most of them near certainty, so the average stays at 1 in 100 while its
last records each carry a 3-in-10 chance of being real. That is what a false-discovery rate
promises, and it does keep its promise overall. But it is spent where duplications are rare: of
the SNPs removed on clean ground at 30x, 13 of 14 are true.

A per-record cut does not help: requiring every flagged record to reach 0.99 removes no true SNP
at 30x, but also only 41 of the 83 false SNPs on duplicated ground (§4).

### 3.4 At 5x, the depth tail is too thin

**Raw read depth in 500-base windows, over the 99,229 such windows of clean confident ground**
(no GC correction; relative to the median, the Normal's spread taken from the interquartile range):

| depth | 5x: windows per 10,000 | Normal | 30x: windows per 10,000 | Normal |
|---|---|---|---|---|
| ≥ 1.6 | 330 | 205 | 6.9 | 2.8 |
| ≥ 1.8 | 84 | 32 | 0 | 0.02 |
| ≥ 2.0 | 18 | 3.3 | 0 | 0 |
| ≥ 2.2 | 3.0 | 0.2 | 0 | 0 |

At 5x a window holds about 17 reads of 150 bases, so its depth is a count whose upper tail is longer than a
Normal's, and one-copy windows reach two copies 5 to 14 times as often as the score allows. Windows
holding a called SNP reach two copies about 1.3 times as often as windows in general, so most of the
excess is not that heterozygotes are easier to call where reads are deep.

There, the two-copy story with an even split (two copies, one carrying the ALT) explains a
heterozygote's reads as well as one copy does, so depth alone decides. 77 of the 79 true SNPs
removed on clean ground at 5x sit at 1.8 copies or more; the same positions read 0.89 to 1.52
copies at 30x.

## 4. The filter on the whole chromosome, by posterior band

SNP calls with a non-reference genotype in the confident regions, true when a called ALT is the
truth's ALT at that position. *Expected real* is the sum over the band of 1 − posterior, which is
what the posterior predicts.

**30x** (cut 0.69):

| posterior | flagged | duplicated ground: true / false | clean ground: true / false | expected real |
|---|---|---|---|---|
| 0.50 – 0.69 | no | 0 / 13 | 2 / 0 | 5.8 |
| 0.69 – 0.80 | yes | 2 / 9 | 8 / 0 | 4.8 |
| 0.80 – 0.90 | yes | 0 / 3 | 2 / 1 | 0.9 |
| 0.90 – 0.99 | yes | 1 / 30 | 3 / 0 | 1.3 |
| 0.99 – 0.999 | yes | 0 / 17 | 0 / 0 | 0.1 |
| ≥ 0.999 | yes | 0 / 24 | 0 / 0 | 0.0 |

**5x** (cut 0.83):

| posterior | flagged | duplicated ground: true / false | clean ground: true / false | expected real |
|---|---|---|---|---|
| 0.50 – 0.69 | no | 5 / 3 | 92 / 2 | 41.6 |
| 0.69 – 0.83 | no | 1 / 8 | 50 / 8 | 16.7 |
| 0.83 – 0.90 | yes | 1 / 6 | 29 / 1 | 5.0 |
| 0.90 – 0.99 | yes | 0 / 9 | 34 / 2 | 2.6 |
| 0.99 – 0.999 | yes | 0 / 6 | 11 / 0 | 0.0 |
| ≥ 0.999 | yes | 0 / 0 | 5 / 0 | 0.0 |

On clean ground the posterior is not calibrated at either depth. At 30x the flagged bands predict
about 7 real variants among the clean-ground SNPs and duplicated-ground SNPs together and hold 16;
at 5x the clean-ground SNPs at 0.99 and above are 16 true and none false.

## 5. Options

Every arm below was recomputed from this run's own likelihood ratios with the windowed rate and
the cut copied from `src/paralog/windowed_prior.rs`. The recomputation reproduces the run: no SNP's
verdict differs and the largest posterior difference is 1.4 × 10⁻⁵. Counts are SNPs removed in
the confident regions, true / false.

| arm | 30x clean | 30x duplicated | 5x clean | 5x duplicated |
|---|---|---|---|---|
| as run: pull 10 records towards 0.145, mean cut | 13 / 1 | 3 / 83 | 79 / 3 | 1 / 21 |
| no pull | 5 / 1 | 3 / 84 | 54 / 1 | 1 / 22 |
| pull 10 records towards 0.01 | 5 / 0 | 3 / 83 | 49 / 1 | 1 / 16 |
| pull 10 records towards 0.02 | 6 / 0 | 3 / 83 | 49 / 1 | 1 / 17 |
| pull 50 records towards 0.145 | 20 / 1 | 2 / 83 | 114 / 4 | 2 / 19 |
| every flagged record ≥ 0.99, pull as run | 0 / 0 | 0 / 41 | 16 / 0 | 0 / 6 |
| one run-wide rate (window 0) | 42 / 1 | 5 / 70 | 194 / 9 | 4 / 19 |

### 5.1 Leave it

At 30x it removes 92 false calls for 17 true ones, which is the trade it was built for. At 5x it
costs 94 true calls for 26 false ones in the confident regions, and an unknown number of each in
the 13,000 records it flags outside them.

### 5.2 Aim the pull at a typical window — recommended

**The pull exists so that a window holding a handful of records does not set its rate from them
alone, and that protection is worth keeping. What it should pull towards is the rate of an
ordinary window, not the run-wide rate**, which is a mixture of clean windows near 0.005 and a
minority of duplicated ones near 0.5. With the anchor anywhere from 0.005 to 0.02 the result
barely moves: clean-ground true SNPs removed fall from 13 to 5 or 6 at 30x and from 79 to 49 at
5x, and false SNPs removed on duplicated ground stay at 83 at 30x and fall from 21 to 16 or 17 at
5x. In code the anchor could be the median of the windows' own rates.

**What it needs before adoption.** The pull of 10 towards the run-wide rate is what the
2026-10-23 measurements were taken at (spec §3.3): the GIAB trio at 5x, 10x, 30x and 300x, the
simulated genome with hidden duplications, and the 20-sample cohort at 3x. The change has to be
re-measured on those, because this chromosome tests clean ground well and duplicated ground only
inside GIAB's confident regions. The median window rate needs no cohort, so it behaves the same at
one sample as at thousands.

### 5.3 The 5x depth tail — a design question, not a knob

Of the 49 true SNPs that the recommended pull still removes on clean ground at 5x, at least 47 sit
at 1.8 copies or more, by chance. Only a depth term with a longer upper tail at low depth stops that: for
example, one-copy depth modelled as a count (a negative binomial on the window's reads) instead
of a Normal on relative depth, or the Normal's upper spread fitted from the histogram's upper side.
Either changes the statistic the spec says is production's, and either has to be checked against
the simulated genome, where at 3x the filter catches 764 of 1,401 duplication artifacts. I
recommend a short design note before any code.

### 5.4 Not recommended

- **A per-record cut** (every flagged record at 0.99 or above). It removes no true SNP at 30x, but
  it also removes only 41 of the 83 false SNPs on duplicated ground at 30x, and 6 of 21 at 5x.
- **One run-wide rate** (`--hidden-dup-window 0`). Worse at both depths: 42 and 194 true SNPs
  removed on clean ground.

## 6. Reproducing

The scripts are kept, untracked, in the main checkout's
`tmp/kept_from_worktrees/dup_filter_misses_2026-11-01/`; the runs' outputs went with the worktree.
§2 to §5 were measured with an investigation-only dump that was never committed: with
`PVC_PARALOG_DUMP` set, pass two printed each record's four scorer inputs and its ratio, and the
windowed rate each record's rate and posterior. `checks.py`, `tail.py` and `counterfactual.py` read
that dump, so rerunning them needs it added back.

- `run_chr20.sh 30x 5x` — the runs, inside the dev container with
  `DEV_EXTRA_MOUNT=/Users/jose/devel`.
- `bands.py` — SNP calls by posterior band, ground and truth; `vcfeval_onoff.sh` — the
  filter on and off over all variant types.
- `tail.py`, `checks.py` — window depth against the Normal; the 5x sites at 30x; rates
  without the pull.
- `counterfactual.py <depth>` — the §5 arms.
- `windows.py`, `wintail.py` — raw window depth on clean ground (`samtools bedcov`).
- The ground: GIAB stratifications v3.3, `GRCh38_alllowmapandsegdupregions.bed.gz`,
  downloaded from GIAB and not kept.

## 7. Built and measured (2026-11-01, the owner's go-ahead)

**The change.** `windowed_verdicts` (`src/paralog/windowed_prior.rs`) estimates every bin's rate
twice. The first round pulls each bin towards the run-wide rate, as before; the median of those
rates is the typical window's rate; the second round pulls each bin towards that median, and its
rates give the posteriors. The median is written on the `##paralogFilter` line as `pull_to=` and in
the run report. Where windows hold too few records to move far from the run-wide rate, the median
stays near it, and the filter behaves as before. Three new tests cover the clean-window case, the
sparse case and the empty case; the 255 paralog tests and the 25 header tests pass. The spec's §3.3
and §3.5 describe it.

**What it was measured on, and what not.** The 2026-10-23 sets — the chromosome-22 trio at 5x to
300x, the simulated genome with planted duplications, the 20-sample simulated cohort — are not on
this laptop and were not re-run. The simulated genome is the only set that counts duplications
caught, so **how many real duplications the change lets through is unmeasured**. Measured instead,
one thread a run, `--defaults`:

**HG002 chromosome 20**, vcfeval inside the confident regions, filter applied:

| depth | true calls, old → new | false calls, old → new | true SNPs removed on clean ground | false SNPs removed on duplicated ground |
|---|---|---|---|---|
| 5x | 57,975 → 58,008 (+33) | 4,616 → 4,620 (+4) | 79 → 53 | 21 → 18 |
| 10x | 76,381 → 76,396 (+15) | 2,090 → 2,091 (+1) | 55 → 41 | 17 → 16 |
| 30x | 81,832 → 81,839 (+7) | 552 → 554 (+2) | 13 → 6 | 83 → 83 |

The anchor is 0.018 at 30x against a run-wide rate of 0.145; the records flagged fall from 18,237
to 18,213. The 5x and 30x SNP counts match §5's offline recomputation exactly.

**GIAB panels**: HG002, HG003 and HG004, each alone over its 100 benchmark regions, at 5x, 10x, 30x
and 300x. The run-wide rate is zero in 20 of the 24 runs and the filter flags at most 5 records. In
all 12 pairs every record, genotype and verdict is the same; only `PARALOG_POST` moves. This shows
the rule leaves sparse, clean ground alone; it cannot show a gain.

**Tomato**, 57 accessions at about 3x (BAQ-free CRAMs, 80 regions of 100 kb). Run-wide rate 0.107,
anchor 0.016. The change unflags 13 of 22,864 flagged records and flags none new. With no truth, the
13 are compared on the two signals a duplication leaves in an inbred cohort:

| records | how many | heterozygous share of called samples, median | depth over the sample's median, median | at 1.5x depth or more |
|---|---|---|---|---|
| flagged under both | 22,851 | 0.10 | 1.20 | 21% |
| flagged under the old rule only | 13 | 0.07 | 1.09 | 0% |
| flagged under neither | 196,418 | 0.02 | 1.00 | 0% |

The 13 sit between the two groups on both signals, so they are neither clearly duplications nor
clearly real variants.

**Cost.** Two rounds of bin rates double the windowed step's time. What share of a run that is
was not measured: the machine ran at a load average near 30 throughout, so no timing from these
runs means anything.

**Reproducing**: in the kept scripts, `anchor/run_one.sh`, `anchor/jobs.txt` (one line a run; strip
trailing spaces before `xargs -L 1`), `anchor/score_one.sh`, `bands.py` and `anchor/tomato_flips.py`;
the old rule is main 8578ce52.
