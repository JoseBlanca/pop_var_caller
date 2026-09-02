# HG002's real slippage against the constants ng calls with

Every benchmark run of ng today scores repeat tracts under HipSTR's shipped constants —
10 reads in 100 misreport the tract length, half of those report a shorter tract, 5 in
100 of the misreports move by more than one repeat — one triple for a 40-base
homopolymer and a 5-copy tetranucleotide alike. This measures what HG002's own 30x reads
actually do, over the same 20,204 repeat tracts ng calls on that sample.

**The short version.** The slip share is not one number and the shipped one is too high:
over the tract ground as a whole the reads slip **27 times in 1,000** where the caller
assumes 100, and on the 8,091 tracts that make up two fifths of the ground the shipped
number is more than ten times too high. The other two constants are wrong in the other
direction: reads lose repeats about **1.8 times as often as they gain them** (3,926 reads
against 2,145, not the assumed 1:1), and about **1 slip in 6 moves by more than one
repeat** (1,053 of 6,071) where the constant says 1 in 20.

| | shipped | measured | how far off |
|---|---|---|---|
| `share_of_reads_that_slip` | 0.10 flat | 0.0039 at an 8-base homopolymer, 0.088 at a 30-base one, 0.0027-0.068 for period 2 | 25x too high at the short end, roughly right only at the longest tracts; 3.7x too high averaged over the tract ground |
| `shorter_share` | 0.50 | 0.636 (period 1), 0.726 (period 2), 0.582 (period 3-6) | 1.3x too low, one direction |
| `fall_off` | 0.05 | 0.161 (period 1), 0.182 (period 2), 0.453 (period 3-6) | 3-9x too low |

The measured numbers are written as parameters-file rows in `slippage_rows.toml`
alongside this file, one row per (period, repeat count) for every stratum in the tract
ground and 15 counts past its longest tract, so a run cannot fall back to the default
anywhere.


## How it was measured

At a tract where HG002's two haplotypes are the **same length**, every read is a read of
one known repeat count, so the spread of observed counts around it is the stutter
distribution with no model in the way. The reads come from
`tmp/attrib/tier_30x_candidates.tsv` — ng's own per-tract candidate dump at 30x, one row
per distinct sequence a read showed, with its read count. The true length comes from
GIAB's phased tandem-repeat truth VCF applied to GRCh38.

Two things had to be fixed before any number meant anything.

### The truth window was too narrow, by 1,652 tracts

`benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py` reconstructs a haplotype from
truth records inside `[tract start - 1, tract end]`. A repeat-length change is
left-aligned to the base before the **repeat run**, and where ng's catalog tract starts
inside a run that begins earlier, that record sits outside the window and the tract reads
as homozygous reference. The symptom: of 1,887 tracts called homozygous-reference whose
reads put fewer than half their weight on the reference length, **1,652 have a truth
record within 30 bases of the tract**.

So the true length is taken here as the reference tract length plus the **net length
change** the truth records make over a window padded 15 bases either side. That is stable
— padding 15 and padding 30 give 1,993 and 1,995 usable homozygous-alternate tracts — and
it cuts the tracts whose reads contradict the truth from 1,943 of 5,837 to 402 of 4,542.

### The merge only shows us tracts where something slipped — at homozygous-reference ones

ng's cohort merge builds a locus only where non-reference reads reach `max(2, 2% of
compared reads)`; at 30x the floor of 2 binds. Checked directly: **every** one of the
2,703 built homozygous-reference tracts has 2 or more non-reference reads and none has
fewer, so that rule is exactly the inclusion rule.

At a tract HG002 carries **at the reference length**, "non-reference reads" *is* "reads
that slipped", so the tracts in the dump are the slippy ones: 2,703 of 13,911 are built
and 11,208 refused. Averaging over the built ones alone gives 140 slipped reads in 1,000,
which is not slippage, it is the condition for being looked at.

At a tract HG002 carries **away from** the reference, every read is already
non-reference, so the rule is met whatever the reads did: 1,993 of 2,034 such tracts were
built, and the 41 that were not have almost no reads at all. **That class is unconditioned
and is the direct measurement.**

For the homozygous-reference class the condition is divided out with a beta-binomial over
the built tracts' non-reference counts *plus* the refused tracts, whose non-reference
count the merge's own rule pins below two — so the sample is complete and there is no
truncation left to correct. A plain binomial does not work here and the data say so:
fitted to the built tracts and asked how many tracts should have been refused, it predicts
4 in 100 where 77 in 100 were refused. Tracts differ from one another far more than
binomial noise allows, which is exactly the spread the beta-binomial's second parameter
carries. The refused tracts' depth is not observed, so it is taken from the
homozygous-alternate tracts, which are built whatever their reads did (median 20 reads,
quartiles 16 and 26).

### Tracts whose reads contradict the truth are dropped

**The rule: fewer than half a tract's reads at the truth's length.** No stratum measured
here puts a third of its reads off the true length, so a tract with a *majority* off it is
telling us its truth label is wrong, not that slippage there is 80 in 100. It removes 292
of 2,625 homozygous-reference tracts and 110 of 1,917 homozygous-alternate ones — 9 tracts
in 100 — and those few carry **51 and 42 reads in 100 of all the apparent slips**.

Keeping them roughly doubles every level (homozygous-alternate 0.049 -> 0.081) and moves
`fall_off` from 0.161 to 0.277. So this rule is the single biggest analyst choice in the
whole measurement, and it is stated rather than buried: everything below is with the
contradicted tracts dropped.


## The measurement

`level` is the combined estimate — the direct one and the beta-binomial one pooled by
reads. `alt lvl` is the direct estimate on its own and `ref lvl` the corrected one;
where they disagree the reader can see it.

### Period 1 — homopolymers, 12,441 of the 20,204 tracts

| repeats | tracts away from ref | tracts at ref | refused | reads | slipped reads | alt lvl | ref lvl | **level** | shorter | fall_off |
|---|---|---|---|---|---|---|---|---|---|---|
| 4-8 | 33 | 75 | 1,721 | 39,216 | 143 | 0.033 | 0.0033 | **0.0039** | 0.40 | 0.29 |
| 9 | 49 | 83 | 1,138 | 27,650 | 127 | 0.025 | 0.0042 | **0.0052** | 0.69 | 0.22 |
| 10 | 61 | 136 | 809 | 22,535 | 267 | 0.031 | 0.012 | **0.013** | 0.69 | 0.18 |
| 11 | 69 | 190 | 727 | 22,746 | 351 | 0.024 | 0.019 | **0.019** | 0.74 | 0.08 |
| 12 | 79 | 160 | 526 | 17,591 | 359 | 0.030 | 0.023 | **0.024** | 0.77 | 0.11 |
| 13 | 102 | 160 | 431 | 15,714 | 459 | 0.044 | 0.035 | **0.036** | 0.67 | 0.14 |
| 14 | 102 | 165 | 282 | 12,666 | 468 | 0.051 | 0.041 | **0.043** | 0.66 | 0.16 |
| 15 | 84 | 138 | 260 | 10,640 | 379 | 0.045 | 0.043 | **0.044** | 0.74 | 0.12 |
| 16 | 102 | 115 | 196 | 9,082 | 459 | 0.065 | 0.055 | **0.057** | 0.63 | 0.14 |
| 17 | 74 | 87 | 122 | 6,013 | 362 | 0.085 | 0.064 | **0.069** | 0.56 | 0.20 |
| 18 | 76 | 77 | 98 | 5,150 | 333 | 0.083 | 0.068 | **0.073** | 0.60 | 0.14 |
| 19 | 48 | 69 | 93 | 4,177 | 246 | 0.078 | 0.066 | **0.069** | 0.55 | 0.19 |
| 20 | 40 | 50 | 67 | 3,057 | 170 | 0.077 | 0.062 | **0.066** | 0.60 | 0.20 |
| 21 | 50 | 37 | 61 | 2,786 | 202 | 0.104 | 0.073 | **0.082** | 0.51 | 0.19 |
| 22 | 46 | 41 | 61 | 2,831 | 165 | 0.104 | 0.054 | **0.067** | 0.46 | 0.18 |
| 23 | 27 | 43 | 50 | 2,202 | 170 | 0.131 | 0.078 | **0.088** | 0.65 | 0.24 |
| 24 | 28 | 25 | 39 | 1,549 | 114 | 0.164 | 0.069 | **0.090** | 0.49 | 0.16 |
| 25 | 18 | 21 | 23 | 1,060 | 89 | 0.126 | 0.091 | **0.099** | 0.70 | 0.28 |
| 26-27 | 25 | 28 | 32 | 1,506 | 97 | 0.099 | 0.070 | **0.077** | 0.72 | 0.12 |
| 28-38 | 20 | 25 | 28 | 1,235 | 88 | 0.114 | 0.080 | **0.088** | 0.50 | 0.21 |

A homopolymer's slip rate rises **twenty-five-fold from 8 repeats to 25**, from 4 reads
in 1,000 to 99. The shipped 0.10 is the value at the very top of that range, applied
everywhere. The median homopolymer in the tract ground is 12 repeats, where the reads say
24 in 1,000.

### Period 2 — 5,359 tracts, median 9 repeats

| repeats | tracts away from ref | tracts at ref | refused | reads | slipped reads | alt lvl | ref lvl | **level** |
|---|---|---|---|---|---|---|---|---|
| 4-6 | 19 | 91 | 1,364 | 31,816 | 76 | 0.0090 | 0.0026 | **0.0027** |
| 7 | 22 | 56 | 611 | 14,964 | 78 | 0.0103 | 0.0050 | **0.0052** |
| 8 | 19 | 41 | 270 | 7,237 | 43 | 0.0087 | 0.0064 | **0.0066** |
| 9 | 23 | 24 | 181 | 4,991 | 39 | 0.0160 | 0.0074 | **0.0084** |
| 10 | 27 | 27 | 132 | 4,098 | 27 | 0.0226 | 0.0042 | **0.0072** |
| 11 | 28 | 23 | 106 | 3,486 | 44 | 0.0193 | 0.0130 | **0.0142** |
| 12 | 16 | 19 | 93 | 2,758 | 17 | 0.0190 | 0.0057 | **0.0075** |
| 13 | 22 | 19 | 73 | 2,509 | 31 | 0.0201 | 0.0135 | **0.0148** |
| 14 | 26 | 12 | 48 | 1,752 | 34 | 0.0197 | 0.0260 | **0.0242** |
| 15 | 28 | 14 | 49 | 1,867 | 24 | 0.0148 | 0.0167 | **0.0161** |
| 16 | 40 | 8 | 24 | 1,448 | 63 | 0.0720 | 0.0126 | **0.0445** |
| 17 | 30 | 13 | 27 | 1,406 | 45 | 0.0450 | 0.0341 | **0.0384** |
| 18 | 22 | 10 | 20 | 1,008 | 40 | 0.0686 | 0.0298 | **0.0455** |
| 19 | 29 | 10 | 12 | 933 | 49 | 0.0860 | 0.0176 | **0.0543** |
| 20-21 | 42 | 15 | 16 | 1,276 | 85 | 0.0732 | 0.0671 | **0.0703** |
| 22-23 | 33 | 13 | 16 | 1,058 | 74 | 0.1123 | 0.0381 | **0.0762** |
| 24-47 | 32 | 16 | 15 | 999 | 53 | 0.0915 | 0.0285 | **0.0568** |

The emitted rows follow the pooling that makes these non-decreasing along the
repeat count: 0.0027 at 6 repeats, 0.0052 at 7, 0.0066 at 8, 0.0078 at 9, 0.0112 at 11,
0.0148 at 13, 0.020 at 15, 0.042 at 16, 0.046 at 18, 0.054 at 19, 0.068 at 23 and above.
Even at the longest period-2 tracts the reads never reach the shipped 0.10.

### Periods 3 to 6 pooled — 2,404 tracts, 12 in 100 of the ground

| repeats | tracts away from ref | tracts at ref | refused | reads | slipped reads | alt lvl | ref lvl | **level** |
|---|---|---|---|---|---|---|---|---|
| 3-4 | 7 | 7 | 33 | 1,025 | 7 | 0.0055 | 0.0102 | **0.0094** |
| 5 | 15 | 14 | 117 | 3,104 | 3 | 0.0000 | 0.0015 | **0.0014** |
| 6 | 24 | 65 | 575 | 14,179 | 22 | 0.0019 | 0.0021 | **0.0021** |
| 7 | 32 | 35 | 275 | 7,285 | 18 | 0.0030 | 0.0034 | **0.0034** |
| 8 | 28 | 17 | 131 | 3,655 | 16 | 0.0086 | 0.0049 | **0.0055** |
| 9 | 19 | 22 | 103 | 2,980 | 25 | 0.0188 | 0.0102 | **0.0113** |
| 10 | 32 | 7 | 53 | 1,850 | 36 | 0.0319 | 0.0191 | **0.0234** |
| 11 | 15 | 6 | 45 | 1,350 | 8 | 0.0197 | 0.0024 | **0.0063** |
| 12-13 | 19 | 10 | 30 | 1,136 | 20 | 0.0054 | 0.0312 | **0.0229** |
| 14-22 | 25 | 14 | 25 | 1,115 | 46 | 0.0394 | 0.0577 | **0.0514** |

Four periods are pooled because separately they are too thin to say anything:
periods 5 and 6 contribute 3 and 7 slipped reads. The emitted rows, after pooling to
non-decreasing, are 0.0024 at 6 repeats rising to 0.051 at 16 and above.

### The two shares, and the size distribution behind `fall_off`

| period pool | slipped reads | shorter | shorter_share | moved > 1 repeat | fall_off | tracts contributing a >1 slip |
|---|---|---|---|---|---|---|
| 1 | 5,048 | 3,212 | 0.636 | 812 | 0.161 | 494 |
| 2 | 822 | 597 | 0.726 | 150 | 0.182 | 72 |
| 3-6 | 201 | 117 | 0.582 | 91 | 0.453 | 50 |
| all | 6,071 | 3,926 | 0.647 | 1,053 | 0.173 | 616 |

Neither share shows a usable trend along the repeat count once the tract-to-tract noise
is allowed for, so one pair per period pool is what is emitted.

**The step-size distribution is not geometric, and the emission model assumes it is.**
Over homopolymers the slips are 4,236 at one repeat, 427 at two, 132 at three, and 253 at
four or more. A geometric with `fall_off = 0.161` — which is what the model builds from
this number (`stutter_rates.rs`: step *n* weighs `(1 - fall_off) * fall_off^(n-1)`) —
would put 683 reads at two steps and 21 at four or more. The one-step share is matched by
construction; the far tail is twelve times heavier than the model can express. That tail
is not one bad tract: the five biggest contributors carry 56 of the 812 far slips.

**A note on the part-repeat placeholder.** `PART_REPEAT_SHARE_OF_WHOLE = 0.05` in
`stutter_rates.rs` is documented as a placeholder, never estimated. Measured here: 36
part-repeat reads against 822 whole-repeat slips at period 2 (0.044 — the placeholder is
close) and 41 against 201 at period 3-6 (0.20 — four times the placeholder). At period 1
a part-repeat change cannot exist. This is thin evidence, offered as a pointer rather than
a replacement.


## What this cannot see

- **A read that ran out inside the tract names no length**, so it is not in the dump at
  all and nothing here can count it. The size of that hole is not measurable from the
  dump: the median tract carries 20 reads on a run nominally at 30x, and the gap is some
  mix of partial reads, MAPQ and flank filters, and real coverage variation, which these
  columns cannot separate.
- **Homozygous tracts are not a random sample.** A tract is usable here only if HG002's
  two haplotypes are the same length: 4,158 of 20,204 tracts are dropped for being
  heterozygous in length, and heterozygosity is highest exactly where slippage is —
  the long, unstable tracts. So the long strata rest on the *stabler* long tracts.
- **The two homozygous classes disagree at short tracts, by a factor of ten.** At period 1
  with 4 to 8 repeats the direct estimate is 0.033 and the corrected reference-carrying
  one 0.0033. Those 33 direct tracts are short homopolymers where HG002 differs from
  GRCh38 — a rare and unusual thing, enriched for hard regions — against 1,796 tracts on
  the other side. The combined number follows the reference-carrying class because it
  holds 39,000 of the 39,216 reads, and the direct estimate is shown beside it so the
  disagreement is visible rather than averaged away.
- **The two classes bracket the truth for a different reason too.** A tract HG002 carries
  at the reference length has not mutated in the human lineage and a tract it carries away
  from it has; germline instability and polymerase slippage are driven by the same tract
  properties. So the reference-carrying estimate is the lower side and the direct one the
  upper side of what a randomly chosen tract of that stratum would give. They differ by
  about 1.3-fold over most of period 1 and by more at the ends.
- **The refused tracts' depth is assumed, not observed.** Every correction for the merge's
  inclusion rule uses the depth spread of the tracts HG002 carries away from the
  reference. If refused tracts are systematically shallower — which is plausible, since a
  shallow tract needs a higher slip rate to reach two slipped reads — the corrected levels
  here are too low. The hard bound the merge's rule gives has no such assumption: a
  refused tract has at most one non-reference read, which puts period 1 as a whole between
  0.020 and 0.056, and the fitted 0.025 sits inside it.
- **One sample, one library, one aligner, one depth.** Everything here is HG002 at 30x on
  GRCh38 inside GIAB's confident intervals. Slippage is a property of the library
  chemistry and the aligner as much as of the tract, so these numbers are a fact about
  this benchmark corner, not about the caller — and the parameters file is per read group
  and per slippage group precisely because a second library would need its own.
- **Substitutions inside a tract are invisible to this measure and to the model.** A read
  of the right length with the wrong spelling counts as not-slipped here, which is what
  the length-keyed emission model does too.


## Will re-running with these move genotype accuracy?

**Expect a move; do not assume its sign from the parameters alone.** The three errors do
not push the same way.

The level error is the big one and it is one-directional: over the tract ground the
caller assumes 100 slipped reads in 1,000 where the reads give 27, and on the 8,091
tracts that are two fifths of the ground it assumes more than ten times what happens. A
caller that thinks a tenth of its reads are stutter can explain away a real second allele
at almost any tract; correcting it downward should sharpen exactly the calls that
0.886/0.903 genotype accuracy is losing. That is the same lever the simulator showed,
turned the other way: there, a caller assuming 0.10 against reads slipping 0.25 gave
0.932 and the true model gave 0.990.

Against that, `fall_off` is 3x too low and `shorter_share` 1.3x too low, and correcting
those makes the caller *more* willing to attribute a two-repeat-away read to stutter.
That partly offsets the level correction at long tracts, which are where the two-step
slips are.

Two things temper the expectation. The model cannot express the measured step-size
distribution — the four-or-more tail is twelve times heavier than a geometric with the
fitted `fall_off` — so a share of the remaining error is structural, not a parameter
value. And the numbers are estimated on the tracts the same caller's merge kept, which is
not the tract set a run scores.

The recommendation is to run it: append `slippage_rows.toml` to a run's
`.parameters.toml` (delete its `slippage_by_stratum_and_group = []` line first and put
these rows **at the end** of the file) and re-score the two genotype-accuracy numbers.
It costs one run and it is the only way to settle the sign.


## Files

| file | what it is |
|---|---|
| `extract_observations.py` | the first join, truth window `[start-1, end]` — kept because the padded one imports its FASTA, truth and dump readers, and because the diagnosis of the window bug is only reproducible against it |
| `extract_observations2.py` | the join used, truth window padded 15 bases; writes `observations_pad15.tsv` |
| `summarise.py`, `diagnose.py`, `diag2.py`, `sum2.py`, `why_homref_is_wrong.py`, `checks.py`, `thin_check.py` | the diagnoses quoted above, in the order they were run |
| `make_tables.py` | the two period tables above, generated from `fit_report.txt` rather than typed |
| `estimate.py` | the plain truncated-binomial attempt, kept because its refusal-prediction failure (4 in 100 predicted against 77 in 100 observed) is what forced the beta-binomial |
| `estimate2.py` | the beta-binomial fit and the hard bounds from the merge's rule |
| `fit_and_write.py` | the per-stratum numbers and `slippage_rows.toml`; output in `fit_report.txt` |
| `ground_weighted.py`, `validate.py` | the tract-ground comparison against 0.10, and the check that the emitted file parses and covers every ground stratum |
