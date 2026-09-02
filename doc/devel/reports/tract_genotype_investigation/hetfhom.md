# The heterozygotes ng calls where the truth says homozygous

**Ground:** GIAB's HG002 tandem-repeat benchmark, 50,000 Tier intervals, one
sample, at 30x and 50x. **Caller:** ng's `.raw.vcf` (ungated), run with
`--defaults` — no parameter was fitted, which turns out to be the finding.
**Scorer:** `benchmarks/lib/tract_qual_experiment.py`'s own
`compare_genotypes` rule, reused rather than rewritten; `classify.py` here
writes one row a tract instead of one row a period class.

Every count in this file is over the **comparable** tracts — the ones both the
truth set and ng call and whose two sequences the scorer can rebuild. That is
3,515 homopolymer tracts and 2,543 tracts of period 2 or more at 30x. ng's
genotype is right at 3,113 of the 3,515 (0.886) and 2,297 of the 2,543 (0.903),
which reproduces the standing figures.

Of the 402 homopolymer errors, **86 are a heterozygote called where the truth
is homozygous**; of the 246 errors at period 2 and above, **51 are**. Together
137 tracts, about 2 in every 100 comparable tracts of either class.

---

## 1. What these tracts look like

| | homopolymer het-for-hom (86) | homopolymer right (3,113) | period 2+ het-for-hom (51) | period 2+ right (2,297) |
|---|---|---|---|---|
| tract length, median repeat units | **19** | 14 (het) / 15 (hom) | **17** | 13 (het) / 12 (hom) |
| DP, median | 30 | 32 / 31 | 31 | 34 / 34 |
| GQ, median | **34.5** | 72 (het) / 57 (hom) | **38** | 94 (het) / 57 (hom) |
| QUAL, median | 165 | 91 (het) / 276 (hom) | 192 | 200 (het) / 333 (hom) |
| reads on the weaker called allele, median | **3** | 9 (het) | **5** | 7 (het) |
| weaker allele's share of the two, median | **0.31** | 0.42 (het) | 0.39 | 0.41 (het) |

Read the table by column pairs: the tracts ng gets wrong this way are **longer**
and their two called alleles are **more lopsided**, and they are not deeper or
shallower than the rest.

**Genotype shape.** Of the 86 homopolymer cases, 38 are `0/1`, 35 are `1/2`
(two different non-reference alleles), and 13 are tracts holding two records.
Of the 51 period-2+ cases, 17 are `0/1` and 34 are `1/2`. So in half of them ng
is not calling "reference plus a slip product" — it is calling two different
non-reference lengths.

**Period.** All 86 homopolymer cases are period 1 by construction. The 51
period-2+ cases are 38 dinucleotide, 6 trinucleotide and 7 tetranucleotide,
against a correct-heterozygote background of 1,136 / 145 / 308 / 39 / 10 for
periods 2 to 6 — so dinucleotides carry the class about as heavily as they carry
the ground.

### Tract length is the one thing that separates them

Rate of this error per 100 comparable tracts, 30x:

| tract length (repeat units) | homopolymer | period 2+ |
|---|---|---|
| under 10 | 1.1 per 100 (3 of 266) | 0.3 per 100 (2 of 737) |
| 10–13 | 1.1 per 100 (9 of 806) | 1.0 per 100 (5 of 499) |
| 13–16 | 1.5 per 100 (13 of 861) | 2.7 per 100 (10 of 370) |
| 16–20 | 2.5 per 100 (19 of 768) | 4.0 per 100 (17 of 422) |
| 20–25 | 4.0 per 100 (22 of 544) | 3.9 per 100 (14 of 357) |
| 25 and over | 7.4 per 100 (20 of 270) | 1.9 per 100 (3 of 158) |

A homopolymer of 25 A's or more is wrongly called heterozygous about **7 times
in 100**; one under 13 A's, about **1 time in 100**. That is a 6.6-fold rise
across the length range, and it is the strongest gradient anywhere in this data.

---

## 2. What the spurious second allele is

Taking the truth's one sequence and the other member of ng's called pair
(`mechanisms.py`, 30x):

| | homopolymer (86) | period 2+ (51) |
|---|---|---|
| a **different tract length** | 66 | 44 |
| a substitution **inside** the tract | 7 | 7 |
| a substitution in the **single flanking base**, tract length right | 12 | 0 |
| neither called allele is the truth's | 1 | 0 |

**110 of the 137 are a length difference** — 8 in every 10. Of those 110:

| distance | homopolymer | period 2+ |
|---|---|---|
| one repeat **shorter** | 27 | 19 |
| one repeat **longer** | 19 | 13 |
| two or more repeats shorter | 13 | 10 |
| two or more repeats longer | 7 | 2 |

So **78 of the 110 are exactly one repeat unit away** (7 in 10), and short
products outnumber long ones 69 to 41. Not one of the 110 is a length change
that is not a whole number of repeat units. That is the signature of a slip
product being read as a second allele, and the short-over-long imbalance is what
a `shorter_share` of 0.50 — dead even — cannot express.

**Two subclasses a stutter model would not touch.** The 14 substitutions inside
the tract are interruptions or base errors, not lengths. The 12 flanking-base
substitutions are not the tract caller's error at all: ng's tract call is
correct and a SNP-path record sitting on the one base beside the tract carries
the heterozygote — e.g. `chr1:9,955,403`, where the tract is called `1/1` on
14 reads and a neighbouring `0/1` record splits 22 reads to 7 on the base after
a 22-A run. The scorer counts that base because of its one-base anchor pad, and
counting it is right — but the fix lives in the SNP path.

---

## 3. Would a stricter candidate support bar remove them? No.

ng admits a tract allele on `max(2 reads, 10% of that sample's reads at the
locus)`. `bar_sweep.py` re-applies a ladder of bars to the read-level candidate
dump (`tmp/attrib/tier_30x_candidates.tsv` — 20,204 tracts at 30x, of which the
6,630 both sides call are the ones scored here) and counts two things: het-for-hom tracts that lose the spurious allele, and tracts
ng currently gets right that lose an allele the truth carries. The second is
the worse error — a lost true allele cannot be recovered by any genotyper.

At the current setting the sweep reproduces today exactly (0 removed, 0 lost),
which is the check that the link between the dump and the outcome table holds.

| bar | homopolymer: spurious removed / true alleles lost | period 2+: spurious removed / true alleles lost |
|---|---|---|
| 2 reads, 10% (today) | 0 of 86 / 0 of 3,113 | 0 of 51 / 0 of 2,297 |
| 2 reads, 15% | 7 of 86 / **10** of 3,113 | 5 of 51 / 4 of 2,297 |
| 3 reads, 10% | 11 of 86 / **50** of 3,113 | 10 of 51 / **57** of 2,297 |
| 3 reads, 20% | 19 of 86 / **64** of 3,113 | 11 of 51 / **72** of 2,297 |
| 4 reads, 10% | 17 of 86 / **127** of 3,113 | 12 of 51 / **150** of 2,297 |
| 5 reads, 30% | 23 of 86 / **284** of 3,113 | 14 of 51 / **333** of 2,297 |

**There is no setting worth taking.** The mildest one that removes anything
(2 reads, 15%) removes 7 homopolymer errors and destroys 10 correct calls. The
next rung up removes 11 and destroys 50. And the ladder saturates: even at
5 reads and 30% — a bar that costs 617 correct calls across both classes — only
37 of the 137 errors go away, because the spurious allele usually has real
support. Measured in the dump, the spurious allele carries a median of 4 reads
and 32% of the locus's reads at homopolymers, against 12 reads and 54% for a
true non-reference allele. The two distributions overlap heavily.

### The same sweep run downwards prices a discovery round

Loosening the bar from 10% to 5% of the locus's reads would supply the truth's
sequence to **2 more** of the 242 homopolymer tracts where it is currently
missing, and would hand **137 more** correct tracts a candidate the truth does
not carry (51 → 188). About seventy new wrong candidates for every truth
recovered. A discovery feature that works by lowering this bar is not worth
building; one that works by scoring candidates better might be.

*(A lead, not a conclusion: at the current bar, the truth's sequence is already
among the candidates ng kept at 93 of the 242 homopolymer `never_offered`
tracts and 45 of the 165 period-2+ ones — so a substantial part of that larger
class may be lost after candidate selection rather than at it. `loosen_out.txt`
holds the counts. This was not the assigned question and the comparison
between a dump sequence and a VCF record's sequence is approximate; it needs
its own check before anyone acts on it.)*

---

## 4. Does GQ separate them? Not enough, and a gate is the wrong instrument

GQ does separate: median 34.5 for the class against 72 for the heterozygotes ng
gets right (homopolymer). But a GQ gate does not correct a call, it withdraws
it — the wrong genotype becomes a no-call and correct calls are withdrawn with
it. The trade at 30x:

| GQ floor | homopolymer: het-for-hom withdrawn / correct withdrawn | period 2+: same |
|---|---|---|
| 10 | 18 of 86 / 77 of 3,113 | 4 of 51 / 45 of 2,297 |
| 20 | 33 of 86 / 203 of 3,113 | 9 of 51 / 122 of 2,297 |
| 30 | 43 of 86 / 413 of 3,113 | 20 of 51 / 242 of 2,297 |
| 40 | 51 of 86 / 709 of 3,113 | 26 of 51 / 385 of 2,297 |
| 50 | 65 of 86 / 1,051 of 3,113 | 33 of 51 / 614 of 2,297 |

Half the class needs a floor near GQ 30, which withdraws **413 correct
homopolymer calls to withdraw 43 wrong ones** — nearly ten good calls per bad
one, and the bad ones become no-calls rather than right ones.

### An allele-balance rule is no better

The obvious correction — where the two called alleles are too lopsided, call the
sample homozygous for the better-supported one — is a correction rather than a
withdrawal, so it deserves its own curve. `net_balance.py` scores it on
outcomes: a case counts as fixed only when collapsing lands on the truth's
sequence.

| minor share below | homopolymer fixed / correct heterozygotes broken | net |
|---|---|---|
| 0.05 | 2 / 0 | +2 |
| 0.10 | 2 / 1 | +1 |
| 0.15 | 6 / 7 | −1 |
| 0.20 | 13 / 16 | −3 |
| 0.25 | 16 / 49 | −33 |
| 0.30 | 24 / 171 | −147 |

The best net effect anywhere on the curve is **+2 tracts of 3,515**. Part of the
reason is that the truth is the better-supported of ng's two alleles only 42
times in 86 (19 times it is the *worse*-supported, 25 undecidable), so
collapsing onto the majority allele is not even reliably the right move.

---

## 5. Does depth matter? Almost not at all — and that is the point

Rate per 100 comparable tracts by DP, 30x homopolymers: 5.7 (3 of 53) below
DP 15, 3.2 (14 of 444) at DP 20–25, 2.6 (23 of 897) at DP 30–35, 1.1 (2 of 188)
at DP 45 and over. A gentle decline, no concentration.

More telling is what happens when the whole run goes from 30x to 50x:

| class | 30x | 50x |
|---|---|---|
| homopolymer, truth allele never offered | 242 | 219 |
| homopolymer, called hom where truth is het | 27 | 11 |
| **homopolymer, called het where truth is hom** | **86** | **84** |
| period 2+, truth allele never offered | 165 | 146 |
| period 2+, called hom where truth is het | 14 | 8 |
| **period 2+, called het where truth is hom** | **51** | **51** |

Every other error class shrinks when reads are added. **This one does not
move.** An error that more data does not buy down is not a sampling error; it
is a model error, and it will be there at 100x as well.

---

## 6. What would remove them, and at what cost

ng ran with `--defaults`. Its own parameters file says what that means, at
`benchmarks/ssr_hg002/results/ng/HG002_30x.raw.parameters.toml`:

> This table is empty, which is not the same as a missing row: **no stratum was
> fitted at all** … `share_of_reads_that_slip` = 0.10, `shorter_share` = 0.50,
> `fall_off` = 0.05 … **One pair of numbers stands in for every stratum.** A
> 20-base mononucleotide run and a 5-copy tetranucleotide are scored identically
> here, where real slippage rises steeply as the period falls and as the tract
> lengthens — short-period long tracts are where this is furthest wrong.

The measurements above are that prediction coming true. One flat slip rate of
10 reads in 100 is applied to a 30-A homopolymer and an 8-copy dinucleotide
alike; where the true rate is higher, the extra slip reads are too many for the
model to call noise, so it prices them as a second allele — one repeat away,
usually shorter, on 3 to 5 reads.

**The simulator says how much of the class this accounts for.** Re-scored with
the current scorer (`sim/genotype.tsv` here; the copy in `tmp/tract_qual/`
predates the scorer's last fix and should not be quoted):

| reads slipping | model ng used | homopolymer het-for-hom | period 2+ het-for-hom |
|---|---|---|---|
| 10 in 100 | the shipped 10 in 100 — correct | 0 of 333 | 1 of 1,682 |
| 25 in 100 | the shipped 10 in 100 — 2.5x too low | 25 of 334 | 106 of 1,683 |
| 25 in 100 | the true 25 in 100 | 2 of 332 | 7 of 1,677 |

Handing the caller the slippage its reads were actually drawn under removes
**122 of 131** of this error class — 93 in 100. When the shipped model happens
to be right, the class is 1 tract in 2,015.

**Headroom on the real data is smaller than that, and can be bounded.** Only
the 110 length-change cases are prices a slippage model sets. If every tract
length carried the rate of the shortest bin, those 110 would be about 50
(`headroom.py`: homopolymer 66 → 36, period 2+ 44 → 14). So a length-aware
slippage fit should remove somewhere between **60 and 110 of the 137** — the
lower bound if it only flattens the length gradient, the upper if it does on
real reads what it does on simulated ones. The remaining 27 are the 14
in-tract substitutions and the 12 flanking-base SNP-path calls, which no
slippage number touches.

**The cost.** The machinery is already in the caller:
`slippage_by_stratum_and_group` is keyed by stratum, and
`src/ng/parameter_estimation/joint/slippage_curve.rs` already holds how the
level rises with repeat count (`census.rs`'s repeat-count band was widened from
±4 to ±8 on 2026-08-20 with the note that ±4 "was losing real slippage at long
tracts"). What is missing is the command that fits a parameters file and feeds
it back — the benchmark runner says so in as many words ("no command fits a
parameters file yet, so every cell on either GIAB ground is `Defaulted`"). So
the cost is building and validating that fitting pass, plus one extra pass over
the reads per run. It is not a threshold change and cannot be had by tuning one.

### What is not worth doing, plainly

- **A stricter candidate support bar.** Every setting loses more true alleles
  than it removes spurious heterozygotes, and the whole ladder saturates at 37
  of 137 removed for 617 correct calls destroyed. Do not raise it.
- **A GQ floor.** Withdraws about ten correct calls per wrong one, and turns
  the wrong ones into no-calls rather than right ones.
- **An allele-balance collapse rule.** Best net effect on the curve is +2
  tracts in 3,515, and the majority allele is the truth only 42 times in 86.
- **A discovery round that works by lowering the support bar.** It buys 2 truth
  sequences and hands out 137 new wrong candidates at homopolymers. A discovery
  round is only safe once the slippage model can tell a slip product from an
  allele — which is the same fix as above, and is why it should come first.

---

## Scripts

All under `tmp/agent_hetfhom/`, run with `uv run --no-project python` from the
worktree root.

| file | what it does |
|---|---|
| `classify.py <depth>` | re-runs the scorer's comparison, one row a tract → `tracts_<depth>.tsv` |
| `analyse.py` | tract shape, DP/GQ/QUAL/AD quantiles, first cut at the offset distribution |
| `ad_shape.py` | read support behind each called allele, read through the GT indices |
| `mechanisms.py` | the three mechanisms and the repeat-unit distance distribution |
| `bar_sweep.py` | the candidate support bar swept upward, against the read-level dump |
| `loosen.py` | the same bar swept downward — output kept in `loosen_out.txt` |
| `support_gq_depth.py` | spurious-allele support, the GQ trade curve, the depth bins |
| `balance.py`, `net_balance.py` | the allele-balance rule, scored on collapses and on outcomes |
| `length_and_major.py` | is the truth the better-supported allele; the length gradient |
| `headroom.py` | the length gradient restricted to length-change errors, and the counterfactual |
| `sim/genotype.tsv` | the simulator arms re-scored with the current scorer |
