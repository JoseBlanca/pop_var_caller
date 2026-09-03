# What a fitted stutter model buys on real reads: HipSTR against ng on GIAB HG002

ng scores every repeat tract under one pair of shipped stutter constants
(HipSTR's: 10 reads in 100 slip, half of those short, 5 in 100 of the slips move
more than one repeat). HipSTR fits a stutter model **per locus** by
expectation-maximization over the read lengths before it genotypes. Both were
run on the same HG002 reads, and both callsets are on disk, so the question
"what would fitting be worth" has an outside answer on real data rather than on
a simulator.

This scores HipSTR with `benchmarks/lib/tract_qual_experiment.py` — the same
instrument, the same truth set, the same confident regions and the same tract
ground ng was scored on.

**The short answer.** Where both callers put a record at the same tract, at 30x
and at period 2 or more, they are within a third of a percentage point of each
other: ng 1,832 right of 1,974 (0.928), HipSTR 1,733 right of 1,874 (0.925). The
gap is not in candidate selection — ng's allele sets miss the truth's sequence
*less* often than HipSTR's (4.0% of tracts against 4.8%). It is in the
genotyper: over sets that did contain the right sequences, ng picks wrong at
3.19% of tracts and HipSTR at 2.72%, and that gap **widens** with depth rather
than closing. So the fitted model is worth something, it is small, and it is a
model story rather than a statistical-power story.

**And the reach is not close, the other way.** HipSTR was handed a region set
with no homopolymer in it at all, and it reaches 2,116 of the truth's 2,827
period-2-or-more tracts against ng's 2,660. Counting the tracts each caller
gets *right* on the whole ground: at period 2 and above, ng 2,297 and HipSTR
1,778; at period 1, ng 3,113 and HipSTR nothing, because it never attempts one.

---

## 1. HipSTR has no homopolymer records here, so ng's homopolymer numbers have no comparator

The `--regions` BED HipSTR ran on
(`results/hipstr/HG002_Tier.hipstr_regions.bed`, built from this project's own
`ssr-catalog` by `catalog_to_hipstr_bed.py`) holds 13,272 loci:

| period | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| loci | **0** | 3,689 | 1,890 | 4,772 | 2,429 | 492 |

HipSTR was never asked to genotype a homopolymer. The scorer nevertheless files
about a hundred HipSTR records under `homopolymer`, because it charges a record
to whichever ground tract its span touches and ng's ground has homopolymers
abutting HipSTR's dinucleotide and tetranucleotide regions. Every one of those
records carries its own `PERIOD` of 2 to 6 (60 at period 2, 63 at 3, 252 at 4,
134 at 5, 30 at 6) — not one is a period-1 record.

**Every `homopolymer` row in HipSTR's tables below is therefore a mislabel and
must not be read as homopolymer accuracy.** ng's homopolymer figures — 3,113
right of 3,515 comparable tracts at 30x, 0.886 — stand without a comparator on
this benchmark. Everything that follows is period 2 and above.

## 2. Reach: how many of the sample's repeat genotypes each caller even attempts

The tier tract ground holds 12,441 homopolymer and 7,763 period-2-or-more
tracts. The truth set calls a variant at 3,817 and 2,827 of them.

| | homopolymer | period 2+ |
|---|---|---|
| tracts the truth calls | 3,817 | 2,827 |
| ng 30x writes a record there | 3,653 (95.7%) | 2,660 (94.1%) |
| HipSTR 30x writes a record there | 193 — all spillover | 2,116 (74.9%) |
| HipSTR 30x with a non-reference genotype | 115 | 2,042 |

(These counts are of tracts holding a truth record. The scorer's own
`tracts_truth_calls` is 5 lower in each class — 3,812 and 2,818 — because it
also requires the truth's genotype to be readable, and 10 tracts across the
ground carry a truth record whose genotype is not. Percentages here use the
row's own denominator.)

HipSTR reaches three quarters of the period-2-or-more tracts the truth calls,
and of the 12,947 records it writes at 30x only 3,163 give the sample a
non-reference allele; the rest are hom-ref (9,444) or no-call (340).

## 3. Genotype accuracy, whole tract ground, 30x and 50x

`comparable` is the scorer's `tracts_both_call` less `not_comparable`, and
`accuracy` is `genotype_right / comparable` — the definition ng's standing
numbers use. `hipstr` is every HipSTR record; `hipstr_var` is only the records
giving the sample a non-reference allele, which is ng's emission policy.

Period 2 and above:

| arm | depth | comparable | right | accuracy | no-call |
|---|---|---|---|---|---|
| ng | 30x | 2,543 | 2,297 | 0.903 | 0 |
| hipstr_var | 30x | 1,920 | 1,778 | 0.926 | 0 |
| hipstr | 30x | 1,947 | 1,778 | 0.913 | 44 |
| ng | 50x | 2,570 | 2,347 | 0.913 | 0 |
| hipstr_var | 50x | 1,961 | 1,855 | 0.946 | 0 |
| hipstr | 50x | 1,973 | 1,855 | 0.940 | 51 |

Read alone this says HipSTR is 2.3 points better at 30x — but over 1,920 tracts
against ng's 2,543. Section 4 removes that.

Homopolymer, ng only (HipSTR has no comparator, §1):

| arm | depth | comparable | right | accuracy |
|---|---|---|---|---|
| ng | 30x | 3,515 | 3,113 | 0.886 |
| ng | 50x | 3,542 | 3,170 | 0.895 |

## 4. Head to head on the tracts both callers reach

`ground/shared_fixed.bed` holds the 2,235 tracts every arm reaches at both
depths — the intersection of the tracts ng writes a record at (30x and 50x) and
the tracts HipSTR writes a record at (30x and 50x). It is fixed across depths,
so a change from 30x to 50x is a change in accuracy and not a change in which
tracts were attempted. 2,044 of those tracts are period 2 or above and the
truth calls a variant at every one.

| arm | depth | comparable | right | accuracy |
|---|---|---|---|---|
| ng | 30x | 1,974 | 1,832 | **0.9281** |
| hipstr_var | 30x | 1,874 | 1,733 | **0.9248** |
| ng | 50x | 1,973 | 1,852 | **0.9387** |
| hipstr_var | 50x | 1,890 | 1,792 | **0.9481** |

At 30x ng is ahead by 0.33 percentage points — about 6 tracts in 1,900. At 50x
HipSTR is ahead by 0.94 points, about 18 tracts. Neither margin is a different
class of caller.

(The 174 homopolymer tracts of this shared set are homopolymers abutting
HipSTR's regions, and ng scores 124 of 161 right at 30x, 0.770 — well below its
0.886 over all homopolymers. Homopolymers sitting inside compound repeats are
harder; HipSTR still has nothing to compare against there.)

## 5. Where the errors are: candidate selection or the genotyper

The scorer's four error counters partition every wrong genotype by what would
have to change to fix it. `truth_allele_never_offered` means a sequence the
truth carries was not among the ones the caller's records name — **no genotype
over that allele set could have been right, so it is candidate selection's**.
The other three are errors made over a set that did hold the right sequences,
so they are the genotyper's: likelihood and prior picked the wrong pair.

Period 2 and above, on the fixed shared ground of §4:

| arm | depth | comparable | errors | never offered | called hom, truth het | called het, truth hom | other |
|---|---|---|---|---|---|---|---|
| ng | 30x | 1,974 | 142 | 79 | 10 | 42 | 11 |
| hipstr_var | 30x | 1,874 | 141 | 90 | 14 | 36 | 1 |
| ng | 50x | 1,973 | 121 | 64 | 6 | 40 | 11 |
| hipstr_var | 50x | 1,890 | 98 | 61 | 6 | 30 | 1 |
| hipstr_var | 300x | 1,893 | 80 | 46 | 1 | 32 | 1 |

As rates over the comparable tracts:

| | ng 30x | HipSTR 30x | ng 50x | HipSTR 50x | HipSTR 300x |
|---|---|---|---|---|---|
| candidate misses | **4.00%** | 4.80% | **3.24%** | 3.23% | 2.43% |
| genotyper errors | 3.19% | **2.72%** | 2.89% | **1.96%** | 1.80% |

**ng's candidate sets are not the problem relative to HipSTR's.** At 30x ng
fails to offer the truth's sequence at 79 tracts in 1,974 and HipSTR at 90 in
1,874 — ng is better by 0.8 points. At 50x the two are level (3.24% against
3.23%). Whatever else fitting buys, it does not buy a wider allele set here.

**The gap is in the genotyper, and it does not close with depth.** ng's error
rate over sets that held the right sequences is 1.17× HipSTR's at 30x (3.19%
against 2.72%) and 1.47× at 50x (2.89% against 1.96%). The absolute gap grows
from 0.47 points to 0.93 points as depth doubles. A gap that closes with depth
is statistical power; this one widens, which says ng is not converting extra
reads into genotype accuracy as well as a fitted model does.

**In absolute terms the prize is small.** Closing the whole genotyper gap at
50x on this shared ground is 18 tracts of 1,973.

**Two shapes inside the genotyper class are worth naming.**

- *Called heterozygous where the truth is homozygous* — a slip product admitted
  as a real allele — is the largest single error for both callers and does not
  shrink with depth: ng 42 at 30x and 40 at 50x, HipSTR 36 and 30 (and still 32
  at 300x). This is the failure a stutter model exists to prevent, and neither
  caller prevents it; HipSTR is about a quarter better at it.
- *Wrong some other way* — both sides heterozygous, the wrong pair of alleles —
  is ng 11 at both depths against HipSTR's 1. Small, but it is entirely ng's,
  and it is not a stutter-pricing failure.

## 6. Depth trend

HipSTR alone across every depth on disk, period 2 and above, on the fixed
shared ground of §4 (variant records only):

| depth | comparable | right | accuracy |
|---|---|---|---|
| 5x | 357 | 275 | 0.770 |
| 10x | 1,153 | 968 | 0.840 |
| 15x | 1,579 | 1,401 | 0.887 |
| 20x | 1,767 | 1,604 | 0.908 |
| 30x | 1,874 | 1,733 | 0.925 |
| 50x | 1,890 | 1,792 | 0.948 |
| 300x | 1,893 | 1,813 | 0.958 |

ng has runs at 30x and 50x only, so the pair can be read at two depths.
Between them ng gains 1.06 points (0.9281 → 0.9387) and HipSTR gains 2.33
(0.9248 → 0.9481). HipSTR then gains a further 0.96 points from 50x to 300x, so
it is close to its own ceiling at 50x; ng's 50x accuracy sits where HipSTR's
is at about 35x.

**Reading it against §5**: the two callers' *candidate* miss rates converge as
depth rises (4.00% against 4.80% at 30x, 3.24% against 3.23% at 50x) while
their *genotyper* error rates diverge (0.47 points apart at 30x, 0.93 at 50x).
So the part of the gap that closes with depth is the candidate part, and the
part that does not is the model part.

## 7. What HipSTR's fit actually looks like on this data

HipSTR writes the parameters it fitted into each record's INFO.
`INFRAME_UP` — the chance a read gains one repeat unit — over the 12,947
records at 30x:

| min | p10 | median | p90 | max | distinct values |
|---|---|---|---|---|---|
| 0.020 | 0.030 | 0.040 | 0.070 | 0.330 | 29 |

ng's shipped constant is 0.05 for a gain (10 in 100 slip, half of them short).
**The median locus's fitted gain rate, 0.04, is within a fifth of the constant
ng assumes** — which is the mechanical reason fitting buys so little on
average. What fitting captures is the tail: at the 90th percentile the fitted
rate is 0.07 and at the extreme 0.33, six times what ng assumes, and those are
the loci where a constant should cost genotypes.

## 8. What is NOT controlled — read the numbers with these

1. **Different tract catalogs.** ng's ground is `ng_typed_region_dump`'s
   `ssr_locus` rows over the Tier intervals: 20,204 tracts. HipSTR's regions come
   from this project's `ssr-catalog`: 13,272 loci, none of them period 1.
   HipSTR's records land on 4,660 of the 20,204 ground tracts. The shared ground
   of §4 is 2,235 tracts — 11% of the ground and 34% of the tracts the truth
   calls. **Nothing here measures either caller outside that slice.**
2. **HipSTR writes no QUAL.** Its QUAL column is `.` at every record and the
   scorer drops a query record without one, so the whole callset read as empty.
   `hipstr_add_qual.py` fills the column from the per-sample posterior `Q`, as
   `-10 log10(1-Q)` capped at 60. The genotype comparison never reads QUAL, so
   §§3–6 are unaffected; but **no calibration or threshold-sweep number can be
   produced for HipSTR this way** and none is reported — those outputs went to
   `/dev/null`.
3. **Different emission policies.** HipSTR writes a record at every region it
   attempts, hom-ref and no-call included; ng writes only variant records. Scored
   on all HipSTR records, HipSTR is charged for hom-ref calls at tracts the truth
   calls variant while ng is simply absent there. Both arms are reported
   (`hipstr` and `hipstr_var`) and they differ by 1.3 points at 30x.
4. **HipSTR ran with two non-default flags** — `--use-unpaired` (the BAMs were
   sliced to the Tier intervals, so mates are gone) and `--min-reads 5` (its
   default of 100 is a cohort threshold). `--min-reads 5` means HipSTR declines
   tracts below 5 reads that ng attempts, which flatters HipSTR's accuracy and
   costs it reach.
5. **More than the stutter model differs.** The two callers use different
   candidate rules, different priors, different aligners, and HipSTR
   marginalizes over haplotype alignments. **A gap of 18 tracts at 50x cannot be
   attributed to the stutter fit alone** — it is the whole genotyping stack.
   One confound is ruled out: HipSTR did no SNP-based read phasing here (`DSNP`
   is 0 at every record, since it ran single-sample with no `--snp-vcf`).
6. **The shared ground is defined by where both callers emit**, which is
   enriched for tracts ng called variant — 2,222 of its 2,235 tracts are truly
   variant. It is not a random slice of repeat tracts.
7. **One sample, one species, one depth ladder.** Everything here is HG002 on
   GRCh38 at a single sample. Nothing in it speaks to cohorts, and per the
   project's own range commitment a figure measured on one high-coverage human
   is a fact about that corner.

## 9. Verdict

**A fitted stutter model is worth roughly one percentage point of genotype
accuracy on real HG002 reads at 50x, and nothing measurable at 30x, on the
tracts both callers reach.** ng at 30x is 0.33 points *ahead* of HipSTR on the
shared ground; at 50x it is 0.94 points behind. That is 18 tracts of 1,973.

**The gap that exists is a model problem, not a candidate problem.** ng's
allele sets contain the truth's sequences as often as HipSTR's (better at 30x,
level at 50x); ng's genotyper then picks wrong 1.2 to 1.5 times as often over
sets that held the right answer. The gap in that class widens as depth doubles,
so it will not be fixed by more reads.

**But the size of the prize argues against fitting as the next move.** Two
things are larger by an order of magnitude on the same tables:

- **Homopolymers.** ng is right at 3,113 of 3,515 comparable homopolymer tracts
  at 30x (0.886) against 0.903 at period 2 and above, and 402 wrong homopolymer
  tracts is 22 times the 18-tract genotyper gap HipSTR demonstrates at period 2
  and above. HipSTR offers no evidence on homopolymers at all here, because it
  was handed no period-1 locus.
- **ng's own candidate misses in absolute terms.** On the whole ground at 30x,
  242 homopolymer and 165 period-2-plus tracts fail because the truth's sequence
  was never offered. That is 407 tracts against the 18 a perfect stutter fit
  would buy at period 2 and above.

The one piece of evidence that a per-locus fit would matter where the constants
are worst is §7's spread: HipSTR's fitted gain rate runs from 0.02 to 0.33 with
a median of 0.04 against ng's fixed 0.05. The median locus does not need
fitting; the top decile (0.07 and up) is where a constant is wrong by enough to
cost genotypes, and that is a targeted fit at high-slip tracts rather than a
per-locus EM everywhere.

---

## Files

| file | what it is |
|---|---|
| `hipstr_add_qual.py` | fills HipSTR's empty QUAL column from the per-sample posterior; writes an all-records and a variant-only VCF |
| `build_shared_ground.py` | the tracts every named callset puts a record at, as a ground BED; imports the scorer's own record-to-tract rule |
| `run_all.sh` | every arm on the whole tract ground, and on a per-depth shared ground |
| `run_fixed.sh` | every arm and depth on one fixed shared ground |
| `reach.py` | tracts reached per caller, against the tracts the truth calls |
| `homopolymer_check.py` | the period-1 question of §1 |
| `out/genotype_full.tsv` | §3 — whole tract ground |
| `out/genotype_shared.tsv` | per-depth shared ground (superseded by the fixed one) |
| `out/genotype_fixed.tsv` | §§4–6 — one fixed shared ground |
| `ground/shared_fixed.bed` | the 2,235 tracts of §4 |
| `vcf/HG002_*.vcf`, `vcf/HG002_*.var.vcf` | HipSTR with a QUAL column |
