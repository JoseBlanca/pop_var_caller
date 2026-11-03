# pvc's errors inside STR segments: what they are and which to fix first

*Report, 2026-11-02, branch `str-inside-errors` from main e285095f. No caller code was changed. Scratch,
not committed: `tmp/` in the worktree `../pop_var_caller-str-inside-errors` (`tmp/classify/` holds every
script named below). Handoff that set the task: main checkout `tmp/handoff/str_inside_errors/context.md`.*

---

## 1. The answer

**At 30× on HG002, most errors inside STR segments come from the stutter model: pvc confuses a real
allele with stutter, in both directions.** That is 123 of the 254 errors on chromosome 20, at 76 of
the 162 segments that hold an error. The next largest class is false alleles that replace an
interrupting base inside a repeat (39 errors). The separate context-error work already targets that
class. Depth accounts for only 28 of the 254.

**At 5×, 92% of the errors inside segments are depth.** In a random sample of 478 errors, 440 sit at
segments where at most one read carries the true allele, or where a call was made on three reads or
fewer.

**The stutter class has two causes, and only one is the model's form:**

- **`--defaults` scores every repeat with one stutter rate**: 5% of reads one unit longer and 5% one
  unit shorter, whatever the motif or the length. These are HipSTR's shipped constants. HG002 is
  PCR-free, and its reads show 0–4 such reads in 100 at dinucleotides and short homopolymers, but 17–28
  in 100 at homopolymers of 15 copies or more. One rate is wrong in both directions at once.
- **Fitting the rates on HG002 does not repair it.** The fit gets homopolymers and dinucleotides
  right, but it says that 23% to 65% of reads slip at repeats with a 4- or 5-base motif. Real
  alleles are then lost there instead. The same fit, run on the stuttering simulation, recovers the
  true rates to within about 10%. So something in HG002's real repeats misleads it, and that something
  is not yet known.

**Recommendation (§8): find the root cause of the fit's inflated rates at 4- to 6-base motifs first,
before building anything.** A correct per-motif, per-length stutter is what the largest class needs. On
the simulation it already cuts one sample's errors at 30× from 133 to 23. On HG002 it cannot be judged
until the fit is right.

## 2. What was measured, and how

**Data.**

- **HG002 chromosome 20** against GIAB v4.2.1, at 30× and thinned to 5×.
- **GIAB's tandem-repeat set** (HG002 TR v1.0.1, 50,000 regions) at the same two depths, for the
  overall numbers only.
- **The stuttering simulation** (`examples/ng_str_stutter_benchmark.rs`), with PCR stutter, at the
  two ends of the range:
  - one sample at 30×;
  - a cohort of 50 inbred samples at 3×, drawn from 8 founder haplotypes, with 80% of each genome in
    runs where both copies are the same.

**Two runs of pvc, both from main e285095f.**

- **`--defaults`**: the configuration the benchmarks use. Nothing is fitted; the parameters are
  pvc's stated defaults.
- **fitted**: `generate-psps`, then `estimate-parameters` on the same sample(s), then
  `call-from-psps` with that fit.

**Scoring.**

- **HG002**: vcfeval with genotype matching. An error is a *missed* variant (in the truth, not
  called right) or a *false* one (called, not in the truth). It is placed against pvc's own STR
  segments (`ssr_locus` or `ssr_bundle` in `ng_typed_region_dump`) by `by_distance.py`.
- **Simulation**: `scripts/ng_str_stutter_score.py`. One unit is one sample at one STR segment, and
  an error is any outcome other than *right* or *ref_right*.

**How each error was classified.** pvc was rerun with `--explain-loci` over the error segments; the
calls were checked byte for byte against the run being classified. The explanation lists:

- the candidate alleles;
- the reads that spell each allele in full (*complete reads*: reads crossing the whole segment);
- the likelihood of every genotype;
- the call and its genotype quality (GQ);
- what became of the record.

The truth's two sequences over the segment were rebuilt from GIAB's records and compared with those
facts. Each segment then got exactly one class, by the first rule below that applies:

| class | rule |
|---|---|
| depth | the missed truth allele is on fewer than 2 complete reads and fewer than 8 reads cross (a heterozygote shows on at most one of 8 reads about 1 time in 30), or a wrong allele was called on 3 crossing reads or fewer; or no locus was built and at most one spanning read in the BAM carries the variant |
| allele absent from the reads | enough reads cross, but at most one spells the truth allele |
| real allele taken as stutter | the truth is heterozygous, the truth allele is a candidate, and the call is homozygous for its other allele |
| stutter taken as an allele | the truth is homozygous and the call adds an allele a whole number of units away |
| length off | a heterozygote called as another heterozygote of whole-unit lengths |
| proposed too late | the truth allele is in reads but was not made a candidate (only whole-unit cases occur) |
| false substitution allele | an allele of the reference length with a base changed inside the repeat, not in the truth |
| non-unit indel or real substitution | an allele whose length change is not a whole number of units, or a true substitution called as something else |
| seam or edge | a truth record crosses the segment's edge or overlaps another one (the unfinished design of `str_flank_ownership.md`; counted, not studied) |
| right genotype | the locus's call is the truth; the error is in scoring or in another record |

The 5× classification is of a random sample: 330 errors were drawn, and every error at the segments
they fell in was kept, giving 478 errors at 328 segments.

## 3. The baseline on current main

`--defaults`, HG002 chromosome 20. Each cell is SNP / indel. "Overlapping" means the variant's
changed bases touch a segment.

| where | 30× missed | 30× false | 5× missed | 5× false |
|---|---|---|---|---|
| overlapping a segment | 16 / 99 | 0 / 139 | 326 / 2,758 | 0 / 788 |
| 1–8 bases from one | 8 / 6 | 1 / 14 | 257 / 85 | 58 / 17 |
| 9–20 bases from one | 1 / 4 | 2 / 0 | 271 / 48 | 75 / 6 |
| farther | 469 / 79 | 294 / 23 | 19,208 / 1,495 | 3,496 / 179 |
| inside a satellite (not called by default) | 14 / 57 | 0 / 0 | 14 / 57 | 0 / 0 |

These are within a handful of the handoff's numbers from 0e635332. The tandem-repeat set stands at
14,245 missed and 2,185 false at 30×, and 24,028 / 2,537 at 5×. It was not split by segment.

## 4. The classes on HG002

`--defaults`. 30× covers every error overlapping a segment; 5× is the sample of §2.

| class | 30× segments | 30× missed / false | median GQ | 5× segments (sample) | 5× missed / false |
|---|---|---|---|---|---|
| **depth** | 13 | 16 / 12 | 9 | 303 | 340 / 100 |
| **real allele taken as stutter** | 23 | 23 / 16 | 7 | — | — |
| **stutter taken as an allele** | 38 | 17 / 38 | 10 | 10 | 6 / 10 |
| **length off** | 8 | 8 / 8 | 12 | — | — |
| **proposed too late** | 7 | 8 / 5 | 13 | — | — |
| *whole-unit total* | *76* | *56 / 67* | | *10* | *6 / 10* |
| **false substitution allele** | 36 | 3 / 36 | 18 | 1 | 0 / 1 |
| **non-unit indel or real substitution** | 13 | 11 / 12 | 26–34 | — | — |
| **allele absent from the reads** | 15 | 20 / 4 | 44–63 | 5 | 5 / 0 |
| **seam or edge** | 7 | 7 / 6 | 19 | 7 | 10 / 3 |
| **right genotype** | 2 | 2 / 2 | 4 | 2 | 1 / 2 |
| total | 162 | 115 / 139 | | 328 | 362 / 116 |

At 5×, "allele absent from the reads" means a locus was never built even though at least two spanning
reads carry the variant. That rule needs 2 reads disagreeing with the reference. Rebuilt from the BAM
for every unbuilt locus at 5×: at 94, no spanning read carries the variant; at 92, one does; at 5, two
or more do. A one-read rule was tried on 2026-10-10 and rejected: it made as many false calls as right
ones at 5×.

### 4.1 Whole-unit confusion between a real allele and stutter

Seventy-six segments. Most are homopolymers (45) or dinucleotides (18); 13 have a 3- to 6-base motif.
Of the 22 dropped alleles traced for direction, 20 are one unit from the called allele (11 longer, 9 shorter) and 2 are two units shorter.

**At chr20:734,611**, a GT repeat of 14 copies, the truth is a heterozygote of +3 and +4 units. Nine
reads show +3 and two show +4. pvc calls +3 on both copies, with GQ 5.

**At chr20:3,383,802**, a homopolymer of 15 T, the truth is the reference plus one T on one copy.
Twenty-two reads show the longer allele and 4 the reference. pvc calls the longer allele on both
copies, with GQ 20.

**The reads favour the wrong genotype at all 23 dropped alleles under `--defaults`.** The fixed 5%
rate makes 2 reads in 11 cheaper to explain as stutter than as a heterozygote. It is also wrong for
HG002. Counting, at 395 segments where GIAB has a homozygote and pvc called right (30×, at least 8
complete reads), the reads one unit off the true length:

| motif length | copies | reads | one unit longer, per 100 | one unit shorter, per 100 |
|---|---|---|---|---|
| 1 | under 10 | 1,426 | 1.1 | 0.6 |
| 1 | 10–14 | 3,275 | 3.2 | 2.5 |
| 1 | 15–19 | 2,054 | 11.3 | 6.2 |
| 1 | 20 or more | 739 | 9.5 | 18.1 |
| 2 | any | 983 | 0.3 | 1.9 |
| 3–6 | any | 1,091 | 7.2 | 7.2 |

At 3- to 6-base motifs HG002's reads do show about 7 in 100 one unit off, more than at dinucleotides;
that is still a third or less of what the fit puts there (below).

**A minority allele's share of reads does not separate the two errors.** At the dropped real alleles
it is 12–21% of complete reads; at the stutter called as an allele, 8–33%. Real heterozygotes with a
share that low are not rare: of 351 whole-unit heterozygotes pvc called right at 30× (a random
sample), 23 have the rarer allele on fewer than 20% of reads.

**Fitted on HG002** (chromosome 20, 30×), the class shrinks only slightly, from 76 segments to 73 and
from 123 errors to 116:

- stutter taken as an allele falls from 38 segments to 30;
- real allele taken as stutter rises from 23 to 30.

Inside segments, missed variants go 115 → 118 and false ones 139 → 123. Over all of chromosome 20,
missed go 765 → 851 and false 450 → 378. At 5×, the fit is slightly worse everywhere.

Why the dropped alleles rise: the fit puts 23% of reads slipping at a 9-copy repeat of a 5-base motif,
and 51% at 14 copies. At a 4-base motif of 25 copies it puts 65%. A PCR-free library does nothing of
the kind. Homopolymers and dinucleotides come out plausible: 3 reads in 1,000 at 8 copies, 46 at 16
copies, about 130 at 28 copies or more.

**chr20:41,015,365**, (GTTTT)×9, shows the result. Ten reads carry the reference length and eight
carry one more unit. Even so, "reference on both copies" scores only 2.6 log units below the true
heterozygote, and with the prior pvc calls 0/0. At 19 of the 30 dropped alleles in the fitted run, the
reads favour the true genotype and the call goes the other way.

**The genotype prior is not what is wrong.** The fitted prior is the share of chromosomes at each length
offset from the reference, for each motif length and copy number. Set beside GIAB's own genotypes at
chr20's segments, it agrees at short repeats: homopolymers of 8 copies have 98.0% of chromosomes at
the reference length in the truth and 96.8% in the prior; dinucleotides of 8 copies, 91.9% and 90.3%.
At long homopolymers the prior is *more* generous to variants than the truth. Script:
`tmp/classify/prior_check.py`.

**On the simulation the same fit is right.** One sample at 30×, PCR stutter: at a 5-base motif of
9 copies the true rate is 0.76% times that library's multiplier of 1.24, which is 0.94%, and the fit
says 1.07%. What makes HG002's 4- to 6-base repeats different is not known yet. Real heterozygotes being
counted as slips in strata with few loci in one sample is a guess, not a measurement.

### 4.2 False substitution alleles inside a repeat

Thirty-six segments: 16 bundles, 9 homopolymers, 11 others. In 24 of the 37 false alleles, the
reference repeat carries one interrupting base and the false allele replaces it with the repeat's own
base: `AAAAAA[G>A]AAAAAA`, `TTTTTT[C>T]TTTTTT`, `GTGTGT[C>G]TCTGTG`. Two to five of 10–25 reads
carry it. **This is the mechanism the context-error work measures**
(`doc/devel/ng/spec/str_context_error.md`, worktree `../pop_var_caller-context-error`, not in main): a
read sequenced through a long run of one base misreads the next base as one more copy of the run. Its
prototype removed 20 of the 138 false calls inside segments at 30× on the 2026-11-01 build. Periods 2
to 4 are its next step, by the owner's ruling. At 5× the class is 1 of 478 errors: the false allele
needs about three reads.

### 4.3 The smaller classes

- **Non-unit indels and real substitutions** (13 segments, 23 errors): 7 of the 9 false non-unit
  alleles are at bundles. Not traced further.
- **Allele absent from the reads** (15 segments, 24 errors at 30×): reads cross,
  but none or one spells the truth allele. GQ is high (44–63), so pvc is confident. The candidates are
  a truth error, mapping (the allele's reads aligned elsewhere or soft-clipped), or a long insertion
  the reads cannot span. It needs tracing read by read before it can be called a caller class.
- **Seam or edge** (7 segments, 13 errors): counted; it belongs to the unfinished flank design.

## 5. The range

The simulation, both ends, fixed stutter defaults against fitted. The unit is one sample at one STR
segment.

**One sample at 30×** (3,000 cells):

| class | `--defaults` | fitted |
|---|---|---|
| stutter taken as an allele | 106 | 4 |
| real allele taken as stutter | 7 | 10 |
| length off | 7 | 2 |
| depth, or the allele not proposed | 6 | 6 |
| flank change in the truth allele (SNP beside, deletion across the edge) | 7 | 1 |
| total | 133 | 23 |

**Fifty inbred samples at 3×** (150,000 cells):

| class | `--defaults` | fitted |
|---|---|---|
| no call (./.) | 3,969 | 3,618 |
| depth: the sample's reads lack its allele | 9,401 | 1,515 |
| whole-unit confusion (all at 3 reads or fewer, but 115 fitted) | 13,706 | 583 |
| flank change in the truth allele | 2,026 | 503 |
| allele in reads, not proposed | 61 | 58 |
| substitution allele involved | 192 | 2 |
| total | 29,355 | 6,279 |

What this says:

- **At 3× in a cohort, 82% of the fitted run's errors are depth**: the no-calls (3,719 of all 3,803
  no-calls, flank cases included, have no complete read in that sample) and the samples whose reads
  lack their allele. Whole-unit confusion is 9%.
- **Most of the fit's gain at 3× is not the stutter.** The run with the simulator's true stutter and
  every other parameter at its default still makes 24,829 errors. The rest comes from fitting the
  inbreeding and the allele prior: with `--defaults`, an 80% inbred cohort is scored as outbred.
- **At one sample at 30×, the fitted stutter removes 102 of 106 "stutter taken as an allele" errors**
  and adds 3 dropped real alleles.

The simulator had to be repaired first; uncommitted, in the worktree. Its copy of the repeat scanner was
out of date since the scanner change of 2026-10-30, and it laid 19 of 500 STR bundles that pvc types
as something else. A subagent ported the split into the copy: all 24 of its tests pass, and every
segment and bundle now types as planted. The driver
(`scripts/ng_str_stutter_benchmark.sh`) and the suite's `stutter_sim.py` call `ng_str_bundle_runs`,
an example deleted in 0f9281c5, so the stutter step is broken on main. It is restored, untracked, in
the worktree.

The tomato cohort's rick results were not split by segment; that waits on the owner.

## 6. Which errors are the caller's

| | 30× HG002 (254) | 5× HG002 (sample of 478) | 3× sim cohort, fitted (6,279) | 30× sim, fitted (23) |
|---|---|---|---|---|
| depth (no-calls included in the cohort) | 28 | 440 | 5,133 | 3 |
| whole-unit confusion (stutter model) | 123 | 16 | 583 | 16 |
| false substitution allele | 39 | 1 | 2 | 0 |
| non-unit or substitution | 23 | 0 | — | — |
| allele absent from reads, or not proposed | 24 | 5 | 58 | 3 |
| seam or edge / flank change | 13 | 13 | 503 | 1 |
| right genotype | 4 | 3 | — | — |

## 7. What fixing the two largest caller classes would take

**Whole-unit confusion: a stutter rate per motif length and copy number that is right for the
library.**

- *What a fix looks like.* Score STRs with stutter by motif length and copy number, fitted from the
  run's own reads, instead of one shipped rate. Before that, find why the one-sample fit on HG002
  gives 23–65% slips at 4- to 5-base motifs. Is it real heterozygotes counted as slips, interrupted
  repeats, few loci, or something else?
- *What it touches that was ruled out before.* The tract-accuracy programme of 2026-09-03 discarded
  "fitted per-stratum stutter rows as default": they bought heterozygotes at motifs of 2 bases or more
  by selling homopolymers. The new data: on the simulation the fit is now right at every motif length;
  on HG002 its error is concentrated at 4- to 6-base motifs; and the fixed default (10 in 100 one unit
  off) is about 4 times what HG002's dinucleotides show (2.2 in 100).
- *Cost at the other end.* In a cohort at 3×, the fit already gets the stutter right (fitted
  pentanucleotide rates 0.1–1.7%; the simulator's PCR rates 0.1–1.5% before each library's own
  multiplier),
  and depth dominates. A change made for one deep sample must not loosen it there; the cohort grid
  above is the check. At 5× on HG002 the class is 16 errors of 478.

**False substitution alleles: the context-error term.** It already has a design and a prototype, and the owner has
set its next steps. Nothing new is proposed here; this report only sizes it — 39 of 254 errors at 30×,
none at 5×.

## 8. Recommendation

**First, find the root cause of the inflated slip rates in the one-sample fit at 4- to 6-base motifs
on HG002. This is diagnosis, not a fix.** It is the step between the largest caller class and any
fix for it.

- If the fit is wrong at those motifs, repairing it serves every input.
- Once the fit is right, the choice is between running the fit by default and shipping per-motif,
  per-length defaults. That choice can be measured on both benchmarks. Today it cannot: on HG002 the
  fit trades one error for another.

The second-largest class already has its own line of work.
