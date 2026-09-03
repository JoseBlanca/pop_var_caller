# ng at repeat tracts: how we measure it, what moves it, and what to try next

**Date:** 2026-09-03. **Written to be read by somebody starting fresh**, so it repeats things the
last conversation took for granted. It covers three questions the owner asked for: how the numbers
are produced (§2, and **§3 is the trap** — read it before trusting any genotype figure), what has
been tried and what it was worth (§4, §5), and what to do next (§6).

**Branch** `ng-ssr-calling-loop`, worktree `/Users/jose/devel/pop_var_caller-ssr-calling-loop`,
six commits ahead of `main`, working tree clean, 6,026 library tests green.

---

## 1. What ng is doing at a repeat tract, in one paragraph

ng is the from-scratch SNP/indel/STR caller in `src/ng/`. A **repeat tract** is a microsatellite —
a stretch like `ATATATATAT` — and ng routes those to their own path: candidate alleles chosen by
`select_ssr`, evidence shaped by `shape_ssr_locus`, and genotypes from the same calling loop the
SNP path uses but under a **stutter model** instead of a base-error model. Stutter (also called
slippage) is a read reporting a different number of repeat copies than the allele it was copied
from; it is the dominant error at a tract and everything here turns on how it is priced.

Every benchmark run of ng today is `--defaults`, which means **nothing is fitted**: the stutter
model is HipSTR's shipped constants, one pair of numbers for every class of tract, and the
genotype prior's length spectrum is flat. The run's own report says so — *"0 of 7 groups the file
says were fitted"*.

---

## 2. How the numbers are produced

### 2.1 The three grounds

| ground | what it is | why it exists |
|---|---|---|
| **`benchmarks/ssr_hg002`** — GIAB's HG002 tandem-repeat benchmark, 50,000 Tier intervals over 6.09 Mb, 36,497 assembly-based truth records, BAMs at 5×–300× | ng types **20,204 repeat tracts** in it and writes **6,351 tract records at 30×** | **This is the ground that can answer anything.** One sample only, so it says nothing about cohorts |
| **`benchmarks/giab/per_sample`** — the GIAB trio's 100 random confident intervals a sample, at 5×–300× | ng writes **149 tract records pooled over three samples at 30×** | Where every standing ng number was measured, so a figure here is comparable with C4's. **Far too small for anything about tracts** — its calibration table is a column of zeroes |
| **`examples/ng_tract_simulator`** — tracts whose genotypes we chose, sequenced under a slippage we set | exact truth, settable slippage, and the only place the fitted-versus-defaulted split has two sides | **It cannot measure candidate-selection losses**: its "truth allele never offered" count is 0 at every tract, because its alleles are drawn within three repeats of the reference and its reads carry nothing else |

**A trap in the third one.** The simulator says a wrong stutter model costs a lot (period-2+
genotype accuracy 0.932 at a true slippage of 0.25 against an assumed 0.10, rising to 0.990 when
the true model is supplied). **On real reads that effect is an order of magnitude smaller** (§4).
The simulator over-states the stutter lever because it contains no other failure — no mapping
error, no interruption, no allele the aligner mis-reads. A recommendation made from it alone was
wrong once already.

### 2.2 The instrument

[`benchmarks/lib/tract_qual_experiment.py`](../../../../benchmarks/lib/tract_qual_experiment.py)
takes a truth VCF, a caller's VCF, a confident-region BED and a **tract-ground BED**, and writes
three tables. Read its module docstring — it records why each rule is what it is.

- **`calibration.tsv`** — emitted records binned by QUAL against the share that really are at a
  variant tract. **Unit: the tract.** A record counts as truly variant when any truth record falls
  inside the tract it sits at.
- **`sweep.tsv`** — precision and recall as a QUAL threshold sweeps. **Unit: the allele**, on
  `benchmarks/giab/src/score_ng_recall.sh`'s rule (contig, position, REF and ALT equal after both
  sides are left-aligned and split), so its numbers are comparable with the standing ones. Carries
  `fp_with_no_called_copy` — false alleles listed in ALT that no sample was given a copy of, which
  is a different failure from a sample genotyped wrong.
- **`genotype.tsv`** — **the one that matters most, and §3 is about how to read it.** Where truth
  and caller both call a tract, are the caller's two alleles the truth's two. Reported twice: as
  the two tract **sequences** and as the two **repeat lengths**.

Everything is split by motif period — homopolymer against period 2 and above — because slippage
rises steeply as the period falls.

### 2.3 The tract ground, and which floors define it

The BED is built from `examples/ng_typed_region_dump <reference> <regions.bed> calling`, keeping
the `ssr_locus` rows as `chrom start end period`. **Pass `calling`, not `catalog`.** `catalog` is
the floors the catalog file was *stored* at; `calling` is what `call-from-alignments` actually
routes on, and the two differ by a lot. The dump's own doc comment still claims `catalog` is what
the run uses, which stopped being true when the routing floors changed.

Built ground lives at `tmp/tract_qual/ground/` — `tier.bed` (20,204 tracts) and `tier_sorted.bed`
(the confident regions). `benchmarks/lib/run_tract_qual_experiment.sh tandem_repeat_tier` rebuilds
both.

### 2.4 Running things

The benchmark data is gitignored and exists **only in the main checkout** at
`/Users/jose/devel/pop_var_caller/benchmarks`. From the worktree:

```
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks ./scripts/dev.sh <command>
```

mounts it read-only. `dev.sh` forwards only `CARGO_*`/`RUST*` environment, so anything else must
be set **inside** the container command: `./scripts/dev.sh bash -lc 'FOO=bar ...'`. Every ng run
needs `NG_REFERENCE_CHECK=skip`. The container has `python3`, `bcftools` and `samtools` and no
`uv`; the host is the other way round.

- ng over the tandem-repeat ground: `benchmarks/ssr_hg002/src/run_ng_coverages.sh 30x 50x`
  (51 s and 61 s; writes `results/ng/HG002_<cov>.raw.vcf` and its `.parameters.toml`).
- One parameter setting, run and scored: `benchmarks/ssr_hg002/src/sweep_tract_parameters.sh`.

---

## 3. ⚠ The genotype-coding trap, and it cost two wrong reports

**This is the most important section here.** Repeat-tract genotype accuracy was reported wrong
twice before it was reported right, and every wrong version looked plausible. If you change the
comparison, re-read this first.

### 3.1 What GIAB's truth does that ng does not

**A heterozygote carrying two different non-reference alleles is written as TWO records**, at the
same position, phased, one per haplotype:

```
chr1 4416244 . C CGTGT     ... GT=0|1
chr1 4416244 . C CGTGTGT   ... GT=1|0
```

ng writes the same genotype as **one multi-allelic record**:

```
chr1 4416244 . C CGTGT,CGTGTGT ... GT=1/2
```

These agree perfectly. **1,412 of 6,303 tracts on this benchmark are that shape.** Any comparison
that keeps one record a side, or that treats two records at a tract as two independent edits to
compose, gets them all wrong.

### 3.2 Four ways to get it wrong, each measured

1. **Compare allele strings as written, after `bcftools norm`.** Two records describing one event
   over different spans — `AGT → A` against `AGTGTGT → A,AGTGT` — share no string, because
   normalisation trims each record against *its own* ALT column, so a record with two ALTs trims
   less than one with a single ALT. **324 of 6,303 tracts read as genotype errors that were only a
   difference of spelling.** Reported accuracy 0.771 and 0.628 against the true 0.886 and 0.903.
2. **Pad both records to their union span.** Fixes that pair, and still fails wherever the two
   sides put their records at different places in the tract, because the span between them is in
   neither record's REF.
3. **Keep one record a tract a side.** Throws away 1,711 tracts of 6,303, 1,412 of them §3.1's
   shape — the class most worth measuring.
4. **Reconstruct haplotypes but refuse overlapping records.** Two truth records at the same
   position overlap as *edits* and do not as *haplotypes*: each copy takes a non-reference allele
   from exactly one of them. Refusing on spans alone scored 3,354 of 3,648 homopolymer tracts
   incomparable.

**And one more, in the plumbing:** applying the confident-region masks **before** left-alignment.
A truth record already written on an anchor base one base outside its interval is dropped, while
the query record describing the same event has not moved there yet, is kept, and then moves onto
the anchor with nothing to match. `chr1:69,233,430` is such a site — the truth carries
`TATAATAATA → T`, the Tier interval starts at 69,233,431, and ng's identical call scored a false
positive at QUAL 922. **Left-align first, mask after**, so both sides are treated alike. This is a
deliberate departure from `score_ng_recall.sh`, which masks first.

### 3.3 The rule that works

For each tract: take **all** of a side's records, read the tract's reference bases (one batched
`samtools faidx` for the whole ground), and for each haplotype copy lay on it only the records
whose allele for that copy is **non-reference**. Compare the resulting sequences as a multiset. A
side whose phase is stated (GIAB writes `0|1`) yields one pair; a side whose phase is open (ng
writes `0/1`) yields every assignment, and the two agree if any pair is shared.

**The window reaches one base either side of the tract**, because a left-aligned insertion at a
repeat's first base is anchored on the base *before* it.

**3 to 4 tracts in 100 still cannot be laid out** — a copy needing two overlapping edits at once —
and are counted apart and excluded from the denominator, because a tract the instrument cannot
compare is neither right nor wrong.

Three hand-checkable cases pin the implementation and should be kept: two phased truth records
against one multi-allelic call agree; a genuinely different call does not; one event written over
two spans agrees.

### 3.4 Sequence or repeat length — an open question that moves the headline by 3 points

Scored on the two tract **sequences**, ng is at 0.886/0.903 at 30×. Scored on the two **repeat
lengths** — which is how the STR field scores, and what ng, HipSTR and the existing caller all
emit as `REPCN` — it is at **0.915/0.913**. The two differ on 126 tracts of 6,058: both lengths
right, one spelling wrong.

**This is not a cosmetic choice.** 174 of the 268 "a truth sequence no read carried" cases are at
the *reference length* — the truth record inside the tract is a substitution, not a repeat-count
change — so scored on length, three quarters of that bucket is not a bucket. **The owner should
rule on which is the headline.** Both are in `genotype.tsv` (`genotype_accuracy` and
`length_accuracy`).

### 3.5 Two defects in existing tooling, found and not fixed

- `benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py` keeps a truth record only where its REF
  span reaches the tract's start, so an insertion anchored at `start − 1` — where left alignment
  puts every repeat-length gain — is dropped. 35 of the 268 cases have such a record, and 2,737 of
  the 19,613 clean tracts do. Widening the window takes the missing-sequence total from 434 to
  564, so the window and the overlap rule need reworking together.
- The same script's truth reconstruction misses records left-aligned before the tract; the
  slippage investigation had to work around it and reports 1,652 tracts called homozygous-reference
  that have a truth record within 30 bases.

---

## 4. Where ng stands, and what the errors are

**Genotype accuracy, comparable tracts, tandem-repeat benchmark:**

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| shipped (`--defaults`) | 0.8856 | 0.9033 | 0.8950 | 0.9132 |
| **fitted slippage + outlier weight 0.10** | **0.8907** | **0.9095** | **0.8997** | **0.9160** |
| shipped, scored on repeat length | 0.9147 | 0.9127 | 0.9294 | 0.9237 |
| both changes, scored on repeat length | **0.9185** | **0.9182** | **0.9325** | **0.9249** |

**The 648 errors at 30×, partitioned by what would have to change to fix them:**

| | homopolymer (402) | period 2+ (246) |
|---|---:|---:|
| a truth sequence was **never offered** as a candidate | **242** | **165** |
| called heterozygous, truth homozygous | 86 | 51 |
| called homozygous, truth heterozygous | 27 | 14 |
| wrong some other way, over a set that held the right alleles | 47 | 16 |

**Six errors in ten are decided before the model is consulted.** Where those missing sequences go,
over the same ground at 30× — 434 missing true sequences:

| | |
|---:|---|
| **268** | no read carried it, so the merge's allele table never held it |
| 61 | it cleared the support bar; the per-sample top-`ploidy` cut dropped it — **the only class a discovery round is aimed at** |
| 59 | the merge refused the tract, so no locus was built |
| 46 | the merge's table held it and the support bar refused it |

and the 268 split:

| | |
|---:|---|
| 121 | absent from the reads even at 300× with a median 115 spanning reads — not recoverable |
| **67** | **an alignment loss: reads carry the sequence and ng's table does not** (46 unambiguous at 30×) |
| 66 | every base the truth needs is in the reads, but no read spells the tract that way |
| 14 | the allele is longer than a 150 bp read can span with 20 bp of flank |

---

## 5. The levers: what is settable, what was tried, what it was worth

### 5.1 What can be changed today without touching code

A run's own `.parameters.toml` is written beside every output and can be edited and fed back with
`--parameters`. Two things in it matter here:

- **`repeat_tracts.slippage_by_stratum_and_group`** — per `(period, reference_repeats,
  slippage_group)`: `share_of_reads_that_slip`, `shorter_share`, `fall_off`. Empty at `--defaults`,
  which means HipSTR's shipped 0.10 / 0.50 / 0.05 everywhere.
  **Rows must be appended at the END of the file** with the empty
  `slippage_by_stratum_and_group = []` line deleted — an array-of-tables closes the table it sits
  in. Write them with `benchmarks/ssr_hg002/src/tract_slippage_rows.py`.
  **A row is needed for every repeat count a *candidate* can reach**, not just every reference
  length, because the lookup is by the candidate's count.
- **`stated_constants.repeat_tract_outlier_weight`** — 0.01, inherited from the existing caller
  and never measured. A value somebody typed must carry `warrant = "supplied"`.

`benchmarks/ssr_hg002/src/sweep_tract_parameters.sh` does all of this and scores the result. **Its
control is byte-identical to the `--defaults` run**; run it after any change to the script.

### 5.2 The sweep, on real reads, 30×

| setting | homopolymer | period 2+ | het called for a hom truth | hom called for a het truth |
|---|---:|---:|---:|---:|
| slip share 0.02 | 0.8805 | 0.9034 | 120 | 13 |
| slip share 0.05 | 0.8842 | 0.9048 | 98 | 20 |
| slip share 0.10 (shipped, supplied as rows) | 0.8851 | 0.9037 | 88 | 27 |
| slip share 0.15 | 0.8829 | 0.8976 | 82 | 42 |
| slip share 0.20 | 0.8781 | 0.8941 | 80 | 58 |
| slip share 0.30 | 0.8681 | 0.8850 | 71 | 96 |
| slip share 0.40 | 0.8384 | 0.8527 | 61 | 195 |
| shorter share 0.65 (shipped 0.50) | 0.8849 | 0.9035 | 87 | 30 |
| shorter share 0.80 | 0.8841 | 0.9031 | 89 | 31 |
| fall-off 0.15 (shipped 0.05) | 0.8857 | 0.9037 | 87 | 27 |
| fall-off 0.30 | 0.8854 | 0.9029 | 91 | 25 |
| slip rising with tract length, best of three shapes | 0.8859 | 0.9048 | 83 | 28 |
| **outlier weight 0.10 (shipped 0.01)** | **0.8892** | **0.9059** | **73** | **31** |
| outlier weight 0.20 | 0.8887 | 0.9051 | 70 | 37 |
| fitted per-stratum slippage | 0.8868 | 0.9053 | 90 | 21 |
| **fitted slippage + outlier weight 0.10** | **0.8907** | **0.9095** | **78** | **23** |

**What worked:**

- **The outlier weight, 0.01 → 0.10** — the largest single move, +0.41 points at homopolymers.
  It is the share of reads at a tract the stutter model cannot explain; `λ · U` is a floor under
  every read's emission, so the number is really **a bound on how far one read may pull a
  genotype** — the job freebayes does with a read-dependence factor and GATK with a Phred-45 cap,
  and which ng has nothing else doing. It removes 15 spurious heterozygotes of 88 and costs 4
  collapsed ones.
- **A per-stratum slippage fit**, +0.17 alone, and it fixes a *different* error (collapsed
  heterozygotes 27 → 21), which is why the two compose to +0.56.

**What did not work, and this contradicts the recommendation that started the work:**

- **Any flat change to any of the three stutter numbers.** Over a twenty-fold range of the slip
  share, accuracy moves half a point while the two error classes swing nine-fold and three-fold in
  **opposite** directions. It is a dial that trades one error for the other. The shipped 0.10 sits
  at the peak. `shorter_share` and `fall_off` do nothing at all — every setting is within 0.1
  points of the shipped one.
- **A slippage that rises with tract length**, hand-set: +0.08 at best.

### 5.3 What the six investigations found

Raw reports in [`../../reports/tract_genotype_investigation/`](../../reports/tract_genotype_investigation/),
with a README saying which claims a run corroborated.

| investigation | verdict |
|---|---|
| **slippage** — fit it from HG002's reads | Real slip share is **0.0039 at an 8-repeat homopolymer rising to 0.088 at 30**, 0.027 tract-weighted: the shipped 0.10 is **too high by 3.7-fold**, and is the value at the very top of the measured range. `shorter_share` is 0.58–0.73 against the shipped 0.50 and `fall_off` 0.16–0.45 against 0.05 — both too **low**. Produced the rows the winning run used. |
| **readmodel** — every inherited constant | The outlier weight is the only one that changes calls. The part-repeat shares can be set to anything, zero included, and change nothing. The slip-size cutoffs change at most 1 call in 8,686. The per-base substitution rate inside a tract is measurably wrong (reads give 1 mismatch in 756 at homopolymers against the model's 1 in 3,300) and worth 2–3 calls in 8,686. **Its doc comment's premise is false**: base quality inside a tract is 35.3 against 34.0 outside — better, not worse. |
| **hetfhom** — the spurious heterozygotes | 110 of 137 are a length change, **78 of those exactly one repeat**. Concentrated at long tracts: **7 in 100 at homopolymers of 25+ repeats, 1 in 100 below 13**. **The class does not shrink from 30× to 50×** (86→84, 51→51) while every other class does — an error more reads do not buy down is a model error. Every proposed filter costs more than it buys (below). |
| **prior** — the flat length spectrum | The truth is nothing like flat — **79 chromosomes in 100 sit at the reference length** against a flat shape's 11 in 100, and it is strongly stratified (0.97 at a 6–8-base homopolymer, 0.51 at 21+). **But fitting it reaches 10 of 648 errors and risks 77 correct calls**, because the dominant error is a homozygote called heterozygous and a reference-peaked prior makes that worse. **Do not fit it.** Half A answers `population_diversity.md` §4.4's open question: a per-period pooled spectrum is wrong at both ends; it would have to be per stratum. |
| **hipstr** — the outside bar | On the 2,044 period-2+ tracts both reach: ng **0.928** against HipSTR **0.925** at 30×, ng 0.939 against HipSTR **0.948** at 50×. HipSTR's own median fitted slip level is **0.04 against ng's fixed 0.05** — the median locus does not need fitting; only its top decile does. Uncontrolled: different catalogs, candidate rules, priors, `--min-reads 5`. **HipSTR's region file holds no period-1 loci**, so ng's homopolymer numbers have no comparator. HipSTR writes `.` in QUAL, so no calibration or sweep number for it is producible. |
| **unseen** — the 268 no read carried | The split in §4, three alignment losses verifiable by hand, and the observation that 174 of 268 are at the reference length. |

**Ruled out with numbers** (from `hetfhom`):

- **A stricter candidate bar** — 2 reads and 15% removes 7 homopolymer errors and destroys 10
  correct calls; 3 reads and 10% removes 11 and destroys 50; the ladder saturates at 37 of 137
  removed for **617 true alleles lost**.
- **A GQ floor** — GQ 30 withdraws 43 of 86 wrong homopolymer calls and 413 of 3,113 correct ones,
  and the bad become no-calls rather than right calls.
- **An allele-balance collapse rule** — best net anywhere on the curve is +2 tracts in 3,515.
- **Lowering the support bar to feed discovery** — 10% to 5% supplies the truth to 2 more of the
  242 never-offered tracts and hands 137 more tracts a candidate the truth does not carry.

---

## 6. What I would do next, and why

Ordered by expected value, with what I am confident about and what I am guessing at marked.

### 6.1 Settle the scoring convention first — it is cheap and it re-ranks everything else

**Confident this matters.** §3.4: sequence-scored 0.886, length-scored 0.915. If the headline is
length, then 174 of the 268 "unseen sequences" stop being errors and the candidate-selection bucket
shrinks by more than half — which changes whether §6.3 is worth doing at all. Both numbers are
already produced. **This is a decision, not work.**

### 6.2 The tract aligner: 46 tracts where the reads plainly carry the allele

**Confident this is a real defect.** Three cases with coordinates:

- `chr3:33,877,690` — an 11 bp poly-A, heterozygous. 10 of 23 reads at 30× carry `CAAAAAAAAAA`;
  ng's table holds `AAAAAAAAAAA` and a 10-base run. **The leading `C` is dropped and the tract is
  reported one repeat short.**
- `chr3:37,126,860` — homozygous non-reference. 10 of 11 reads at 30×, and 120 of 138 at 300×,
  spell `AAAAAAAAAAAAAGAAA`; ng's table holds a bare 13-base A-run. **The `GAAA` tail inside the
  tract span is discarded.**
- `chr11:37,147,255` — same mechanism, and the only 17-base allele ng holds is the one-read
  sequencing-error spelling rather than the true one carried by 12 reads of 14.

That is 2 tracts in 1,000, and the mechanism is legible: something at the tract's edge or an
interruption inside it is being discarded by the tract realigner
(`src/ng/alignment/ssr_best_path_*`). It is the only large class in §4 that is a defect rather than
a limit of the data. **Start by reproducing one of the three by hand** — pull the reads with
`samtools view`, run the tract locus generator over that interval, and see where the sequence is
lost.

### 6.3 Set the outlier weight deliberately, per period

**Confident in the direction, not the value.** 0.10 beat 0.01 on both period classes at both
depths, so it is a candidate for the shipped constant — but it was picked from a four-point sweep
on one sample. Two things say it needs care: read literally as *the share of reads nothing
explains*, the measured value is **1 in 2,300 at homopolymers and 1 in 209 at period 2+** — the
opposite ordering to what the sweep wants — and a fifth of the tracts the change touches get
worse. **So sweep it per period, and check it at 3 reads a position on the tomato cohort before
changing a shipped default**, because a floor under every emission behaves very differently when
there are three reads rather than thirty. That is the project's stated range commitment and
nothing here has tested it.

### 6.4 The fit-mode command

**Confident it is worth building; less confident it is next.** §5.2 says a per-stratum fit is worth
+0.17 points, and it is a fit for *this sample and this chemistry*, so it cannot ship as a default —
it argues for the deferred fit-mode command (spec `calling_loop_ssr.md` §3.4) that produces a
parameters file per run. The machinery to read one back works. **But +0.17 points is small**, and
HipSTR — which does exactly this per locus — is level with ng, so do not expect more from fitting
than the sweep measured.

### 6.5 Milestone E (allele discovery): I would leave it parked

**Confident.** The plan's Milestone E builds a discovery round that admits tract lengths hiding
under stutter. Its decision half is built, tested and committed (`src/ng/calling/inference/discovery.rs`,
12 tests, five mutations each predicted and each caught). **But it is aimed at 61 of 434 missing
sequences**, and the class it can only enlarge — a heterozygote called for a homozygote — is
already 137. Two further things came out of building it and should be recorded before anyone
resumes:

- **ng's retrace needs no posteriors.** The tract locus generator already realigns every read, so
  an observation's bases *are* the sequence HipSTR's alignment retrace would imply. The eligible
  set is a function of the observations and the candidate table alone, and a test asserts that a
  second round over the same evidence admits nothing. **So discovery is a pre-pass, not a round
  wrapped around the loop** — no second convergence, and no append-only emission store is needed.
  That contradicts the premise of spec §4.1's decision to "look against the converged posteriors",
  and the spec should be amended rather than the code bent to it.
- **The remaining wiring belongs inside `select_ssr`, not the loop.** A discovered allele's reads
  have to move out of the record's "no written allele explains these" column and into its `AD`,
  and that bookkeeping — the merge-index-to-candidate map and the per-sample leftover — is built
  inside selection.

### 6.6 Things I suspect and have not tested

**These are guesses. Treat them as hypotheses with a stated test, not as findings.**

- **The one-repeat spurious heterozygote at long tracts may be a per-tract problem, not a
  per-stratum one.** 78 of 137 are exactly one repeat away and they concentrate at long tracts,
  which is where a *locus's own* slippage differs most from its stratum's. The spec already has a
  per-locus slippage re-fit designed and unbuilt (`calling_em_loop.md` §5.1, three pull-back
  settings). **Test:** score those 137 tracts alone under a re-fit and see whether the class moves;
  if the per-stratum fit already got most of it, drop the idea.
- **The outlier term's uniform shape may matter more than its weight.** `λ · U` spreads the junk
  mass evenly over every reachable length; real junk is probably concentrated near the called
  allele. **Test:** replace `U` with a shape falling off with distance and re-run the same sweep.
  If the weight's optimum moves toward the measured 1-in-209, the shape was the problem.
- **`bcftools norm` may be doing damage nobody has checked.** Every number in §4 passes both sides
  through it. **Test:** score once with a comparison that never normalises (§3.3's reconstruction
  does not need it) and see whether anything moves.
- **The 66 "ambiguous" cases in §4 — every base the truth needs is in the reads, but no read spells
  the tract that way — smell like the same aligner defect as the 46.** If they are, §6.2's prize is
  112 tracts rather than 46.

### 6.7 What nothing here has touched

**One sample, one chemistry, two depths.** Every number is HG002 at 30× and 50×. Nothing says what
any of it does on the 63-accession tomato cohort at three reads a position, where the outlier
weight's floor, the candidate bar and the prior's concentration all bind differently — and that is
half of the range the caller is committed to (`CLAUDE.md`, `design_principles.md` §0).

---

## 7. State of the plan and the tree

**Plan:** [`doc/devel/ng/impl_plan/calling_loop_ssr.md`](../impl_plan/calling_loop_ssr.md).
A1–D2 ✅, E1–E3 ☐. **Checkpoint D is a hard pause the owner has not yet cleared** — the QUAL report
exists and `calling_quality_ssr.md`, which takes the decision from it, is unwritten.

**Commits on `ng-ssr-calling-loop`, six ahead of `main`:**

| | |
|---|---|
| `314cb3da` | D1 — the QUAL experiment's instrument |
| `8eef9887` | D2 — the QUAL experiment and Checkpoint D |
| `aafc7d54` | the genotypes at a tract (**its numbers were wrong; superseded**) |
| `e20abef0` | E1 — what a discovery round would admit (decision half only) |
| `2efa47da` | the genotype measurement corrected — 0.886, not 0.771 |
| `65df5c24` | two parameter values lift the tract genotypes by 0.6 points |

**Reports:**
[`ng_tract_qual_experiment_2026-09-02.md`](../../reports/ng_tract_qual_experiment_2026-09-02.md) —
the QUAL answer, and §5 the genotypes;
[`ng_tract_genotype_improvement_2026-09-02.md`](../../reports/ng_tract_genotype_improvement_2026-09-02.md) —
the sweep and what moved;
[`tract_genotype_investigation/`](../../reports/tract_genotype_investigation/) — the six raw
investigations.

**The QUAL answer, for completeness**, since it is the milestone that is formally open: ng's tract
QUAL is as well calibrated as its SNP QUAL on the same benchmark (above QUAL 200, wrong 5 times in
2,882 against ordinary sequence's 19 in 10,272) and is **not** a usable gate — at a homopolymer
tract precision peaks at 0.850 at QUAL 50 and falls back to 0.831 by QUAL 200 having shed more
than half the recall.
