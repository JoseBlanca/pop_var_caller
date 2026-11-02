# ng — a misread rate that depends on the run a read has just crossed

**DRAFT, 2026-11-01. A prototype with a hard-coded table is built (uncommitted, worktree
`pop_var_caller-context-error`) and measured in §6.1; the owner kept it and set the next two steps
(§6.2).**

Reads on: [`read_likelihoods.md`](read_likelihoods.md) §4.3 (the STR substitution term),
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1 and §4.1 (what the pre-pass fits and
counts), [`calling_quality.md`](calling_quality.md) §8.2 (the strand test that was withdrawn), and
[`../../specs/design_principles.md`](../../specs/design_principles.md) principle 0.

## 1. The problem

**A read that has just been sequenced through a long run of one base misreads the bases after it,
and most often reads the next base as one more copy of the run.** The STR path cannot see this: it
scores every base of a read against a candidate with one substitution rate, fitted per read group
and stratum (`read_likelihoods.md` §4.3), so a base the sequencer read right after sixteen `A`s costs
the same as a base read anywhere else.

**What it costs.** On HG002 chromosome 20 at 30×, 264 of pvc's 1,310 errors lie inside a tract. Of
30 of those traced read by read (`tmp/intract/` on branch `flank-indel-junction`), 9 are false
substitutions inside the tract, and in each the reads carrying the false allele come from one
strand: the strand whose reads were sequenced through a long run before reaching the misread base.
Two examples, from the `--explain-loci` rows of that tracing:

- **chr20:44,076,996**, reference `A`×16 `G` `A`×10. Three reads show the `G` as an `A`; all three
  are forward reads, which reach the `G` after reading the sixteen `A`s. Called `0/1`, GQ 31; GIAB
  has no variant.
- **chr20:36,176,622**, reference `T`×8 `C` `T`×16. Four reads show the `C` as a `T`; all four are
  reverse reads, which reach the `C` after reading the sixteen `T`s. Called `0/1`, GQ 46; GIAB has
  no variant.

Four more of the 30 are the same misread one step outside the tract: a flank base read as the run's
base, which the mapper writes as a deletion of that flank base (`str_flank_ownership.md` §1.3).

**A strand test does not catch these** (`calling_quality.md` §8.2). True calls at tracts are
strand-skewed about as often as false ones, and most of these false alleles rest on three reads,
where the strand test charges nothing by design. The evidence against them is where the misread
base sits, not the strand split.

## 2. What the reads show

**Measured on HG002, chromosome 20 from 10 to 20 Mb, 30×** (`tmp/homopolymer/run.sh`; reads of
mapping quality 20 or more, every base within 10 of a GIAB variant skipped). For each run of one base
in the reference, the base just past each end of it is looked at, in two groups of reads: those
sequenced **through** the run before reaching that base, and those that reach it **before** the run.
The second group is the control: a true variant shows in both, a misread caused by the run only in
the first.

**The next base, read as one more copy of the run** (reads in 10,000):

| run length | reads through the run | reads reaching the base first |
|---|---|---|
| 4 | 27 | 18 |
| 8 | 51 | 41 |
| 10 | 139 | 71 |
| 12 | 287 | 94 |
| 14 | 474 | 32 |
| 16 | 971 | 69 |
| 18 | 1,436 | 156 |
| 20 | 1,582 | 60 |
| 23 | 2,491 | 30 |

So below about 8 bases the run barely matters, and from 12 bases on it dominates: after a 16-base
run, about 1 read in 10 sequenced through it reads the next base as the run's base, against fewer
than 1 in 100 from the other side. Misreads to a base other than the run's are rare at every length
(after a 16-base run, 30 reads against 285).

**The bases after the next one are affected too.** Reads sequenced through a run of 16 bases or
more mismatch the reference about 6 times in 100 at every base for at least the next 20
(1,633 in 10,000 at the next base, 500 to 710 at bases 2 to 20), against about 3 in 1,000 for reads
reaching the same bases first. After runs of 12 to 15 bases it is about 2 in 100 against 2 in 1,000;
after runs of 8 to 11, about 4 in 1,000 against 2 in 1,000.

**The sequencer already marks these bases.** The misread bases carry a mean base quality of about 10
— an error probability of 1 in 10, close to the rate measured — against 35 for bases read as the
reference. The generic path, which scores reads by their own qualities (`read_likelihoods.md` §3.2),
already discounts them; the STR path does not, because it uses the fitted rate instead.

The same pattern appears after a short tandem repeat: at chr20:64,021,935, reverse reads crossing
`GT`×7 read the next `C` as a `T`, continuing the repeat. This draft measures homopolymers only; §7
leaves the other periods open.

### 2.1 After repeats of period 2 to 4 (measured 2026-11-01, night)

The same measurement, over repeats of period 2 to 4 (`tmp/homopolymer/measure_str.py`, same reads
and region): the base right after a repeat, read as **one more base of the repeat** — the base
the repeat would put there next — in reads sequenced through the repeat against reads reaching the
base first. Reads in 10,000:

| period | repeat length | through the repeat | reaching the base first | mean base quality of those misreads |
|---|---|---|---|---|
| 2 | 12–15 | 68 | 58 | 24 |
| 2 | 16–19 | 220 | 171 | 29 |
| 2 | 20–29 | 303 | 179 | 27 |
| 2 | 30+ | 420 | 187 | 17 |
| 3 | 12–15 | 55 | 45 | 28 |
| 3 | 20–29 | 177 | 187 | 25 |
| 3 | 30+ | 647 | 739 | 27 |
| 4 | 12–15 | 26 | 20 | 12 |
| 4 | 20–29 | 104 | 66 | 18 |
| 4 | 30+ | 282 | 243 | 25 |

**Unlike after a homopolymer, the two directions are close and the misread bases carry good
quality.** After a 16-base homopolymer, reads sequenced through it misread the next base 14 times
as often as reads from the other side, at a base quality of about 10; after these repeats the ratio
is between 0.9 and 1.3 except for period-2 repeats of 20 bases or more (1.7 and 2.2), and the base
quality is 17 to 29. A read showing the repeat continuing from both sides with a good base quality
is most likely a read one unit longer that the mapper wrote as mismatches — stutter, which Model A
already prices — not a misread caused by the run. The `GT`×7 case of §1 (chr20:64,021,935, 14
bases) sits in the 12–15 row, where the two directions differ by 10 reads in 10,000.

## 3. The decision

**Give the STR path's substitution term a per-base error rate that depends on the run the read has
just crossed, in the read's own sequencing direction, fitted from the data per read group.**

1. **The context of a base is the run it follows in sequencing direction**: the length of the
   longest run of one base that ends within the last `d` bases before it, read in the direction the
   read was sequenced — left to right for a forward read, right to left for a reverse one. It is
   computed on the **candidate** sequence, with the flank bases beside it, not on the read: the
   question is how likely this read is to have misread the candidate.
2. **The rate at a base is the read group's fitted rate for that context**: a small table, by run
   length class (under 8, 8–11, 12–15, 16 or more) and distance from the run (the next base, then
   2–20), with two entries per cell — the chance of reading the run's base, and the chance of
   reading any other wrong base. Contexts under 8 use today's fitted rate unchanged.
3. **An observation is scored once per direction.** Every stored observation carries how many of
   its reads were forward (`SequenceObservation::num_fwd`); the row adds its forward reads' term and
   its reverse reads' term, each computed with its own direction's rates. The sequence is the same,
   so this is at most two evaluations per distinct sequence, not one per read.
4. **The table is fitted in the parameter pre-pass from reads against the reference**, with the
   control of §2: at each context cell, the rate in reads sequenced through the run minus the rate
   in reads reaching the base first. A true variant raises both and cancels; the misread raises only
   the first. No truth set and no called genotype is needed.

**Why a fitted table and not the reads' own base qualities.** The qualities carry the signal
(§2), but the STR evidence does not keep them per base: an observation is one sequence and a count
of reads, with a per-read error sum (`read_likelihoods.md` §4.3), so two reads of one sequence with
different qualities are one row. Keeping per-base qualities would split every observation by its
reads' qualities, and the psp file and the merge would grow with depth. The context is a property of
the candidate sequence and the read's direction, both already stored, so the table costs no storage.
The qualities are used to check the table: a fitted rate far from what the misread bases' qualities
say would point to a fault in the fit.

## 4. One sample and thousands, three reads and three hundred

- **One sample, three reads a position.** The table is pooled over the whole genome, not fitted per
  locus. The measurement of §2 saw 2,935 reads through 16-base runs in 10 Mb at 30×; a whole
  genome at 3× gives about thirty times that, enough for every cell. Applied at a locus, it
  changes how much three one-strand reads count against the reference: at chr20:44,076,996 each of
  the three forward reads showing `A` after sixteen `A`s costs the reference allele about
  `ln(1/10)` instead of about `ln(ε/3)` at the fitted rate, a difference of about 3 to 6 nats a read
  for fitted rates of 20 to 1 in 1,000. That is the case where the strand test can do nothing.
- **Thousands of samples.** The table is per read group and counts only; its cost grows with reads,
  as the pre-pass already does. The scoring cost is at most twice today's per distinct sequence.
- **Three hundred reads.** Unchanged in kind: two evaluations per distinct sequence. A true allele
  supported by both strands loses nothing, since only the through-strand reads are discounted, and
  only at bases after a long run.
- **Where it can hurt.** A true variant just past a long run whose reads come mostly from one
  direction — a heterozygous base right after a 20-base run at the end of a capture target — is
  discounted on those reads. The reads from the other direction still count at full weight, so
  it costs depth, not the call, wherever both directions are present.

## 5. What it touches

- **The STR emission** (`src/calling/likelihood/ssr_emission.rs`): the substitution term takes a
  per-base rate instead of one rate, and a row is summed per direction. Bundles use the same per-base
  rate in place of their per-run rates where a context applies.
- **Partial observations** keep no strand count today (`PartialObservation`; `calling_quality.md`
  §8.2, item 4). Until they do, a partial read is scored with the average of the two directions'
  rates, weighted by the locus's complete reads' forward share. Storing the count is a psp change
  and is left to the measurement (§6).
- **The parameter pre-pass** counts, per read group and context cell, reads through and reads
  reaching first, and mismatches in each. It counts on the generic path's reference positions, so
  the STR loci it is applied to are not the ones it is fitted on.
- **The parameters file** gains the table, per read group.
- **The generic path does not change**: it already scores by base quality.

## 6. Measuring it

A prototype first, with the table fitted once from the HG002 chr20 counts of §2 and hard-coded,
before any pre-pass work:

1. HG002 chromosome 20 and the tandem-repeat benchmark, at 5× and 30× (`tmp/hg002/run2.sh`), against
   GIAB with vcfeval: the bar is fewer false calls inside tracts with no rise in missed ones. Report
   how many of the 9 traced one-strand cases go.
2. The stuttering benchmark and the flank simulation, where the simulator plants no context error:
   nothing should change beyond noise. A change there means the term is acting where it should not.
3. Tomato (63 samples, about 3×): records that change, by run length beside them.

Then the pre-pass fit, checked against the hard-coded table on HG002 and against the base qualities
as §3 says.

### 6.1 What the prototype measured (2026-11-01)

The table hard-coded from §2's counts, homopolymers only. Rule: the fitted rate plus the excess
(reads through the run minus reads reaching the base first), not the larger of the two. Partial reads
carry no strand count and are scored at the flat rate, as today.

HG002 against GIAB, vcfeval, genotype match, before → after (`tmp/hg002/run_ctx.sh`):

| run | missed | false calls |
|---|---|---|
| chr20, 30× | 748 → 741 | 549 → 529 |
| chr20, 5× | 24,611 → 24,611 | 4,617 → 4,617 |
| tandem-repeat benchmark, 30× | 14,212 → 14,179 | 2,393 → 2,359 |
| tandem-repeat benchmark, 5× | 24,042 → 24,044 | 2,544 → 2,547 |

On chr20 30× every change is inside a tract: false calls 138 → 118, missed 116 → 109. The flank
simulation's right, no-call and wrong counts are unchanged. Of the 9 traced one-strand cases, 4 go
(misreads right after a homopolymer), 1 drops from quality 35 to 5 and stays, and 4 are untouched:
three after a `GT`, `CA` or `GATG` repeat, one inside a `G` run. At 5× nothing moves: the misreads
that fool the caller come about three reads at a time, and at 5× a site rarely has the reads for the
term to decide anything. A read matching the reference still pays a small charge after a long run
(about 0.65 nats over the ten bases after the run at chr20:44,076,996).

### 6.2 The owner's decisions (2026-11-01)

1. **Keep the term and build the fitted version**: the table estimated in the parameter pre-pass
   (§3, point 4) instead of hard-coded.
2. **Before the fitted version, extend the table to repeats of period 2 to 4** (§7, item 1), which
   covers three of the four untouched cases. Measured as §6 says before it is kept.
3. **2026-11-02, after §2.1 and §6.3: do not extend the hard-coded table by hand; the fitted table
   counts repeats of period 2 to 4 the same way as runs of one base, and each read group's data
   decide what they charge.** The fitted version is specified in
   [`str_context_error_fit.md`](str_context_error_fit.md).

### 6.3 The prototype on main, re-measured (2026-11-02)

The same prototype, ported unchanged onto main 92b00dd8 (branch `str-context-error`), against main's
own build. HG002 against GIAB, vcfeval, genotype match, before → after (`tmp/hg002/run.sh`,
`summary.tsv`). "Fitted" is generate-psps → estimate-parameters → call-from-psps; the two builds fit
identical parameters, since the fit does not score STR reads.

| run | missed | false calls |
|---|---|---|
| chr20, 30×, `--defaults` | 753 → 746 | 473 → 453 |
| chr20, 30×, fitted | 841 → 837 | 373 → 357 |
| chr20, 5×, `--defaults` | 24,519 → 24,519 | 4,619 → 4,619 |
| chr20, 5×, fitted | 24,913 → 24,908 | 4,346 → 4,341 |
| tandem-repeat benchmark, 30×, `--defaults` | 14,245 → 14,212 | 2,185 → 2,151 |
| tandem-repeat benchmark, 30×, fitted | 14,512 → 14,478 | 1,919 → 1,890 |
| tandem-repeat benchmark, 5×, `--defaults` | 24,028 → 24,030 | 2,537 → 2,540 |
| tandem-repeat benchmark, 5×, fitted | 24,300 → 24,306 | 2,389 → 2,392 |

On chr20 every change is inside an STR segment: at 30× false calls there go 139 → 119 (`--defaults`)
and 118 → 102 (fitted), missed 115 → 108 and 108 → 104.

**The 5× losses on the tandem-repeat benchmark are the cost §4 names.** Each is a heterozygous
variant inside a long run of one base seen by five to eight reads, of which one carries the reference
base and has crossed ten or more bases of the run to reach the variant. At chr10:15,111,520 (reference `A`×22, true allele `A`×10 `G` `A`×11; one reference read, four `G`
reads) the term makes that one read's `A` about 2.4 nats cheaper to explain as a misread of the `G`
allele, main's `0/1` led by 0.7, and the call becomes `1/1`. chr18:40,078,132 and
chr1:92,216,311 are the same.

**The simulations do not move.** The flank simulation's scores are identical in all four cases (one
and 20 samples, 3× and 30×). The stuttering simulation is identical at one sample, 30×, in every
route; in its 50-sample 3× cohort, of about 150,000 sample-segment cells, false calls go
7,388 → 7,395 with `--defaults` and 387 → 386 fitted.

## 7. Open

1. **Short tandem repeats beyond homopolymers** (§2's `GT` case): settled by §6.2 item 3 — the
   fitted table counts periods 2 to 4 and the data decide.
2. **The distance cut-off of 20 bases** is where the measurement stopped, not where the effect does.
3. **Whether the table should also enter the delimiter's measurement** of where a tract ends, which
   is where the four one-strand flank misreads of §1 arise. That belongs with the flank redesign
   (`str_flank_ownership.md`), which measures each read's window from its sequence.
