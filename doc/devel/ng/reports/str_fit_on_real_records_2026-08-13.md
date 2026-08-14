# The repeat-tract half of the joint fit, built and run on real reads

*Research report, 2026-08-13. Written for a reader who has read none of the specifications.*

*What is new: `src/ng/parameter_estimation/joint/ssr_fit.rs` — the estimator, lifted out of
`examples/ng_joint_str_harness.rs` rather than written again; three new modes in that harness
(`library`, `borrow`, `borrow-edge`); and the repeat-tract half of
`examples/ng_joint_records_walk.rs`, which drives it from records built off real alignments. Raw
output under `tmp/records/` and `tmp/borrow_crossover_*.log`.*

---

## 1. What this is

While a polymerase copies a repeat tract it sometimes **slips**, adding or dropping a whole repeat
unit, so a read reports a tract one unit longer or shorter than the DNA it came from. Before calling
any variant, ng estimates how often that happens — once, over the whole cohort. Tracts are grouped
into **strata**, every tract sharing a motif length and a reference repeat count, because slippage
depends on repeat count more than on anything else.

Per stratum the fit produces five numbers:

| number | what it is |
|---|---|
| slippage level | how often a read reports a length other than its allele's |
| shorter-share | of the reads that slip, the share showing a **shorter** tract |
| fall-off | how fast two-unit slips fall off against one-unit slips |
| concentration | how monomorphic the stratum's tracts are — small means most tracts carry one length |
| length spectrum | how the stratum's chromosomes are spread over the tract lengths |

The first three are per (read group × stratum); the last two are per stratum. **None of them is per
tract**: a tract's own length frequencies are a latent vector drawn from the stratum's Dirichlet and
integrated away on a fixed 256-point rule, which is what
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.1 and §4.2 settle.

**Until today this half did not exist as code.** The ordinary-position half has been running on real
records since 2026-08-13 morning
([`joint_fit_against_truth_2026-08-13.md`](joint_fit_against_truth_2026-08-13.md)); the repeat-tract
half existed only inside a research harness that fed itself drawn data.

### What the real records gave

**On tomato it recovered the rise in slippage with repeat count, which is the thing the whole
stratification exists for.** A read misreads a homopolymer tract of 8 repeats **2 times in 1,000** and
one of 12 repeats **9.4 times in 1,000** — 4.7 times as often across four repeat counts, rising at
every step, about **1.47-fold a repeat count** against a specification that predicted roughly 1.3-fold
from a different measurement. Five separate fits, no constraint tying them together (§6.1).

**It also answers a question the specification records as unanswerable until now**: how monomorphic a
real stratum's tracts are. Tomato homopolymers come back at 0.52 to 1.56 and its dinucleotides at six
repeats at 5.25 (§6.2).

**Two findings are uncomfortable and both are about coverage, not about the model.** Half the reads
that reach a tomato repeat tract never cross it, so they report no length and the fit drops them
(1,538,186 against 1,587,703); and 65 of tomato's 71 strata hold too few tracts to say anything on
their own. On the GIAB trio, whose region set holds 216 tracts in 32 strata, that second problem is
total: all fifteen homopolymer strata borrowed from each other and came back with **one identical
answer**, the repeat-count axis flattened to a constant with nothing in the output to say so (§5.3).

---

## 2. What was built

`ssr_fit.rs` is the harness's estimator moved into the library. **It was lifted rather than
rewritten**, because two implementations of one model are two things to keep agreeing. Three things
are generalised in the move, and each is the specification catching up with the harness rather than a
new idea:

1. **Alleles reach further than the read buckets.** The records store offsets `−4 … +4` with the ends
   saturating; the lengths the fit may place allele mass on reach `±6`
   ([`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §3.2), which is
   what lets an end bucket be attributed to a distant allele rather than to a far slip. The end bucket
   is scored by its marginal — *at least four repeats short* gets the sum over every offset it
   absorbs, including a slip that started outside the recorded range. With the two spans set equal the
   arithmetic reduces to the harness's exactly, which is what makes §3's comparison possible.
2. **The homozygote excess is per sample**, arriving from the ordinary-position half, where the
   harness supplied one number for the whole panel. It enters as two dot products a sample rather than
   a rebuilt genotype prior, so per-sample inbreeding costs nothing.
3. **Slippage is per read group**, as spec §4 has it. Read groups are named into **slippage groups**
   so that a run may pool them; one group per read group is the specified default.

Everything else is the harness's code: the same 256-point Halton rule pushed through stick-breaking
Beta quantiles, the same read-likelihood cache held across parameter moves, the same coordinate ascent
from three starting points.

**What is not in it**, and each is stated in the module's own documentation rather than left to be
discovered:

- **Reads that reached a tract without crossing it are counted and then dropped.** They report no
  length, so they carry nothing about slippage; but the censoring is not random — a tract longer than a
  read is never crossed, in every sample at every depth. §5 shows this is not a small number on real
  data.
- **The guard's reads are counted and dropped**, and a tract over the guard's threshold — more than
  one read in ten of those differing from the reference length differing by a *non-whole* number of
  repeat units — is left out of the fit entirely, as the records specification asks.
- **The mismatch list is unread.** The substitution rate inside a tract is a separate parameter this
  work does not fit.

---

## 3. Does the library agree with the harness?

**Yes, to every digit.** The harness's new `library` mode draws one stratum, fits it with the
harness's own estimator, converts the identical draw into the library's types and fits it again.

| number | harness | library | relative difference | drawn truth |
|---|---:|---:|---:|---:|
| slippage level | 0.080624 | 0.080624 | 0 | 0.0800 |
| shorter-share | 0.827487 | 0.827487 | 0 | 0.830 |
| fall-off | 0.256167 | 0.256167 | 0 | 0.250 |
| concentration | 0.492193 | 0.492193 | 0 | 0.500 |
| spectrum, three classes | 0.271220 / 0.472834 / 0.255947 | identical | 0 | 0.262 / 0.476 / 0.262 |

*20 samples, 6 reads a tract, 4,000 tracts, three length classes.* Repeated at **three reads a tract
and five length classes** — tomato's depth, and enough classes to exercise four stick-breaking
dimensions — the two agree to every digit again (slippage level 0.081126 both ways, concentration
0.571768 both ways). The agreement survives the two savings §7 describes.

```sh
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness library 20 6 4000 3 0.5
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness library 20 3 4000 5 0.5
```

**The positive control is in the library's own tests**, not only in the harness: a drawn stratum of
1,500 tracts comes back with the slippage level within 10% of the truth, the shorter-share within
0.05 and the concentration within 25%. Without it a clean-looking answer on real data could not be
told from an estimator with no power in it — which is the failure this project has recorded before.

---

## 4. What a thin stratum should do — the one design decision here

68 of tomato's 141 strata hold fewer than a hundred tracts, which is far below any count that can
carry an answer: at 50 tracts and three reads a site the twelve draws of the earlier sweep put the
fall-off anywhere between 0.139 and 0.338 against a truth of 0.250
([`str_stratum_size_sweep_2026-08-13.md`](str_stratum_size_sweep_2026-08-13.md)). The design says such
a stratum **borrows from its neighbouring repeat counts**. Where borrowing has to start had never been
measured and nobody had built it.

### 4.1 Recommendation

**Borrow below 1,000 tracts; take both sides of a repeat count together; refuse below 50 even after
borrowing.** In the code these are `DEFAULT_BORROWING_FLOOR = 1_000` and
`DEFAULT_REFUSAL_FLOOR = 50`, both named constants with the measurement in their documentation.

### 4.2 Why borrowing is nearly free, and it is the symmetry that makes it so

Borrowing trades one error for another. A thin stratum's own answer is centred on its own truth and
moves a long way from draw to draw; a borrowed answer is steady and sits on its **neighbours'** value
instead of its own, because slippage genuinely rises with repeat count — roughly **1.3-fold a count**
([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3).

So the neighbours in this measurement are drawn at 1.3 times and 1/1.3 times the thin stratum's own
slippage level, and the statistic is the **root-mean-square error against the truth** over the draws —
one number that adds the two failures together, where a mean error alone would call the noisy arm
perfect and a spread alone would call the displaced arm perfect.

*Its two neighbours bring 600 tracts each; 20 samples, three reads a tract, three length classes,
concentration 0.5, eight draws at each size. Each cell is the mean error, the spread between draws,
and the root-mean-square error, all as a percentage of the truth.*

| at 50 tracts | fitted on its own | with both neighbours pooled in |
|---|---|---|
| slippage level | +1.0% ± 10.2% (**rmse 9.6%**) | +1.2% ± 1.5% (**rmse 1.8%**) |
| shorter-share | −1.4% ± 5.3% (**5.1%**) | +0.2% ± 0.7% (**0.7%**) |
| fall-off | −1.3% ± 37.6% (**35.2%**) | +1.7% ± 4.2% (**4.2%**) |
| concentration | −2.4% ± 17.7% (**16.7%**) | +1.4% ± 3.2% (**3.4%**) |
| commonest length | −3.5% ± 12.9% (**12.6%**) | +0.1% ± 1.6% (**1.5%**) |

| at 250 tracts | fitted on its own | with both neighbours pooled in |
|---|---|---|
| slippage level | +0.1% ± 2.8% (**rmse 2.7%**) | +1.0% ± 1.4% (**rmse 1.6%**) |
| shorter-share | −0.3% ± 2.9% (**2.8%**) | +0.2% ± 0.7% (**0.7%**) |
| fall-off | −3.8% ± 8.8% (**9.1%**) | +0.8% ± 3.6% (**3.4%**) |
| concentration | −3.7% ± 8.3% (**8.6%**) | +0.7% ± 2.7% (**2.6%**) |
| commonest length | +0.0% ± 4.2% (**3.9%**) | +0.1% ± 1.6% (**1.5%**) |

| at 1,000 tracts (five draws) | fitted on its own | with both neighbours pooled in |
|---|---|---|
| slippage level | +0.3% ± 1.7% (**rmse 1.5%**) | +0.7% ± 0.7% (**rmse 0.9%**) |
| shorter-share | +0.3% ± 1.3% (**1.2%**) | +0.3% ± 0.7% (**0.7%**) |
| fall-off | −0.7% ± 12.8% (**11.5%**) | +0.7% ± 7.7% (**7.0%**) |
| concentration | −1.2% ± 2.9% (**2.9%**) | +1.0% ± 0.7% (**1.2%**) |
| commonest length | −0.8% ± 1.5% (**1.5%**) | −0.1% ± 1.7% (**1.5%**) |

```sh
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness borrow 20 3 50 3 0.5
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness borrow 20 3 250 3 0.5
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness borrow 20 3 1000 3 0.5
```

**Borrowing wins at every size measured, and the margin closes as tracts accumulate**: on the
fall-off it is eight times better at 50 tracts, under three times at 250, and 1.6 times at 1,000. On
the commonest-length share the two are equal at 1,000 — 1.5% either way. That is the shape a
crossover has, and **the two have not crossed anywhere in the measured range.** The floor is
therefore not the place where borrowing stops paying in error; §4.3 says what it is instead.

**And the displacement everyone expects — the neighbours are 30% away in slippage level — does not
appear: the borrowed level comes back +1.2%, not +30%.** That is not luck. The two neighbours sit
either side, so their displacements point opposite ways and the tracts-weighted mean of 0.0615,
0.0800 and 0.1040 lands at 0.0827, three percent above the middle stratum's own truth.

**That is why the library takes both sides of a repeat count together rather than one neighbour and
then a test.** A rule that stopped as soon as the floor was reached would keep whichever neighbour it
happened to reach first and carry that neighbour's whole 30% displacement. Taking the ring is a
one-line difference in the borrowing loop and it is the difference between +1.2% and −23%.

*A cross-check that the library behaves like the program it came from:* its own-tracts-only arm at 50
tracts gives spreads of 10.2%, 5.3%, 37.6% and 17.7% on the level, shorter-share, fall-off and
concentration, against the harness's 9.0%, 5.0%, 31.3% and 17.1% at the same tract count and depth.

### 4.3 Why 1,000, when borrowing never actually lost

**The measurement does not name the floor, and it is worth saying so plainly.** Borrowing was better
or equal on every number at every tract count run — 50, 250 and 1,000. Taken alone the tables argue
for borrowing everywhere, and that answer is wrong for a reason no single stratum's error can show.

**What borrowing everywhere costs is the axis.** Each stratum borrows its own neighbours, so each
answer is centred on its own repeat count, which is what keeps the error small. But two neighbouring
strata that pool the same tracts get the same answer, and a run in which every stratum borrows is a
run in which the fitted slippage barely changes from one repeat count to the next. Slippage genuinely
rises about 1.3-fold a repeat count and reaches a twenty-two-fold spread across the range
([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4, §4.3). **Smoothing that away
would remove the reason repeat count is a stratification axis at all** — and §5.3 is what it looks
like when it happens on real data.

So the floor is where a stratum's own answer becomes good enough to be worth keeping, not where it
becomes better than a borrowed one:

- **At 1,000 tracts the own-tracts answer is already close.** Its errors are 1.5%, 1.2%, 11.5%, 2.9%
  and 1.5% against the borrowed arm's 0.9%, 0.7%, 7.0%, 1.2% and 1.5% — a factor of 1.6 on the worst
  of them and equal on one. Keeping the stratum's own answer costs about a percentage point and buys
  a repeat-count axis that varies.
- **At 250 it is not close**: 2.7%, 2.8%, 9.1%, 8.6% and 3.9% against 1.6%, 0.7%, 3.4%, 2.6% and
  1.5% — three times worse on the fall-off and the concentration.
- **A floor of 5,000 would make every stratum a borrower.** The per-stratum selection cap is itself
  5,000 tracts ([`str_stratum_size_sweep_2026-08-13.md`](str_stratum_size_sweep_2026-08-13.md)), so
  no stratum ever arrives at the fit holding more than that.

### 4.4 What this recommendation rests on, and where it does not hold

- **The symmetric ring is the whole of the argument, and a stratum at the end of the repeat-count
  range does not get one.** The longest tracts have no longer neighbour, so their borrowing is
  one-sided and the displacement no longer cancels: the earlier work prices that at **15 to 25% of the
  level per repeat count borrowed** ([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md)
  §4.3). Set that against the own-tracts arm at 50 tracts — 9.6% on the level, 35.2% on the fall-off —
  and **a one-sided borrow costs the level about twice what keeping it would, while still saving the
  fall-off eight times over.** Which way that trade should go is not measured; `borrow-edge` mode
  exists in the harness to measure it and was not run, because the machine was carrying four other
  jobs. Until it is, the recommendation above is a recommendation about strata with neighbours on both
  sides.
- **One stratum shape, one depth, one panel size.** Three length classes, concentration 0.5, 20
  samples, three reads a tract. A stratum whose tracts are more nearly monomorphic carries less signal
  per tract for the concentration.
- **The neighbours here are the same shape as the receiving stratum apart from their slippage level.**
  A real neighbouring repeat count also has a different length spectrum and a different concentration,
  and pooling those was not measured.

---

## 5. The GIAB trio — the first repeat-tract numbers this route has produced from real reads

Three deeply sequenced human samples over the 100 regions all three benchmark sets share: 452,288
bases, ~300 reads a position. The records were rebuilt by the same walk that fits them
(`tmp/run_records_trio.sh`), because the records themselves are held in memory and never written to
disk.

**The region set holds 216 repeat tracts in 32 strata.** That is the first thing the run says and it
governs everything after it: 216 tracts is a fortieth of the 5,000 the size sweep calls the floor at
three reads a site, and a fifth of the 1,000 it calls the floor at six.

### 5.1 What the reads did at those tracts

| | |
|---|---:|
| reads that crossed a whole tract, over the three samples | 171,930 |
| reads that reached a tract and crossed no whole copy of it | **102,308** |
| reads whose tract differed by a non-whole number of repeat units | 9 |
| tracts over the guard's threshold, and so left out of the fit | 6 of 216 |

**More than a third of the reads that reach a repeat tract here never cross it** — 102,308 against
171,930, so 37 reads in 100 arrive and report no length at all. Nothing about slippage can be
estimated from them and the fit drops them, which is correct; what is not obviously correct is that
**the censoring runs along the very axis the parameters are fitted within.** A tract longer than a
read is never crossed in any sample at any depth, so the reads that go missing are concentrated in the
long-repeat strata — exactly the strata that are already thinnest. This number had not been measured
before; the records specification names the mechanism (§1.1's four states) and prices nothing.

*What this measurement cannot say:* which tract lengths the uncrossed reads sat at. The count is
accumulated per stratum in the code and printed only as a cohort total, so the shape of the loss along
repeat count is one print statement away and was not taken.

### 5.2 What was fitted

`SSR_ALLELE_SPAN=6` — the specification's span — borrowing below 1,000 tracts, refusing below 50. The
whole fit took 98.8 s after a 467 s run.

| motif length | repeat counts | tracts with reads | fitted |
|---|---|---:|---|
| 1 (homopolymer) | 8 – 23, fifteen strata | 187 | slippage level **0.0120**, shorter-share **0.723**, fall-off **0.131**, concentration **0.666** |
| 2 | 6 – 38, eight strata | 17 | refused — below the floor of 50 even with every period-2 stratum pooled |
| 3 | 6 – 7 | 1 | refused |
| 4 | 6 – 13, six strata | 5 | refused |
| 5 | 5 | 0 | refused — no spanning reads |

**The one fitted number is plausible and it is the first of its kind from this route.** A slippage
level of 0.0120 means 12 reads in 1,000 report a homopolymer length other than their allele's, which
sits inside the range the per-sample route measured from tomato — 9 reads in 10,000 below four repeats
against 2 in 100 at six repeats and above
([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §5) — for tracts that are all eight
repeats and longer. The shorter-share of 0.723 says a slipped read shows a **shorter** tract 2.6 times
as often as a longer one, the same direction as tomato's dinucleotides and weaker than their 4.9-fold.

### 5.3 And here is what is wrong with it

**All fifteen homopolymer strata carry the identical answer.** A homopolymer of 8 repeats and one of
23 are given the same slippage level, the same direction split, the same fall-off and the same
concentration — 0.0120, 0.723, 0.131, 0.666 down every row of the table.

That is borrowing doing exactly what it was told. Every one of the fifteen is far below the floor, so
every one reaches for its neighbours; each of them ends up pooling all fifteen; the pooled sets are
identical, so the library fits them once and hands the same answer to all fifteen. **The stratification
has been flattened to nothing**, and with it the rise along repeat count that is the reason repeat
count is a stratification axis at all — the same rise §6.1 measures on tomato at 1.47-fold a count.

**This is not a defect in the borrowing rule; it is the rule's cost made visible, and the region set
is what makes it total.** On 452 kb of genome there are 187 homopolymer tracts. On the whole human
genome there would be hundreds of thousands, most strata would clear 1,000 on their own, and borrowing
would touch only the tails. **What this run proves is that the answer a run gets depends on how much
genome it walked, in a way nothing in the specifications currently says.** A user who walks a small
region set and reads the emitted per-stratum table will see a repeat-count axis that does not vary,
with nothing in the output to say the axis was averaged away.

**What the emitted parameters must therefore carry, and today do not:** how many tracts stood behind
each stratum's own answer, and which strata were pooled to produce it. The library records both
([`StratumFit::borrowed`](../../../src/ng/parameter_estimation/joint/ssr_fit.rs) and `tracts_fitted`)
and the walk prints them; nothing in the specification requires them of the emitted parameters. §8
proposes that.

---

## 6. Sixty-three tomato accessions — the run that says the estimator works

63 accessions at 2.4 to 30.6 reads a position, over the bench region set's 80 spans and 8 Mb: 1.99
million ordinary positions and **4,164 repeat tracts in 71 strata**. The ordinary-position half ran
first and returned exactly what it returned this morning — a read misreading at 0.00333 at an ordinary
position and 0.0239 at a mismapped one, 1 position in 30 mismapped, the population's expected
heterozygosity 4.886 per kilobase — so the repeat-tract half sits on an unchanged base. The whole run
took 4,539 s, of which the repeat-tract fit was 1,690 s.

**This arm was run with borrowing switched off** (`SSR_BORROWING_FLOOR=0`), so every stratum speaks
from its own tracts alone and nothing is smoothed. Six of the 71 clear the refusal floor of 50 tracts.
Those six hold 3,661 of the 4,164 tracts, so **8% of the strata carry 88% of the tracts.**

### 6.1 Slippage rises with repeat count, and the rise is measured here for the first time from real reads

| motif | repeats | tracts | reads crossing | slippage level | shorter-share | fall-off | concentration |
|---|---:|---:|---:|---:|---:|---:|---:|
| homopolymer | 8 | 2,082 | 937,557 | **0.0020** | 0.595 | 0.636 | 0.61 |
| homopolymer | 9 | 887 | 363,500 | **0.0027** | 0.637 | 0.631 | 0.88 |
| homopolymer | 10 | 350 | 121,244 | **0.0037** | 0.623 | 0.629 | 1.56 |
| homopolymer | 11 | 153 | 34,904 | **0.0059** | 0.713 | 0.540 | 0.52 |
| homopolymer | 12 | 83 | 15,976 | **0.0094** | 0.754 | 0.758 | 0.60 |
| dinucleotide | 6 | 106 | 45,723 | 0.0017 | 0.799 | 0.520 | 5.25 |

**A read misreads a homopolymer of 8 repeats 2 times in 1,000 and one of 12 repeats 9.4 times in
1,000 — 4.7 times as often over four repeat counts, and it rises at every step.** Step by step the
ratio is 1.35, 1.37, 1.59 and 1.59, so **about 1.47-fold a repeat count**. The specification predicts
roughly **1.3-fold a repeat count**, inferred from stutter rates the per-sample route measured a
different way ([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3). **Nothing told
this estimator to expect a rise**: each of the five strata was fitted on its own tracts, in a separate
optimisation, with no monotonicity constraint anywhere in the code.

**The direction of the slip rises with the tract too.** The share of slipped reads showing a
*shorter* tract goes 0.595, 0.637, 0.623, 0.713, 0.754 — from 1.5 shorter for every longer at eight
repeats to 3.1 at twelve. Tomato's dinucleotides at six repeats sit at 0.799, four shorter for every
longer, against the 0.83 the per-sample route measured on the same crop from a per-read tally.

### 6.2 The concentration a real stratum carries — an open question in the specification, now answered

[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.1 closes with *"what is
not measured is which `κ` a real stratum has, and it cannot be until the STR records exist and a
genome walk fills them"*. It is measured here:

- **Tomato homopolymers: 0.52 to 1.56**, four of the five between 0.52 and 0.88.
- **Tomato dinucleotides at six repeats: 5.25**, ten times more spread out.

A small concentration means most tracts in the stratum are fixed at one length while the stratum as a
whole spans many. So tomato's homopolymers are largely monomorphic tract by tract and its short
dinucleotides are not — which is the direction a geneticist would predict and had not been measured.

**It also retro-validates the drawn work.** Every measurement of how few tracts a stratum can hold was
made at a concentration of 0.5
([`str_stratum_size_sweep_2026-08-13.md`](str_stratum_size_sweep_2026-08-13.md)), and tomato's
homopolymers come back at 0.52 to 0.88. **The dinucleotide stratum at 5.25 is outside the range any
drawn measurement covered**, so the floor of 1,000 tracts is not established for strata like it.

### 6.3 What the run also says, and it is less comfortable

- **Half the reads that reach a repeat tract never cross it**: 1,538,186 against 1,587,703 that do.
  That is a worse ratio than the trio's 37 in 100 (§5.1), on a crop with shorter tracts, and it says
  the effective depth at a tract is about half the depth at an ordinary position.
- **125 of the 4,164 tracts are over the guard's threshold** — one in 33 carries so many reads
  differing from the reference by a *non-whole* number of repeat units that the noise model does not
  describe them. They are left out and counted, which is what the records specification asks for.
- **65 of the 71 strata got no answer at all from their own tracts.** Every dinucleotide stratum above
  six repeats, every trinucleotide stratum, and every homopolymer above twelve repeats is below 50
  tracts. That is the population §4's borrowing exists for, and on this region set it is most of the
  table.
- **Three of the six fitted strata are themselves below the borrowing floor** — 350, 153, 83 and 106
  tracts — so their numbers carry the scatter §4 measures: at 100 tracts and three reads a site the
  fall-off moves ±21.5% between draws. **The monotone rise in §6.1 is stronger evidence than any one
  of those cells**, because five independent fits landing in order is not what scatter does.

---

## 7. What this costs to run, and the two things that made it affordable

A fit's time goes in two places, and which one dominates depends on the stratum:

- **The integral.** Building the 256-point rule over a stratum's length frequencies costs **19 ms at
  thirteen allele classes** and is rebuilt every time the climb moves the spectrum or the
  concentration — a few thousand times a stratum. On a thin stratum this is most of the run.
- **The likelihood sweep**, which is tracts × samples × genotypes × 256 points. On tomato's 63 samples
  this is the larger half.

Two exact savings are in the library. Neither changes an answer — §3's agreement was re-checked after
both.

1. **The integral is held while only slippage moves.** The climb asks about dozens of parameter sets a
   round that leave the tract's length frequencies alone.
2. **Strata that borrow the identical set are fitted once.** Where a whole motif length is thinner
   than the floor, every stratum in it pools every other, so the answers are one object computed as
   many times as there are strata. On the trio that collapsed 32 fits to four.

**The allele span is the expensive knob.** The specification's `±6` is thirteen classes, 91 genotypes
and twelve stick-breaking dimensions; the records' own `±4` is nine, 45 and eight — about a third of
the time. The trio was run at `±6`; the tomato cohort at `±4`, and the run says which in its own
output. What the wider span buys on real records was not measured.

**What the two real runs actually cost**, both on a machine carrying other work at the time:

| | strata | tracts | samples | repeat-tract fit | whole run |
|---|---:|---:|---:|---:|---:|
| GIAB trio, `±6`, borrowing below 1,000 | 32 | 216 | 3 | 99 s | 467 s |
| tomato cohort, `±4`, no borrowing | 71 | 4,164 | 63 | 1,690 s | 4,539 s |

**Borrowing is what makes a large cohort expensive, not the tracts.** With borrowing off, each of
tomato's 71 strata reads only its own tracts, and the fit is 1,690 s. With the floor at 1,000 every
thin stratum would pool its neighbours up to a thousand tracts, so the same 4,164 tracts would be read
tens of thousands of times over — the arithmetic says hours, and that arm was not run.

---

## 8. Changes I would propose to the specifications

None of these were made — the specification and architecture documents were being edited elsewhere
while this was written.

1. **[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.2** — where it says
   "a thin stratum borrows from its neighbours", give the two numbers §4 of this report measures:
   **borrow below 1,000 tracts**, and **take both sides of a repeat count together rather than one
   neighbour and then a test**. The second is not a detail: it is the difference between a borrowed
   slippage level 1.2% from the truth and one 23% from it.

2. **[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4** — require the
   emitted per-stratum parameters to carry **how many tracts stood behind each answer and which strata
   were pooled into it**. §5.3 is what happens without it: fifteen homopolymer strata emitting the
   identical four numbers, with nothing in the output to say the repeat-count axis had been averaged
   away. The library records both; no specification asks for them.

3. **[`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §1.1** — the
   *reached but did not cross* state is specified as a stored field and priced at nothing. On the GIAB
   trio it is **102,308 reads against 171,930 that crossed — 37 in every 100 that arrive**. Say that
   the count must be **reported per stratum**, because the censoring runs along the repeat-count axis
   the parameters are fitted within, and a stratum unreadable at this read length must not look like
   one that was unlucky with coverage.

4. **[`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md)** — add the check
   this work found missing: **how many tracts a stratum holds is a property of the analysed regions,
   not of the genome**. A 452 kb region set holds 216 tracts in 32 strata and can fit exactly one
   thing; the same reference walked whole holds 462,701 in 141. A run should say, before fitting,
   how many strata clear the floor on their own.

5. **[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.1** — carry over
   the size sweep's proposal that has not been made yet: the concentration is counted in **tracts and
   not reads**, so it, and not the slippage numbers, is what sets how many tracts a stratum needs, and
   a deeper cohort does not relieve it.

6. **[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.1, last line** —
   close the open question. It says which concentration a real stratum has "cannot be measured until
   the STR records exist and a genome walk fills them". They exist and it walked: **tomato
   homopolymers carry 0.52 to 1.56 and its dinucleotides at six repeats 5.25** (§6.2). Record with it
   that the drawn work was done at 0.5, which covers the homopolymers and not the dinucleotides.

7. **[`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3** — the monotonicity rule
   there merges strata whose fitted levels dip along the repeat-count axis, and its own text notes
   nobody had checked how often it fires when the truth is monotone. On tomato this estimator returned
   a **monotone rise across five consecutive homopolymer strata with no constraint applied at all**
   (§6.1), so on that cohort the rule would never have fired. That is one cohort and one motif length,
   but it is the first evidence either way.

---

## 9. What is still open

- **Borrowing was never run on the tomato cohort.** The trio ran with it and the tomato cohort
  without, so the two real arms differ in more than their data. The tomato borrowing arm is the one
  that would show what borrowing does to a repeat-count axis that genuinely varies — §6.1 measures the
  axis, §5.3 shows borrowing erasing one, and no run has yet done both at once. It is the next thing
  to run, and §7 says why it is expensive.
- **The one-sided borrow is unmeasured.** A stratum at either end of the repeat-count range has
  neighbours on one side only, and §4's argument does not reach it. `borrow-edge` mode exists and was
  not run. On tomato this is not a corner case: every homopolymer stratum above twelve repeats has
  neighbours only below it.
- **The crossover was measured at three tract counts** — 50, 250 and 1,000 — and borrowing won at all
  three, so no crossing point was found. The remaining rows (100, 500, 2,500) were dropped to give the
  real-data runs the machine.
- **Nothing here fits slippage per read group.** Both cohorts were run with every read group pooled
  into one slippage group, because 63 single-read-group samples would otherwise ask 189 slippage
  numbers of a stratum holding a few dozen tracts. The per-read-group grain the specification names is
  implemented and was not exercised.
- **The substitution rate inside a tract is not fitted.** The records carry a mismatch list and a
  denominator for it; this estimator reads neither.
- **Contamination from the STR loci** stays where spec §4.3 leaves it.
