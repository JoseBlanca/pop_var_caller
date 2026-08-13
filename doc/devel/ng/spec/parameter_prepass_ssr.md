# ng — the parameter pre-pass: the STR path

*Design spec, 2026-08-03. **Revised 2026-08-06 against measurement**, and the measurements changed
the accumulator's key rather than only its constants. **Revised again 2026-08-07**: §5.1 records
that the copy floors cannot be swept downward far by re-typing. **Revised 2026-08-08**: the archive
survey of 2,457 tomato libraries delivers the per-period floors table, and it is a measurement rather
than a decision — §5.1 records what adopting it would cost. **Revised 2026-08-09**: §4.5 adds what
has to hold where a stratum barely stutters — a second floor, on slipped reads rather than on loci,
because the level and the two shares starve at rates 20,000-fold apart; and agreement with the
generic path's substitution rate, which is obligatory there because the two noise models coincide
when nothing slips. **Revised 2026-08-11** in three places, none of which moves the noise model:
nothing detects a repeat while a sample is being read any more — the loci arrive from a catalog
built once per reference (§5.3), which is also why §5.1's sweep-and-re-type instrument no longer
exists; the copy floors §5.1.1 measured are now the ones the code applies, so §4.4's stratum count
follows them; and §4.6 asks a question nobody has asked of this path — what a locus that is not the
tract the stratum thinks it is does to the four numbers.*

*Source for the 2026-08-06 numbers:*
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§6. **No code yet — this settles the design.** One of five documents covering ng step 4. The shared
framing — why parameters are estimated without calling genotypes, and how the maximum-likelihood fit
works — is in [`parameter_prepass.md`](parameter_prepass.md), which this assumes and does not repeat.
**Scope: the noise model for repeat tracts, what it measures, and what the walk must accumulate for
it.** The empirical figures are measured on 51 tomato read groups (8.1M observations) and on HG002
(whole genome), both with ng's default delimiter, and are against each unit's own modal allele.
`src/ssr/` and `src/pileup/` are frozen production: everything said about them here is a record, not
a change.*

> **Every statistical claim in this document carries one of three marks**, so that a reader can tell
> what has been decided from what has been argued:
>
> - **Measured** — an exact-bias measurement stands behind it, with the number and where to find it.
> - **Measurable now** — nothing blocks the measurement; the experiment is named and the harness
>   exists.
> - **Deferred** — it needs something that does not exist yet, and the home is named.
>
> The empirical sections (§2, §3, §5, §7) are measurements on real data and carry their own numbers.

---

## 1. What this path fits, and what it borrows

**The core procedure is not this document's**, and that is the point of the split. Summing over the
genotype instead of choosing one, and maximising the result, are
[`parameter_prepass.md`](parameter_prepass.md) §3 and apply here unchanged. What differs is the
**noise model** — what a read is assumed to get wrong — and, because this path stratifies its loci
and the generic one does not, two implementation choices that only make sense where strata exist
(§4.3).

**The generic path's noise model has one parameter and cannot express what happens here.** Its `p_j`
says a read shows the wrong allele with probability `ε`, a per-base substitution rate, and nothing
else. That is right where the alternatives to a reference base are three other bases. A repeat tract
can also **slip**, showing a whole copy more or fewer than the allele carries, which `ε` has nowhere
to put.

**A read at an STR locus is noisy in two independent ways, and keeping them apart is what makes both
estimable.** The tract can slip, changing its **length** by whole copies; and bases can be misread,
changing its **composition** at fixed length. Mark-2 relies on exactly that separation —
*"`θ` changes length in whole units, `ε` changes composition via substitutions"*
([`../../specs/ssr_cohort_mark2.md`](../../specs/ssr_cohort_mark2.md) line 289). So this path's noise
model **contains** the generic one and adds slippage:

| | generic path | STR path |
|---|---|---|
| observation | reads supporting the alternative, out of depth | reads at each whole-repeat offset from the reference tract length, out of depth |
| what one accumulator entry holds | one site | **one locus** (§4.1) |
| genotype summed over | one genotype per site, dosage 0…P | the pair of allele lengths |
| noise parameters | `ε` (substitution) | `ε` **and** slippage: a level, a direction split and a distance decay, per stratum |
| stratified by | nothing | motif period × repeat count (§4) |

### 1.1 Four numbers, not three — and an earlier version of this section counted wrong

**This path fits four things per (read group × stratum):**

- **how often a read slips at all** — the *level*, `P(a read shows a length other than its allele's)`;
- **which way it slips**, shorter or longer than the allele, which is strongly asymmetric (§3);
- **how far it slips when it does**, which decays, and decays the same way in both directions (§3);
- **the substitution error**, which changes composition at fixed length.

**An earlier version of this document listed three and left out the level**, as did
[`parameter_prepass.md`](parameter_prepass.md) §1's summary table and §3.1's "scan its three noise
parameters". The level is not an optional fourth: it is the quantity §4 stratifies by, §4.3 holds
monotonic from one stratum to the next, and §5 tabulates at 0.091%, 0.170% and 2.006% across
repeat-count bands. A model with no level cannot express any of those sentences. **Measured** — the
correction was forced by writing the model down for the harness
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.1).

**Decision: this path's substitution error is fitted separately from the generic path's, and the two
are never tied.** Both are called a per-base substitution rate, and it is tempting to say a sequencer
misreads a base at one rate and therefore one number should serve — the same read group, the same
chemistry, the same machine.

**The reason not to is that each number is the error parameter of *its own model*, not a measurement
of the machine.** A model is a description, not the thing described. The generic path's rate absorbs
everything its model cannot otherwise explain at an ordinary site; this path's absorbs everything the
slippage model cannot explain inside a tract. Those are different residuals, so forcing one number to
carry both would make each model wrong in a way neither could report. Whichever model is the more
constrained would quietly wear the difference.

**What it costs is one parameter, fitted from counts already accumulated (§4.1).** What it buys is
that the two numbers can be *compared*: if the error rate inside repeat tracts comes out higher than
outside — which is plausible and, as far as these specs record, unmeasured — that is a finding.
Tying them destroys the observation before it can be made.

**And the comparison has a place where it is not merely interesting but obligatory.** Where a
stratum barely slips, this path's noise model degenerates into the generic one, so the two rates have
to meet or one of them is wrong. §4.5 makes that a requirement and fixes what "meet" means.

*Naming:* both are written `ε` in their own formulas, because each is its model's error parameter and
no formula contains both. In prose this document says **the STR path's substitution error** and
never assumes it equals the generic one.

### 1.2 What GATK's DRAGstr gives us, and what it does not

**DRAGstr is the closest existing thing to this document**, and it is worth being precise about which
half of it transfers. `DragstrParametersEstimator.java` (vendored) fits an STR error model without
ever calling a genotype — the precedent [`parameter_prepass.md`](parameter_prepass.md) §3 rests on —
and its search is worth seeing concretely. It fixes a grid of candidate parameter values: 41 error
values crossed with 41 values of the variant prior. For **each** of those 1,681 candidates it
evaluates the three-term mixture **at every locus in the stratum**, adds the logarithms, and records
the sum. Then it keeps the candidate whose sum came out largest. Nothing is being made larger — the
sums are what they are, and the search is a choice **among** candidates.

**"At every locus" is the half that matters most, and it is the half this document got wrong.** A
locus's reads share a genotype, and DRAGstr's unit of evidence is the locus for that reason. §4.1
below adopts it, and the research note measures what pooling reads across loci instead costs.

**Two of its implementation choices are copied here**, and both are specific to a stratified path,
which is why they live in this document rather than the general one: thin strata borrow from their
neighbours, and the fitted sequence is held monotonic along the repeat-count axis. §4.3 carries them,
with what each costs.

**What is not copied is DRAGstr's noise model.** It collapses every kind of slippage into a single
indel-error rate per stratum and describes neither its direction nor its size. §3 measures both and
finds them worth fitting, so a model with nowhere to put them is not enough here.

### 1.3 HipSTR fits a different model, not a different estimator

**HipSTR marginalises over the genotype exactly as this design does.** It reaches the answer by a
different route — `em_stutter_genotyper.cpp` (vendored) alternates: start from a guess at the stutter
parameters; work out, for each sample, how probable each genotype is given that sample's reads; then
re-fit the parameters, letting every read count towards *every* genotype in proportion to how
probable that genotype turned out to be; repeat until nothing moves.

**That is expectation-maximization, and it is an algorithm rather than an alternative.** Its two steps
provably increase the same marginal likelihood that
[`parameter_prepass.md`](parameter_prepass.md) §3 evaluates directly. **Choosing between them cannot
change the answer beyond optimisation error** — but which optimiser is affordable *is* now a live
question on this path, because there are three slippage parameters rather than one and a flat scan
over all three is 4.2 million scores per stratum (§4.2).

**What differs is the model.** Working out how probable each genotype is requires knowing which
alleles are common, and HipSTR takes that from **the locus itself**, fitting each locus's allele
frequencies inside the same loop. This design marginalises against **per-stratum** parameters
pooled across many loci. Those are different models:

| | per-stratum (here) | per-locus (HipSTR) |
|---|---|---|
| what the genotype is weighed against | one set of allele frequencies per *(period, repeat count)* | this locus's own allele frequencies |
| can describe a locus that behaves unlike its stratum | no | yes |
| data behind each parameter | every locus in the stratum | one locus, across samples |
| needs several samples at one locus | no | **yes** |

**The last row is why the comparison was never run**: a per-sample walk had no cross-sample evidence
at a locus, so only the per-stratum model was available. **Nothing here claims it is more accurate**
— it claims it fitted the shape of the pass, which is a weaker and more honest thing to say.

**The STR census removes that obstacle** ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)
§1.1), supplying several samples at one repeat tract. One caveat before anyone reads that as a plan:
a census holds far fewer loci than the per-stratum histogram holds loci, so a per-locus fit over it
is not being handed equivalent evidence. §8.5 carries the comparison, with no leaning.

---

## 2. Slippage moves whole repeats

Out of every 100 reads that differ from the allele, 98 differ by a whole number of motif copies at
dinucleotides, 95 at trinucleotides, 93 at hexamers. **Homopolymers are 100% by arithmetic — every
integer is a multiple of one — and are therefore no evidence either way.**

The remainder is not slippage: it is one- and two-base indels, which this noise model has no way to
describe. That residue is small where tracts are long (3 in 1,000 reads at ≥6 repeats) and large
where they are short (§5).

---

## 3. Which way a read slips, and how far — one of these is asymmetric, the other is not

**A read that has slipped differs from the true allele in two separate ways, and the model needs a
parameter for each.** It slipped in some **direction** — the tract it shows is shorter or longer than
the allele — and by some **distance**, usually one repeat, occasionally two. Whether the up and down
cases behave alike is a different question for each, and the answers differ.

**Direction: strongly asymmetric, and the model must carry it.** Reads showing a *shorter* tract than
the allele are far more common than reads showing a longer one, and the imbalance grows with the
motif period:

| | reads showing a shorter tract | reads showing a longer one | ratio |
|---|---:|---:|---:|
| tomato, homopolymer | 5,072 | 3,592 | 1.4× |
| tomato, dinucleotide | 2,438 | 501 | **4.9×** |
| human, homopolymer | 3,037 | 1,640 | 1.9× |
| human, dinucleotide | 545 | 162 | **3.4×** |

At tomato dinucleotides a read is nearly **five times** as likely to have lost a repeat as to have
gained one. **This is a fitted parameter, one per stratum** — the "direction split" of §1.1 — and
[`parameter_prepass.md`](parameter_prepass.md) §2.2 is the measurement of what happens when it is
fitted on uncontrolled data: it inverts, and reports gains as marginally *more* common than losses.
**§4.1 now shows the same inversion arising from a keying choice rather than from thresholding**, at
the same size.

**Distance: the decay is the same in both directions.** Of the reads that slipped by one repeat, some
smaller number slipped by two instead. That ratio is the **fall-off** — read `0.065 (5,072 → 329)` as
*5,072 reads slipped by one, 329 slipped by two, so about 7 in every 100 took the second step*. The
question this table answers is whether a read that is losing repeats decays differently from one
gaining them:

| | fall-off when losing | fall-off when gaining | gap vs its counting error |
|---|---|---|---|
| tomato, homopolymer | 0.065 (5,072 → 329) | 0.074 (3,592 → 265) | 1.5 SE |
| tomato, dinucleotide | 0.087 (2,438 → 211) | 0.102 (501 → 51) | 0.9 SE |
| human, homopolymer | 0.097 (3,037 → 296) | 0.115 (1,640 → 188) | 1.6 SE |
| human, dinucleotide | 0.106 (545 → 58) | 0.123 (162 → 20) | 0.5 SE |

**Decision: one fall-off parameter shared by both directions.** The reason is the thinness of the
gaining arm: above dinucleotides it rests on 3 to 13 reads, so giving it a free parameter would fit
counting noise rather than capture a difference.

**So the two decisions together.** A stratum carries **an asymmetric direction split** and **a
symmetric distance decay**. In plain terms: yes, the model expects to see more reads short of the
allele than long of it — several times more at dinucleotides — but once a read has slipped, the
chance of it having slipped twice rather than once is taken to be the same either way.

**Read the fall-off table carefully, though — it does not say there is no difference.** Each gap sits
inside its own counting error, at 1.5, 0.9, 1.6 and 0.5 standard errors. But the gaining side decays
faster in **all four** rows, and four consistent signs pool to roughly 2 SE. The finding is "no
difference we can afford to fit", not "no difference". If the fall-off ever becomes a function of
level (§8.2), revisit this with the pooled test rather than the per-row one.

**The fall-off value does not transfer between datasets.** About 10 reads in 100 take a second step in
human against about 7 in tomato. The structure is portable; the number has to be fitted.

---

## 4. Stratify by repeat count, not by base length

**How much a tract slips depends on how many repeats it has, more than on anything else the data
offers** — 9 reads in 10,000 below four repeats against 2 in 100 at six or more, a twenty-two-fold
spread within one dataset (§5). A model that does not condition on repeat count is averaging across
that range. So repeat count is a stratification axis; what this section settles is whether it beats
the obvious alternative.

The Mark-2 spec fits the stutter level as linear in tract **length in bases**
([`../../specs/ssr_cohort_mark2.md`](../../specs/ssr_cohort_mark2.md) §4.4). The data says repeat
**count** is the better axis, which is what a per-copy slippage mechanism predicts: ten copies offer
ten chances to slip whether the copy is 2 bp or 6 bp.

At 12–15 repeats, tomato homopolymers, dinucleotides and trinucleotides stutter at 14.3%, 15.0% and
8.6% — within a factor of two of each other. On a base-length axis the same three periods at
20–29 bp were 12.9%, 12.6% and 1.3%, an order of magnitude apart. Homopolymers also become monotonic
on the copy axis and are not on the length axis.

**Decision: stratify by (period, repeat count).** Same parameter count as today, better behaved.

**The monotonicity that §4.3's second borrowed choice depends on is a property of this axis, not of
the other one.** On the copy axis homopolymer slippage rises with every step; on the base-length axis
it does not. So the constraint "a longer tract must not come out less slippery than a shorter one" is
something the data supports here and would misfire if the strata were keyed by length.

### 4.1 What the walk accumulates

**Still a histogram, and still one that forgets which locus was which. What changed is what it
counts: loci instead of reads.**

One table per `(read group, period, repeat count)`. Each locus is reduced to its **shape** — how
many of that locus's reads fell at each whole-repeat offset from the reference tract length — and
the table counts *how many loci in this stratum had each shape*. **Two loci that looked alike
collapse into one entry with a count of two**, exactly as the generic path's histogram counts sites
that looked alike; which loci they were is never asked again, and neither their coordinates nor
their sequence ever enters.

Three loci at four reads each make the difference concrete. Two of them show all four reads at the
reference length; the third shows three there and one a repeat short:

| | what the table holds |
|---|---|
| **counting reads** — the earlier version | one row for the whole stratum: *11 reads at the reference length, 1 read one repeat short* |
| **counting loci** — this version | two entries: *"all four reads at the reference length" — 2 loci*, and *"three at the reference, one a repeat short" — 1 locus* |

The read row is equally consistent with two stories the fit has to tell apart: every locus is
reference-length and one read slipped, or one locus's own allele is a repeat short and it happened
to get one read. **Keeping a locus's reads together is what distinguishes them**, and the collapse
is real rather than nominal — at HG002's 300× the table holds 0.43 entries per locus, so most loci
share a shape with another ([research note](../research/parameter_estimator_experiments_2026-08-06.md)
§6.8, and the memory paragraph at the end of this section).

Each entry therefore holds:

- **how many of that locus's reads showed each whole-repeat offset**, over a bounded range with the
  end buckets saturating (below);
- **how many showed something that is not a whole number of copies** — the guard bucket (below);
- and, pooled across the stratum rather than per locus, **how many bases were compared and how many
  of them mismatched**. Two running counts, and they are what `ε` is fitted from.

**Why the locus and not the read, and this is the finding that reopened the section.** A read carries
no genotype: it drew one of the locus's two alleles and then slipped. Pool reads across loci and what
the table holds is the **allele spectrum convolved with the slippage kernel** — and recovering the
kernel from that means undoing a convolution with both halves unknown. **Measured:** with the allele
spectrum handed to the fit, a per-read tally is exactly unbiased; with the spectrum fitted from the
same tally, the answer for the slippage level moves **333-fold depending only on where the search
starts**, on a stratum whose alleles span three repeats either side of the reference
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.2). That is the
signature of a quantity that is *not identified* rather than badly estimated — the same signature
§2.2 of the research note records for the generic path's pooled multi-library key. Keying by locus
returns the level to 0.000% bias with the same four starting points agreeing to 1.000×.

*The other way out of it is closed for a reason worth recording.* The same measurement shows a
per-read tally is **exactly unbiased when the allele spectrum is handed to it** rather than fitted
from it, so a design that got the spectrum from somewhere else could keep the cheaper object. There
is nowhere else to get it during the walk: the STR census is the only object that holds allele
lengths across samples, and it stays raw until the gather
([`parameter_prepass.md`](parameter_prepass.md) §1.3), which is after every stutter parameter has
been emitted. **Deferred rather than rejected**, and it is the same shape as §8.5's per-locus
question — both become askable only once something upstream of the gather can supply a spectrum.

**This is the STR twin of a rule the generic path already carries**, and the correspondence is exact:
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §1 rejects a windowed histogram keyed
per read group because splitting one site between two entries lets each entry draw its own genotype
independently. Pooling reads across loci is the same mistake with the site dissolved entirely.

**A locus covered by two read groups still makes two entries, and that one is safe.** It looks like
the same mistake — a locus has one genotype and its reads are being split — so it is worth saying
exactly why it is not. **Measured:** each read-group entry's own distribution is correctly specified,
because the genotype is drawn once for the locus and enters both entries through the same mixture.
The product over them is what the literature calls a *composite likelihood*: every factor is a true
marginal, so the estimator stays consistent and what the split throws away is the dependence between
a locus's entries, which is precision. The generic path measured exactly that split — fitting
everything from its read-group table alone — and found it unbiased in all 25 worlds tried
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §2.6). **What must not be
split is a locus's reads within one read group**, which is the paragraph above; and 1,550 of the
1,707 samples in the tomato archive survey carry one read group anyway.

**Decision: the offsets are measured from the reference tract length.** This closes the `OPEN` that
this section and §8 carried before this revision, and it closes it the other way from the leaning
they recorded.

*What the old leaning was, and what it cost.* Every measurement in §2, §3 and §5 uses each unit's own
**modal observed length** as the origin, and the accumulator used to be written to match, on the
grounds that centring each locus on its own mode keeps the offset range tight. The doubt was stated
in the same paragraph — *the mode of a heterozygous locus is not the allele* — and answered with
*"the origin is a binning choice and not a genotype call: the fit marginalises over the genotype, so
a heterozygous locus's second allele is explained by the genotype term rather than charged to
slippage"*, followed by *"that is the claim, and it needs checking rather than assuming."*

**Measured, and the claim is false as stated** ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.3).
The modal origin makes the origin a function of the reads, so a fit that treats it as a fixed
property of the locus is misspecified — and marginalising over the genotype does not repair it,
because the genotype term is now being asked to explain a quantity that moved with the data:

| how the accumulator was keyed / what the fit believed | slippage level | direction split | heterozygosity |
|---|---:|---:|---:|
| reference origin, scored as such — at every depth and heterozygosity tried | **0.0%** | **0.000** | **0.0%** |
| modal origin, scored as a fixed origin — 3.2 reads a locus, 10 loci in 100 heterozygous | **+50%** | +0.166 | −62% |
| modal origin, scored as a fixed origin — 3.2 reads a locus, 46 in 100 heterozygous | **+408%** | +0.311 | −75% |
| modal origin, scored as a fixed origin — 9.6 reads a locus, 46 in 100 heterozygous | **+90%** | +0.225 | −12% |

**The direction split is the row to read twice.** The truth in every world is a split of 0.17 — a
read is 4.9 times as likely to lose a repeat as gain one, which is exactly what §3 measures at tomato
dinucleotides. Centring on the mode and scoring it as a fixed origin returns 0.48 at three reads with
46 loci in 100 heterozygous: **a 1.1-fold asymmetry where the truth is 4.9-fold.** That is
[`parameter_prepass.md`](parameter_prepass.md) §2.2's real-data confound — 0.9× on all loci against
3.4× on known-homozygous ones — reproduced exactly, from a keying choice rather than from
thresholding, and it is the same size.

**The bias shrinks with depth and does not go away**, because a heterozygous locus's mode is the
wrong origin however many reads confirm it: still +90% at 9.6 reads a locus.

*A fit that models the centring is unbiased, and still does not pay for itself.* Scoring the
mode-centred table by summing over everything the centring forgot returns the level to within 1.1%
and the direction split to within 0.004, and the fall-off is recoverable too — a residual that
looked like bias turned out to be the measuring harness's own climb stopping short
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.3.1). What it costs
is two things the reference origin does not. It needs an enumeration over every way a locus's reads
could have landed, where the reference origin scores a locus directly. And **it makes the climb over
the genotype frequencies converge a thousand times more slowly** — mode-centring leaves a locus's
shape saying much less about which alleles it carries, and that overlap is exactly what sets an
expectation-maximization scheme's rate. The reference-origin climb converges in under 200 passes on
the same data; the mode-centred one has not converged in 200,000. That climb sits inside the inner
loop of a search that runs once per stratum. **There is no case left for the modal origin**, and
one further reason to be glad of it: the STR census already keys by the reference tract length
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §2.1), so the two objects
now share an origin and §4.2 of [`parameter_prepass.md`](parameter_prepass.md)'s comparison between
them is a comparison of the same quantity.

**The offset range saturates at its ends, and the end buckets are scored by their marginal.** "At
least four repeats short" is one bucket, and its probability is the sum over every offset it absorbs
— never the probability of sitting exactly on the edge. **Measured:** the marginal is exactly
unbiased at every range tried, ±1, ±2 and ±3, and costs a sum over a handful of kernel terms. Plugging
in the edge instead fails the algebraic gate outright — the bucket probabilities sum to 0.9488 at ±1
rather than to one — and, rescaled so they do sum to one, costs **+33% of the slippage level** where
30 in 100 slipped reads take a second step
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.4). That is the regime
long tracts sit in, so the error is largest exactly where slippage matters most.

**Against the reference origin the end buckets absorb whole alleles, and the marginal rule survives
that too.** A read's offset is now *its allele's distance from the reference length* plus the slip,
so a narrow range folds alleles into its ends rather than only far slips — which looks like the one
thing the reference origin costs. **Measured, and it does not cost it:** on a stratum whose alleles
reach three repeats either side and where 30 loci in 100 carry one of them, a range of **±1** still
returns the slippage level to within 0.05% and both shares to within 0.002
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.4). The same rows
score a plug-in at **−52% of the level**, so this is where the scoring rule earns its keep rather
than where the range does.

**Two widths, and only one of them is load-bearing.** The *recorded offset range* can be narrow; the
*allele lengths the fit is allowed to place mass on* must reach beyond it, because that is what lets
the marginal rule attribute an end bucket to an allele rather than to a far slip. What a narrow
recorded range costs is the heterozygosity that falls out of the fitted genotype frequencies — 1.5%
at ±1, 0.7% at ±2 — and heterozygosity is a by-product here, not something this path emits.
**A recorded range of ±4 is the working value**, chosen so that the ends absorb little on ordinary
strata rather than because anything forces it. The fitted allele support is the other width, it is
the load-bearing one, and §8.1 settles it at **±6** —
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) carries both as
`OFFSET_HALF_RANGE` and `ALLELE_OFFSET_LIMIT`.

**A locus's reads are capped, and the cap is a subsample.** The number of distinct entry shapes grows
with a locus's depth, so a deep locus is entered from a random subsample of its reads down to the
cap. Losing reads costs precision and no correctness: a uniform subsample of a locus's reads is a
thinning of the same multinomial, so the entry is distributed exactly as it would be at the lower
depth. Seed the draw from the locus's position, so a region-sharded walk and a single-threaded one
keep the same reads and merging stays exact — the same rule and the same reason as
[`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §2.2. **The cap's
value is not an open design question** (§8.8): both things that would have made it one are measured
away below — the table's size and the scoring rule's depth limit — leaving a precision trade against
the width of the entry's counters. [`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md)
§2.1 carries the working value.

**The guard bucket factorises, which is why it can be a diagnostic rather than a parameter.** A read
that differs from the allele by something that is not a whole number of copies is modelled as an
independent per-read outcome, so the likelihood splits exactly into *how many reads were
non-whole-repeat* times *how the rest fell across the offsets*. Nothing about the slippage parameters
is estimated from it, and nothing about it disturbs them. **Measurable now, and it is the assumption
under the split, not the split itself:** if a non-whole-repeat outcome is more likely at some allele
lengths than at others, the factorisation stops being exact — which the harness can test by making
the guard rate depend on the allele.

**And it needs a threshold, which an earlier version left as "a stratum where it is large" — a field
nobody reads.** The threshold is **one non-whole-repeat read in ten of the reads that differ from the
reference tract length**, and §5 is where that number comes from rather than being chosen for
looking round.

**Which reads the ten are counted out of is worth pinning, because the model and the diagnostic use
two different denominators and only one of them is computable here.** The *model* statement above is
about reads that differ from the **allele** — that is what a slip is. The *diagnostic* cannot be:
the accumulator does not know the allele, only the reference tract length, so the share it reports
is over reads differing from the **reference**. The two coincide at a locus whose alleles are the
reference length and diverge otherwise, always in the same direction: a real non-reference allele
contributes whole-repeat differences to the denominator and none to the numerator, so **the reported
guard share is diluted relative to the model's**, never inflated. §5's table is measured on
known-homozygous loci at the reference length, where the two are the same number, and the research
note's real-data figures — 1.37% on HG002, 1.97% on tomato — are over the reference denominator
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.8). A stratum that
crosses the threshold on the diluted number has crossed it on the model's too.

**The composition channel is a division, not a search.** Each read is compared against the tract at
**the length that read shows**, so a mismatch is a substitution and not a slip, and a read's mismatch
count is binomial at `ε` whatever length it showed. **Measured:** the two pooled counters are a
sufficient statistic and the maximum-likelihood rate is mismatches over bases compared — recovered to
four decimal places by a search that had no need to run
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.7). Where a stratum
holds reads of two different true rates the pooled counters return their base-weighted mean, which is
the right answer for a model carrying one rate and the same behaviour the generic path's `ε` has.

**A subtlety worth stating, because it looks like the thing this spec forbids.** Counting mismatches
means comparing each read against the tract sequence at the length that read shows — not at a
genotype the fit has chosen. It is an alignment, not a call, so it does not reintroduce the
threshold-then-count bias of [`parameter_prepass.md`](parameter_prepass.md) §2.

**A constraint this path imposes on read admission upstream, and the reason is here rather than
there.** Reads are admitted on a quality test before any of this
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2), and **base quality is
systematically worse inside repeat tracts** — so a test sensitive to a run of bad bases would drop
reads at STR loci preferentially, and among those the long tracts and the slipped reads first. That
would bias every number in this document downward, slippage most of all. It is why the admission test
is a mean rather than a quantile, and why that document carries a check comparing the admitted
fraction at STR loci against elsewhere.

**Keyed by read group because slippage is chemistry**
([`parameter_prepass.md`](parameter_prepass.md) §5), and by `(period, repeat count)` because §4 is
the argument that this is the axis slippage varies on.

**What it costs in memory: measured, and it is not the concern the revision first raised.** Making
an entry a locus rather than a read takes the table's size with it — a tally of reads is a handful
of buckets a stratum, while a table of locus shapes grows with how many *distinct* shapes the loci
take — so this section first said the size ran between two bounds, with nearly every locus its own
entry at 300×. **Measured on HG002 at 300× over the 50,000-interval GIAB tandem-repeat set, the
uncapped table is 12,727 entries for 29,811 loci — 0.43 entries a locus, 0.36 MB**
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.8). Deep data
deduplicates because most loci at a clean tract are "every read at the reference length", so what
separates two entries is mostly their depth, and depths repeat. **Entries saturate rather than
scaling with the genome**: a whole tomato genome's 1.73 million STR loci make 70,305 entries —
2.01 MB — and a whole human genome extrapolates to about 120,000 entries, **3.5 MB**, beside a
windowed histogram [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §9 prices at
115 MB per human sample. **This object is not where step 4's memory goes**, and the read cap is not
what keeps it from being.

**Which leaves the read cap less to do than it was given.** It was introduced to hold the table
down, which it does not need to do; and it looked like it also marked where the evidence for the
scoring rule stopped, at 12 reads a locus. **It no longer does: the rule is exactly unbiased at
every depth to 45 reads a locus** — 0.00% on the level, 0.0000 on both shares, four starting points
agreeing to 1.00× ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.8).
Reaching those depths cost only a narrower recorded range, which §4.1's saturation measurement had
already shown loses nothing.

**So the cap is a precision trade and a counter width, not a correctness limit.** An entry holds its
bucket counts in single bytes, which a 300× locus overflows, so the cap and the counter width are
one decision; and reads dropped above the cap are precision the fit does not get. Neither is a
reason to keep it low, and a design that raises it is not stepping outside what has been measured.

### 4.2 How the four numbers are fitted

**The substitution rate first, because it does not need the others.** Mismatched bases over bases
compared, per stratum. One division (§4.1).

**Then the three slippage parameters, together, from several starting points.** The genotype is
summed over, so the fit alternates in the shape [`parameter_prepass.md`](parameter_prepass.md) §3.1
settles: hold the slippage parameters, climb to the allele-pair frequencies that fit the stratum's
loci best — the half with a concavity proof behind it, so it cannot stop on a false summit — then
move the slippage parameters and repeat.

**What the genotype is on this path.** An unordered pair of allele lengths, expressed as whole-repeat
offsets from the reference tract length, so a stratum whose alleles span `A` distinct lengths has
`A(A+1)/2` genotype frequencies — 45 at nine lengths. *Above diploidy it is a `P`-tuple and the
count is `C(A + P − 1, P)`, which is 495 at nine lengths and four copies; the accumulator does not
change and the fit's inner loop does, the same trade
[`parameter_prepass.md`](parameter_prepass.md) §3 makes for the generic path.* **They are fitted
freely**, matching the generic path's decision to
fit its genotype frequencies freely rather than tie them through one allele frequency
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §11.4), and for the same reason: a
Hardy-Weinberg tie presumes the inbreeding coefficient is zero, and the inbreeding coefficient is a
quantity this run measures rather than assumes.

*A cheaper parameterisation exists and is worth measuring before a thin stratum needs it.*
`A(A+1)/2` free numbers is 45 at nine allele lengths, against 8 for an allele spectrum plus the
sample's own `F` supplied from the generic path (`π(a,a) = F·p_a + (1−F)·p_a²`). That is not the tie
§11.4 rejects — `F` is supplied rather than assumed zero — and it would let a thin stratum be fitted
where a free set cannot. **Measurable now:** the harness fits both against the same truth and reads
the bias off. **Leaning:** free where the stratum can afford it, tied where it cannot, with the
choice recorded per stratum as provenance.

**Not a flat scan over the three, and this is a change from
[`parameter_prepass.md`](parameter_prepass.md) §3.1's decision.** That section steps through the
noise parameters end to end, and prices the STR path at "4.2 million, a single flat scan, end to
end". Two things are wrong with carrying it over unchanged:

- **The arithmetic prices the wrong three parameters.** `ε` is a division (§4.1), so what would be
  scanned is the level, the direction split and the fall-off. `161³` happens to be the same 4.2
  million, but a quarter-Phred ladder is the wrong spacing for two parameters that are shares in
  `(0, 1)` rather than rates spanning orders of magnitude.
- **It is per (read group × stratum).** The generic path scans 161 rungs once per read group. This
  path would scan millions of combinations several hundred times, each combination costing a climb
  over the stratum's entries rather than a single pass. Whether that is affordable is not a question
  the shared spec's arithmetic answers.

**Decision: search from several starting points, spread over all three slippage parameters, and keep
the best-scoring.** The reason the shared spec gave for a flat scan is that nobody has shown the
profile curve has a single hump; **measured, and it does** — profiling the level with the other two
parameters and the genotype frequencies maximised out gives exactly one interior maximum on 41 rungs
from 0.0001 to 0.3, in both worlds tried
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.5). **That is two
worlds and one axis profiled**, so it is evidence against a second hump in the level and says nothing
about the two shares — which is why the answer is several starts rather than one.

**The starting points must disagree about all three parameters, not only the level.** This is the
trap the generic path's inbreeding fit fell into: five starts that disagreed about the headline
number while sharing one guess at a nuisance axis returned a confident zero on a genome 29% covered
by runs ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §3.4). Every
number in §4.1's tables comes from four starts spread over the level, the direction and the fall-off
at once, and the spread across them is reported beside each fit.

**How sharply the level is determined varies enormously with the data, and that is precision rather
than bias.** At 3.2 reads a locus with 46 loci in 100 heterozygous the profile is nearly flat over a
three-fold range of the level; at 9.6 reads with 10 in 100 heterozygous the same span costs eight
times as much score. The observation count emitted beside each fit
([`parameter_prepass.md`](parameter_prepass.md) §6) is what a consumer has to tell those apart.

### 4.3 Borrowing and merging across strata, and what each costs

Two procedures copied from DRAGstr, both of which **change the estimate**:

1. **A stratum too thin to fit takes its neighbours' value** — adjacent repeat counts at the same
   period — rather than fitting noise. *"Too thin" is two tests and not one:* a floor on loci, which
   is this rule; and a floor on the reads that actually slipped, which protects the direction split
   and the fall-off and which a stratum can fail while clearing the first by four orders of
   magnitude (§4.5). Where only the second fails, only those two parameters are borrowed.
2. **The fitted levels are held monotonic along the repeat-count axis.** Slippage genuinely rises
   with repeat count (§4), so a fitted sequence that dips in the middle — tracts of 7 repeats coming
   out *less* slippery than tracts of 6 — is reporting the noise in one stratum, not a fact about
   repeats. Fit the strata in order, check each against the one before, and where a fit breaks the
   expected direction, **merge the two strata and refit them together**, repeating until the sequence
   behaves.

**Measured: what a merge costs is the merged stratum's distance from the pooled mean.** Two strata
pooled and fitted as one return close to the loci-weighted mean of their levels, so each then carries
its own distance from it: a **1.5-fold** difference between the two costs about **a quarter** of the
level, a two-fold difference about **half**, a four-fold difference up to **141%**
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.6). With two identical
strata a merge costs exactly nothing, which is the control.

**What that is worth on real strata.** §5's levels run 0.091% below four repeats, 0.170% at four to
five and 2.006% at six or more, and §4's dinucleotides reach 15.0% at 12–15 repeats — so slippage
rises roughly **1.3-fold per repeat count** over that range. Merging or borrowing across one repeat
count therefore costs on the order of **15 to 25% of the level**. That is a price worth paying for a
stratum that would otherwise be fitted on noise, and it is not a price to pay silently: every
borrowed or merged stratum carries `Provenance::Borrowed` and the strata it was fitted with
([`parameter_prepass.md`](parameter_prepass.md) §6).

**What is not measured: how often the monotonicity constraint fires when the truth is monotone.** A
merge triggered by sampling noise pools two strata that did not need pooling and charges both the
cost above. The exact method cannot see it — with cells weighted by their exact probabilities there
is no sampling noise to trigger a spurious merge — so this one needs draws. **Measurable now:**
simulate a monotone sequence of strata at each stratum's real locus count and count how often the
fitted sequence dips. **Leaning:** the constraint is worth keeping, because the alternative is a
stratum whose value is noise, and the failure it prevents is invisible while the failure it causes is
recorded.

### 4.4 What comes out beside the numbers, and why it has to be aggregated

**This path runs one fit per (read group × stratum), against the generic path's four fits in total.**
How many strata there are is bounded by arithmetic rather than guessed: six periods, each running
from its copy floor (`[8, 6, 6, 6, 5, 4]`, §5.1.1) up to as many repeats as a read can span, which at
150 bp reads is 150 for homopolymers and 25 for hexamers. That is 143 + 70 + 45 + 32 + 26 + 22 =
**338 strata**, of which the ones holding loci are fewer. *(Two earlier versions counted differently:
"about 370" counted each period from one repeat rather than from its floor, and 350 counted from the
unmeasured floors `[6, 4, 4, 3, 3, 3]` the code carried until 2026-08-10.)*

The generic path emits four things beside each fit: whether the answer sat on the edge
of its ladder, the starting points tried, the estimator's resolution, and how the fit terminated. At
four fits those are readable records; at several hundred they are a file nobody opens, and
[`parameter_prepass.md`](parameter_prepass.md) §3.1 already warns that a flag nobody reads is how a
badly-fitted parameter reaches a caller.

**So the walk emits a summary over strata, not a record per stratum.** The summary must answer, for
each read group:

- **how many strata were fitted in place, how many borrowed, how many merged**, and which — the
  merged sets by name, because a merge is a claim about two strata at once;
- **how many fits disagreed across their starting points**, with the worst offender named. This is
  the diagnostic §4.2's four starts exist to produce, and the one that separates a fitted number
  from a stopped search. *"Disagreed" needs a size, and this path has no ladder to read one off* —
  it searches rather than scanning, so there are no rungs. The size is borrowed from the generic
  path's ladder spacing, a quarter-Phred or 6%, which is that spec's argued estimate of the finest
  difference a caller can feel ([`parameter_prepass.md`](parameter_prepass.md) §3);
  [`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) fixes it as
  `START_AGREEMENT_LIMIT`;
- **how many strata carry a large guard-bucket share** (§5's threshold), with the worst named;
- **the observation count distribution across strata** — how many loci stood behind the thinnest fit
  and the thickest.

**A per-stratum record is still written**, because a fit that looks wrong has to be traceable; what
changes is that nothing downstream is expected to read it, and the summary above is what a person
sees. **This is an architecture decision as much as a spec one**, and
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) fixes the shape.

### 4.5 What has to hold where there is almost no slippage

**Nearly one locus in two arrives in a stratum where slippage barely happens, and the fit has to be
as trustworthy there as it is at six repeats.** §5 measures it: tracts below four repeats are 19.9%
of the loci and tracts of four to five another 28.2%, and their slippage levels are 0.091% and
0.170% against 2.006% at six or more. So in half the loci this path sees, fewer than 2 reads in
1,000 move at all. **Two requirements follow and they are different in kind** — the first is about
which of the four numbers survives thin evidence, the second is about agreeing with the other half of
step 4.

**First: the level and the two shares starve at very different rates, so one floor cannot protect
both.** §4.3's `MIN_LOCI_TO_FIT` counts loci, and loci are the wrong unit for the direction split and
the fall-off: those are measured only by the reads that actually moved. Take a stratum of 100,000
loci at 5 reads each — half a million reads — at the 0.091% level §5 measures below four repeats.
**455 reads slipped. Of those about 77 gained a repeat rather than losing one (§3's split of 0.17),
and about 5 of those 77 gained two (§3's fall-off of 0.065).** The locus count and the count that
measures the fall-off's gaining arm differ by a factor of 20,000, and §3's own tables already show
that arm resting on 3 to 13 reads above dinucleotides on real data.

**The level is unharmed by the same thinness, which is what makes the split worth making.** It is a
proportion over every read, not over the slipped ones: 455 out of 500,000 measures 0.091% to about
5% of itself, inside the 6% of `START_AGREEMENT_LIMIT`. The fall-off, on 5 reads, is measured to
about 45% of itself. **The same stratum measures one of its parameters well and another barely at
all.**

**So a second floor is needed and it counts slipped reads.** How large follows from the precision
wanted rather than from a preference: at §3's measured values, holding the direction split to 6% of
itself takes about **1,400** slipped reads and holding the fall-off to the same takes about
**4,000** — so the floor is on the order of a few thousand, and
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) fixes the value. **That the
fall-off is the parameter which starves first is a second reason for §3's decision to share one
fall-off between the two directions**, arrived at from precision where §3 arrived at it from the
thinness of the gaining arm.

**Below that floor the fit is split rather than borrowed whole: the level is fitted in place, the
direction split and the fall-off are borrowed from the nearest stratum at the same period that
clears the floor.** Borrowing the level as well would cost 15 to 25% of it per repeat count (§4.3),
and the level at the bottom of the range is exactly the number §5.1's copy-floor decision reads — so
it is the one thing not to throw away. The provenance records the two halves separately, because a
level fitted here with shares from two repeat counts up is a different claim from either a fit or a
borrow.

**At the bottom of the range the shares will mostly be borrowed, and that is the design working
rather than failing.** At a level of 0.091%, reaching 4,000 slipped reads takes about **880,000 loci
at 5 reads each** — against tomato's 1.73 million STR loci in total, no stratum below four repeats
comes near it. The alternative is not a better-measured share; it is a share fitted on 5 reads and
reported as though it were measured.

**Second: where nothing slips, this path's noise model *is* the generic path's, so the two
substitution rates must meet.** With the level at zero a read differs from its allele by substitution
and nothing else, which is exactly the generic model of
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §5.1. §1.1 fits the two rates
separately and refuses to tie them, and **that decision stands** — but fitted separately is not the
same as free to disagree. A stratum at 0.091% slippage is a place where the two models describe
nearly the same thing, so **they must agree there or one of them is wrong**, and §1.1's "the two
numbers can be *compared*" stops being an aspiration and becomes a test. It runs at step 4's own
surface, which is where a sample's two halves meet.

**The size is borrowed rather than derived: a quarter-Phred, 6% of the rate** — the generic path's
ladder spacing, which that spec argues is the finest difference a caller can feel and which §4.4
already borrows for `START_AGREEMENT_LIMIT`. **Soft, and honestly so:** nothing has yet measured how
close the two actually come, so the first run of this comparison is what sets the number.

**Three things move the two rates apart legitimately, and each is measurable rather than a tolerance
to widen:**

- **Tract impurity.** §4.1 counts a mismatch by comparing a read against the motif tiled to that
  read's length, so an interruption inside the tract is charged to this path's substitution rate and
  to nothing in the generic path's. `SsrSegment::purity_fraction()`
  ([`segment_criteria.rs:209`](../../../../src/ng/region_typing/segment_criteria.rs)) measures it per
  stratum, so the comparison runs on strata at full purity or states the offset it expects.
- **Which reads survived admission.** §4.1 records that base quality is systematically worse inside
  repeat tracts, so the two rates describe the same chemistry only if the same share of reads was
  admitted in both places — which is the check §4.1 already asks read admission to carry.
- **A real difference, which is the finding §1.1 exists to protect.** If the error rate inside tracts
  genuinely is higher, this comparison is how it is seen. A gap that persists across strata and read
  groups once the two above are accounted for is a result, not a failure — and the thing that must
  not happen is it being absorbed by widening the limit.

**What the exact-bias method cannot answer here, and at a low level it is the whole question.** The
harness weights every entry by its exact probability under the truth, so what it reports is bias with
no sampling noise in it — and near zero, sampling is precisely what is in doubt. The level is bounded
below by zero, so a finite stratum's estimate piles up against that boundary and its spread is not
symmetric about the truth: an estimator can be exactly unbiased in the harness's sense and still
return zero from a large share of real strata. **Measurable now, and it needs draws rather than the
exact method**, for the same reason §4.3's spurious-merge question does. **The experiment:** draw
strata at §5's measured levels and at the locus counts real strata have, fit each, and report how
often the level comes back at zero and how wide the spread is. **Leaning:** emit the slipped-read
count beside the level, the way §4.4's summary already carries the locus count — a level of 0.0003
standing on 4 slipped reads and one standing on 4,000 are different claims, and nothing downstream
can currently tell them apart.

### 4.6 Loci that are not the tract the stratum thinks they are — ADDED 2026-08-11

**A stratum's four numbers are pooled over hundreds of thousands of loci, and the pooling assumes
each of them is one tract of that period and that repeat count, read with that chemistry. Some are
not, for reasons that are biological rather than statistical**, and this section is the question of
what those do to the fit. It is asked because the sibling path asked it and the answer was expensive:
one substitution rate per read group described the body of the generic path's distribution and not
its tail — 818 loci carrying no benchmark variant showed three or more alternative reads where the
model predicted 29 — and because the only class that could explain such a locus was the heterozygote,
**the surplus came back as heterozygosity at 1.41 times the benchmark's count**
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1).

**Four mechanisms produce a repeat tract that behaves unlike its stratum, and they do not belong to
the same object:**

- **A duplication the reference does not carry.** Two paralogous tracts, of different lengths,
  collect their reads at one locus. A diploid genotype offers two allele lengths and the reads then
  show three or four, so no genotype in the support explains the shape. This is a property of the
  **genome**, so every library made from it shows it — the same first cause the generic path lists.
- **Mismapping from a similar tract elsewhere**, which puts a minority of reads at that other
  tract's length. A property partly of the **library**, since mapping difficulty follows read length
  and insert size.
- **Instability within the individual.** A long tract that mutates in the soma gives a spread of
  lengths wider than any slip kernel, and the spread is a property of the **locus**.
- **An interruption or a nearby indel** that makes the delimiter anchor a subset of the reads
  differently, so those reads' lengths move together by more than a copy.

**Which of the four numbers pays for them is not known, and this path has two absorbers where the
generic path had one.** A stratum's genotype frequencies are fitted freely over allele pairs
reaching ±6 repeats (§4.2), so a locus with reads at two distant lengths can be explained as a
locus carrying two distant alleles — which would corrupt the emitted allele spectrum and leave the
slippage parameters alone. Or the reads land where no pair of alleles reaches and the only
explanation left is slippage, which is what would move the level and flatten the direction split.
**Measurable now**, and the harness that measures it exists.

**One class of it is already measured, and it says this path is not robust by construction.** A
locus whose true allele sits outside the fit's support has its reads explained the only way left,
as slippage: leaving 2.5% of loci outside costs 0.1% of the level, 7.9% costs 2.5%, and **19.3%
costs +499% of the level with the direction asymmetry destroyed** — a split of 0.17 becoming 0.47
([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.4.1, and §8.1 of this
document, which reads it as the reason to treat the allele support as a threshold to clear rather
than a number to tune). That is the same collapse the modal origin produced, from a different cause.
A locus that is a collapsed duplication is not the same thing as a locus carrying a long allele, but
it fails the fit in the same place: **the model has nowhere to put it.**

**What must not be done about it, whatever the measurement says: dropping the loci that fit badly.**
A rule that scores each locus against the fitted model and discards the worst is thresholding on
the data and then counting what survives, which is the exact bias this whole step exists to remove
([`parameter_prepass.md`](parameter_prepass.md) §2). It would also remove the wrong loci first: a
real long allele and a collapsed duplication both sit where the model has least mass, so a trimming
rule would take genuine STR variation out before it took an artefact. **So if the check fails, the
fix is a class the model carries, exactly as the generic path's was** — a locus is *ordinary* with
probability `1 − w` and *aberrant* with probability `w`, an aberrant locus's reads drawn from a
much flatter distribution over the offsets, and `w` fitted rather than chosen. That costs this path
one further parameter and touches nothing shared: the noise parameters are an associated type of the
fitting seam, which is what let the generic path add its second site class without disturbing this
one.

**The diagnostic it asks of the fit is free, and it is worth having before any of this is decided.**
Every candidate already evaluates each entry's likelihood, and an entry is one locus's shape — so
the fit can report **what share of a stratum's loci its own fitted model calls very unlikely**
without keeping a coordinate, changing the entry key, or walking anything twice. §4.4's summary is
where it belongs, beside the guard share. **Measured per read rather than per locus**: a shape's
likelihood is a product over its reads, so a fixed floor on the total would flag every stratum's
deepest loci by arithmetic and say nothing about any of them.

**And the guard share does not already cover this**, which is the reason a second diagnostic is
needed rather than a threshold on the first. The guard share counts *reads* that moved by a
non-whole number of copies, pooled over a stratum (§4.1); a collapsed duplication moves its reads by
whole copies, and one locus in a thousand behaving strangely is invisible in a stratum-wide ratio.
The two diagnostics answer different questions: *is this model right for these tracts* against *are
these all the tracts this model was told they were*.

---

## 5. Short tracts carry almost no slippage, and most of what they carry is the wrong kind

**This section is evidence, not a decision.** Which tracts count as STR loci is settled upstream by
region typing; this step fits parameters from whatever loci it is handed. What follows is a
measurement whoever sets those boundaries will want, plus one thing it implies about how to read a
fitted value.

| repeats | share of loci | share of all slippage | slippage rate | of that, **not** a whole repeat |
|---|---:|---:|---:|---:|
| < 4 | 19.9% | 1.7% | 0.091% | **58.5%** |
| 4–5 | 28.2% | 4.5% | 0.170% | 33.8% |
| ≥ 6 | 51.9% | 93.7% | 2.006% | 0.9% |

*(HG002, known-homozygous loci. Tomato agrees in shape: 11.2% of loci, 6.0% of slippage, 19.8%
not-whole-repeat below four repeats.)*

Short tracts are a fifth of the loci, produce under 2 in every 100 slippage reads, and **nearly six in
ten of even that is not a whole-repeat change** — it is an ordinary indel, which §2's noise model has
no way to express.

**What this implies for the fitter is a diagnostic, not a rule.** A stratum whose slipped reads are
mostly *not* whole-repeat changes is one where the model does not describe the data, and a slippage
rate fitted there is mostly mis-modelled indel however much data stands behind it. That is
distinguishable from an ordinary thin stratum and should not be silently treated as one, so **the
fitter reports the non-whole-repeat fraction per stratum alongside the fitted values**.

**The threshold is one non-whole-repeat read in ten of the reads that differ from the reference tract
length** (§4.1 pins that denominator against the model's, which is the allele), and
the table above is where it comes from rather than a round number chosen for looking like one: the
strata this model describes well sit at 0.9%, the strata it describes badly at 33.8% and 58.5%, and
there is nothing in between. Ten percent separates them by a factor of three either way, so a stratum
crossing it is unambiguous. **Soft, and honestly so:** three bands of one dataset is a coarse basis
for a boundary, and the number should move if the per-stratum distribution turns out to be
continuous rather than the two clumps this table suggests. **Measurable now**, from the same tooling
that produced the table, at per-stratum grain rather than in three bands.

No repeat-count threshold is hardcoded: thin strata already borrow under §4.3's rule, and a stratum
whose fit is meaningless for this other reason is now visible rather than merely averaged in.

**These are also the strata §4.5 places two requirements on**, and the table above is what makes both
necessary: at 0.091% the direction split and the fall-off have almost no reads behind them however
many loci the stratum holds, and at 0.091% this path's noise model is close enough to the generic
path's that their two substitution rates have to agree.

### 5.1 What makes a locus an STR locus, and where the copy floors go

**An STR locus, on this path, is not a locus that contains a short tandem repeat. It is a locus that
is likely to stutter.** The repeat detector answers the first question and answers it from the
reference; this section answers the second. They are different questions and they have different
answers: a hexamer at three copies is unambiguously a tandem repeat and, by this definition, not an
STR locus.

**ng's code already said this and had no evidence for its numbers; it has them now.** The floors
were `[6, 4, 4, 3, 3, 3]` for periods 1–6, documented as *"the copy number at which a repeat starts
to stutter — below it, the generic SNP/indel caller handles the tract fine and only a stuttering one
needs the STR route"*, with every value marked *"a starting value, soft and swept"*. This section is
the sweep, and **the floors it measured — `[8, 6, 6, 6, 5, 4]` — are what the code applies since
2026-08-10** ([`segment_criteria.rs:444`](../../../../src/ng/region_typing/segment_criteria.rs), with
the two surveys' reasoning at `:402-443`). §5.1.1 is the derivation, and the paragraphs below are the
evidence that got there; where one of them says "ng's `[6, 4, 4, 3, 3, 3]`" it is describing what was
being argued against, not what runs.

**Two things a floor has to clear, and the second is the one that bites.** There has to be enough
stutter to be worth a stratum at all — §5's 0.091% below four repeats against 2.006% at six or more.
But the sharper test is what *kind* of difference the reads show. Below four repeats **58.5% of the
reads that differ from the allele differ by something that is not a whole number of copies**, which
this noise model has no way to express. A thin stratum is recoverable: it borrows from a neighbour
and the provenance records it (§4.3). **A mis-modelled stratum is not.** It returns a confident
slippage rate that is mostly ordinary indel wearing the wrong model, and — worse than being wrong on
its own loci — it enters the monotonicity walk of §4.3 and can drag a neighbour into a merge with
it. So the floor goes where the guard share crosses §5's one-in-ten threshold.

**Below the floor nothing is lost; it is re-described.** Those tracts go to the generic path, where
an ordinary indel is exactly what the model expects. That is the argument for moving them rather
than a cost of moving them.

**Whether a locus stutters is a property of the library, not of the locus — and the reference cannot
know it.** A PCR-amplified library stutters more than a PCR-free one, so a tract not worth the STR
path in one is worth it in the other. The ideal routing is therefore per (locus × library) and
region typing, which sees only the reference, cannot deliver it.

**Settling for the reference approximation is not only a practical compromise — something else in
this design depends on it.** Region typing's output being a pure function of the reference is what
makes it *identical in every sample*, and that is what the STR census rests on:
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §3 selects its loci from
that output precisely because "region typing delimits STR loci from the reference alone, never from
the reads, so its output is identical for every sample". Route per library and two samples no longer
share a locus set, so the census stops holding the same questions for everyone and the cohort
comparisons it exists for become meaningless rather than merely noisier. **The reference
approximation is load-bearing, not just convenient.**

**What this step can still do about it, and it is not nothing: measure whether the floor was right
for this library.** The guard-bucket share is fitted per (read group × stratum), so a library whose
low-repeat strata come back above the threshold is a library whose floor should have been higher —
and §4.4's summary already reports exactly that, per read group, with the worst stratum named. **So
the floor is set from the reference and checked against the data**, and a mismatch is visible rather
than silent. What a run cannot do is act on it; that would be a second pass.

**Which decides how to choose the default, because the two errors are not symmetric.** A floor set
**too low** sends tracts to the STR path that do not stutter in this library. They arrive as strata
with a high guard share, and the per-read-group audit above names them — **the error announces
itself**. A floor set **too high** sends tracts that *do* stutter to the generic path, where they
are called as ordinary indels and nothing anywhere reports that a stuttering tract was mis-routed —
**a silent error**, and silent errors are what this whole step exists to remove.

**Decision: set the default from the most stuttering library available, and let the audit catch the
over-inclusion.** That is the direction where being wrong is visible.

**It also decides which of our two datasets sets the numbers, and it is not the one with the better
truth set.** GIAB's HG002 is prepared with unusual care and sits at the *low*-stutter end, so floors
derived from it would be too high for an ordinary library — the silent direction. The tomato
accession is an ordinary short-read library (`@RG PL:ilumina LB:PRJDA59759_DRR000741`, no PCR-free
marker) and stutters more: 1.97% of its differing reads are not whole-repeat against HG002's 1.37%,
with 28 of 148 strata above the guard threshold against 10 of 132. **So the floors below are set
from tomato and checked against HG002**, which inverts the usual order of preference between those
two datasets and is deliberate.

*The honest limit on both, and it is the one that matters rather than the species:* **each is a
single library**, so neither varies the axis that actually drives stutter. What would settle the
defaults properly is the same measurement across libraries of known preparation — PCR-amplified
against PCR-free — which the tomato archive can supply and this document cannot.

**Deferred, with a home:** per-library re-routing, as a *downstream refinement* rather than a change
to region typing — a locus that region typing sent to the STR path can be handled by the generic
model at calling time without disturbing the locus set the censuses share. **Home:** whichever spec
takes on the calling step's marker routing.

**The floors are a default, not a constant.** `MinCopies` is already a per-period knob
([`segment_criteria.rs:355`](../../../../src/ng/region_typing/segment_criteria.rs)), so what follows
changes a default rather than adding a mechanism. Someone cataloguing repeats rather than genotyping
them wants a lower floor and should set one — and since 2026-08-11 that is exactly what the catalog
does, at floors of its own that sit under every routing floor (§5.3).

**One thing the measurement has to check, because the current defaults appear to be keyed on the
axis §4 rejected.** In copies the floors are `[6, 4, 4, 3, 3, 3]`; in **bases** they are
`[6, 8, 12, 12, 15, 18]` — close to a constant tract length, and the comment's reasoning
(*"shorter motifs stutter more, so periods 4–6 sit at 3"*) is a base-length argument. §4 measured
that the copy axis is the one slippage varies on and the base-length axis is not: at 12–15 repeats
three periods stutter within a factor of two of each other, while at 20–29 bp the same three are an
order of magnitude apart. **If that holds at the low end too, the floors should be far closer to
uniform in copies than they are** — which would raise periods 4–6 substantially. That is a
prediction, and the per-period measurement below is what confirms or refutes it.

**One thing the guard share cannot do, and it changes how it is read.** Its absolute level is not
comparable *between* periods: a hexamer has five ways to differ from the reference by a non-whole
number of copies and a mononucleotide has none, so **the guard share for period 1 is identically
zero at every repeat count** — which is §2's "homopolymers are 100% by arithmetic … and are
therefore no evidence either way", arriving from the other direction. So the guard share says where
a period crosses **its own** curve, and monos need the other criterion: how often a read differs
from the reference length at all.

**The per-period evidence — measured 2026-08-07 over the tomato archive.** 2,457 libraries walked
across chromosome 1 at copy floors one step under ng's, one setting throughout.

**What the table below holds, in one sentence: of the reads whose tract came out a different length
from the reference, how many differed by something that is *not* a whole number of motif copies.**
That is the guard share of §4.1, and it is the sharp criterion because it asks whether the movement
is **the kind this noise model can describe**. A dinucleotide read that came out 1 bp short did not
lose a `CA`; it took an ordinary deletion, and a slippage rate fitted from such reads is
mis-modelled indel however many of them there are. Counted out of every 100 reads that differ from
the reference, so the one-in-ten threshold is the number 10:

| period | at 3 repeats | 4 | 5 | 6 | 7 | floor |
|---|---:|---:|---:|---:|---:|---:|
| 2 — dinucleotide | 66 | 43 | 26 | **9** | 6 | **6** |
| 3 — trinucleotide | 64 | 39 | 10 | 20 | **6** | **7** |
| 4 — tetranucleotide | 76 | 32 | 27 | **5** | 2 | **6** |
| 5 — pentanucleotide | 52 | 13 | **7** | 12 | — | **5** |
| 6 — hexanucleotide | 12 | **4** | — | — | — | **4** |

So a three-copy dinucleotide is a tract where two reads in three that move at all move in a way the
STR path cannot express, and a six-copy one is a tract where fewer than one in ten do. The floor is
where a period crosses ten and stays across.

**Split into its two parts, the criterion says something sharper than the ratio does.** Per 100
reads at a dinucleotide tract:

| dinucleotide | 3 repeats | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| moved by whole copies | 0.19 | 0.72 | 2.58 | **9.07** | 15.30 | 23.17 |
| moved by something else | 0.38 | 0.54 | 0.90 | 0.90 | 0.95 | 0.39 |

**The second line barely moves and the first climbs 122-fold.** Ordinary indel error is a property
of the sequencer rather than of the tract, so it sits at roughly half a read in 100 at every repeat
count — inside a factor of 2.5 across the whole range. Repeat slippage is a property of the tract
and rises steeply with it. **So a period's floor is not where tracts stop taking indels; it is where
whole-copy movement finally outruns a constant indel floor.** At four copies the two are comparable,
0.72 against 0.54, so a slippage rate fitted there is fitting about as much indel as slippage even
though the majority is nominally the right kind. At six copies it is ten to one.

**One thing the whole-copy line is not: a slippage rate.** It counts reads differing from the
**reference** by whole copies, and that is slippage *plus* reads that correctly measured an allele
the sample genuinely carries. Nothing here separates the two, and nothing here needs to — telling
them apart is what §4's fit does by summing over the genotype, and it is why the fit exists. The
criterion above asks only which *kind* of movement a tract produces, which is answerable from the
reference alone.

*(Mononucleotides are absent because the quantity is vacuous for them — every integer is a multiple
of one, so no read can differ by a non-whole number of copies and the share is always zero. They
need the other criterion, how often a read differs from the reference length at all. The walk also
censors every period below the floor it was swept to, so nothing here sees a period-1 tract under 5
repeats.)*

The evidence behind each cell is large: the dinucleotide 3-repeat entry pools 292 million loci and
3.4 billion reads.

**Measured floors: `[6, 6, 7, 6, 5, 4]` against ng's `[6, 4, 4, 3, 3, 3]`.** Every period except
mononucleotides is currently set too low, periods 4 and 5 by three copies. **This confirms the
prediction this section recorded**: keyed in copies rather than bases, the floors come out far
closer to uniform than they are, and periods 4 to 6 rise substantially.

**Period 1 is not moved by this survey.** Its only criterion is how often a read differs from the
reference length at all, and both 5 and 6 repeats clear the one-in-a-hundred rate the survey uses
(1.13% and 2.54%) — a rate that is chosen rather than derived. The floor of 6 keeps the reason it
already had, the Illumina read-artifact onset, which this measurement neither supports nor
contradicts.

**Read the pooled curve, not the per-library one, and the reason is a bias worth recording.** The
survey also reports a floor per library, and those sit systematically **below** the pooled crossing
— pentanucleotides call 4 in three libraries out of four where the pooled curve says 5. A thin
stratum in one library often observes *no* non-whole-repeat read at all, so its guard share is
exactly zero and passes. The per-library floor is therefore optimistic wherever data is thin, which
is precisely the low-repeat strata a floor is placed on. The pooled curve has no such problem.

**Two of the crossings rest on a single stratum and should be read as softer than the others.**
Trinucleotides pass at 5 repeats (0.098), fail at 6 (0.203) and pass again at 7; the floor of 7
follows from requiring the crossing to hold, and a rule accepting the first crossing would say 5.
Pentanucleotides pass at 5 and fail at 6 on 5,333 loci, which is thin. Dinucleotides, trinucleotides
at the low end, and tetranucleotides are not in doubt.

**How far the floor moves between libraries — the axis nothing had varied, and the reason for the
survey.** It moves by two to four repeats. Dinucleotide floors run 5 to 8 across libraries with 7%,
38%, 34% and 19% of them at each; tetranucleotides run 4 to 6. So a single number is a compromise
across real variation rather than a constant of the species.

**What that variation is has *not* been established, and the obvious reading is not available.**
Grouping libraries by their project accounts for 45% of the variance in the dinucleotide floor — but
a project is only a proxy for how a library was prepared, and a poor one in both directions:
libraries within a project need not share a preparation, and a project also bundles the instrument,
the submission batch and the **read length**. Read length is the confound this survey documents and
cannot see, and it would produce exactly this pattern for purely geometric reasons, since a 100 bp
library cannot span the tracts a 150 bp one spans. Separating preparation from read length needs the
join to `rick_sample_manifest.sh` that the survey's own header calls for, and until that is done
"PCR amplification stutters more" remains the motivating hypothesis rather than a finding.

**What is ruled out is depth**, which is the alternative that would have made the whole axis an
artifact. The median dinucleotide floor is flat — 6, 7, 7, 6, 7 — across depth quintiles spanning
2,606 to 55.9 million reads a library, a twenty-thousand-fold range.

**And 28 libraries produced nothing**, each because no read witnessed a locus on chromosome 1 while
region typing found the same 644,194 as everywhere else. They are reported rather than counted as
surveyed, which is the failure this section's instrument was rebuilt to make loud.

**And the instrument built to supply it cannot reach below the floors — measured 2026-08-07, and
this is why the table was pending.** The obvious way to see under a floor is to lower it and
re-type, which is what `examples/ng_str_stutter_by_library.rs` and
`scripts/ng_str_library_survey.sh` were built to do. It does not work, and the reason is structural
rather than a bug in either. *(Both tools survive and both now read the catalog rather than typing
the genome themselves, so "re-type" below means "read the file back at a lower floor". The
structural failure is unchanged — it is a property of what a low floor does to bundling, not of
where the repeats came from — and a second bound now sits under it, §5.3.)*

**`MinCopies` decides two different things and only one of them is the floor.** It is read by
`prefilter`, which runs **before** bundling
([`segment_criteria.rs:648`](../../../../src/ng/region_typing/segment_criteria.rs)), and again by
`classify` after it (`:820`). The first reading decides **what counts as a neighbouring repeat**;
the second decides **what is admitted as a locus**. Folding them into one value was deliberate —
the type's own note calls a swept floor moving both "the whole point" (`:346-353`) — and it is
correct at the defaults, where the two questions have the same answer. Under a sweep they come
apart, because two copies of a mononucleotide is any `AA` and occurs every few bases. Admit those
and every real tract acquires a neighbour inside the bundle threshold, so no tract has clean flanks,
so region typing emits `SsrBundle` — which names no locus at all.

Over a 2 Mb slice of tomato SL4.0ch01, one library at 80×:

| copy floors | STR loci typed | bundles | bundle bp |
|---|---:|---:|---:|
| `[6, 4, 4, 3, 3, 3]` — ng's defaults | 6,237 | 1,943 | 74,289 |
| uniform 4 | 7,434 | 11,720 | 675,372 |
| uniform 3 | 848 | 7,950 | 1,623,636 |
| uniform 2 | **0** | 1 | 1,177,849 |
| `[2, 4, 4, 3, 3, 3]` — period 1 alone | **0** | 18 | 1,599,157 |
| `[6, 2, 4, 3, 3, 3]` — period 2 alone | 1,618 | 13,913 | 1,017,906 |

At a uniform floor of two the whole 2 Mb span comes back as **one bundle covering 1.18 Mb** and
zero loci. That is the failure the archive survey hit on 2,475 tomato CRAMs: every walk succeeded,
and every walk measured nothing.

**Lowering one period at a time does not rescue it, and that is the finding that closes this
route.** The loss is not confined to the period that moved: dropping period 2 alone takes
**period-1 tracts at six copies from 2,678 loci to 225** — 92% gone, at a period whose own floor
never changed, and with the off-reference share moving with it (0.0263 to 0.0210 at six repeats).
So two settings do not give two ends of one curve. They give two different populations of loci, and
the second is a selected subset: the tracts that happened not to acquire a neighbour.

**What follows for §5.1's table: one step down is usable, and the criterion this section places the
floors on survives it.** A sweep of `[5, 3, 3, 3, 3, 3]` — one step under ng's defaults at the three
periods that can move — still costs about half the loci of every stratum it shares with a default
walk. But the two criteria are not equally damaged, and the sharper one holds:

- **The guard share is essentially unmoved.** The same stratum, dinucleotides at 4 repeats, reads
  0.345 in a default walk against 0.344 in a swept one; at a 20 bp radius, 0.431 against 0.408. So
  the criterion that decides periods 2 to 6 can be read off a swept walk directly.
- **The off-reference share is not**, coming back about 30% low, because the surviving tracts are
  the isolated ones and isolated tracts sit in cleaner context. That is the only criterion
  **mononucleotides** have (their guard share being zero by arithmetic, §2), so the period-1 floor
  rests on softer evidence than the rest. The bias understates stutter, which pushes a floor *up* —
  the direction this section's decision rule calls silent.

**And the answer the section turns on is robust to all of it.** Dinucleotides at 3 repeats read a
guard share of **63%** against the one-in-ten threshold, on 2,031 to 3,259 loci, unchanged across
every bundle radius and both copy-floor settings tried. Trinucleotides at 3 repeats read 51%, and
dinucleotides one repeat higher read **0%**. There is a clean step between 3 repeats and 4, and it
is six-fold clear of the threshold, so no correction of the size measured above can reach it.
**Periods 2 and 3 keep their floor of 4.**

**What the 2 Mb slice cannot settle is periods 4, 5 and 6**, which hold 67, 20 and 7 loci at 3
repeats. Those are locus-starved rather than biased, and the archive walk is what fixes them.

**The separation of the two roles is therefore not needed** — a bundling floor held at ng's defaults
while a classification floor moves. It remains the clean way to extend a curve below its floor if a
period ever comes back near the threshold rather than far from it.

**The direction of this section's decision survives intact.** A floor set too low announces itself
through the guard-share audit; a floor set too high is silent. So the default still comes from the
most stuttering library available, and the crossing per period per library — which the archive
survey delivers — is what places it.

### 5.1.1 The floors as set, 2026-08-10 — `[8, 6, 6, 6, 5, 4]`

**Two measurements stand behind them, asking different questions, and where both can answer they
agree.**

- **Where the model applies** — the guard share, over 2,457 libraries (§5.1 above). Crossing one in
  ten puts it at `[·, 6, 7, 6, 5, 4]`.
- **Where stutter grows large enough to matter** — the per-read slippage rate, fitted by summing
  over the genotype, over 181 libraries
  ([`ng_str_stutter_rate.rs`](../../../../examples/ng_str_stutter_rate.rs)). The repeat count at
  which the rate first reaches **5%**:

| tract | libraries | earliest | median | latest |
|---|---:|---:|---:|---:|
| mononucleotide | 165 | 9 | 11 | 13 |
| dinucleotide | 46 | 6 | 7 | 9 |

**The floor tracks the earliest library, not the median**, and that follows from this section's
asymmetry rather than from caution: a stuttering tract sent to the generic path is genotyped by a
model that cannot express what it does, and nothing reports it; a quiet tract sent to the STR path
is genotyped correctly, only more slowly. **So the library axis this survey exists to measure moves
the answer by three to four repeats**, and the floor has to sit under all of it.

**Mononucleotides at 8 rather than the measured 9**, for one repeat of margin: the most-stuttering
library reaches 3.6% at 8 against 1.0% at 7 and 0.9% at 6, so 8 is where the margin stops being
free.

**Trinucleotides at 6 rather than the model-fit 7** — the one place the two criteria disagree, and
the same asymmetry decides it. The worst library reaches 3.4% at 6 repeats and the rate climbs about
2.5-fold per repeat, so 6 catches a library that 7 would miss. **This is a deliberate step below
where the noise model describes the data well**, and it is not free: at 6 repeats about one in five
of the reads that differ do so by a non-whole number of copies, so a slippage rate fitted there is
part ordinary indel. The trade is taken knowingly — a mis-modelled rate on a tract that is genotyped
beats a correct rate on a tract that is not.

**Periods 4 to 6 have no 5% crossing in this data at all** — no library reached it at any repeat
count holding enough loci to fit — so those floors are the guard share's alone. What the rate survey
adds there is only a bound: the most-stuttering library reaches 1.1%, 0.55% and 2.4%.

*What is thin, said plainly.* The dinucleotide range rests on 46 libraries against the
mononucleotide's 165, and the trinucleotide extrapolation on 37 libraries at 6 repeats. The run was
cut short at 181 of 1,400 libraries; more would tighten the tails, and nothing seen so far suggests
the medians would move.

**What adopting the measured floors cost, because it is larger than "changing a default" suggests.**
Against the loci ng emitted at the old floors:

| period | current | measured | STR loci today | kept | re-routed |
|---|---:|---:|---:|---:|---:|
| 1 | 6 | 6 | 272,444,044 | 272,444,044 | 0.0% |
| 2 | 4 | 6 | 38,651,283 | 3,716,800 | 90.4% |
| 3 | 4 | 7 | 5,168,746 | 355,932 | 93.1% |
| 4 | 3 | 6 | 7,976,568 | 83,722 | 99.0% |
| 5 | 3 | 5 | 1,338,183 | 25,771 | 98.1% |
| 6 | 3 | 4 | 523,635 | 58,836 | 88.8% |

**15.2% of ng's STR loci in total, and roughly nine in ten of every non-mononucleotide one.** The
total is small only because mononucleotides are 84% of the loci and do not move. **Nothing is lost
by this** — §5.1's own argument is that those tracts are re-described rather than dropped, and they
go to the generic path where an ordinary indel is what the model expects, which is exactly what
their 30% to 76% guard shares say they are producing.

**The change was taken on 2026-08-10** and `MinCopies::default` is `[8, 6, 6, 6, 5, 4]`
([`segment_criteria.rs:444`](../../../../src/ng/region_typing/segment_criteria.rs)). **Two things
this document says about measurements are now dated by it**, and both matter to anyone rerunning
them: every empirical figure in §2, §3, §5 and §5.1 was taken over the loci the *old* floors
admitted, and so was the per-stratum accounting the research note reports (the entry counts, the
guard-share tallies, the modal-offset distribution). Those numbers stand as records of what was
measured; **a run today sees a different, smaller population of non-mononucleotide loci, and a test
asserting one of them will fail on arrival** rather than catch a regression.

**And the mononucleotide floor of 6 carries a reason that survives this framing rather than being
overridden by it.** It was chosen deliberately over ~9 — "the Illumina read-artifact onset, not the
higher ~9-unit germline-slippage threshold"
([`segment_criteria.rs`](../../../../src/ng/region_typing/segment_criteria.rs)). Under the
definition above that is right: the stutter model is a **noise** model, so a tract that stutters
because of the instrument needs the STR route exactly as much as one that stutters in the germline.
Raising that floor means overturning a considered choice, not filling a gap. What the survey adds is
the other direction: mononucleotides at 5 repeats read 1.34% off-reference against the survey's 1%
criterion, and the swept walk's bias understates that, so 5 is worth testing on the archive.

### 5.2 The bundle radius, and why it moved to 15 bp

**A tract is only nameable as a locus if it has clean sequence either side**, and how much is the
bundle radius (`DEFAULT_BUNDLE_THRESHOLD`). It is a second lever on the same problem as the copy
floors: narrowing it means fewer neighbours spoil a tract, so more tracts survive as loci. Measured
on the same 2 Mb slice at ng's copy floors:

| radius | loci | clusters | bases in clusters |
|---:|---:|---:|---:|
| 30 bp | 6,237 | 1,943 | 74,289 |
| 20 bp | 7,245 | 1,582 | 47,623 |
| **15 bp** | **7,820** | **1,334** | **36,053** |
| 10 bp | 8,467 | 1,045 | 24,906 |

**Decision: 15 bp, changed 2026-08-07 from 30.** The 30 was itself unmeasured, chosen as "more than
enough unique sequence to anchor a short read".

**What had to be checked first, because no real-data measurement can check it.** The radius caps
`flank_bp`, the anchor the STR aligner places a read against, so narrowing it shortens every read's
anchor. A mis-anchored read placed a whole repeat off lands in an offset bucket and reads as
**slippage** — the very parameter this step fits — so a measurement with no truth to compare against
would report the failure as a finding. The answer comes from the synthetic bake-off
([`ng_ssr_synthetic_bakeoff.rs`](../../../../examples/ng_ssr_synthetic_bakeoff.rs)), which builds
each read from a chosen allele and scores the recovered length against it. For the shipped delimiter
across 166 scenarios — including flanks carrying a repeat of their own below every copy floor, which
is what the radius permits and what a marginal locus looks like:

| flank | exact on clean reads | exact on noisy | 1 bp flank indel | composite |
|---:|---:|---:|---:|---:|
| 25 bp | 1.000 | 0.999 | 1.000 | 0.9994 |
| 20 bp | 1.000 | 0.999 | 1.000 | 0.9993 |
| **15 bp** | **1.000** | **0.999** | **1.000** | **0.9992** |
| 10 bp | 1.000 | 0.999 | 1.000 | 0.9984 |
| 6 bp | 1.000 | 0.999 | 1.000 | 0.9985 |
| 4 bp | 1.000 | 0.992 | **0.500** | 0.8699 |
| 2 bp | 1.000 | 0.927 | **0.250** | 0.7854 |

**Nothing moves until 6 bp, and it breaks between 6 and 4** — so 15 sits at two and a half times the
width where the first number budges. What gives out first is tolerating a 1 bp indel *inside* the
anchor, which four bases cannot absorb.

**10 was defensible on the same evidence and not taken.** It buys 8% more loci, and more loci is the
only thing narrowing buys, so it would spend a third of the margin on something that is not scarce.
The dimension the harness holds favourable is base composition: its flanks are deliberately
aperiodic, where real tomato flanks are AT-rich, and 10 bp of AT-rich anchor carries less than 10 bp
of aperiodic anchor. At 15 that shortfall sits inside the margin.

**The radius does not move this section's answer**, which is worth recording because it means the
archive survey's floors do not have to be re-measured when it changes: dinucleotides at 3 repeats
read a guard share of 63%, 62%, 65% and 64% at radii of 30, 25, 20 and 15.

### 5.3 The loci arrive from a file built once, not from a scan run per sample — ADDED 2026-08-11

**Nothing looks for a tandem repeat while a sample is being read any more.** The reference is
scanned once, by a command whose only job is that, and the repeats it found are written to a
catalog file; every run afterwards reads that file
([`repeat_catalog.md`](repeat_catalog.md)). The walk that used to scan as it went is deleted. On
tomato chromosome 1 the same 1,829,315 loci over the same 60,580 regions now take 0.618 s against
2.101 s, and the difference is the scan.

**This changes nothing in this document above §5.** A locus still arrives carrying its reference
bases, its motif and its reads; the stratum key and every read's offset are still computed from
those; the noise model, the accumulator and the four fits never knew where the tract came from.
**Three practical things do change**, and each is a property of a run rather than of the estimator:

- **A run needs the catalog as an input**, beside the alignments and the reference. Every check in
  §10 that touches real data takes one, and a catalog is bound to the reference it was built from —
  it carries the contig table and a content digest, and pairing it with a different reference is
  refused rather than silently served.
- **The file's own floors bound what any run can ask for, and only downward.** It is built at
  `[5, 5, 4, 4, 4, 3]` copies, 15 bp of flank either side and periods 1 to 6, deliberately below
  every floor a caller routes on (§5.1.1's `[8, 6, 6, 6, 5, 4]`). A reader asking beneath any of
  those three is refused **by name**, saying which axis and both numbers. Above them, the routing
  floor, the purity floor, the score floor, the satellite cap and the bundle radius are filters over
  stored fields, so moving one costs a re-read of the file rather than a re-scan of the genome.
- **So §5.1's per-period evidence cannot be re-run as a flag on a run.** Its tables report guard
  shares at 3 repeats for periods 2 to 6, which is under the file's own floors at four of them.
  Reproducing that evidence means building a catalog at lower floors — a build of the reference,
  taken deliberately, not a sweep inside a walk.

**One thing the file makes possible that this path does not use, said here so nobody wires it in by
reflex.** The catalog can count and randomly sample loci per `(period, repeat count)` before a read
is touched, which is the reason it exists at all: the second route to these same parameters keeps
evidence at a **fixed sample** of loci per stratum and needs the stratum enumerated first
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.2). The per-sample route
described here walks every STR locus the sample covers and needs neither. Its stratum populations
are therefore available as a denominator — how many loci the genome holds against how many this
sample witnessed — and that is a diagnostic worth reporting rather than an input to any fit:
**what a stratum's parameters stand on is the loci this sample actually observed**, and the two
floors of §4.3 and §4.5 count those.

---

## 6. A different STR number this path does not fit: how diverse repeat tracts are

A repeat tract mutates orders of magnitude faster than a base does, so the population's diversity at
STR loci is a different quantity from its diversity at ordinary sites — not a correction to it. The
STR path already depends on such a number: `SFS_THETA = 0.01`
([`src/ssr/cohort/freebayes_emit.rs:42`](../../../../src/ssr/cohort/freebayes_emit.rs)), described in
its own comment as *"freebayes' default `-T`"*, which is a **SNP-scale** value applied to repeat
tracts. It is not incidental — the same comment explains that each distinct allele pays a factor of
`θ`, so this constant is what decides how much read evidence a rare STR allele must produce before it
is believed. Too small a value suppresses real STR variation.

**Nothing measures it today**, and measuring it is why the **STR census** exists
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1.1) — accumulated in this
same walk, holding tract lengths rather than base counts. The number itself is computed at the gather
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3), kept distinct from the generic
diversity.

**It is not fitted here** because it is a property of the population rather than of the chemistry,
and because no single sample contains it.

**One thing §4.1's revision changes about it.** The per-stratum genotype frequencies this path now
fits (§4.2) are an allele-length distribution per stratum, in one sample. That is not STR diversity —
it is one individual's allele spectrum, not the population's — but it is the same kind of object, and
the gather is where the two should be reconciled. **Deferred, with a home:**
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3.

---

## 7. Alleles longer than a read — the approach we are not taking, and what it would cost

This path reads an allele's length off reads that contain the whole tract, so an allele longer than a
read is invisible to it. **There is a well-worked alternative that is not, and it is worth recording
why it is deferred rather than rejected**, because the reason is a property of our data and not of the
method.

**How GangSTR does it.** An allele longer than a read produces no containing read at all — the code
makes this a hard zero, not a soft penalty
([`GangSTR/src/enclosing_class.cpp:41-42`](../../../../GangSTR/src/enclosing_class.cpp)). So the
length is inferred from three indirect signals instead
([`read_pair.h:27-36`](../../../../GangSTR/src/read_pair.h)):

| evidence | what it observes | what it says about a long allele |
|---|---|---|
| flanking read | a read that runs into the tract and off the end | a **lower bound** only — nonzero likelihood only above what was seen ([`flanking_class.cpp:59`](../../../../GangSTR/src/flanking_class.cpp)) |
| spanning pair | both mates in the flanks, tract between them | the fragment length shifted by the allele's excess over the reference ([`spanning_class.cpp:122-135`](../../../../GangSTR/src/spanning_class.cpp)); pairs beyond mean + 3 SD are discarded, so it runs out quickly |
| **fully-repetitive read** | a read lying entirely inside the tract, its mate anchored in a flank | **the length measurement itself** — impossible unless the tract is at least one read long ([`frr_class.cpp:38-41`](../../../../GangSTR/src/frr_class.cpp)) |

**The measurement is a count, not an observation.** A read can only fall wholly inside the tract if it
starts within the first `allele_bp − read_length` bases of it, so that excess is a target window whose
size is read off how many reads land in it — a Poisson term
([`frr_class.cpp:187-190`](../../../../GangSTR/src/frr_class.cpp)):

```text
expected fully-repetitive reads  =  coverage / 2 / read_length  ×  (allele_bp − read_length)
```

**This is why such a method needs a parameter pre-pass of its own, and a different one from ours.**
Every term above is an integral over the library's fragment-length distribution, and the count term
needs coverage; GangSTR therefore profiles fragment length, coverage, GC and read length before
genotyping anything ([`main_gangstr.cpp:424`](../../../../GangSTR/src/main_gangstr.cpp), ahead of the
genotyper at `:437`), keeping the fragment-length distribution as a non-parametric 2,000-bin histogram
([`bam_info_extract.cpp:348-372`](../../../../GangSTR/src/bam_info_extract.cpp)). The profile is not an
optimisation for the long-allele case — **it is the measuring instrument**, since the direct evidence
is exactly zero. The code is explicit: profiling aborts below 500 usable fragments
(`bam_info_extract.cpp:321-326`), `--nonuniform` deletes the count term outright
(`likelihood_maximizer.cpp:416`), and the tool documents its quality score as confidence *"in short
allele calls (shorter than read length)"* ([`README.md:247`](../../../../GangSTR/README.md)).

**Why it does not fit our data: the count is Poisson, so precision is bought with coverage.** Below is
that formula evaluated at 150 bp reads and a 3 bp motif — **arithmetic from the code, not a
measurement**, and the count term in isolation (the fragment-length terms add information, so real
precision is somewhat better):

| coverage | expected reads, 100-copy allele | expected reads, 200-copy allele | 1 SD on a 100-copy allele |
|---|---:|---:|---:|
| 30× | 15 | 45 | ± ~13 copies |
| **3× — the tomato cohort** | **1.5** | **4.5** | **± ~41 copies** |

At 30× the count triples when the allele doubles and the method separates them. At 3× the whole
discrimination lives between Poisson means of 1.5 and 4.5, and the modal observation for a 100-copy
allele is **one read or none**. GangSTR also assumes paired-end throughout and uniform coverage for
the count term; neither is safe to assume across the tomato archive.

**Decision: deferred, not rejected — and if it is ever taken up, the tract-length error must be
measured first.** Note that even at 30× the table gives ± ~13 copies on a 100-copy allele. A
count-based length estimate is intrinsically imprecise, and how imprecise is a function of coverage,
motif period, read length and the true length together. **The measurement that would settle it:**
simulate tracts of known length across that grid and plot recovered length against truth — the
synthetic harness ([`synthetic_validation.md`](synthetic_validation.md)) already sweeps depth, period
and tract length, and would need only the long-allele cases added. Until that error is known, "we can
genotype long alleles" is a claim with no number attached to it.

---

## 8. Open questions

1. **How far should the fit be allowed to place an allele from the reference tract length?** — the
   question §4.1's decision opened, and now **half measured**. The *recorded* offset range turns out
   not to matter much — ±1 costs 0.05% of the level with the marginal rule (§4.1) — but the allele
   lengths the fit may place mass on do, because that is what lets an end bucket be attributed to an
   allele rather than to a far slip.
   **The distribution is measured** ([research note](../research/parameter_estimator_experiments_2026-08-06.md)
   §6.8): on HG002, **88.9 loci in 100 sit exactly at the reference length**, ±4 holds 99%, ±12
   holds 99.9% and ±19 holds 99.99%; tomato is tighter, **95.7 in 100 at zero** and ±1 holding 99%,
   which is what a within-species reference should give against a human sample carrying two
   haplotypes' divergence from GRCh38. So a limit of about ±6 covers all but roughly one human locus
   in 200. *(An earlier version of this paragraph wrote ±5 and ±18 for two of those rungs; the
   research note's ±4 and ±19 are the measurement.)*
   **CLOSED, and what the remaining loci cost turns out to be a threshold rather than a slope.**
   A locus outside the support has its reads explained the only way left, as slippage. Measured
   ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.4.1): leaving
   **2.5% of loci** outside costs 0.1% of the slippage level, **7.9%** costs 2.5%, and **19.3%**
   costs **+499% with the direction asymmetry destroyed** — 0.17 becoming 0.47, the same collapse
   the modal origin produced. The cost is nothing and then everything, with the transition between
   about 8% and 19%.
   **Decision: ±6**, which leaves about 0.5% of HG002 loci outside — a fifth of the way to the row
   that is already free. The limit is not a number to tune but a threshold to clear, and nothing in
   the measured distribution comes near it.
2. **Should the fall-off depend on the level?** — OPEN. Read groups that stutter more also decay more
   slowly: at tomato dinucleotides, level against the one-step share ran ρ = −0.69, and that survived
   removing real alleles. If real, the fall-off is a function of level rather than a free parameter,
   which keeps one number per group while fitting better. *Leaning:* model it as a function of level.
   **Soft:** one cohort, one correlation. **Settled by:** the exact-bias harness rather than the
   synthetic validation an earlier version named — generate from a level-dependent fall-off, fit the
   level-independent model, and read the bias off. That is cheaper than simulating reads and it
   answers the question the synthetic run was being asked. Replicating ρ = −0.69 on human read groups
   would confirm the *biology* and cannot be done: HG002's 300× BAM declares a single `@RG`.
3. **Why do hexamers put more mass at −3 than at −1?** — OPEN, and unexplained. Either a real
   long-tract behaviour or an artefact of tract delimitation. It breaks the geometric assumption for
   that period alone. **Settled by:** the synthetic validation
   ([`synthetic_validation.md`](synthetic_validation.md)), which can inject known hexamer alleles and
   see whether the delimiter reproduces them.
4. **Is the low whole-repeat fraction at tetra and penta real?** 62% and 53% in tomato, on 464 and 131
   reads. Both are thin and both are dominated by 3-copy tracts, which sit at the bottom of §5's
   table. *Leaning:* it disappears if region typing raises its copy floors. **Settled by:** re-running
   §2 over loci restricted to higher copy counts — a question this step can answer from the loci it
   already has, without anything upstream changing first.
5. **Can a per-locus stutter model beat the per-stratum one?** — OPEN, and note it is a comparison of
   **models**: both marginalise over the genotype, and how each is maximised is an implementation
   detail that cannot decide it (§1.3). The STR census supplies the cross-sample evidence a per-locus
   model needs ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1.1), which
   is what makes the question askable at all; it is the STR half of
   [`parameter_prepass.md`](parameter_prepass.md) §4.2. *Leaning:* none, but with a specific doubt
   worth recording — **slippage is a property of the chemistry, and a locus is not a chemistry.** The
   per-stratum model pools across loci precisely because the thing being measured does not vary from
   one locus to the next, so a per-locus model may mostly be fitting noise. Mark-2 already resolves
   this by making the per-locus value a *refinement* shrunk toward the stratum prior rather than a
   free fit, which is likely the shape the answer takes.
6. **Free genotype frequencies, or an allele spectrum plus the sample's inbreeding coefficient?** —
   OPEN and **measurable now** (§4.2). Free costs `A(A+1)/2` numbers per stratum where the tied form
   costs `A`, which is what decides whether a thin stratum can be fitted at all. *Leaning:* free
   where the stratum can afford it, tied where it cannot, recorded as provenance either way.
7. **How often does the monotonicity constraint fire on a truly monotone sequence?** — OPEN and
   **measurable now**, but it needs draws rather than the exact method (§4.3): a spurious merge is
   triggered by sampling noise, and the exact method has none.
8. **What does the per-locus entry cost at 300×, and what should the read cap be?** — **CLOSED.**
   The table is not where this step's memory goes: 0.43 entries a locus uncapped at 300×, 0.36 MB
   over 29,811 loci. And the cap is not a correctness limit either: the scoring rule is exactly
   unbiased at every depth to 45 reads a locus (§4.1). What is left is a precision trade and the
   width of the entry's counters, which is an implementation choice
   ([`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) §2.1) rather than a
   question about the design.
9. **Where do the per-period copy floors go?** — **CLOSED, and the change is made** (§5.1.1). The
   archive survey walked 2,457 tomato libraries and a slippage-rate survey walked 181; together
   they put the floors at `[8, 6, 6, 6, 5, 4]`, and that is what the code applies since 2026-08-10.
   It re-routed 15% of ng's STR loci to the generic path, and nine in ten of every
   non-mononucleotide one. **What this leaves behind is not an open question but a dated set of
   measurements**: every empirical figure in this document was taken over the loci the old floors
   admitted (§5.1.1).
   *Two things that stay soft:* the trinucleotide and pentanucleotide crossings each rest on one
   noisy stratum, and the two-to-four-repeat spread between libraries is **not** yet attributable to
   preparation — grouping by project explains 45% of it, but a project bundles preparation with the
   instrument, the batch and the read length, and read length alone would produce the same pattern
   geometrically. Separating them needs the read-length join the survey's header calls for.
10. **How does a fitted level behave when the truth is near zero?** — OPEN and **measurable now**,
    and like §8.7 it needs draws rather than the exact method (§4.5): the level is bounded below by
    zero, so a finite stratum's estimate piles against that boundary and the exact method — which
    has no sampling noise — cannot see it. **The experiment:** draw strata at §5's measured levels
    and at real locus counts, and report how often the level returns zero and how wide the spread
    is. *Leaning:* none on the estimator; the reporting change §4.5 names — the slipped-read count
    emitted beside the level — is worth making whatever the answer.
11. **How closely do this path's substitution rate and the generic path's actually agree at a
    low-slippage stratum?** — OPEN, and **the number that decides §4.5's limit**. That limit is
    currently the generic path's ladder spacing borrowed for want of a measurement, so the first run
    of the comparison sets it. **Blocked only on the generic half existing**, not on any design
    question. *Leaning:* none, and deliberately — the interesting outcome is a persistent gap, which
    is §1.1's unmeasured finding rather than a defect, and pre-committing to a number invites
    reading it as one.
12. **What does a locus that is not the tract its stratum thinks it is do to the four numbers?** —
    OPEN and **measurable now** (§4.6). A collapsed duplication, a mismapped minority, a somatically
    unstable tract or a mis-anchored subset of reads all put reads where no allele pair in the
    support reaches, and this path has two places for the surplus to go — the genotype frequencies
    or the slippage parameters — where the generic path had one. **The experiment:** contaminate a
    stratum in the exact-bias harness with a share `q` of loci drawn from another process and read
    the bias in all four numbers against `q`. *Leaning:* none on the size, and one on the shape of
    any fix — a class the model carries, never a rule that drops loci scoring badly, which is
    threshold-then-count with the long alleles taken first. The one measurement that already exists
    says robustness cannot be assumed: 19.3% of loci placed outside the fit's allele support costs
    **+499% of the level** with the direction split collapsing from 0.17 to 0.47.

---

## 9. Deferred, with a recommended home

- **Long-allele recovery, and the sample profiling it needs** (§7). The two are one item: the profile
  *is* the instrument that measures a long allele. Deferred on coverage, not on merit. **Home:**
  whichever spec takes on long-allele recovery, and its first job is the tract-length error §7 says
  is unmeasured.
- **Per-locus refinement of the stutter parameters.** This step produces the prior; Mark-2 §5 already
  specifies refining it per locus inside the EM. **Home:** stays there.
- **STR diversity** (§6) — accumulated by the STR census here, computed at the gather. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3, together with reconciling it
  against the per-stratum allele spectra §4.2 now fits.
- **A guard-bucket rate that depends on the allele length** (§4.1). The factorisation that makes the
  guard bucket free assumes it does not. **Home:** this document, once the per-stratum guard shares
  of §5's threshold have been measured and show whether it matters.
- **A second class of locus, for tracts the stratum's model does not describe** (§4.6, §8.12).
  Deferred on a measurement rather than on merit: nobody has yet asked what such loci cost, and the
  answer decides whether one parameter is worth adding. **Home:** this document, and if it lands it
  lands the way the generic path's second site class did — its own milestone, on its own evidence,
  between the end-to-end run and the anchors it would move.

---

## 10. How we know it works

**Every check below that reads real alignments also takes a repeat catalog** (§5.3), and the numbers
each quotes were measured before the copy floors moved on 2026-08-10 (§5.1.1). They are records of
what was measured, not thresholds to assert: a run today sees a smaller population of
non-mononucleotide loci, so the shape is what a test should hold — sharded accumulation is an
equality, the guard share separates some strata from others, the table is tens of thousands of
entries rather than millions — with the current numbers recorded when the real accumulator first
walks each cohort.

1. **The fit recovers known stutter behaviour, exactly rather than by simulation.**
   [`../../../../examples/ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs)
   weights every entry by its exact probability under a known truth, so what it reports is bias with
   no sampling noise in it. **Three algebraic checks run before any fit, and each rejects a broken
   scoring rule in one line**: the rule sums to one over the entry space at any parameter values; no
   bucket is charged a negative number of reads; and with the slippage level at zero every locus's
   reads land on its own alleles. **Any change to §4.1's scoring re-runs these first.**
2. **The control that the whole method rests on.** Key the accumulator by locus with the reference
   origin, generate and fit under the same key, and the bias must be **exactly zero** on all four
   numbers. It is: 0.000% on the level, 0.0000 on both shares, and the four starting points agree to
   1.000×. A number other than zero there is the harness's, not the estimator's.
3. **It agrees with truth where truth exists — run 2026-08-08, and it does.** On HG002, the stutter
   parameters fitted by the marginal likelihood must match those measured directly on
   known-homozygous loci. **This is the test production's estimator fails**
   ([`parameter_prepass.md`](parameter_prepass.md) §2.2), and the reason for the whole design. **It
   is also the only check in this document that does not come from the model itself**: every
   recovery test above generates its data from the model it then fits, so a shared misspecification
   cancels and passes.

   Fitted by [`ng_str_stutter_rate.rs`](../../../../examples/ng_str_stutter_rate.rs) over the
   50,000-interval tandem-repeat set at 30×, against the direct measurements in §3 and §5:

   | | fitted, summing over the genotype | measured on known-homozygous loci |
   |---|---:|---:|
   | slippage at 6 or more repeats | **1.851%** | 2.006% |
   | dinucleotide, share of slipped reads gaining | **0.250** | 0.229 |
   | mononucleotide, the same | **0.280** | 0.351 |

   **The direction split is the row that matters**, because it is the one that *inverts* under the
   biased estimator: §2.2 of the shared spec records 0.9× measured over all loci against 3.4× over
   known-homozygous ones — gains appearing more common than losses. The fit returns 3.0× **from all
   loci, with no truth genotypes**, which is the design's central claim demonstrated on real reads
   rather than on constructed worlds.

   *That run predates the copy floors of §5.1.1*, so it was fitted over strata reaching down to 4
   repeats at periods 2 and 3 and 3 repeats at periods 4 to 6. Repeating it today fits a shorter
   ladder of strata; what must survive is the agreement and the direction, not the third decimal.

   *One trap this run walked into first, worth recording because it is a property of the pooling and
   not of the estimator.* A 500-locus floor on which strata are fitted returned 1.077% for the level
   — half the target — because it dropped every stratum beyond 15 repeats, and those carry 14.5% of
   the reads at a mean rate of 6.494%. **A summary pooled over "6 or more repeats" is dominated by
   where its strata stop**, so the floor has to be low enough to keep the tail before the pooled
   number means anything. The per-stratum rates were unaffected.
4. **The non-whole-repeat fraction is reported per stratum** (§5), and a stratum above the
   one-in-ten threshold is distinguishable in the output from a stratum that merely had few loci.
5. **Monotonicity holds on the copy axis and is not silently enforced where it fails.** §4.3's
   merge-and-refit rule assumes slippage rises with repeat count (§4); assert that the fitted
   sequence is monotone after fitting, and that a deliberately non-monotone synthetic input triggers
   the merge rather than being accepted.
6. **The diagnostics aggregate.** A run over a cohort must produce a summary a person reads (§4.4),
   and the test is that a deliberately unfittable stratum — one whose loci are generated with the
   level at zero and the alleles spread — appears in it. A per-stratum record that only a debugger
   would open does not satisfy this.
7. **Sharded accumulation is exact.** The same sample walked in one region and in many must give
   identical tables, which integer entry counts make an equality rather than a tolerance. The read
   cap is seeded from the locus's position for exactly this reason (§4.1).
8. **A barely-stuttering stratum is well behaved, and it is checked in three ways** (§4.5). *The
   level survives thinness:* a stratum generated at 0.091% with the locus count §5's low-repeat
   bands really have returns the level within its own sampling error, and does not return zero more
   often than that error predicts. *The shares do not pretend to:* the same stratum reports its
   direction split and fall-off as borrowed, with the strata they came from named, rather than
   reporting the value its handful of slipped reads landed on. *And the two halves of step 4 agree
   there:* on one real sample, this path's substitution rate at its lowest-slippage strata and the
   generic path's for the same read group must sit within a quarter-Phred of each other, on strata
   at full tract purity. **This is the second check in this document that does not generate its data
   from the model it tests** — the other is §10.3's — because the two rates are fitted by different
   models from different sites, so a misspecification shared with the harness cannot make them
   agree.
9. **The four numbers survive a minority of loci the model does not describe** (§4.6). Fit a stratum
   in which a share `q` of the loci were generated by another process — reads split between two
   tracts several repeats apart, as a duplication the reference does not carry would give — and
   report the bias in the level, the two shares and the substitution rate against `q`. **What this
   check must report rather than pass:** the size of the bias at the shares such loci plausibly
   reach, and **which parameter absorbed it**, because this path can absorb an aberrant locus into
   its fitted allele spectrum instead of into slippage and the two failures need different fixes.
   A run that reports "no bias at q = 0.001" without saying what the emitted allele frequencies did
   has answered half the question.
