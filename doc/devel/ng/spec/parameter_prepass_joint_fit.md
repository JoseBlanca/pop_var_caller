# ng — the joint parameters fit: every parameter estimated once, over every sample at the same loci

*Design spec, 2026-08-10. **No code yet — this settles the design.** **Read this one first.** It says
what the route is, what it produces, why it exists, and how the estimate is computed. Two companion
documents settle the two things it stands on, and each maps to a module someone can build on its own:*

| document | what it settles |
|---|---|
| **this one** | what the route is and what it produces (§1), what having every sample at one locus changes (§2), the estimator (§3–§4), inbreeding (§5), what it cannot reach (§6), and the comparison it exists for (§8) |
| [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) | **which loci every sample keeps evidence at** — the rule, the stratified STR variant and the reference catalog it selects from, the size knobs |
| [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) | **what is recorded at each kept locus** — the two record shapes, the depth ladder, the encoding |

*Types and interfaces: [`../arch/parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md).*

*The shared framing for all of step 4 — the parameters and their grains, why production's numbers are
biased, and the decision to sum over the genotype rather than choose one — is
[`parameter_prepass.md`](parameter_prepass.md), which this assumes.*

***This route is a rival to the per-sample fits of***
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) *and*
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md), ***built beside them and not replacing them***
— [`parameter_prepass.md`](parameter_prepass.md) §4.1 sets up the comparison and these three documents
are what make it runnable. `src/ssr/` and `src/pileup/` are frozen production: everything said about
them here is a record, not a change.*

---

## 1. What this is, and what it is for

**Step 4 has to tell the caller what to expect from this data before any calling starts** — how often
a read shows the wrong base, how often a repeat tract gains or loses a copy, how variable the
population is. There are two ways to get those numbers out of the same reads, and this document is the
second one.

**The two routes differ in one thing: what the genotype is weighted against while a parameter is
fitted.** Both sum over the unknown genotype rather than choosing one
([`parameter_prepass.md`](parameter_prepass.md) §3) — that decision is settled and neither route
reopens it.

- **The per-sample route** walks each sample, folds its loci into histograms, and fits from those. A
  histogram has forgotten which locus each observation came from, so the genotype can only be weighted
  by **one pooled set of genotype frequencies per sample**. It is the only thing available.
- **This route** keeps raw evidence at a bounded set of loci, **the same loci in every sample**, and
  fits everything once when they are all in. Because the loci keep their identity, the genotype can be
  weighted by **that locus's own allele frequency in the cohort** — a quantity the other route cannot
  see at all.

Everything else follows from that one difference. Fitting against a per-locus frequency means every
sample's evidence at a locus enters one likelihood, so the parameters fit cannot be split by sample, so it runs
once — after every sample has been walked and before any calling begins.

**What this route reads, in one line: the records of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md), at the kept loci, and
nothing else.** One entry per kept position holding that position's allele counts and its binned
depth, per read group; one entry per kept STR locus holding its offsets, guard and differences. It
accumulates no summary of its own — a summary has forgotten which locus it saw, which is the property
§2 exists to exploit. **The one object of the other route that survives beside it is not read by
anything here**: the autozygosity `F_autozygosity` needs runs of homozygosity and scattered loci cannot give them
at any budget, so that estimator stays in the genome walk whatever the comparison decides (§6). The two
routes' objects are listed side by side in
[`parameter_prepass.md`](parameter_prepass.md) §5.1, which is where that comparison belongs.

### 1.1 What it produces

**Every parameter [`parameter_prepass.md`](parameter_prepass.md) §1 lists, and one of them is a
different quantity under this route.** The last six rows are already the cohort gather's
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)) and this document changes only how one
of them is reached; the first eight are what it adds.

| parameter | per-sample route fits it from | this route fits it from | same quantity? |
|---|---|---|---|
| per-base error rate, generic | the read-group histogram | the generic kept loci, all samples | yes |
| per-base error rate, STR | the STR table's composition counts | the STR kept loci's, likewise | yes |
| how often an STR read slips | the STR table, per stratum | the same strata, weighted **per locus** (§4) | yes |
| which way it slips | the same | the same | yes |
| how far it slips | the same | the same | yes |
| observed heterozygosity `Hobs` | the windowed histogram, summed | the kept loci, as a sum of genotype posteriors (§3.2) | yes |
| homozygous-non-reference rate `π_hom_alt` | the same | the same | yes |
| **inbreeding** — `F_autozygosity` there, `F_hom_excess` here | runs of homozygosity, per sample | **homozygote excess** against the panel (§5) | **no — two quantities, §5** |
| the cohort's diversity `Hexp` | — | the fitted frequency density **directly**, with no division by `1 − F` (§5.3) | yes, and better conditioned |
| STR diversity | — | the STR kept loci, **reweighted by stratum** ([loci](parameter_prepass_joint_loci.md) §3.3) | yes, once reweighted |
| the frequency spectrum | — | **derived from this parameters fit's own parameter**, the frequency density (§2.1.2) | yes |
| **contamination `α`** | — **it cannot** ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2) | the generic kept loci, against individual-specific allele frequencies this parameters fit derives (§3.4) | it is only produced here |
| relatedness, read-group grouping | — | unchanged: [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6–§7 | yes |

**One row answers *no*, and it is inbreeding** — `F_autozygosity` and `F_hom_excess` are two different quantities, which §5 is about. **One row has
no comparison to make at all**: contamination is produced by this route and not by the other, because
identifying it needs the locus and the allele and both histograms have kept neither
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2). §3.4 is the estimator.

### 1.2 Goals

1. Produce every parameter in that table from the kept loci alone, so that both routes can be run on
   the same data and **their estimates compared** — which is what this route exists for, and the
   reason it is built rather than argued about.
2. **Use what the per-sample route cannot reach**: the locus's own allele frequency (§2.1), and the
   fact that a badly-behaved locus is badly behaved in every sample (§2.2).

**Non-goals.**

- **Replacing the per-sample fits.** Both are built and both run. §8 is the comparison; nothing is
  deleted before it.
- **Designing the caller's priors**, and **calling anything** — inherited from
  [`parameter_prepass.md`](parameter_prepass.md) §1.2.

**It does not:**

- estimate anything **local** — a rate per window, a run of homozygosity, a haplotype (§6);
- change the read groups' role: chemistry is still fitted per read group and biology per sample
  ([`parameter_prepass.md`](parameter_prepass.md) §1.1);
- assume diploidy. The likelihood is written for any ploidy; what is diploid-only is the population
  genetics on top (§10).

---

## 2. What having every sample at one locus changes

### 2.1 The genotype is weighted against the locus's own allele frequency — and that frequency is summed over too

**The tempting version of this design is wrong, and it is wrong for the reason
[`parameter_prepass.md`](parameter_prepass.md) §3 already gives about genotypes.** That section's
table has three things one can do with an unknown: choose it, maximise over it, or marginalise it
away, and the middle one fails because each new locus brings its own new parameter, so the bias does
not shrink as data accumulates (the incidental-parameters problem, Neyman & Scott 1948). **A
per-locus allele frequency is exactly such a parameter.** Two million of them, each backed by fifty
samples at three reads, maximised alongside a handful of noise rates, would reproduce that failure one
level up.

**Decision: two nested sums, and only the outer objects are free parameters.**

1. **Inside a locus, sum over each sample's genotype**, weighted by what the locus's allele frequency
   implies for that sample — this is [`parameter_prepass.md`](parameter_prepass.md) §3's likelihood,
   with the pooled genotype frequencies replaced by the locus's.
2. **Sum over the locus's allele frequency itself**, weighted by how common that frequency is across
   the cohort — that weighting is the **frequency density**, and it is fitted.

So the free parameters are: the noise rates (per read group, and per stratum on the STR path), the
four numbers that describe the frequency density (§2.1.2), and two numbers per sample (§5, §3.4).
**None of them grows with the number of loci.** A locus contributes evidence and holds no parameter
of its own.

*In the standard vocabulary the per-locus frequency is a latent variable given an estimated prior —
empirical Bayes.*

**What this answers.** [`parameter_prepass.md`](parameter_prepass.md) §4.2 asks whether a prior fitted
against per-locus weights beats one fitted against pooled weights, and records that the question is
askable only where the loci keep their identity. This route *is* that arm of the question, on both
paths at once.

#### 2.1.1 A count in the panel or a frequency in the population — and one of them cannot carry inbreeding

**There are two ways to be uncertain about a locus's allele frequency, they are different models
rather than two notations, and an earlier draft of this document wrote one and named the other.**

- **A count in the panel.** The panel's `2N` chromosomes carry `c` copies of the allele; `c` is
  unknown and drawn from a **spectrum of counts**. Given `c` the samples' genotypes are *not*
  independent — they have to add up to `c` — so the sum over genotype configurations is a
  convolution across samples. This is ANGSD's `realSFS` (Nielsen et al. 2012; *Bioinformatics* 2015),
  and its virtue is that the conditional distribution of the genotypes given `c` **contains no
  frequency at all**, so the spectrum is estimated without committing to any population.
- **A frequency in the population.** The population the panel came from carries the allele at
  frequency `f`; `f` is unknown and drawn from a **density**. Given `f` the samples' genotypes *are*
  independent, and the panel's own count is a consequence rather than a parameter.

**Decision: the frequency, and inbreeding is what forces it.** The cancellation that makes the count
form work needs each individual's genotype to be a pair of independent draws at `f`, because that is
what makes a sample's weight `f^j (1−f)^(P−j)` times a combinatorial constant — and only then does
`f` divide out of the conditional. **An inbreeding coefficient breaks that factorisation**: under
`F_hom_excess` a diploid heterozygote's weight is `2f(1−f)(1−F)` and each homozygote's carries an
extra `F·f(1−f)`, which is not of that form, so `f` survives in the conditional and the convolution
is no longer free of it. That is why the estimators which fit inbreeding from genotype likelihoods
work with a per-site frequency rather than a per-panel count (ngsF; Vieira et al., *Genome Research*
2013). §5 makes one inbreeding coefficient per sample a required output of this route, and §3.4 adds
one contamination fraction per sample which needs to know **which allele** the population carries and
not only how many copies of it there are. Neither is expressible in the count form.

**What is *not* the reason, and an earlier version of this document said it was.** Nothing is thrown
away by treating the samples as independent given `f`. The correlation between two samples at a locus
is *induced* by the shared unknown frequency, and integrating over it keeps all of it — the frequency
form is an exact likelihood under its own generative story, not an approximation to the count form.
The two place different priors on the same latent, and the count form is the more flexible of the
pair: every frequency density induces a spectrum of counts, and not every spectrum of counts comes
from a frequency density.

#### 2.1.2 The density is fitted with four numbers, not with one weight per allele count

**The real hazard in the frequency form, and it is a hazard the count form does not have.** A weight
per allele count is a distribution over an almost-observable quantity. A weight per *frequency* on
the same `2N + 1` grid is not: a frequency reaches the data only through the genotypes it draws, so
recovering its density means undoing a binomial blur. That is a deconvolution, and an unregularised
maximum-likelihood deconvolution onto a grid as fine as the sample size collapses onto spikes — the
classic behaviour of a nonparametric mixing distribution, whose maximiser is discrete with at most as
many support points as there are distinct data patterns (Lindsay 1983). Fifty samples at three reads
would not support a hundred and one free weights.

**Decision: a shape with four free numbers, and the grid is quadrature rather than parameters.**

```text
π(f)  =  p_invariant · [f = 0]
       + p_fixed_alt · [f = 1]
       + (1 − p_invariant − p_fixed_alt) · Beta(f; a, b)
```

- **`p_invariant`** — the share of positions where the population carries only the reference base.
  It is most of the genome, and giving it a mass of its own rather than leaving the Beta to imitate
  it is what keeps the Beta describing the sites that actually vary.
- **`p_fixed_alt`** — the share where the population carries only a non-reference base, which is the
  reference accession's own private alleles. On a crop reference that is not a rounding term, and it
  is where `π_hom_alt` (§3.2) mostly comes from.
- **`a`, `b`** — the Beta's shape over the segregating sites. `a < 1` reproduces the rare-allele
  pile-up a neutral population has (`π(f) ∝ 1/f`), and letting it be fitted rather than assumed is
  what allows for a bottlenecked, selfing crop being nothing like neutral.

**The integral over `f` is a fixed Gauss–Jacobi quadrature over `(0, 1)`**, whose node count is a
numerical accuracy choice and **not a count of parameters** — doubling the nodes costs time and adds
no freedom. The two masses are exact terms beside it.

**What is emitted.** A consumer asking for "the frequency spectrum" usually means counts, so both are
emitted and they are named apart: the fitted density `π` — four numbers — and the **panel allele-count
spectrum it implies**, `spectrum(c) = ∫ π(f) · P(c copies among 2N | f, the fitted F's) df`, which is
the object [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4 describes and the one a
caller's prior is built from.

**The alternative is one weight per allele count, fitted in the count form with inbreeding held at
zero, and it is not adopted.** It would be better conditioned and it cannot produce two of this
route's required outputs. *Open, and §12 measures it:* whether the four-number shape is flexible
enough for the tomato panel, checked by fitting a drawn cohort whose true density is deliberately not
a Beta — a two-subpopulation panel, whose frequency density is bimodal — and reporting what the
fitted `Hexp` and `spectrum(c)` cost. **Settled by:** §12.9.

### 2.2 A noisy locus is noisy in every sample, and no per-sample marginal can see that

**This is the strongest statistical argument for the route, and it rests on a measurement made on
2026-08-10.** The generic path's noise model was one substitution rate per read group until
[`../research/noise_model_overdispersion_2026-08-10.md`](../research/noise_model_overdispersion_2026-08-10.md)
measured its tail on HG002: at the 550,976 loci where the GIAB benchmark records no variant of any
kind, **818 carry three or more alternative reads where one rate predicts 29**. The three-genotype
mixture has exactly one class that can absorb the surplus, and fitted heterozygosity came out **1.41
times the benchmark's count** — 776 heterozygous sites where the truth has 550.

The fix adopted there is a second class of site: a locus is *clean* with probability `1 − w` and
*noisy* with probability `w`, and at 30× that is about **one locus in 110 disagreeing with the
reference at 5% instead of 0.19%**. It cuts the heterozygosity excess from 1.41× to 1.09×.

**What it cannot do is say *which* loci they are**, because a histogram has forgotten. So `w` is a
mixture weight applied blindly, and every locus pays a share of it.

**Under this route the locus keeps its identity, and mismapping is a property of the locus rather
than of the sample.** A collapsed paralog raises the alternative-read fraction in *every* sample at
*that* position; a genuine heterozygote raises it in the samples that carry the allele and not in the
others. Fifty samples at one locus separate those two patterns; one sample at that locus does not.

**Decision: carry the same two site classes, and let the class be a per-locus latent variable rather
than a blind mixture weight.** `w` and `ε_noisy` stay cohort-level free parameters, fitted as they are
today; what changes is that each locus's posterior probability of being noisy is computed from all
fifty samples' evidence at it. *Soft, and it is the headline thing to measure:* whether that closes
the residual 1.09× is unknown, and §8 is where it gets tested.

**And the two classes are not enough.** `ε_noisy` lives on the error-rate ladder, which stops at 10%,
while the loci that ask for a class of their own sit at about **half** — and
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1 is right to refuse to widen the
ladder for them: as a noisy rate approaches a half, a noisy site and a heterozygous site are the same
distribution, so a class stretched that far would take real heterozygotes with it. **What those loci
are is a duplication the reference does not carry**, collecting two copies' reads at one position, so
every position where the copies differ shows alternative reads at every depth. Two of five real
alignments ask for such a population — tomato SRR7279482 and SRR7279483, at **0.42% and 0.49% of
sites** — and the model refuses, leaving the only class it has for them: *heterozygous*.

**How many loci that is — measured, and it is thirty times smaller than the fitted share implied.**
The estimate here was that the population is about **8,400 of two million** kept positions if the
fitted share is read directly, or about **1,700** if the share is the alternative-read mass a 10%
rate had to carry. [`../reports/duplicated_locus_probe_2026-08-12.md`](../reports/duplicated_locus_probe_2026-08-12.md)
walked eight tomato alignments and found **150 to 590 of two million**, and it says where the factor
of thirty went. Two quantities were being conflated:

- **0.6% to 3.2% of positions sit in a window carrying about twice the sample's normal coverage.**
  That is the same order as the fitted share, and it is what the mixture is picking up.
- **Only 0.3% to 1.2% of those positions read near half**, because a duplication is silent wherever
  its two copies agree — which is 99% of their length, that 1% being the divergence between the two
  copies rather than anything the caller sets.

The artefact is the product of the two, and the product is small. **Against the same sample's
near-half positions in ordinary-coverage windows — 668 per two million on SRR7279482, which is the
right order for that accession's heterozygosity — the artefact is about a third**, where this
paragraph previously had it six to thirty times larger.

**The class is still worth carrying and it is no longer the largest term.** It is concentrated 24
times over where coverage says it should be, it has the sign the model cannot otherwise produce, and
the parameters fit's only alternative for those positions is to call every sample heterozygous at a
mid-frequency variant. But the 6,000 error positions and 2,500 noisy-class positions
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2 counts are each an order of
magnitude larger than it, and that document is corrected to say so.

**Refusing costs more here than there.** At tomato's
three reads, calling such a locus homozygous-reference under a 10% rate costs about 180 nats across
fifty samples while calling every sample heterozygous at an allele frequency of one half costs about
143 — so the maximum takes the second, and the locus enters the **frequency spectrum** as a
mid-frequency variant with every sample heterozygous. A real half-frequency variant produces that
pattern with probability 0.5⁵⁰. Three things move together as a result, and none of them is a
parameter the locus could absorb locally: the spectrum gains mass it should not have, `Hexp` is read
straight off that spectrum (§5.3), and `Hobs` is a mean of genotype posteriors that are heterozygous
at every one of those loci. The homozygote excess cannot vent it either, being bounded below by zero
where these loci demand the opposite sign.

**Decision, revised 2026-08-14: a third class of site, and what tells it from a heterozygote is not a
rule but two terms in its own likelihood.** The class is a component of the mixture like the other two:
every position gets a posterior over the three, and nothing assigns a position to one. What the
duplicated component predicts, and what its likelihood therefore multiplies at every position, is

- **a read composition near a half** — which alone cannot separate it from a heterozygote, and never
  could; and
- **a depth about twice that sample's own median**, two copies' reads having landed there.

Beside them, the cohort supplies the third piece: a duplication's carriers are all at a half and **none
is homozygous for the non-reference allele**, where a real variant at the same frequency puts some
samples there. That is not a separate mechanism either — it is what the likelihood of the whole
cohort's genotypes at one locus already says.

**Two measured floors say where each term has power, and neither is a threshold anywhere in the code.**
With three samples the missing homozygote means little and with fifty it is nearly decisive: the
inbreeding coefficient comes back 0.5807 against a truth of 0.5942, where a fit with no third class
returns 0.4471 and heterozygosity 50.6% high. At three reads a position a depth of 3 against 6 barely
separates, and at 25 the position's own depth reaches 10.8-fold enrichment with 2.8 wrong calls in 100
against 44 in 100 at 2.5 reads. **The likelihood consults neither number**; it multiplies the terms,
and a flat term contributes nothing.

**Where the floors are used is the output.** Below about twenty-five samples *and* about twenty-five
reads a position neither term carries information, so the fitted weight is not identified and is
emitted as **not identified** rather than as a number — the treatment §6.1 already gives the homozygote
excess at one sample. A fitted zero there is what `CLAUDE.md`'s range principle forbids while *absent*
is what it allows.

**The per-sample coverage-by-window summary is removed** (owner, 2026-08-14;
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4). It discriminated better
than the position's own depth at the shallow end — 14-fold against 1.14 at 2.5 reads a position, read
at the width that depth requires — and it cost a genome-wide accumulator, a GC curve fitted per sample,
denominators from the reference and the analysed regions, and a second full pass over every pileup in a
two-phase run. **What it bought was a corner where calling is marginal anyway**: above twenty-five
samples the cohort pattern does the work for nothing, and above twenty-five reads a position the
position's own depth does. It also never worked end to end — as specified it could not apply its own GC
correction — and the one real-data test with a truth set returned a class weight of zero.

**What made the second discriminator viable is a change to the ladder, not to the model.** At twenty
bins topping out at 124 reads, a doubled position stopped sitting above an ordinary one from 76 reads a
position and was written identically from 98 — inside the range this caller commits to. The ladder now
runs to about 1,500 in the same five bits, and the depth recorded is the position's true depth rather
than the subsampled one
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §2.2).

The class's weight is fitted as `w` is, and a site the run believes is in the class is drawn from a
component at about half alternative reads.

**Its grain is the (locus, sample) pair, where the other two classes' is the locus, and that
difference is not a detail.** A collapsed paralog is a property of the reference and of the aligner,
so it mismaps in *every* sample, and *clean* against *noisy* is rightly a property of the locus alone.
A duplication carried by an individual is a property of **that individual** — it is what a copy-number
variant is, and it is why the discriminator below is *this sample's* window coverage rather than the
locus's. So **class membership is per sample and only the class's weight and its alternative-read
fraction are cohort-level.** Writing it as a third per-locus class beside the other two would make
every sample share one sample's amplification.

**Measured, and both components are there with the per-sample one the larger.** Eight tomato samples
on one window grid hold 84 windows that at least one of them reads near two copies; **40 of the 84
are read that way by exactly one sample and 11 by seven or eight**
([`../reports/duplicated_locus_probe_2026-08-12.md`](../reports/duplicated_locus_probe_2026-08-12.md)
§5). The threshold-free form is sharper: at SRR7279482's ten such windows every other sample also
reads between 1.36 and 2.10 — the reference's own collapse — while at SRR7279501's fifty-two,
SRR7279481 and SRR7279484 read 2.20 and 2.27 and SRR7279482 and SRR7279540 read 0.82 and 0.76, which
is one copy. **That is copy number segregating in the panel, and a per-locus class would force one
accession's amplification onto samples that do not carry it.**

**It must not be conditioned on one position's alternative-read fraction**, which is the constraint
that has survived every revision: a duplicated position and a heterozygous one both read about half
alternative. What changed on 2026-08-14 is what it *is* conditioned on. A per-sample coverage-by-window
summary was specified here and is now removed; the two discriminators are the cohort's genotype
composition and the position's own depth, and the decision above says which applies where.

**What the class is worth, and what it costs where there is nothing to find — MEASURED 2026-08-13
inside the estimator itself**
([`../reports/contamination_floor_and_duplicated_class_2026-08-13.md`](../reports/contamination_floor_and_duplicated_class_2026-08-13.md)).
With the class **off**, a drawn panel carrying duplications returns heterozygosity **60.8% above** the
truth and a homozygote excess of 0.4209 against a drawn 0.600. With it **on**, −1.2% and 0.5948. On a
panel with **no duplications at all** the class is not free: it takes heterozygosity 0.5% low at fifty
samples and **6.4% low at ten**. That last figure is the one a small cohort needs in order to decide
to switch it off, and it is another face of the same threshold — below about twenty-five samples the
cohort pattern cannot tell the class from a real variant, so it costs more than it recovers.

**And it does not explain the human benchmark trio's heterozygosity excess.** The trio comes back 1.23
to 1.28 times its benchmark VCF's count; with the class fitted and a coverage summary supplied, the
three rates are unchanged to three decimals and **the class's weight is zero**. That excess is
something else, and this route has now ruled out the explanation that looked most likely.

**Two consequences for what the class may be used for.** Its fitted weight **must not be emitted as a
measurement of how much duplication a sample carries** — both discriminators recover it about twice
too large while sorting the positions correctly. And the class's *grain* is unchanged and now has a
second reason: it is the (position, sample) pair, because the coverage reading is per sample and works
at any panel size while the cohort pattern is per position and fails on exactly the duplications a
single accession carries.

**Measured, 2026-08-12, and the class is kept**
([`../reports/duplicated_locus_probe_2026-08-12.md`](../reports/duplicated_locus_probe_2026-08-12.md)).
On tomato SRR7279482 at 25× depth, **1 position in 8,600** is both in a window carrying about twice
the sample's normal coverage and reading between 35% and 65% alternative. Inside those windows the
near-half rate is **1.26%** against **0.033%** in ordinary-coverage windows, a factor of 38 at matched
read depth, and the joint cell holds **24.8 times** what independence between the two would predict.
The alternative-fraction distribution inside the two-copy band carries a bump centred on a half that
the one-copy band does not have — 44 times the mass in 0.4–0.5. All eight samples walked, from 2.5× to
28.7×, separate the same way.

**One constraint the measurement adds, and it belongs to the summary rather than to the class**: the
window has to collect about **12,000 aligned bases** before its mean depth tells one copy from two.
At 25× a 500 bp window does; at 3.6× it does not, and the enrichment falls to 1.3 — no separation at
all — until the window reaches 5 kb.
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4 carries what that means
for the stored grid.

*This class was first asked for by the histogram route and is recorded here because this is the route
being planned; that route is implemented and is not being changed for it.*

---

## 3. The generic path

### 3.1 What is fitted

Free parameters, all cohort-level except the last:

| parameter | grain | count |
|---|---|---|
| clean error rate `ε_clean`, noisy error rate `ε_noisy`, noisy-locus fraction `w` | read group | 3 per read group |
| the frequency density `π` — `p_invariant`, `p_fixed_alt`, `a`, `b` (§2.1.2) | cohort | 4 |
| homozygote excess `F_hom_excess` | sample | 1 per sample (§5) |
| contamination `α` | sample | 1 per sample (§3.4) |

**A locus's likelihood.** The reference base is a property of the position; **which non-reference base
segregates is not known and is summed over**, for the reason §3.1.1 gives.

```text
                                   ⎧                                                           ⎫
P(locus) =  Σ    class_weight   ·  ⎪  p_invariant  ·   Π    P(reads | every copy reference)    ⎪
           class                   ⎪                 sample                                    ⎪
                                   ⎪                                                           ⎪
        (class ∈ {clean, noisy})   ⎪ + p_fixed_alt  ·  ⅓ Σ    Π    P(reads | every copy alt)   ⎪
                                   ⎪                     alt sample                            ⎪
                                   ⎪                                                           ⎪
                                   ⎪ + p_segregating · ⅓ Σ   ∫ Beta(f;a,b) ·                   ⎪
                                   ⎪                      alt                                  ⎪
                                   ⎪                        Π    Σ  genotype_freq(j | f, F)    ⎪
                                   ⎪                      sample j    · P(reads | j, ε_class)  ⎪
                                   ⎩                        df                                 ⎭
```

Read outward: a sample's reads at the locus are scored under each genotype `j` and added, weighted by
what a population frequency of `f` and that sample's `F_hom_excess` imply about how common that
genotype is; the samples multiply because **given the frequency they are independent**; the frequency
is integrated over the fitted density (§2.1.2); which non-reference base is the segregating one is
summed over with an equal prior over the three; and the locus's error rate is drawn from one of two
classes. **The innermost sum is [`parameter_prepass.md`](parameter_prepass.md) §3's likelihood
unchanged** — the same `p_j`, the same `ε/3`, the same ploidy-generic loop bound.

### 3.1.1 Which base is the alternative one is summed over, never chosen from the data

[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2 gives the gather the job of deciding
each site's allele set: *"it unions what the samples observed, works out which allele is major"*. That
rule was written for an estimator that reads counts, and **importing it here would bias the rare end
of the density, which is the end everything downstream reads.**

**Why.** At a position where the population carries only the reference base, the non-reference reads
are errors, spread over the three other bases. Choosing the observed largest of the three and scoring
it as one allele's evidence is conditioning on a maximum: with 150 read observations across a
fifty-sample cohort at `ε = 0.002` a site expects about 0.3 error reads, and the sites that show any
show them as *the winner of three draws*. The excess is small per site and there are two million
sites, and it lands entirely on the classes `p_invariant` has to be told apart from — the singleton
and near-singleton frequencies that `Beta(a, b)` describes, that contamination is measured from
(§3.4) and that `Hobs` is a sum over (§3.2).

**Decision: sum over the three, with an equal prior.** It costs a factor of three on the segregating
term alone — the invariant term, which is nearly every locus, carries no allele at all and is
evaluated once. Nothing in the record changes: the four allele counts are already stored per position
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §2.1), and it is the parameters fit
that stops picking one of them.

*Two things this does not do.* It does not model **two** segregating non-reference alleles at one
position; a triallelic site is scored under whichever of the three the sum favours, and the residue
falls into the noisy class. And it does not remove the same hazard from the STR path, where the
alleles are lengths and the candidate set is bounded differently (§4.2).

### 3.2 What is derived rather than fitted

**A sample's heterozygosity and its homozygous-non-reference rate are not free parameters here, and
that is a real difference between the routes.** Once the parameters fit has converged, every sample has a
posterior over its genotype at every kept locus, so:

```text
Hobs(sample)       =  mean over kept loci of  P(genotype is heterozygous | that locus's reads, fitted parameters)
π_hom_alt(sample)  =  mean over kept loci of  P(every copy non-reference | the same)
```

**These are rates over the kept loci, and they estimate the genome-wide rates because the loci were
chosen by a rule that never looks at the data**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §1.2). A subsample drawn without
reference to what is at the position is unbiased for a rate; that is the property the whole route
rests on, and §12.2 checks it.

**A posterior mean is a rate only where there is data, and at three reads a fifth of the kept loci
have none worth the name.** A locus at which a sample has no read has a posterior equal to its prior,
and that prior is `Hexp · (1 − F_hom_excess)` — the model's own prediction of that sample's
heterozygosity. Such a locus therefore contributes the answer rather than evidence for it. At
tomato's three reads a site, Poisson arithmetic puts about **5 kept loci in 100 with no read at all
and 15 in 100 with exactly one**, and a single read barely moves a genotype posterior. **So `Hobs` is
in part a measurement and in part the model repeating itself, and how much of each depends on the
sample's depth.**

**What follows, and it is a reporting requirement rather than a fix.** Emit beside each sample's
`Hobs` and `π_hom_alt` **how many kept loci carried at least one read and at least two**, so that a
consumer can see how much of the number is data. Nothing about the estimator changes: dropping the
empty loci would trade a self-consistent estimate for a truncated one, and the loci are the same in
every sample, so a comparison between samples is unaffected either way. What must not happen is the
number being read as a count of observed heterozygotes when a fifth of its support saw nothing.

**And those empty loci are where a quarter of one known failure lands — MEASURED 2026-08-13.** When
the third class of §2.2 is missing, about a quarter of the resulting excess in `Hobs` comes from
positions with **no read at all**: the posterior there is the prior, the prior is read off a frequency
density the same missing class has already inflated, and the inflation is charged a second time. So
the two counts this paragraph requires are not only a caveat on how much of `Hobs` is data — they are
the size of the surface on which a modelling error compounds.
*This is also why §5.2's opening of the circularity is only half an opening — the same paragraph
there says which half.*

**Trap: the per-sample route fits these two as free parameters and this one does not.** They are
therefore not two estimates produced the same way, and a disagreement between the routes is not
automatically an error in either. It has a specific likely cause: under this route the two rates are
constrained by the spectrum, so a sample whose real genotype distribution departs from what the
panel's spectrum and its own `F_hom_excess` predict cannot express that departure. §8 says how to tell that
apart from noise.

**On a nearly homozygous sample they are a residual under a much larger background, and that is not a
precision problem.** Measured on a proportion's own scale the estimate gets *tighter* as heterozygosity
falls — `√(p/n)` shrinks with `p` — so a low-heterozygosity sample needs no more loci
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2). What changes is the
mixture: at tomato's least heterozygous sample, **at least 97 in every 100 positions carrying an
alternative read are sequencing error or a noisy locus** rather than a real heterozygote, so `Hobs` is
what survives subtracting the background, and a few percent of error in the background's fitted weight
is the whole of the answer. More loci do not change that ratio.

**This is where the route has something to offer, which is why it belongs in this document too.** The
background is fitted from every sample and every locus at once rather than from the one inbred sample,
and §2.2's per-locus identification names the offending loci instead of averaging over them.

### 3.3 How the maximum is found

**Alternate, as the per-sample route already does.** [`parameter_prepass_generic.md`](parameter_prepass_generic.md)
§5.1 settled this for two coupled fits and measured it: starting every error rate at three times the
truth and every frequency at half, the loop converged to the truth in all 25 worlds tried. The same
shape applies here, with three blocks instead of two:

1. hold the frequency density and each sample's `F_hom_excess`, and fit each read group's three noise
   numbers;
2. hold the noise numbers, and climb to the density's four numbers (§2.1.2);
3. hold both, and fit each sample's `F_hom_excess`;
4. repeat until the fitted values stop moving.

**Two of those blocks can trade against one another, and the measurement has to say by how much.**
The density's shape and the per-sample inbreeding both control how often a genotype comes out
heterozygous. They are separately identified in principle — the density is pinned by the distribution
*across* loci and `F_hom_excess` by which samples are heterozygous *within* a locus of given
frequency — but at three reads that separation is thin, and at one sample it does not exist at all
(§6.1). §12.9 profiles the likelihood in the two together and reports the correlation rather than
assuming it away.

**A flat scan over the noise parameters is not affordable here, and that is the one place the
procedure genuinely departs from [`parameter_prepass.md`](parameter_prepass.md) §3.1.** That section
prices the generic scan at 161 steps over a few hundred binned cells. Here one score is a pass over
two million loci × fifty samples, so 161 of them per read group per outer iteration is not a scan over
a table any more — it is a pass over the data. **Decision: climb rather than scan, from several
starting points**, which is the choice the STR path already made for its own reasons
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2). The starting points must span the
separation between the clean and noisy classes, for the reason
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5 records for its own two-state
model: a start that puts the two classes close together empties one of them and reports convergence.

*Open, and it is the one that could change this: whether the profile curve over `ε` has one hump*
([`parameter_prepass.md`](parameter_prepass.md) §9.3). A scan does not care and a climb does. §10
carries it.

### 3.4 Contamination — a per-sample number this path fits, and only this path can

**`α` is the fraction of a sample's reads that come from another individual** — a second plant in the
tube, a neighbouring library on the same run. It is not a rate per base: the contaminant's reads are
ordinary reads of another genome, so wherever the two genomes agree it is invisible, and only where
they differ does it show.

**Three signatures identify it, and they answer three different questions** (owner, 2026-08-12):

- **It appears in *few* reads, which is what tells it from a real heterozygote.** A heterozygote's two
  alleles are balanced, near half the reads each; a contaminant allele sits at about `α` of them, or
  half that where the contaminant is heterozygous.
- **It appears only at loci that are variable in the population.** The contaminant is another
  individual of the same species, so where the population is monomorphic it carries the same allele as
  the sample and nothing shows. Sequencing error has no such preference.
- **The allele it shows is the allele the population carries there**, where an error's wrong base is
  one of three at random.

**All three need what only this route has.** The first needs a locus's reads to stay together rather
than being folded into a genome-wide count — *which is not the same as deciding, locus by locus,
whether a locus is contaminated; nothing here does that, and §3.4.4 is about what the identification
does rest on.* The second and third need to know **which loci the population varies at and which
allele it carries** — and this route fits exactly that, per locus, as the spectrum of §2.1.
So the frequencies are not an input to be supplied from outside: they come out of the same fit, over
the same loci, in the same pass — **though not as the pooled spectrum, which §3.4.2 is about.** The
histogram route cannot do any of it, having kept neither the allele nor the locus
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2), which is why contamination was a
cohort-gather parameter before this route existed
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5, whose criterion this is).

**One number per sample, fitted alongside the rest**, as a fourth block in §3.3's alternation or after
it — the choice is the estimator's and nothing here turns on it. The evidence is the generic kept
loci: for each sample, how much of its low-fraction alternative-read mass falls on loci the parameters fit says
are segregating, carrying the alleles it says segregate there.

**Two exclusions are requirements of this estimator, not tuning — MEASURED 2026-08-13 on real reads**
([`../reports/contamination_floor_and_duplicated_class_2026-08-13.md`](../reports/contamination_floor_and_duplicated_class_2026-08-13.md)).
Without them the 63 tomato accessions came back with a **median accession 6.5% contaminated**, which
is not a property of the archive.

- **A position the fit judges more likely mismapped than not is not a contamination marker.** Two
  stretches of genome the reference holds once, piling reads onto one position, put a few unexpected
  reads into **every** sample — which is the contamination signature exactly. **Two tomato markers in
  five were such positions**, 20,767 of 52,525. Leaving them in, the median accession reads 0.0684.
- **A sample's reads are scored against every depth its stored code could stand for**, not against the
  middle of the range. The count of disagreeing reads is exact while the depth is a five-bit code
  standing for a *range* above nine reads (§2.2 of the records document), so a heterozygote's read
  share lands away from a half for a reason that is not the sample. Skipping this, a drawn panel with
  **nothing in it** reads 0.025 at ten reads a position and 0.0013 at three.

With both, tomato's median goes to **0.0000** and its worst accession to **0.0090**, while a drawn
panel holding one genuinely 3%-contaminated sample still returns 0.0102 for it against 0.0003 for the
worst clean one.

**And the per-position probability of being mismapped becomes an output of this route, not an
internal.** It is what having every sample at one position buys and nothing else in step 4 can
produce, and it now has two consumers: this estimator, and whatever calls variants afterwards. Four
bytes a position.

**Correcting the sample's own ancestry coordinates for the contamination was measured and is NOT
adopted.** `verifyBamID2` maximises over the fraction and the sample's coordinates together, on the
reasoning that reads from the contaminant drag those coordinates toward the panel average. Undoing
that drag moves a drawn 3% sample from 0.0115 to 0.0166 against a truth of 0.030 — closer — and moves
the worst clean sample from 0.0008 to 0.0046, so **the separation between them falls from 14-fold to
3.6-fold**; searching each axis freely is worse again. **The attenuation and the floor move together**,
so correcting the coordinates buys accuracy in the value at the cost of the thing a threshold actually
needs.

### 3.4.1 The estimator is a two-genotype mixture, and it is standard

**This is `verifyBamID`'s model** (Jun et al., *AJHG* 2012), and the sequence-only form is the one we
need, since ng has no array genotypes for anyone. At each marker, a read's base comes from the
intended sample's genotype with probability `1 − α` and from the contaminant's with probability `α`;
**both genotypes are unknown and both are summed over**, the sample's against its allele frequency and
the contaminant's against the frequency of the population it was drawn from. That is one more mixture
layer on §3.1's likelihood, over the same records, and it fits alongside as one number per sample.

**Two properties of that form to inherit deliberately.** The sequence-only likelihood is *symmetric*
in `α`, so it cannot tell a 20% contaminated sample from an 80% one and the search is restricted to
`α ≤ ½` — a sample swap is invisible to it by construction. And the paper's own list of what biases
`α` **upward** is reference bias, poorly aligned bases and wrong base qualities; the first is already
an open question here ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §8) and the
third does not apply, since ng fits one error rate and does not read base qualities at all.

### 3.4.2 The contaminant's allele frequencies must be individual-specific, not the pooled spectrum

**This is what `verifyBamID2` changed** (Zhang et al., *Genome Research* 2020), and the failure it
fixes is one this cohort is already known to have. The original assumed the contaminant was drawn from
a population whose allele frequencies you supply, and the estimate collapses when that is wrong:
**using African frequencies on East Asian samples returned 2.9% for a true 10% contamination** — under
any 1–3% flagging threshold, so contaminated samples pass as clean. Their fix is to stop using one
frequency per marker. Each individual gets **its own** frequency at each marker, a linear function of
that individual's coordinates in a principal-component space, and the parameters fit **jointly maximises over
`α`, the intended sample's coordinates and the contaminant's** — four components by default, with an
equal-ancestry and an unequal-ancestry model compared by AIC. It removes about 80% of the bias.

**The tomato panel has exactly the structure that breaks the pooled version.** It is landraces from
several regions, which is why §5 rejects homozygote excess as a measure of autozygosity — the Wahlund
effect.

**Measured, 2026-08-12, and the direction is the opposite of what this section used to claim**
([`../reports/joint_contamination_2026-08-12.md`](../reports/joint_contamination_2026-08-12.md),
`examples/ng_joint_contamination_harness.rs`). An earlier version said a sample from a diverged
subpopulation carries alleles the pooled spectrum calls rare, that rare alleles turning up in a sample
are the contamination signature, and so **structure would be read as contamination**. It is not. Fifty
samples, four subpopulations, three reads a site:

| | true `α` = 0.010 | 0.030 | 0.100 |
|---|---:|---:|---:|
| no structure, pooled frequency | 0.0103 | 0.0325 | 0.1065 |
| `F_st` = 0.10, pooled | 0.0037 | 0.0209 | 0.0829 |
| **`F_st` = 0.20, pooled** | **0.0000** | **0.0050** | 0.0584 |
| `F_st` = 0.20, each sample's own subpopulation frequency | 0.0131 | 0.0347 | 0.1044 |

**Structure does not invent contamination; it hides it.** At `F_st` 0.20 a sample contaminated at 3%
comes back at 0.5% and one at 1% comes back at exactly zero — both under any 1–3% flagging threshold,
so **contaminated samples pass as clean**. That is the same direction `verifyBamID2`'s own number goes
— 2.9% returned for a true 10% — so the claim being corrected contradicted the evidence quoted two
paragraphs above it.

**The decision is unchanged and the reason for it is now the measured one.** Individual-specific
frequencies are necessary because without them contamination becomes invisible on a structured panel,
not because a clean panel would be flagged. **What follows for the pipeline is the opposite of what a
false-positive story implies**: the thing to watch is a panel that flags nothing, and §3.4.5's
post-hoc re-estimate is the check that would catch it.

**One warning the measurement adds, and it is about *how* an individual frequency is obtained.** A
per-subpopulation frequency estimated from that subpopulation's own twelve samples adds about **+0.015
to every sample's `α`** and puts **41 to 47 of 50 clean samples** over a 1% threshold — worse than the
pooled frequency on both counts. `verifyBamID2`'s frequencies are a smooth function of an individual's
principal-component coordinates, fitted across the whole panel, and that is the property to preserve:
**borrow strength across the panel, never partition it**. A correct individual frequency needs no
apology — handed the frequency the genotypes were drawn at, a clean panel returns **0 of 50** above 1%
at every divergence tested — so the whole difficulty is in the estimate.

**Measured 2026-08-13: what fitting them actually buys, and the one number that says whose to
trust** ([`../reports/joint_contamination_2026-08-12.md`](../reports/joint_contamination_2026-08-12.md)
§5a). Fifty samples, four subpopulations at `F_st` 0.20, three reads a site, 80,000 loci, one sample
contaminated at 3%:

| | that sample | worst of the other 49 |
|---|---:|---:|
| pooled frequency | **0.0022** | 0.0000 |
| fitted per sample, four axes, shrunk | **0.0320** | 0.0141 |
| each sample's true subpopulation frequency (unattainable) | 0.0305 | 0.0055 |

**It turns a blind estimate into a detectable one** — 2.3 times its own noise floor, where the pooled
frequency cannot separate the contaminated sample from the clean ones at all. It reaches the
unattainable estimate and pays in floor, and **that floor is not a budget knob**: it falls only from
0.0154 to 0.0141 between 20,000 and 80,000 loci, so it is the cost of fitting the frequencies rather
than of having too few markers.

**Shrink each locus's slopes.** A locus whose slopes are indistinguishable from noise must keep only
its intercept — the pooled frequency — so that modelling structure is never worse than not modelling
it. Unshrunk, the same fit returns 0.0443 for that true 0.030, and on an unbalanced panel it
degenerates to the search boundary.

**And the unbalanced panel is where this needs a guard.** With subpopulations of 40, 5, 3 and 2 and
**nobody contaminated**, the worst spurious `α` runs 0.0136, 0.0078, 0.0133, **0.0311** across the four
groups — and 0.2346 in the group of two without shrinkage. **The failure is the opposite of the
intuitive one**: a small group does not fail to get an axis and fall back to the panel average, it
sits at an axis's *extreme*, where a straight line is most sensitive to it, so its own noisy dosages
bend the line towards themselves and its frequency becomes its own echo. By the mechanism above, a
noisy frequency manufactures contamination.

**Decision: emit each sample's leverage, and refuse an estimate above it.** How much of its own fitted
frequency a sample supplies depends only on the coordinates, so it is **one number per sample for the
whole run, computable before a single locus is fitted**. It tracks the damage exactly — 0.027, 0.307,
0.429, **0.857** across those four groups, against a fair share of `(components + 1) / samples` = 0.100
— so a sample supplying more than about half of its own frequency gets `NotIdentified` rather than a
number. That is the same refusal this route makes everywhere else, it costs nothing, and it turns a
silently wrong estimate into an absent one.

**Decision: derive the individual-specific frequencies from this cohort, not from an external panel.**
There is no tomato HGDP; what there is, is the kept loci in every sample, which is the matrix a
principal-component decomposition wants — and the same matrix
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6 already builds for relatedness. So the
components come free with an object this route holds anyway. *Soft:* how many components a plant panel
needs is not four by inheritance, and it is set by the panel rather than by us.

**Built 2026-08-13** — `src/ng/parameter_estimation/joint/contamination.rs`, reading the records rather
than a harness's own panel, with the shrinkage and the leverage refusal above. **Two things the
implementation settled that this section had left implicit.**

**The two genotypes are drawn against two different frequencies, and it is not a detail.** The
sample's own genotype is drawn at *its* frequency — the line fitted at its own coordinates — because
that is a statement about its ancestry. **The contaminant's is drawn at the frequency of whoever was
sequenced beside it**, which by default is the whole panel (§3.4.3), because a neighbouring library on
a plate is not chosen for ancestry. Scoring both against the sample's own frequency is the obvious
reading of §3.4.1's "two-genotype mixture" and it is wrong in the expensive direction, on forty samples
in four subpopulations at `F_st` 0.20 with one contaminated at 3%:

| the contaminant's genotype drawn at | that sample | worst of the 39 clean | a clean panel's mean |
|---|---:|---:|---:|
| **the panel's frequency** — correct | 0.0166 | 0.0032 | **0.0004** |
| the sample's own frequency | 0.0481 | 0.0195 | 0.0099 |

A contaminant from a different subpopulation carries alleles the *sample's own* frequency calls rare,
and rare alleles turning up is the contamination signature — so the wrong reading manufactures about
1% of contamination in every clean sample, which is the threshold itself.

**It finds the sample and it understates the fraction, and the second half is a bias.** With the
correct prior the contaminated sample comes back at 0.0166 against a truth of 0.030, and the value
**does not move with more positions** — 0.0163 at 60,000 against 0.0166 at 12,000 — while the noise
floor on the clean samples does, from 0.0032 to 0.0004. So the separation a threshold needs improves
with the budget and the magnitude does not. **The cause is the one thing of `verifyBamID2`'s that is
not yet built**: it maximises over `α` *and the intended sample's own coordinates together*, and here
those coordinates are estimated from the sample's own reads, which the contamination has already
pulled towards the panel average — so its fitted frequency sits closer to the contaminant's than it
should and the difference the estimator lives on shrinks. **Until that lands, `α` is to be read as
*this sample stands out from the panel* rather than as a fraction.**

### 3.4.3 The contaminant is a sample sequenced alongside this one, and who that was is the user's to say

**The two allele frequencies in this estimator are different questions with different sources, and an
earlier version of this section treated them as one** (owner, 2026-08-13).

- **The sample's own genotype prior** asks what alleles *this* plant is likely to carry. That is
  ancestry, and §3.4.2's individual-specific frequencies are for it.
- **The contaminant's genotype prior** asks what alleles the *stray reads* are likely to carry — and
  the stray reads did not come from a random member of the species. They came from a second plant in
  the tube or a neighbouring library on the same run. **So the population that matters is the set of
  samples sequenced together, not the biological population**, and that set is a list rather than
  something to be inferred.

**That removes the harder half of the problem.** §3.4.2's warning is about how an individual-specific
frequency is obtained without partitioning the panel; nothing of the kind arises here, because the
grouping is stated rather than estimated. It also makes the answer stronger: a sequencing batch is a
few dozen samples, every one of which this parameters fit already holds genotype posteriors for, so the estimator
can sum over **named candidates** and report which batch-mate the stray reads resemble. A real event
favours one donor consistently across loci; a spurious estimate favours none. **That is a check the
scalar `α` cannot make, and it is the form a laboratory can act on.**

**Decision: the batching is a run input with a default, and is never guessed** (owner, 2026-08-13).

- **By default every read group in the run is one batch.** That is the honest statement of what a
  pipeline knows when nobody has told it otherwise, and it makes the contaminant's prior the whole
  cohort's allele distribution.
- **The CLI accepts groups of read groups that were sequenced together**, and the parameters fit uses each
  sample's own group as the contaminant's population. **Read groups rather than samples**, because one
  sample's libraries may have run on different flowcells and the read group is the grain the header
  gives.
- **A run that names any group must name every read group.** An unlisted read group is refused rather
  than swept into a default batch: a user who lists three plates and forgets four samples would
  otherwise get a wrong contaminant prior for those four with nothing said.
- **The batching used travels with `α`**, because two runs under different batchings produce different
  numbers and neither is comparable to the other.

**Why it cannot be inferred, and this is measured rather than assumed** (2026-08-12). The grouping is
not in the alignments of either cohort this project uses:

| | what the header says | what the read names say |
|---|---|---|
| tomato archive | `@RG ID SM PL LB` only — no `PU`, and `LB` synthesised per run accession, so every sample is its own library | rewritten by SRA to `SRR7279481.37559618:TTAGGC:37559618` — the barcode survived, the flowcell and lane did not |
| GIAB HG002 | `PU:unknown` | intact: `HISEQ1:23:H9UD5ADXX:2:1210:8315:21713` — instrument, run, **flowcell `H9UD5ADXX`, lane 2** |

**So one cohort has it in the reads and the other has lost it**, and neither has it where the SAM
specification puts it. A pipeline that inferred a batching from what survives — shared barcodes, say —
would be guessing at exactly the point where guessing wrong is silent. *`ReadGroup` carries no
platform unit today ([`../arch/parameter_prepass_joint_fit.md`](../arch/parameter_prepass_joint_fit.md)
§1.6); reading `PU` when a file declares one is worth doing, but as a **default the user can
override**, never as the answer.*

**What this costs the tomato archive is stated rather than hidden.** With every sample in one batch,
the contaminant's prior is the whole cohort's, which is the pooled frequency — and §3.4.2 measures
what that does on a structured panel: at an `F_st` of 0.20 a sample truly contaminated at 3% comes
back at 0.5%. **So contamination on a cohort with no batching information is a weaker number, and it
is emitted as one** rather than as something comparable to a run where the batches are known.

---

### 3.4.4 Depth, marker count, and what that does to the budget

**Correction to an earlier draft of this section, which said `α` would not be identified at three reads
a site.** That was the wrong axis. No single site is classified: the information is pooled across
markers, so what governs is **how many markers segregate**, with depth entering only because a marker
with two or more reads separates a heterozygote from a homozygote better than one with a single read.
`verifyBamID2` runs on **4× genomes**, well below the depth at which any one site is legible.

**How the accuracy moves with the marker count**, as they report it — their error figure for `α` at
1%, 2% and 5% contamination. **The units are theirs and are not reproduced here**: a mean squared
error of 0.69 at a true `α` of 1% would be a root-mean-square error of 83%, which is not a
measurement anyone would publish, so the column is a relative or rescaled quantity whose definition
we have not confirmed. **Only the ratios below are used**, and a ratio is safe under any monotone
rescaling of a squared error.

| markers | 1% | 2% | 5% |
|---:|---:|---:|---:|
| 1,000 | 0.69 | 0.25 | 0.11 |
| 10,000 | 0.11 | 0.04 | 0.01 |
| 100,000 | 0.02 | 0.01 | 0.01 |

**The census budget is sized to yield about ten thousand segregating sites**
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1), the middle row. Going
to the top row needs about **twenty million** kept positions at a segregating rate near 1 in 200 bp —
ten times the target and ten times the memory — and buys a factor of about **2.3** in the error of `α`
(the square root of 0.11 against 0.02; the ratio is safe to read even where the units are not).

**Decision: do not size the budget for contamination.** These are **priors** for a caller, not a
laboratory measurement of a mixture, and a couple of percent of error in a prior changes nothing a
caller does — the same standard [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md)
§4.1 sets for every parameter here: precision keeps improving as loci are added and *usefulness*
plateaus, so the question is never "how precise" but "precise enough to feel". **Ten times the memory
for 2.3 times the accuracy of one prior fails that test**, and the middle row is already where their
own recommendation sits.

**Measured on this route's own estimator, 2026-08-12, and it confirms the ten-thousand figure while
showing it is tight rather than comfortable**
([`../reports/joint_contamination_2026-08-12.md`](../reports/joint_contamination_2026-08-12.md)).
Fifty samples, three reads a site, each sample scored against its own subpopulation's frequency:

| segregating markers | a clean panel's **worst** fitted `α` | a sample truly at 3% |
|---:|---:|---:|
| 3,434 | **1.85%** | 2.90% |
| 13,748 | **0.86%** | 3.00% |
| 54,735 | 0.53% | 3.20% |
| 218,943 | 0.21% | 3.08% |

**A contaminated sample's estimate is right from 3,400 markers up.** What more markers buy is the
**noise floor on the clean samples** — and that floor, not the estimate, is what a flagging threshold
has to clear. Ten thousand segregating markers puts it at about **1%**, which is the threshold itself.

**So the number to emit beside `α` is the panel's own floor.** Reporting the distribution of fitted
`α` across the panel costs nothing — every sample is fitted anyway — and it turns a constant threshold
into a comparison a reader can make: a sample is contaminated when it stands out from that
distribution, not when it crosses 1%. That is cheaper than the five-times budget increase 55,000
markers would need, and it is the same reasoning §3.2 uses when it requires the two evidence counts to
travel beside the rates.

**SUPERSEDED 2026-08-13 — the floor those figures describe was two defects, and both are fixed**
([`../reports/contamination_floor_and_duplicated_class_2026-08-13.md`](../reports/contamination_floor_and_duplicated_class_2026-08-13.md)).
The table above measures the worst spurious estimate in a clean panel *before* mismapped positions were
excluded from the markers and before a sample's reads were scored across the range its stored depth
code stands for (§3.4). With both, **the drawn floor is 0.0000 for the median accession and 0.0014 for
the worst, at 8,879 markers** — where this table's 13,700-marker row reported 0.86%.

**So the marker budget is not set by the noise floor, and the recommendation to raise it is withdrawn.**
An earlier revision of this section argued for a census of about twenty-seven million positions, on the
grounds that the floor fell from 1% to 0.21% between 10,000 and 219,000 markers and that a file on disk
made the extra positions cheap. **There is no longer a floor to buy down.** The two-million-position
budget yields about ten thousand usable markers, which is where a contaminated sample's estimate has
been right since the first measurement, and the sizing question returns to
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3 — where three hundred and
twenty thousand positions is enough for everything else.

**What does not change is that the value is attenuated.** A drawn sample truly at 3% returns 0.0102,
and correcting for that attenuation costs more separation than it buys (§3.4). So a sample is still
judged against the panel's own spread of fitted values rather than against a constant, for the reason
§3.4.4 gave before the floor was fixed: it costs nothing, since every sample is fitted anyway.

### 3.4.5 Measure it again from the final calls, and report the two

**The better estimate comes free after calling, and it is not a workaround** (owner, 2026-08-12).
`verifyBamID`'s strongest mode is the one where the intended sample's genotypes are **known** —
originally from a genotyping array — because then only the contaminant's genotype has to be summed
over. **Called genotypes are that input.** So once the caller has run, contamination can be
re-estimated from every called site rather than from the kept subsample, with the sample's own side of
the mixture no longer inferred. It is a different and stronger measurement of the same quantity, and
it costs one pass over output that already exists.

**Report both, and let the user decide what to do about a gap.** Two courses, and they belong to the
user rather than to the pipeline: raise the memory budget for the pre-pass and run again, or keep the
prior-stage calls and take the post-hoc number as the contamination estimate the study reports. **What
the pipeline must not do is loop** — re-running the caller under the corrected `α` is a burn-in, and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §1 forbids a second traversal of the
reads as the one constraint not to be given up quietly. Reporting a disagreement respects that;
iterating on it does not.

**And it is a check on the whole route, not only on `α`.** The pre-pass estimate is fitted from a
subsample under a model; the post-hoc one uses every site and a genotype call. **A large gap says the
pre-pass model is wrong somewhere**, and contamination is the only parameter here that can be
re-measured this cheaply after the fact — which makes it worth emitting even on a run where nobody
suspects a contaminated library.

---

## 4. The STR path

**Here the route is HipSTR's model, and the obstacle that kept it out of the per-sample design is
gone.** [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3 sets out the two models side by
side and its last row is the whole story: a per-locus stutter model **needs several samples at one
locus**, and a per-sample genome walk never has them. That document is explicit that the per-stratum choice
"fitted the shape of the pass" rather than being the more accurate model.

**And there is a second thing the per-sample route cannot do that this one can.**
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 records that a per-read tally is exactly
unbiased **when the allele spectrum is handed to it** rather than fitted from the same tally — and
that there is nowhere to get one during the genome walk, so the cheaper object was deferred rather than
rejected. **The kept STR loci are that spectrum**, fitted per locus. This route is where that deferral
comes due.

**Decision: fit the four slippage numbers per (read group × stratum), as today, and weight the
genotype by the locus's own length frequencies — which are a latent drawn from a fitted per-stratum
prior, never a parameter of the locus.** The stratification stays — slippage depends on repeat count
more than on anything else, which is also why the loci are chosen per stratum
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3) — and what changes is the
weighting inside the sum, **in exactly the shape §2.1 uses on the generic path**.

**This path runs after the generic one and depends on it in one direction only — stated here because
it is what lets a run hold one half at a time.** What it takes from the generic path is each sample's
homozygote excess, which weights the genotype drawn from a locus's length frequencies (§4.1). What it
gives back is nothing: contamination is fitted on the generic loci (§4.3), and so are the noise
classes, the frequency density and the homozygote excess itself (§3.3). **And within this path a
stratum is fitted on its own** — slippage is per (read group × stratum), the concentration is per
stratum, and the length spectrum is that stratum's own, so nothing but sums of counts crosses between
them. Both facts are load-bearing for memory rather than only tidy:
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.2 shapes the file around
them, and §11 question 10 is where the cohort-scale consequence sits.

### 4.1 The locus's length frequencies are summed over, not fitted — and this is one design with §2.1, not two

**HipSTR fits each locus's allele frequencies inside its loop**
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3), and taking that literally would
reintroduce the failure §2.1 exists to avoid: each locus brings its own new parameter, so the bias
does not shrink as loci accumulate. **The STR path has already measured what that costs on its own
data.** With the allele spectrum handed to the parameters fit, a per-read tally is exactly unbiased; with the
spectrum fitted from the same tally, the slippage level moves **333-fold depending only on where the
search starts** ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1) — the signature of a
quantity that is not identified rather than badly estimated.

**So the two paths carry one design.** On the generic path a locus's allele frequency is a latent
drawn from a fitted density and integrated away (§2.1.2). Here a locus's **length frequencies** are a
latent vector drawn from a fitted **Dirichlet** whose mean is the stratum's own length spectrum and
whose concentration is one further fitted number per stratum:

```text
locus's length frequencies   ~   Dirichlet( κ_stratum · stratum_spectrum )
each sample's genotype       ~   a draw of P lengths from those frequencies, with F_hom_excess
```

**The concentration is the parameter that says how monomorphic loci are, and it is the whole point.**
A small `κ` makes each locus nearly fixed for one length while the *stratum* still spans many — which
is what a repeat tract actually looks like across a cohort, and what no per-stratum marginal can
express, because drawing every chromosome independently from the stratum's spectrum would make one
locus's samples carry different lengths at random. A large `κ` returns the per-stratum model exactly,
so **the per-stratum route is a special case of this one** and the comparison between them is a
comparison of one fitted number rather than of two designs. In the biallelic case a Dirichlet is a
Beta, so this is the same object §2.1.2 fits on the generic path with the mean and the concentration
named separately.

**The concentration has a second consequence a run has to act on: it is what sets how many tracts a
stratum needs — MEASURED 2026-08-13**
([`../reports/str_stratum_size_sweep_2026-08-13.md`](../reports/str_stratum_size_sweep_2026-08-13.md)).
The four slippage numbers are measured from *reads*, so a deeper cohort gets them from fewer tracts.
**The concentration is measured from tracts** — it says how far each tract's own lengths depart from
the stratum's, which is a comparison between tracts and cannot be sharpened by reading any one of them
harder. Doubling the depth from three reads a site to six halves the scatter of the read-driven
numbers and moves the concentration's not at all: 14.3% against 14.2% at 100 tracts, 7.9% against 8.4%
at 250. **So the per-stratum cap cannot be relaxed because a cohort was sequenced deeply**, and 5,000
tracts is the floor at three reads a site
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6). The number that breaks
*first* is not this one — it is how fast two-repeat slips fall off against one-repeat slips, which
rests on the fifth of the slipped reads that slipped by two — but that one depth does buy back, and
this one it does not.

**Free parameters per stratum: the four slippage numbers, the length spectrum's shape, and `κ`.**
None of them is per locus.

**Measured, 2026-08-12, and it is the strongest result on this path**
([`../reports/joint_str_estimator_2026-08-12.md`](../reports/joint_str_estimator_2026-08-12.md),
`examples/ng_joint_str_harness.rs`). Twenty samples, drawn truths.

- **`κ` is identified**: fitted 0.487 against a truth of 0.500 at three length classes.
- **Per-stratum really is the large-`κ` limit**, so the choice between the two is one fitted number
  as this section claims. The per-stratum model's error on the slippage level runs from **−37.1%**
  where 87% of loci carry one length to **−0.9%** where none does, and one Dirichlet fit is within a
  percent of the truth at every point of that thousand-fold range.
- **At tomato's three reads a site the per-stratum model does not lose accuracy, it loses the
  parameters**: the slippage level comes back **70.9% low** (0.0233 against 0.0800), the direction
  split pins at 1.000 and the fall-off collapses to zero, where the per-locus fit returns +0.3%,
  −0.2% and +1.1%. **That is this route's case on this path, on the cohort the caller is aimed at.**

**Which concentration a real stratum carries — MEASURED 2026-08-13, and this closes the question that
stood here** ([`../reports/str_fit_on_real_records_2026-08-13.md`](../reports/str_fit_on_real_records_2026-08-13.md)
§6.2). The records exist and a genome walk filled them: on 63 tomato accessions the **homopolymer
strata carry 0.52 to 1.56, and the dinucleotides at six repeats carry 5.25**.

*Read that against where the drawn work sits.* Every drawn panel behind §4.1 and behind the
per-stratum cap was drawn at a concentration of **0.5** — which covers the homopolymers and does not
reach the dinucleotides. So the cap of 5,000 tracts
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6) is measured on a stratum
shape tomato's homopolymers have and its dinucleotides do not, and the size sweep's own caveat — that
a differently-shaped stratum could have a different floor — now has a real number attached to it
rather than a hypothetical.

### 4.2 Which lengths a locus may carry, and what bounds the sum

A Dirichlet over thirteen length classes cannot be integrated by a **grid** the way a Beta over one
frequency can: nested quantile quadrature at 24 nodes a dimension needs `24^(classes − 1)` points,
which is 576 at three length classes, 331,776 at five and 1.1 × 10¹¹ at the ±4 the record stores.

**An earlier version of this section concluded from that arithmetic that the sum had to be over a
bounded set of configurations — a locus fixed for one length, or segregating two — and that decision
is withdrawn. It was measured on 2026-08-12 and it costs the slippage level up to eight times over**
([`../reports/joint_str_estimator_2026-08-12.md`](../reports/joint_str_estimator_2026-08-12.md)). Two
things were wrong with it.

- **The support cannot hold a locus carrying three or more lengths**, so the parameters fit's only way to
  explain three lengths among a locus's reads is to say the reads slipped. The fitted slippage level
  tracks that population and nothing else: **+0.9% where no locus carries three lengths, +23.7% where
  18% do, +722% where 99.9% do.** The bias is six to twelve times the spread between draws, so it is
  what the description converges to rather than the luck of one data set. Narrowing the candidate
  lengths to those the reads reported makes the support smaller, not wider, so it does not help.
- **It destroys the property §4.1 was adopted for.** A large concentration *is* the per-stratum
  model, and it is also the regime in which every locus carries every length — which is where this
  support is worst. Under it, "per locus or per stratum is one fitted number" cannot be tested at
  all.

**Decision: integrate the Dirichlet over a fixed low-discrepancy point set** — a Halton sequence in
`classes − 1` dimensions pushed through the stick-breaking Beta quantiles. **256 points, whatever the
class count**, which at the record's nine classes is smaller than the withdrawn support's own 441.
Measured against the grid at three classes it returns the same answers to within 0.3 percentage
points on the concentration and 0.2 on every slippage number, in half the time; at five classes the
grid cannot be run and the point set fits in 96 seconds.

**The points are fixed and the quantile map is continuous in the concentration**, so the objective is
a smooth function of it rather than a jittery one — which is what makes this quadrature rather than
Monte Carlo, and what stops the search chasing sampling noise.

**And it is robust to the model being wrong.** With the truth drawn from the withdrawn support's own
family instead — 71.8% of loci fixed for one length, none carrying three — the point set still
returns the slippage level to **−1.0%**, against that family's own −0.7%. The reverse does not hold.

*Withdrawn with it: the claim that a locus needing more than two lengths is one the guard bucket
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §3.3) is already
watching. The guard catches reads differing by a non-whole number of motif copies, which is a
different thing from a locus segregating three lengths; nothing was watching the second.*

**What survives from the withdrawn decision is the candidate set's reasoning, and it is worth keeping
even though the support no longer needs it** — a
length no read reported in the cohort's 150 reads at that locus is genuinely absent, where a base no
read showed is merely a base that no error happened to produce. *Open, and cheap to settle:* what the
candidate set's size distribution actually is on tomato, which prices the sum. **Settled by:** one
pass over the STR records once they exist.

**Two of the per-stratum route's mechanisms are still needed and are unchanged**: a thin stratum
borrows from its neighbours, and the fitted level is held monotonic along the repeat-count axis
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3). Nothing about per-locus weighting makes
a stratum with eleven loci in it fittable.

**How a thin stratum borrows — MEASURED 2026-08-13**
([`../reports/str_fit_on_real_records_2026-08-13.md`](../reports/str_fit_on_real_records_2026-08-13.md)
§4). **Borrow below 1,000 tracts. Take both neighbouring repeat counts together, never one and then a
test. Refuse below 50 tracts even after borrowing.**

*The both-sides rule is the part that is not a detail.* Slippage rises along the repeat count, so each
neighbour's level is displaced from the thin stratum's — by about 30% a count on the drawn panel.
Taking one neighbour and then testing whether the floor is cleared keeps whichever side was reached
first and carries that displacement whole: the borrowed slippage level comes back **23%** from the
truth. Taking both sides of a distance together leaves the two displacements pointing opposite ways
and it comes back **1.2%** from the truth. The obvious implementation is the wrong one.

*And borrowing never lost in the measured range*, so the floor is not where it stops paying: at 50
tracts, fitting alone against borrowing, the slippage level is 9.6% against 1.8%, the fall-off 35.2%
against 4.2%, the concentration 16.7% against 3.4%; borrowing also won at 250 and at 1,000. **The
floor is set at 1,000 because that is where a stratum's own answer costs about a percentage point**,
and keeping its own answer is what preserves a repeat-count axis that varies.

**And borrowing can erase the axis it is defending, which is why the emitted parameters must say what
stood behind them.** On the human benchmark trio — 216 tracts in 32 strata, from a 452 kb region set —
every one of the fifteen homopolymer strata fell below the floor, so each reached for its neighbours,
each ended up pooling all fifteen, and all fifteen were handed the identical four numbers for repeat
counts 8 through 23. **The stratification was flattened to nothing and no field in the output said
so.** That is the rule working as specified, not a defect in it — but it means *the answer a run gets
depends on how much genome it walked*, and a reader of the emitted table cannot see it.

**So each emitted per-stratum parameter carries two more values: how many tracts stood behind its own
answer, and which strata were pooled to produce it.** They cost nothing — the estimator holds both
already — and without them a flattened axis is indistinguishable from a flat one.

**How much of the genome's STR loci this route holds decided its standing, and the answer is: all of
them.** Measured on tomato SL4.00 at the calling floors `[8, 6, 6, 6, 5, 4]`
(`examples/ng_joint_loci_probe.rs`, 2026-08-12): **141 strata holding 462,701 loci**, so a cap set
above the largest stratum keeps every one of them
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.5). **This route therefore
holds the same STR loci as the per-stratum histogram *and* remembers which was which** — a far
stronger position than the generic path's one position in three hundred and ninety-one, and it makes
§8's second measurement a like-for-like comparison on this path where it is not one on the other.

### 4.3 Contamination is not fitted from the STR loci, and the reason is sharper than "noisier"

**Every signature §3.4 uses has an STR analogue on paper**: the alleles are lengths, this route fits a
per-locus length spectrum, and a contaminant would show a length the sample does not carry, in few
reads, at a locus the population varies at. **Two things break it, and only the first is about noise
being large.**

- **At long tracts the noise is bigger than the signal.** At six repeats and above, **2 reads in 100**
  sit at a length the sample does not carry ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
  §5), so a 1% contamination arrives buried under twice as much of the noise model's own output. *This
  half of the objection is about size and it does not hold across the whole range — see below.*
- **The second and third signatures lose their power, and this is the decisive one.** They work on the
  generic path because sequencing error is *indifferent* to the population. An error falls everywhere,
  including at the 199 positions in 200 the population is monomorphic at, where contamination cannot
  show at all because the contaminating individual carries the same allele; so stray reads
  concentrated on the variable positions mean contamination and not error. And when an error does
  land on a variable position, its wrong base matches the population's alternative allele about one
  time in three, where a contaminant read matches it always. **At a repeat tract the noise imitates
  the variation instead of being indifferent to it.** A genuine STR allele differs from its neighbour
  by one repeat unit, and one repeat unit is exactly the step stutter takes, so the wrong length
  stutter produces *is* a length the population segregates at rather than one of many wrong lengths.
  **Even if stutter were a hundred times rarer than contamination those two tests would return the
  same answer for both**, which is why this objection is not about size.

**But both objections weaken at short tracts, and they weaken together** (owner, 2026-08-12). The
"everything segregates" half is a property of **long** tracts: slippage varies twenty-two-fold across
repeat counts, and so does mutation, so below four repeats a great many loci are monomorphic across a
cohort. That restores the contrast the second signature needs — a population of quiet loci where
stutter shows and contamination cannot. And it is the same end of the range where the noise is small:
**9 reads in 10,000 below four repeats against 2 in 100 at six and above**
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5), so a 1% contamination is ten times the
stutter background there and half of it at the other end.

**Decision: fit `α` on the generic loci, and settle the STR question per stratum rather than in
general.** The generic set supplies `α` either way, so nothing waits on this. **Two numbers decide
it, both per stratum and both cheap**: what fraction of the stratum's loci segregate across the
cohort, and what fraction of stutter products land on a length the cohort segregates at. *Leaning:*
the long-tract strata are useless for this and the short-tract ones may be usable — the opposite of
their standing everywhere else in the STR path, where the long tracts carry the information.

---

## 5. Inbreeding — two quantities, and from here they have two names

**This route produces an inbreeding coefficient, and it is not the same quantity the per-sample route
produces.** That has to be stated before anything else, because both are one number per sample between
0 and 1, both are called *the inbreeding coefficient* in the literature, and a consumer handed the
wrong one gets a plausible answer.

**So neither is written `F` from here on.** A bare `F` in any step-4 document means whichever the
author had in mind, which is the failure itself; the two names below are used in prose, in the emitted
parameters and in code. *(The per-sample route's own spec still says plainly `F` throughout —
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6 — because that path is built and is
not being reopened for a rename. It means `F_autozygosity` everywhere.)*

| | `F_autozygosity` | `F_hom_excess` |
|---|---|---|
| what it measures | the fraction of **this genome** where both copies descend from one ancestral copy | how much less heterozygous this individual is than random mating in the panel would predict |
| how | walk 100 kb windows, ask each whether it sits inside a long stretch nearly free of heterozygotes | `1 − Hobs/Hexp`, with `Hexp` from the fitted spectrum |
| needs | one sample | the whole panel |
| produced by | the per-sample genome walk ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6) | this route, §5.1 |
| what the caller's prior asks for | **this one** | — |
| the literature's name for it | `F_ROH`, the inbreeding coefficient from runs of homozygosity | `F_IS`, the within-population fixation index |

**The caller's genotype prior mixes `F_autozygosity·π_i + (1−F_autozygosity)·π_i^ploidy`**
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6) — that is literally the question
*did these two copies come from the same ancestral copy?*, so realized autozygosity is what it wants.
Homozygote excess coincides with it only when nothing else suppresses heterozygosity.

**Three things make them come apart, and all three are live on the tomato cohort.**

- **Population structure.** A panel that is really landraces from several regions shows a homozygote
  excess in every individual with no individual's parents being related. Homozygote excess counts it as
  inbreeding; autozygosity does not. This is the Wahlund effect, and it is what
  [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.2 rejected the ratio for.
- **False heterozygotes.** Collapsed paralogs and mismapping add heterozygous sites roughly uniformly.
  Adding them at five times tomato's real rate of one per kilobase **moves the runs estimate not at
  all** — both states' heterozygote rates lift together and `F_autozygosity` reads only the gap — while the
  whole-genome heterozygosity the ratio reads **inflates eight-fold**
  ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
  §3.3).
- **Old inbreeding.** 100 kb windows resolve runs of about 300 kb and longer. Autozygosity from far
  enough back to have been broken into shorter tracts is invisible to runs and still suppresses
  heterozygosity, so the ratio sees it and the runs do not.

### 5.1 Decision: fit it, emit it under its own name, and never as a substitute

Fit one `F_hom_excess` per sample here, in step 3 of §3.3's alternation, as the departure from
Hardy-Weinberg proportions the fitted density predicts. **Emit it as a distinct, differently-named
parameter from the autozygosity `F_autozygosity`** — never as an alternative value for the same
field. A caller must not be able to receive one where it expected the other.

**It is constrained to `[0, 1]`, and the constraint is load-bearing rather than cosmetic.** `F_IS` is
negative wherever an individual is *more* heterozygous than random mating predicts, and an
unconstrained fit will go there — which is precisely how the duplicated loci of §2.2 would escape.
A duplication the reference does not carry shows about half its reads disagreeing at every position
where the copies differ, in every sample; §2.2's whole argument that they cannot be absorbed
elsewhere rests on the homozygote excess being unable to move in that direction. A negative
`F_hom_excess` would let one sample's mismapping be booked as biology, silently and with a plausible
number. **The constructor refuses a value outside `[0, 1]`** — a fit that wants to leave the interval
is reporting a modelling failure, and the right response is the diagnostic §5.1 lists rather than the
number.

**What that constraint catches, and what it does not — MEASURED 2026-08-13**
([`../reports/duplicated_class_identification_2026-08-13.md`](../reports/duplicated_class_identification_2026-08-13.md)).
On an **outbred** panel a missing third class drives the coefficient to **−0.09** and the constraint
refuses it, so the failure is loud and the run stops. On a **selfing** panel the same missing class
lands it at **0.4471 against a true 0.5942** — a perfectly legal value, arrived at silently, with the
two error rates and the noisy share all within a percentage point of the truth. **So the constraint is
not a safeguard against this, and nothing else in the fit is either**: on the panels this caller is
aimed at, a quarter of the inbreeding coefficient can go missing with every other number looking
right.

**Why fit it at all, given the caller wants the other one:**

- **It is the only external check the autozygosity estimate has.** Nothing else in step 4 measures
  inbreeding twice, and the runs estimator's known failure mode is a confident zero
  ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5: a start that guesses the
  inside state's heterozygote rate far below the truth returns `F_autozygosity` = 0.0000 on a genome with 29% of
  its length in runs, and reports convergence).
- **A disagreement between the two is diagnosable rather than confusing**, because the three causes
  above push in known directions: structure raises the excess above the autozygosity, artifacts push it
  below, and old inbreeding raises it.
- **It costs one number per sample inside a fit that is running anyway.**

### 5.2 Half the circularity opens here, and it is worth being exact about which half

[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.3 and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3 both carry the same trap: the panel's
expected heterozygosity is computed as `Hobs/(1 − F_hom_excess)`, so taking `F_hom_excess` from the ratio and then computing
the diversity from it returns whatever was assumed. **That loop is what makes the ratio a diagnostic
rather than an estimate in the per-sample design.**

**Under this route half of that loop opens, and it is worth being exact about which half.** Expected
heterozygosity comes from the fitted frequency density — measured across the panel — and never from
`Hobs`. So the *diversity* is no longer whatever inbreeding was assumed, which is what §5.3 is about
and it is a real gain.

**What does not open: `Hobs` is not an independent measurement of the same thing.** It is derived from
genotype posteriors (§3.2) whose prior is the fitted density *and this sample's own
`F_hom_excess`* — so at a locus with no reads the posterior is the prior, whose heterozygosity is
exactly `Hexp · (1 − F_hom_excess)`, and the ratio returns the fitted value by construction. At
tomato's three reads a site that is about a fifth of the kept loci (§3.2).

**So `1 − Hobs/Hexp` is a restatement of the fitted parameter, not a check on it.** It is emitted
because a consumer asks for the ratio in those terms, and it is emitted with the two evidence counts
§3.2 requires. **The check on `F_hom_excess` is the other estimator**, `F_autozygosity` from the
per-sample genome walk, which reads a different feature of the data entirely — the ratio and the runs
disagree in known directions (§5, three causes) and that is what makes the disagreement diagnosable.
Nothing inside this route checks it against itself.

### 5.3 What this does to the cohort's diversity

**The diversity stops needing an inbreeding coefficient at all.** [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)
§3 computes `Hexp = mean over samples of Hobs/(1 − F)`, which inherits every sample's `F_autozygosity` and its
uncertainty, and breaks entirely if the runs estimator returned a confident zero. Here `Hexp` is read
off the fitted frequency density directly:

```text
Hexp  =  ∫ π(f) · 2 f (1 − f) df
```

— the chance that two copies drawn at random **from the population** differ, which is what expected
heterozygosity means. *It needs no finite-sample correction, and that is a property of having fitted
a population frequency rather than a panel count: the `2N/(2N − 1)` factor a sample allele-count
spectrum would need (§2.1.1) is the correction from panel to population, and here there is nothing to
correct.* **That is a genuine improvement independent of everything else in
this document**, and it is worth flagging to whoever finishes the comparison: even a route that lost
on every other axis would still be the better source of this one number.

---

## 6. What this route cannot produce

**Anything local.** A 100 kb window holds about a hundred thousand sites but only a couple of hundred
kept ones ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §2), so at one
heterozygote per kilobase it carries a **fraction of one** expected heterozygote instead of a hundred
— far too thin to tell a run of homozygosity from ordinary sequence. *The reason is thinness, not
independence: the kept positions are a uniform thinning, which leaves most of them with a neighbour
within a kilobase and does not make them independent of one another
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1).*

So **runs of homozygosity, and the `F_autozygosity` the caller's prior actually reads, are not
available by this route and no budget within reach changes that.** The windowed histogram and its
estimator stay in the genome walk whatever the comparison decides — which means this route can never make the
genome walk a pure accumulation pass, and the genome walk's largest object (37 MB per sample on tomato, 145 MB on
human) survives regardless. **That is the honest ceiling on what deleting the per-sample fits could
ever save.**

### 6.1 What it does at one sample, and why saying so matters

**This caller must also run on one sample**, and the route has a clean answer that no section stated
before: at `N = 1` the frequency density is still fitted, the integral over `f` is still taken, and
the likelihood becomes [`parameter_prepass.md`](parameter_prepass.md) §3's per-sample estimator with
the pooled genotype frequencies replaced by what `π` implies for one individual. **Nothing breaks and
nothing is a special case in the code.** So this route is a *generalisation* of the per-sample fit
rather than a cohort-only tool.

**What is lost at one sample is exactly the two things a second sample buys**, and both should be
emitted as absent rather than as numbers:

- **The per-locus site class.** §2.2's whole mechanism is several samples at one locus disagreeing;
  with one, the class posterior is the mixture weight and the parameters fit *is* the blind two-class model
  [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1 already ships. §8's third
  measurement is the curve that says how many samples it takes to matter.

  **The duplicated class is the exception, and what carries it at one sample is depth.** Its cohort
  discriminator — no sample homozygous for the non-reference allele — has no power at all here. Its
  other one does: the position's own depth against the sample's median is a comparison inside the one
  sample's own genome, and it discriminates from about twenty-five reads a position (§2.2). **So a
  single sample at 30× or better keeps the class, and a single shallow sample emits it as absent
  rather than as zero.**
- **`F_hom_excess` separately from the density.** One individual's heterozygote deficit and a density
  concentrated near zero are the same observation, so the pair is not identified. Emit
  `F_hom_excess` as *not identified* below the sample floor §12.5 measures, never as a fitted zero.

Contamination survives at one sample in principle — its evidence is *within-sample*, the low-fraction
alternative reads concentrated on segregating positions — but the positions it needs to be told are
segregating come from the panel, so at `N = 1` there is no panel and it is `NotIdentified` too.

**Linkage, haplotypes, and anything else reading a stretch of genome are out of reach for a different
reason, and it is a matter of shape rather than of budget.** A uniform thinning throws away exactly
the close pairs where linkage information lives and keeps unlimited distant ones, which carry none;
the instrument for those is a small clustered budget
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §8), which nothing in step 4
asks for today.

---

## 7. Cross-cutting concerns

**The scheduling cost is the one a manager needs, and it is not memory.** This route puts a **barrier**
in the pipeline: no locus can be called until every sample has been walked, because the parameters fit needs every
sample at once. The per-sample route has no such barrier — a sample's parameters are ready when that
sample's genome walk ends. **Nothing about the kept loci makes this avoidable**; it is what "fitting against
the locus's frequency in the cohort" means. A run that adds one sample later must refit, and the
parameters every earlier call was made under have changed.

**Memory** is [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6, and
**this route adds no per-sample object beyond the records** — a revision until 2026-08-14 added a
coverage-by-window summary at 1.6 to 6.2 MB a sample, and removing it is most of what that change
bought. What the *parameters fit* adds is working memory only: one posterior over the frequency
quadrature per locus, one number per node, held for one locus at a time.

**The resident bill scales with the cohort, and a cohort of thousands breaks it — OPEN, 2026-08-13.**
At two million positions a sample's records are about 6 MB, so a thousand samples is 6 GB and five
thousand is 30 GB. **No budget knob fixes that**, because the cost is per sample rather than per locus.
The way out is that **nothing here needs every sample**. The population-level parameters — the
per-locus allele frequencies, the diversity, the inbreeding coefficient, the error rates, the share of
noisy loci — are properties of a population measured from a panel, and a panel of a few hundred is
enough: the same sweep run at fifty, two hundred and a thousand samples shows the read-driven
parameters improving between fifty and two hundred, and the count of usable markers rising only about
14% between fifty and a thousand
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3.1, §4.3.2). Contamination is
the opposite kind of number — per sample, one at a time, given the frequencies — so it streams (§3.4.4).

**Proposed, and not measured: bound the parameters fit's resident set by a subsample of *samples*, chosen the same
data-blind way the loci are, and compute every remaining sample's per-sample parameters in a pass
afterwards.** Memory then stops depending on cohort size along either axis — loci bounded by the
budget, samples bounded by the subsample. **What nobody has measured is how large the subsample has to
be**, and it is the same shape of experiment as §4.3's site-budget sweep in the loci document, run
along the sample axis instead. §11 question 8 carries it.

**And a subsample is a statistical answer, not a budget.** A run handed five thousand samples and a
memory ceiling has to stay under it whatever the sweep says the parameters fit would like, which is a second
question with its own levers — reading the per-sample files locus-major so peak memory is set by how
many loci are held rather than how many exist, and keeping the iterative part on a small resident set
of loci while the full census is swept once. **§11 question 10 carries that one**, and it is the one a
user of this caller feels.

**Compute, honestly.** One evaluation of §3.1's likelihood over the whole generic set is on the order
of `loci × samples × genotypes × (1 + 3 · quadrature nodes)` operations — two million by fifty by
three by about fifty, roughly 10¹⁰ — and the parameters fit needs tens of iterations of that. **Nearly all of it
is spent on positions where nothing happened**, and that is the lever: at a locus where no sample
showed a non-reference read, the three-allele sum collapses to one term and every sample's inner sum
depends only on its depth bin, so the per-sample factors come out of a `(depth bin × genotype)` table
of 21 × 3 entries built once per candidate. **It is embarrassingly parallel over loci** and needs no
communication between them within an iteration. **This is arithmetic, not a measurement**, and §8's
last item replaces it. If it dominates, the budget is the knob and
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4 is how to set it.

**Determinism.** The parameters fit sums over loci in a fixed order and over samples in a fixed order, so no
parameter varies with thread count. Multiple starting points are enumerated, not sampled.

**Errors.** The twelve recording terms of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §5 must all match across
samples, and **a mismatch must refuse, not average**. Seven say the samples were asked for the same
loci; one is a digest of the loci each one actually kept, and it is the only one that would catch a
hash function or a threshold's arithmetic changing underneath two otherwise identical runs; the last
five say in what **units** the evidence was recorded — the depth ladder, the two caps, the stratum
counts and the window size — and a mismatch in any of those leaves every other value agreeing while
two samples' rows mean different things.

---

## 8. The comparison this route exists for

**Both routes run on the same data and their estimates are compared, on three axes: accuracy, memory
and speed.** The reference is the per-sample whole-genome histogram route
([`parameter_prepass_generic.md`](parameter_prepass_generic.md),
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)), and the medium is **synthetic data**, because
only there is the truth known — draw genotypes at known allele frequencies, draw depths, draw reads at
a known error rate, and fill both routes' accumulators from the same draw.

**Two things that framing does not cover, and both are cheap to add rather than reasons to change
it.**

- **Synthetic data answers "does each estimator recover what it was given", not "is the model right
  about real reads".** The noisy-locus class exists because 818 non-variant HG002 loci carried three
  or more alternative reads where one error rate predicts 29 (§2.2) — that is a property of real
  mismapping and error-prone context, and a generator reproduces it only if someone builds it in, in
  which case the parameters fit is being graded against an assumption. **So measurement 3 below is deliberately
  not synthetic**, and it is the one that grades the model rather than the arithmetic.
- **Memory and speed are only comparable at realistic scale.** The numbers that decide this route —
  the records held for every sample at once, and the per-sample windowed histogram the other route
  holds in flight — do not appear on a small fixture. Run the cost measurements at the shape of a real
  cohort: a genome-sized position budget, fifty samples, and the concurrency the genome walk would actually
  use.

Five measurements, in the order they can be made.

1. **Bias, against synthetic truth.** Draw genotypes at known allele frequencies, draw depths, draw
   reads at a known error rate; fill both routes' accumulators from the same draw; fit both. **Repeat
   the simulation** so that bias separates from imprecision — imprecision shrinks with the budget and
   bias does not, and bias is what a subtly broken selection rule produces
   ([`parameter_prepass.md`](parameter_prepass.md) §4.1). Run it at more than one coverage, since this
   route loses relatively more where reads are scarce, and at `P = 2` and `P = 4`.
2. **How well each uses the data it is given.** Report each fit's error against truth **with the count
   of observations behind it**, because the two routes are not handed the same evidence: the generic
   set holds one position in a few hundred, while the STR set may hold every locus (§4). A route that
   is less precise on a four-hundredth of the data is not thereby worse; a route that is *no more*
   precise on all of it is.
3. **The two-class residual, and it needs a cohort — an earlier version of this item asked for it on
   one sample.** §2.2's claim is that **several samples at one locus** separate a mismapped locus from
   a heterozygous one: a collapsed paralog raises the alternative-read fraction in every sample at
   that position, a real heterozygote raises it only in the samples carrying the allele. **At one
   sample there is nothing to disagree with** — the per-locus class posterior is the blind mixture
   weight, and the parameters fit reduces to the model it is supposed to beat (§6.1). So *"refit HG002's
   551,843 confident-region loci and see whether heterozygosity comes in below 1.09×"* cannot move the
   number it is aimed at, and this item is three arms instead:

   | arm | samples | truth | what it can show |
   |---|---:|---|---|
   | **the GIAB trio at 30×** — `benchmarks/giab/per_sample/bam/30x/`, HG002/HG003/HG004, a v4.2.1 benchmark VCF for each | 3 | yes | real mismapping and real error-prone context, graded against truth. Three samples is the mechanism at its weakest useful strength, which is exactly why it is the arm that cannot be argued with |
   | **a drawn cohort with planted noisy loci**, refitted at 2, 3, 5, 10, 25 and 50 samples on the same drawn loci | 2 → 50 | yes | **where the mechanism starts paying**, which is the number the trio alone cannot give and the one a user needs |
   | **the tomato cohort** | 63 | no | that it behaves at the real shape — the depth, the sample count and the mismapping a crop reference has |

   **What is reported is the trend in sample count, not a single ratio.** The claim being tested is
   that the excess falls as samples are added, and a claim of that shape is settled by a curve. The
   drawn arm is what makes the curve affordable and the trio is what keeps it honest, since a
   generator only reproduces the mismapping someone built into it (see the framing above).

   **Two things to hold fixed across the arms, or the curve measures something else.** Depth, because
   the trio is 30× and tomato is three reads and §2.2's excess is itself depth-dependent; and the
   locus set, because HG002's confident regions are not a uniform thinning of the genome. Report the
   drawn arm at both depths.
4. **The least heterozygous samples, separately.** Report every measurement above **split by the
   sample's own heterozygosity**, never pooled over the cohort. The reason is not that the estimate
   gets noisier there — on a proportion's own scale it gets tighter (§3.2) — but that **the failure
   mode changes**: at the cohort median a shortfall shows up as scatter, and at tomato's floor of
   0.149 per kilobase, where at least 97 in 100 positions with an alternative read are artefact, a
   mis-fitted background shows up as a **confident wrong number** that no cohort mean would reveal and
   no extra loci would cure. Include a drawn sample an order of magnitude below that floor, for a
   selfing line. **It is also the axis on which the two routes are least interchangeable**, since the
   per-sample route reads every site and this one reads a subsample of them.
5. **Memory and wall clock, measured rather than computed.** §7's figures are arithmetic. Report the
   records at rest, the parameters fit's working set, and the parameters fit's wall clock at several core counts, on the
   whole tomato cohort — the run that stresses sample count, which is the axis this route's cost
   scales on. **Report it at each budget of
   [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3's downward sweep**, since
   the two questions share one experiment.

**What must be recorded so that nobody re-runs the genome walk to finish this**: both routes' fitted values,
the evidence count behind each, the memory each cost, and the wall clock, written up in
`doc/devel/ng/reports/`.

**Two comparisons that are not like-for-like, and reporting them as if they were is the mistake to
avoid.** Inbreeding is two different quantities (§5) and must be reported as two rows, never as one row with
two values. `Hobs` and `π_hom_alt` are free parameters in one route and derived quantities in the
other (§3.2), so their disagreement has a specific candidate cause and should be read against it.
**Contamination has no counterpart to compare against** — only this route produces it (§3.4) — so what
the report carries is the value **and the depth it was fitted at**, with *not identified at this
depth* surviving as that rather than as a blank.

---

## 9. Reuse over rewrite

| what | existing code | how it is reused |
|---|---|---|
| the sum over genotypes at one sample's locus | `src/ng/parameter_estimation/generic/noise_model.rs` (`NoiseModel`) | used as the innermost term of §3.1, unchanged. It is already the seam both paths implement |
| the climb over mixture weights | `src/ng/parameter_estimation/fitting/mixture_weights.rs` | the three masses of the frequency density (§2.1.2) are mixture weights, so the same climb applies to them. **The Beta's two shape numbers are not weights** and need their own climb |
| the two site classes | `src/ng/parameter_estimation/generic/` (`SiteNoise`) | the parameters are the same three; what changes is that the class becomes a per-locus latent variable (§2.2) |
| alternating between coupled blocks | `src/ng/parameter_estimation/generic/coupled_fit.rs` | the same loop shape and the same termination handling, with three blocks instead of two (§3.3) |
| provenance and evidence count | `src/ng/parameter_estimation/mod.rs` (`Provenance`, `Estimate`) | used as-is; every parameter here carries both, as in the per-sample route |
| the frequency spectrum | [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4 | **not a separate step under this route** — it is §2.1's outer integral, and the panel allele-count spectrum falls out of the same fit (§2.1.2). That section's `realSFS` reference belongs to the count form this route does not use, for the reason §2.1.1 gives |

*The two companion documents carry their own reuse: the depth ladder
([`records`](parameter_prepass_joint_records.md) §2.2) and the even spread across a stratum
([`loci`](parameter_prepass_joint_loci.md) §3.1), both already written.*

**No parity oracle.** Neither route is a port of the other, and agreement is the thing being measured
rather than the thing being asserted (§8).

---

## 10. Deferred, with a recommended home

- **Runs of homozygosity, and `F_autozygosity`.** Not deferred so much as unreachable here (§6).
  **Home:** [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6, where it already lives
  and where it stays.
- **Relatedness and read-group grouping.** Already specified and unchanged by this route, which reads
  the same records they do. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6–§7.
- **Contamination's arithmetic.** §3.4 settles the evidence, the three signatures and the states `α`
  may return; how the mass is turned into a fraction is not written. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5, which defers the same arithmetic
  and already states the criterion.
- **The population genetics above diploidy.** The likelihood here is ploidy-generic as
  [`parameter_prepass.md`](parameter_prepass.md) §3 requires, but "allele count in the panel" and
  "homozygote excess" both need restating above `P = 2`. **Home:** the same spec that owes the
  diploid-only definitions, [`parameter_prepass.md`](parameter_prepass.md) §8.
- **Adding a sample to a finished run.** This route's barrier (§7) makes it a refit rather than an
  increment, and what a pipeline should do about that is the pipeline's decision. **Home:** whichever
  spec settles incremental cohorts; nothing here forecloses it.

---

## 11. Open questions

*The two that belong to the companion documents — how many STR strata there are and how many loci each
holds, and where the generic budget starts to matter — are
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6.*

1. **Does per-locus weighting close the two-class residual?** — OPEN. §2.2 argues it should and nobody
   has measured it. *Leaning:* partly — mismapping is a per-locus property and this route can see it,
   but the 1.09× also contains error-prone sequence contexts that a single locus's fifty samples
   describe no better than a stratum's thousand loci do. **Settled by:** §8's third measurement.
2. **Does the profile curve over the noise rates have one hump?** — OPEN, inherited from
   [`parameter_prepass.md`](parameter_prepass.md) §9.3 and **sharper here**, because §3.3 climbs where
   the per-sample route scans, and a climb can be trapped where a scan cannot. *Leaning:* one hump,
   with the two-class model the main reason for doubt — a mixture's component parameters are where
   multimodality classically appears. **Settled by:** plotting the curve on synthetic data before the
   starting points are fixed.
3. **How many starting points, and spanning what?** — **PARTLY CLOSED 2026-08-12: spanning the
   separation between the clean and the noisy class, and it is worth more than any other choice in
   the parameters fit.** Three starts, ten samples, 40,000 drawn loci, the same data throughout:

   | where the search began | clean error rate | `Hexp` |
   |---|---:|---:|
   | the classes far apart — `ε_clean` 5 × 10⁻⁴, `ε_noisy` 2 × 10⁻¹ | **−1.3%** | −0.8% |
   | a middling start — 6 × 10⁻³ and 2 × 10⁻² | −8.9% | +2.4% |
   | **the classes close together** — 2 × 10⁻² and 6 × 10⁻² | **−45.8%** | +10.6% |

   **A start that puts the two classes near each other collapses them into one and reports
   convergence**, which is the failure
   [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5 records for its own two-state
   model — measured here rather than inherited, at 46% of the clean rate. Take the best-scoring of
   several starts. *Still open:* **how many**. Three is what the measurements above used and it was
   enough for one of them to land well; whether nine is needed, as the per-sample route settled for
   its own model, wants the profile curve of question 2.
4. **How many principal components does a tomato panel need for the individual-specific allele
   frequencies?** — OPEN (§3.4.2). `verifyBamID2` defaults to four on human reference panels, and a
   crop panel's structure is not human structure: a landrace collection may need more, and an
   inbred-line panel behaves differently again. *Leaning:* decide it from the panel by the usual
   scree-plot criterion rather than by inheriting four, and record what was used beside `α`. **A
   number fitted under too few components is biased *downward*** — §3.4.2's measurement: at an `F_st`
   of 0.20 a pooled frequency, which is what zero components gives, returns 0.5% for a true 3%. So too
   few components hides contamination rather than inventing it, and a panel that flags nothing is the
   symptom. **Settled by:** the components come out of the same matrix relatedness uses
   ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6), so this is a plot rather than an
   experiment. **One constraint the measurement adds**: whatever produces the individual frequency
   must borrow strength across the whole panel. A frequency estimated from a subpopulation's own
   twelve samples puts 41 to 47 of 50 *clean* samples over a 1% threshold, which is worse than using
   no structure at all.
5. **Do the short-tract STR strata contribute to contamination after all?** — OPEN (§4.3). *Leaning:*
   possibly, and it does not matter much, since the generic loci supply `α` either way. **Settled by:**
   two per-stratum counts — what fraction of a stratum's loci segregate across the cohort, and what
   fraction of its stutter products land on a length the cohort segregates at.
6. **Is a monomorphic mass, a fixed-alternative mass and a Beta enough shape for a domesticated
   selfing panel?** — **CLOSED 2026-08-12: yes, and a second Beta buys nothing.** Fitted against a
   drawn cohort whose true density has two bumps — the shape two diverged subpopulations leave, and
   one a single Beta cannot make — the four-number shape returns `Hexp` **4.9% high**. Giving it a
   second bump, three more fitted numbers, returns **5.8% high**: no better on the truth the second
   bump exists for, and the likelihood gain is 10⁻⁶ nats a locus, which an extra parameter always
   buys. Against a one-bump truth the two are −0.8% and +0.4%. **The noise rates barely feel the
   misspecification at all** — they read the density as a background, and they move by about two
   percentage points between the two truths.

   *What this does not license* is calling the emitted density a site-frequency spectrum. Five percent
   is small for a prior, and on an autogamous panel `1 − F_hom_excess` **is** `Hobs/Hexp`, so it
   passes into the inbreeding coefficient at full size
   ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2). The number to watch is
   there, not in the caller's genotype prior.
7. **How far do the density's shape and the per-sample inbreeding trade against each other at three
   reads?** — OPEN (§3.3). Both control how often a genotype comes out heterozygous, and they are
   separated only by *where* the variation sits: across loci for one, within a locus for the other.
   *Leaning:* separable at fifty samples and not at five. **Settled by:** §12.9's profile.
8. **How many samples does the parameters fit actually need at once?** — OPEN, raised 2026-08-13 (§7). A thousand
   samples' records are 6 GB and five thousand are 30 GB, so a cohort of that size cannot be held, and
   the proposal is to fit the population-level parameters on a data-blind subsample of samples and
   compute the per-sample ones in a streaming pass afterwards. *Leaning:* a few hundred is enough —
   the measured sweeps show the read-driven parameters settling between fifty and two hundred samples
   and the usable-marker count rising only 14% between fifty and a thousand
   ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3.1, §4.3.2). **Settled
   by:** the site-budget sweep of that document's §4.3, run along the sample axis — refit at 1,000,
   200, 100 and 50 samples drawn from one cohort and report each parameter's error against the drawn
   truth, one row per parameter, beside what the subsample cost at rest.
9. **Does re-fitting `α` alone on a larger census keep the noise floor it promises?** — **CLOSED
    2026-08-13 as moot.** The question existed to price a two-step scheme for buying down a noise floor
    of about 1%. That floor was two defects rather than a sampling limit — mismapped positions used as
    markers, and a depth code read as a point rather than a range — and with both fixed the drawn floor
    is 0.0000 for the median accession at 8,879 markers (§3.4, §3.4.4). **There is nothing left to buy
    down, so there is no large census to re-fit on.**
10. **How do these estimates run inside a memory budget when the cohort is thousands of samples?** —
    OPEN, raised 2026-08-13 (owner). Question 8 asks how few samples the parameters fit *needs*; this asks what it
    *does* when it is handed five thousand and a ceiling it must not cross. The two have different
    answers because one is a statistical choice and the other is an engineering one, and a run must
    survive even where the statistical answer is "use them all". **The size in play:** a sample's
    records are about 6 MB at two million loci, so five thousand samples are 30 GB, and no locus budget
    reduces it — the cost is per sample.

    **Two of the levers are free, in the sense that they cost no accuracy at all** (owner,
    2026-08-13). They come from the estimator's own shape rather than from any approximation:

    - **The generic and repeat-tract halves are never resident together.** The repeat-tract half takes
      one number per sample from the generic half — that sample's homozygote excess, which weights the
      genotype drawn from a locus's length frequencies (§4.1) — and returns nothing to it: `α`, the
      noise classes, the frequency density and the homozygote excess are all fitted on the generic loci
      (§3.3, §4.3). So the order is generic, then repeat tracts, and the generic records are dropped in
      between. At two million positions and 462,701 tracts that is the difference between holding
      6 MB a read group and holding the larger of 1.25 MB and 4–5 MB.
    - **A stratum is fitted on its own.** Slippage is per (read group × stratum), the concentration is
      per stratum, and a locus's length frequencies are drawn from that stratum's own spectrum. Only
      sums of counts cross strata, and a sum is accumulated as the strata go past.
      [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.2 makes the file
      readable that way. **The measured limit:** one stratum holds 47% of tomato's kept tracts, so
      per-stratum reading halves the repeat-tract peak rather than dividing it by 141, and the
      per-stratum cap is what turns it into a bound
      ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.5).

      *The reading order is free of accuracy cost; the cap that turns it into a bound is not, and
      2026-08-13 priced it.* At **5,000 tracts** it costs nothing measurable — every fitted number
      within 2.4% of the truth and moving no more than 2.3% between draws at three reads a site. Below
      1,000 it costs the fall-off of two-repeat slips first and the concentration second
      ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6). So the bound is real
      and it is not free at every setting: **50 kB a sample per stratum section, and the accuracy is
      paid for only if the cap goes below the floor.**

    **Three more levers bound the three axes, and these do trade something:**

    - **Read locus-major rather than sample-major.** Every sample's records file is in genome order
      over the same loci, so a merge across files yields one locus's evidence from every sample with
      nothing else resident. Peak memory becomes `loci held × samples × 5 bits` instead of
      `all loci × samples`: a chunk of a hundred thousand loci at five thousand samples is about
      310 MB, and **the chunk size is the budget knob**. What it costs is one pass over those 30 GB per
      iteration, and §3.3's alternation takes tens of iterations. This is the same shape the cohort
      caller already reads its per-sample files in, so it is a reuse rather than a new mechanism.
    - **Iterate on a resident subsample of loci and sweep the whole set once at the end.** The pooled
      rates are finished at about 320,000 positions
      ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3), so the iterative
      part need never see the rest; the full census is read once, for contamination (§3.4.4).
    - **Subsample samples** — question 8, which bounds the other axis.

    *Leaning:* **take the two free levers always** — nothing is given up by fitting the halves in
    order and the strata one at a time, so they belong in the design rather than in a low-memory mode.
    Of the other three, iterate over a few hundred samples and a few hundred thousand loci held
    resident, with locus-major reading as the mechanism that makes either expressible, then one
    streaming pass over everything for the per-sample numbers. **A cohort of fifty needs none of the
    five and can take the run that never writes a file.**

    **Settled by:** measuring one iteration's wall time and peak resident memory at a thousand and at
    five thousand samples, held against streamed, on records written for the tomato cohort. **Report
    the two separately**: a design that meets the ceiling by re-reading 30 GB thirty times has not
    solved the problem, it has moved it from memory into wall time, and a single "it fits" would hide
    that. **Report peak resident for the generic half and the largest stratum apart from each other**,
    since the free levers make the total meaningless — a run that never holds both has no total to
    report.
11. **Can the duplicated class be identified from the cohort pattern alone, with no coverage summary?**
    — **CLOSED 2026-08-13, and superseded 2026-08-14.** Yes from about twenty-five samples up: the
    pattern returns the inbreeding coefficient at 0.5807 against a truth of 0.5942 where a fit with no
    third class returns 0.4471, and observed heterozygosity 3.0% high against 50.6%. Below that it
    fades — 6.8% high at twenty-five samples, 21.2% at ten — and at one sample it has no power at all.
    **The summary it was weighed against has since been removed entirely** (§2.2, and
    [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4), so what covers the
    small-cohort and single-sample cases is the position's own depth, not a window.
12. **Does the locus's own read count do what the coverage-by-window summary does?** — **CLOSED
    2026-08-13: no, and the gap is widest at tomato's own depth**
    ([`../reports/locus_depth_vs_window_2026-08-13.md`](../reports/locus_depth_vs_window_2026-08-13.md)).
    Read at the width the sample's depth requires, the window gives 14-fold enrichment at 2.5 reads a
    position where the position's own read count gives 1.1-fold, and the cheap arm calls 44 in every
    100 positions it scores "about two copies" against the window's 2.8 (§2.2). **The coverage-by-window
    summary is kept.**

    *Three things the measurement settled beyond the question asked.* The GC correction the cheap arm
    needs does **not** drag the window's machinery in with it — a depth-against-GC curve fitted from one
    position in 300 matches one fitted from all 7.5 million, 10.72 against 10.77 — so the cheap arm
    stayed cheap and simply does not work, which is the clean form of the answer. The **per-position
    depth cap puts a ceiling on any per-position companion**: from 76 reads a position a doubled
    position's stored code no longer sits above an ordinary one, and from 98 the two are written
    identically. That never fires on tomato, whose deepest accession has a median position of 31 reads,
    but it is inside the range this caller commits to, and the window has no such ceiling. And **the
    five-bit encoding costs the discriminator 11% of its enrichment at 13.3× and 37% at 28.7×**, where
    below five reads a position it costs nothing — the first measurement of what per-record binning
    costs a consumer, as against the pooled-cell figure
    ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4).

    *Not adopted:* using the position's own read count as a **second** condition beside the window. It
    would take 9.6× to 16.7× at 9.9 reads a position and 21.3× to 25.7× at 25.2×, and cost 14.4× → 13.0×
    and 6.6× → 5.5× at the shallow end, so it needs a depth-conditional rule the fit does not have.
    **Settled by**, if it is ever wanted: the change in fitted heterozygosity on a drawn panel, not the
    enrichment — enrichment does not say whether a discriminator reaches anything the caller emits.

---

## 12. How we know it works

*The selection rule's tests are [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md)
§7 and the records' are [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §7.
These are the estimator's.*

1. **The parameters fit recovers known parameters.** Draw allele frequencies from a known density, draw genotypes,
   draw reads at known clean and noisy error rates with a known noisy-locus fraction, and fill the
   records directly — no reads, no alignments. The parameters fit must return every drawn value: the three noise
   numbers per read group, the density's four, and each sample's homozygote excess. **At `P = 2` and
   `P = 4`**, since the likelihood is written for any ploidy and an untested loop bound is an
   assumption. **Draw from the fitted family first and from outside it second** — a fit graded only on
   data drawn from its own shape reports that its arithmetic is right and nothing about its model.
2. **The derived rates are unbiased.** §3.2 derives `Hobs` and `π_hom_alt` from posteriors rather than
   fitting them. On the same synthetic draw, both must match the drawn values, and must match what the
   per-sample route fits from the identical genotypes. **This is the test that would catch the kept set
   being a biased subsample**, and it is the criterion [`parameter_prepass.md`](parameter_prepass.md)
   §4.1 says cannot be waived.
3. **The two inbreeding coefficients behave as §5 says.** Three synthetic cohorts: one outbred and unstructured (both must
   return ~0), one with two subpopulations and no autozygosity (**homozygote excess up, autozygosity
   ~0**), and one with a known autozygous fraction and no structure (both up, and agreeing). A route
   that returns the same number in all three has not implemented two estimators.
4. **Adding a false-heterozygote floor moves the two in opposite directions.** Add spurious heterozygous
   sites uniformly at up to five times the real rate. The autozygosity estimate must not move (it is
   measured not to, at that floor) and the homozygote excess must fall.
5. **The panel is not secretly required to be large — MEASURED 2026-08-12, and it is required to be
   larger than anything here had said.** Refit at 1, 2, 5, 10, 25 and 50 samples on the same drawn
   genomes and report where each parameter stops being estimable.

   | samples | clean rate | noisy rate | noisy-locus share | `F_hom_excess`, truth 0.6 |
   |---:|---:|---:|---:|---:|
   | 1 | −35.0% | −62.6% | +474% | 0.530 |
   | 5 | −21.7% | −59.0% | +383% | 0.607 |
   | 10 | −8.9% | −29.6% | +73% | 0.640 |
   | 25 | **−1.4%** | **−0.9%** | **+5.4%** | 0.613 |
   | 50 | −0.7% | +3.2% | −2.3% | 0.601 |

   **The two-class noise model is what needs the panel, and it needs about twenty-five samples.**
   Below ten, the noisy class's share comes out four to five times the truth in the same direction at
   every size — the flat ridge §3.3 warns about, with too little data to pin it, and not scatter.
   **Inbreeding needs about five**, and at one sample it is not identified at all: 0.097 where the
   truth is zero and 0.530 where the truth is 0.6, pulled towards the middle from both sides exactly
   as §6.1 says it must be. *(`Hexp` is deliberately not in this table: at the 40,000 loci the sweep
   used, its scatter is set by the locus budget rather than by the panel, and the budget is
   [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3's experiment.)*

   **What it obliges the route to do**: emit `F_hom_excess` as *not identified* at one sample, and
   emit the noisy class's two numbers with a provenance saying the panel was too small below whatever
   floor a further sweep settles. A fitted `w` five times the truth is not a degraded prior — it is a
   background subtraction, and on a sample at tomato's heterozygosity floor the background *is* the
   answer ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2).
6. **Contamination is recovered, and the test that catches a broken frequency is a *spiked* structured
   panel — RUN 2026-08-12, and it turned this item round.** Three cohorts.

   - **Mix a known fraction of one drawn sample's reads into another's** at 1%, 2% and 5% and require
     `α` back. On an unstructured panel a pooled frequency already passes this — 0.0103, 0.0325 and
     0.1065 for truths of 0.010, 0.030 and 0.100 — so it checks the arithmetic and nothing else.
   - **Draw subpopulations with different allele frequencies, contaminate nothing, and require `α` ≈ 0
     in every sample.** An earlier version of this item called this the test that matters more, on the
     grounds that a pooled frequency would inflate `α`. **It does not**: measured, a pooled frequency on
     a clean structured panel returns `α` = 0.0000 in every sample at `F_st` 0.10 and 0.20 — it passes
     this test perfectly, because it has lost its power rather than its calibration
     ([`../reports/joint_contamination_2026-08-12.md`](../reports/joint_contamination_2026-08-12.md) §3).
     Keep it, because a *noisy* frequency does fail it — one estimated from a twelfth of the panel puts
     41 to 47 of 50 clean samples over 1% — but it is no longer the one that catches the failure §3.4.2
     is about.
   - **The one that does: a structured panel with one sample truly contaminated.** At `F_st` 0.20 a
     pooled frequency returns **0.0050 for a true 0.030** and **0.0000 for a true 0.010**, where the
     correct per-individual frequency returns 0.0347 and 0.0131. **Run it at more than one degree of
     divergence**, since the loss grows with it and a barely-structured panel lets a broken
     implementation through.

   **And a fourth, which is not statistical**: a batching that leaves a read group out must be refused
   and must name it (§3.4.3), since the alternative is a wrong contaminant prior for exactly the
   samples the user forgot.
7. **The post-hoc estimate agrees with the pre-pass one where nothing is wrong** (§3.4.5). On the
   synthetic mixture, re-estimating from the called genotypes must land within the pre-pass estimate's
   own error rather than somewhere else. **This is the only step-4 parameter with a second, stronger
   measurement available**, so it is also the only place a silent modelling error in the route would
   show up without a truth set — which is the reason to build it even when no library is suspect.
8. **The parameters fit is deterministic** — same records, same parameters, independent of thread count and of the
   order samples were walked in.
9. **The frequency density's shape is enough, or it is measured not to be** (§11 questions 6 and 7).
   **Run 2026-08-12 for the first half and it closed question 6**: fitted against a drawn cohort whose
   true density has two bumps, the four-number shape returns `Hexp` 4.9% high and a seven-number shape
   with a second bump returns 5.8% high, while the noise rates move about two percentage points
   between the one-bump and two-bump truths. **Keep it as a regression**: it is the only test here
   that grades the *model* rather than the arithmetic, and it must be re-run whenever the density's
   shape is touched.

   **The second half is still to run**: profile the likelihood over `F_hom_excess` and the density's
   mean together at 5, 10 and 50 samples and report the correlation between them; a ridge rather than
   a peak is what says the pair is not separable at that panel size.

   **And a third thing the first run turned up, which belongs in this list in its own right.**
   Report every fit **from at least three starting points, with each one's score**, not just the
   winner. On the identical data the clean error rate came back −1.3%, −8.9% and −45.8% from three
   starts, and the worst was the one that put the two noise classes close together (§11 question 3).
   A suite that reports one number per parameter cannot tell a converged fit from a trapped one.
10. **The alternative allele is summed over and not picked** (§3.1.1). On a drawn cohort that is
    entirely monomorphic — no real variant anywhere, only sequencing error — the fitted
    `p_invariant` must come back at 1 within its own error. **A fit that chooses each site's
    alternative allele as the observed largest of three fails this and passes every other test in
    this list**, because the bias it introduces is small per site, one-directional, and lands only on
    the classes this test isolates.
