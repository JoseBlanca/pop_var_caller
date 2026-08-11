# ng — the joint fit: every parameter estimated once, over every sample at the same loci

*Design spec, 2026-08-10. **No code yet — this settles the design.** **Read this one first.** It says
what the route is, what it produces, why it exists, and how the estimate is computed. Two companion
documents settle the two things it stands on, and each maps to a module someone can build on its own:*

| document | what it settles |
|---|---|
| **this one** | what the route is and what it produces (§1), what having every sample at one locus changes (§2), the estimator (§3–§4), inbreeding (§5), what it cannot reach (§6), and the comparison it exists for (§8) |
| [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) | **which loci every sample keeps evidence at** — the rule, the stratified STR variant and the reference catalog it selects from, the size knobs |
| [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) | **what is recorded at each kept locus** — the two record shapes, the depth ladder, the encoding |

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
sample's evidence at a locus enters one likelihood, so the fit cannot be split by sample, so it runs
once — after every sample has been walked and before any calling begins.

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
| **inbreeding `F`** | runs of homozygosity, per sample | **homozygote excess** against the panel (§5) | **no — §5** |
| the cohort's diversity `Hexp` | — | the fitted per-locus frequencies **directly**, with no division by `1 − F` (§5.3) | yes, and better conditioned |
| STR diversity | — | the STR kept loci, **reweighted by stratum** ([loci](parameter_prepass_joint_loci.md) §3.3) | yes, once reweighted |
| the frequency spectrum | — | **it is one of this fit's own parameters** (§2.1) | yes |
| contamination, relatedness, read-group grouping | — | unchanged: [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5–§7 | yes |

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
   the cohort — that weighting is the **frequency spectrum**, and it is fitted.

So the free parameters are: the noise rates (per read group, and per stratum on the STR path), the
spectrum (one weight per possible allele count in the panel), and one number per sample (§5). **None
of them grows with the number of loci.** A locus contributes evidence and holds no parameter of its
own.

*In the standard vocabulary the per-locus frequency is a latent variable given an estimated prior —
empirical Bayes. This is also, exactly, how the frequency spectrum is estimated from genotype
likelihoods without calling genotypes (ANGSD's `realSFS`,* Bioinformatics *2015), which
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4 already adopts for the spectrum
alone. **Under this route that estimator is not a separate step: it is the inner loop of every fit
here**, and the spectrum comes out of it as a by-product rather than being computed afterwards.*

**What this answers.** [`parameter_prepass.md`](parameter_prepass.md) §4.2 asks whether a prior fitted
against per-locus weights beats one fitted against pooled weights, and records that the question is
askable only where the loci keep their identity. This route *is* that arm of the question, on both
paths at once.

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

---

## 3. The generic path

### 3.1 What is fitted

Free parameters, all cohort-level except the last:

| parameter | grain | count |
|---|---|---|
| clean error rate `ε_clean`, noisy error rate `ε_noisy`, noisy-locus fraction `w` | read group | 3 per read group |
| the frequency spectrum | cohort | one weight per allele count, `2N + 1` for `N` diploid samples |
| homozygote excess `F` | sample | 1 per sample (§5) |

**A locus's likelihood**, for a locus whose reference base is known and whose alternative allele the
gather has chosen ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2):

```text
                        2N                              class ∈ {clean, noisy}
P(locus) =  Σ  spectrum(c) ·  Σ    class_weight ·   Π    Σ   genotype_freq(j | c/2N, F) · P(reads | j, ε_class)
            c                class                sample  j
```

Read outward: a sample's reads at the locus are scored under each genotype `j` and added, weighted by
what an allele frequency of `c/2N` and that sample's `F` imply about how common that genotype is; the
samples multiply because given the frequency they are independent; the locus's error rate is drawn
from one of two classes; and the frequency itself is summed over the fitted spectrum. **The innermost
sum is [`parameter_prepass.md`](parameter_prepass.md) §3's likelihood unchanged** — the same `p_j`,
the same `ε/3`, the same ploidy-generic loop bound.

### 3.2 What is derived rather than fitted

**A sample's heterozygosity and its homozygous-non-reference rate are not free parameters here, and
that is a real difference between the routes.** Once the fit has converged, every sample has a
posterior over its genotype at every kept locus, so:

```text
Hobs(sample)       =  mean over kept loci of  P(genotype is heterozygous | that locus's reads, fitted parameters)
π_hom_alt(sample)  =  mean over kept loci of  P(every copy non-reference | the same)
```

**These are rates over the kept loci, and they estimate the genome-wide rates because the loci were
chosen by a rule that never looks at the data**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §1.2). A subsample drawn without
reference to what is at the position is unbiased for a rate; that is the property the whole route
rests on, and §11.2 checks it.

**Trap: the per-sample route fits these two as free parameters and this one does not.** They are
therefore not two estimates produced the same way, and a disagreement between the routes is not
automatically an error in either. It has a specific likely cause: under this route the two rates are
constrained by the spectrum, so a sample whose real genotype distribution departs from what the
panel's spectrum and its own `F` predict cannot express that departure. §8 says how to tell that
apart from noise.

### 3.3 How the maximum is found

**Alternate, as the per-sample route already does.** [`parameter_prepass_generic.md`](parameter_prepass_generic.md)
§5.1 settled this for two coupled fits and measured it: starting every error rate at three times the
truth and every frequency at half, the loop converged to the truth in all 25 worlds tried. The same
shape applies here, with three blocks instead of two:

1. hold the spectrum and each sample's `F`, and fit each read group's three noise numbers;
2. hold the noise numbers, and climb to the spectrum (the `realSFS`-style expectation-maximization of §2.1);
3. hold both, and fit each sample's `F`;
4. repeat until the fitted values stop moving.

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

---

## 4. The STR path

**Here the route is HipSTR's model, and the obstacle that kept it out of the per-sample design is
gone.** [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3 sets out the two models side by
side and its last row is the whole story: a per-locus stutter model **needs several samples at one
locus**, and a per-sample walk never has them. That document is explicit that the per-stratum choice
"fitted the shape of the pass" rather than being the more accurate model.

**And there is a second thing the per-sample route cannot do that this one can.**
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 records that a per-read tally is exactly
unbiased **when the allele spectrum is handed to it** rather than fitted from the same tally — and
that there is nowhere to get one during the walk, so the cheaper object was deferred rather than
rejected. **The kept STR loci are that spectrum**, fitted per locus. This route is where that deferral
comes due.

**Decision: fit the four slippage numbers per (read group × stratum), as today, but weight the
genotype by the locus's own length spectrum rather than by the stratum's.** The stratification stays —
slippage depends on repeat count more than on anything else, which is also why the loci are chosen
per stratum ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3) — and what
changes is the weighting inside the sum, exactly as on the generic path. The genotype at an STR locus
is a pair of tract lengths, so the locus's "allele frequency" is a distribution over lengths and the
spectrum being fitted is over length classes rather than over allele counts.

**Two of the per-stratum route's mechanisms are still needed and are unchanged**: a thin stratum
borrows from its neighbours, and the fitted level is held monotonic along the repeat-count axis
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3). Nothing about per-locus weighting makes
a stratum with eleven loci in it fittable.

**How much of the genome's STR loci this route holds is the open question that decides its standing**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6). If the per-stratum cap
keeps every locus in most strata, this route holds the same loci as the per-stratum histogram *and*
remembers which was which — a much stronger position than the generic path's one position in a few
hundred.

---

## 5. Inbreeding — two quantities, both called `F`

**This route produces an `F`, and it is not the same quantity the per-sample route produces.** That
has to be stated before anything else, because both are one number per sample between 0 and 1, both
are called the inbreeding coefficient in the literature, and a consumer handed the wrong one gets a
plausible answer.

| | runs of homozygosity (the walk) | homozygote excess (here) |
|---|---|---|
| what it measures | the fraction of **this genome** where both copies descend from one ancestral copy | how much less heterozygous this individual is than random mating in the panel would predict |
| how | walk 100 kb windows, ask each whether it sits inside a long stretch nearly free of heterozygotes | `1 − Hobs/Hexp`, with `Hexp` from the fitted spectrum |
| needs | one sample | the whole panel |
| what the caller's prior asks for | **this one** | — |

**The caller's genotype prior mixes `F·π_i + (1−F)·π_i^ploidy`**
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
  all** — both states' heterozygote rates lift together and `F` reads only the gap — while the
  whole-genome heterozygosity the ratio reads **inflates eight-fold**
  ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
  §3.3).
- **Old inbreeding.** 100 kb windows resolve runs of about 300 kb and longer. Autozygosity from far
  enough back to have been broken into shorter tracts is invisible to runs and still suppresses
  heterozygosity, so the ratio sees it and the runs do not.

### 5.1 Decision: fit it, emit it under its own name, and never as a substitute

Fit one homozygote-excess `F` per sample here, in step 3 of §3.3's alternation, as the departure from
Hardy-Weinberg proportions the spectrum predicts. **Emit it as a distinct, differently-named parameter
from the autozygosity `F`** — never as an alternative value for the same field. A caller must not be
able to receive one where it expected the other.

**Why fit it at all, given the caller wants the other one:**

- **It is the only external check the autozygosity estimate has.** Nothing else in step 4 measures
  inbreeding twice, and the runs estimator's known failure mode is a confident zero
  ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5: a start that guesses the
  inside state's heterozygote rate far below the truth returns `F` = 0.0000 on a genome with 29% of
  its length in runs, and reports convergence).
- **A disagreement between the two is diagnosable rather than confusing**, because the three causes
  above push in known directions: structure raises the excess above the autozygosity, artifacts push it
  below, and old inbreeding raises it.
- **It costs one number per sample inside a fit that is running anyway.**

### 5.2 It is not circular here, and that is new

[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.3 and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3 both carry the same trap: the panel's
expected heterozygosity is computed as `Hobs/(1 − F)`, so taking `F` from the ratio and then computing
the diversity from it returns whatever was assumed. **That loop is what makes the ratio a diagnostic
rather than an estimate in the per-sample design.**

**Under this route the loop opens.** Expected heterozygosity comes from the fitted spectrum — the
per-locus allele frequencies, measured across the panel — and never from `Hobs`. So `Hexp` and `Hobs`
are two independent quantities and their ratio is a measurement.

### 5.3 What this does to the cohort's diversity

**The diversity stops needing `F` at all.** [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)
§3 computes `Hexp = mean over samples of Hobs/(1 − F)`, which inherits every sample's `F` and its
uncertainty, and breaks entirely if the runs estimator returned a confident zero. Here `Hexp` is read
off the fitted spectrum directly. **That is a genuine improvement independent of everything else in
this document**, and it is worth flagging to whoever finishes the comparison: even a route that lost
on every other axis would still be the better source of this one number.

---

## 6. What this route cannot produce

**Anything local.** The kept positions are scattered on purpose so that no two are inherited together
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §2). A 100 kb window holds about
a hundred thousand sites but only a couple of hundred kept ones, so at one heterozygote per kilobase it
carries a **fraction of one** expected heterozygote instead of a hundred — far too thin to tell a run
of homozygosity from ordinary sequence.

So **runs of homozygosity, and the autozygosity `F` the caller's prior actually reads, are not
available by this route and no budget within reach changes that.** The windowed histogram and its
estimator stay in the walk whatever the comparison decides — which means this route can never make the
walk a pure accumulation pass, and the walk's largest object (37 MB per sample on tomato, 145 MB on
human) survives regardless. **That is the honest ceiling on what deleting the per-sample fits could
ever save.**

Linkage, haplotypes, and anything else reading a stretch of genome are out of reach for the same
reason.

---

## 7. Cross-cutting concerns

**The scheduling cost is the one a manager needs, and it is not memory.** This route puts a **barrier**
in the pipeline: no locus can be called until every sample has been walked, because the fit needs every
sample at once. The per-sample route has no such barrier — a sample's parameters are ready when that
sample's walk ends. **Nothing about the kept loci makes this avoidable**; it is what "fitting against
the locus's frequency in the cohort" means. A run that adds one sample later must refit, and the
parameters every earlier call was made under have changed.

**Memory** is [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §5 — roughly
60–110 MB of records for a fifty-sample cohort. **This route adds no per-sample object.** What it adds
is working memory for the fit itself: one posterior over allele counts per locus, `2N + 1` numbers,
which at fifty diploid samples is 101 values and need not be held for more than one locus at a time.

**Compute, honestly.** One evaluation of §3.1's likelihood over the whole generic set is on the order
of `loci × samples × (2N+1)` operations — two million by fifty by a hundred and one, about 10¹⁰ — and
the fit needs tens of iterations of that. **It is embarrassingly parallel over loci** and needs no
communication between them within an iteration, so it is a matter of cores rather than of
feasibility. **This is arithmetic, not a measurement**, and §8's last item replaces it. If it
dominates, the budget is the knob and
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4 is how to set it.

**Determinism.** The fit sums over loci in a fixed order and over samples in a fixed order, so no
parameter varies with thread count. Multiple starting points are enumerated, not sampled.

**Errors.** The eight identity values of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4 must all match across
samples, and **a mismatch must refuse, not average**.

---

## 8. The comparison this route exists for

**Both routes run on the same data and their estimates are compared.** Four measurements, in the order
they can be made.

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
3. **The two-class residual.** §2.2's claim is that fifty samples at one locus separate a mismapped
   locus from a heterozygous one, and the measurement is already set up: refit HG002's 551,843
   confident-region loci and see whether heterozygosity comes in below the 1.09× the blind two-class
   mixture reaches. **This is the measurement most likely to decide the comparison**, because it is the
   one thing this route can do that no amount of extra data can buy the other.
4. **Memory and wall clock, measured rather than computed.** §7's figures are arithmetic. Report the
   records at rest, the fit's working set, and the fit's wall clock at several core counts, on the
   whole tomato cohort — the run that stresses sample count, which is the axis this route's cost
   scales on. **Report it at each budget of
   [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.1's downward sweep**, since
   the two questions share one experiment.

**What must be recorded so that nobody re-runs the walk to finish this**: both routes' fitted values,
the evidence count behind each, the memory each cost, and the wall clock, written up in
`doc/devel/ng/reports/`.

**Two comparisons that are not like-for-like, and reporting them as if they were is the mistake to
avoid.** `F` is two different quantities (§5) and must be reported as two rows, never as one row with
two values. `Hobs` and `π_hom_alt` are free parameters in one route and derived quantities in the
other (§3.2), so their disagreement has a specific candidate cause and should be read against it.

---

## 9. Reuse over rewrite

| what | existing code | how it is reused |
|---|---|---|
| the sum over genotypes at one sample's locus | `src/ng/parameter_estimation/generic/noise_model.rs` (`NoiseModel`) | used as the innermost term of §3.1, unchanged. It is already the seam both paths implement |
| the climb over mixture weights | `src/ng/parameter_estimation/fitting/mixture_weights.rs` | the spectrum is a mixture over allele counts, so the same climb applies with a different component set |
| the two site classes | `src/ng/parameter_estimation/generic/` (`SiteNoise`) | the parameters are the same three; what changes is that the class becomes a per-locus latent variable (§2.2) |
| alternating between coupled blocks | `src/ng/parameter_estimation/generic/coupled_fit.rs` | the same loop shape and the same termination handling, with three blocks instead of two (§3.3) |
| provenance and evidence count | `src/ng/parameter_estimation/mod.rs` (`Provenance`, `Estimate`) | used as-is; every parameter here carries both, as in the per-sample route |
| the frequency spectrum by expectation-maximization | [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4 | **not a separate step under this route** — it is §2.1's inner sum, and the spectrum falls out of the same fit |

*The two companion documents carry their own reuse: the depth ladder
([`records`](parameter_prepass_joint_records.md) §2.2) and the even spread across a stratum
([`loci`](parameter_prepass_joint_loci.md) §3.1), both already written.*

**No parity oracle.** Neither route is a port of the other, and agreement is the thing being measured
rather than the thing being asserted (§8).

---

## 10. Deferred, with a recommended home

- **Runs of homozygosity, and the autozygosity `F`.** Not deferred so much as unreachable here (§6).
  **Home:** [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6, where it already lives
  and where it stays.
- **Contamination, relatedness, and read-group grouping.** Already specified and unchanged by this
  route, which reads the same records they do. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5–§7.
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
3. **How many starting points, and spanning what?** — OPEN. The per-sample route's answer for its own
   two-state model is nine starts spanning the separation between the states
   ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5), chosen after a start that
   guessed wrong returned a confident zero. The analogous separation here is between the clean and
   noisy classes. *Leaning:* the same shape of answer; the number needs the curve from question 2.

---

## 12. How we know it works

*The selection rule's tests are [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md)
§7 and the records' are [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.
These are the estimator's.*

1. **The fit recovers known parameters.** Draw allele frequencies from a known spectrum, draw genotypes,
   draw reads at known clean and noisy error rates with a known noisy-locus fraction, and fill the
   records directly — no reads, no alignments. The fit must return every drawn value: the three noise
   numbers per read group, the spectrum, and each sample's homozygote excess. **At `P = 2` and
   `P = 4`**, since the likelihood is written for any ploidy and an untested loop bound is an
   assumption.
2. **The derived rates are unbiased.** §3.2 derives `Hobs` and `π_hom_alt` from posteriors rather than
   fitting them. On the same synthetic draw, both must match the drawn values, and must match what the
   per-sample route fits from the identical genotypes. **This is the test that would catch the kept set
   being a biased subsample**, and it is the criterion [`parameter_prepass.md`](parameter_prepass.md)
   §4.1 says cannot be waived.
3. **The two `F`s behave as §5 says.** Three synthetic cohorts: one outbred and unstructured (both must
   return ~0), one with two subpopulations and no autozygosity (**homozygote excess up, autozygosity
   ~0**), and one with a known autozygous fraction and no structure (both up, and agreeing). A route
   that returns the same number in all three has not implemented two estimators.
4. **Adding a false-heterozygote floor moves the two in opposite directions.** Add spurious heterozygous
   sites uniformly at up to five times the real rate. The autozygosity estimate must not move (it is
   measured not to, at that floor) and the homozygote excess must fall.
5. **The panel is not secretly required to be large.** Run the fit at 2, 5, 10 and 50 samples on the
   same drawn genomes and report where each parameter stops being estimable. The spectrum has `2N + 1`
   weights, so at two samples it has five, and a route whose parameters are fitted against it will
   degrade somewhere. **Where** is what a user needs to be told, and nothing currently says.
6. **The fit is deterministic** — same records, same parameters, independent of thread count and of the
   order samples were walked in.
