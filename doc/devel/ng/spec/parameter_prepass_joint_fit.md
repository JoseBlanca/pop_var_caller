# ng — the joint fit: every parameter estimated once, over every sample at the same loci

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
sample's evidence at a locus enters one likelihood, so the fit cannot be split by sample, so it runs
once — after every sample has been walked and before any calling begins.

**What this route reads, in one line: the records of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md), at the kept loci, and
nothing else.** One entry per kept position holding that position's allele counts and its binned
depth, per read group; one entry per kept STR locus holding its offsets, guard and differences. It
accumulates no summary of its own — a summary has forgotten which locus it saw, which is the property
§2 exists to exploit. **The one object of the other route that survives beside it is not read by
anything here**: the autozygosity `F_autozygosity` needs runs of homozygosity and scattered loci cannot give them
at any budget, so that estimator stays in the walk whatever the comparison decides (§6). The two
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
| the cohort's diversity `Hexp` | — | the fitted per-locus frequencies **directly**, with no division by `1 − F` (§5.3) | yes, and better conditioned |
| STR diversity | — | the STR kept loci, **reweighted by stratum** ([loci](parameter_prepass_joint_loci.md) §3.3) | yes, once reweighted |
| the frequency spectrum | — | **it is one of this fit's own parameters** (§2.1) | yes |
| **contamination `α`** | — **it cannot** ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2) | the generic kept loci, against individual-specific allele frequencies this fit derives (§3.4) | it is only produced here |
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

**And the two classes are not enough.** `ε_noisy` lives on the error-rate ladder, which stops at 10%,
while the loci that ask for a class of their own sit at about **half** — and
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1 is right to refuse to widen the
ladder for them: as a noisy rate approaches a half, a noisy site and a heterozygous site are the same
distribution, so a class stretched that far would take real heterozygotes with it. **What those loci
are is a duplication the reference does not carry**, collecting two copies' reads at one position, so
every position where the copies differ shows alternative reads at every depth. Two of five real
alignments ask for such a population — tomato SRR7279482 and SRR7279483, at **0.42% and 0.49% of
sites** — and the model refuses, leaving the only class it has for them: *heterozygous*.

**How many loci that is, on the sample where it matters most.** If the fitted share is right the
population is about **8,400 of two million** kept positions; if instead the share is the
alternative-read mass a 10% rate had to carry and the true rate is 50%, it is about **1,700**. Against
either, tomato's least heterozygous sample has **300** genuinely heterozygous positions in two million
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2) — so the artefact
population is between six and thirty times the quantity being estimated, on exactly the samples this
caller is aimed at.

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

**Decision: a third class of site, and its discriminator is local relative coverage — not the
alternative-read fraction.** A duplicated locus and a heterozygous one both read about half
alternative, so the read counts at one locus cannot separate them; what differs is that a duplication
collects two copies' reads. Production settled the same question for the same artefact:
[`../../specs/hidden_paralog_filter.md`](../../specs/hidden_paralog_filter.md) §2 makes coverage the
primary signal, and it is also the only one that leaves an introgression alone, being single-copy and
so of normal depth. **A site whose window sits near two copies is drawn from a class at about half
alternative reads**, and the class's weight is fitted as `w` is.

**It must be conditioned on the window, not on the site**, which is the same document's measured
constraint: per-base coverage at 6× has no power to tell a two-copy carrier at ~12 reads from a
single-copy sample reading high, and tomato's three reads a site is half that depth. **The cost is a
per-sample coverage-by-window summary, GC-corrected** — 500 bp windows over tomato's 800 Mb is 1.6 M
windows at a few bytes, single-digit megabytes per sample, plus a small GC curve. Production computes
both in its Stage-1 pileup; `src/pileup/` is frozen, so ng builds its own. **Nothing in it needs a
cohort**, which matters because this caller must also run on one sample.

**What many samples at one locus add is identification, and it is a bonus rather than the mechanism.**
A real half-frequency variant leaves about a quarter of the samples in each homozygous class; a
duplication leaves every sample at a half. So the cohort names *which* loci they are instead of merely
weighing how many there are.

**Open, and the measurement comes before any more design**: on one tomato sample, the fraction of
positions that sit in windows near two copies **and** show a near-half alternative fraction. One walk
over one alignment, no cohort and no truth set, and it is what says whether the class is worth an
accumulator. *This finding was made on the histogram route and is recorded here because this is the
route being planned; that route is implemented and is not being changed for it.*

---

## 3. The generic path

### 3.1 What is fitted

Free parameters, all cohort-level except the last:

| parameter | grain | count |
|---|---|---|
| clean error rate `ε_clean`, noisy error rate `ε_noisy`, noisy-locus fraction `w` | read group | 3 per read group |
| the frequency spectrum | cohort | one weight per allele count, `2N + 1` for `N` diploid samples |
| homozygote excess `F_hom_excess` | sample | 1 per sample (§5) |

**A locus's likelihood**, for a locus whose reference base is known and whose alternative allele the
gather has chosen ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2):

```text
                        2N                              class ∈ {clean, noisy}
P(locus) =  Σ  spectrum(c) ·  Σ    class_weight ·   Π    Σ   genotype_freq(j | c/2N, F) · P(reads | j, ε_class)
            c                class                sample  j
```

Read outward: a sample's reads at the locus are scored under each genotype `j` and added, weighted by
what an allele frequency of `c/2N` and that sample's `F_hom_excess` imply about how common that genotype is; the
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
rests on, and §12.2 checks it.

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

1. hold the spectrum and each sample's `F_hom_excess`, and fit each read group's three noise numbers;
2. hold the noise numbers, and climb to the spectrum (the `realSFS`-style expectation-maximization of §2.1);
3. hold both, and fit each sample's `F_hom_excess`;
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
whether a locus is contaminated; nothing here does that, and §3.4.3 is about what the identification
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
loci: for each sample, how much of its low-fraction alternative-read mass falls on loci the fit says
are segregating, carrying the alleles it says segregate there.

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
that individual's coordinates in a principal-component space, and the fit **jointly maximises over
`α`, the intended sample's coordinates and the contaminant's** — four components by default, with an
equal-ancestry and an unequal-ancestry model compared by AIC. It removes about 80% of the bias.

**The tomato panel has exactly the structure that breaks the pooled version.** It is landraces from
several regions, which is why §5 rejects homozygote excess as a measure of autozygosity — the Wahlund
effect. **A pooled spectrum used as the contaminant's frequency makes the same error in a different
place**: a sample from a diverged subpopulation carries alleles the pooled spectrum calls rare, and
rare alleles turning up in a sample is exactly the contamination signature. Structure would be read as
contamination.

**Decision: derive the individual-specific frequencies from this cohort, not from an external panel.**
There is no tomato HGDP; what there is, is the kept loci in every sample, which is the matrix a
principal-component decomposition wants — and the same matrix
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6 already builds for relatedness. So the
components come free with an object this route holds anyway. *Soft:* how many components a plant panel
needs is not four by inheritance, and it is set by the panel rather than by us.

### 3.4.3 Depth, marker count, and what that does to the budget

**Correction to an earlier draft of this section, which said `α` would not be identified at three reads
a site.** That was the wrong axis. No single site is classified: the information is pooled across
markers, so what governs is **how many markers segregate**, with depth entering only because a marker
with two or more reads separates a heterozygote from a homozygote better than one with a single read.
`verifyBamID2` runs on **4× genomes**, well below the depth at which any one site is legible.

**How the accuracy moves with the marker count**, as they report it — mean squared error of `α` at
1%, 2% and 5% contamination:

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

### 3.4.4 Measure it again from the final calls, and report the two

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
hundred. **That question is now one call rather than an experiment**, `count_loci_per_stratum` on the
repeat catalog at the STR path's calling floors, and it is the cheapest thing outstanding in these
three documents.

### 4.1 Contamination is not fitted from the STR loci, and the reason is sharper than "noisier"

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
| produced by | the per-sample walk ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6) | this route, §5.1 |
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
Hardy-Weinberg proportions the spectrum predicts. **Emit it as a distinct, differently-named parameter
from the autozygosity `F_autozygosity`** — never as an alternative value for the same field. A caller must not be
able to receive one where it expected the other.

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

### 5.2 It is not circular here, and that is new

[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.3 and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3 both carry the same trap: the panel's
expected heterozygosity is computed as `Hobs/(1 − F_hom_excess)`, so taking `F_hom_excess` from the ratio and then computing
the diversity from it returns whatever was assumed. **That loop is what makes the ratio a diagnostic
rather than an estimate in the per-sample design.**

**Under this route the loop opens.** Expected heterozygosity comes from the fitted spectrum — the
per-locus allele frequencies, measured across the panel — and never from `Hobs`. So `Hexp` and `Hobs`
are two independent quantities and their ratio is a measurement.

### 5.3 What this does to the cohort's diversity

**The diversity stops needing an inbreeding coefficient at all.** [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)
§3 computes `Hexp = mean over samples of Hobs/(1 − F)`, which inherits every sample's `F_autozygosity` and its
uncertainty, and breaks entirely if the runs estimator returned a confident zero. Here `Hexp` is read
off the fitted spectrum directly. **That is a genuine improvement independent of everything else in
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
estimator stay in the walk whatever the comparison decides — which means this route can never make the
walk a pure accumulation pass, and the walk's largest object (37 MB per sample on tomato, 145 MB on
human) survives regardless. **That is the honest ceiling on what deleting the per-sample fits could
ever save.**

**Linkage, haplotypes, and anything else reading a stretch of genome are out of reach for a different
reason, and it is a matter of shape rather than of budget.** A uniform thinning throws away exactly
the close pairs where linkage information lives and keeps unlimited distant ones, which carry none;
the instrument for those is a small clustered budget
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §8), which nothing in step 4
asks for today.

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

**Errors.** The ten identity values of
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4 must all match across
samples, and **a mismatch must refuse, not average**. Nine of them say the samples were asked for the
same loci; the tenth is a digest of the loci each one actually kept, and it is the only one that would
catch a hash function or a threshold's arithmetic changing underneath two otherwise identical runs.

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
  which case the fit is being graded against an assumption. **So measurement 3 below is deliberately
  not synthetic**, and it is the one that grades the model rather than the arithmetic.
- **Memory and speed are only comparable at realistic scale.** The numbers that decide this route —
  the records held for every sample at once, and the per-sample windowed histogram the other route
  holds in flight — do not appear on a small fixture. Run the cost measurements at the shape of a real
  cohort: a genome-sized position budget, fifty samples, and the concurrency the walk would actually
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
3. **The two-class residual.** §2.2's claim is that fifty samples at one locus separate a mismapped
   locus from a heterozygous one, and the measurement is already set up: refit HG002's 551,843
   confident-region loci and see whether heterozygosity comes in below the 1.09× the blind two-class
   mixture reaches. **This is the measurement most likely to decide the comparison**, because it is the
   one thing this route can do that no amount of extra data can buy the other.
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
   records at rest, the fit's working set, and the fit's wall clock at several core counts, on the
   whole tomato cohort — the run that stresses sample count, which is the axis this route's cost
   scales on. **Report it at each budget of
   [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3's downward sweep**, since
   the two questions share one experiment.

**What must be recorded so that nobody re-runs the walk to finish this**: both routes' fitted values,
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
3. **How many starting points, and spanning what?** — OPEN. The per-sample route's answer for its own
   two-state model is nine starts spanning the separation between the states
   ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.5), chosen after a start that
   guessed wrong returned a confident zero. The analogous separation here is between the clean and
   noisy classes. *Leaning:* the same shape of answer; the number needs the curve from question 2.
4. **How many principal components does a tomato panel need for the individual-specific allele
   frequencies?** — OPEN (§3.4.2). `verifyBamID2` defaults to four on human reference panels, and a
   crop panel's structure is not human structure: a landrace collection may need more, and an
   inbred-line panel behaves differently again. *Leaning:* decide it from the panel by the usual
   scree-plot criterion rather than by inheriting four, and record what was used beside `α`, since a
   number fitted under too few components is biased in the direction §3.4.2 describes — structure read
   as contamination. **Settled by:** the components come out of the same matrix relatedness uses
   ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6), so this is a plot rather than an
   experiment.
5. **Do the short-tract STR strata contribute to contamination after all?** — OPEN (§4.1). *Leaning:*
   possibly, and it does not matter much, since the generic loci supply `α` either way. **Settled by:**
   two per-stratum counts — what fraction of a stratum's loci segregate across the cohort, and what
   fraction of its stutter products land on a length the cohort segregates at.

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
3. **The two inbreeding coefficients behave as §5 says.** Three synthetic cohorts: one outbred and unstructured (both must
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
6. **Contamination is recovered, and — the test that matters more — an uncontaminated structured panel
   returns zero.** Two synthetic cohorts. In the first, mix a known fraction of one drawn sample's
   reads into another's at 1%, 2% and 5% and require `α` back; that only checks the arithmetic. In the
   second, **draw two subpopulations with different allele frequencies, contaminate nothing, and
   require `α` ≈ 0 in every sample.** A fit using the pooled spectrum as the contaminant's frequency
   fails this and passes the first, which is exactly the failure §3.4.2 exists to prevent and the one
   `verifyBamID2` measured as 2.9% reported for a true 10%. **Run the second at more than one degree of
   divergence**, since the bias grows with it and a barely-structured panel would let a broken
   implementation through.
7. **The post-hoc estimate agrees with the pre-pass one where nothing is wrong** (§3.4.4). On the
   synthetic mixture, re-estimating from the called genotypes must land within the pre-pass estimate's
   own error rather than somewhere else. **This is the only step-4 parameter with a second, stronger
   measurement available**, so it is also the only place a silent modelling error in the route would
   show up without a truth set — which is the reason to build it even when no library is suspect.
8. **The fit is deterministic** — same records, same parameters, independent of thread count and of the
   order samples were walked in.
