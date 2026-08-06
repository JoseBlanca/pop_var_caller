# ng — the parameter pre-pass: what it produces, and the machinery both paths share

*Design spec, 2026-08-03. **No code yet — this settles the design.** Companion arch and plan docs
do not exist yet. This is ng step 4 of [`ng_proposal.md`](ng_proposal.md), and it **reverses that
step's central mechanism**: no rough caller, no confident-genotype gate, so `SampleSummarizer` is
fed sufficient statistics rather than calls. `CohortEstimator` survives unchanged. §2.3 says what
that means for the two documents that assume otherwise.*

*ng step 4 is specified in **five documents**, and the split is about **length and modularity, not
about passes**. One document holding all of it would be unreadable; and there is a lot of code here,
which — although every piece of it reads the same data — divides cleanly. Each document below covers
**one accumulator and everything fitted from it**, so each maps to a module with its own interface.
Splitting the prose splits neither the traversal nor the data. **This one says what the step produces
and carries the machinery the others share**, and is the one to read first.*

| document | what it settles |
|---|---|
| **this one** | the parameters and their grains (§1), why production is biased (§2), the estimator (§3), what is still open (§4), and the map of what the walk accumulates (§5) |
| [`parameter_prepass_generic.md`](parameter_prepass_generic.md) | the SNP/indel path, end to end: both histograms, the error rate, the sample's rates, and the inbreeding coefficient `F` |
| [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) | the STR path, end to end: its noise model, what stutter actually looks like, its accumulator and how the strata are cut |
| [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) | the two **censuses** — small sets of loci, one over ordinary sites and one over STR tracts, kept raw and identical in every sample so that answers can be compared across them |
| [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) | the gather: diversity, frequency spectrum, contamination, relatedness, sample grouping — run later, inside cohort variant calling |

*The first four all describe that **one** walk over the reads — four accumulators filled in a single
traversal, not four passes. The fifth is the exception and the only one: the cohort gather runs
later, reads the summaries the walk produced, and never opens an alignment file (§1.3).*

*Grounded in the research note*
[`rough_caller_alternatives_2026-07-23.md`](../../reports/research/rough_caller_alternatives_2026-07-23.md)
*and in measurements made on the tomato and HG002 cohorts during the STR stutter study.
`src/ssr/` and `src/pileup/` are frozen production: everything said about them here is a record,
not a change.*

---

## 1. What this step produces

Before the real caller runs, something has to tell it what to expect from this data — both **the
noise in it** and **the variation in it**: how often a base is read wrong, how often a repeat tract
gains or loses whole copies, how often an individual's two copies of a site differ, and how often an
individual differs from the reference at all.

**Here is the whole output of step 4.** Everything else in these five documents is how one of these
rows is produced, or why it is produced that way. **The "specified in" column also says which pass
computes it:** the first seven come out of the walk, one sample at a time; the last five need every
sample and are computed later, at the gather (§1.3). This document's walk *accumulates* what those
five need and computes none of them (§8).

| parameter | grain | fitted from | specified in |
|---|---|---|---|
| **per-base error rate**, generic path | read group | the read-group histogram | [generic](parameter_prepass_generic.md) §3 |
| **per-base error rate**, STR path | read group × stratum | the STR table's composition counts | [ssr](parameter_prepass_ssr.md) §1.1, §4.1 |
| **how often an STR read slips at all** | read group × stratum | the STR table | [ssr](parameter_prepass_ssr.md) §4 |
| **which way an STR read slips** | read group × stratum | the STR table | [ssr](parameter_prepass_ssr.md) §3 |
| **how far it slips when it does** | read group × stratum | the STR table | [ssr](parameter_prepass_ssr.md) §3 |
| **observed heterozygosity** `Hobs` | sample | the windowed histogram, summed over windows | [generic](parameter_prepass_generic.md) §5 |
| **homozygous-non-reference rate** `π_hom_alt` | sample *and reference* | the same, summed over windows | [generic](parameter_prepass_generic.md) §5 |
| **inbreeding coefficient** `F` | sample | the same, window by window | [generic](parameter_prepass_generic.md) §6 |
| **the cohort's diversity** `Hexp` | cohort | the census sites, with every sample's `Hobs` and `F` | [cohort](parameter_prepass_cohort.md) §3 |
| **STR diversity** | cohort | the STR census | [cohort](parameter_prepass_cohort.md) §3 |
| **the frequency spectrum** | cohort | the census sites | [cohort](parameter_prepass_cohort.md) §4 |
| **contamination** | read group | the census sites | [cohort](parameter_prepass_cohort.md) §5 |
| **relatedness** | sample pairs | the census sites | [cohort](parameter_prepass_cohort.md) §6 |
| **which read groups share a chemistry** | read groups | the fitted parameters, with the evidence behind each | [cohort](parameter_prepass_cohort.md) §7 |

**Every parameter also carries where it came from** — fitted here, borrowed from a neighbour, or
defaulted — and **how much data stood behind it** (§6). A downstream model that cannot tell a stutter
rate measured on 80,000 reads from one measured on 8 will treat them alike, and one of them deserves
it. **No uncertainty interval is emitted**; §6 says why not.

### 1.1 Three grains, and why the table has three of them

The grain column is not bookkeeping. It is the claim that decides which accumulator each parameter
comes out of, and getting it wrong is the failure §5 is built to prevent.

- **Noise is chemistry.** The per-base error rate and the STR slippage behaviour are properties of
  how the DNA was prepared and sequenced, not of the genome — so their unit is the **read group**.
- **Variation is biology, and there are two kinds of it.** *Within* an individual: how often its two
  copies differ, and the inbreeding coefficient `F` that says how much of the genome is homozygous
  by descent. *Against the reference*: how often the individual carries a non-reference allele at
  all — a different thing, because a selfing landrace can be almost entirely homozygous and still
  differ from the reference accession at a great many sites. Both take the **sample** as their unit,
  though the second is really a distance between the sample and whichever genome was chosen as the
  reference.
- **The cohort's variation is a third thing, and no single sample holds it.** How variable the
  population is, whether a library is contaminated, which samples are relatives: each is a
  *comparison between* samples, so no amount of data about any one of them contains the answer.

**A sample sequenced from two libraries therefore has two error rates and one heterozygosity.** That
sentence is the whole content of §5, and the reason there are two generic histograms rather than one.

### 1.2 Goals

1. Estimate every row of §1's table from Stage-1 data alone, each at the grain it belongs to.
2. Do it **without first calling genotypes and keeping the confident ones**, which is what production
   does and what biases it (§2).
3. Share the **estimation machinery** across both paths, SNP/indel and STR — marginalise the genotype
   away, then maximise what is left, both defined in §3 — while keeping **two noise models**, because
   the two paths are noisy in different ways (§3.2). *How* the maximisation is done is §3.1.
4. **Compute the inbreeding coefficient `F`** here rather than deferring it, because the cohort's
   diversity divides by `1 − F` and nothing downstream can start without it.
5. **Gather in the same pass** what the cross-sample parameters are computed from later, so that none
   of them needs a second traversal of the reads.
6. Emit parameters the caller can consume as frozen inputs, so genotyping stays a pure function of
   (reads, parameters).

**Non-goals.**

- **Genotyping.** Nothing here calls a variant. The estimator sums over genotypes precisely so it
  never has to pick one.
- **Designing the caller's priors.** This step estimates quantities. What a caller builds out of
  them is the caller's design.
- **Long-allele recovery, and the sample profiling it needs.** GangSTR genotypes alleles longer than
  a read by counting reads that fall entirely inside the tract, which requires profiling the
  library's fragment-length distribution and coverage first. It is a real alternative, deferred on a
  property of our data rather than rejected on merit: the count is Poisson, so its precision is
  bought with coverage, and the tomato cohort has about 3 reads per plant.
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §7 works the numbers and names the error
  that would have to be measured before adopting it. Deferred (§8).

**It does not:**

- decide which loci exist (that is region typing, step 3);
- change anything in `src/ssr/` or `src/pileup/` — production is frozen;
- **assume diploidy.** The likelihood and the accumulators are written for any ploidy (§3), which
  costs a loop bound now and would cost a rewrite later. What is still diploid-only is the
  *population genetics* on top, and that is deferred with a home (§8), not baked in;
- estimate **chemistry** per sample, or **biology** per read group (§1.1, §5).

### 1.3 Two passes, and where each fit runs

**The step runs as two passes over data, and only the first touches reads.**

1. **The walk, per sample.** Read that sample once, accumulate the five objects of §5, and **fit
   everything that sample's own accumulators can answer before moving on** — the error rate, its
   two genotype frequencies, `F`, the stutter parameters. Samples need nothing from each other, so
   this is parallel across the cohort.
2. **The gather, once, over every sample's summary**
   ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)) — **invoked later, inside cohort
   variant calling**, not at the end of the walk. It reads summaries and never reads a read.

**The histograms do not outlive the sample that filled them**. Each is
reduced to its parameters at the end of that sample's walk, then dropped. **Only
two things travel onward**: those estimates, and the census sites — which stay raw, because they hold
the cross-sample correspondence and reducing them per sample is exactly what would destroy it. So the
peak cost of the histograms is *one sample's, times however many samples are walked at once*, not
fifty samples' held together (§6).

The interfaces have this shape — `SampleSummarizer` then `CohortEstimator`
([`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) lines 343 and 351) — and §2.3 says
what changes about the first of them.

**Where each fit runs follows from which object it reads**, so none of it is in doubt. What *is* open
is which of two estimates of the middle four survives — and until that is measured, both are computed:

| parameter | computed | open? |
|---|---|---|
| inbreeding `F` | after that sample's walk — it needs one sample's windowed histogram and nothing else | no |
| diversity, frequency spectrum, contamination, relatedness, read-group grouping | only at the gather — each is a comparison between samples | no |
| error rate, heterozygosity, distance from the reference | **twice, for now**: after the walk from the histograms, and again at the gather from the censuses | **which estimate is kept — §4.1** |
| stutter | once, after the walk, from the STR histogram | no. Whether the STR census could *also* fit it is a different question (§4.2), and §10.2 excludes stutter from §4.1's comparison for that reason |

**What this means for the coder: both routes are built, and both run.** §4.1 is a comparison, and a
comparison needs two numbers. So the histogram fits run after each sample's walk, over that sample's
own accumulators; the census fits run at the gather, over data that exists nowhere else. They are
separate code reading separate objects, and neither waits on the other. **When the measurements are
in, one of the two is deleted** — that is a decision about which code to keep, not about how to
arrange it.

**Where the walk's output lives between here and its consumers is not this step's decision.** It may
be held in memory, serialized per sample, or folded into whatever the pipeline already writes. What
this spec requires is a property, not a mechanism: **the gather must be able to read it without
walking the reads again.** Everything §5 accumulates is shaped by that one requirement.

**The one cross-document constraint, worth knowing in advance.** The gather needs nothing from the
walk but its output. In the other direction there is a single dependency: the gather divides by
`1 − F`, and the simplest estimator of `F` would need that same quantity to start from, so the two
cannot be paired. [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6 explains it where
the reader first needs it, and it is why that document builds the runs-of-homozygosity estimator
rather than the ratio.

**All of it is gathered in one pass over the reads, because the pass is the cost.** Reading and
decoding dominates; the fitting afterwards is arithmetic over accumulators small enough to hold in
memory (§6). A second traversal to pick up `F`, or a third for the cohort's diversity, would multiply
the expensive half to gain nothing. That is why the **accumulators** are settled here even where the
**estimators** that read them are still open (§4).

**Separately, and for a different reason, the noise and the variation have to be *fitted* together.**
A read carrying the alternative allele is either an error or a real second allele, and nothing in
the data says which until both rates are estimated together (§3). **That is a statement about the
model, and it does not settle the procedure.** The two are read from two different tables — the error
rate from the read-group histogram, the frequencies from the windowed one — because they are counted
over different things (§5.1). How the coupling between them is resolved across those two tables is
open, with alternating fits as the leaning:
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §5.1.

---

## 2. Why production's numbers are biased

Two problems — what gets thresholded away, and what never gets looked at — then what they cost, and
what fixing them changes elsewhere in ng.

### 2.1 Both rough callers threshold, then count

The generic path classifies each site as het, hom-alt or ambiguous, then forms
`Hobs = n_het / (n_het + n_hom_alt)` from two of those three counters
([`src/sample_summary/het.rs:266`](../../../../src/sample_summary/het.rs), `observe_site`). The STR
path pools reads from loci that passed a confident-genotype gate
([`src/ssr/cohort/prepass.rs`](../../../../src/ssr/cohort/prepass.rs)). Both estimate a parameter
from **the sites that were easy to call**, and inherit whatever the threshold selected for — which
sites those are depends on depth, so the bias is depth-dependent.

**And the sites with no alternative allele never arrive at all.** The classifier's input builder
returns `None` for a pure-reference column
([`het.rs:147-148`](../../../../src/sample_summary/het.rs)), so a site with depth 30 and no alt
reads is never seen. **The record itself is still staged and written** — the skip is the het
accumulator's alone
([`pileup_to_psp.rs:90-108`](../../../../src/pileup/per_sample/pileup_to_psp.rs)) — so this is a loss
in the estimator, not in the data. Those columns are the majority of the genome and the strongest
evidence there is about `ε`: a clean 30-read column says far more about the error rate than a
marginal 2-alt one does.
This is the loss that shapes ng's accumulator rather than merely explaining production's number:
§5 keeps them, and §3's likelihood is what makes them count.

### 2.2 What that costs, measured

The STR stutter study put a number on it. Slippage was measured against each sample's own modal
allele, on HG002, whose truth genotypes are known from the GIAB assembly-based benchmark. **Eleven
in every hundred loci inside the truth regions are heterozygous**, and a heterozygous locus has a
second real allele whose reads are indistinguishable from slippage unless the genotype is handled.

| fitted on | slippage rate at ≥6 repeats | how much more often reads lose a repeat than gain one (dinucleotide) |
|---|---:|---:|
| all loci | 4.9% | 0.9× — *gains marginally ahead* |
| known-homozygous loci only | 2.0% | 3.4× — *losses well ahead* |

Ignoring the genotype inflates the rate **2.4-fold and reverses the direction**. A stutter model
fitted on the uncontrolled numbers would be wrong in size and in sign.

**The catch that shapes this whole design:** that clean column came from assembly truth (human) and
from cohort recurrence — a length that is some other sample's consensus (tomato). Neither exists
when the pre-pass runs. So the answer cannot be "fit on homozygous loci". It has to be to sum over
the genotype, which is §3.

### 2.3 What this changes elsewhere in ng

**Two documents describe this step the way ng first proposed it, and neither has been brought into
line.** They are handled differently, and both are settled:

- [`ng_proposal.md`](ng_proposal.md) is **kept as a historical document** — the record of what ng set
  out to do — and carries a note at its head saying so and pointing here. It is not being revised.
- [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) will be revised as part of a
  separate pass over the interfaces, before any code is written. **The table below is the input to
  that discussion**, not a change request against the file.

A coder who builds against either today builds the thing §2.1 says is biased, which is why the
differences are enumerated rather than left to be noticed.

| where | what it says today | what replaces it |
|---|---|---|
| [`ng_proposal.md`](ng_proposal.md) §4 (lines 326-330, 354-357, 368-370) | a rough first-pass caller runs first, and only its confident calls teach the parameters | **No rough caller.** This step never chooses a genotype, so there is nothing to bootstrap from and no confident-subset gate (§3) |
| the same section, line 334 | lists GATK's DRAGstr calibration among the "rough passes" | `CalibrateDragstrModel` fits its parameters without ever calling a genotype (§3) |
| the same section, line 338 | calls HipSTR's stutter pre-pass "a rough length-genotyper" | It calls nothing. `em_stutter_genotyper.cpp` alternates between soft per-genotype posteriors and a refit, never taking a hard call ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3) |
| [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) lines 343-344, 360 | `fn summarize(&self, confident: &[ConfidentGenotype]) -> SampleSummary` | The input is the accumulated statistics of §5, not a list of calls. `CohortEstimator` (line 351) is unchanged |

**Decided: ng does not build a rough caller, here or anywhere.** Its only job was to break the
chicken-and-egg — you need genotypes to estimate parameters — and §3 breaks it instead by summing
over the genotype. **The decision survives however §4's open questions resolve**: every model still
in contention marginalises the genotype away rather than choosing one, whether it is fitted per
stratum or per locus, and by grid search or by EM. No branch of this design brings back a rough
caller.

**The durable constraint is narrower than the decision, and it is about accuracy rather than
permission.** A parameter *can* be estimated from thresholded calls — production estimates several
that way, and the arithmetic is perfectly well defined. It comes out **biased**, by the amount §2.2
measured: the slippage rate inflated 2.4-fold with its direction reversed. **The bias belongs to the
thresholding rather than to the step doing it**, so it would follow the practice anywhere it was
reintroduced, which is why the constraint outlives the particular decision above. A later step
wanting a called genotype for something other than teaching parameters (a display, a cheap
pre-filter) is a separate question this does not settle.

---

## 3. The estimator — sum over the genotype, do not choose one

At a site, an individual carries some number of copies of the alternative allele — none of them,
some of them, or all of them — and we do not know which. Rather than decide, score the site under
every possibility and add, each weighted by how common that kind of site is. Those weights — the
**genotype frequencies** — are unknown too, so they are fitted alongside the error rate.

**Written for any ploidy, because we want to support more than diploids.** For an individual with
`P` copies of the genome, write `j` for how many of them carry the alternative allele, `n` for the
site's depth and `k` for the reads supporting the alternative:

```text
                    P
     P(site)   =   Σ    genotype_frequency(j) · p_j^k · (1 − p_j)^(n−k)
                   j=0

where   p_j  =  (j/P)·(1 − ε/3)  +  (1 − j/P)·ε
```

`p_j` is *the chance one read shows something other than the reference base*: it draws one of the `P`
copies, hits an alternative copy with probability `j/P`, and either reading can be misread.

**The `3` is not decoration, and getting it wrong biases a headline output.** The two corruptions are
not the same event. A read over a *reference* copy shows a non-reference base if it is misread into
**any** of the three others — rate `ε`. A read over an *alternative* copy stops showing non-reference
only if it is misread into the **one** reference base — rate `ε/3` under a symmetric substitution
model, so it still counts as alternative with probability `1 − ε/3`. Writing `1 − ε` there charges the
alternative copy three times too much reversion, which biases `π_hom_alt` — the very output §5 argues
matters most for a landrace far from the reference. GATK's own code carries the same caveat as a
comment it never acts on (`DragstrParametersEstimator.java:229`). **Soft: the symmetric-substitution
assumption behind the `3`** — real substitution spectra are not flat, and if that matters the factor
becomes a fitted quantity rather than a constant.

**Two things are bundled in that formula and they come apart later.** The outer sum — score every
genotype, weight each by how common it is, add — is the **procedure**, and both paths use it. The
inner `p_j`, which says the only thing that can go wrong is a misread base at rate `ε`, is a
**noise model**, and it is the generic path's alone. A repeat tract can also slip a whole copy,
which `p_j` has nowhere to express. §3.2 separates the two properly.

**The diploid case is the one to keep in your head**, and it is where every other section's language
comes from — three terms, each obtained by putting `P = 2` into the `p_j` above:

| `j` | `p_j` | term | the name used elsewhere in this spec |
|---:|---|---|---|
| 0 | `ε` | `ε^k (1−ε)^(n−k)` | homozygous-reference rate |
| 1 | `½ + ε/3` | `(½+ε/3)^k (½−ε/3)^(n−k)` | observed heterozygosity |
| 2 | `1 − ε/3` | `(1−ε/3)^k (ε/3)^(n−k)` | homozygous-non-reference rate |

**A heterozygote is not exactly a half, and an earlier version of this table said it was.** It
printed `½` and `1−ε`, which is the very model the paragraph above denounces — `1−ε` charges the
alternative copy three times too much reversion. The middle row is `½ + ε/3` because errors are not
symmetric between the two copies: at a heterozygote, a misread of the reference copy *adds* a
non-reference read at rate `ε/2`, while a misread of the alternative copy *removes* one at only
`ε/6`, so the expected non-reference fraction sits above a half. The offset is small — `ε/3` is
3 in 10,000 at `ε = 0.001` — but it points the other way from the reference-bias term
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §8 is prepared to spend a parameter
on — which is 0.01 to 0.03 against a half, so thirty to ninety times larger. Small, then, and
recorded because it is free to get right rather than because it competes.

> **Names, once — because two symbols are overloaded in this field and both appear here.**
> This document uses **words** for these quantities and keeps the symbols for formulas.
>
> - The three weights above are **genotype frequencies**. In the formulas they are `π_hom_ref`,
>   `π_het`, `π_hom_alt`. **They are not nucleotide diversity**, which is what `π` usually denotes
>   in population genetics — that quantity appears here as *expected heterozygosity* ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)).
> - **Observed heterozygosity** (`Hobs`) is how often an individual's two copies of a site differ.
>   **Expected heterozygosity** (`Hexp`) is how often they would differ under random mating; it is
>   the population's **diversity**, written `θ` in that literature and `SFS_THETA` in our code
>   ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md)).
> - `θ` also means the **repeat-slippage rate** in the STR specs ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)). Same
>   letter, unrelated quantity. This document says *slippage* for that one and never writes `θ` for
>   it outside a quotation.
> - `ε` appears in **both** noise models and is **two parameters, not one**. Each is the error term of
>   its own model and they are fitted separately, deliberately
>   ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1). No formula holds both, so the symbol
>   is unambiguous inside any one of them; in prose this document names the path.
>
> And three words that are easy to confuse, because each of them names a kind of bucket:
>
> - a **cell** is one entry of an accumulated histogram — one `(depth, alt-count)` pair and how many
>   sites fell in it (§5). It belongs to the **data**;
> - a **stratum** is one group of loci that gets its own fitted parameter — on the STR path, one
>   *(motif period, repeat count)* combination ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
>   §4). **The generic path's rates are not stratified**: the error rate, the heterozygosity and the
>   distance from the reference are each fitted over all sites at once. (Its *windows* look like
>   strata and are not — a window is a unit of evidence the runs model reads in order, not a group
>   with its own fitted rate; [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4.) It
>   belongs to the **loci**;
> - a **step** is one value of a parameter being scanned during the fit — one error rate out of the
>   several hundred §3.1 tries. It belongs to the **search**. (This document also calls ng's stages
>   "step 3", "step 4"; those are numbered and always carry the number, so the two never collide in
>   practice.)
>
> How they fit together: one stratum's evidence is many cells, and the fit scores **every cell at
> every step**. So "a few hundred cells" and "a few hundred steps" are different few hundreds, and
> multiplying them is the cost of a fit. None of the three is a read or an observation.

**How the estimate comes out of this.** The formula above gives the probability of **one site's
reads**, for a given error rate and a given set of genotype frequencies. Multiply it across every
site and you have the probability of **everything this sample's reads showed**, still for that one
set of rates. Now try other sets. Whichever set makes the observed data most probable is the
maximum-likelihood fit, and **those rates are the estimate** — the error rate, the heterozygosity and
the homozygous-non-reference rate all fall out of the same search.

**"Now try other sets" is doing a lot of work in that paragraph, so here is what it means: lay the
parameters on a grid and score every combination.** A **candidate** is one complete set of values for
the parameters being fitted — an error rate of 0.001 paired with a heterozygosity of 0.002, say — and
*scoring* it means evaluating the formula above at every cell of the accumulated histogram and adding
the logarithms. Note what that is not: it is a pass over **cells**, not over reads. A dense
`(depth ≤ 100, alt-count)` histogram holds 5,151 of them; binned by depth under the rule in
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §4 it holds a few hundred.
**That binning applies to every histogram in step 4, the read-group one included** — it is the same
rule, and no accumulator is exempt. Take the plausible range of each parameter, step across it, and score every
combination; the best-scoring one is the fit. That is a **grid search**, and it is what this document
means wherever it says so.

**"Step across it" is where a continuous parameter becomes a finite list, and that needs saying.**
An error rate is a real number, so there are infinitely many candidates and no search enumerates
them. What a grid search does is **replace the continuum with a ladder of values** — a range and a
step — and accept that the answer lands on one of its steps. Two choices follow, and neither is arbitrary.

**Space the steps on a log scale, not a linear one.** These are rates spanning orders of magnitude,
and the distance that matters is a ratio: 0.0001 against 0.0002 is the same size of difference as
0.01 against 0.02. A linear grid would spend nearly all its steps at the high end and skip past the
region the data actually occupies. DRAGstr does this by laying its grid out in **Phred** units, which
are log units — its slippage grid is Phred 10 to 50 in steps of 1
([`DragstrHyperParameters.java:12`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/dragstr/DragstrHyperParameters.java),
declared `start:step:limit`), so 41 steps running from a rate of 0.1 down to 0.00001, each 1.26× the
one below it. *Two neighbouring declarations are easy to misread as part of the same grid and are
not:* line 15 is the variant prior, also scanned, and line 17 is the gap-open penalty, which
**DRAGstr does not scan at all** — it sets it analytically once the other two are chosen
(`DragstrParametersEstimator.java:221`).

**Make the step fine enough that the answer stops mattering, and no finer.** The question is not how
precisely the data *could* pin the parameter down — it is how precisely the **caller** needs it. These
are priors. A genotype likelihood is driven by read counts, and shifting `ε` by a few percent moves it
by far less than one more read would. **So a few percent between neighbouring steps is enough.**

*Soft, and the one number here worth checking before it is trusted:* "a few percent" is an argument
from what a prior does, not a measurement of what a caller tolerates. §10.1's synthetic fits can test
it directly — refit at two spacings and see whether anything downstream moves.

**At that resolution the scan is small on this path, and needs no refinement stage.** A quarter-Phred
step is about 6% apart in probability, so Phred 10 to 50 is 161 steps — which is DRAGstr's own
spacing for its gap-open parameter. One parameter is 161 scores, once per read group. **That is a
single flat scan, end to end, and nothing is skipped.**

**The STR path does not scan, and an earlier version of this paragraph said it did.** It priced
"the STR path's three, scanned together" at 4.2 million — three being the substitution rate, the
direction split and the distance decay. Three things are wrong with that and
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2 carries the replacement: the
substitution rate is fitted in closed form from two counters and is not an axis at all; the
parameter that *is* missing from the three is how often a read slips; and 4.2 million is per
(read group × **stratum**), of which there are several hundred, each combination costing a climb
rather than a pass. That path searches from several starting points instead, which the same
document's measurement of the profile curve's shape is what permits.

### 3.1 How the maximum is found, and why the shape of the likelihood decides it

**What decides the choice: whether a method finds the true maximum. Not how fast it runs.** All the
methods below hunt for the highest point of the same score. Any two of them that find it give the
same answer. They differ only in whether they can miss it. **Speed is not a criterion here** — if
fitting turns out to be slow, that is something to optimise when it happens.

**Say plainly what was and was not established.** That GATK's DRAGstr grid-searches is *precedent* —
evidence the approach survives contact with real STR data at scale — and precedent is not a
measurement. An earlier draft of this section presented the two as if they were the same thing.

**The candidates.** A grid with successive refinement scans coarsely, then rescans around the winner.
Quasi-Newton (L-BFGS) follows the gradient. Nelder–Mead needs no gradient. **Expectation-maximization**
is the textbook choice for a mixture like this one, and is what HipSTR uses
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3).

**What the literature says about EM specifically, because it is the one with a reputation for
reliability.** Dempster, Laird & Rubin (1977) proves **monotone ascent** — an iteration never
decreases the likelihood — and that is often misread as convergence to the maximum. It is not: Wu
(1983, *Annals of Statistics*) established that EM reaches a **stationary point**, which is the
global maximum only where the likelihood is unimodal. For general finite mixtures it frequently is
not, and multiple random restarts are standard practice for that reason (McLachlan & Peel, *Finite
Mixture Models*, 2000). Convergence is also only linear, at a rate set by how much information is
missing (Redner & Walker, 1984, *SIAM Review*) — **slowest where the mixture components overlap**,
which for us is low coverage, where a heterozygous site and a homozygous-reference one look alike.
**So monotone ascent buys less than its reputation suggests, and buys least in our regime.**

**But our likelihood is not a general mixture, and that is what settles this.** Hold the noise
parameters still — the error rate, or the STR path's three slippage numbers — and ask only which
genotype frequencies fit the data best. **That question has exactly one answer, with no false summit
anywhere on the way to it.** It is a property of the formula's shape, not a happy accident, so it
holds for every sample at every depth. Ploidy does not disturb it either: a tetraploid has five
frequencies where a diploid has three, and the surface still has a single summit.

*The precise statement, for anyone who wants to check it.* With the noise parameters fixed, the
log-likelihood in the genotype frequencies is

```text
log L (frequencies)  =  Σ over cells   count(cell) · log ( Σ over j   π_j · f_j(cell) )
```

and every `f_j(cell)` is now a **constant**. So the inner sum is linear in the frequencies, the
logarithm of a linear function is concave, and a non-negative sum of concave functions is concave.
The frequencies live on a simplex, which is a convex set. Maximising a concave function over a convex
set has no local maximum that is not also global. (*Concave* means the surface never bends upward —
so it has no dip in which a second peak could hide.) These are the standard composition rules — Boyd
& Vandenberghe, *Convex Optimization* (2004) §3.2 — and that a mixture likelihood is concave in the
mixing weights when the components are fixed is the basis of Lindsay (1983), *The geometry of mixture
likelihoods*, Annals of Statistics.

*One caveat, so the claim is not stronger than it is.* "No local maximum that is not global" is not
the same as "one maximum". The maximising frequencies are unique only if the component profiles are
linearly independent across cells, which fails in degenerate corners — the genotypes become
indistinguishable where `p_0 = p_2`, that is at `ε = 0.75`, and the summit is then a flat ridge.
(An earlier version said `ε = 0.5`, which reads off the retired `1−ε` form of `p_2`. Either way it
is far above the ladder, whose highest error rate is 0.1, so it can never be reached in a fit.) The climb still cannot get **stuck**, which
is what the decision below relies on.

**About the noise parameters we have no such result, and the difference is in what is proved rather
than in what is known to happen.** The argument above does not extend to them: profiling the
frequencies out can leave a curve in `ε` with more than one hump. **Nobody has shown that it does.**
What can be said is weaker and worth saying exactly:

- for the frequencies there is a **proof** of good behaviour;
- for the noise parameters there is **no proof either way**, and the analogy points the wrong way —
  `ε` is a *component* parameter, fixing where each component sits, and component parameters are
  where multimodality classically shows up in mixtures, which is why multiple restarts are standard
  practice there (McLachlan & Peel, *Finite Mixture Models*, 2000). The mixing weights are the
  well-behaved half in that literature too.

**The decision below does not depend on resolving this**, and that is the point of it: stepping
through one or three parameters end to end is correct whether the curve has one hump or three. The
missing proof is the *reason* for the scan. If someone establishes unimodality, the payoff is that a
one-dimensional optimiser would do instead. §9.3 records it as an open question, with plotting the
curve as what settles it — worth doing before the refinement spacing is chosen, since a second hump
would mean the coarse pass has to be dense enough to find it.

**Decision: step through the noise parameters; at each step, climb to the best genotype frequencies.**
In full:

1. **Scan the error rate across its whole plausible range**, coarsely at first, with the steps
   spaced on a log scale for the reason §3 gives. **The coarse pass is what looks everywhere**; how
   coarse it may be is bounded by §9.3, since a hump narrower than the spacing can be stepped over.
   *On the STR path this step is replaced by a search from several starting points*
   ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2), for the two reasons §3 above gives
   — three parameters times several hundred strata, and a measurement of the curve's shape that this
   path has and that one does not.
2. **At each step, hold those values fixed and find the genotype frequencies that score highest** —
   scoring, as always, over every cell of the accumulated histogram. This climb cannot get stuck,
   because with the noise parameters fixed the surface has no local summit that is not also the
   global one. **Any method that climbs will do**: EM is a reasonable default and is what HipSTR
   uses, but nothing here depends on it being EM.
3. **Take the highest-scoring step.** That value is the estimate, and the genotype frequencies found
   at it come with it. There is no refinement pass: §3's spacing is already finer than a caller can
   feel, so a second, narrower scan would move the answer by less than the answer matters.

**Why this beats either method alone.** EM by itself would never step along the error rate, so it
could stall on a ridge there. A grid over *everything* would spend its resolution on the frequencies,
which never needed searching at all. Splitting the two puts the effort exactly where the difficulty
is.

*In the standard vocabulary this is a **profile likelihood** over the noise parameters: the
frequencies are maximised out at each value of the parameter being scanned, leaving a curve in that
parameter alone.*

**One more thing follows, and it matters most where the fits are many.** Scanning leaves no
convergence flag to check per fit. The STR path runs one fit per (read group × stratum) — thousands
of them — and a flag nobody reads is how a badly-fitted parameter reaches a caller. Step 2's climb
does have one, on the half that provably cannot get stuck; **treat a climb that fails to converge as
a bug rather than a data condition**, since on a concave surface it has no legitimate reason to, and
assert convergence in the tests (§10.1) rather than propagating a flag no consumer would read.

*One argument for scanning is lost by dropping intervals, and it was the strongest of them.* A scan
leaves a score at every step, so an interval could be read off it directly; an optimiser returns a
point. §6 no longer wants that, so what remains for the scan is the reason above and the one that
matters most — **we have no proof the curve has a single hump** (§9.3), and a scan does not care.

**Cost is not a criterion here, but one arithmetic is worth having.** The STR path fits per
(read group × stratum), and those fits are near-independent — which is why the per-stratum count is
not alarming, and why the neighbour information that does exist is used deliberately rather than
implicitly (thin strata borrow, and the fitted sequence is held monotonic:
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3). If fitting ever does dominate, that is
the moment to measure it, not now.

**The approach has a name, and the name is worth having because two rivals sit next to it.**
Summing an unknown out of a likelihood — here the genotype, each possibility weighted by how common
it is — is called **marginalising** over it, and what remains is the **marginal likelihood** of the
reads given the rates alone. Maximising that over the rates is **maximum marginal likelihood**, and
it is what the rest of these documents mean by "the fit".

What makes the name earn its keep is that the genotype is a **nuisance variable** — not a parameter
anyone wants, not observed either, and in the way. There are three things to do with one, and two of
them are wrong:

| what to do with the unknown genotype | what it costs | who does it |
|---|---|---|
| **choose** one, and keep the confident choices | the estimate is drawn from the sites that were easy to call. §2.2 measures the damage: 2.4-fold, with the sign reversed | production (§2) |
| **maximise** over it, as one more free parameter per site | one new parameter for every site, so the parameter count grows with the data — and the bias **does not shrink as data accumulates** | nobody here; it is the tempting wrong answer |
| **marginalise** it away | nothing. The genotype leaves the expression entirely, and what remains is a proper likelihood in the parameters we want | this document, GATK's DRAGstr, and HipSTR |

**The middle row is worth understanding, because it is not obvious that it fails.** It is the
incidental-parameters problem (Neyman & Scott, 1948), and the sting is that more data does not cure
it, because each new site brings its own new parameter with it. Work through a 3-read site showing
one alternative read. Its best-fitting genotype is **heterozygous** — `(½)³` scores far above
`ε·(1−ε)²` at any sensible error rate — so choosing the best genotype books that site as a het,
including where that read really was an error. **The error rate is then squeezed toward zero and
heterozygosity inflates**, and sequencing ten more genomes only adds more sites doing the same thing.

**This is a low-coverage argument, and it is worth being exact about where it stops**, because an
earlier version said "so does every other site with at least one alternative read" and that is
false. A single alternative read favours het only while `(½)^n` beats `ε(1−ε)^(n−1)` — **up to about
nine reads** at `ε = 0.001`, and six at `ε = 0.01`. At depth 30 with one alternative read the
best-fitting genotype is homozygous-reference by six orders of magnitude. So the incidental-parameters
problem bites hard at tomato's three reads per plant and not at all at HG002's 300×; what makes
marginalising the right choice at every depth is that it needs no such argument. Marginalising books the same
site as a *fraction* of a heterozygote and a fraction of a homozygous-reference site, weighted by how
common each kind is, so an error and a real allele each get their share and neither is lost.

**Nothing is counted, and that is the whole point.** No site is classified, no threshold is applied,
nothing ambiguous is set aside, and the per-base error rate comes from the reads rather than from
the qualities the instrument reported about itself. A 4-read site with two alternative reads does not
have to be called anything; it contributes its share of evidence to every rate at once.

**What ploidy costs, and what it does not.** It does not touch **what Stage-1 accumulates**: the
site enters the likelihood only through its depth and alternative count, so the `(depth, alt-count)`
histogram of §5 is a sufficient statistic at *any* ploidy, and nothing about the expensive half of
this step needs to know `P`. Ploidy enters only at the fit, where it costs `P + 1` terms instead of
3 and `P` free genotype frequencies instead of 2. **Decision: make the likelihood and the
accumulator ploidy-generic from the first line of code**, since doing so costs a loop bound and
retro-fitting it later would not. Ploidy is per sample and per region — sex chromosomes and
aneuploidy make it vary within one genome — so it is an input to the fit, not a global constant.

**None of this is a novel construction, which is worth saying because it is the precedent §2.3
corrects `ng_proposal.md` about.** GATK's STR calibration (`DragstrParametersEstimator.java`,
vendored) does exactly the above: the same three-term mixture, marginalised the same way, maximised
by a grid search over 41 error values crossed with 41 values of a variant prior. **There is no rough
genotyper anywhere in DRAGEN's calibration** — the tool that most resembles this step already refuses
to call a genotype in order to estimate one.

**Two more of its choices are copied, and both live in
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3 because both need strata, which only the
STR path has**: pooling a stratum that holds too little data into its neighbours, and requiring the
fitted values to change monotonically from one stratum to the next, merging and refitting where they
do not. Its **noise model** is not copied at all — it collapses every kind of slippage into one
indel-error rate, which describes neither the direction a read slips nor how far — the two things §3
of that document measures and finds worth fitting separately.

**Expectation-maximization is not a fourth option, and treating it as one confuses two decisions.**
EM is an *algorithm* for maximising exactly this marginal likelihood: its expectation and
maximisation steps provably increase it at every iteration. You reach for it when the sum over the
nuisance is intractable, or when the maximisation has no closed form. Ours is `P + 1` terms summed
over a few hundred binned cells, so it can simply be evaluated and the parameter space scanned.
**Grid search versus EM is therefore an optimisation choice and cannot change the answer beyond
optimisation error** — both are climbing the same function. What *can* change the answer is the model
being marginalised, and that is §4.

### 3.2 One procedure, two noise models

**Everything above is the *procedure*: marginalise the genotype away, then maximise what is left. It
applies unchanged to both paths.** What differs is the **noise model** — what each path assumes can
go wrong with a read.

**The formula in §3 quietly carried one.** Its `p_j` says a read shows the wrong allele with
probability `ε`, a per-base substitution rate, and nothing else can go wrong. **That is the generic
path's noise model**, appropriate where the alternatives to a reference base are three other bases,
and [`parameter_prepass_generic.md`](parameter_prepass_generic.md) is what it produces. A repeat
tract can also **slip**, showing a whole copy more or fewer than the allele carries, which `p_j` has
nowhere to put — so the STR path keeps this procedure and replaces the noise model, in
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md), together with the measurements that decide
its shape.

**What the two paths share is everything around the likelihood**, and it is the part that carries the
bias this spec exists to remove:

- the genotype is summed over, never chosen;
- a sufficient statistic is accumulated per read group, and per stratum where the path has strata;
- where there are strata, thin ones borrow from their neighbours instead of being fitted on noise;
- the maximisation runs over a small parameter space, by whichever search §3.1 settles on;
- the output is frozen before genotyping.

So step 4 has **one implementation of the procedure and two noise models behind it**, chosen by
marker type — not one estimator, and not two estimators either.

---

## 4. What is not settled: two data objects, two models

Three separate choices sit behind this step, and they are easy to run together into a single "which
approach?" — which is wrong, because one of them is settled and the other two are not, and the
settled one holds whichever way the others go.

| the choice | the options | status |
|---|---|---|
| **how we estimate** | marginalise the genotype, or call it and count | **settled** — marginalise (§3), whatever the other two rows do. §9.1 measures *how much* it wins by; it does not reopen the choice |
| **what we accumulate** | genome-wide histograms, or census sites | **open** — build both, measure, decide (§4.1) |
| **what the genotype is weighted by while we fit** | one pooled set of genotype frequencies, or each locus's own | **open on both paths**, and only askable where a census is the data object (§4.2) |

**GangSTR is not one of these.** Its distinguishing move — four classes of read evidence, to reach
alleles longer than a read — is an STR matter with no generic counterpart, and it is worked through
in [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §7.

### 4.1 Histograms or census sites — build both and measure

**The two objects differ in what they can reach, not only in shape.**

- **The whole-genome histogram** pools every analysed site within one sample. Its advantage is
  precision: it uses every site there is and throws nothing away. Its limitation is absolute — **it
  cannot yield the population's diversity**, because that needs the same loci in every sample and a
  histogram has forgotten which loci it saw.
- **The census sites** keep raw evidence at a bounded set of positions, identical in every sample. An
  error rate estimated from two million sites is still an error rate, so this object can do
  everything the histogram can — **and** the cross-sample work the histogram cannot. It could cover
  the whole genome in principle, but only with a second pass over the reads or a per-sample store, so
  what is actually on the table is whatever fits a fixed memory budget.

**Nothing stops the walk recording every position**: it visits them all anyway. What stops us is the
bill — even at the compact encoding of
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5 that is roughly 400 MB
per sample on tomato and 20 GB across a fifty-sample cohort, which is not a sensible thing to carry
beside a variant caller. So we fix a budget and keep the sites it buys: about **1 MB per read group**. *Soft, and chosen
rather than derived* — it is the round number at which the census stops being worth arguing about
beside a caller that already holds gigabytes, and it buys the two million positions the paragraph
below shows are ample.
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5 works the other direction,
deriving the count needed from an assumed segregating rate, and marks *that* as the measurement to do
first; if the two disagree, its derivation wins and this number moves.
**The subsample is a budget, not a limitation.**

**The budget costs remarkably little, because these are distributions.** The precision of a
distribution estimated from a random sample of sites depends on **how many sites it holds, not on
what fraction of the genome they are**: ten thousand variable sites characterise a frequency spectrum
whether the genome is 800 Mb or 3 Gb. Two million selected positions is about a four-hundredth of the
storage and nearly all of the information.

**If a bounded sample of loci gives good enough estimates, it is the better object, because it is
strictly more capable.** That is worth testing rather than assuming, so: **implement both, estimate
the per-sample parameters both ways on synthetic data where truth is known, measure the memory each
costs, and decide afterwards what to keep.** §10 carries the measurements.

**What "good enough" should mean.** The histogram wins on precision by construction — four hundred
times the sites is one twentieth of the standard error — so the question is never *which is more
precise*, and a comparison framed that way answers itself. It is worth agreeing roughly in advance
what would count as enough, not to bind the decision — the numbers may show a structure nobody
predicted, and the final call should weigh them — but because "precise enough" has no meaning until
someone says what the precision is *for*. Three criteria, in increasing usefulness:

- **Is it biased?** Look at the mean error over repeated simulations, not one run's error.
  Imprecision shrinks as the site budget grows; bias does not, and bias is exactly what a subtly
  broken selection rule would produce. **This is the criterion that must not be waived**, because it
  is the one that says the object is wrong rather than merely small.
- **Is the difference wide enough to feel?** Take the gap between the two routes' estimates, perturb
  each parameter by that much, and see whether anything downstream moves. This replaces an earlier
  version that perturbed by each parameter's fitted interval — §6 no longer emits one, and the gap
  between the two routes is the more direct quantity anyway, since it is the thing being decided.
- **Do the calls change?** Fit both ways, genotype the same synthetic reads with each parameter set,
  and compare the calls. If none changes, the difference does not matter *by definition*. This is the
  criterion that is not arbitrary, because it measures the thing the parameters exist for rather than
  the parameters themselves.

**The expected answer, written down so that a surprise is recognisable.** Two million sites at three
reads is six million observations; at an error rate near 0.001 that is roughly 6,000 errors, pinning
`ε` to about one part in eighty. Heterozygosity at one per kilobase leaves ~2,000 heterozygous sites
among the census sites, so about 2% relative — **and that one is optimistic**, being a plain count
error. At 3 reads a het site shows no alternative read about an eighth of the time and is invisible,
and the single-alternative cell is mostly error background, which pushes the real figure nearer 3%
even before the two rates are fitted together. **Both sit far inside anything a genotype
likelihood can feel.** So the precision half of this comparison ought to come out comfortably for the
census sites, and if it does not, suspect the selection rule before suspecting the arithmetic.

**Only the first criterion can be run today, and that decides when the question can be closed rather
than merely postponed.** Perturbing a parameter to see whether anything moves, and asking whether
any call changes, both need the downstream caller, which is not built. The bias check needs only
synthetic reads and the two fits. So:

- **If the two estimates agree to within rounding**, at a memory cost that is not absurd, **the
  question is settled here.** A difference that small cannot propagate to any caller, so there is
  nothing a downstream measurement could later overturn — waiting would buy no information.
- **Otherwise the decision stays open**, and deliberately. An appreciable difference is exactly the
  case where only the downstream can say whether it matters, so choosing now would mean inventing an
  answer. Keep both implementations, record the numbers, and revisit when there is a caller to ask.

**What to record now, so that whoever finishes the comparison need not re-run the walk:** the fitted
parameters from both routes, the evidence behind each, and the memory each cost, written up in
`devel/ng/reports/`. Criteria two and three are then a day's work against those numbers rather than a
repeat of the experiment.

**One thing the census sites cannot replace, whatever the comparison says: the windowed statistic
inbreeding needs** ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4). It estimates
a rate *per window*, not per genome. A 100 kb window holds about a hundred thousand sites but only
~250 selected ones, so at one heterozygote per kilobase it would carry about **0.25** expected
heterozygotes instead of a hundred — far too thin to tell a run of homozygosity from ordinary
sequence. A subsample is the wrong instrument for a local quantity, and no site budget within reach
changes that.

**STR loci are not such a case: they get a census of their own.** Built by the same rule over
the delimitation region typing already produces
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1.1), and holding tract
lengths rather than base counts. It exists chiefly to supply **STR diversity** — a different quantity
from the generic one, because repeat tracts mutate orders of magnitude faster, and one the STR path
currently takes from a SNP-scale constant
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §6). Whether it also improves the slippage
priors is an open question the set makes askable (§4.2).

**One honest deflation of this comparison, so nobody expects too much of it.** Deleting the
histogram route saves less than it sounds. The windowed histogram is built regardless — the runs
estimator needs it and no census can replace it — and the whole-sample histogram is a free fold of
it. So what a census victory actually removes is the **read-group** histogram, which
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §9 prices at kilobytes. The comparison
is therefore about which estimate to *trust*, not about which object to stop building; only one of
the two was ever chargeable.

**And the two are endpoints of one axis rather than rivals.** Raise the site budget far enough to
make the local statistics work and the census grows into the histogram — keeping every
site while forgetting its identity is precisely what a histogram *is*. The choice is where to sit on
that axis: identity for some sites, or all sites without identity.

**What it cannot do is anything that needs a site's neighbours.** Two million positions across 800 Mb
is one every four hundred bases, and they are spread out on purpose so that no two are close enough
to be inherited together. Linkage between variants, haplotypes, and anything else that reads a
stretch of genome rather than a pile of separate sites are therefore out of reach — by design.

**Which loci — and the requirement is stronger than "at random".** For a given genome and a given set
of analysed regions, **every sample must select the identical positions**, even though the samples are
walked independently — on different machines, at different times, with no sample able to see what any
other chose. The set therefore cannot be negotiated or handed round; it has to be **computed from the
run's inputs alone** and arrived at separately by each sample. Those inputs are the reference,
intersected with the `--regions` BED when one is given.

**The BED belongs to that identity; it is not a detail of how the run happened to be invoked.** Two
things go wrong when it is treated as incidental. Selecting across the whole genome in a run
restricted to a BED leaves almost every chosen position unvisited and empty. And two samples analysed
over different region sets share no loci at all — not a degraded estimate but a meaningless one.
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) gives the rule that delivers
this, and what has to be recorded so a later gather can check it held.

### 4.2 What the genotype is weighted by while a prior is fitted

**First, what is not in question. This step produces priors, and priors are pooled by definition.**
Nothing here emits a per-locus number, and nothing here decides what the caller does at a locus —
per-locus allele frequencies, per-locus genotypes, and everything else that reads one locus at a time
belong to the calling phase and are not this step's business.

**The live question is narrower, and it is about the fitting rather than the output.** Marginalising
the genotype away means weighting each possible genotype by how common it is (§3), so some weighting
has to be supplied *while the prior is being fitted*. This step supplies one pooled set. HipSTR
instead lets each locus supply its own, fitting that locus's allele frequencies in the same loop
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.3). **Both produce a pooled parameter at
the end**; they differ in what the sum inside was weighted by.

**Which data object is in hand decides whether the question can be asked at all.** A histogram has
forgotten which site each observation came from, so there is no locus to take a weighting from and
the pooled set is the only thing available — by construction of the accumulator, on **both** paths.
A census has not forgotten: it holds several samples at the same locus, which is the whole reason it
exists. **There are two censuses, one per path (§4.1), so the question is askable on both** — it
rides on the census, not on the marker type.

**What differs between the paths is where the per-locus version already lives, and what it is
expected to buy.**

- **Generic.** Letting every site carry its own allele frequency is not a new model to build — it is
  the **frequency spectrum**, which [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4
  computes from the same census. So the question here is not "build a per-site model?" but "should
  the pre-pass fit be weighted by the spectrum rather than by one heterozygosity?", and the object it
  needs is already specified.
- **STR.** The per-locus version is HipSTR's, and it is
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §8.5 — *can a per-locus stutter model beat
  the per-stratum one?*

**Soft: the payoff is expected to be larger on the STR path, and that is a guess.** A repeat tract is
multi-allelic and its allele spectrum varies far more from locus to locus than a substitution site's
does, which is the motivation HipSTR gives. **Nobody has measured how much either is worth**, so the
asymmetry is a reason to look at the STR case first, not a reason to drop the generic one.

**Both carry the same data-volume caveat:** a census holds far fewer loci than a histogram holds
sites, so a per-locus fit is not being handed the same evidence as a pooled one, and any comparison
has to say so.

**One thing this replaced, because the wrong version is memorable.** An earlier draft asked "marginal
likelihood or expectation-maximization". That has the answer *identically*: EM is an algorithm for
maximising a marginal likelihood, not an alternative to one (§3).

---

## 5. What Stage-1 accumulates

**Chemistry belongs to the library preparation, not to the individual.** A sample sequenced from
two libraries has two error rates and two stutter behaviours; averaging them describes neither. How
often the two come apart is measured and recorded next door
([`read_groups.md`](read_groups.md) §1, from a 2,085-file, 68-project tomato archive survey): in 157
of 1,707 samples — 133 carry two libraries, 20 carry three, and four carry 7, 16, 16 and 42.

**A trap that comes with the read-group grain: read length is chemistry too.** Suppose one library
was sequenced with 100 bp reads and another with 250 bp. The shorter reads span fewer of the long
repeat tracts, so that library's slippage rate is measured over a different mix of loci — and comes
out lower, without either library having slipped any more or less than the other. Stratifying by
`(period, repeat count)` removes the confound: inside one stratum every locus has the same tract
length, so the two libraries are compared like with like. That stratification exists for a stronger
reason anyway — slippage depends far more on repeat count than on anything else
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4) — but this is a second reason to keep it.
What is left is not quite bookkeeping, and there are two residues rather than one:

- **Empty strata.** The short-read library's high-repeat strata hold nothing, so they borrow (§6) and
  must be marked as having borrowed.
- **A biased subset of reads inside the strata that are not empty**, which stratification does *not*
  remove. Within one stratum the two libraries cover the same loci, but a 100 bp read spans a 40 bp
  tract only when the tract sits near its centre, so that library's spanning reads carry less flank
  and are anchored worse; the 250 bp library spans the same loci freely. **Like-with-like in loci is
  not like-with-like in reads.** *Soft and unquantified* — nothing here measures how much slippage
  depends on flank length. **Settled by:** §10.3's third case, which manufactures two read groups of
  known chemistry, run with two read lengths instead of two error rates.

**So: never compare two read groups on an overall slippage rate** — that number mixes chemistry with
read length. Compare them stratum by stratum.

**Soft: neither cohort we measure on exercises this.** The tomato benchmark declares one read group
per sample and *synthesizes* the library name (`tmp/tomato_stutter/rg_check.tsv`); the HG002 300×
BAM declares one `@RG` for the whole file (`ID:HG002 SM:HG002 LB:HG002`). On both, read group,
library and sample coincide. The grain comes from the archive survey and from first principles, and
is unvalidated on the data this spec's numbers come from. §10.3 says how to test it anyway — by
manufacturing the contrast a simulator can supply and our cohorts cannot.

ng already carries the evidence at this grain: `ObservedSequence` includes its read group in its
identity ([`src/ng/locus_generation/mod.rs:178`](../../../../src/ng/locus_generation/mod.rs)), so an
allele seen from two groups is two rows with their own quality moments, and the read-group table
([`src/ng/read/input/read_groups.rs:43`](../../../../src/ng/read/input/read_groups.rs), the
`ReadGroup` record) resolves that identity to sample, library and experiment — recording whether each grouping name was
**declared** by the file or **synthesized** because the file gave none.

**Decision: estimate the *chemistry* parameters per read group; expose the fold.** Read group is the
finest grain available and the safest default for anything that describes how the DNA was prepared
and read — the error rate and the stutter behaviour. (The sample's own quantities are fitted at
sample grain from their own accumulator.) Which grain a run fits at is a knob, because a library
sequenced across four lanes usually shares one chemistry and pooling its lanes buys precision for
free, while two libraries of one sample usually do not. The estimator does not guess: it fits at
whatever grain it is handed.

### 5.1 The five objects, and who owns each

| accumulator | key | accumulated | specified in |
|---|---|---|---|
| generic noise | **read group** | a histogram of `(depth, alt-count)` cells, **including `k = 0` cells**, which is what production discards. Base qualities are neither carried nor used to filter reads: the model has one error rate ([generic](parameter_prepass_generic.md) §2) | [generic](parameter_prepass_generic.md) §2, §3 |
| windowed heterozygosity | **sample × genomic window** (100 kb) | the same histogram shape, but each site entered once with its **total** depth and alternative count | [generic](parameter_prepass_generic.md) §4 |
| STR noise | read group × `(period, repeat count)` | a histogram **of loci, not of reads**: each locus is reduced to how many of *its own* reads fell at each whole-repeat offset from the reference tract length, and the table counts how many loci had each such shape. Plus a bucket for reads differing by something that is not a whole number of copies, and **bases compared and bases mismatched** pooled over the stratum — the length buckets alone cannot yield an error rate | [ssr](parameter_prepass_ssr.md) §4.1 |
| generic census | **read group × selected position** | reads supporting each allele (A/C/G/T + other) at each position of a fixed set drawn from the analysed regions, **identical in every sample** | [census sites](parameter_prepass_census_sites.md) §2 |
| STR census | **read group × selected locus** | reads at each whole-repeat offset from the locus's reference tract length, a non-whole-repeat bucket, and the same two composition counts, over a fixed set of STR loci | [census sites](parameter_prepass_census_sites.md) §2.1 |

**The first two rows look like the same histogram twice, and the difference between them is why both
exist.** One is keyed by read group, so a site covered by two libraries enters it **twice, at two
depths**; the other is keyed by sample, so that site enters **once, at its full depth**. Which one is
right depends on what is being counted, and step 4 counts both kinds:

- a **rate per read** — how often a read shows an allele the individual does not carry — loses
  nothing when a site's reads are split between two tables, and *needs* the read-group key, because
  chemistry differs between libraries;
- a **rate per site** — how often a genome is heterozygous, how often it differs from the reference
  at both copies — cannot be counted at all in a table where one site has become two.

**Neither table can do the other's job, and no fold repairs that.** What *is* free is the fold from
the windowed histogram to a whole-sample one: summing a site's windows is addition, so the sample's
`(depth, alt-count)` histogram costs nothing and is exact. **There are two accumulated generic
objects, not three.**

**Why the censuses take the read group and not the sample.** Two consumers want two different
grains, and unlike the histograms above they can both be served from one object. The cross-sample
work — diversity, the frequency spectrum, relatedness — wants **the sample's** counts at a position;
the error rate wants **the read group's**. Read group sits *below* the site, so summing a position's
read groups to get the sample's is addition of raw counts at one place, exact and free. **Keying at
the finer grain therefore costs a memory multiplier and loses nothing**, while keying at the sample
grain would throw away the chemistry axis irrecoverably — and it is the chemistry axis that lets the
censuses stand as an alternative to the histograms at all (§4.1). The multiplier is the read groups
per sample: **1 for 1,550 of the 1,707 samples in the survey above**, 2 or 3 for nearly all the rest.

**Decision: accumulate all five directly, in one walk.** The rule is narrower than "never fold", and
worth stating exactly, because two of the folds this document relies on are perfectly sound: **no
fold may cross a key that has already split a site.** Below the site, summing is addition of raw
counts at one place and is exact — a position's read groups (§5.1), a sample's lanes when a run
chooses to pool them (§5), a window's cells summed to the whole sample
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §1). Across sites, a histogram has
forgotten which site was which, and no fold puts that back. The walk over the reads is the cost; a second and third increment into a
second and third small histogram is not, and every object is then exact for its own consumer.

*Rejected: one generic histogram keyed by `(read group × window)`, with the other two derived by
summing.* It is the tempting shape — accumulate at the finest grain, let each consumer fold to what
it wants — and it does not work here, for a reason worth recording because the same temptation will
recur:

- **It is not exact**, for the reason just given: read group and window are not nested, and the thing
  that must stay whole is neither — it is the **site**. Summing a joint key's cells back over read
  groups counts each site once per group.
- **And it does not save memory where it matters.** Per sample a joint key costs
  `read groups × windows × cells`, against `(read groups + windows) × cells` for the two separate
  objects — so the separate pair is smaller once a sample has **two or more read groups**, by roughly
  the read-group count. At exactly one read group the joint key is smaller, by one histogram's worth
  of cells out of some 8,000 windows: about a hundredth of a percent, on 1,550 of the 1,707 samples
  in the survey above. **So memory is not the argument** — the argument is the previous bullet, and
  it holds at every read-group count.
- **The exactness is free precisely because the site is the unit.** Summing a site's read groups
  *before* binning it into the windowed histogram is not an approximation at all: these are raw
  counts at one position, and 20 + 10 really is 30.

**"Accumulate finely and let the consumer fold" is still the right principle — it just stops at the
site.** Below a site, read groups sum exactly. Across sites, a histogram has already forgotten which
site was which, so no fold can put back what the key threw away.

**Three of the five are reduced at the end of their own sample's walk; the two censuses are not,
and the difference is not a preference.** A histogram can be reduced to parameters as soon as that sample is
done, because everything it can ever say concerns that sample alone. The census sites cannot: their
whole content is *which samples carried an allele together*, and a per-sample reduction destroys
exactly that. So the census stays raw all the way to the gather and the histograms do not survive
their sample (§1.3).

---

## 6. Cross-cutting concerns

**Memory — the generic read-group histogram is small; the STR table is no longer certainly so.** The
per-read-group generic histogram is a few hundred cells, kilobytes per read group, and does not
register against the gigabytes the observation stream itself costs.

**The STR table's key changed, its size did not, and both halves of that are measured.** An earlier
version of this paragraph priced it at "a handful of offset buckets and two counters, likewise
kilobytes", which was true of a table of *reads*; the key is now one entry per **locus**, for the
identification reason [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 gives. That was
expected to cost memory, since the number of distinct locus shapes grows with a locus's depth. **It
does not: on HG002 at 300× the uncapped table is 0.43 entries a locus — 12,727 entries over 29,811
loci, 0.36 MB** ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.8).
Deep data deduplicates, because most loci at a clean tract are "every read at the reference length"
and what separates two entries is mostly their depth. **This object is not where step 4's memory
goes**; the windowed histogram still is.

| object | per sample | lives until | what a 50-sample run holds |
|---|---:|---|---|
| generic read-group histogram | kilobytes | end of that sample's walk | kilobytes × samples in flight |
| STR locus table | 0.36 MB over 29,811 HG002 loci at 300×, uncapped — **measured**, not arithmetic | end of that sample's walk | likewise × samples in flight |
| windowed heterozygosity | 30 MB tomato / 115 MB human, binned | end of that sample's walk | 30 MB × samples in flight |
| generic census | ~1–1.3 MB per read group | **the gather** | ~50–65 MB at one read group per sample |
| STR census | ~1–2 MB per read group | **the gather** | ~50–100 MB, likewise |
| the fitted parameters | bytes | the gather | negligible |

**The "lives until" column decides peak memory, not the "per sample" column** (§1.3). A histogram is
reduced to its parameters when its sample finishes, so a fifty-sample run never holds fifty of them —
it holds as many as it is walking concurrently. **The census sites are the only object that
accumulates across the cohort**, at ~50 MB for fifty samples. They are no longer the only one worth
watching, though: at 30 MB per sample the windowed histogram reaches that at two samples in
flight and passes it beyond, so **the concurrency multiplier is now the number that decides this step's peak
memory**, and §10.6 measures it rather than assuming it. Sizes are priced in
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5 for the censuses, and in
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §9 for the windowed histogram.

**Every one of those figures is arithmetic rather than measurement**, and §10.6 replaces them all —
including the "samples in flight" multiplier, which is a scheduling decision nobody has made yet.

**The windowed histogram is not chargeable to inbreeding.** The sample's own rates are fitted from
that same object, so it is accumulated on every run. What `F` adds to it is the **window** key:
without inbreeding — or with `F` supplied outright
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.4) — it collapses to one histogram
per sample, 3.7 kB rather than 30 MB. That is the line to hold when someone asks what inbreeding
costs: a factor of 8,000 on this object, and nothing at all on the other four.

**Errors, and how well determined each parameter is.** Two different things, and every parameter
carries both.

*Where it came from.* Too little data to fit is not an error, it is a provenance. **There are two
different borrows and they happen at different levels**, which is worth keeping straight:

- **within a sample** — a thin stratum takes its neighbouring strata's value, adjacent repeat
  counts at the same period;
- **at the cohort gather** — a read group too thin to fit *at all* takes the panel-pooled value,
  which is the only borrow that needs other samples.

Either way the parameter is **marked as having borrowed**, because one that came from a neighbour is
softer than one fitted in place and the consumer should be able to tell. Four states — **fitted
here**, **borrowed**, **defaulted**, and **supplied** by the person who ran the tool. The last is
neither a fit nor a fallback: it is more authoritative than a default and less checkable than either,
and `F` is the parameter that can arrive that way
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.4).

*How much data stood behind it.* Emit the count of observations the fit used — reads for a noise
parameter, sites for a rate. One number, already to hand, and enough to tell a stutter rate measured
on 80,000 reads from one measured on 8.

**Decision: no uncertainty intervals. These are priors, and a prior does not need one.** Computing one
honestly would cost an extra scan per parameter, because §3.1's climb keeps a best value at each step
rather than a curve — and nothing downstream would read it, since a caller consumes a prior as a
number. *Rejected: emit a likelihood-based interval from the score curve.* It is free for the scanned
parameters and manufactured for the rest, and that asymmetry is what makes it not worth having: half
the outputs would carry a real interval and half a differently-derived one, inviting a comparison
between them that means nothing.

**An interval would not have caught the failure this spec exists to fix, which is the honest reason
not to regret dropping it.** §2.2's slippage rate was inflated 2.4-fold with its direction reversed,
from a whole-genome measurement — so its interval would have been **narrow** and its provenance would
have read "fitted here". A bias is a confident wrong answer. What catches that is external:
recovering known parameters from synthetic data (§10.1), and agreeing with assembly truth on HG002,
which is §9.1's measurement.

**Concurrency.** Every accumulator is keyed — by read group, by stratum, by window — and keyed
accumulation is associative, so a region-sharded walk merges by summing and needs no communication
between shards. **Almost every fit is then per sample and so embarrassingly parallel across
samples** — the error rate, the stutter parameters, heterozygosity, the homozygous-non-reference
rate and `F` all come out of one sample's own accumulators. Only the panel-level work needs every
sample: the cohort's diversity, the sample-group clustering, and the fallback for a read group too
thin to fit at all. That gather is one small single-threaded pass over the per-sample summaries, not
over reads — and it runs inside cohort variant calling, not at the end of this step (§1.3), so it is
not on this step's critical path at all. The only other sequential step is the runs-of-homozygosity
pass, which walks 8,000 windows of one tomato sample in order and is irrelevant to the walk's cost.

**Determinism.** Every fit is a sum over cells; fix the summation order so no parameter varies with
thread count.

---

## 7. Reuse over rewrite

| what | existing code | how it is reused |
|---|---|---|
| read-group identity and grouping | `src/ng/read/input/read_groups.rs` | used as-is; the fold to library/experiment is already there |
| read group on each observation | `src/ng/locus_generation/mod.rs:166` | the read-group half of the two *chemistry* accumulators' keys (§5.1). The other two are keyed by sample |
| the three-genotype likelihood | `src/sample_summary/het.rs:266` (`observe_site`) | the three binomials are already computed there; they are **added** instead of compared |
| the site's sufficient statistic | `src/sample_summary/het.rs:146` (`SiteCounts::from_record`) | the `(depth, alt-count)` reduction is reusable — **minus** its `None` return for a pure-REF column (§2.1) |
| grid search over accumulated counts | GATK `DragstrParametersEstimator.java` (vendored) | algorithm copied, not code. Its stratification choices are the STR path's ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3) |
| the cohort gather | `../arch/ng_step_interfaces.md:351` (`CohortEstimator`) | **reused unchanged** — its signature is agnostic about how a `SampleSummary` was built; what it computes is [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) |
| the per-sample summary | `../arch/ng_step_interfaces.md:343` (`SampleSummarizer`) | **not reused: re-specified.** It takes the accumulated statistics of §5, not `&[ConfidentGenotype]` (§2.3) |

**No parity oracle.** This is not a port: production's estimator is the thing being replaced, so
agreeing with it would be failure. §10 says how correctness is shown instead.

---

## 8. Deferred, with a recommended home

- **Everything cross-sample** — the cohort's diversity, the frequency spectrum, contamination,
  relatedness, and the grouping of read groups into shared chemistries. This document's walk
  **accumulates** what all of them need and computes none of them. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md), which is a spec rather than a
  promise: none of it needs a second pass over the reads, which is the whole reason the census sites
  exist.
- **Separating contamination from sequencing error.** The `ε` fitted here is the sum of both, plus
  mismapping and everything else that makes a read show an allele the individual does not carry
  ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2). Splitting them
  needs population allele frequencies, since contaminant reads carry segregating alleles and
  concentrate at polymorphic sites while sequencing errors do not. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5. Nothing in the walk changes when
  it lands.
- **Alleles longer than a read, and the sample profiling that would be needed.** The two are one
  item: the profile *is* the instrument. Deferred on coverage rather than merit, and worked through
  in [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §7. **Home:** whichever spec takes on
  long-allele recovery; its first job is the tract-length error that document says is unmeasured.
- **The population genetics above diploidy** — several identity-by-descent coefficients instead of
  one `F`, and heterozygosity replaced by gene diversity. The likelihood and the accumulators are
  ploidy-generic (§3); only the population genetics on top is not.
  [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §7 states what breaks. **Home:** a
  spec of its own, which should pick definitions that degrade to the diploid ones at `P = 2` rather
  than keeping two code paths.
- **Overdispersion — a beta-binomial in place of the binomial.** One extra parameter for the
  correlated errors the independence assumption ignores.
  [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §8 states the problem and adopts the
  *other* of the two binomial assumptions; it is not itself the home. **Home:** a spec of its own,
  taken together with the reference-bias term, since both change the same `½`.

---

## 9. Open questions that belong to no single path

1. **Does marginalising beat production's threshold-then-count on our data?** — OPEN. The research
   note proposes the experiment: fit both on the same cohort and compare the recovered `ε` against a
   held-out truth set. *Leaning:* yes, on the strength of §2.2 — but the measurement that settles
   it is HG002, where truth genotypes exist, comparing the fitted stutter rate against the rate
   measured on known-homozygous loci. **This does not gate the design** — §4 settles the choice on
   §2.2's measurement and on the incidental-parameters argument in §3. What §9.1 buys is the size of
   the win on our own data, which is worth knowing and worth publishing.
2. **Histograms or census sites?** — OPEN by construction; §4.1 is the plan and §10.2 the experiment.
   §10.2 needs only sites drawn at known frequencies, so it is buildable now — which matters, because
   §4.1 says the bias criterion is the one that cannot be waived and it is also the only one runnable
   before a caller exists.
3. **Does the score curve over the error rate have one hump or more?** — OPEN, and it is the only
   part of the search still in doubt. §3.1 *proves* good behaviour for the genotype frequencies and
   proves nothing either way for the error rate, so a second hump could only appear there (or along
   the STR path's three slippage parameters). *Leaning:* one hump, for a model this constrained —
   but that is a guess, and the analogy to general mixtures points the other way. **Settled by:**
   emitting the curve from §10.1's synthetic fits and looking at it. No new machinery; it decides how
   closely spaced the steps have to be. **Note this question is not blocking**: the design steps
   through the parameter end to end precisely because the answer is unknown, so it is correct either
   way. What the answer buys is the option of a cheaper search, not a correction.

   **Partially answered on the STR path, and it is worth being exact about how little that
   transfers.** Profiling the *slippage level* — the other two slippage parameters and the genotype
   frequencies maximised out at each value — gives exactly one interior maximum on 41 rungs, in both
   worlds tried ([research note](../research/parameter_estimator_experiments_2026-08-06.md) §6.5).
   That is a different model, a different parameter and two worlds, so it is encouragement rather
   than an answer here; what it did settle is that path's own search
   ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2), where a flat scan was not
   affordable and the question was therefore blocking.

**Open questions belonging to the other four documents** — everything about stutter, including the
STR half of §4.2's weighting question
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §8); the `ε`/genotype-frequency coupling,
whether the two `F` estimators agree, and reference bias
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §11); whether the
censuses could replace these histograms
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §6); and the cohort gather
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §10).

---

## 10. How we know it works

*Tests that span both paths. Each path document adds its own.*

**What may be simulated, and where the line is.** **Reads are in scope.** Take a reference stretch,
put known alleles on it, emit reads at a chosen error rate, read length and read group, and the
estimator can be checked end to end against truth you set. Two read groups with deliberately
different chemistry are two different error rates in the generator. Contamination is reads from a
second individual mixed in. A collapsed paralog is reads from two similar loci placed at one. Errors
in runs rather than singly are a loop change. **All of that is a generator, not a project**, and the
repo already carries two of this shape (`src/ssr/cohort/sim.rs`,
`examples/ssr_delimiter_comparison.rs`, both seeded so a failure replays).

**Out of scope: simulating a population.** Coalescent genealogies, recombination, demographic
histories, pedigrees — the `msprime`/SLiM stack. That is a serious undertaking and this step does not
need it.

**And the reason it does not need it is a property of what this step accumulates, not a
compromise.** Every object here is either a marginal — a histogram that has forgotten which site each
observation came from — or a census of positions **deliberately chosen far enough apart that no two
are inherited together** (§4.1). Nothing this step estimates can see linkage, so a generative model
that produced correct linkage would be producing something no estimator here could read. Sites drawn
independently at chosen allele frequencies are therefore not an approximation *for this purpose*:
they are the right generator. Even the cohort's frequency spectrum can be tested that way — draw the
frequencies from a chosen spectrum, draw genotypes, and check what comes back.

**Where that stops being true is any claim about linkage or ancestry itself** — a relatedness estimate
tested against a real pedigree, or anything reading a stretch of genome rather than a pile of separate
sites. Those are out of reach here by the same design decision, and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) defers them on those grounds.

The cheapest test of the lot needs no simulator at all — §10.3's first case relabels half of one read
group's reads in the tomato CRAMs, which is a metamorphic check requiring no truth.

1. **The estimator recovers known parameters.** Draw genotypes from known frequencies,
   draw a depth per site, draw the alternative count from the model at a known error rate, and fill
   the accumulator directly — no reads, no alignments, no reference. Then fit. The fit must return the
   values drawn from. **This is the primary test, and it is a test module rather than a project.** **Run it at more than one ploidy** — at minimum
   `P = 2` and `P = 4` — because §3's likelihood is written for any ploidy and an untested loop
   bound is an assumption, not a generalisation. A tetraploid case also exercises the intermediate
   dosages (one and three alternative copies of four) that have no diploid counterpart.
2. **The two data objects are compared head to head, on the same synthetic data.** Fit the per-sample
   parameters from the whole-genome histograms and from the census sites, against the same known
   truth (§4.1). Report each fit's error against truth, and the spread of that error across repeats,
   at more than one coverage, since the census sites lose relatively more where reads are scarce. **Repeat each
   simulation** so that bias is separable from imprecision — §4.1 says why that is the criterion that
   cannot be waived, and why it is also the only one of the three runnable before a caller exists.
   The comparison is meaningful only for the *global* rates; the windowed statistic and stutter are
   excluded for the reasons §4.1 gives. **This test closes the question only if the two agree to
   within rounding**; otherwise its job is to record the numbers well enough that the remaining
   criteria can be applied later without repeating it.
3. **The read-group grain does what §5 claims it does.** No cohort we hold exercises this — tomato declares one read group per sample, HG002
   one for the whole file (§5) — so the strongest cases need a simulator to **make two read groups
   whose chemistry differs by a known amount**. Three cases, in increasing strength:

   - **Split a real read group in two** (a metamorphic relation — no truth needed, and it runs on
     the tomato CRAMs today). Relabel half of one read group's reads as a second group. The two
     per-group fits must agree with each other, and with the fit on the unsplit group, within
     error. This exercises the whole per-group path — keying, thin-stratum borrowing, the fold —
     on real reads, and a failure names its own bug.
   - **Two synthetic read groups, one chemistry.** Generate both from a single `(ε, stutter)`
     parameter set. Fitting them separately and pooling must equal fitting the pool. This is the
     arithmetic of the fold in §5's "expose the fold" decision.
   - **Two synthetic read groups, deliberately different chemistry.** Generate one at `ε = 0.001`
     and one at `ε = 0.01`, with different stutter levels. The per-group fit must recover both
     values; forcing a single pooled fit must be visibly worse than the two-group fit, by a
     margin the test asserts. **This is the case that tests §5's claim rather than its
     arithmetic** — that the read group is where chemistry lives, and that merging groups destroys
     something real. Nothing else in this spec can check that.

   **This third case carries one more assertion, and it is the one that would have caught the design
   [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §5 replaced: the sample's
   heterozygosity must come out right.** Generate both libraries from **one genome**, so there is a
   single true heterozygosity however the reads are split, and vary only the two error rates. The fit
   must return that one value. A design that fits a heterozygosity per read group and averages will
   fail here, because the per-read-group histogram has split every site covered by both libraries
   into two shallower ones — and it will fail *worse* the more unevenly the reads are divided, so run
   at least one lopsided split (90/10) alongside an even one.

   **All three are buildable**: the first on the tomato CRAMs today, the other two from the read
   generator above, since two chemistries is two error rates in the generator.

   **The read-group axis belongs in
   [`synthetic_validation.md`](synthetic_validation.md) §4**, whose matrix has none today — the
   nearest is `samples: 1 · 2 · 50`. One axis, `read groups per sample: 1 · 2 (same chemistry) ·
   2 (different chemistry)`, and this step's acceptance is expressible there rather than in a
   private harness.
4. **The scan is fine enough that its spacing does not decide the answer.** Refit test 1's synthetic
   data at two spacings — §3's quarter-Phred and one twice as fine — and the recovered parameters
   must agree to well within what a caller could feel. This is the measurement §3 marks soft, and it
   is the only check the spacing rule needs now that no interval is emitted.
5. **Provenance is not decorative.** A parameter that borrowed from a neighbour, and one that was
   defaulted, must be distinguishable in the output from one fitted in place (§6).
6. **Memory is measured, not assumed.** Every size in §6 and in the other four documents is
   arithmetic — dense cells, a depth cap of 100, an assumed segregating rate, an assumed error rate.
   Replace all of it with a measurement of peak resident memory **and of each of §5.1's five objects
   separately**, on two runs that stress different axes:
   - **HG002 at 300×**, which stresses *depth*. §6's cell counts assume depth ≤ 100; at 300× the
     binning runs past the end of that table and the `(depth, alt-count)` cell count grows with it.
     This is where the arithmetic is likeliest to be wrong, and the depth cap is the knob if it is.
     It also inflates the census sites' sparse list, which grows with depth × error rate. **The STR
     locus table is already done**: 0.36 MB uncapped over 29,811 loci at 300× (§6), so it is the
     one object of the five that this run need not re-price.
   - **The whole tomato cohort, every sample**, which stresses *sample count*. What grows here is the
     census sites, held for every sample until the gather. The windowed histograms should **not**
     grow with the cohort — they are dropped per sample (§1.3) — so this run doubles as the check
     that they actually are: peak memory rising with sample count is a leak, and it is the cheapest
     place to catch one. Vary the number of samples walked concurrently, since that is the real
     multiplier.

   Report the objects separately rather than one total, because §9.2 turns on which of them costs
   anything. Report the windowed histogram **with and without its window key**, since §6 claims the
   windows are what `F` costs and the cells are what step 4 costs, and that split has never been
   measured.
7. **The fit is deterministic** — same accumulators, same parameters, independent of thread count.
8. **Sharded accumulation is exact.** The same sample walked in one region and in many must give
   byte-identical accumulators, for all five objects.
