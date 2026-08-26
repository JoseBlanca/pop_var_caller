# ng — how the caller learns how variable the population is

**Status:** 2026-08-26. **No code yet — this settles the design.** It supersedes part of
[`calling_priors.md`](calling_priors.md) §5, which specifies the repeat tract's prior seed as a
constructed geometric shape scaled by a measured diversity; §4.2 below replaces both halves of that
construction. Nothing has been built against §5's version — `fill_ssr_seed` exists and its two
run-level inputs have no producer — so this is a change of design, not a supersession anyone has to
do archaeology on.

**Companions:** [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3 (the ordinary-site
diversity, and the STR one it delegates), [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §6
(why the repeat number is not fitted there),
[`calling_priors.md`](calling_priors.md) §4–§5 (what the caller does with either),
[`str_slippage_level_curve.md`](str_slippage_level_curve.md) §5.1 (the curve-and-ladder device this
document borrows its discipline from).

---

## 1. What this is

**Before a caller looks at a read, it needs some idea of how variable the population is.** Without
one it cannot judge whether a rare allele is surprising. Every cohort caller has such a number; ours
has two, because a repeat tract and an ordinary site are not the same population question, and this
document settles **where each comes from, what happens when there is too little data to measure it,
and how it reaches the caller**.

### Goals

1. **Both numbers are measured on the run's own data** wherever the data supports it, rather than
   taken from a constant.
2. **An answer at one sample and at a thousand.** Neither number may simply be absent at the small
   end; what changes is which rung of a stated ladder the answer came from.
3. **The rung travels.** A call resting on a measured number and one resting on a stated constant
   must be distinguishable in the run's output without re-running anything.
4. **The two numbers cannot be confused for one another.** Today's repeat path applies a SNP-scale
   constant to tracts (`SFS_THETA = 0.01`, `src/ssr/cohort/freebayes_emit.rs:42`); that is the
   failure this separation exists to stop recurring.

### Non-goals

- **Fitting the ordinary-site frequency spectrum.** That is the cohort gather's own content and is
  much the largest piece here; §3.4 says what it buys and defers it with a home. This document
  specifies the *diversity*, which is the rung below it, and how both reach the caller.
- **Changing what the caller does with either number once it has it.** The projection from a
  spectrum to a seed pair is built and tested (`project_spectrum_seed`,
  `src/ng/calling/genotype_prior/seed_generic.rs:636`); this document does not touch it.
- **Repeat-tract candidate selection.** Which lengths a tract is called over is
  [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)'s.

### What it does not do

- It does not introduce a per-locus diversity. Both numbers are per run or per stratum; nothing
  here is fitted per site or per tract.
- It does not add a knob. Every fallback is a named constant with its source, not a setting.
- It does not decide the *encoding* of anything on disk.

---

## 2. Why one number will not do

**A repeat tract mutates orders of magnitude faster than a base does**, so the population's
variability at tracts is a different quantity from its variability at ordinary sites — not a
correction to it (`parameter_prepass_cohort.md` §3). A consumer that applies the ordinary-site
number to a tract badly understates how many alleles to expect there, and the current caller does
exactly that.

The two also differ in **shape**, and that is what makes them different design problems:

| | ordinary site | repeat tract |
|---|---|---|
| the allele set | two, usually | many, indexed by length |
| what "diversity" is | how often two copies differ | how often two copies carry different **lengths** |
| what it is fitted **across** | samples, at each site | **tracts, within a stratum** |
| what one sample gives | that genome's own heterozygosity ÷ (1 − F) | thousands of tracts per stratum |

**That third row is the whole of why the small-cohort answers differ**, and it is not obvious: the
ordinary-site spectrum needs a panel because a frequency spectrum has no shape without one, while a
stratum's length spectrum is estimated from every tract of that shape in the genome. Tomato holds
462,701 kept tracts in 141 strata (`src/ng/parameter_estimation/joint/ssr_fit.rs:18`), so a single
genome still puts thousands of tracts behind each stratum.

---

## 3. The ordinary-site path

### 3.1 What is measured

**Expected heterozygosity: how often two copies of an ordinary site, drawn at random from the panel,
differ.** It is a property of the population, so an individual's inbreeding is divided out before it
is a population statement (`parameter_prepass_cohort.md` §3):

```text
Hobs = Hexp · (1 − F)          so      Hexp = mean over samples of  Hobs(sample) / (1 − F(sample))
```

**Both inputs are already fitted, per sample.** Each sample's genotype frequencies at each ploidy
are `SampleRates` (`src/ng/parameter_estimation/generic/mod.rs:343`), carried on
`GenericSampleParameters` beside that sample's inbreeding coefficient (`:446`, `:450`). What does
not exist is the fold across samples, because the step that owns it — the cohort gather — has a
specification and no code.

**Build it from observed heterozygosity, not from the non-reference rate.** A site where every
sample is homozygous for the alternative is not polymorphism; it is a place where the reference
accession carries the odd allele. Estimating diversity from "how often we see a non-reference
allele" counts every quirk of the reference as cohort variation
(`parameter_prepass_cohort.md` §3).

**⚠ A constraint this places on an open choice elsewhere.** The inbreeding coefficient must *not*
come from the ratio estimator `F = 1 − Hobs/Hexp`, because that estimator needs an expected
heterozygosity to produce its answer, and feeding it back returns whatever was assumed. The
runs-of-homozygosity estimator has no such problem — it reads inbreeding off the genomic
distribution of heterozygosity and never needs a population expectation. Which estimator ships is
open in [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §11; **this section is a
constraint on that choice**, not a re-opening of it.

### 3.2 One sample

**The formula returns that genome's own observed heterozygosity divided by its own `(1 − F)`, and
that is a measurement rather than a fallback.** One diploid genome carries two copies of every site,
and how often those two differ is exactly what the quantity asks; the inbreeding correction is what
turns a statement about the individual into one about the population it was drawn from. So **this
number is fitted at every cohort size down to one**, which the frequency spectrum is not.

### 3.3 The ladder, and what each rung means

The consumer already implements three regimes
(`src/ng/calling/genotype_prior/seed_generic.rs:597`), and this document changes none of them — it
supplies the inputs they were written for:

| rung | when | what the seed is |
|---|---|---|
| **fitted spectrum** | a panel large enough for one | shape and scale both from the spectrum; the diversity is not read |
| **fitted diversity** | no spectrum | a neutral shape at the measured diversity |
| **stated constant** | neither | a neutral shape at `ExpectedHeterozygosity::SPECIES_FALLBACK` |

**The middle rung is the one sample's rung**, and it is why the diversity is worth fitting even
after the spectrum exists: at one sample there is never a spectrum. **This corrects a reading that
the spectrum makes the diversity redundant** — true only where a spectrum exists.

**⚠ One sentence in the consumer's own documentation contradicts this** and should be re-derived
rather than repeated when that file is next touched: `seed_generic.rs:604` says *"a cohort of five
arrives here without one while a single sample arrives with one"*, which has the two cohort sizes
the wrong way round against `parameter_prepass_cohort.md` §3. Nothing depends on it — the code
branches on whether a spectrum arrived, never on cohort size — so it is a wrong sentence rather than
a wrong behaviour.

### 3.4 Deferred: the frequency spectrum

**What it buys is the *shape* of variation, on panels large enough to have one**, and it is the
larger half of this subject by a wide margin: fitting allele-count classes across a panel is the
cohort gather's real content. **Deferred, with a home:**
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4, which specifies it.

**What deferring it costs, stated plainly:** with the diversity alone, the ordinary-site prior moves
off a species-range constant onto a number measured on this cohort and keeps a neutral shape. That
is most of the benefit for a small fraction of the work, which is why the two are separated here
rather than built together.

---

## 4. The repeat-tract path

### 4.1 What is already measured, and it is not a scalar

**The joint repeat fit already estimates, per stratum, exactly the object a tract's prior needs**
(`src/ng/parameter_estimation/joint/ssr_fit.rs:281`):

- **the length spectrum** — how that stratum's chromosomes are spread over tract lengths, indexed in
  whole repeat units either side of the reference tract length (`:289`);
- **the concentration** — how monomorphic the stratum's tracts are; small means most tracts sit at
  one length while the stratum as a whole spans many (`:291`).

Together they are a Dirichlet over a tract's length frequencies: **a shape and a strength, both
fitted**. The fit conditions on each sample's homozygote excess as it goes (`:59`), so inbreeding is
divided out inside the estimator rather than afterwards.

**This has been run end to end on both benchmark cohorts**
([`../reports/str_slippage_curves_on_both_cohorts_2026-08-21.md`](../reports/str_slippage_curves_on_both_cohorts_2026-08-21.md)).

### 4.2 Decision: seed the tract prior from the fit, not from a constructed shape

**Decided.** The tract's prior seed is the stratum's fitted length spectrum and concentration,
mapped onto the locus's candidate lengths.

**What it replaces.** `calling_priors.md` §5 specifies a shape *constructed* as a geometric decay
away from the cohort's commonest length at the tract, scaled to reproduce a separately measured
cohort-wide repeat gene diversity (`fill_ssr_seed`,
`src/ng/calling/genotype_prior/seed_ssr.rs:200`). Both halves go.

**Why the fitted pair wins, and these are facts about the two constructions rather than a
preference:**

- **The constructed shape has one free parameter and the fitted one has none.** The decay is a
  single number per group of loci with a coded fallback of 0.5 (`src/ng/types.rs:810`); the fitted
  spectrum is a distribution estimated from that stratum's own tracts.
- **The constructed version has a failure mode the fitted one does not have.** Scaling a shape to
  reproduce a measured diversity is only possible below a ceiling the shape itself sets, and
  `SsrSeedOutcome::DiversityUnreachable` is what happens above it
  (`src/ng/calling/genotype_prior/seed_ssr.rs:85`). **At one outbred sample that is every tract** —
  a single diploid shows at most three lengths, whose shape can imply at most 0.625, against the
  ~0.72 repeat diversity HG002 actually has. A fitted Dirichlet asserts no such scaling and cannot
  fail this way.
- **It removes a per-locus input that has no source.** The constructed shape needs the cohort's
  commonest length *at this tract*, which is cohort-derived and would come from repeat-tract
  candidate selection — unwritten. The fitted spectrum is indexed by offset from the **reference**
  tract length, which every locus already knows.
- **It removes the need for a cohort-wide repeat diversity number entirely**, which nothing emits.

**What it costs.** The prior's belief becomes per stratum rather than per locus: a tract whose own
commonest length differs from its stratum's centre is not distinguished from one that sits on it.
**Unmeasured**, and §9's open question 1 says what would settle it.

**The alternative that was live** was to build the cohort-wide repeat diversity in the cohort
gather and keep §5's construction. It loses on all four points above and costs an unbuilt subsystem
first.

### 4.3 One sample is not the thin case here — a thin stratum is

**The stratum's spectrum is fitted across tracts**, so cohort size is not what makes it thin. The
fit's own refusal floor is measured in **tracts, and is 8** of them
(`src/ng/parameter_estimation/joint/ssr_fit.rs:650`), chosen from draws run deliberately at both
ends of this caller's range: at 8 tracts, 3 fits in 100 collapse on a single deep sample and none on
a 63-sample cohort.

**So the small-cohort answer on this path is: nothing special happens.** A single genome carries the
same tracts as a panel does; what it does not carry is many *samples* per tract, and the estimator
does not need them.

### 4.4 The fallback, and why it is not a curve

**A stratum too thin to fit is furnished from its motif period's slippage curves, and such a stratum
carries no length spectrum and no concentration at all** (`ssr_fit.rs:420`). So the ladder needs a
rung below the stratum's own fit.

**The obvious move is the device this project already chose for slippage** — a curve per motif
period through every stratum, each weighted by how precisely it holds its own answer, with a
stratum that has nothing taking the curve whole, and three further rungs below that so a curve
always comes back (`src/ng/parameter_estimation/joint/share_curve.rs:1-37`). **The principle
transfers and is adopted: always answer, from data where there is data, and say which rung.** The
machinery does not, and the reason is a measurement.

**How much this rung actually carries**, from the two cohorts' real tables:

| | strata with their own fit | strata furnished from curves (no spectrum) | strata with nothing |
|---|---:|---:|---:|
| HG002 (deep, one sample) | 79 of 117 | **38** | 15, holding 36 tracts |
| tomato (3×, 63 accessions) | 17 of 39 | **22** | 10, holding 24 loci |

The stratum counts look alarming and the **loci counts do not**: on HG002 the strata that gained a
parameter set at all hold 529 tracts, 2% of the run's 27,399; on tomato the whole gain from
furnishing was 280 loci of 3,965, 7%. **A thin stratum is thin because few loci sit in it.**

**Contrast with the case the curve was built for.** The two slippage shares had a 4,000-slipped-read
floor that only one motif period in twelve ever cleared, so **69 of HG002's strata and every one of
tomato's got nothing at all** (`share_curve.rs:10-15`). That is why a curve was worth its machinery
there. Here most strata have their own spectrum and the gap is a few percent of loci.

**Decided: the rung below a stratum's own fit is its motif period's pooled tracts** — one fit over
every tract of that period, giving a spectrum and a concentration in the same form. Three reasons,
and the first is the one that matters:

- **proportionate to what it carries** (2–7% of loci, above);
- **it stays a coherent distribution.** A curve through a *distribution* means one curve per length
  class, refitted and renormalised — the classes are not independent, so the renormalisation is an
  approximation. A pooled fit is a real distribution by construction;
- **it is one mechanism rather than two**, since the concentration comes out of the same pooled fit.

**What it gives up:** the repeat-count trend within a period. A longer tract spreads over more
lengths, and pooling flattens that. Bounded by the loci it applies to; §9's open question 2 says
what would revisit it.

**The ladder, in full:**

| rung | when | provenance |
|---|---|---|
| the stratum's own fitted spectrum and concentration | ≥ 8 tracts in the stratum | fitted here |
| its motif period's pooled spectrum and concentration | the stratum has no own fit | borrowed |
| a stated flat spectrum over the reachable lengths at a named concentration | the period has no fitted stratum either | defaulted |

**The bottom rung's two values are soft and are marked so**: a flat shape asserts no belief about
which length is likelier, and the concentration below it is a stated constant with no measurement
behind it. Naming it is what makes it movable.

---

## 5. How both reach the caller

**Decision: the run-wide bundle the caller already assembles, present or absent as a whole.**

The caller takes one `FrozenParameters` per run (`src/ng/calling/mod.rs:567`), gathered once and
handed to every locus. It already carries two run-level repeat-tract parameters — the slippage
lookup and the substitution-rate map — so this follows a precedent rather than setting one.

**Absent-or-present as a whole**, exactly as the contamination views are: a run whose fit produced no
repeat-tract parameters carries none, and **a tract in such a run is refused by name** rather than
seeded from the ordinary-site numbers, which are indexed by allele rather than by length and would
be meaningless there. *Absent* and *a fitted zero* are different claims about a run.

**The alternative was to attach the tract numbers to each locus's evidence.** It loses because they
are run-level facts: every locus of the run would repeat them, and `FrozenParameters`' own
documentation is explicit that mixing per-run with per-locus and per-sample values is the confusion
it exists to prevent.

**What has to change at the seam.** `StratumFits` is the one lookup that crosses into calling
(`src/ng/parameter_estimation/joint/stratum_fits.rs`). It gathers each stratum's slippage numbers
and their provenance and **drops the length spectrum and the concentration**, which `fit_strata`
produced. They must be carried, keyed the same way, with the rung recorded beside them.

**The ordinary-site side needs no seam work.** `RunParameters::project_seed`
(`src/ng/calling/run_parameters.rs:97`) already takes the spectrum and the diversity as arguments;
only the producer is missing.

---

## 6. Cross-cutting concerns

**Cost.** Both numbers are computed once per run. The ordinary-site diversity is a fold over
per-sample values already in memory. The repeat numbers are already computed; carrying them across
the seam adds one spectrum and one scalar per stratum — tens to a couple of hundred strata per run.

**Errors.** Neither number can fail a run. Every absence has a rung below it and the rung is
reported; the only refusal in this document is a repeat tract in a run with no repeat-tract
parameters at all, which is a run that fitted no tracts.

**Concurrency.** Nothing here is per locus or per worker; both are frozen before calling starts and
read-only thereafter.

---

## 7. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the stratum's length spectrum and concentration | `joint/ssr_fit.rs:289`, `:291` | read as-is; nothing re-fitted |
| the seam into calling | `joint/stratum_fits.rs` | widened to carry two more values and their rung |
| per-sample genotype frequencies and inbreeding | `generic/mod.rs:343`, `:450` | the two inputs of the ordinary-site fold |
| the three-regime seed projection | `seed_generic.rs:636` | unchanged; this supplies its inputs |
| the always-answers ladder discipline | `joint/share_curve.rs` | the principle, not the curve machinery (§4.4) |
| the run-wide bundle | `calling/mod.rs:567` | two more fields, absent-or-present |

**Parity oracle.** None: neither number exists in production to compare against. What replaces it is
§8's first two checks.

---

## 8. How we know it works

1. **The repeat prior reproduces the fit.** Seeded from a stratum's own spectrum and concentration,
   the prior's implied length distribution matches the fitted one — a property of the mapping, not
   of any cohort.
2. **One sample no longer fails at every tract.** The current construction refuses every tract at
   one outbred sample (§4.2); the replacement is checked on a single-sample fixture and refuses
   none.
3. **The rung reaches the output.** A locus seeded from a stratum's own fit, one from its period's
   pool, and one from the stated constant are distinguishable in the run's record.
4. **The ordinary-site fold is checked against a hand-computed cohort** at one sample and at
   several, including a sample with a non-zero inbreeding coefficient, where the uncorrected and
   corrected answers differ.
5. **The two numbers cannot be crossed.** They are separate types; a test that hands one where the
   other belongs must not compile.
6. **The end-to-end check** is a repeat tract called from real evidence, which is what this unblocks.

---

## 9. Resolved decisions & open questions

**Resolved.**

- *Where does the tract prior's belief come from?* **The fitted per-stratum Dirichlet**, not a
  constructed geometric shape scaled by a cohort-wide diversity (§4.2). Rejected because the
  construction has a free parameter with a coded fallback, fails at every tract at one sample, and
  needs two inputs nothing produces.
- *Is the ordinary-site diversity redundant once the spectrum exists?* **No** — it is the
  one-sample rung, and at one sample there is never a spectrum (§3.3).
- *Should a thin stratum borrow a neighbouring stratum?* **No.** That mechanism was deleted from
  this project on 2026-08-20 along with the floor and the copy rule; the replacement is a rung that
  always answers and says which rung it was (§4.4).
- *Should the thin rung be a curve, like slippage?* **No — the principle yes, the machinery no**,
  because the gap is 2–7% of loci where slippage's was most strata on both cohorts (§4.4).

**Open.**

1. **Does a per-stratum belief cost anything against a per-locus one?** The prior's shape is now the
   stratum's, so a tract whose own commonest length sits off its stratum's centre gets a belief
   centred elsewhere. **Leaning: it does not matter enough to act on**, because the shape is a prior
   the reads move. **What would settle it:** the census holds per-locus length evidence, so the
   spread of per-locus modes within a stratum is measurable directly. Confirm before relying on the
   per-stratum grain for anything but the prior.
2. **Should the thin rung keep the repeat-count trend?** Pooling a period flattens it, and longer
   tracts genuinely spread over more lengths. **Leaning: leave it** until the loci it applies to are
   worth more than 7%. **What would settle it:** the spread of the pooled spectrum against repeat
   count within one period, on tomato, where the thin rung carries the larger share.
3. **What concentration does the bottom rung state?** A flat shape has an obvious answer; its
   strength does not. **Leaning:** the run's own median fitted concentration where the run fitted
   any stratum, and a stated constant only where it fitted none. **Confirm before code.**

**Deferred, with a recommended home.**

- **The ordinary-site frequency spectrum** —
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4 (§3.4 above).
- **Where the ordinary-site diversity fold is computed.** It needs a home, and the cohort gather is
  the specified one — but that step has no architecture document and no code, and the fold itself is
  small. **Recommended home:** the cohort gather when it is built; until then, wherever the run
  assembles `RunParameters`, marked as a temporary lodging.
- **Reconciling the per-stratum length spectra with the per-sample allele spectra** the per-sample
  STR route fits — they are the same kind of object at different grains
  (`parameter_prepass_ssr.md` §6). **Home:** `parameter_prepass_cohort.md` §3.
