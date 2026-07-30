# ng — the parameter pre-pass: what each read group tells us before genotyping

*Design spec, 2026-07-30. **No code yet — this settles the design.** Companion arch and plan docs
do not exist yet. This is ng step 4 of [`ng_proposal.md`](ng_proposal.md); the interfaces it fills
are `SampleSummarizer` and `CohortEstimator` in
[`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) (lines 343, 351). Grounded in two
research notes — [`rough_caller_alternatives_2026-07-23.md`](../../reports/research/rough_caller_alternatives_2026-07-23.md)
and [`pileup_partial_coverage_ref_fill_2026-07-27.md`](../../reports/research/pileup_partial_coverage_ref_fill_2026-07-27.md)
— and in measurements made on the tomato and HG002 cohorts during the STR stutter study.
`src/ssr/` and `src/pileup/` are frozen production: everything said about them here is a record,
not a change.*

---

## 1. Scope — goals, non-goals, and what it does not do

Before the real caller runs, something has to tell it **how noisy this data is**. Two numbers on
the generic path — how often a base is read wrong, and how often a site is heterozygous — and on
the STR path a third: how often a repeat tract gains or loses whole copies. Those numbers are the
priors the cohort model leans on, and they are properties of **how the DNA was prepared and
sequenced**, not of the genome. This step measures them.

**Goals.**

1. Estimate the per-base error rate, the variant rate, and the STR slippage behaviour **per read
   group**, from Stage-1 data alone.
2. Do it **without first calling genotypes and keeping the confident ones**, which is what
   production does and what biases it (§2).
3. Use **one estimator for both paths**, differing only in how sites are stratified. The generic
   path stratifies by nothing; the STR path by motif period and repeat count (§6).
4. Emit parameters the cohort caller can consume as frozen inputs, so genotyping stays a pure
   function of (reads, parameters).

**Non-goals.**

- **Genotyping.** Nothing here calls a variant. The estimator sums over genotypes precisely so it
  never has to pick one.
- **Replacing `F`.** The inbreeding coefficient is estimated today as `1 − Hobs/Hexp`. The research
  note argues for a runs-of-homozygosity estimator instead; that is a separate piece of work and is
  deferred with a home (§9).
- **Long-allele recovery.** GangSTR profiles insert size, coverage, GC and read length before
  genotyping. Those pay off only for alleles longer than a read, which our STR path does not
  attempt. Deferred (§9).

**It does not:**

- decide which loci exist (that is region typing, step 3);
- change anything in `src/ssr/` or `src/pileup/` — production is frozen;
- estimate anything per *sample* where a sample holds several read groups. The unit is the read
  group (§4).

---

## 2. Why production's numbers are biased

Two independent problems. Both are recorded in the research notes; the second is now measured.

### 2.1 Both rough callers threshold, then count

The generic path classifies each site as het, hom-alt or ambiguous, vetoes some hets on strand
bias, and increments one of three counters
([`src/sample_summary/het.rs:266`](../../../../src/sample_summary/het.rs), `observe_site`). The STR
path pools reads from loci that passed a confident-genotype gate
([`src/ssr/cohort/prepass.rs`](../../../../src/ssr/cohort/prepass.rs)). Both therefore estimate a
parameter from **the sites that were easy to call**, and inherit whatever the threshold selected
for.

Two further consequences of the same design, both already noted in `het.rs`'s own comments:
ambiguous sites are dropped rather than counted (`het.rs:48`), and the strand-bias veto removes
heterozygotes "broadly rather than targeting artifacts" (`het.rs:83`).

There is also a quieter loss upstream: `pileup_to_psp.rs:91` skips pure-reference columns
altogether. Those zero-alt sites are exactly the observations that pin the error rate down — a site
with depth 30 and no alt reads is strong evidence about `ε`, and it never reaches the accumulator.

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

---

## 3. The estimator — sum over the genotype, do not choose one

For a site with depth `n` and `k` reads supporting the alternative, and unknown genotype
frequencies `π`:

```text
P(site | ε, π) = π_hom_ref · ε^k (1−ε)^(n−k)
               + π_het     · (½)^n
               + π_hom_alt · (1−ε)^k ε^(n−k)
```

Multiply over sites, maximise over `(ε, π_het)`. The maximiser **is** the variant-rate estimate.
No threshold, no counting, no discarded ambiguity, and the error rate comes from the data rather
than from base qualities.

This is not a research proposal. GATK's STR calibration
(`DragstrParametersEstimator.java`, vendored) computes exactly these three terms and grid-searches
the error parameter to maximise the total; there is no rough genotyper anywhere in DRAGEN's
calibration. Three of its implementation choices are worth copying:

1. **Stratify, then pool the thin strata.** Fit per cell; where a cell has too little data, borrow
   from its neighbours rather than fitting noise.
2. **Constrain the shape across strata** instead of fitting a free parameter per cell.
3. **Grid-search rather than solve.** The likelihood is cheap and the parameter space is small.

**Why this and not a soft-EM gate.** HipSTR's EM (`em_stutter_genotyper.cpp`, vendored) also avoids
the hard gate: genotypes are latent and every read contributes weighted by its responsibility. It
is a legitimate alternative and strictly better than production's gate. It lost on scope: it
carries a per-locus allele-frequency model we do not need at this step, whereas the marginal
likelihood above gives us the two numbers we came for and nothing else. If §10's first open
question resolves against the marginal estimator, HipSTR's EM is the fallback.

**The architectural payoff.** One estimator serves both paths, differing only in the stratification
axis. That is the "SNP, indel and STR at one level" property the proposal asks for, arriving at the
step where production is most duplicated.

---

## 4. The unit is the read group, and what Stage-1 must accumulate

**Chemistry belongs to the library preparation, not to the individual.** A sample sequenced from
two libraries has two error rates and two stutter behaviours; averaging them describes neither. In
the tomato archive this is not hypothetical — one BioSample holds sixteen libraries, and one holds
three whose runs differ in read length.

ng already carries the evidence at this grain: `ObservedSequence` includes its read group as part
of its identity ([`../../../../src/ng/locus_generation/mod.rs:166`](../../../../src/ng/locus_generation/mod.rs)),
so an allele seen from two groups is two rows, each with its own quality moments. The read-group
table ([`src/ng/read/input/read_groups.rs:43`](../../../../src/ng/read/input/read_groups.rs))
resolves that identity to sample, library and experiment, and records whether each grouping name
was **declared** by the file or **synthesized** because the file gave none.

**Decision: estimate per read group; expose the fold.** Read group is the finest grain available
and the safest default. Which grain a run actually fits at — read group, library, or experiment —
is a knob, because a library sequenced across four lanes usually shares one chemistry and pooling
its lanes buys precision for free, while two libraries of one sample usually do not. The estimator
does not guess: it fits at whatever grain it is handed.

**What Stage-1 accumulates, per read group:**

| path | accumulated | why this shape |
|---|---|---|
| generic | a histogram of `(depth, alt-count)` cells | the sufficient statistic for §3; **including `k = 0` cells**, which is what production discards |
| STR | per `(period, repeat count)`: a histogram of `(reads on the modal length, reads at each whole-repeat offset, reads at a non-repeat offset)` | the sufficient statistic for the stutter model, stratified on the axis §6 shows is the right one |

Both are smaller objects than what production keeps today: three counters become a histogram, but
the histogram is sparse and bounded by depth.

**Trap.** The `k = 0` cells are the majority of the data and the reason the error rate is
identifiable at all. Whatever accumulates them must not treat "no alternative allele" as "nothing
happened here" — that is precisely the bug in `pileup_to_psp.rs:91`, inherited if the port is
mechanical.

---

## 5. The generic path — what is estimated

Two parameters per read group: the per-base error rate `ε` and the variant rate `π_het`
(with `π_hom_alt` following from it under Hardy-Weinberg, or fitted freely — §10).

**Base qualities are not the source of `ε`.** They are the instrument's own claim about itself,
recalibrated by a tool that assumed a variant catalog; the research note lists this as one of the
two refinements the literature is clearest about. `ε` is estimated from the data by §3 and the
qualities become a covariate, not the answer.

**The het model's `½` is optimistic.** The `(½)^n` term assumes both alleles are sampled equally.
Reference bias means they are not. The note recommends Bryc et al.'s reference-bias term, which
replaces `½` with a fitted per-read-group constant. **Leaning: adopt it**, since it costs one
parameter and this step already fits a grid — but it is not measured on our data (§10).

---

## 6. The STR path — what stutter actually looks like

Measured on 51 tomato read groups (8.1M observations) and on HG002 (whole genome), both with ng's
default delimiter. All figures are against each unit's own modal allele.

### 6.1 Slippage moves whole repeats

Out of every 100 reads that differ from the allele, 98 differ by a whole number of motif copies at
dinucleotides, 95 at trinucleotides, 93 at hexamers. **Homopolymers are 100% by arithmetic — every
integer is a multiple of one — and are therefore no evidence either way.**

The remainder is not slippage: it is one- and two-base indels, which the STR model has no way to
describe. That residue is small where tracts are long (3 in 1,000 reads at ≥6 repeats) and large
where they are short (§6.4).

### 6.2 One fall-off number, both directions

Reads that move usually move one repeat. Of those that moved one, the fraction that moved two
instead is the **fall-off**, and it is the same going up as going down:

| | lose (down) | gain (up) | difference vs its counting error |
|---|---|---|---|
| tomato, homopolymer | 0.065 (5,072 → 329) | 0.074 (3,592 → 265) | 1.5 SE |
| tomato, dinucleotide | 0.087 (2,438 → 211) | 0.102 (501 → 51) | 0.9 SE |
| human, homopolymer | 0.097 (3,037 → 296) | 0.115 (1,640 → 188) | 1.6 SE |
| human, dinucleotide | 0.106 (545 → 58) | 0.123 (162 → 20) | 0.5 SE |

In all four the gap is smaller than the noise on the counts. **Decision: one fall-off parameter
shared by both directions**, rejecting a separate up and down decay. Above dinucleotides the
expansion arm rests on 3 to 13 reads, so a free second parameter would be fitting noise, not
capturing a difference.

The asymmetry is real but it lives elsewhere — in **how often reads go down rather than up**, which
grows with period: 1.4× at tomato homopolymers to 4.9× at dinucleotides; 1.9× and 3.4× in human.
That is a second parameter, per period.

**The fall-off value does not transfer between datasets.** About 10 reads in 100 take a second step
in human against about 7 in tomato. The structure is portable; the number must be fitted.

### 6.3 Stratify by repeat count, not by base length

The Mark-2 spec fits the stutter level as linear in tract **length in bases**
([`../../specs/ssr_cohort_mark2.md`](../../specs/ssr_cohort_mark2.md) §4.4). The data says repeat
**count** is the better axis, which is what a per-copy slippage mechanism predicts: ten copies offer
ten chances to slip whether the copy is 2 bp or 6 bp.

At 12–15 repeats, tomato homopolymers, dinucleotides and trinucleotides stutter at 14.3%, 15.0% and
8.6% — within a factor of two of each other. On a base-length axis the same three periods at
20–29 bp were 12.9%, 12.6% and 1.3%, an order of magnitude apart. Homopolymers also become
monotonic on the copy axis and are not on the length axis.

**Decision: stratify by (period, repeat count).** Same parameter count as today, better behaved.

### 6.4 Below four repeats there is nothing to model

| repeats | share of loci | share of all slippage | slippage rate | of that, **not** a whole repeat |
|---|---:|---:|---:|---:|
| < 4 | 19.9% | 1.7% | 0.091% | **58.5%** |
| 4–5 | 28.2% | 4.5% | 0.170% | 33.8% |
| ≥ 6 | 51.9% | 93.7% | 2.006% | 0.9% |

*(HG002, known-homozygous loci. Tomato agrees in shape: 11.2% of loci, 6.0% of slippage, 19.8%
not-whole-repeat below four repeats.)*

Short tracts are a fifth of the loci, produce under 2 in every 100 slippage reads, and **nearly six
in ten of even that is not a whole-repeat change** — it is an ordinary indel, which is the generic
path's business. The STR machinery is being handed a problem it does not model.

**This bears on step 3, not on this step.** Where the boundary between the two paths sits is region
typing's decision; what this section provides is the evidence for moving it. ng's current copy
floors are `[6, 4, 4, 3, 3, 3]` for periods 1–6
([`src/ng/region_typing/segment_criteria.rs:368`](../../../../src/ng/region_typing/segment_criteria.rs)).
The measurement suggests roughly `[9, 5, 5, 4, 4, 3]`: homopolymers stay under 1 in 100 all the way
to 9 repeats in human, so their floor should rise, not fall. Recorded here, decided there.

---

## 7. Cross-cutting concerns

**Memory.** The accumulators are per read group and sparse. The generic histogram is bounded by
maximum depth; the STR histogram by (period × repeat count × offset), all small. A cohort of 50
read groups holds 50 of them — kilobytes, not the gigabytes the observation stream itself costs.

**Errors.** A read group with too little data to fit is not an error: it borrows from the pooled
stratum (§3.1) and is marked as having done so, because a parameter that came from a neighbour is
softer than one fitted in place and the consumer should be able to tell.

**Concurrency.** Accumulation is per read group and therefore embarrassingly parallel; the fit is a
cohort-wide gather, run once, single-threaded. No shared mutable state on the hot path.

**Determinism.** The grid search must be deterministic given the same accumulated histograms —
same grid, same order, no floating-point reduction whose result depends on thread count. `ε` is
frozen after this step precisely so that genotyping is reproducible.

---

## 8. Reuse over rewrite

| what | existing code | how it is reused |
|---|---|---|
| read-group identity and grouping | `src/ng/read/input/read_groups.rs` | used as-is; the fold to library/experiment is already there |
| read group on each observation | `src/ng/locus_generation/mod.rs:166` | the input to the STR accumulator |
| the three-genotype likelihood | `src/sample_summary/het.rs:146` (`SiteCounts::from_record`) | the three binomials are already computed; they are **added** instead of compared |
| stratify-and-pool, shape constraint | GATK `DragstrParametersEstimator.java` (vendored) | algorithm copied, not code |
| step-4 interfaces | `../arch/ng_step_interfaces.md` (343, 351) | `SampleSummarizer` per read group; `CohortEstimator` fits |

**No parity oracle.** This is not a port: production's estimator is the thing being replaced, so
agreeing with it would be failure. §11 says how correctness is shown instead.

---

## 9. Deferred, with a recommended home

- **`F` from runs of homozygosity.** `1 − Hobs/Hexp` applies one cohort-wide `Hexp` to every
  sample regardless of ancestry, and a uniform floor of false hets biases every sample downward.
  F_ROH is closer to the definition for a selfing crop. **Home:** its own spec — it is a
  per-genome-structure estimator, not a per-read-group chemistry one, and mixing them would make
  both harder to read.
- **Insert size, coverage, GC and read-length profiling.** Needed only for alleles longer than a
  read. **Home:** whichever spec takes on long-allele recovery.
- **HDplot and other cohort artifact signals.** A cohort-level signal for collapsed paralogs.
  **Home:** the artifact-filter step (step 10).
- **Per-locus refinement of the stutter parameters.** This step produces the prior; Mark-2 §5
  already specifies refining it per locus inside the EM. **Home:** stays there.

---

## 10. Open questions

1. **Does the marginal estimator beat production's gate on our data?** — OPEN. The research note
   proposes the experiment: fit both on the same cohort and compare the recovered `ε` against a
   held-out truth set. *Leaning:* yes, on the strength of §2.2 — but the measurement that settles
   it is HG002, where truth genotypes exist, comparing the fitted stutter rate against the rate
   measured on known-homozygous loci. **Confirm before code.**
2. **Should the fall-off depend on the level?** — OPEN. Read groups that stutter more also decay
   more slowly: at tomato dinucleotides, level against the one-step share ran ρ = −0.69, and that
   survived removing real alleles. If real, the fall-off is a function of level rather than a free
   parameter, which keeps one number per group while fitting better. *Leaning:* model it as a
   function of level. **Settled by:** repeating that correlation on human read groups — which needs
   a multi-read-group human cohort we do not currently have.
3. **Free `π_hom_alt`, or Hardy-Weinberg from `π_het`?** — OPEN. HW is one fewer parameter and is
   wrong for a selfing crop. *Leaning:* fit both and let the cohort decide, since the grid is cheap.
4. **Adopt the reference-bias term in place of `½`?** — OPEN. *Leaning:* yes; costs one parameter.
   **Settled by:** whether the fitted value departs from ½ by more than its standard error on real
   data.
5. **Why do hexamers put more mass at −3 than at −1?** — OPEN, and unexplained. Either a real
   long-tract behaviour or an artefact of tract delimitation. It breaks the geometric assumption
   for that period alone. **Settled by:** the synthetic validation
   ([`synthetic_validation.md`](synthetic_validation.md)), which can inject known hexamer alleles
   and see whether the delimiter reproduces them.
6. **Is the low whole-repeat fraction at tetra and penta real?** 62% and 53% in tomato, on 464 and
   131 reads. Both are thin and both are dominated by 3-copy tracts, which §6.4 would route
   elsewhere anyway. *Leaning:* it disappears once the copy floor rises. **Settled by:** re-running
   §6.1 with the floors of §6.4 applied.

---

## 11. How we know it works

1. **The estimator recovers known parameters from synthetic data.** Simulate reads at a known `ε`,
   known het rate and known stutter behaviour; the fit must return them. This is the primary test
   and it does not need real data.
2. **It agrees with truth where truth exists.** On HG002, the stutter parameters fitted by §3 must
   match those measured directly on known-homozygous loci (§2.2) — the 2.0% at ≥6 repeats and the
   3.4× direction split, within the fit's own error. **This is the test production's estimator
   fails**, and the reason for the whole design.
3. **It is stable across grains.** Fitting per read group and then pooling lanes of one library
   must give the same answer as fitting that library directly, within error. If it does not, the
   grain is carrying something other than chemistry.
4. **Thin strata do not produce wild parameters.** A cell with 20 reads must borrow, and be marked
   as having borrowed, rather than return an extreme value.
5. **The fit is deterministic** — same histograms, same parameters, independent of thread count.
