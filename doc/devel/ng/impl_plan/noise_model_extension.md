# The generic noise model's second site class — implementation plan

**Status:** approved by the owner 2026-08-10, on the evidence in
[`../research/noise_model_overdispersion_2026-08-10.md`](../research/noise_model_overdispersion_2026-08-10.md).
**Position:** a new milestone in
[`parameter_prepass_generic.md`](parameter_prepass_generic.md), **between Milestone F and
Milestone G**. F is complete; G1 is held until this lands, because the anchors G1 asserts are
the ones this changes.

## Why

`spec/parameter_prepass_generic.md` §2 gives the generic path one error rate per read group and
nothing else. Measured on HG002's confident regions, that distribution's **body is right and
its tail is wrong**: 818 loci with no benchmark variant carry three or more alternative reads
where the model predicts 29. The three-genotype mixture has one class that can explain such a
site, so **fitted heterozygosity comes back 1.41 times the benchmark count** — 776 sites where
the truth has 550.

A second site class — clean with probability `1 − w`, noisy with probability `w`, the genotype
emission using that site's rate — cuts that to **1.09**, beats a beta-binomial by 425 nats for
one further parameter, and is found independently at 30× and at 300×. It also halves the error
rate's drift across depth, from +6.1% to +3.0%, which is the property Milestone G's coverage
sweep exists to test.

**The bias is an absolute number of sites, not a percentage** — about 230 on this depth
distribution, whether the sample is 0.1% or 1% heterozygous. It is therefore worst on the
low-heterozygosity samples this caller is aimed at, which is why it is worth a milestone.

## What this is not

- **Not a change to the STR path.** The site class lives in `generic/noise_model.rs`. The
  `NoiseModel` trait's `NoiseParams` is an associated type, which is the seam that makes this
  possible without touching `fitting/`'s signature or the STR path's own model (arch §4).
- **Not a fix for the homozygous-non-reference rate**, which stays about 7% below the benchmark
  and whose cause is unknown.
- **Not a claim the tail is now right.** A residual of 9% at 30× and 14% at 300× remains.

## The decision this plan takes, and why it is stated here rather than assumed

**What a sample emits as its error rate stays one number, and that number is the
share-weighted marginal** `(1 − w)·ε_clean + w·ε_noisy` — the probability that a read
disagrees with the reference at a site drawn at random. Three reasons. It keeps
`Estimate<ErrorRate>` and every consumer of it unchanged. It is the quantity the model-free
count measures, so the G1 anchor still applies — measured, it is 2.344 × 10⁻³ against a
model-free 2.263 × 10⁻³, **3.6% high and inside one ladder rung**. And emitting `ε_clean`
instead would put the reported rate **16% below** the model-free count, which is the one
comparison arch §9 calls an unambiguous bug.

The pair travels beside it as a diagnostic, so a later consumer that wants to score a read
against its own site class can, without another fit.

---

## The steps

### N1. The types, and the emitted rate.  ✅

`SiteNoise { noisy_fraction, noisy_error_rate }` as constrained newtypes beside the existing
ones; `SampleLibraryNoise` gains it. `GenericSampleParameters` keeps `error_rate` as the
share-weighted marginal and gains `site_noise` as the diagnostic pair.

**One question this step must not answer silently.** `w` and `ε_noisy` are fitted **per sample,
shared across read groups**, while `ε` stays per read group. Both cohorts hold single-library
samples, so no data distinguishes a per-sample from a per-library noisy class, and the
per-sample choice is the one that keeps the ladder one-dimensional. **Recorded as an assumption
to revisit when a multi-library alignment exists**, in the type's own doc.

*Tests:* the newtypes reject out-of-range; the marginal is computed once and agrees with a
hand-worked example.

### N2. The scoring rule.  ✅ **Own commit, do not bundle.**

`ln L(cell | θ)` becomes a convex combination of the existing rule evaluated at `ε_clean` and at
`ε_noisy`. The multi-library closed form of spec §5.1 factors cleanly, because the site class is
a property of the **site** and the library split is a property of the **reads**: the sum over
"which library produced each alternative read" happens inside each branch.

**The silent failure this isolates:** a convex combination of two correct rules is still a
probability over the cell space, so none of the identities *by itself* can see the branches
being weighted wrongly, or `ε_noisy` being applied to the wrong branch.

*Oracle:* all four of D2's identities re-run — the rule sums to one over the cell space at any
parameters, no cell is charged a negative count of reference reads, every library's rate equal
reproduces the exact per-library likelihood, and agreement with `ng_multilib_key_harness.rs` to
floating point at `w = 0`. Plus a fifth that is new and is what catches a mis-weighted branch:
**at `ε_noisy = ε_clean` the rule must reproduce the one-class rule exactly, at every `w`.**

### N3a. The oracle's world: a measured depth distribution, generated in closed form.  ✅

 and  in
. Landed ahead of the fit because the histogram came from a walk
over an alignment no worktree carries and could not be re-derived from anything in the
repository; everything else in N3 can.

### N3b. Fitting the two new parameters.  ☐ **Own commit, do not bundle.**

Expectation-maximisation for `w` and `ε_noisy` inside the coupled loop, beside the ladder scan
that settles each library's `ε`. Multi-start over the separation between the two classes, not
over `w` alone — **the same trap the inbreeding fit hit**, where starts that disagreed only
about how much of the genome was inside a run returned `F` = 0.0000, converged and silent
(research note 2026-08-06 §3.4).

**The silent failure this isolates:** two classes that collapse into one, or swap, both return
a plausible pair and report convergence.

*Oracle:* the two exact controls of the research note, both computed in closed form over the
real depth distribution with no sampling noise in them.
- **On a world with one error rate the genotype frequencies must be unchanged**: measured, the
  extension returns heterozygosity **−0.0006%** and the homozygous-non-reference rate
  **+0.0000%** from the truth, while splitting off a spurious `w` = 0.48% and reporting a clean
  rate **1.10% low**. That last number is the price of the extension and the test pins it, so a
  later change that makes it worse is visible.
- **On a world that has a noisy class it must be recovered**: `w` to four decimal places,
  `ε_noisy` to 0.02%, both frequencies to 0.000%.

### N4. The harnesses.  ☐

`ng_multilib_key_harness.rs` and `ng_inbreeding_harness.rs` extended to the new model, and
**E2's 25 worlds re-run**. The gate is that the coupled fixed point is still the truth in all 25
from a start at three times the true rates, and that `F` still recovers a drawn genome's
realised autozygous fraction to four decimal places under a false-heterozygote floor. Any world
that moves is a finding, not a tolerance to widen.

### N5. The anchors, re-measured.  ☐

F3's four real-alignment tests re-run on all five alignments, and G1's model-free comparison
computed before and after, on HG002 at 30× and 300× and on three tomato CRAMs. **The deliverable
is the before-and-after table**, including the tomato samples, which have no truth set and
therefore test only that the fit does not degenerate — `w` between 0 and 1, the two classes
separated, the fit converging.

> **Checkpoint N:** the scoring rule passes five identities, the fit is proven unbiased on a
> world that does not need it and exact on one that does, both harnesses agree, and
> heterozygosity on HG002 has moved from 1.41 to about 1.09 times the benchmark. Pause for
> review. **G1 resumes after this.**

---

## What the owner still owns

The spec and architecture describe the one-rate model and this plan does not edit them:

- `spec/parameter_prepass_generic.md` **§2** — *"a per-base substitution rate and nothing
  else"* is what this milestone changes.
- `spec/parameter_prepass_generic.md` **§5.1** and `arch/parameter_prepass_generic.md`
  **§5.1** — the multi-library closed form gains its convex combination.
- `arch/parameter_prepass_generic.md` **§2.4** — `GenericSampleParameters` gains a field.
- `arch/parameter_prepass_generic.md` **§8** — the reference-bias `OPEN:` item is adjacent to
  this and may be subsumed by it; worth re-reading once N3 lands, because a noisy site class
  and a reference-bias term compete to explain some of the same reads.

Recorded here rather than edited, per the project rule.
