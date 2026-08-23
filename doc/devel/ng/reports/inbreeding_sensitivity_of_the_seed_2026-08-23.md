# How much does the caller's starting belief depend on getting inbreeding right?

*2026-08-23. Harness: [`examples/ng_inbreeding_sensitivity.rs`](../../../../examples/ng_inbreeding_sensitivity.rs),
run in release in the dev container. Prompted by the owner's question at Checkpoint D of
[`calling_prior.md`](../impl_plan/calling_prior.md): the pre-pass measures a panel's diversity
through its heterozygotes, and in a selfer most of the alternative allele is not in heterozygotes.*

## The question

The genotype prior starts every locus from two numbers: one attached to the reference base, one
shared among the alternatives. The panel's inbreeding coefficient `F` reaches those two numbers by
two separate routes, and the worry is about the second:

1. **Through the fit.** Working out which pair of numbers produced the panel's allele counts needs
   a model of how a panel's chromosomes come about, and in a selfer an individual's two copies are
   usually one ancestral copy inherited twice. That model carries `F`.
2. **Through the diversity.** Where the pre-pass emits no allele-count distribution — one sample,
   or a cohort below its panel-size floor — the pair falls back to `(1, θ)`, and `θ` is each
   sample's observed heterozygosity divided by `(1 − F)`
   ([`calling_priors.md`](../spec/calling_priors.md) §4). **At an inbreeding coefficient of 0.85
   about 15 alternative copies in 100 sit in heterozygotes**, so diversity is measured through the
   thinnest channel the panel has and then multiplied by 6.7 to recover it.

## The answer, in one line

**At cohort scale the fit barely notices; at one sample the exposure is real, and it is entirely in
the diversity rather than in the fit.**

## What a wrong `F` does to the fit

The panel really is at `F true`; the fit is told `F used`. The truth in every row is a reference
number of 1 and an alternative total of `θ = 6 in 10,000`, tomato's fitted diversity. The last
column is the prior odds on a heterozygote against a homozygous-variant call before any read,
computed by running the shipped prior — 2:1 is what a reference number of 1 gives.

| individuals | `F` true | `F` used | reference number | alternative / θ | het : hom-alt |
|---|---|---|---|---|---|
| 26 | 0.85 | 0.75 | 0.9775 | 0.972 | 1.954 : 1 |
| 26 | 0.85 | 0.80 | 0.9874 | 0.984 | 1.974 : 1 |
| 26 | 0.85 | **0.85** | **1.0013** | **1.001** | **2.001 : 1** |
| 26 | 0.85 | 0.90 | 1.0115 | 1.015 | 2.022 : 1 |
| 26 | 0.85 | 0.95 | 1.0243 | 1.032 | 2.047 : 1 |

**A coefficient wrong by 0.10 moves the reference number by 2.4% and the prior odds by 2.3%.** At
63 individuals it is smaller still — 1.6% for the same error. At `F = 0.6` rather than 0.85, 2.0%.

**That is the number that matters for a cohort**, because the projection is handed *one*
coefficient standing for samples that each have their own, and no document yet says how that one
is arrived at (raised separately at Checkpoint D). Whatever aggregation is chosen, a sample
departing from it by 0.10 costs the panel's starting belief about one part in forty.

## What a wrong `F` does at one sample

The same table at one individual looks alarming and must not be read that way:

| individuals | `F` true | `F` used | reference number | het : hom-alt |
|---|---|---|---|---|
| 1 | 0.85 | 0.80 | 0.600 | 1.199 : 1 |
| 1 | 0.85 | **0.85** | **1.001** | **2.001 : 1** |
| 1 | 0.85 | 0.90 | 3.010 | 6.009 : 1 |
| 1 | 0.85 | 0.95 | 996.3 | 1237 : 1 |

**Those rows hold the coefficient wrong at one end only, and at one sample that cannot happen.**
No site can vary *across* a panel of one, so the distribution the fit is handed is the pre-pass's
own neutral prior, built from the same coefficient the fit then uses. A wrong coefficient is wrong
at both ends and the two errors cancel exactly in the reference number, which comes back at 1
whatever `F` is.

What survives at one sample is route 2, and it is exact arithmetic rather than a fit:

| `F` | alternative copies in heterozygotes, of 100 | error in `θ` from `ΔF = 0.05` |
|---|---|---|
| 0 | 100 | 5% |
| 0.6 | 40 | 12% |
| 0.85 | 15 | **33%** |
| 0.9 | 10 | 50% |

## What I conclude

- **No change to the diversity estimator is warranted for a cohort.** Wherever the pre-pass can
  emit an allele-count distribution, the fit reads both numbers off it — counting every
  alternative copy, homozygous ones included — and `θ` enters only through the neutral shape the
  estimate is held toward, diluted by however much the real sites outweigh it (3,100 to 1 in
  aggregate on tomato, 39 to 1 in the thinnest class).
- **The single-sample case is where the sensitivity lives, and dividing by `(1 − F)` is not
  avoidable there.** A single genome carries no information about how common an allele is in a
  population; its heterozygosity, corrected for its own inbreeding, is the whole signal. The
  obvious alternative — measuring diversity from non-reference homozygous sites as well — is the
  non-reference rate spec §4 rejects, because on a crop reference it counts the reference
  accession's own private alleles as population variation.
- **So the honest response is to report it, not to repair it.** A single-sample run at a high
  inbreeding coefficient should carry the coefficient it used and the fact that its diversity was
  divided by `(1 − F)`, so a reader can size the uncertainty themselves. That belongs with the
  pre-pass's own output rather than here.

## What this does not measure

The pre-pass's cohort gather does not exist in the code, so route 2's upstream half — how much a
wrong `θ` moves the *fitted* spectrum through the regularizer — is estimated from the ratios spec
§4.1 records rather than measured. It can be measured the day the gather lands.
