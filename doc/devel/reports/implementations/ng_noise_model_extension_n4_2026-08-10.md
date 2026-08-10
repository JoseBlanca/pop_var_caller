# N4 — the two measurement programs, extended to a second class of site

**Step:** N4 of `impl_plan/noise_model_extension.md`. Run **after** N5 on the owner's call: the
two steps share nothing, and the real-alignment numbers are the milestone's point.
**Date:** 2026-08-10.

## What these programs are, and what a result from them means

`examples/ng_multilib_key_harness.rs` and `examples/ng_inbreeding_harness.rs` build **simulated
samples** — a made-up sample with a known error rate per library, known heterozygosity and
known read depth — and ask the estimator to recover what went in.

The first one computes each cell's probability under the truth rather than drawing reads, so
what the fit converges to is what an infinite genome would give: **a departure is bias, with no
sampling noise in it, and a zero is consistency rather than luck.** The second cannot do that —
a hidden Markov model's likelihood does not factor over windows — so it draws genomes and reads
the error against the fraction of the drawn genome that really lay inside a run.

## The question only the first one can ask

Every piece of evidence behind this milestone came from **single-library** samples: both cohorts
carry one read group, and the research note's two synthetic cases have one error rate. The
second class of site is **one share and one rate for the whole sample**, while `ε` is **one rate
per read group** — so nothing anywhere had asked whether the pair survives a cell key that has
thrown away which library produced each alternative read.

Six simulated samples, at the shares and rates HG002 and tomato returned on real alignments,
each built once as one library and once as two, fitted through both keys:

| simulated sample | key | `ε₁` rungs | `ε_noisy` rungs | `w` | `π_het` | `π_hom_alt` |
|---|---|---|---|---|---|---|
| hg002-30x, one library | pooled | −0.000 | 0.000 | 0.00% | −0.00% | 0.00% |
| hg002-30x, one library | attributed | 0.000 | 0.000 | −0.00% | 0.00% | −0.00% |
| hg002-30x, two libraries | pooled | −0.000 | 0.000 | 0.00% | −0.00% | 0.00% |
| hg002-30x, two libraries | attributed | −0.000 | −0.000 | 0.00% | 0.00% | −0.00% |
| hg002-300x, one library | pooled | 0.000 | 0.000 | −0.00% | −0.00% | 0.00% |
| hg002-300x, one library | attributed | 0.000 | 0.000 | −0.00% | −0.00% | 0.00% |
| hg002-300x, two libraries | pooled | 0.000 | 0.000 | −0.00% | −0.00% | 0.00% |
| hg002-300x, two libraries | attributed | 0.000 | 0.000 | −0.00% | −0.00% | 0.00% |
| tomato, one library | pooled | −0.000 | −0.000 | 0.00% | 0.00% | −0.00% |
| tomato, one library | attributed | 0.000 | 0.000 | −0.00% | −0.00% | 0.00% |
| tomato, two libraries | pooled | −0.000 | −0.000 | 0.00% | 0.00% | −0.00% |
| tomato, two libraries | attributed | −0.000 | −0.000 | 0.00% | 0.00% | −0.00% |

**All five quantities, all twelve fits, exactly zero.** The estimator is consistent for the pair
under a key that keeps the library attribution and under one that does not.

**The control, and it must find nothing:** the same six with the second class taken out of the
truth. Every row returns the generating rate to −0.000 rungs, both frequencies to 0.000%, and
declines the second class.

## Two defects, both in the measuring program rather than the estimator

**The likelihood floor was in the wrong units, and it made every simulated sample decline its
own second class.** The caller refuses a second class that buys less than three units of
likelihood — a test on a real table, whose weights are counts of sites, χ²(2) at p ≈ 0.05 over
half a million of them. Copied into a program whose weights are **probabilities summing to
one**, three units is a threshold five orders of magnitude above any gain a real second class
buys. What this program measures is the argmax at infinite data, where such a floor has no
meaning; only arithmetic noise is worth rejecting.

**The first table was the search's resolution, not bias.** Twenty-five points spanning 10⁻³ to
0.3 step by a factor of 1.27 — **4.3 rungs** — so the noisy rate could only ever land within
about two rungs, and it did: 1.0 to 2.3 rungs, which reads as bias and is nothing but spacing.
A golden-section refinement between the winner's neighbours removes the grid from the answer.
**A residue of 0.015 rungs then survived on two of the twelve rows and was also mine**: three
times the inner passes and a convergence tolerance a thousand times finer took every row to
exactly zero.

## The second program: does a noisy class move the inbreeding coefficient?

It does not. Across 21 fits — seven genomes of 8,004 windows at a true inbreeding coefficient
of 0.30, shares up to 3% and both measured noisy rates — the fitted value differs from the
fraction of the drawn genome that really lay inside a run by **at most 1 part in 10,000**, and
told the pair it rounds to zero at four decimal places in all seven. The same estimator under
the evenly-spread floor of false heterozygotes the program already tested misses by up to
9 parts in 10,000.

**What the noisy class does move is the pair of heterozygote rates the model reads a run from,
and this measures a decision made earlier on reasoning alone.** N3c argued that the two rates
must reach the runs model as a pair rather than averaged into one number first. Measured, at a
share of 0.88%:

| what the fit was told | inside a run | outside a run | ratio |
|---|---|---|---|
| generating truth | 5.0 × 10⁻⁵ | 1.0 × 10⁻³ | 20× |
| the pair | 5.0 × 10⁻⁵ | 1.002 × 10⁻³ | 20× |
| the share-weighted marginal | 2.63 × 10⁻⁴ | 1.201 × 10⁻³ | 4.6× |
| the clean rate alone | 4.30 × 10⁻⁴ | 1.593 × 10⁻³ | 3.7× |

At a share of 3% the clean-rate-alone row collapses to 1.6×. **The averaged rate recovers about
half of what is lost and no more**, which is the argument N3c made, now with numbers.

The inbreeding coefficient survives all three because a noisy class is patchy site by site and
even window by window: one site in a hundred is noisy, but a 100 kb window holds 100,000 of
them, so a window's count varies by about 3% around a thousand and every window is lifted by
nearly the same amount. Both states rise together and the contrast the chain reads is
unchanged. The only trace is the chain's confidence — the share of windows it could not decide
goes from nothing to 5 in 1,000 at the most extreme setting tested.

## The regression the plan asks for

The 25 simulated samples the coupled fit was originally proved on all carry **one** error rate,
so the extended model has no second class to find in any of them. Re-fitted through it, from a
start at three times the true rates and half the true frequencies, every one returns **0.000
rungs on both libraries' rates and 0.000% on both genotype frequencies** — the plan's gate, and
the same answer the research note recorded before the model gained a second class.

That is a check on the extension, not evidence it works: a sample with nothing to find is a
sample where the only possible result is *unchanged*. It is stated here as what it is.

**A byte-for-byte diff of the whole report was not run for this program** — it was for the
other one. What stands in its place is structural: with no second class the mixture is never
formed at all, the truth builder evaluates the original expression and the scoring rule returns
the one-class components directly, so nothing else can move.

## What did not change, which is most of both programs

The second program's other nine sections are **byte-identical**: the full report before and
after differs by 50 lines, every one an addition. `None` takes a branch that calls the original
code path rather than a share of zero, so this is true by construction as well as by
measurement.

## Validation

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --lib --bins --tests --all-features` and `cargo doc --no-deps --lib` clean, the last
at the 12-unresolved-link pre-existing baseline with none in this module.
