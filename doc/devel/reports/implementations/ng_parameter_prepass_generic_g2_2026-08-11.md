# G2 — the coverage sweep, and the bias it makes visible as a curve

**Step:** G2 of `impl_plan/parameter_prepass_generic.md`, Milestone G.
**Date:** 2026-08-11. **Runner:** seven arms of G1's anchor. There is no separate program — the loop is in
`generic/truth_anchors.rs`'s own doc comment, beside the single-arm invocation it varies.

## What it is, and why it needed no new estimator code

The same genome, the same confident regions, seven depths. An error rate is per read and the
other three parameters are properties of a genome, so **none of them has any business depending
on how deeply the sample was sequenced**. Any slope is bias, and its sign names the mechanism.

G1's anchor already fits the sample, counts the benchmark over the loci walked, and checks the
two against each other. A sweep is that measurement at seven depths — so **every rung is bounded
as well as plotted**, rather than being a row in a table nothing checks.

**Seven rungs, where the plan asked for four.** It specifies 300×, 30×, 10× and 3×; the
repository holds **5×, 10×, 15×, 20×, 30×, 50× and 300×**. The plan's 3× does not exist and was
not created: 5× is the shallowest real rung, and what this measures is a slope, which seven
points describe better than four.

## The table

| depth | fitted `ε` | model-free `ε` | apart | heterozygosity / truth | every-copy-non-ref / truth |
|---|---|---|---|---|---|
| 5× | 2.3083 × 10⁻³ | 2.2940 × 10⁻³ | +0.6% | **0.960** | **0.994** |
| 10× | 2.3441 × 10⁻³ | 2.3538 × 10⁻³ | −0.4% | 1.055 | 0.939 |
| 15× | 2.3638 × 10⁻³ | 2.3504 × 10⁻³ | +0.6% | 1.038 | 0.942 |
| 20× | 2.3349 × 10⁻³ | 2.3553 × 10⁻³ | −0.9% | 1.071 | 0.939 |
| 30× | 2.3327 × 10⁻³ | 2.3300 × 10⁻³ | +0.1% | 1.085 | 0.937 |
| 50× | 2.3489 × 10⁻³ | 2.3452 × 10⁻³ | +0.2% | 1.111 | 0.949 |
| 300× | 2.5039 × 10⁻³ | 2.4963 × 10⁻³ | +0.3% | **1.142** | 0.947 |

All seven pass the anchor's check.

## The error rate is flat, and it tracks a number owing nothing to the model

Seven values scattered from −0.9% to +0.6% of the model-free count, with no trend in depth. The
count itself moves little across the range — 2.294 to 2.355 × 10⁻³ from 5× to 50×, with 300× at
2.496 × 10⁻³ — and the fit follows it rung for rung. This is the parameter step 4a changed most,
and across a sixty-fold range in depth it stays on a count with no model in it.

## Heterozygosity is not flat, and the shape rules out the obvious culprit

It rises from **0.960** of the benchmark at 5× to **1.142** at 300×.

**It is not monotonic** — 10× sits at 1.055 and 15× at 1.038, out of order by 1.7%. The trend
across the range is unmistakable and the local ordering is not; two adjacent rungs differing by
less than two percent are not evidence of anything on their own.

**It crosses the truth between 5× and 10×.** Only the shallowest rung undercounts heterozygotes;
every rung from 10× up invents them.

**The depth cap is not the cause, and that is what the sweep settles.** The cap draws a locus
down to 124 reads and bites at 300× and nowhere else — 545,863 of that arm's 550,049 sites
against zero at 30×. A mechanism that fires on one rung cannot bend the five below it, and the
five below it are already bent. The slope is present from 10× to 50×, where no site is capped at
all.

**What it looks like instead** is the mechanism the research note described from the other end: a
shallow site cannot tell an error from a variant, and the direction of the confusion is not
symmetric. Below the crossing point the fit loses real heterozygotes it has too few reads to
confirm; above it, the remaining tail of bad positions is read as variation. Step 4a removed most
of the second half — the same span was +9.7% between 30× and 300× before the second class of site
and is +6.0% after — and the residue is what this curve shows.

## A third thing, which nobody was looking for

**The every-copy-non-reference shortfall is absent at 5×.** It sits at 0.994 of the benchmark
there, and at 0.937 to 0.949 at every deeper rung — appearing between 5× and 10× and flat
thereafter. That rate has been 5 to 6% low since this milestone began and its cause is
unexplained; the sweep says whatever causes it **needs depth to express itself**, which is the
opposite of what a missing-evidence explanation would predict.

**One point is not a finding.** At 5× the fit has three reads a site and the least information of
any rung, so 0.994 could be luck. What is worth recording is that the shortfall is not constant
across depth, because that is a testable handle on a question that has had none.

## What this step corrected in G1

**G1's one assertion was wrong, and the sweep is what exposed it.** It demanded the fitted error
rate be no lower than the counted one, and two rungs — 10× and 20× — failed it by 0.4% and 0.9%.

The estimator was not at fault. `arch` §9 argues the count is a floor because the confident
regions are the *easy* ones: a fit that saw the whole genome should land above a count taken only
from the easy parts. **G1 removed that premise without noticing** — it counts and fits over
exactly the same loci, all inside the confident regions, so there is no easy-against-hard
asymmetry left. And the fitted rate is a **rung**: the ladder steps by a quarter of a Phred, 5.9%,
so asking a quantised number to stay reliably on one side of a continuous one is asking for a
resolution it does not have.

The assertion is now **agreement to within half a rung**, expressed in the ladder's own units.
That is three times the worst deviation in the sweep, and it still catches what the inequality
caught: the wrong sample's benchmark misses by 1.1 rungs, reading `POS` as 0-based by about five,
a classifier that returns nothing by about six.

## Validation

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --lib --bins --tests --all-features` and `cargo doc --no-deps --lib` at the
12-unresolved-link pre-existing baseline, none in this module. All seven arms of the sweep green.
