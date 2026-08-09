# ng step 4, the SNP/indel path — F2: recovery from a directly-filled table

**Date:** 2026-08-09. **Plan:** F2. **Design:** arch §9, spec `parameter_prepass.md` §10.1.

Tests only — no production behaviour changes, and one **defect fixed in the shared fixture
generator** that F2 is the first step to reach.

## What was built

`generic/recovery.rs`, a `#[cfg(test)]` module: a table filled cell by cell from a known truth,
handed to `fit_coupled_from_tables`, and asked to find its way back. Seven tests.

**The arms, and why each is a condition the others hide:**

| ploidy | depth | what only this arm can see |
|---|---|---|
| 2 | 3 | tomato's regime — the one every fit before F2 was proven in |
| 2 | 124 | the widest geometric bin, `98..=124`, where the bin index and row offset both have to be right |
| 4 | 8 | dosages 1, 2 and 3 — an answer returned off by one is a wrong number, not a compile error |
| 4 | 124 | both at once, the only arm where a dosage mix-up and a binning fault could cancel |

Plus three tests that pin the arms' own premises: that the shallow and deep depths really do
sit on opposite sides of the binning rule, and that a tetraploid at three **and at four** reads
is not identified.

## The defect this found, and it is in the file every fit's tests share

**Two arms failed on the underflow**, with the top dosage coming back at **0.0000** against a
generating 0.020 — in *both* ploidies at depth 124. (Three failed on the first run overall; the
third was the tetraploid shallow arm, and its cause was identifiability rather than underflow.
The commit message for `7a207bac` says "three of the four failed" without separating the two
causes, which a reviewer flagged and this corrects.)

`table_generated_at` built its binomial term up from `(1 − p)^depth` by repeated
multiplication, on the stated grounds that this "keeps every intermediate inside `f64` at the
depths used here". It did, until F2 raised the depth. For a homozygous non-reference genotype
`p` is `1 − ε/3`, so at 124 reads the starting term is `0.00033^124` ≈ **10⁻⁴³¹**, which
underflows to exactly zero; every later multiplication keeps it there. The genotype then
contributes nothing to any cell, **the fixture silently omits its homozygous non-reference
sites, and the fit correctly reports 0.0000 for a class the table never contained.**

That is precisely the failure `expected_counts.rs` exists to prevent — "two copies of a fixture
generator are two chances for one of them to drift" — in the file shared by every fit's tests.
Fixed by computing the term in logs (`ln` of the same quantity is −992, an ordinary number).
**Every existing user of the generator was re-run and none moved**, including E2's 25 worlds;
the previous uses top out at depth 40, below where the cliff is.

## The measured finding worth carrying forward

**The counting argument gets the boundary wrong on its own, and the review is what exposed
that.** A depth-`d` table has `d + 1` cells whose probabilities sum to one, so it carries `d`
independent numbers against a ploidy-`P` truth's `P` free frequencies — hence `d ≥ P`. **That is
the condition for the frequencies given the error rate**, and this fit must find the rate too.
Measured at `P = 4`: at `d = 4` the rate never leaves the rung the fit started from and the
frequencies come back **7.2%** away; at `d = 3`, 17.9% away; both are identified by `d = 8`.

**An earlier version of this file asserted that four reads recovers a tetraploid, and it
passed** — because the tables were generated at Phred 30, which is `DEFAULT_ERROR_RATE`, which
is exactly where the coupled fit begins. The rate looked recovered by never moving, and the
frequencies were right because they were conditional on a rate that happened to be correct.
Generating ten rungs away turned both halves into claims with content. Two tests now pin three
and four reads as *not* identified, because shallow-and-polyploid is where a later reader would
blame the fit.

Tomato is unaffected: it is diploid, and three reads identifies a diploid.

## Tolerances, and why they differ by arm

Asserted per arm rather than once, because the arms differ by two orders of magnitude:

| arm | rung found | worst frequency error | asserted |
|---|---|---|---|
| 2 × 3 | one off | 0.334% | 0.5% |
| 2 × 124 | exact | 0.0009% | 0.01% |
| 4 × 8 | exact | 0.0019% | 0.01% |
| 4 × 124 | exact | 0.0024% | 0.01% |

**The shallow diploid arm's 0.334% is the coupling, not noise.** Its error rate lands one rung
from the generating one — inside the tolerance the design argues for, and about 6% in the rate —
and the frequencies are conditional on the rate, so they absorb it. That is the coupled fit
doing what it is for.

An earlier version asserted 1% everywhere, justified by the ladder's 0.3% binning bias (research
note §4.3). **That figure does not apply here**: it was measured on a mixed-depth world, and
every table in this file puts all its sites at one exact depth — which is exactly what
`expected_counts`' own doc says makes the binning rule contribute nothing. The looseness had a
measured cost: at 1%, the fixture underflow this milestone exists to have fixed passes the whole
suite.

## Recorded deviation: the plan's "300×"

F2 is specified at "3 reads a site and at 300×". **A directly-filled table cannot hold a site
above the ladder's cap** — `add_site` refuses it — and a 300× sample does not reach a table at
300 either: C2's cap subsamples it to 124 first, and 124 is the number the cell records. So
**124 is the deep arm**, and what happens to the reads above it is C2's own oracle, which proves
the kept alternative count hypergeometric in mean and variance. "300×" and "the deepest a cell
can be" are the same intent and different numbers, so this is recorded rather than absorbed.

The plan's stated *reason* for the deep arm survives intact: at three reads every site sits in a
one-per-depth bin, so a binning fault cannot show, and at 124 the row is 27 cells wide and
shared.

## What F2 does not establish

- **Recovery cannot catch a scoring rule that is wrong self-consistently.** On a table of
  expected counts the score is `N · Σ_c p_c(θ₀) · ln p_c(θ)`, which Gibbs' inequality maximises
  at `θ = θ₀` for *any* rule whose cell probabilities sum to one. That is D2's four identities'
  job and F2 adds nothing to it. What recovery does catch is a fit that mislabels, misgathers,
  misreports or misranks.
- **Nothing here reads a locus.** F2 fills the tables directly; F3 is where the two cohorts'
  alignments are walked.
- **The runs model is not exercised.** These arms fit error rates and genotype frequencies. `F`
  needs 3,000 windows and belongs to the harnesses.
