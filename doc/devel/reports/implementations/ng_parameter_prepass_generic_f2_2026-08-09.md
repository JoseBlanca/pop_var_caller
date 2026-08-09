# ng step 4, the SNP/indel path — F2: recovery from a directly-filled table

**Date:** 2026-08-09. **Plan:** F2. **Design:** arch §9, spec `parameter_prepass.md` §10.1.

Tests only — no production behaviour changes, and one **defect fixed in the shared fixture
generator** that F2 is the first step to reach.

## What was built

`generic/recovery.rs`, a `#[cfg(test)]` module: a table filled cell by cell from a known truth,
handed to `fit_coupled_from_tables`, and asked to find its way back. Six tests.

**The arms, and why each is a condition the others hide:**

| ploidy | depth | what only this arm can see |
|---|---|---|
| 2 | 3 | tomato's regime — the one every fit before F2 was proven in |
| 2 | 124 | the widest geometric bin, `98..=124`, where the bin index and row offset both have to be right |
| 4 | 4 | dosages 1, 2 and 3 — an answer returned off by one is a wrong number, not a compile error |
| 4 | 124 | both at once, the only arm where a dosage mix-up and a binning fault could cancel |

Plus two tests that pin the arms' own premises: that 3 reads and 124 really do sit on opposite
sides of the binning rule, and that tetraploid-at-3-reads is **not identified** (below).

## The defect this found, and it is in the file every fit's tests share

**Three of the four arms failed on first run**, with the top dosage coming back at **0.0000**
against a generating 0.020 — in *both* ploidies at depth 124.

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

**A ploidy-`P` sample needs depth ≥ `P` before its genotype frequencies are identified at all.**
A table at depth `d` has `d + 1` cells whose probabilities sum to one, so it carries `d`
independent numbers against `P` free frequencies. At `d = 3, P = 4` that is three equations for
four unknowns: the likelihood has a ridge, and the climb lands somewhere on it — dosage 1 comes
back at 0.158 against a generating 0.150, 5.6% away, **with the error rate on the right rung
throughout**. Pinned by `a_tetraploid_sample_at_three_reads_is_not_identified`, because the
shallow-and-polyploid corner is exactly where a later reader would blame the fit.

Tomato satisfies the rule at three reads a site only because tomato is diploid.

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
