# ng calling loop — B1: the E-step for one sample

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step B1
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §2, §7, §8;
[arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §1, §2;
[spec/calling_priors.md](../../ng/spec/calling_priors.md) §3, §6
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Five category reviews raised three
> Blockers and seven Majors, and the fixes changed what §3, §4 and §5 describe: the
> total-weight check is now release-held and is the function's only `NaN` detector, a
> release-held check on the sample's own expected copies was added, the copy-table check went
> the other way to `debug_assert`, all three new items are `pub(crate)`, and the six tests
> became twelve. **One finding is a real bug this report's §3 got wrong**: the finiteness check
> it describes cannot see a `NaN`, because a `NaN` never wins the comparison that picks the
> maximum. [The review](../reviews/ng_calling_loop_b1_2026-08-25.md), and
> [what was done about it](../reviews/fixes_applied_2026-08-25_v3.md).

---

## 1. Plan

The first arithmetic in the calling loop: **one sample's share of the E-step**. Given that
sample's filled genotype-likelihood row, work out how probable each candidate genotype is and
turn the answer into expected allele copies.

Four steps, each the next one's input (spec §2's pseudocode, the `for each sample s` block):

1. the sample's **concentration** — the locus's seed plus what the *other* samples showed,
   its own expected copies subtracted from the cohort total;
2. its **log-prior** over candidate genotypes, through the `GenotypePriorModel` seam;
3. its **posterior** — likelihood plus prior, normalised;
4. its **expected allele copies**, folded out of that posterior.

Pure, scratch-backed, no allocation. Nothing calls it yet: B2 adds the M-step beside it and
D1 assembles the loop around both.

## 2. Assumptions and deviations

Two, and the first is the substantial one.

### 2.1 `SampleScoringBuffers` — a borrow bundle the architecture did not need

Architecture §2 sketches `CallingScratch` with **public fields**, so a loop body would write
`scratch.posterior_row` and `scratch.prior_row` side by side and the borrow checker would
split the struct for it. **A1's review made every field private behind a per-buffer
accessor**, and each accessor borrows the whole scratch — so the four buffers step 1 needs
live at once cannot be reached one at a time. The three ways out are: make the fields public
again (undoing a review finding), copy a buffer to break the borrow (an allocation, in the
one function whose whole shape exists to have none), or hand the disjoint borrows out
together. This step takes the third.

`SampleScoringBuffers<'a>` is eight borrows of eight different fields, made by
`CallingScratch::sample_scoring_buffers(sample)`, which is also the one place the two flat
sample-major tables are sliced. Nothing about the buffers themselves changed.

**This is the third divergence from architecture §2's `CallingScratch` sketch**, after A1's
two (private fields with three per-allele buffers where the sketch has one
`concentration: Vec<f64>`, and `SampleGenotypeCall` as an enum). Whether the architecture is
amended is the owner's call; the divergences are recorded, not applied there.

### 2.2 The prior arrives as `&dyn GenotypePriorModel`

Not a generic parameter. The seam exists to compare two priors that differ by 11 points of
genotype accuracy on GIAB at 5× (`spec/calling_priors.md` §2.2), and a run selects between
them at run time, so the trait object is what a caller will already be holding. The prior's
own module pins object safety for the same reason.

## 3. Changes made

- **[src/ng/calling/inference/summarise_condition.rs](../../../../src/ng/calling/inference/summarise_condition.rs)**
  — new. The module doc names arm A and defines the E-step and the M-step in the reader's
  language before either does work in a sentence. `score_one_sample` is the whole of the
  step's arithmetic.
- **[src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs)** — `SampleScoringBuffers<'a>`
  and `CallingScratch::sample_scoring_buffers`.
- **[src/ng/calling/inference/mod.rs](../../../../src/ng/calling/inference/mod.rs)** — one
  `pub mod` line.

**Six checks are held in release and one in debug**, per spec §8's rule that caller bugs are
assertions — *amended after the review, which found two of the original four untested and one
of them unable to do what its doc claimed.* Three of the six guard a *truncation*: `zip` and
`chunks_exact` stop at the shorter side, so a likelihood row or a posterior row of the wrong
length leaves the genotypes past its end holding the previous sample's numbers, and the call
comes out confident and wrong. The seed check buys the message rather than the catch. The
remaining two are the review's: the total weight, which is the function's **only** `NaN`
detector, and the finiteness of the sample's own expected copies, without which the scratch's
`NaN` sentinel is absorbed as a zero cohort term. The debug-held one is the copy table's
width, which cannot fire while every view comes from `GenotypeTable::build`.

**Each of the six is reached by a test that fails under `--release` when it is downgraded** —
measured in one run with all six downgraded together: `5 passed; 6 failed`, a one-to-one
mapping. Before the review, two of the four were reached by nothing in either profile.

## 4. Tests added

Six, in the file's own `mod tests`.

| test | what it pins |
|---|---|
| `one_samples_e_step_matches_the_arithmetic_done_by_hand` | every intermediate of one pass on numbers a reader can follow: concentration `[3, 0.5]`, posterior `1/7, 2/7, 4/7`, copies `4/7` and `10/7` |
| `the_posterior_is_unchanged_by_a_constant_that_would_underflow_the_exponentials` | the max-subtraction: 1,000 nats added to every genotype leave the posterior where it was |
| `at_one_sample_the_concentration_comes_back_as_the_seed` | spec §7 — the leave-one-out term is exactly zero at one sample, by arithmetic and not by a branch |
| `every_samples_posterior_is_a_distribution_and_its_copies_sum_to_the_ploidy` | the two locus-wide invariants, on a triallelic locus where a mis-shaped fold still writes plausible numbers |
| `a_short_likelihood_row_is_refused` | the release-held truncation check, reached through a hand-built bundle because the scratch cannot produce a short row |
| `a_non_finite_score_is_refused` | a `NaN` is raised where it arrives, rather than reaching the cohort's copies and every other sample's next prior |

**Three mutations were run against them, and each is quoted in the doc comment of the test
that catches it**, so the claim cannot go stale silently:

| mutation | what the suite prints |
|---|---|
| drop `+ prior_of_genotype.get()` | `posterior 0.25 against 0.14285714285714285` |
| drop `- largest_score` before `.exp()` | debug: `the normaliser cannot come out below one: got 0 over 3 genotypes`; `--release`: `posterior NaN against 0.25` |
| fold with `chunks_exact(genotype_count)` | `copies 0 against 1.4285714285714286`, and `expected copies sum to 1.9882062059027987` |

The third is the one worth noticing: folding over the wrong axis lands **0.0118 away from the
ploidy**, close enough that only an exact check finds it — which is why the tolerance there is
`1e-12` and not a loose one.

## 5. Validation

All in the container, from this worktree.

- `cargo fmt --all -- --check` — exit 0, no output.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4540 passed; 0 failed; 14 ignored`. The pre-B1 baseline was **4,528**,
  so B1 adds twelve: the six above, five the review added, and one on the scratch accessor.
- `cargo test --release --lib ng::calling::inference::summarise_condition` — `11 passed; 0
  failed`, which is the run that can tell `assert!` from `debug_assert!`.

**Two aggregate gates are red and neither is this work's**, both verified rather than assumed:

- `cargo test --all-targets --all-features` aborts in `benches/psp_writer_perf.rs:386`
  (`index out of bounds: the len is 3300000 but the index is 3300000`) — the standing item
  `PROJECT_STATUS.md` already records as failing on an unmodified `HEAD`.
- `cargo test --lib --bins --tests --examples --all-features` fails `--example
  ng_generic_loci_dump` with `2 passed; 11 failed`, every failure at that file's line 956 on
  `ReferenceCheck::VerifyAgainstIndex`. **Checked in a detached worktree at `2b6aceae`
  without this patch: the identical `2 passed; 11 failed`.**

## 6. Trade-offs and follow-ups

- **Nothing here is exercised above 3 samples or 3 alleles.** `CLAUDE.md`'s range commitment
  runs to several thousand samples and to loci with more alleles than that; the E-step has no
  cohort-size branch to break, but nothing yet demonstrates it. The natural home is D1, where
  a loop exists to run.
- **The flat prior-free first pass (C1) does not go through this function**, because it uses
  no prior at all and must not build a concentration. **The review measured what happens if it
  tries**: handing this function a flat `GenotypePriorModel` still runs step 1, which reads
  cohort copies nothing has written — a panic in debug, and in release a concentration that
  comes back as the seed exactly, which is the seed-only first pass spec §3 spends four
  paragraphs ruling out. A repair exists and is green in a reviewer's worktree
  (`PassPrior::{Flat, LeaveOneOut}`, one entry point and one `match`); it is **not applied**,
  because the shape of C1's entry point is C1's decision.
- **`score_one_sample` is `pub(crate)` with no caller yet**, and the dead-code lint is
  `expect`ed under `cfg(not(test))` so that D1 adding a real caller turns the line into a
  compile error rather than leaving a stale suppression behind.
