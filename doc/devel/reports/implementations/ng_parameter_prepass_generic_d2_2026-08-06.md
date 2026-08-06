# ng step 4, D2 — the generic `NoiseModel`, §5.1's closed form

**Date:** 2026-08-06. **Branch:** `ng-parameter-estimation`. **Plan:**
`doc/devel/ng/impl_plan/parameter_prepass_generic.md`, Milestone D step D2 — *own commit, do not
bundle*, the fourth of the six steps whose failure is silent. **Design:**
`doc/devel/ng/arch/parameter_prepass_generic.md` §4.2 and §5.1,
`doc/devel/ng/spec/parameter_prepass_generic.md` §1, §2 and §12.8,
`doc/devel/ng/spec/parameter_prepass.md` §3.

## What landed

**`src/ng/parameter_estimation/fitting/mod.rs`** — the `NoiseModel` trait, the one genuine swappable
seam in step 4. Two associated types (`Cell`, `NoiseParams`) and one method.

**`src/ng/parameter_estimation/generic/noise_model.rs`** — a new file, and the SNP/indel path's
implementation:

- `LibraryNoise` — one library's share of the sample's reads and its error rate, named in a struct
  because a share and a rate are both small fractions and transposing them is a plausible wrong
  answer rather than a compile error.
- `SampleLibraryNoise` — every library of one sample, ascending by read group, shares summing to
  one. This is the model's `NoiseParams`: what the profile scan steps through. `single()` is the
  whole of the read-group table, where each entry belongs to one group already.
- `SubstitutionNoiseModel` — stateless, because everything travels in the parameters.
- The scoring rule itself, in two arms.

**`examples/ng_multilib_key_harness.rs`** — one print-only section, reachable only by naming it
(`--only=oracle`). It generates the fixture for the fourth oracle below. It computes nothing the
harness measures and cannot run on any other invocation; the harness's own output was re-checked
unchanged afterwards.

## What the rule is

Each read independently picks a library — library `g` with probability `w_g` — and then shows the
alternative allele with probability `p_j(ε_g)` or the reference with `1 − p_j(ε_g)`, where

```text
p_j(ε)  =  (j/P)·(1 − ε/3)  +  (1 − j/P)·ε
```

So a cell of the **attributed** arm — which records how many alternative reads came from each
library and how many reads showed the reference in total, having forgotten how the *depth* split —
is one multinomial over `G + 1` categories:

```text
                       n!                                                   n−k
L(cell | j)  =  ────────────────  ·  Π (w_g·p_j(ε_g))^{k_g}  ·  ( Σ w_g·(1 − p_j(ε_g)) )
                 Π k_g! (n−k)!       g                           g
```

**The sum in the last factor runs over every library the sample has**, listed at this cell or not: a
library that showed no alternative read here still produced reads that showed the reference. The
product before it runs only over the listed ones, which is the same thing, since a `k_g` of zero
contributes a factor of one. A **pooled** cell is the same expression with the `G` alternative
categories collapsed into one, leaving a binomial at the share-weighted rate `Σ_g w_g·p_j(ε_g)`.

Nothing is approximated in either. What is *not* done is inventing the forgotten split — giving each
library `n̂_g = w_g·n` — which is not a probability over the cell space and reports heterozygosity
68% high and the homozygous-non-reference rate 78% low at three reads a site, on two libraries with
the *same* error rate, without shrinking as data accumulates (research note §2).

## The four oracles, and that each one bites

The plan asks for spec §12.8's three identities plus agreement with the research harness. All four
are unit tests, and none of them needs a fit.

1. **The rule sums to one over the cell space, at any parameters.** Every cell a site of a given
   depth can land in — keyed exactly as the accumulator keys it, attributed at or below four
   alternative reads and pooled above — summed under one genotype, comes to one within 1e-12. Run
   at depths 1, 3, 5 and 9 on two libraries at ploidy 2 and 4; at depths 2 and 6 on **four**
   libraries, where four alternative reads have 35 splits rather than five; and at an error rate of
   0.3, far above anything the ladder reaches, because the identity is algebraic rather than a
   numerical accident near plausible rates.
2. **No cell is charged a negative count of reference reads.** Asserted inside the model, and
   exercised over the deep alternative-rich corner — every depth from 100 to 124 crossed with every
   alternative count from 0 to that depth, 2,900 cells — which is the exact case the per-bin mean
   fails. A cell whose depth sits below its own alternative count is named rather than scored.
3. **With every library's error rate equal, the attributed rule reproduces the pooled one.** In two
   halves, because the claim is exactly true only in one of them. At **one** library the two
   expressions are the *same expression* — `w = 1` collapses the multinomial coefficient back to a
   binomial one — and that is asserted directly, to 1e-12. At two libraries with equal rates the
   attributed rule still keeps a per-cell factor the pooled one does not, the probability of *that*
   split of the alternative reads; it carries no genotype, so what must agree exactly is the
   **differences between genotypes**, which is what a fit is computed from.
4. **Agreement with `examples/ng_multilib_key_harness.rs`'s `ln_component_attributed`.** Eighteen
   cells of the harness's `ratio=4 depth=6 split=skew90` world — two libraries at `ε` = 0.001 and
   0.004 with a 90/10 split — generated by the new `--only=oracle` section and pasted in as a
   fixture with its provenance. **Compared as differences between genotypes, and that is the
   sharper comparison rather than the looser one.** The two implementations reach `ln Γ`
   differently — the harness by a Lanczos series, this file through `libm` via
   `crate::genetics::lgamma` — so their factorial prefactors disagree in the last bits. Those
   prefactors carry no genotype, so they cancel from every difference, and two rules that agree on
   the differences give the same responsibilities, the same climb and the same answer. The absolute
   value of the prefactors is held instead by oracle 1, which cannot pass if they are wrong.

**Six mutations, five killed.** Because "the code is right" is not what these tests are for:

| mutation | killed by |
|---|---|
| the retired plug-in: `Σ_g w_g·ln(1−p_g)` in place of `ln Σ_g w_g·(1−p_g)` — the mean of a log for the log of a mean | all three sum-to-one tests, the pooled-equals-its-splits test, and the harness oracle (5 tests) |
| `1 − ε` in place of `1 − ε/3` on the alternative copy | the spec's three rows, and the harness oracle (2) |
| the attributed arm forgets the library share `w_g` | all three sum-to-one tests and the splits test (4) |
| the reference-read sum runs over the **listed** libraries only | all three sum-to-one tests, the splits test, and the harness oracle (5) |
| the negative-reference-read guard deleted | its refusal test (1) |
| the reference rate written `1 − Σ_g w_g·p_g` instead of `Σ_g w_g·(1 − p_g)` | **nothing — and correctly so** |

That last one is not a gap. The two forms are equal because the shares sum to one:
`Σ_g w_g·(1 − p_g) = Σ_g w_g − Σ_g w_g·p_g = 1 − Σ_g w_g·p_g`. The comment beside the function
originally claimed they were *not* complements, which was wrong; it now states the identity and why
the longhand form is kept anyway — it is the form the spec writes, and the shares' sum is an
invariant of `SampleLibraryNoise` rather than of this expression.

Note what the sum-to-one identity does **not** catch: the `1 − ε/3` mutation, which is still a valid
probability. That is what the harness oracle is for, and the two together are why the plan asks for
four checks rather than three.

## Deviations from the architecture, recorded

1. **`genotype_likelihoods` is `append_genotype_likelihoods`, and it appends.** The architecture
   writes `out: &mut Vec<f64>` without saying whether the callee clears. Appending is what the
   profile scan needs — it clears one flat buffer per rung and calls this once per cell, and what
   comes out is exactly the row-major table `GenotypeLikelihoodTable` borrows, with no per-cell row
   and no copy. The name says so because a model that cleared instead would silently leave the scan
   holding one cell's row. Pinned by `the_model_appends_rather_than_clearing`.
2. **`NoiseParams` carries the libraries' shares as well as their rates.** The architecture calls
   the associated type "the noise parameters being scanned", and a share is data rather than a
   parameter. The alternative is to put the shares on the model and the rates in the parameters —
   two collections indexed by position, which is how a library's rate gets scored against another
   library's share. They are inseparable in the expression, so they travel together.
3. **A fifth file under `generic/`.** `noise_model.rs`, following A4's `depth_bins.rs` precedent.
   `generic/mod.rs` already does four jobs (a recorded Milestone-C deferral) and this is ~330 lines
   of implementation plus its tests.
4. **`ploidy` stays a method argument even though `Cell` carries one**, and the two are asserted
   equal. The architecture's §4.2 says ploidy travels with the cell and its trait signature passes
   it separately; keeping the argument keeps the trait path-independent for the STR path, whose cell
   may not carry a ploidy, and the assertion turns the redundancy into a checked invariant instead
   of a place the two can disagree.
5. **One print-only section added to a research harness.** `--only=oracle`, unreachable without the
   flag, so no measurement moves. It exists because a unit test cannot call into an example binary
   and the alternative — transcribing the harness's expression into the test — would make the
   fourth oracle a comparison of the author with the author.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → the library binary at **3,080
passed**, 0 failed, 5 ignored, up from 3,060; `cargo doc --no-deps --lib` at the 12-unresolved-link
pre-existing baseline. `ng::parameter_estimation` 158 → **178** tests, of which `noise_model` is 20.

`examples/ng_multilib_key_harness.rs` was run to completion both before the step and after the edit
to it, with its measured output unchanged.

## Open, carried forward

- **Nothing here reads a locus, and nothing here is fitted.** The scan that consumes this is D3;
  `SampleLibraryNoise` has no producer until E1 computes each library's share from the read-group
  table's read counts.
- **The `½` a heterozygote's reads are assumed to split at** is the model's one soft assumption
  (`spec/parameter_prepass_generic.md` §8, open). Adopting a reference-bias term later lengthens the
  scan's ladder and changes no signature here, because `NoiseParams` is already an associated type.
