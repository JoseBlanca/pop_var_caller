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
   at depths 1, 3, 5 and 9 on two libraries at ploidy 2 and 4; at the same depths on a **single**
   library at ploidy 1, 2 and 4, which is the configuration every entry of the read-group table and
   1,550 of the 1,707 archive samples are scored in; at depths 2 and 6 on **four**
   libraries, where four alternative reads have 35 splits rather than five; and at an error rate of
   0.3, far above anything the ladder reaches, because the identity is algebraic rather than a
   numerical accident near plausible rates.
2. **No cell is charged a negative count of reference reads.** Asserted inside the model, and
   exercised over the deep alternative-rich corner — every depth from 100 to 124 crossed with every
   alternative count from 0 to that depth, **2,825** cells. What that sweep asserts is that every
   score is a log-probability (`ln L ≤ 0`) at every genotype; **the guard itself is held by its own
   one-line refusal test**, and deleting the guard fails only that one. Two further tests bracket
   `REFERENCE_READ_TOLERANCE` from both sides — a cell a millionth of a read short is refused, one a
   trillionth short is clamped — so widening the constant to swallow the per-bin-mean bug cannot
   pass.
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
   fixture with its provenance, and verified bit-identical against a fresh regeneration.
   **Both the absolute values and the genotype differences are asserted**, and the first is the one
   that took measuring: the two implementations reach `ln Γ` differently — the harness by a Lanczos
   series, this file through `libm` via `crate::genetics::lgamma` — so it is not obvious the
   absolute values can be compared at all. Over these 54 numbers they agree to **1.4 × 10⁻¹⁴**,
   under six units in the last place and seventy times inside the test's tolerance.

   **An earlier version asserted only the differences, and the reasoning was backwards.** It ran:
   the factorial prefactors carry no genotype, so they cancel, so the difference is the sharper
   comparison. The premise is true and the conclusion is not — cancelling a term makes the check
   **weaker**, and it costs a real kill, because a rule that dropped the library share `w_g` adds
   `−Σ_g k_g·ln w_g`, which also carries no genotype and slips straight through. Found by review,
   measured, and now both are asserted.

**Thirty-three mutations, thirty-two killed.** Six before review and twenty-seven during it;
"the code is right" is not what these tests are for. A representative selection:

| mutation | killed by |
|---|---|
| the retired plug-in: `Σ_g w_g·ln(1−p_g)` in place of `ln Σ_g w_g·(1−p_g)` — the mean of a log for the log of a mean, **applied to the attributed arm**, which is where the plug-in was | all three sum-to-one tests, the pooled-equals-its-splits test, and the harness oracle (5 tests) |
| `1 − ε` in place of `1 − ε/3` on the alternative copy | the spec's three rows, and the harness oracle (2) |
| the attributed arm forgets the library share `w_g` | all three sum-to-one tests and the splits test (4) |
| the reference-read sum runs over the **listed** libraries only | all three sum-to-one tests, the splits test, and the harness oracle (5) |
| the negative-reference-read guard deleted | its refusal test (1) |
| `REFERENCE_READ_TOLERANCE` widened from 1e-9 to 1.0 | the millionth-short refusal test (1) |
| `single()` builds its library under read group 0 whatever it was handed | the single-library constructor test (1) |
| the genotype loop runs `0..ploidy` rather than `0..=ploidy` | 11 |
| `count_times_ln`'s two arguments transposed | 11 |
| the two libraries' **shares swapped between their read groups** in the fixture | the harness oracle, **and nothing else** |
| the two libraries' **rates swapped between their read groups** in the fixture | the harness oracle, **and nothing else** |
| the reference rate written `1 − Σ_g w_g·p_g` instead of `Σ_g w_g·(1 − p_g)` | **nothing — and correctly so** |

**The last two kills are the ones to carry forward.** Pairing a share to the read group that
produced it has *no identity behind it*: a rule with the shares swapped is still a probability over
the cell space, so none of the three sum-to-one tests can see it, and only the harness fixture —
which happens to live in a 90/10 world — catches it. That is why `LibraryNoise` holds the read
group, the share and the rate in one struct rather than in parallel collections, and it is a warning
for **E1**, which computes those shares from read counts and will have no identity oracle behind the
pairing at all.

That last one is not a gap. The two forms are equal because the shares sum to one:
`Σ_g w_g·(1 − p_g) = Σ_g w_g − Σ_g w_g·p_g = 1 − Σ_g w_g·p_g`. The comment beside the function
originally claimed they were *not* complements, which was wrong; it now states the identity and why
the longhand form is kept anyway — it is the form the spec writes, and the shares' sum is an
invariant of `SampleLibraryNoise` rather than of this expression.

Note what the sum-to-one identity does **not** catch: the `1 − ε/3` mutation, which is still a valid
probability. That is what the harness oracle is for, and the two together are why the plan asks for
four checks rather than three.

## Deviations from the architecture, recorded

0. **The trait declares its own genotype count, `genotypes(&self, ploidy) -> usize`.** Not in the
   architecture at all, and added on review. The append contract was originally written as
   "`ploidy + 1` entries", which is the SNP/indel path's number and not the seam's: on the STR path
   a genotype is an unordered tuple of allele lengths, so a diploid stratum spanning nine lengths
   has `A(A+1)/2` = **45** (`spec/parameter_prepass_ssr.md` §4.2). D3 hands that width to
   `GenotypeLikelihoodTable::from_natural_logs` as the table's only shape argument, and 45 columns
   read as 3 still divides — so a mis-shaped table would be accepted and the climb would run on
   transposed rows. It also gives back the one thing appending costs against clearing: under a
   clearing contract the row width is `out.len()` for free, under an appending one it is a
   difference nobody records.
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
5. **Not a deviation but recorded here: one print-only section added to a research harness.** `--only=oracle`, unreachable without the
   flag, so no measurement moves. It exists because a unit test cannot call into an example binary
   and the alternative — transcribing the harness's expression into the test — would make the
   fourth oracle a comparison of the author with the author.

## Validation

All in the container. `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -D
warnings` clean; `cargo test --lib --bins --tests --all-features` → the library binary at **3,089
passed**, 0 failed, 5 ignored, up from 3,060, with the nine other binaries unchanged at 69 passing;
`cargo test --doc ng::parameter_estimation` → 1 passed; `cargo doc --no-deps --lib` at the
12-unresolved-link pre-existing baseline, none of the twelve in `parameter_estimation`.
`ng::parameter_estimation` 158 → **187** tests, of which `noise_model` is 29.

`examples/ng_multilib_key_harness.rs` was run to completion both before the step and after the edit
to it, with its measured output unchanged.

## Open, carried forward

- **Nothing here reads a locus, and nothing here is fitted.** The scan that consumes this is D3;
  `SampleLibraryNoise` has no producer until E1 computes each library's share from the read-group
  table's read counts.
- **The `½` a heterozygote's reads are assumed to split at** is the model's one soft assumption
  (`spec/parameter_prepass_generic.md` §8, open). Adopting a reference-bias term later lengthens the
  scan's ladder and changes no signature here, because `NoiseParams` is already an associated type.
