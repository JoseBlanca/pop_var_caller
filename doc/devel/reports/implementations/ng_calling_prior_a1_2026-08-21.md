# ng genotype prior — A1: the `genotype_prior/` folder and the diversity scalar

*Implementation report, 2026-08-21. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step A1 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone A.*

> **This report describes the code as submitted for review.** The review changed part of it —
> most visibly, `DEFAULT_SPECIES_DIVERSITY_FALLBACK` became the associated const
> `ExpectedHeterozygosity::SPECIES_FALLBACK`, and four tests were added or strengthened. What
> landed is in [the review](../reviews/ng_calling_prior_a1_2026-08-21.md) and
> [the fixes](ng_calling_prior_a1_fixes_2026-08-21.md).

## 1. Plan

Two things, neither of them logic:

- **`src/ng/calling/genotype_prior/`** — the folder step 8 lives in, its `mod.rs` declaring the
  four files the plan's later milestones fill, wired into
  [`src/ng/calling/mod.rs`](../../../../src/ng/calling/mod.rs).
- **Two additions to ng's shared vocabulary**,
  [`src/ng/types.rs`](../../../../src/ng/types.rs): `ExpectedHeterozygosity`, the cohort's
  expected heterozygosity at ordinary sites — the spec's `θ` — constrained to `[0, 1]` with the
  file's `try_new`/`get` shape; and `DEFAULT_SPECIES_DIVERSITY_FALLBACK = 1e-3`, the value a run
  with no fitted diversity falls back to.

Design authority: [`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §2.1 (both names,
and the doc comments' obligations) and [`spec/calling_priors.md`](../../ng/spec/calling_priors.md)
§4 (what `θ` is, what it is not, and why the fallback must reach the run's output).

## 2. Assumptions — what the design left open and what was chosen

**The four sub-module files are created empty, carrying only a doc comment that names what lands
there and at which plan step.** The plan says `mod.rs` "declares the four files", and a `mod`
declaration does not compile without a file behind it. The alternative — add each `mod` line as
its step lands — was rejected because it makes the folder's shape invisible until Milestone F, and
because the plan asked for the opposite. Three of the four (`plug_in.rs`, `seed_spectrum.rs`,
`seed_ssr.rs`) stay empty past this session's Checkpoint B; `dirichlet_multinomial.rs` is filled by
the next two steps.

**`ExpectedHeterozygosity` is placed beside `InbreedingF`, not appended at the end of the file.**
Both are population-genetics quantities the pre-pass fits and the genotype prior consumes, and the
section they now share opens with a comment enumerating its members — which is edited in the same
commit, from four types to five. The parallel branch this plan runs beside
(`ng-calling-read-likelihoods`) adds nothing to `types.rs`, so the append discipline A1 of the
foundations plan followed — insert nowhere, to avoid a conflict — does not bind here.

**Its `DomainError` variant is its own rather than shared with `GenotypeFrequency`.** A
heterozygosity averaged over sites and the share of sites carrying one genotype are different
quantities; the file's own reasoning for giving each scalar its own variant is that a message
naming the wrong one sends the reader to the wrong fit.

## 3. Changes made

Seven files, `git diff --cached --numstat`:

| file | + | − |
|---|---|---|
| `src/ng/types.rs` | 104 | 14 |
| `src/ng/calling/genotype_prior/mod.rs` | 55 | 0 |
| `src/ng/calling/mod.rs` | 6 | 3 |
| `src/ng/calling/genotype_prior/plug_in.rs` | 9 | 0 |
| `src/ng/calling/genotype_prior/dirichlet_multinomial.rs` | 8 | 0 |
| `src/ng/calling/genotype_prior/seed_spectrum.rs` | 8 | 0 |
| `src/ng/calling/genotype_prior/seed_ssr.rs` | 8 | 0 |

### `src/ng/types.rs`

- **`ExpectedHeterozygosity(f64)`** — private field, `try_new` through the shared
  `checked_probability` predicate, `#[inline] get()`. The doc comment carries the three things
  the architecture asks it to: that it is the chance two chromosomes drawn from the cohort differ
  at an ordinary site; that it is **not** the non-reference rate, which books the reference
  accession's own quirks as cohort polymorphism; and that the STR path's diversity is a separate
  number this one must never stand in for — the substitution production's STR path made with a
  fixed SNP-scale constant.
- **`DomainError::ExpectedHeterozygosity(f64)`**, inserted beside `DomainError::InbreedingF` so
  the variant order follows the type order.
- **`DEFAULT_SPECIES_DIVERSITY_FALLBACK: f64 = 1e-3`** — port of production's
  `DEFAULT_DIVERSITY_PRIOR` ([`diversity.rs`](../../../../src/var_calling/diversity.rs)), value
  and reasoning: weakly informative, overridable, and a run that lands on it must say so, because
  a run on a species-range guess and a run on a measured diversity are otherwise
  indistinguishable in the output.
- The section-header comment above the four parameter scalars now says five types and four
  fractions, and its sentence claiming the consuming steps do not exist yet is dropped — step 8
  is the folder this same commit creates.

### `src/ng/calling/genotype_prior/mod.rs`

Module documentation only. It states what the step produces (one log-probability per candidate
genotype, before any read is looked at), the two sources of that belief, the two functions the
module is made of, and what a *concentration* is in the one reading that makes the cohort term
obvious — chromosomes the prior behaves as though it had already seen, so observed allele copies
add straight onto it. Then the four files and which plan step fills each.

### `src/ng/calling/mod.rs`

`pub mod genotype_prior;` beside `pub mod genotype_table;`, and the paragraph listing the four
sub-modules now says which one has arrived and points at its plan.

## 4. Tests added

No new test function for the scaffold — there is nothing to call. Four assertions and one new
test in `src/ng/types.rs`'s `mod tests`:

| test | what it pins |
|---|---|
| `the_constrained_rates_accept_both_endpoints` (extended) | `0.0` and `1.0` both construct. Zero is the fully invariant cohort, a real answer, so a half-open check would reject valid data. |
| `each_constrained_rate_rejects_out_of_range_in_both_directions` (extended) | `-0.5` and `1.5` are refused, each carrying its own `DomainError` variant with the offending value. |
| `the_constrained_rates_reject_nan_and_the_infinities` (extended) | `NaN`, `+∞` and `-∞` all fail to construct. |
| `the_species_diversity_fallback_is_a_constructible_heterozygosity` (new) | The fallback constant is inside the range of the type it seeds. It is the one constant in the file that is a *diversity* rather than a bound, so an edit that made it a percentage or a per-kilobase rate would be caught here rather than at the first run with no fitted diversity. |

## 5. Validation

Run in the dev container (`./scripts/dev.sh`, Apple `container` on macOS), from this worktree's
own copy of the script:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 13.35s` |
| `cargo test --lib` | 0 | `test result: ok. 4005 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 639.43s` |

4,005 against the 4,004 the branch inherited from `ng-calling-foundations` — the one new test
below. The four in-scope tests, from `cargo test --lib ng::types::tests`:

```
test ng::types::tests::each_constrained_rate_rejects_out_of_range_in_both_directions ... ok
test ng::types::tests::the_constrained_rates_accept_both_endpoints ... ok
test ng::types::tests::the_constrained_rates_reject_nan_and_the_infinities ... ok
test ng::types::tests::the_species_diversity_fallback_is_a_constructible_heterozygosity ... ok
```

The three aggregate gates that were already red before this branch existed — 18 clippy errors in
two benches and one example, a bench panic in `benches/psp_writer_perf.rs`, and 17 unresolved
intra-doc links — are untouched by it and were not re-run.

## 6. Trade-offs and follow-ups

- **Three empty files** until Milestones E and F. See §2.
- **Nothing consumes `ExpectedHeterozygosity` yet.** Step D2 is where the projection returns one;
  until then it is exercised only by its own tests. In particular **no code path yet converts the
  pre-pass's `JointFit::expected_heterozygosity`, a bare `f64`, into the newtype** — that crossing
  lands with D2 and is where the range check first meets real data.
- **The fallback constant is a bare `f64`, not an `ExpectedHeterozygosity`.** The architecture
  writes it that way, and a `const` of a type with a private field cannot be built through
  `try_new`. The new test is what keeps the two in step.
