# ng genotype prior — A2: the local types and the step-8 seam

*Implementation report, 2026-08-21. Branch `ng-calling-prior`, worktree
`../pop_var_caller-calling-prior`. Step A2 of
[`calling_prior.md`](../../ng/impl_plan/calling_prior.md), Milestone A.*

> **This report describes the code as submitted for review.** The review changed a great deal of
> it — most visibly, the seam's eight parameters became one checked bundle, `PriorRow`, whose
> constructor runs the shape checks. What landed is in
> [the review](../reviews/ng_calling_prior_a2_2026-08-21.md) and
> [the fixes](ng_calling_prior_a2_fixes_2026-08-21.md).

## 1. Plan

Four declarations in
[`src/ng/calling/genotype_prior/mod.rs`](../../../../src/ng/calling/genotype_prior/mod.rs), and
no logic behind any of them:

- **`Concentration<'a>`** — a borrow of the caller's buffer holding one strictly positive number
  per allele; invariant checked in debug.
- **`SeedRegime`** — which of three kinds of information produced the run's starting point. **A
  branch on absence, never on cohort size.**
- **`SpectrumSeed`** — the SNP/indel starting point: two numbers for the whole run, plus the
  regime.
- **`GenotypePriorModel`** — the seam the marginalized default and the plug-in comparator both sit
  behind.

No `Result` anywhere in the module: a mis-shaped input is a caller bug, so it is an assertion, and
the structural ones are held in release.

Design authority: [`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §1.1, §2.2, §2.3,
§3.2; [`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §1 (what a concentration is), §4
and §4.1 (what the regimes mean), §8 (numerics, determinism, no allocation).

## 2. Assumptions and deviations — what the design left open and what was chosen

Three, and the second is the one that changes a signature the architecture writes down.

### The trait takes `Concentration<'_>`, not a bare `&[f64]`

Arch §3.2's sketch spells the first parameter `concentration: &[f64]`. It is the newtype instead,
because otherwise **`Concentration` has no consumer anywhere** — arch §2.2 mints it and §3.2 does
not use it — and because the newtype is the only place its documented invariant ("every entry
`≥ MIN_ALT_CONCENTRATION`, checked in debug") can live. Arch §0's rule is that the prior takes
**flat slices rather than the loop's types**; `Concentration` is this document's own type and is
itself a borrow of a flat slice, so the rule is kept.

### The trait gained a per-allele scratch parameter, and this is a real addition

Arch §3.2's sketch has six parameters. The implementation has seven, the extra being
`per_allele_scratch: &mut [f64]`.

**The reason is spec §8**, which the arch sketch did not spell: *"Nothing may allocate inside the
per-sample loop — the caller hands in scratch **sized by allele count and genotype count**, and
the prior fills it."* `out` is the genotype-count half; this is the allele-count half, and without
it the first thing B1's ported primitive does is allocate a `Vec` of `lgamma(α_a)` per call —
which is exactly the allocation production lifted out of its own loop after a profile put the
allocator's self-time at about one cycle in six.

**What it costs to do without it**, so the trade is on the record rather than asserted: the
per-allele `lgamma` baseline would have to be recomputed at every (genotype, non-zero allele)
pair. At a diploid locus with six alleles — the shipping candidate cap — that is **72 `lgamma`
calls per sample per pass against 42**. Counted off the table itself rather than by hand: a probe
over `GenotypeTable::build(Ploidy(2), 6)` reports `alleles=6 genotypes=21 nonzero_pairs=36`, so
the cached form is 6 baselines plus 36 and the uncached form is two per pair. It is a count of
calls, not a profile.

**Blast radius:** the calling loop's `CallingScratch` needs an allele-count-sized slot. The plan
already records those slot sizes as that plan's impl-time confirmation, so this is inside what it
anticipated. Raised at Checkpoint A so arch §3.2 can be brought into line.

### The parameter is `homozygous_alleles`, not `homozygous_allele_for`

The genotype table's field keeps the architecture's name and its accessor does not
(`GenotypeTableView::homozygous_alleles`). The trait's parameter follows the accessor, so a call
site reads `view.homozygous_alleles()` into `homozygous_alleles`.

### Not done: the review's open question about the fallback's provenance

The A1 review asked whether the species-range fallback should move down beside `SeedRegime`, or
whether the only public door should be a constructor returning value and regime together. **Neither
is applied, and the reason is that `SpectrumSeed` already does the job.** A seed is a struct of
three fields with no other constructor, so **no seed can be built without a regime**, and the
regime is what travels to the output. What the pairing function would additionally prevent — a
consumer reading `ExpectedHeterozygosity::SPECIES_FALLBACK` directly and never building a seed —
is prevented instead by `project_spectrum_seed` (step D2) being the only sanctioned way to obtain
a run's starting point. Recorded rather than decided quietly; it is on the Checkpoint A list.

## 3. Changes made

One file, purely additive: `src/ng/calling/genotype_prior/mod.rs`, **+525 / −0**
(`git diff --numstat`), of which the test module is roughly the second half.

- **`Concentration<'a>(&'a [f64])`** — `new` (asserting non-empty in release, and every entry
  finite and at least `MIN_ALT_CONCENTRATION` in debug), `get`, `allele_count`. `Copy`, because it
  is a borrow.
- **`SeedRegime`** — `FittedSpectrum { regularizer_site_weight: f64, data_dominated: bool }`,
  `NeutralShape`, `FallbackDiversity`. The last variant's doc links
  `ExpectedHeterozygosity::SPECIES_FALLBACK`, so the constant and the thing that must report it
  name each other.
- **`SpectrumSeed`** — `alpha_ref`, `alpha_alt_total`, `regime`, all public: the invariant is that
  all three are present, which three public fields already enforce.
- **`GenotypePriorModel`** — one method, `genotype_log_priors`, taking the concentration, the
  genotype table's three flat views, the sample's inbreeding coefficient, the per-allele scratch
  and the output row.
- **`assert_row_shapes`** — the four structural checks, in one function so the two implementations
  cannot disagree about them, held in release.

`MIN_ALT_CONCENTRATION` is imported from `crate::genetics` with its reasoning, as arch §6
directs. **`PROBABILITY_FLOOR` is not**, though the same row of that table names it: nothing in
A2 takes a logarithm, and an unused import is a clippy error. Its consumers are the Wright test
oracle (B2) and the plug-in comparator (F1).

### One lint is silenced, and it is load-bearing

`genotype_log_priors` has eight arguments counting `&self`, against clippy's ceiling of seven.
**Verified rather than assumed** — before the attribute was added,
`cargo clippy --lib --all-features -- -D warnings` gave:

```
error: this function has too many arguments (8/7)
```

The justification in the doc comment is that four of the eight are caller-owned buffers, and they
are parameters precisely because spec §8 forbids this seam from allocating. Bundling them into a
struct moves the count without removing anything; bundling the three table views would take the
loop's own type, which arch §7 decides against.

## 4. Tests added

Nine, all in `src/ng/calling/genotype_prior/mod.rs`.

| test | what it pins |
|---|---|
| `the_seam_takes_the_genotype_tables_views_unadapted` | **The one thing A2 can get wrong that nothing else would catch until the loop is written**: the trait names four flat views by type, and a real `GenotypeTable`'s accessors are passed straight into it with nothing adapting in between. If `homozygous_alleles()` yielded bare ids, or the counts were `u16`, this stops compiling. |
| `the_seam_is_reachable_through_a_trait_object` | `&dyn GenotypePriorModel` works, which is what lets a run select between the two implementations without the calling loop being generic over the choice. |
| `a_concentration_borrows_the_buffer_it_is_given` | The wrapper copies nothing — the pointer that comes out is the pointer that went in. |
| `an_empty_concentration_is_refused` | Refused in release: every locus has a reference allele, so a zero-length concentration is a wiring bug, and `Σα` over no entries is a division by zero waiting to happen. |
| `a_non_positive_concentration_entry_is_refused_in_debug` | `lgamma` is defined only for a positive argument. |
| `a_concentration_entry_under_the_floor_is_refused_in_debug` | The floor is `MIN_ALT_CONCENTRATION` and not merely "above zero" — a smaller value means a seed builder's flooring step was skipped. |
| `every_row_shape_check_names_the_length_that_was_wrong` | Each of the four structural checks fires on **its own** length and says which. Four cases, each shortening one slice. |
| `the_row_shape_checks_pass_on_a_well_formed_call` | The companion to the above: the fixture is well formed to begin with, so the four failures are the lengths the test changed. |
| `a_spectrum_seed_carries_the_regime_that_produced_it` | Two seeds differing only in regime are different values — the regime is part of what a seed *is*, not an annotation on it. |

## 5. Mutation results

**Six mutations run, five killed, one survived, none changed no behaviour.** Each was applied to a
copy of the reviewed file and reverted; the baseline (`9 passed`) was re-run after the last.

| mutation | outcome |
|---|---|
| the per-allele scratch length check deleted | killed — `8 passed; 1 failed` |
| the allele-counts table check relaxed from `==` to `>=` | killed — `8 passed; 1 failed` |
| the empty-concentration refusal deleted | killed — `8 passed; 1 failed` |
| the debug floor relaxed from `>= MIN_ALT_CONCENTRATION` to `> 0.0` | killed — `8 passed; 1 failed` |
| the homozygous-lookup length check deleted | killed — `8 passed; 1 failed` |
| **an implementation that never calls `assert_row_shapes`** | **survived — `9 passed; 0 failed`** |

The survivor is real and is not fixable by a test in this file. A Rust trait cannot require that a
method's body opens with a particular call, so the shape checks are an obligation the trait
documents and each implementation owes. The doc comment says so, in those words, and names the
measurement. **What closes it is a test per implementation that a mis-shaped call panics** — B2's
to write for the marginalized prior, F1's for the comparator.

## 6. Validation

Run in the dev container from this worktree's own `scripts/dev.sh`:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 2.89s` |
| `cargo test --lib ng::calling::genotype_prior` | 0 | `test result: ok. 9 passed; 0 failed; 0 ignored; 0 measured; 4018 filtered out` |
| `cargo test --lib` | 0 | `test result: ok. 4016 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 669.93s` |

The three aggregate gates already red on `main` are untouched and were not re-run.

## 7. Trade-offs and follow-ups

- **The shape checks are an obligation, not an enforcement** (§5). Closed per implementation at B2
  and F1.
- **Nothing implements the trait outside the test module.** `MarginalizedDirichletPrior` arrives at
  B2 and is the first real implementation; the stand-in in the tests writes the multinomial
  coefficient and is deliberately not a prior, so nothing in this file can be mistaken for a check
  of one.
- **`SpectrumSeed` has no builder yet.** `project_spectrum_seed` (D2) is the only thing that will
  produce one, and `seed_for_locus` (D3) the only thing that expands it onto a locus.
- **Two arch signatures now differ from the code** — §3.2's parameter list (the scratch) and
  §2.1's spelling of the fallback constant, from A1. Both are on the Checkpoint A list.
