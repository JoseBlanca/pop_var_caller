# ng calling foundations — B3: `LocusInference` and `SampleGenotypeCall`

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, on top of B1+B2 (`56b65ae0`).
Step B3 of [`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md), the last step of
Milestone B.*

## 1. Plan

Add the outcome type the calling loop produces at each locus, as
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 writes it: the region, the final
allele table, one call per sample in run order, the cohort's expected allele copies, and four
pieces of evidence for how the answer was reached — `converged`, `passes`, `weakest_provenance`,
`seed_diversity_unreachable`. Plus `SampleGenotypeCall`, the pair of a genotype and its quality.
Plain data; the loop plan fills it.

## 2. Assumptions and deviations

**Public fields, plus a checked `new`.** The plan says "plain data, no logic", and this file's
house style for record types is public fields — `SsrDetail` and `Estimate<T>` are the same shape.
But two things about a `LocusInference` are checkable and worth checking where both halves are in
scope, so `new` asserts them:

- **The copies are one entry per allele.** `ExpectedAlleleCopies` is already built against *an*
  allele table, so what remains is the residue: copies built against a **different** table that
  happened to be the same width. This is the one place both tables are in scope.
- **`passes > 0`.** Every locus takes at least one pass, so zero is a counter that was never
  incremented rather than a locus that converged instantly — and it would distort the very
  distribution the pass cap is to be set from.

**A struct literal bypasses both**, which is the cost of keeping the fields public. Flagged to the
reviewers as the decision most worth challenging: the alternative is eight private fields and eight
accessors on a record whose whole purpose is to be read.

**No relation is checked between `converged` and `passes`.** A capped locus should have `passes`
equal to the cap, but the cap is a run-level configuration this type does not see
(`CallingLoopConfig`, `arch/calling_em_loop.md` §2.1, built by the loop plan). Checking it here
would mean passing the config in.

**`seed_diversity_unreachable` is not forced to `false` on the SNP/indel path**, though the doc
says it is always false there and `alleles.kind()` is in scope. It is an assertion the type could
make; flagged to the reviewers rather than made, because the STR-bundle kind's behaviour is
deferred (spec §11) and I would be guessing which side of the line it falls on.

## 3. Changes made

One file, **+287 / −2** (`git diff --stat`): `src/ng/calling/mod.rs`.

- **`SampleGenotypeCall`** — `genotype: Genotype`, `genotype_quality: Phred`, both public.
- **`LocusInference`** — the eight fields above, all public, with doc comments carrying why each
  cannot be reconstructed downstream.
- **`LocusInference::new`** — eight positional arguments, the two asserts above.
- Two import lines: `Provenance`, and `GenomeRegion`/`Genotype`/`Phred` added to the existing
  `types` import.

## 4. Tests added

Five, plus three helpers, in `src/ng/calling/mod.rs`'s `mod tests`:

| test | what it pins |
|---|---|
| `a_locus_carries_its_calls_and_the_evidence_for_how_they_were_reached` | Every field survives construction, and `per_sample` is a *sequence* — the second sample is the heterozygote, so a reversal would show. |
| `a_locus_that_ran_out_of_passes_is_emitted_with_the_flag_set` | The spec's rule: `converged: false` is emitted with its call intact, never dropped and never fatal. Also that the other two warrants travel independently of it. |
| `a_locus_cannot_carry_copies_of_a_width_its_alleles_do_not_have` | The residue check — copies built against a two-allele table, paired with a three-allele one. |
| `a_locus_cannot_report_no_passes_at_all` | The `passes > 0` assert, with its message. |
| `a_cohort_of_one_sample_is_a_locus_like_any_other` | The hardest case in this caller's stated range is representable: one sample, one call, and the cohort's copies are that sample's own. |

## 5. Validation

Run in the dev container, verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 3.83s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 17 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |

*(Those are the numbers as submitted for review; the post-review state is in the fix-application
report beside this one.)*

The three aggregate gates red on `main` are unchanged; recorded in
[A1's review](../reviews/ng_calling_a1_2026-08-21.md) §7.

## 6. Trade-offs and follow-ups

- **A struct literal bypasses `new`** — see §2.
- **`converged` against `passes`, and the SNP/indel path's `seed_diversity_unreachable`**, are both
  unchecked — see §2.
- **Nothing consumes the type yet.** The loop plan fills it; site filtering, emission and the
  writer read it.
