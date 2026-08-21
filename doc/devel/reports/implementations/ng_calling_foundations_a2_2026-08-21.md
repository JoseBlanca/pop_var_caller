# ng calling foundations — A2: `Genotype`

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, on top of A1 (`2cf8be6e`).
Step A2 of [`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md), Milestone A.*

## 1. Plan

Add the opaque output multiset the calling loop mints at its final pass:
`Genotype(Box<[AlleleId]>)`, alleles held sorted so that equal genotypes compare equal, with a
constructor that sorts and an `alleles()` accessor. Design authority:
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 and §Module home, and
[`arch/ng_step_interfaces.md`](../../ng/arch/ng_step_interfaces.md) §2.

## 2. Assumptions — what the design left open and what was chosen

**`new` is infallible and does not refuse an empty vector.** An empty genotype is meaningless —
every [`Ploidy`](../../../../src/ng/types.rs) in this codebase is at least one, and the type
refuses zero. The choice made is to state the contract in prose and let the one minter uphold it:
a `Genotype` is built by expanding a genotype-table row, and a row holds exactly `ploidy` copies,
so the empty case has no producer. Refusing it would cost a `Result`, a third `DomainError`
variant, and a fallible call at the final pass, for a failure nothing can reach. **This is the
assumption most worth challenging at review**, and it is the one the reliability agent was asked
about directly.

**Two accessors the interfaces sketch has were deliberately not built**, and both omissions are
argued in the type's doc comment where a caller will meet them:

- **No `ploidy()`.** `arch/ng_step_interfaces.md` §2 sketches `pub fn ploidy(&self) -> u8`. One
  returning a bare `u8` would re-open exactly the transposition hole the `Ploidy` newtype exists to
  close — §1's own argument for having the newtype at all. One returning `Ploidy` would have to be
  fallible, because `Ploidy` refuses zero and this type does not refuse an empty multiset.
  `genotype.alleles().len()` says the same thing without either problem.
- **No `is_homozygous()`.** The same sketch has one, but
  [`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §3.2 makes
  `GenotypeTable::homozygous_allele_for` the **one** homozygous test — the plan's step C1 repeats
  it: "nothing else may decide homozygosity, so the above-diploidy spec has one place to change".
  A second test on this type would be the place the two silently diverge once ploidy above two is
  specified.

**Where the type sits in the file.** A new section at the end, banner and all, parallel to the
existing "The motif — STR domain vocabulary, shared across steps" section, rather than at the end
of the "Scalar newtypes" section where A1's two scalars went. A genotype is not a scalar: it owns a
heap buffer and is not `Copy`. Appending at the end of the file is also the region furthest from
anything the parallel `ng-calling-prerequisites` branch touches.

## 3. Changes made

One file, purely additive: `src/ng/types.rs`, **+128 / −0** (`git diff --stat`).

- **`Genotype(Box<[AlleleId]>)`**, private field, deriving
  `Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug` — the file's integer-newtype derive set
  minus `Copy`, which a heap-owning type cannot have.
- **`Genotype::new(mut alleles: Vec<AlleleId>) -> Self`** — sorts with `sort_unstable`, then
  `into_boxed_slice`. Owned `Vec` because the sort happens in place and the buffer is handed
  straight over; `sort_unstable` because these are plain indices, so no two equal entries are
  distinguishable and stability buys nothing.
- **`Genotype::alleles(&self) -> &[AlleleId]`**.

## 4. Tests added

Three, in `src/ng/types.rs`'s `mod tests`:

| test | what it pins |
|---|---|
| `a_genotypes_alleles_are_a_multiset_not_a_sequence` | The property the design rests on: `[1, 0]` and `[0, 1]` are the same genotype. Checked four ways — `PartialEq`, the sorted slice, a `HashSet` collapsing to one entry, and `Ord` returning `Equal` — because the derived `Hash` and `Ord` must agree with the `PartialEq` a reader assumes, and only the sorting makes them. |
| `a_genotype_keeps_repeated_alleles` | A multiset, not a set. Two copies of one allele is what homozygous means, so a `dedup` after the sort, or a `HashSet` reached for because "order does not matter", would turn every homozygote into a haploid call — and the test asserts a diploid homozygote is not equal to a haploid call. |
| `a_genotype_holds_one_allele_per_genome_copy_at_any_ploidy` | Haploid and tetraploid, so the type imposes no ceiling — a hexaploid crop region is in scope, as `Ploidy`'s own doc says. The tetraploid case also exercises sorting past the two entries every other test uses. |

## 5. Validation

Run in the dev container (`./scripts/dev.sh`), verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 4.46s` |
| `cargo test --lib ng::types` | 0 | `test result: ok. 44 passed; 0 failed; 0 ignored; 0 measured; 3915 filtered out; finished in 0.02s` |

*(Those are the numbers as submitted for review. The review added three tests and changed one
fixture; the post-fix state is in [the fix-application report](ng_calling_foundations_a2_fixes_2026-08-21.md).)*

```
test ng::types::tests::a_genotype_keeps_repeated_alleles ... ok
test ng::types::tests::a_genotype_holds_one_allele_per_genome_copy_at_any_ploidy ... ok
test ng::types::tests::a_genotypes_alleles_are_a_multiset_not_a_sequence ... ok
```

The two aggregate gates that are red on `main` — `cargo clippy --all-targets` and
`cargo test --all-targets`, in benches and examples this branch does not touch — are unchanged by
this step and are recorded in
[A1's report](ng_calling_foundations_a1_2026-08-21.md) and in `PROJECT_STATUS.md`.

## 6. Trade-offs and follow-ups

- **The empty-genotype contract is prose, not code** — see §2. If a second minter ever appears, the
  contract needs a constructor that enforces it.
- **Nothing consumes `Genotype` yet.** `SampleGenotypeCall` carries it, and that type arrives at
  step B3.
- **No `Display`.** A genotype renders as `0/1` in a VCF, but that spelling belongs to the writer
  (step 12), which also needs the allele table to turn ids into sequence. Deliberately not here.
