# Applying the A2 review — ng calling foundations, step A2

*2026-08-21. Branch `ng-calling-foundations`. Input:
[`ng_calling_a2_2026-08-21.md`](../reviews/ng_calling_a2_2026-08-21.md). Every finding in that
report is accounted for below.*

## Findings table

| id | severity | decision | status |
|---|---|---|---|
| M1 — two wrong sorts survive every fixture; no property test | Major | Apply | **Applied** (both halves: the property test and the fixture swap) |
| M2 — `new(vec![])` succeeds; the invariant is prose only | Major | Apply | **Applied with adaptation** — `assert!`, not `try_new` (§ below) |
| Mi1 — a test name promising a universal two fixtures cannot show | Minor | Apply | **Applied** (renamed and the limit stated) |
| Mi2 — the `into_boxed_slice` "no copy" claim is not quite true | Minor | Apply | **Applied** (documentation) |
| Mi3 — the derived `Ord` has no stated meaning | Minor | Apply | **Applied** (documentation) |
| Missing test — `a_genotype_new_is_idempotent_on_its_own_alleles` | — | Apply | **Applied** |
| Nit — the `HashSet`/`cmp` assertions are derived consequences, not coverage | Nit | Apply | **Applied** |
| Nit — `hashes` names a set of genotypes | Nit | Apply | **Applied** (`distinct_genotypes`) |
| Nit — *pass* undefined in the sentence that justifies the type's shape | Nit | Apply | **Applied** |
| Nit — "the transposition hole" names no pair at risk | Nit | Apply | **Applied** (rewritten to name it) |
| Nit — `std::cmp::Ordering::Equal` fully qualified against the file's habit | Nit | Apply | **Applied** |
| Nit — `new` does not carry the reorder in its name | Nit | Dispute | **Won't fix** |
| Nit — `impl IntoIterator` in place of `Vec` | Nit | Dispute | **Won't fix** |
| Cross-cat — `Genotype` absent from `ng/mod.rs`'s `pub use` list | — | Defer | **Deferred** to step B1 |
| Cross-cat — `arch/ng_step_interfaces.md` §2 and §5 are stale on this type | — | Defer | **Deferred** (a design-doc edit) |
| Cross-cat — the `#[inline]` asymmetry | — | Dispute | **Won't fix** |

## M1 — the property test and the fixture swap

Both halves applied.

`a_genotype_sorts_its_alleles_and_keeps_the_multiset` went into the file's existing
`proptest::proptest! { … }` block, beside the rates', `Ploidy`'s, `SsrPeriod`'s and `Phred`'s. It
draws ids from the whole `u16` range with a dense `0..4` arm (because sampling the full range alone
essentially never repeats an id, and a repeated id is what homozygous means) and one to eight
copies, then asserts three things: the result is sorted, it is the same multiset that went in, and
a rotation of the input gives an equal value.

The tetraploid fixture in `a_genotype_holds_one_allele_per_genome_copy` changed from
`[3, 0, 2, 0]` to `[0, 3, 2, 0]`. The expected value did not change. The point of the swap is in
the test's doc comment now: written largest-first and smallest-last, the fixture trips an "is this
already reversed?" fast path into sorting anyway, so it would pass against a `new` that skipped the
sort — which is exactly why that mutation survived it.

**Verified by re-running both survivors** with `./scripts/dev.sh cargo test --lib ng::types`.

The `first() > last()` fast path — now killed by **two** tests, the proptest and the swapped
fixture:

```
test ng::types::tests::a_genotype_holds_one_allele_per_genome_copy ... FAILED
test ng::types::tests::a_genotype_sorts_its_alleles_and_keeps_the_multiset ... FAILED
Test failed: not sorted: [AlleleId(0), AlleleId(64772), AlleleId(0)] at src/ng/types.rs:1323.
test result: FAILED. 45 passed; 2 failed
```

The narrowed sort key `|a| a.0 as u8` — killed by the proptest, which is the only test reaching ids
above 255:

```
test ng::types::tests::a_genotype_sorts_its_alleles_and_keeps_the_multiset ... FAILED
Test failed: not sorted: [AlleleId(64768), AlleleId(3)] at src/ng/types.rs:1323.
test result: FAILED. 46 passed; 1 failed
```

The file was restored from a copy after each run, and the `proptest-regressions` seed lines the
mutant runs left behind were reverted with `git checkout`.

## M2 — applied with adaptation: `assert!`, not `try_new`

**Two agents recommended opposite fixes.** `idiomatic` wanted `try_new -> Result<Self,
DomainError>` with an `EmptyGenotype` variant. `reliability` argued against a `Result` and for an
`assert!`. The `assert!` is what was applied, and the deciding argument is `reliability`'s, made
from this file's own conventions rather than from taste:

`types.rs` draws a deliberate line. *Constrained* newtypes wrap scalars that arrive from outside
the program — `Phred`, `Ploidy`, `ErrorRate`, `GenotypeFrequency`, `InbreedingF` — and every one of
them has `try_new` and a `DomainError` variant. *Unconstrained* ones wrap internal indices —
`ContigId`, `Position`, `AlleleId` — and `AlleleId`'s doc states the policy outright: "any `u16` is
a legal index at the type level, and an out-of-range id is caught when the table is read — so the
field is public and there is no checked constructor." A `Genotype` is a container of those internal
indices. A `Result` would put it on the wrong side of that line, and on a path that cannot fail it
gets discharged with `.expect()` at the one call site — which relocates the panic without adding a
guarantee, and buys a `DomainError` variant whose message nothing will ever read.

The `assert!` answers `idiomatic`'s actual complaint, which was that the file rules ploidy zero
illegal in `Ploidy` and legal in `Genotype`. It now refuses in both. The price is one length
comparison per sample per locus, on the loop's last pass only — the type's own doc says a
`Genotype` is minted there and nowhere else — so it is not on the EM inner loop.

```rust
assert!(
    !alleles.is_empty(),
    "a genotype holds one allele per genome copy, and the smallest genome has one \
     copy — an empty multiset is not a haploid call, it is a sample with no genome"
);
```

with `a_genotype_cannot_be_built_from_no_alleles_at_all`, a `#[should_panic(expected = "one allele
per genome copy")]` test. **This is a real `assert!`, not a `debug_assert!`** — it fires in release
too, so the `should_panic` test pins behaviour a run will actually have, which is the distinction
this project has recorded against itself before.

`new` gained a `# Panics` section saying what the check stops: a `GT` field naming no allele at
all, written to a VCF with nothing between the genotype table and the writer having objected.

## Mi1 to Mi3 and the nits

`a_genotype_holds_one_allele_per_genome_copy_at_any_ploidy` lost its `_at_any_ploidy`, and its doc
comment now states what a point fixture cannot show — that "no ceiling on ploidy" is a claim about
every length, and unlike `Ploidy`, whose domain is a finite `u8` and so can be enumerated, a
genotype's length domain is not. It names the property test as reaching eight copies and says the
guarantee beyond that rests on `Box<[AlleleId]>`, not on a test.

The constructor's `Vec` rationale is now accurate: `into_boxed_slice` reuses the buffer when the
vector is exactly sized and reallocates when it is not, so a caller pushing one id per genome copy
into a fresh `Vec` pays one copy of a handful of `u16`s — the signature is chosen for clarity, not
for that. The `sort_unstable` sentence now gives the real reason no fixture could tell it from a
stable sort: two entries that compare equal are the same bit pattern.

The type doc gained a paragraph on the derived `Ord`: lexicographic over the sorted entries, so
`0/0` before `0/1` before `1/1`, and at mixed ploidy a shorter genotype before a longer one sharing
its prefix. It exists for a deterministic output order, not to rank genotypes by anything genetic.

"Minted from a row only at the final pass" became "on the loop's last pass, when the locus's calls
are written out". The "transposition hole" sentence now names the pair at risk — a bare `u8` ploidy
is interchangeable at the type level with a bare `u8` mapping quality or base quality, and this
file holds three such types — instead of gesturing at a term whose other uses in the crate all name
their pair. The empty-multiset half of that paragraph now points at `new`'s outright refusal rather
than at an argument about the minter.

`hashes` became `distinct_genotypes`; the two assertions beside it carry a comment saying they
follow from the `alleles()` comparison above and are there to fail if `Genotype` ever grows a
second field or a hand-written impl; and `Ordering` is imported in the test as the file's other
ordering test does.

## The three not applied

- **`from_unsorted` in place of `new`** — disputed, on `naming`'s own argument: `new` takes the
  `Vec` by value, so no caller can observe their own vector reordered, and `from_unsorted` would
  imply a `from_sorted` sibling that will never exist. The file also reserves `try_new` for
  constructors that *reject*, and this one does not.
- **`impl IntoIterator<Item = AlleleId>` in place of `Vec<AlleleId>`** — disputed. It is the more
  general bound and it buys no allocation, only characters. The agent that raised it filed it as a
  nit and noted `&[AlleleId]` would be the worst of the three.
- **`#[inline]` on `new`** — disputed. `new` sorts and allocates and runs once per sample per
  locus; `alleles()` is a field borrow. The asymmetry is the point.

## Validation

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.80s` |
| `cargo test --lib ng::types` | 0 | `test result: ok. 47 passed; 0 failed; 0 ignored; 0 measured; 3915 filtered out` |
| `cargo test --all-targets --all-features` | 101 | lib `test result: ok. 3951 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 690.15s`; every integration-test binary ok; the run then hits the same **pre-existing** `benches/psp_writer_perf.rs:386` panic A1's review records |

The step's own surface, six tests:

```
test ng::types::tests::a_genotype_keeps_repeated_alleles ... ok
test ng::types::tests::a_genotype_holds_one_allele_per_genome_copy ... ok
test ng::types::tests::a_genotype_new_is_idempotent_on_its_own_alleles ... ok
test ng::types::tests::a_genotype_cannot_be_built_from_no_alleles_at_all - should panic ... ok
test ng::types::tests::a_genotypes_alleles_are_a_multiset_not_a_sequence ... ok
test ng::types::tests::a_genotype_sorts_its_alleles_and_keeps_the_multiset ... ok
```

The change is now **+256 / −0** in `src/ng/types.rs` (`git diff --stat`), against the **+128 / −0**
the review was given.

## Follow-ups this run created

1. **`ng/mod.rs`'s curated `pub use types::{…}` list** names none of `AlleleId`, `Phred` or
   `Genotype` — nor `LogProb`, `Position`, `GenomeRegion` or `ReadGroupId`, so it is a partial
   convenience list rather than a contract. Step B1 wires `calling/` into that file and is where to
   decide.
2. **`arch/ng_step_interfaces.md` §5 still records that a `Genotype` is "reached via `alleles()` /
   `ploidy()` / `is_homozygous()`"**, and two of those three do not exist. A design-doc edit, out of this loop's remit.
