# ng calling foundations — C1: what the review changed

*Fix-application report, 2026-08-21. Branch `ng-calling-foundations`. Applies
[`ng_calling_c1_2026-08-21.md`](../reviews/ng_calling_c1_2026-08-21.md) to the step-C1 diff
described in [`ng_calling_foundations_c1_2026-08-21.md`](ng_calling_foundations_c1_2026-08-21.md).*

## 1. What was applied

Everything the review filed except two items held for the owner, both named in §3. The Blocker and
all five Majors are fixed; eleven of the thirteen Minors and all the actionable nits are applied.

**Six defects were provable only by mutation, and each fix was re-verified by re-running its
mutation here** — not by trusting the reviewer's own run. Each mutation below was applied to the
fixed file, the suite run, the mutation reverted:

| what was mutated | before the fix | after the fix |
|---|---|---|
| swap rows 0 and 1 whenever ploidy ≥ 3 and alleles ≥ 3 | 22 tests pass | `every_consecutive_pair_of_rows_is_in_vcf_order` FAILS; the two hand-written order fixtures still pass, so it is the new test that catches it |
| cache slot index `p * (MAX_CACHED_PLOIDY + 1) + allele_count`, **with the new `debug_assert`s deleted** so only a test could catch it | 22 tests pass | `every_cached_shape_gets_back_a_table_of_its_own_shape` FAILS |
| zero the coefficients and blank the homozygous lookup **on the uncached branch only** | 22 tests pass | `a_shape_past_the_cache_bounds_holds_the_right_values` FAILS |
| drop `checked_mul` from the running binomial in `count_genotypes` | 22 tests pass | `a_shape_that_overflows_the_running_binomial_is_refused` FAILS |
| drop the flat-table overflow check (i.e. restore the submitted code) | — | `a_shape_whose_flat_table_outgrows_a_usize_is_refused` FAILS |
| make `count_genotypes` return one more than the true count | — | the new `assert_eq!` fires: `the genotype count and the enumeration disagree at ploidy 1 over 4 alleles` |

Each of the first four is a wrong answer with nothing crashing — the class this module's own header
names.

### The Blocker and the Majors

- **B1 — the enumeration order is now pinned as a law, not at two shapes.**
  `every_consecutive_pair_of_rows_is_in_vcf_order` asserts, over ploidy {1, 2, 3, 4, 6, 8} × alleles
  1–6, that consecutive rows read as sorted allele lists are strictly increasing compared from the
  highest copy downwards — the comparator production's `genotype_order` sorts with. The two
  hand-written fixtures stay: they are what a reader checks by eye.
- **M1 — the cache can no longer hand a locus another shape's table.** Two `debug_assert_eq!`s on
  the cache hit compare the stored `ploidy` and `allele_count` against what was asked for, and
  `every_cached_shape_gets_back_a_table_of_its_own_shape` walks all 128 cached shapes in one thread
  and asks each table which shape it is. Either alone kills the mutation; both are kept because the
  assertion catches it in any test and the test catches it in a release build.
- **M2 — the table's own size is checked where the genotype count is.** `build_uncached` now
  refuses when `genotype_count × allele_count` does not fit a `usize`, in the same expression and
  with the same shape-naming message, and `enumerate_allele_counts` takes that checked total as its
  capacity rather than recomputing it.
- **M3 — the module doc no longer claims a test that is C2's.** Both places now say what this
  commit actually pins — the two hand-written orders and the law — and name step C2 as where the
  value-for-value comparison against production arrives.
- **M4 — the uncached branch is checked by value.** `a_shape_past_the_cache_bounds_holds_the_right_values`
  asserts the genotype count, first row, first coefficient and three homozygous entries at ploidy 9
  over 2 alleles and at ploidy 2 over 17. The old test stays: it pins the *bound*, this pins the
  *values*.
- **M5 — both failure routes through `count_genotypes` are tested.** The shipped fixture (ploidy 40
  over 64 alleles) exits through the final `usize::try_from`; the new one (ploidy 255 over 65,536
  alleles) overflows the running binomial itself.

### The Minors

- **Mi1** — `calling/mod.rs` re-exports `GenotypeIdx`, `GenotypeTable` and `GenotypeTableView`, so
  all of calling's shared vocabulary is named at one depth. The submitted report's `dead_code`
  argument was wrong and is corrected in §4.
- **Mi2** — the field and both accessors are now `genotype_count`; the free function that computes
  it is `count_genotypes`, a verb.
- **Mi3** — the whole-column accessor is `homozygous_alleles()` on both types, beside
  `genotype_allele_counts()` and `log_multinomial_coeffs()`; the row lookup keeps `_of`. The
  **field** stays `homozygous_allele_for`, because that is the parameter name the step-8 seam uses
  (`arch/calling_priors.md` §3.2) and it should stay greppable.
- **Mi5** — `row_of` on the view is the one place an index becomes a row; the three lookups read it
  and `GenotypeTable`'s three delegate to the view. The loop, which takes the view, now has the
  bounds-checked lookups `GenotypeIdx`'s own documentation promises.
- **Mi6** — both bounds are `pub`, and the cost sentence carries the number: the slot array is
  about 1.2 kB a thread, one filled slot at ploidy 8 over 16 alleles is 37,263,864 bytes
  (490,314 × (16 × 4 + 8 + 4)), and a thread meeting every in-bounds shape would hold about 136 MB.
  Both figures were re-derived here rather than taken from the review. The step off the bound is
  stated in genotypes rather than in a machine-specific time: 490,314 at 16 alleles against 735,471
  at 17, rebuilt per locus.
- **Mi7** — the table is built outside the `RefCell` borrow, with two short borrows around it. The
  comment says why: a future call reaching back into `build` would double-borrow, and the release
  profile aborts on panic.
- **Mi8** — `view()` destructures `self` exhaustively, so adding a seventh table (the architecture
  names `nonzero_pairs`) fails to compile there instead of reaching consumers as silence.
- **Mi10** — the width refusal's reason is corrected: an `AlleleId` names exactly 65,536 alleles,
  so that width is the widest one it can name and 65,537 is the first refused.
- **Mi11** — `MAX_ALLELE_COUNT` is one public constant used by the assertion and by the test, and
  the `expect` carries a `// PANIC-FREE:` comment naming the invariant.
- **Mi12** — `build_uncached` asserts that the closed-form count and the enumeration agree, so two
  independent computations of one quantity check each other at every build.
- **Mi13** — the ploidy loops now include 8, the deepest shape the cache keeps.

### Nits

`GenotypeIdx::get` has a doc comment; `only` is `only_allele_seen`; the doc says `allele_count`
where it said `n_alleles`; `count_genotypes` documents that it assumes a non-zero allele count, so
the ordering of the two checks is not accidental; the `lgamma` claim is softened to what is known
("at most a handful of `ln` calls"), since nothing measured it.

## 2. Not applied, and why

- **A `proptest` over the whole `(ploidy, allele count)` domain.** The four laws it would assert
  are each now asserted over an explicit grid, and a grid failure names the shape. A random failure
  would need a seed to reproduce.
- **A 64-thread agreement test.** The concurrency agent measured it — 64 of 64 threads agree — and
  nothing in the type can make it fail: the cache is per-thread and the value is a pure function of
  its key. `each_thread_keeps_its_own_cached_tables` already pins the property that matters.
- **`#[must_use]` on `build`, a doc example, a `criterion` bench.** The first matches the immediate
  neighbours (`types.rs` and `calling/mod.rs` have none); the other two are conventions for
  `src/ng/calling/` as a whole, not for one file, and the bench needs a caller to be worth writing.
- **The stale tree in `arch/module_layout.md`.** It is a design document, and this loop does not
  edit those. Recorded as a follow-up.

## 3. Held for the owner

Two naming and API questions the reviewers raised that are not the implementer's to settle:

1. **`GenotypeIdx` against `GenotypeIndex`.** The crate spells index newtypes `Id` or `Index` and
   this is its only `Idx`; 18 references, all in one file, so it is at its cheapest now. But the
   name is written into the plan's step C1 and into all three calling architecture documents, and
   the four plans that follow will each name it. Left as `GenotypeIdx`.
2. **A `for_locus(ploidy, &CandidateAlleles)` constructor**, which would make the loop's own call
   total and retire the zero-allele assertion for it. Whether the loop wants a table before its
   allele set is final is the calling-loop plan's question.

## 4. A correction to the C1 implementation report

Its §2 gives two reasons for `pub mod genotype_table;`, and one is wrong: it says a private module
of `pub` items nothing calls yet is a `dead_code` warning per item, presenting that as forcing the
choice. A private module **plus a re-export** compiles clean under `-D warnings` — the re-export is
what makes the items reachable. The `examples/` reason stands on its own, and the re-export is now
present regardless. Measured by the reviewer and not disputed: without a re-export, `mod
genotype_table;` raises 13 dead-code warnings.

## 5. Validation

Run in the dev container after every fix, verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 3.09s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 49 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --lib` | 0 | `test result: ok. 4000 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out` |

Twenty-eight tests in `ng::calling::genotype_table`, six more than the 22 submitted; 49 in
`ng::calling` against the 43 submitted and the 21 the branch had before C1.
