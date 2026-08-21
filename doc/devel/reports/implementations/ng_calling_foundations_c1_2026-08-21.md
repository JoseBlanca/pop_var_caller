# ng calling foundations — C1: `GenotypeTable`, `GenotypeIdx` and the flat views

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, on top of B3 (`af458f44`).
Step C1 of [`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md), the first of
Milestone C's two steps.*

## 1. Plan

Add the table that says which genotypes a locus's candidate alleles make, and the three flat
quantities every step of the calling loop reads for each of them — allele counts, log multinomial
coefficients, and the homozygous lookup — built once per `(ploidy, allele count)` shape and cached.
The port of production's `GenotypeShape`
([`shape.rs:42`](../../../../src/var_calling/posterior_engine/shape.rs),
[`:76`](../../../../src/var_calling/posterior_engine/shape.rs)), as
[`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2 and §8 write it, with
`homozygous_allele_for` typed `Vec<Option<AlleleId>>` and named as the caller's **one** homozygous
test ([`arch/calling_priors.md`](../../ng/arch/calling_priors.md) §3.2,
[`spec/calling_priors.md`](../../ng/spec/calling_priors.md) §3.3). Plus `GenotypeIdx`, the loop's
working currency, and `GenotypeTableView<'_>`, the flat borrow the prior and the likelihood both
take. Value parity against production is step C2 and is not in this commit.

## 2. Assumptions and deviations

**The module is `pub mod genotype_table;`, not `mod`.** The plan's B1 text said `calling/mod.rs`
declares it; it could not, because the file did not exist. `pub` rather than private for two
reasons: the bake-off harnesses the last plan of the six builds live in `examples/`, outside the
crate, and a private module would have to be re-exported for them; and a private module holding
only `pub` items that nothing in the crate calls yet is a `dead_code` warning per item, which would
have to be silenced with an attribute that then outlives the reason for it. No re-export was added
to `calling/mod.rs` — one path to each type, `crate::ng::calling::genotype_table::GenotypeTable`.

**The enumeration is built directly in VCF order rather than generated and sorted.** Production
enumerates every non-decreasing tuple and then sorts them with a comparator that walks the tuple
backwards
([`per_group_merger.rs:522`](../../../../src/var_calling/per_group_merger.rs)). The recursion here
picks each genotype's **highest** allele at the outer level and recurses on the copies below it,
which emits the same order with no comparator. Same values, and one fewer thing to get subtly
wrong; C2 pins the claim against production's own function rather than against this reasoning.

**Three refusals the architecture does not mention**, each because the alternative is a wrong
answer rather than a crash:

- **No alleles.** `CandidateAlleles` cannot be built without its reference, so a zero-allele shape
  did not come from a locus. Production does not check it and panics anyway, inside
  `chunks_exact(0)`, with a message about chunk sizes.
- **More alleles than an `AlleleId` can name.** The homozygous lookup holds allele ids, so a table
  above 65,536 alleles could not say which allele its last rows are homozygous for. This mirrors
  the check `CandidateAlleles::admit` already makes.
- **A genotype count that does not fit in a `usize`.** The count is `C(alleles + ploidy − 1,
  ploidy)`, which at ploidy 40 over 64 alleles is about 6.1 × 10²⁸. Refused by name rather than
  overflowing inside a capacity calculation.

**Row lookups (`allele_counts_of`, `log_multinomial_coeff_of`, `homozygous_allele_of`) were added
beyond the plan's list.** Without one, `GenotypeIdx` is a type with nothing to do. They return
`Option` and do not index, for the reason `CandidateAlleles::bases_of` does not: an index minted
against another shape is a legal `u32` here, and indexing would either panic or return a real but
wrong row.

**Minting an owned `Genotype` from a row is deliberately absent.** The architecture says the multiset
is minted from the winning row on the final pass; which pass that is, and what does the minting, is
the calling-loop plan's. Recorded as a follow-up rather than built here.

**The cache bounds are production's — ploidy ≤ 8, alleles ≤ 16 — and the degradation past them is
documented and tested.** Outside the bounds every `build` returns a fresh `Arc` holding the same
values; only the sharing is lost. Two tests pin it, one on each axis, and two more pin the last
shape on each axis that *is* shared. Naming what this caller gives up outside the bound is the
`CLAUDE.md` range rule: a dodecaploid locus, or one with 20 candidates, is correct and uncached.

## 3. Changes made

Two files, **+872 / −0** and **+8 / −0** (`wc -l` on the new file; `git diff --stat` on the other):

- **`src/ng/calling/genotype_table.rs`** (new, 872 lines, of which the test module is 409):
  - `GenotypeIdx(pub u32)` with `get()`.
  - `GenotypeTable` — six private fields (`ploidy`, `allele_count`, `n_genotypes`, and the three
    tables), `build(Ploidy, usize) -> Arc<Self>` with the thread-local slot-array cache, `view()`,
    six accessors, three row lookups.
  - `GenotypeTableView<'a>` — private fields, no public constructor, six accessors.
  - Private helpers: `genotype_count` (checked binomial), `enumerate_allele_counts` +
    `push_genotypes_with_highest_allele_below`, `log_factorial`, `log_multinomial_coefficient`,
    `homozygous_allele`.
- **`src/ng/calling/mod.rs`** — `pub mod genotype_table;` and a paragraph in the module doc saying
  why the file sits beside `mod.rs` rather than inside a sub-module.

`log_factorial` keeps production's summation order (`ln 2 + ln 3 + … + ln n`) on purpose: C2
compares the coefficients value for value, and a different order differs in the last bits.

## 4. Tests added

Twenty-two, in `src/ng/calling/genotype_table.rs`'s `mod tests`. Count and result from
`cargo test --lib ng::calling::genotype_table`: `22 passed; 0 failed`.

| test | what it pins |
|---|---|
| `diploid_triallelic_genotypes_come_out_in_vcf_order` | The six genotypes as allele lists, by hand: `0/0`, `0/1`, `1/1`, `0/2`, `1/2`, `2/2`. |
| `diploid_biallelic_table_holds_the_three_genotypes_with_their_copy_counts` | The flat row-major layout itself, `[2,0, 1,1, 0,2]`. |
| `tetraploid_biallelic_table_runs_from_four_reference_copies_to_four_alternative` | Ploidy above two, by hand — in this caller's range, and not covered by any diploid fixture. |
| `a_haploid_locus_has_one_genotype_per_allele_and_all_of_them_are_homozygous` | Ploidy 1: the other end of the range, where every genotype is homozygous by the prior's definition. |
| `every_genotype_row_accounts_for_exactly_the_ploidy_copies` | Over ploidy {1,2,3,4,6} × alleles 1–6: every row sums to the ploidy. |
| `the_table_holds_one_row_per_multiset_of_alleles` | Nine shapes with the count worked out by hand, up to 490,314 genotypes at ploidy 8 over 16 alleles; and that all three tables have that many rows. |
| `no_genotype_is_enumerated_twice` | The check a count cannot make: a duplicate and a miss would cancel. |
| `diploid_biallelic_coefficients_are_one_two_one_in_logs` | `0`, `ln 2`, `0`, by exact equality. |
| `every_coefficient_is_the_log_of_the_exact_number_of_orderings` | Over ploidy {1,2,3,4,6} × alleles 1–6: each coefficient against `ln` of the exact integer `ploidy!/∏counts!` computed with `u128` factorials — a route that does not repeat the implementation's own summation. |
| `the_homozygous_lookup_names_an_allele_exactly_where_every_copy_is_that_allele` | Over ploidy {1,2,3,4} × alleles 1–6: `Some(a)` exactly where one allele holds every copy, and exactly one homozygote per allele. |
| `diploid_triallelic_homozygous_lookup_marks_the_first_third_and_last_genotypes` | The same rule spelled out where a reader can check it against the VCF order by eye. |
| `a_row_lookup_returns_that_genotypes_counts_coefficient_and_homozygous_allele` | The three lookups agree with the tables, on a homozygous row and a two-allele one. |
| `a_row_lookup_past_the_end_of_the_table_finds_nothing` | An index from a wider shape, and `u32::MAX`, resolve to `None` on all three lookups. |
| `the_view_carries_the_same_shape_and_the_same_three_tables` | The borrow matches the table, and the three slices are parallel. |
| `two_builds_of_one_cached_shape_return_the_same_table` | `Arc::ptr_eq` — the point of the cache. |
| `different_shapes_get_different_tables` | The slot index is a function of both axes, not one. |
| `a_shape_past_the_cache_bounds_is_rebuilt_each_time_with_the_same_values` | Both axes, one step past (ploidy 9, and 17 alleles): fresh `Arc`, equal values. |
| `the_widest_and_deepest_cached_shapes_are_still_shared` | Both axes at the bound (ploidy 8, and 16 alleles): shared. Without it, either comparison could be loosened and nothing would fail. |
| `each_thread_keeps_its_own_cached_tables` | A second thread builds its own, and gets equal values. |
| `a_table_of_no_alleles_is_refused` | With its message. |
| `a_table_wider_than_an_allele_id_can_name_is_refused` | With its message. |
| `a_shape_with_more_genotypes_than_a_usize_can_count_is_refused` | With its message. |

## 5. Validation

Run in the dev container, verbatim:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile [unoptimized + debuginfo] target(s) in 4.34s` |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 43 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --lib ng::calling::genotype_table` | 0 | `test result: ok. 22 passed; 0 failed; 0 ignored; 0 measured; 3983 filtered out` |

The branch's baseline before this step, measured the same way: `cargo test --lib` →
`3972 passed; 0 failed; 11 ignored`, of which 21 were in `ng::calling`.

The three aggregate gates that are red on `main` are unchanged and untouched by this step —
`cargo clippy --all-targets` (18 errors in two benches and one example),
`cargo test --all-targets` (`benches/psp_writer_perf.rs:386` index-out-of-bounds), and
`cargo doc --no-deps --lib` (17 unresolved intra-doc links).

## 6. Tradeoffs and follow-ups

- **The parity test is C2** and is what makes the enumeration-order claim in this file's module doc
  a measurement rather than an argument. Until it lands, `genotypes_as_allele_lists` fixtures are
  the only thing holding the VCF order.
- **`nonzero_pairs` was not ported.** Production carries it to skip zero-count cells in its E-step;
  the architecture says it comes along "only if profiling asks", and nothing has asked
  ([`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §8).
- **Minting a `Genotype` from a `GenotypeIdx`** belongs to the calling-loop plan; see §2.
- **A wide, high-ploidy shape is large**, and the cache holds one per shape per thread: 16 alleles
  at ploidy 8 is 490,314 genotypes, so its allele-count table alone is 490,314 × 16 × 4 bytes ≈
  31 MB. Inside the cache bounds, so it would be kept. No run asks for that shape today — the
  candidate cap ships at 6 alleles — and nothing here bounds it.
