# ng calling foundations — C2: the genotype table against production, value for value

*Implementation report, 2026-08-21. Branch `ng-calling-foundations`, on top of C1 (`cb36127e`).
Step C2 of [`calling_foundations.md`](../../ng/impl_plan/calling_foundations.md), the last step of
the plan. Includes the review of this step and the fixes it produced —
[the review](../reviews/ng_calling_c2_2026-08-21.md).*

## 1. Plan

The value-parity test the plan requires before Milestone C closes: across a grid of
`(ploidy, allele count)` shapes, compare ng's `GenotypeTable` against production's `GenotypeShape`
on the genotype count, every row's allele counts **in production's enumeration order**, every log
multinomial coefficient **to floating-point equality**, and the homozygous lookup — plus the cache
identity, two `build` calls for one shape returning one `Arc`. Its own commit, because its failure
is silent: a wrong genotype index is a wrong call, not a panic.

## 2. Assumptions and deviations

**The plan says this test calls `GenotypeShape` directly. It cannot, and this is the one thing the
next reader of the plan most needs to know.** `src/var_calling/posterior_engine.rs` declares
`mod shape;` privately, so neither `GenotypeShape` nor `shape_for` can be named from `src/ng/`
whatever their own `pub(crate)` visibility, and `:61` imports them privately too, so there is no
re-export route. Compiling a call gives ``error[E0603]: module `shape` is private`` — confirmed
three times independently on this branch. The only way to satisfy the step as written is to widen
that declaration, which would edit a frozen tree.

So the oracle is assembled from two things:

- **`genotype_order`** (`src/var_calling/per_group_merger.rs`), a `pub fn` in a `pub mod` and
  reachable. It is production's *sole* enumeration source — `collect_non_decreasing` is
  production's only enumeration recursion and every consumer, the VCF writer included, reaches it
  through this one function. It is the one thing here that runs production's own code.
- **Four transcriptions** — `GenotypeShape::build`'s fold of that enumeration into the flat count
  table, plus `log_factorial`, `log_multinomial_coefficient` and `homozygous_allele`, copied from
  `shape.rs` including its summation order.

What that gives up: if production changed its fold or one of the three formulas, this test would
keep passing while ng's table stopped matching production's. All four are closed forms rather than
behaviour, and the **order** — the one part that is a recursion with a comparator — is exactly the
part compared against production's live code.

**The grids are wider than the plan's.** The plan names ploidy 2 and 4 over alleles 1–6; that grid
is here as its own test, and two more sit beside it (see §4). The plan's grid is a strict subset of
the second, so its contribution is its total; it is kept because the plan names it.

**`src/ng/mod.rs` gained a paragraph.** It is ng's own tree, not frozen production. The freeze
paragraph now names the test-oracle exception and both instances of it, because this is the first
compile-time edge from `src/ng/` to `src/var_calling/` and the rule a reader consults said nothing
about it.

**One C1 test was strengthened here**, which is a change to a committed file and is recorded as
such: see §3.

## 3. Changes made

Three files:

- **`src/ng/calling/genotype_table_parity.rs`** (new, 315 lines) — the four transcribed oracle
  functions, `compare_against_production` (six comparisons per shape), and four tests.
- **`src/ng/calling/mod.rs`** — `#[cfg(test)] mod genotype_table_parity;` and a clause in the
  module doc saying what it is and why it is its own file.
- **`src/ng/mod.rs`** — the test-oracle exception to the production freeze.
- **`src/ng/calling/genotype_table.rs`** — one C1 test strengthened.
  `a_shape_past_the_cache_bounds_holds_the_right_values` asserted the coefficient of row 0 at
  ploidy 9. Row 0 is `[9, 0]`, whose coefficient is `ln 9! − ln 9!` — zero under the right formula
  and zero under a `log_factorial` capped at 8, at 4, or at 2. **The one row it checked was the one
  row the defect cannot touch.** It now checks row 1, `[8, 1]`, worth `ln 9`.

## 4. Tests added

Four, in `src/ng/calling/genotype_table_parity.rs`. Every shape goes through six comparisons: the
table's own declared ploidy and width, the genotype count, the whole allele-count table as one
slice, the two parallel tables' lengths, then every coefficient by bit pattern and every homozygous
entry.

| test | grid | genotypes | what it pins |
|---|---|---|---|
| `..._over_the_diploid_and_tetraploid_grid` | ploidy {2, 4} × alleles 1–6 | 308 | The plan's own grid, named as such. |
| `..._from_haploid_to_octoploid_up_to_eight_alleles` | ploidy 1–8 × alleles 1–8 | 24,301 | The odd ploidies, the haploid case, widths past the candidate cap. Stops at 8 alleles for cost: the full 8 × 16 grid is 2,042,958 genotypes, eighty-four times as many. |
| `..._past_the_cache_bounds` | ploidy {9, 10} × alleles 1–4, ploidy {2, 3} × alleles {17, 18} | 3,083 | `build`'s **uncached branch** — the one every polyploid or wide locus takes, and which the other two grids never reach. |
| `a_shape_at_the_cache_bound_is_shared_and_one_past_it_is_not` | four shapes on the bound, one past it | — | The cache identity the plan asks for, taken at the boundary and in both directions. |

All three totals were re-derived independently of the test before being written into it:
`sum(comb(a + p − 1, p))` over each grid gives 308, 24,301 and 3,083.

## 5. Validation

Run in the dev container after the review's fixes:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 53 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |

The three aggregate gates red on `main` are unchanged and untouched by this step.

## 6. What the parity test actually catches — measured, not argued

Every one of these was run against the finished code by re-applying the mutation and reading the
failure:

| mutation to the port | caught by |
|---|---|
| `log_factorial` summed downwards — a four-ulp difference | `ploidy 4 over 2 alleles: row 1 [3, 1] coefficient — got 1.3862943611198908, production 1.3862943611198904`. A 1e-12 tolerance would have passed it. |
| the enumeration's outer loop reversed | `ploidy 1 over 2 alleles: allele counts, or the order they are in` |
| the homozygous lookup never fires | `ploidy 2 over 2 alleles: row 2 [0, 2] homozygous lookup` |
| `log_factorial` capped at 8 — silent above the cache bound | the past-the-bounds grid, **and** the strengthened C1 test |
| one extra coefficient row appended | `ploidy 1 over 1 alleles: one coefficient per genotype` |
| one extra homozygous entry appended | `ploidy 1 over 1 alleles: one homozygous entry per genotype` |
| the cache bound turned from `>` to `>=` | `ploidy 8 over 16 alleles is inside the cache bounds` |

## 7. And the port is exact

Over ploidy 1–8 × alleles 1–16 — **128 shapes, 2,042,958 genotype rows** — one of the review's
agents compared ng's table against `genotype_order` and a transcription of production's fold:
enumeration order position for position, every log coefficient **bit-identical** (`to_bits()`
equality, no tolerance), the homozygous lookup entry for entry. Zero differences. The committed
grids are the subset that runs in the suite; that measurement is the reason nothing wider needs to.

## 8. Tradeoffs and follow-ups

- **The plan's step C2 text still instructs the next reader to call `GenotypeShape` directly.** It
  cannot be done, and the plan is a design document this loop does not edit. Raised to the owner;
  the substitution and its justification are recorded in the module's own doc comment, in the
  review, and here.
- **`genotype_order` narrows its allele count with `n_alleles as u8`**
  (`src/var_calling/per_group_merger.rs`), so the oracle silently stops being production above 255
  alleles, where the port reaches 65,536. Frozen production; documented in the module header so a
  future wider grid does not read the resulting failure as the port's fault.
- **The grids stop at 18 alleles.** The 8 × 16 grid is eighty-four times the genotypes of the
  committed one for the same four closed forms; the review measured it at seconds rather than
  hundredths in the debug profile the suite runs in.
