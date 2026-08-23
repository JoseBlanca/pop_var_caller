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

**The oracle is production's own artefact, and reaching it cost production one line.** The plan
says this test calls `GenotypeShape` directly. `posterior_engine.rs` declared `mod shape;`
privately, so neither `GenotypeShape` nor `shape_for` could be named from `src/ng/` whatever their
own `pub(crate)` visibility, and `:61` imports them privately too, so there was no re-export route:
a call from ng gave ``error[E0603]: module `shape` is private``, confirmed three times
independently. **The owner authorised widening that declaration** (2026-08-21), so it is now
`pub(crate) mod shape;` — no behaviour moved, nothing re-exported, and the dependency runs one way.

That is worth stating as the rule it sets rather than as a one-off: **ng may widen a production
item's visibility so a parity test can see it, and may change nothing else.** Anything that would
alter what production computes is still a copy-into-ng.

What it buys is the difference between a strong test and a weak one. The first draft of this module
built its oracle from `genotype_order` plus four verbatim transcriptions of `shape.rs` — its fold
and its three formulas — and its review filed the consequence: if production ever changed any of
those four, the test would keep passing while ng's table quietly stopped matching production's.
`shape_for` closes that gap outright. Every field compared below is the field the shipping
posterior engine reads, and nothing is re-derived.

**The grids are wider than the plan's.** The plan names ploidy 2 and 4 over alleles 1–6; that grid
is here as its own test, and two more sit beside it (see §4). The plan's grid is a strict subset of
the second, so its contribution is its total; it is kept because the plan names it.

**`src/ng/mod.rs` gained a paragraph**, naming both test oracles and the one production edit, since
the freeze rule a reader consults said nothing about either.

**One C1 test was strengthened here**, which is a change to a committed file and is recorded as
such: see §3.

## 3. Changes made

Five files:

- **`src/var_calling/posterior_engine.rs`** — `mod shape;` → `pub(crate) mod shape;`, with a
  comment saying who reads it and why. The only edit ng has made to a frozen tree.
- **`src/ng/calling/genotype_table_parity.rs`** (new, 239 lines) —
  `compare_against_production` (five comparisons per shape, all against `shape_for`'s result) and
  four tests.
- **`src/ng/calling/mod.rs`** — `#[cfg(test)] mod genotype_table_parity;` and a clause in the
  module doc saying what it is and why it is its own file.
- **`src/ng/mod.rs`** — the test-oracle exception to the production freeze.
- **`src/ng/calling/genotype_table.rs`** — one C1 test strengthened.
  `a_shape_past_the_cache_bounds_holds_the_right_values` asserted the coefficient of row 0 at
  ploidy 9. Row 0 is `[9, 0]`, whose coefficient is `ln 9! − ln 9!` — zero under the right formula
  and zero under a `log_factorial` capped at 8, at 4, or at 2. **The one row it checked was the one
  row the defect cannot touch.** It now checks row 1, `[8, 1]`, worth `ln 9`.

## 4. Tests added

Four, in `src/ng/calling/genotype_table_parity.rs`. Every shape goes through five comparisons
against the `GenotypeShape` `shape_for` returns: the table's own declared ploidy and width, the
genotype count, the whole allele-count table as one slice, the coefficients as one slice of bit
patterns, and the homozygous lookup as one slice. Comparing whole slices rather than row by row is
what makes a table with the wrong number of rows fail rather than have its extra rows go unread.

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
| one extra coefficient row appended | all four grids: `log multinomial coefficients as bit patterns, or how many of them` |
| the cache bound turned from `>` to `>=` | `ploidy 8 over 16 alleles is inside the cache bounds` |

## 7. And the port is exact

Over ploidy 1–8 × alleles 1–16 — **128 shapes, 2,042,958 genotype rows** — one of the review's
agents compared ng's table against production over the whole cached region:
enumeration order position for position, every log coefficient **bit-identical** (`to_bits()`
equality, no tolerance), the homozygous lookup entry for entry. Zero differences. The committed
grids are the subset that runs in the suite; that measurement is the reason nothing wider needs to.

## 8. Tradeoffs and follow-ups

- **`genotype_order` narrows its allele count with `n_alleles as u8`**
  (`src/var_calling/per_group_merger.rs`), so the oracle silently stops being production above 255
  alleles, where the port reaches 65,536. Frozen production; documented in the module header so a
  future wider grid does not read the resulting failure as the port's fault.
- **The grids stop at 18 alleles.** The 8 × 16 grid is eighty-four times the genotypes of the
  committed one for the same four closed forms; the review measured it at seconds rather than
  hundredths in the debug profile the suite runs in.
