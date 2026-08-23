# Code Review: ng_calling_c2
**Date:** 2026-08-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step C2 of the calling-foundations plan — the genotype table against production, value for value
**Status:** Request-changes (one Blocker, one Major; both applied — see the fix report)

---

## 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step C2 on branch
  `ng-calling-foundations`, applied to `cb36127e` from `tmp/c2.patch`.
- **In-scope files:**
  [genotype_table_parity.rs](../../../../src/ng/calling/genotype_table_parity.rs) (new, 212 lines
  at review time); [mod.rs](../../../../src/ng/calling/mod.rs), the 2 added lines.
- **Deliberately out of scope:** `genotype_table.rs` (committed at C1 and reviewed there — it is
  the *subject* here, mutated freely to test the oracle but not reviewed); all of
  `src/var_calling/` (production, frozen — read closely as the transcription's source).
- **Categories dispatched**, three agents, each in its own worktree detached at `cb36127e`:

| agent | categories | why |
|---|---|---|
| 1 | `reliability` | the module *is* an oracle; "what wrong port would it let through" is the whole review |
| 2 | `refactor_safety`, `extras`, the diff's numbers | it is a transcription of frozen code, and every claim it makes about that code is checkable |
| 3 | `naming`, `module_structure`, `idiomatic`, `smells`, `errors` | a new file, a new import direction, a 45-line argument in a doc comment |

## 2. Verdict

**Request-changes.** One Blocker, one Major, ten Minors. The oracle itself is sound and the port it
checks is exact; the Blocker is about where the grids stop.

## 3. Execution status

Run by the orchestrator in the dev container on the reviewed diff and handed to all three agents:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.86s` |
| `cargo test --lib ng::calling` | 0 | `52 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --lib ng::calling::genotype_table_parity` | 0 | `3 passed; 0 failed` |

The three aggregate gates red on `main` were excluded by instruction and are untouched.
Findings labelled "Needs verification": **none**.

## 4. Open questions and assumptions

1. **The plan's step C2 says the test calls `GenotypeShape` "directly", and the compiler refuses.**
   Confirmed three ways: by the orchestrator's own compile, and by two agents independently, all
   getting ``error[E0603]: module `shape` is private`` — `posterior_engine.rs:54` declares
   `mod shape;` privately, and `:61` imports `GenotypeShape`/`shape_for` privately too, so there is
   no re-export route. Widening it would edit a frozen tree. **The plan text is not this loop's to
   change**; it is recorded here and raised to the owner. Affects nothing in the code.

## 5. Top 3 priorities

1. **B1** — both grids stop at ploidy 8, exactly where the cache stops, so no coefficient on the
   uncached branch is ever compared with production. A `log_factorial` capped at 8 passes the
   whole suite while understating every heterozygous coefficient at ploidy 9 by `ln 9`.
2. **M1** — the allele axis stops at 8 too, while the cache keeps 16; two defects confined to
   wider tables survive this module.
3. **Mi1** — the coefficient and homozygous tables are read row by row with no length check, so a
   port carrying an extra trailing row passes all three parity tests.

## 6. Findings

### Blocker

**B1: [genotype_table_parity.rs:190](../../../../src/ng/calling/genotype_table_parity.rs) — nothing
compares the uncached branch with production, and a wrong coefficient there is silent**
**Categories:** reliability. **Confidence:** High.

Both grids run `for copies in 1..=8`, which is `MAX_CACHED_PLOIDY`. Above it `build` takes its
uncached branch, documented as supported — the constants' own doc names "a dodecaploid". **Measured:**
capping the port's `log_factorial` at 8 leaves the entire `ng::calling` suite green, 53 passed. The
mutant is not a no-op and not a rounding difference: at ploidy 9 the homozygous rows stay exact
(`ln 9! − ln 9!` is zero however the terms are capped) while **every heterozygous row is understated
by exactly `ln 9 = 2.197` nats** — `[8, 1]` becomes 0 instead of 2.197. That is a genotype prior
tilted toward homozygotes by a factor of nine at every polyploid locus, with nothing crashing, in
the range `CLAUDE.md` commits the caller to.

**There is no second line of defence, and the reason is worth recording.** C1's
`a_shape_past_the_cache_bounds_holds_the_right_values` is the suite's only coefficient assertion
above ploidy 8, and it checks row 0 — `[9, 0]`, the homozygous row, whose coefficient is zero under
the right formula *and* under a `log_factorial` capped at 8, at 4, or at 2. **The one row it checked
was the one row the defect cannot touch.**

**Fix:** a third grid over the uncached branch — ploidy 9 and 10 at 1–4 alleles, ploidy 2 and 3 at
17 and 18 alleles; 3,083 genotypes, widest shape 1,140. Verified green on the clean port and red
under the mutation with production's own number in the message.

### Major

**M1: [genotype_table_parity.rs:194](../../../../src/ng/calling/genotype_table_parity.rs) — the
allele axis stops at 8 while the cache keeps 16**
**Categories:** reliability. **Confidence:** High.

Two mutants confined to wider tables passed this module: capping both the recursion and the
closed-form count at 8 alleles (so the two still agree and C1's internal consistency assert never
fires), and mis-naming the homozygous allele only for the top allele of a table wider than 8. Both
were caught only by C1's own tests. Lower than the Blocker because C1 does catch them — but the
parity module is the artefact whose stated job is agreement with production. The same third grid
closes it: one test kills four of the six survivors.

### Minor

- **Mi1** `:141`, `:155` — the coefficient and homozygous tables are indexed row by row with no
  length assertion, so a port with an extra trailing row passes all three parity tests
  (**measured**: `3 passed`), while the same build fails 6 and 5 of C1's tests respectively. The
  allele-count table is already compared as a whole slice; the other two are not. *(refactor_safety)*
- **Mi2** `:161` — the grid totals are summed from the **oracle's** count, so they certify the grid
  was walked rather than that anything matched: delete all four comparisons and both totals still
  reach 308 and 24,301. *(reliability)*
- **Mi3** `:119` — nothing asks the table what shape it thinks it is. **Measured:** widening the
  stored `allele_count` to `allele_count.max(2)` leaves all three parity tests green. Mitigated by
  C1's grid test. *(reliability)*
- **Mi4** `:206` — the cache test picks an interior point and only the positive direction. Turning
  the bound from `>` to `>=` — which stops the cache at ploidy 8, at 16 alleles, and at the corner
  — leaves it green. It is also a near-duplicate of two C1 tests. *(reliability, smells — convergent)*
- **Mi5** `:30-43` — the doc's argument for the substitution is wrong in two places: **four** things
  are transcriptions, not three (the fold is one), and the *order* is named as the part that could
  drift silently when it is the one part that **cannot**, precisely because it runs production's
  live code. *(refactor_safety)*
- **Mi6** `:30` — the oracle stops being production above 255 alleles: `genotype_order` iterates
  `min_allele..(n_alleles as u8)`, so 256 alleles returns no genotypes at all, where the port
  reaches 65,536. Not a live defect — the count assertion would fire first — but the message would
  read as the *port* being wrong. *(reliability)*
- **Mi7** `:72-76` — the `#[allow(clippy::cast_precision_loss)]` suppresses a lint that cannot fire
  and its `reason` is checkably false. **Measured twice:** with the attribute deleted, clippy under
  `-D warnings` is clean; and forcing pedantic on produces 2,175 warnings crate-wide, **zero** in
  this file. The reason claims `i as f64` is what makes the transcription bit-exact — replacing it
  with `f64::from(i)` leaves all 24,609 coefficient comparisons passing. *(smells)*
- **Mi8** `:190` — `..._over_every_shape_from_haploid_to_octoploid` claims coverage the test does
  not have; it walks 64 of the cache's 128 shapes. *(naming)*
- **Mi9** `:170` — `DEFAULT_MAX_CANDIDATE_ALLELES = 6` is called "the cap the caller ships with",
  and no constant of that name exists in `src/`. The value is right, inherited from production's
  live `DEFAULT_MAX_ALLELES_PER_RECORD`. *(extras)*
- **Mi10** `src/ng/mod.rs` and `calling/mod.rs` — this is the first compile-time ng→production edge
  (`grep -rn crate::var_calling src/ng/` gives four hits, three of them rustdoc links), and neither
  the freeze paragraph a reader consults nor `calling/mod.rs`'s "What is here" list mentions the
  test-oracle exception or the new module. *(module_structure)*

### Nits

`where_` dodges the keyword with a trailing underscore where *shape* is the crate's own word for
the pair; the oracle helpers take `ploidy: u8` in a file where `ploidy` is a newtype; "cheap" is a
placeholder for a number (0.06 s against 6.63 s); "three of its four quantities" never names the
four; the quoted E0603 message uses straight quotes where rustc emits backticks; `genotype_order`
is built twice per shape; "verbatim" is not quite right for the homozygous transcription
(production narrows to `u8`, the oracle keeps `usize`); the `#[cfg(test)] mod` sits mid-`pub mod`
block where `src/ng/mod.rs` puts it first.

## 7. Out of scope observations

- **`per_group_merger.rs:558`** — `genotype_order` takes `n_alleles: usize` and narrows it with
  `n_alleles as u8`, so 256 alleles silently produces the wrong enumeration. Frozen production,
  far outside any grid here; noted because this diff is its first ng caller.
- **The plan's step C2 text** — see open question 1.

## 8. Missing tests to add now

Three, all applied: the past-the-cache-bounds grid (B1/M1), the two length assertions and the two
shape assertions inside `compare_against_production` (Mi1, Mi3), and the cache test moved to the
boundary with its negative direction (Mi4).

## 9. What's good

- **The `to_bits()` comparison earns its keep.** Reversing only the port's summation order fails at
  ploidy 4 over 2 alleles with `1.3862943611198908` against `1.3862943611198904` — four units in
  the last place, which any tolerance of 1e-12 would have passed.
- **The transcription is faithful in all four functions**, checked statement for statement against
  `shape.rs`. The two textual differences — `u32::from(copies)` for `ploidy as u32`, and
  `(i as f64)` against the port's `f64::from(i)` — were both proven value-identical, the second by
  walking all 4,294,967,296 `u32` values and comparing bit patterns: zero mismatches.
- **The grid totals are load-bearing, not decorative.** Of the 4,095 proper subsets of the plan's
  twelve shapes, **zero** sum to 308; of all subsets of the sixty-four, **exactly one** sums to
  24,301 — the whole grid. A silently narrowed loop fails on its total.
- **The oracle is genuinely independent where it matters.** The port's enumeration is a recursion,
  production's is generate-then-sort; a shared enumeration bug would have to be written twice in
  two different shapes, and could not arrive by copying.
- **`genotype_order` really is the sole source of production's enumeration** —
  `collect_non_decreasing` is its only enumeration recursion, and every consumer, VCF writer
  included, reaches it through that one function.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling`

## Author response convention

Address each finding by its identifier (B1, M1, Mi7, …) with `fixed in <commit>` / `disputed
because …` / `deferred to …`. Answer the open question in section 4 first.
