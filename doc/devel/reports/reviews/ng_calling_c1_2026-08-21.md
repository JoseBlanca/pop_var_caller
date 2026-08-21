# Code Review: ng_calling_c1
**Date:** 2026-08-21
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step C1 of the calling-foundations plan — `GenotypeTable`, `GenotypeIdx` and the flat views
**Status:** Request-changes (one Blocker and five Majors; all applied, see the fix report)

---

## 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step C1 on branch
  `ng-calling-foundations`, applied to `af458f44` from `tmp/c1.patch`.
- **In-scope files:** [genotype_table.rs](../../../../src/ng/calling/genotype_table.rs) (new, 872
  lines at review time, of which the test module was 409);
  [mod.rs](../../../../src/ng/calling/mod.rs), the 8 added lines.
- **Deliberately out of scope:** the rest of `calling/mod.rs` (reviewed at B1–B3);
  `src/ng/types.rs`; all of `src/var_calling/` (production, frozen — read as the port's source);
  step C2's parity test.
- **Categories dispatched**, five agents, each in its own git worktree, each detached at
  `af458f44` with the patch applied:

| agent | categories | why |
|---|---|---|
| 1 | `reliability` | always; the module is arithmetic with a silent failure mode |
| 2 | `naming`, `module_structure` | a new module and a new public surface |
| 3 | `idiomatic`, `smells`, `defaults`, `errors` | three panicking refusals, two inherited bounds |
| 4 | `unsafe_concurrency` | a thread-local `RefCell` cache handing out `Arc`s |
| 5 | `refactor_safety`, `extras` + the diff's own numbers | it is a port with a stable-output contract |

## 2. Verdict

**Request-changes.** One Blocker, five Majors, thirteen Minors. The port itself is correct — proved
below — and every finding is about what the tests do not pin, what the prose claims that the code
does not do, or an arithmetic guard that stops one line short.

## 3. Execution status

Run by the orchestrator in the dev container on the reviewed diff, and handed to every agent so
none re-ran them:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | `Finished dev profile … in 2.78s` |
| `cargo test --lib ng::calling` | 0 | `43 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |

The three aggregate gates red on `main` were excluded by instruction and are untouched by this
step: `cargo clippy --all-targets` (18 errors in two benches and one example), `cargo test
--all-targets` (`benches/psp_writer_perf.rs:386`), `cargo doc --no-deps --lib` (17 unresolved
intra-doc links).

Findings labelled "Needs verification": **none**. Every finding below was reproduced by running
something.

## 4. Open questions and assumptions

1. **`GenotypeIdx` or `GenotypeIndex`?** The crate spells index newtypes `Id` (`AlleleId`,
   `ContigId`, `ReadGroupId`) or `Index` (`WindowIndex`); this is its only `Idx`. Affects Mi4.
   Not resolved here: the name is written into the plan and all three calling architecture
   documents, so changing it is the owner's call, not the reviewer's.
2. **Does the calling loop want a table before its allele set is final?** Affects Mi9 (a
   `for_locus(&CandidateAlleles)` constructor). If it does, the constructor is wrong and the
   zero-allele assertion has to stay as the only guard.

## 5. Top 3 priorities

1. **B1** — the VCF genotype order is unpinned at ploidy ≥ 3 with ≥ 3 alleles; a row swap there
   passes all 22 tests and changes which genotype `PL[0]` names.
2. **M1** — one plausible slip in the cache's slot index hands a locus another shape's genotype
   table, silently, and all 22 tests still pass.
3. **M2** — `n_genotypes * allele_count` is the one unchecked multiplication, on a shape the
   `# Panics` documentation says is accepted.

## 6. Findings

### Blocker

**B1: [genotype_table.rs:403-418](../../../../src/ng/calling/genotype_table.rs) — the enumeration
order is unpinned at ploidy ≥ 3 with ≥ 3 alleles**
**Categories:** reliability. **Confidence:** High.

The module's own header names the enumeration order as the contract whose breach is silent. Two
tests assert an order: ploidy 2 over 3 alleles, and ploidy 4 over 2 alleles. Every other test is
invariant under a permutation of rows, because the three tables are derived from the allele-count
table and permute together. **Measured:** a mutant swapping rows 0 and 1 whenever ploidy ≥ 3 *and*
allele count ≥ 3 leaves all 22 tests green while turning ploidy 3 over 3 alleles from
`[3,0,0, 2,1,0, 1,2,0]` into `[2,1,0, 3,0,0, 1,2,0]` — `PL[0]` stops being homozygous-reference.
Ploidy ≥ 3 with ≥ 3 candidates is inside the caller's stated range.

**Fix:** assert the ordering rule as a law over a grid rather than at two shapes — consecutive rows,
read as sorted allele lists, strictly increasing compared from the highest copy downwards, which is
the comparator production's `genotype_order` sorts with.

### Major

**M1: [genotype_table.rs:156](../../../../src/ng/calling/genotype_table.rs) — a one-token slip in
the cache's slot index gives a locus another shape's table**
**Categories:** refactor_safety. **Confidence:** High.

`CACHE_SLOTS` is written from both bounds and the slot index from one of them; the two constants sit
side by side. Substituting the neighbour — `p * (MAX_CACHED_PLOIDY + 1) + allele_count` — keeps every
index inside the 153-slot array, so nothing panics, but ploidy 1 over 16 alleles and ploidy 2 over 7
alleles land in the same slot. **Measured:** all 22 tests pass under that substitution; only a
full-grid comparison against production fails. The four cache tests touch seven shapes, no two of
which collide. C2 will not catch it either — C2 compares one shape at a time, and a collision needs
two shapes in one thread's cache.

**Fix:** a test that walks the whole cached grid and asks each table which shape it is (the table
already stores both), plus a `debug_assert` on the cache hit.

**M2: [genotype_table.rs:383](../../../../src/ng/calling/genotype_table.rs) — the row buffer's
capacity is the one unchecked multiplication, and the docs say it is checked**
**Categories:** errors, reliability, refactor_safety (convergent — three agents). **Confidence:** High.

`build`'s `# Panics` section promised that the genotype-count check "turns an arithmetic overflow
deep in a capacity calculation into a refusal that names the shape". It covered the genotype count
and not the table. **Measured:** `build(Ploidy(4), 65_536)` passes all three refusals —
`count_genotypes` returns `Some(768_684_707_117_285_376)`, which fits a `usize` — and then
`Vec::with_capacity(n_genotypes * allele_count)` overflows. In debug that is `attempt to multiply
with overflow`; in release the multiply wraps and the `Vec` grows until allocation fails.

**Fix:** fold the table size into the same refusal, so the message names the shape.

**M3: [genotype_table.rs:17-19](../../../../src/ng/calling/genotype_table.rs) and `:380-381` — two
doc comments claim a test that this commit does not contain**
**Categories:** reliability, smells, naming (convergent — three agents). **Confidence:** High.

Both say the test module pins the enumeration order against production's `genotype_order`. `grep`
over the file returns only doc-comment mentions; the parity test is step C2, marked "**Own commit,
do not bundle**". The claim is true of the plan and false of the file — and it is how B1 stays
invisible, since a reader takes the contract to be already guarded.

**M4: [genotype_table.rs:813-824](../../../../src/ng/calling/genotype_table.rs) — no test asserts a
value for any shape outside the cache bounds**
**Categories:** reliability. **Confidence:** High.

The only assertion there is `assert_eq!(first, second)` — two runs of one deterministic pure
function, which no wrong-value implementation can violate. **Measured:** a mutant that zeroes the
coefficients and blanks the homozygous lookup *on the uncached branch only* leaves all 22 tests
green. Ploidy above 8 and more than 16 alleles are both inside the caller's stated range, and a
blanked homozygous lookup is a silently wrong inbreeding prior.

**M5: [genotype_table.rs:373](../../../../src/ng/calling/genotype_table.rs) — the running binomial's
`checked_mul` has no test**
**Categories:** reliability. **Confidence:** High.

The shipped refusal fixture, ploidy 40 over 64 alleles, exits through the final `usize::try_from`:
its largest intermediate, C(103, 40) ≈ 6.1e28, fits a `u128`. **Measured:** removing the
`checked_mul` leaves all 22 tests green, while `count_genotypes(255, 65_536)` goes from `None` to an
overflow panic. In release the running binomial would wrap and return a small, plausible count.

### Minor

- **Mi1** `mod.rs:35` — `pub mod` with no re-export puts calling's shared vocabulary at two
  depths; `calling/mod.rs` is the only ng module with a `pub mod` and no matching `pub use`. The
  `dead_code` half of the stated reason is refuted by experiment: `mod` + `pub use` compiles clean
  under `-D warnings`. *(naming/module_structure)*
- **Mi2** `genotype_table.rs:117` — `n_genotypes` the field and `genotype_count` the function are
  one quantity under two names, and `n_` runs against the crate's own convention (71 `allele_count`
  against 29 `n_genotypes`, 27 of them in this file). *(naming)*
- **Mi3** `:259`/`:292` — `homozygous_allele_for()` and `homozygous_allele_of()` differ by one
  preposition while returning different things; the file's other two accessor pairs distinguish
  column from row by number. *(naming)*
- **Mi4** `:74` — `GenotypeIdx` is the crate's only `Idx`. **Not applied** — see open question 1.
- **Mi5** `:270-292` — three row lookups each re-derive the row and decide bounds separately, and
  the view the loop actually consumes has no lookups at all. *(smells, idiomatic)*
- **Mi6** `:38`/`:40` — the two cache bounds are private while `build`'s public doc links them, and
  the stated cost ("the slots") is 1,224 bytes where a filled slot measures 37,263,864 bytes, never
  evicted. *(defaults, unsafe_concurrency — convergent)*
- **Mi7** `:158` — the `RefCell` borrow spans the whole build. Re-entry is unreachable today
  (call tree traced; no call site outside this file's tests), but the release profile aborts on
  panic, so a future re-entrant call would kill the run. Building outside the borrow costs nothing
  measurable: 3.4 ns per cache hit either way. *(unsafe_concurrency, reliability — convergent)*
- **Mi8** `:214-223` — `view()` copies six fields by hand, so a seventh field on `GenotypeTable`
  reaches consumers as silence. **Measured:** adding a field and leaving `view()` untouched compiles
  with no error and no warning. *(refactor_safety)*
- **Mi9** `:151` — `build` takes a bare `usize` where the loop will always hold a
  `CandidateAlleles`. **Not applied** — see open question 2.
- **Mi10** `:141-143` — the width refusal's stated reason is off by one against its own boundary:
  an `AlleleId` names 65,536 alleles, so 65,536 is what it *can* name and 65,537 is one more.
  *(refactor_safety)*
- **Mi11** `:457` — `u16::try_from(allele).expect(…)` relies on an assertion 280 lines away, with
  the bound written twice in two forms and no `// PANIC-FREE:` comment. *(errors)*
- **Mi12** `:181-198` — nothing ties the closed-form genotype count to the number of rows the
  recursion emits; a disagreement would silently shorten two of the three tables. *(reliability)*
- **Mi13** `:559`–`:675` — the ploidy loops stop at 6 though the cache is sized for 8, with no
  recorded reason. **Measured:** widening is safe — worst coefficient error 1.4e-14 at ploidy 30
  against a 1e-12 tolerance. *(reliability)*

### Nits

`GenotypeIdx::get` is the only `pub` item without a doc comment; `only` is half a name for a binding
that lives across a twelve-line loop; the doc says `n_alleles` where the code says `allele_count`
(the only `n_alleles` in `src/ng`); `count_genotypes` and `homozygous_allele` are computations named
as nouns; the `lgamma` remark is a performance claim with no measurement; `GenotypeTableView` is
`Copy` but its accessors take `&self`; no `#[must_use]` on `build` (matches its neighbours); no doc
example (matches `src/ng/calling/`, not the wider tree); no `criterion` bench for a per-locus call.

## 7. Out of scope observations

- **`doc/devel/ng/arch/module_layout.md:93-95`** still lists `GenotypeTable`, `GenotypeIdx`,
  `AlleleId` and `Genotype` as contents of `calling/mod.rs`. Three of the four have moved — two to
  `types.rs` at Milestone A, two to this file — and the block has no `genotype_table.rs` line.
  `calling_em_loop.md:36-49` is the correct tree. A design-document fix, outside this loop's remit.
- **Two `dev.sh` containers outlived the host process that started them** during one agent's
  mutation runs, concurrently rewriting the file under test. That agent stopped them with
  `container stop`, rebuilt from the patch and re-ran everything affected. Worth knowing that
  stopping a `dev.sh` invocation does not stop its container.

## 8. Missing tests to add now

Five, all with bodies supplied by the agents that found them and all applied:
`every_consecutive_pair_of_rows_is_in_vcf_order` (B1),
`every_cached_shape_gets_back_a_table_of_its_own_shape` (M1),
`a_shape_past_the_cache_bounds_holds_the_right_values` (M4),
`a_shape_that_overflows_the_running_binomial_is_refused` (M5),
`a_shape_whose_flat_table_outgrows_a_usize_is_refused` (M2). A sixth,
`build_returns_the_same_values_from_the_cache_and_from_a_fresh_build`, was proposed to stop the two
branches drifting and is also applied.

Proposed and **not** applied: a `proptest` over the whole `(ploidy, allele count)` domain — the four
laws it would assert are each already asserted over a grid, and the grid is what makes a failure
reproducible; a 64-thread agreement test — the concurrency agent measured it (64 of 64 threads agree)
and nothing in the type can make it fail, since the cache is per-thread and the value a pure
function of its key.

## 9. What's good

- **The three refusals in `build_uncached` are the right shape and the right severity.** Production
  does not check a zero-allele shape and panics anyway, inside `chunks_exact(0)`, with a message
  about chunk sizes; this one names the domain fact instead.
- **Building the VCF order directly, rather than enumerating and sorting, leaves no comparator to
  get wrong** — and the two mutations that reorder the enumeration are both caught by the small
  hand-written fixtures, which are carrying more weight than their size suggests.
- **`log_factorial` keeps production's summation order on purpose, and it pays off**: all 2,042,958
  coefficients over ploidy 1–8 × alleles 1–16 are bit-identical to production, `to_bits()` equality,
  no tolerance.
- **The homozygous lookup is one table with one producer**, which is what
  `calling_priors.md` §3.2 asks for and what gives the above-diploidy spec one place to change.
- **`GenotypeTableView` has private fields and no public constructor**, so "a view always describes
  a table that exists" is true rather than aspirational.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling`

## Appendix: the port, measured

The one thing this review can settle that C2 was written to settle: **the port reproduces
production's enumeration and coefficients exactly.** Production enumerates every non-decreasing
allele tuple and sorts them with a comparator walking each tuple backwards
([per_group_merger.rs:522](../../../../src/var_calling/per_group_merger.rs)); the port recurses on
each genotype's highest allele. Over **ploidy 1–8 × allele count 1–16 — 128 shapes, 2,042,958
genotype rows** — three comparisons against `genotype_order` and a transcription of production's
own folding step, all green: the allele-count table position for position, every log coefficient
bit-identical, and the homozygous lookup entry for entry.

**One thing C2's author must know: `GenotypeShape` and `shape_for` are not reachable from `src/ng/`.**
`posterior_engine.rs:54` declares `mod shape;` privately, so the plan's "called directly" is not
possible. Confirmed twice — once by an agent's probe, once by the orchestrator's own compile, which
fails with `error[E0603]: module 'shape' is private`. `genotype_order` **is** reachable
(`pub fn` in `pub mod per_group_merger`), and is what this appendix's oracle was built on.

## Author response convention

Address each finding by its identifier (B1, M1, Mi6, …) with `fixed in <commit>` / `disputed
because …` / `deferred to …`. Answer the open questions in section 4 first.
