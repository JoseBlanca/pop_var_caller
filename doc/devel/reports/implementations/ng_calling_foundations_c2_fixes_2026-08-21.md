# ng calling foundations — C2: what the review changed

*Fix-application report, 2026-08-21. Branch `ng-calling-foundations`. Applies
[`ng_calling_c2_2026-08-21.md`](../reviews/ng_calling_c2_2026-08-21.md) to the step-C2 diff
described in [`ng_calling_foundations_c2_2026-08-21.md`](ng_calling_foundations_c2_2026-08-21.md).*

## 1. What was applied

Everything the review filed. One Blocker, one Major, ten Minors and every actionable nit.

**Four of them were provable only by mutating the subject, and each fix was re-verified here by
re-running its mutation** — not by trusting the reviewer's run. Each was applied to the fixed
tree, the suite run, the mutation reverted:

| what was mutated in `genotype_table.rs` | before the fix | after the fix |
|---|---|---|
| `log_factorial` capped at 8 | the whole `ng::calling` suite green, 53 passed | `the_table_matches_production_past_the_cache_bounds` FAILS **and** the strengthened C1 test FAILS |
| one extra coefficient row appended | 3 parity tests pass | `ploidy 1 over 1 alleles: one coefficient per genotype` |
| one extra homozygous entry appended | 3 parity tests pass | `ploidy 1 over 1 alleles: one homozygous entry per genotype` |
| the cache bound turned from `>` to `>=` | the cache test passes | `ploidy 8 over 16 alleles is inside the cache bounds` |

### The Blocker and the Major

**B1 and M1 — the uncached branch is now compared with production.** A third grid covers
`build`'s uncached branch: ploidy 9 and 10 at 1–4 alleles, and ploidy 2 and 3 at 17 and 18
alleles — 3,083 genotypes, widest shape 1,140, the total re-derived here before it was written
into the assertion. The two grids that existed both stopped at 8 on each axis, which is exactly
where the cache stops, so the branch every polyploid or wide locus takes was never compared with
production at all. What that let through: a `log_factorial` capped at 8 leaves the homozygous rows
exact and understates every heterozygous coefficient at ploidy 9 by `ln 9 = 2.197` nats.

### The Minors

- **Mi1** — the two parallel tables' lengths are asserted before they are read row by row, so a
  port with an extra trailing row now fails here. The per-row loops stay, because their message
  names the genotype; the length assertions are what make the loops sufficient.
- **Mi2** — `compare_against_production` returns the **table's** genotype count rather than the
  oracle's, so the grid totals are evidence about the subject. The two are asserted equal one line
  earlier.
- **Mi3** — the comparison now starts by asking the table its own ploidy and width, since
  everything below slices rows by the width the table declares.
- **Mi4** — the cache test moved to the boundary and gained its negative direction: the three
  shapes on the bound plus one past it, which is where an off-by-one in either comparison shows.
- **Mi5** — the doc's argument is corrected. **Four** things are transcriptions, not three — the
  fold is one — and the *order* is now named as the part that **cannot** drift silently, because
  it is the part that runs production's live code. The previous text said the opposite.
- **Mi6** — the header records that the oracle stops being production above 255 alleles, since
  `genotype_order` narrows its allele count to `u8`. Not reachable from any grid here; recorded so
  a future wider grid does not read the resulting failure as the port's fault.
- **Mi7** — the `#[allow(clippy::cast_precision_loss)]` is deleted. The lint cannot fire (clippy
  does not flag `u32 as f64`, and the lint is not enabled), and its stated reason — that the `as`
  cast is what makes the transcription bit-exact — is false: `f64::from(i)` gives the same bits for
  every `u32`. The cast stays, for transcription fidelity, and the doc now says that instead.
- **Mi8** — the overclaiming test name is now
  `..._from_haploid_to_octoploid_up_to_eight_alleles`, and its doc gives the reason for stopping at
  8 alleles in genotypes: 24,301 against the full grid's 2,042,958, eighty-four times as many.
- **Mi9** — `DEFAULT_MAX_CANDIDATE_ALLELES = 6` is no longer described as a cap the caller ships
  with. It is not a constant in this tree; it is inherited from production's
  `DEFAULT_MAX_ALLELES_PER_RECORD` and recorded in the architecture, and the doc now says so.
- **Mi10** — `src/ng/mod.rs`'s freeze paragraph names the test-oracle exception and both instances
  of it, and `calling/mod.rs`'s "What is here" list names the new module and why it is its own
  file.

### Nits

`where_` is `shape`, the crate's own word for the `(ploidy, allele count)` pair; the oracle helpers
take `copies: u8` rather than `ploidy: u8` in a file where `ploidy` is a newtype; the four
quantities are named where the header first counts them; the E0603 message is quoted with rustc's
backticks; "cheap" is replaced by the genotype counts; the homozygous transcription's one
difference from production (`usize` where production narrows to `u8`) is stated rather than
covered by "verbatim"; the `#[cfg(test)] mod` sits ahead of the `pub mod` block, matching
`src/ng/mod.rs`.

## 2. A C1 test was strengthened, and it is the finding worth remembering

`a_shape_past_the_cache_bounds_holds_the_right_values` was written at C1 to check values on the
uncached branch, because comparing two runs of a pure function cannot catch a wrong one. It
asserted the coefficient of **row 0** at ploidy 9 over two alleles. Row 0 is `[9, 0]`, whose
coefficient is `ln 9! − ln 9!` — zero under the right formula, and zero under a `log_factorial`
capped at 8, at 4, or at 2. **The one row it checked was the one row the defect cannot touch.**

It now checks row 1, `[8, 1]`, worth `ln 9`. The assertion is a tolerance rather than exact
equality, and the comment says why: the coefficient is *summed* as `ln 9! − ln 8!`, so it agrees
with `9.0_f64.ln()` in value and not in every bit — the exact bits are pinned against production
next door. The first attempt at this fix asserted exact equality and failed at
`2.19722457733622` against `2.1972245773362196`, which is how that was found.

## 3. Not applied

- **Replacing the two per-row loops with whole-slice comparisons.** The length assertions close the
  same hole and keep the failure message that names the genotype and its allele counts, which is
  what makes a parity failure diagnosable.
- **A `proptest` over the whole domain, and a 64-thread agreement test.** Same reasons as at C1:
  the laws are asserted over explicit grids that name the failing shape, and the cache is
  per-thread over a pure function of its key.
- **Amending the plan's step C2 text.** See §4.

## 4. Raised, not fixed: the plan asks for something the compiler refuses

Step C2 says the comparison is made "against `GenotypeShape` … called directly". It cannot be:
`posterior_engine.rs` declares `mod shape;` privately and imports `GenotypeShape`/`shape_for`
privately, so there is no route from `src/ng/`, and the only way to satisfy the step as written is
to widen a declaration in the frozen tree. Confirmed three times independently on this branch,
each getting ``error[E0603]: module `shape` is private``.

**This loop does not edit the plan's text**, only its checkboxes. The substitution and its
justification are in the module's own doc comment, in the review, and in the implementation report
— but the plan is what the next reader opens first, and as written it sends them either to
re-derive the refusal or to widen `mod shape;`. Amending that sentence is the owner's call.

## 5. Validation

Run in the dev container after every fix:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --lib --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling` | 0 | `test result: ok. 53 passed; 0 failed; 0 ignored; 0 measured; 3962 filtered out` |
| `cargo test --lib` | 0 | `test result: ok. 4004 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out` |

Four tests in `ng::calling::genotype_table_parity`, one more than the three submitted; 53 in
`ng::calling` against 49 after C1.
