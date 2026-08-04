# Fixes applied — D1 review (`the_first_filter_needs_no_reference` → `reference_free_first_filter`)

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `f5630f8`
**Review:** [`ng_read_filtering_stages_d1_2026-08-04.md`](ng_read_filtering_stages_d1_2026-08-04.md)
— three agents in isolated worktrees (reliability; naming/smells/idiomatic; refactor-safety +
diff-matches-intent). Raw findings under
`tmp/review_2026-08-04_d1-first-filter-no-reference/`.

**Verdict: the step was rebuilt, not patched.** Four Majors across three reviewers, and they agree
on one thing: **the test as submitted had zero unique detection power.** Every reviewer ran the
checkpoint mutation and found it already failing at `cursor.rs:706` — production code, so
`cargo build` breaks — and at `mod tests`'s `pre` helper. One reviewer went further and *broke spec
§5's property outright* with the submitted test green.

---

## 1. What the reviewers established, and what I re-measured

| claim | who | re-measured here |
|---|---|---|
| The checkpoint mutation fails at four sites, two of them pre-existing and one not a test | all three | ✅ confirmed |
| The submitted assertions are a strict subset of `filtering.rs`'s own, so no behavioural mutation reaches them first | reliability | ✅ confirmed by inspection |
| `use super::*` in the module, or merging it into `mod tests`, leaves everything green — "the import list is the assertion" is unenforced | refactor-safety | ✅ accepted, not re-run |
| **`ReadFilterConfig::default()` launders a reference added to the config**: reviewer added a working `&dyn RawRefSeq` field behind `Default`, had the first filter read it, and got **2,867 passed, 0 failed** with the test green | reliability | ✅ confirmed and sharpened — §2 |
| Spec §5's capability is *not* reachable through `AlignmentCursor`, whose `R: RawRefSeq` bound is unconditional | reliability | ✅ confirmed; raised at Checkpoint D, not fixed here |
| The header-wide contig fetch went at **B1**, not C2 | naming, refactor-safety | ✅ spec §5's own amendment says so; prose corrected |
| "this module is the only thing that says so" is false | naming, refactor-safety | ✅ `verdict_on_raw_read`'s doc and arch §3.2 say it too; deleted |
| "one of the two things (spec §8)" — §8 says three | naming, refactor-safety | ✅ the count is the plan's; citation corrected |

## 2. The one mechanism with unique detection power, measured

The reviewers' fix was a function-pointer coercion. It is in, and it is the right statement of the
property — but **it is not an alarm nobody else raises**: the same mutation fails `cursor.rs` and
`filtering.rs`'s `pre` helper with `E0061`. Neither is the sibling-module placement unique: the
visibility mutation breaks `cursor.rs` too.

So I looked for something that *is* unique, and there is exactly one. Reproducing the reliability
reviewer's escape — a reference field on `ReadFilterConfig` — and then repairing it the way its
author would, i.e. supplying it in `Default` **and** in `filtering.rs`'s `post_config`, whose
module already has an `InMemoryRefSeq` in scope:

```text
error[E0063]: missing field `reference` in initializer of `filtering::ReadFilterConfig`
   --> src/ng/read/reference_free_first_filter.rs:100:5
error: could not compile `pop_var_caller` (lib test) due to 1 previous error
```

**One error, in the new file, and nowhere else.** That is why the config is built field by field,
and why the module has to live somewhere a reference cannot be produced to repair it.

## 3. What was applied

**Rebuilt as `src/ng/read/reference_free_first_filter.rs`** — a `#[cfg(test)]` sibling file under
`read/`, declared from `read/mod.rs` exactly as `left_align_parity` is.

| change | finding it answers |
|---|---|
| A `const _FIRST_FILTER_TAKES_NO_REFERENCE: fn(u16, MapQual, &ReadFilterConfig) -> FilterVerdict` coercion | naming Major, refactor-safety Major 1, reliability Minor — verified failing with `E0308` |
| The config built by **exhaustive literal**, no `default()`, no `..` | reliability Major 2 — the only unique detection, §2 |
| A whole reference-free **pass**: reader → narrowing → first filter → tally | reliability Major 3, refactor-safety Major 2 — the shape spec §5's three callers would write |
| Moved out of `filtering.rs` into a sibling module | refactor-safety Major 2(3) — a child of `filtering` sees `pub(in crate::ng::read)` items however private they become |
| Module renamed to a noun phrase, `reference_free_first_filter` | naming Minor — every other named test module in `src/ng` is a noun phrase |
| Narrative moved from `///` to `//!` | both nits — `cargo doc` never renders a `cfg(test)` item's `///` |
| "until C2 … a fetch per contig" → the B1/C2 split as spec §5 states it | naming Minor, refactor-safety Minor |
| "the only thing that says so" deleted; "no test can" deleted | naming Minor, refactor-safety Minor |
| "(spec §8)" → the plan, with a note that §8 lists three | both nits |
| All imports via one route, no `super::`/`crate::` mixing | naming nit |
| The "whole of what a caller must produce" gloss corrected — the config still names filter #8's threshold | naming nit |

## 4. Two things I removed from my own draft, having measured them

Both were written into the rebuilt file and then deleted, because the measurement said they were
padding — the discipline the reviewers applied to the first draft, applied to the second.

- **`the_flag_bits_the_pass_sets_are_the_bits_the_filter_reads`** — asserted noodles' `Flags::*`
  against production's `FLAG_*`. Its stated justification was that a drift between the two crates
  would misroute every drop while the pass still passed. **False:** setting the secondary record's
  flag to `SUPPLEMENTARY` fails the pass test and leaves the flag-bits test green. The pass's
  exhaustive `ReadFilterCounts` equality already catches any misrouting.
- **`a_reference_free_pass_still_knows_which_read_group_a_drop_came_from`** — its justification
  ("if resolving a read group came to need a reference, the pass would still pass and this would
  not") is wrong: the pass drives the same narrowing, so it would not compile either. The one
  assertion worth keeping folded into the pass's loop.

## 5. Not applied, and why

- **Move the `const` beside `verdict_on_raw_read` in non-test code** (reliability Minor). It would
  make `cargo build` carry the check, which is a real gain — but it puts a line in a shipping file,
  and Milestone D is scoped as tests only. Recorded on the const itself as an option for whoever
  wants it.
- **Amend spec §5** to say the capability stops at `AlignedReadsReader` and is not reachable
  through `AlignmentCursor<R: RawRefSeq>` (reliability Major 3, second half). This is a design-doc
  change and this skill does not edit design docs. Written into the new module's doc as an open
  point and raised at **Checkpoint D**.
- **Reword plan D1's "construct it"** (naming cross-category). Predates C2 making the first filter a
  free function; a plan-text matter for the checkpoint.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,867 passed**, 0 failed, 5 ignored |
| `cargo test --examples` | 52 passed, 0 failed |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |
