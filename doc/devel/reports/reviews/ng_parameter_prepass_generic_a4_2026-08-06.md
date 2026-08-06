# Code Review: ng_parameter_prepass_generic_a4
**Date:** 2026-08-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** commit `54378b9` — step 4's parameter pre-pass, plan step A4 (`DepthBin`, `DepthBinEdges`)
**Status:** Request-changes → resolved

---

### 1. Scope

- **What was reviewed:** one commit's diff — one new file, `src/ng/parameter_estimation/generic/depth_bins.rs`, plus a one-line `pub mod` in its parent.
- **Reviewed against:** `54378b9ead8771355451c874c889640da316f9bd`, checked out detached in three isolated worktrees.
- **In-scope files:** [depth_bins.rs](../../../../src/ng/parameter_estimation/generic/depth_bins.rs), [generic/mod.rs](../../../../src/ng/parameter_estimation/generic/mod.rs), and the step's implementation report.
- **Out of scope:** earlier commits on the branch; Milestones B–G, unwritten by design.
- **Categories dispatched:** `reliability` + `errors` (the step's failure mode is silent, so "would a fault be caught, and by what?" is the whole question); `naming` + `defaults` (three behaviourally significant constants, and a docs-heavy file); `idiomatic` + `smells` + `refactor_safety` + `module_structure` (a generator with clamps, and a placement that deviates from the plan's file list). `unsafe_concurrency` skipped — no `unsafe`, `Arc`, lock, atomic or thread. `tooling` skipped — `Cargo.toml` untouched.

### 2. Verdict

**Request-changes**, and the reason is that this is the step whose failure is silent. Two real defects, both found by driving the generator with shapes the adopted constants never produce — not by reading. Both are now fixed.

### 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::parameter_estimation` | 0 | 13 passed |
| `cargo doc --no-deps --lib` | 101 | 12 unresolved links, all pre-existing, none in this file |
| `cargo test --all-targets --all-features` | 101 | one pre-existing failure in `ng::locus_generation::pileup::parity`, out of scope |

Findings labelled "Needs verification": **zero**. Every behavioural claim was produced by mutating the file in an isolated worktree and re-running, with output quoted in the per-category file.

### 4. Open questions and assumptions

1. **Should the architecture's module table gain `depth_bins.rs`?** (affects Mi7) `arch/parameter_prepass_generic.md`'s "Module home" tree names four files under `generic/` and maps §2.2 to `histogram.rs`; the plan's A1 names the same four. All three reviewers agree the fifth file is right, so the documents are now behind the code. Editing a design document is outside this skill's remit — recorded as an owner item.

### 5. Top 3 priorities

1. **M1** — the guard against a degenerate ladder is dead under the adopted constants, and its two clamps are in the wrong order, so when it *does* fire it creates the empty bins it claims to prevent.
2. **M2** — `row_start` answers a plausible wrong number for a bin index one past the end, where its sibling `depth_range` panics on the same input.
3. **M3** — the implementation report overstates the oracle: 2/3/3 failing tests where the truth is 2/2/2, and the five other tests are structurally incapable of catching a constant edit.

### 6. Findings

#### Major

**M1: depth_bins.rs:112-120 — the ladder's degeneracy guard is dead, and its clamps are ordered backwards**
**Categories:** reliability, smells, errors — **three of three, convergent**
**Confidence:** High.
The line `bin_tops.push(top.max(previous + 1).min(MAX_BINNED_DEPTH))` is commented as preventing an empty bin and an overshoot of the cap. Instrumented at the adopted constants, **both clamps fire zero times** — the raw geometric tops already land strictly increasing and hit 124.000000 exactly. Replacing the whole line with a bare `push(top)` leaves all seven tests green, so it was the file's only untested branch, guarding the parameter this step was isolated to protect.

And the order is wrong for the stated purpose: `.min(cap)` runs last and wins, so the guard **creates** the condition it claims to prevent. Driven with other shapes: `exact=8, cap=12, bins=20` gives seven bins the range `13..=12`; `exact=8, cap=10` gives tops `[9, 10, 10, 10, …]`. An inverted `RangeInclusive` returns `false` from `contains` for every depth, so those bins exist in `bin_count()` and own cells in the flat vector while no site can ever land in them — a permanently empty histogram row, which the plan names as the failure "no fit would report".

**M2: depth_bins.rs:151-154 — `row_start` answers where `depth_range` refuses**
**Categories:** errors, idiomatic
**Confidence:** High.
`DepthBin`'s field is public, so an out-of-range bin is constructible. `row_starts` is one longer than `bin_tops`, so `row_start(DepthBin(20))` reads the trailing total and returns **583** — which is `cell_count()`, a number shaped exactly like a row offset, one past the end of the vector the caller is about to index. `depth_range(DepthBin(20))` panics on the identical input. Milestone B indexes a flat cell vector at `row_start(bin) + alt_reads`, so an off-by-one bin index would address another bin's cells rather than failing.

**M3: the implementation report overstates the oracle's breadth**
**Categories:** reliability, naming — convergent
**Confidence:** High.
The report claims 2 / 3 / 3 failing tests for the three one-constant mutations. Measured by two reviewers independently, each fails exactly **2**, and the same two every time — the eleven literal bin tops and the 583-cell count. The reason is the real finding: the other five tests are written *in terms of* the shape constants, so their expectations move with the mutation and they are structurally incapable of detecting a constant edit. The oracle protecting this correctness parameter is **two** tests, not seven. The 2/3/3 figures are also in the commit message, which is now immutable history.

#### Minor

**Mi1: the module doc quotes the cap's cost at the wrong bin count.** **Category:** naming. The sentence says raising the cap 124→300 "at a fixed bin count doubles the error-rate bias and quadruples" the other, with "the same twenty bins" as its subject. Those multipliers are the research note's **sixteen**-bin row (0.545→1.038, 1.83%→7.96%). At twenty bins the note gives 0.054→0.190 rungs and 0.30%→0.88%, i.e. 3.5-fold and 2.9-fold. This was the only in-code statement of the cap's cost at the adopted ladder, and it was wrong in both directions.

**Mi2: `EXACT_DEPTH_LIMIT`'s doc leads with the wrong quantity.** **Category:** naming. "How many depths get a bin to themselves … depths `0..=8`" — that count is nine and the value is 8. It is the deepest depth kept exact, which is what `DEPTH_BIN_COUNT - EXACT_DEPTH_LIMIT - 1` relies on; a reader following the doc computes twelve widening bins instead of eleven.

**Mi3: `max_depth()`'s uniqueness claim is contradicted 130 lines above.** **Categories:** naming, refactor_safety. Its doc says "nothing else in step 4 states this number" while `pub const MAX_BINNED_DEPTH` states it publicly. The architecture flags precisely this as a failure that already happened once in the design — an earlier draft set the subsampling cap to 300 while the ladder ended at 124.

**Mi4: no build-time check on the three shape constants**, where the sibling `generic/mod.rs` guards its three ladder constants with a `const _: () = assert!` and states why. `DEPTH_BIN_COUNT = 8` is caught only by rustc's `arithmetic_overflow` lint const-folding a subtraction — luck, not design — and nothing covers `MAX_BINNED_DEPTH < EXACT_DEPTH_LIMIT` or `EXACT_DEPTH_LIMIT = 0`, the latter making the ratio infinite and every top come out of the `+1` clamp. **Categories:** refactor_safety, errors.

**Mi5: the derived `PartialEq`/`Eq` on `DepthBinEdges` is vacuously true**, and it sits exactly where the design mandates `Arc::ptr_eq`. Measured: two independently built edges compare equal. The type has one constructor taking no arguments, so `==` cannot distinguish any two values — and a future `merge` author writing `assert_eq!(self.edges, other.edges)` gets a check that compiles, reads correctly and cannot fail. **Category:** refactor_safety.

**Mi6: the "330 MB" arithmetic does not work as written.** **Category:** naming. 5,151 cells × 8 bytes is 41 kB; 330 MB needs the ×8,000-windows step the spec supplies. And `MAX_BINNED_DEPTH`'s "losing the extra reads costs nothing" is stated as settled where research §4.6 lists sites deeper than the cap under *what is not measured*.

**Mi7: the design documents were not updated for the fifth file.** **Category:** module_structure. All three reviewers agree the placement is right — the binning rule fits neither side of `generic/`'s data-shaping-vs-mathematics split, the architecture names three consumers for it inside `generic/`, and there are no consumers outside step 4 so lifting it further would be wrong. But `arch`'s module table and the plan's A1 still name four files and map §2.2 to `histogram.rs`, so a reader following them looks in the wrong place.

**Mi8: three `pub const`s with no caller outside the module**, and no `#[must_use]` on any of the seven accessors where the sibling `error_rate_ladder()` has one. **Categories:** refactor_safety, idiomatic.

**Mi9: two `mut` accumulator loops where `scan` states the same rule**, which also removes the `.last().expect(...)` that exists only because `previous` is recovered from the vector being built rather than carried. Verified byte-identical. **Category:** smells.

**Mi10: `Default` and `bin_for`'s documented totality are undertested** — `Default` is never called by any test, and totality over `u32` is checked at three sampled points where a property test is available. **Category:** reliability.

#### Nits

`bin_for`'s exact-region branch is not strictly necessary (a single `partition_point` over the whole ladder gives the same answer, verified over `0..=1000` and at `u32::MAX`), but nothing said why it is there. `EXACT_DEPTH_LIMIT as usize + 1` appears three times and is written a fourth way elsewhere. `(...) as u16` truncates rather than saturates above 65,535 bins. Three `.expect()` sites lack the `// PANIC-FREE:` comment the sibling module uses.

### 7. Out of scope observations

- The pre-existing `ng::locus_generation::pileup::parity` failure, unchanged.
- `pub mod` throughout `ng/` puts step-4-local types on the crate's external surface, where `arch` records them as having no consumer outside this step. A scaffold-wide convention question, not this commit's.

### 8. Missing tests to add now

1. `a_cap_too_low_for_the_bin_count_is_rejected` and `more_bins_than_the_cap_has_room_for_is_rejected` — degenerate shapes through the parameterised generator. These are the tests that would have caught M1.
2. `row_start_rejects_a_bin_past_the_end_of_the_ladder` and its `depth_range` twin — M2.
3. `default_edges_are_the_adopted_ladder` — Mi10.
4. `bin_for_answers_an_existing_bin_for_every_depth` — a proptest over arbitrary `u32`, since the returned index addresses a table.
5. `every_ladder_shape_is_strictly_increasing` — over several shapes, so the clamp that buys it is exercised rather than assumed.

### 9. What's good

- **The generator-checked-against-a-literal-oracle direction is right**, and one reviewer checked this explicitly rather than assuming: `new()` generates the tops from three constants and the test pins the eleven values the research note measured. That is a generator checked against an independent oracle, not two copies of one rule.
- **The fifth-file placement survived scrutiny from the category most likely to object**, on three independent grounds.
- **Field rename, reorder and addition are all compile errors** at the single struct literal; no `..Default::default()`, no non-exhaustive destructure, no manual field-set-dependent impls.
- **`depth_range` returning `RangeInclusive`** rather than a pair was cited as putting inclusivity in the type where its only consumer is a row width.
- **The cast audit came back clean** at the shipped constants — every one either widening, guarded by the branch above it, or saturating by Rust's float-to-int rules.

### 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --all-targets --all-features`
- `./scripts/dev.sh cargo doc --no-deps --lib`

Per-category files kept as an audit trail in `tmp/review_2026-08-06_ng-param-A4/`.
