# ng generic locus generator — prerequisites, Milestone B: the shared locus type

**Date:** 2026-07-28 · **Plan:**
[locus_generation_pileup_prerequisites.md](../../ng/impl_plan/locus_generation_pileup_prerequisites.md)
steps B1–B3 · **Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md) §6, §10
· **Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md) §1.2

Implementation report for Milestone B, the second and last milestone of plan 1 of 3. One report for
the milestone; the three steps are separate commits because the plan requires it, and B1's separation
is load-bearing rather than tidy (below).

## 1. Plan

Four changes to types **both generators fill**, landed together because each done separately would
rebaseline the STR fixtures again. `ReadCoverage` reshaped; `ObservedSequence` given the read group
and `placed_left`; the docs folded in. This is what settles the type plan 2 will fill, which is why
Checkpoint B gates it.

## 2. Assumptions

None that changed direction. Three judgement calls are recorded in §6.

## 3. Changes made

### B1 — `ReadCoverage` → `Complete` + one `Observed` run (`6750d73`)

`PartialLeft(u16)` / `PartialRight(u16)` become
`Observed { offset_in_locus: u16, positions_covered: u16 }`. Two side-tagged variants cannot describe
what a read witnesses once the **events**, not the alignment span, define it — a read can be blind in
the middle of a footprint (an interior `N`, a ref-skip) or at either end, and a widened record can be
wider than the read on both sides.

Prefix-versus-suffix survives as a **derivation**, not a variant:

```rust
pub fn from_left(positions_covered: u16, locus_len: u16) -> Self;
pub fn from_right(positions_covered: u16, locus_len: u16) -> Self;
pub fn is_flush_left(&self) -> bool;
pub fn is_flush_right(&self, locus_len: u16) -> bool;
```

**The constructors clamp the reach before deriving the offset, and the order is the whole point.**
The STR path's reach is in *read* bases, which stutter makes diverge from locus positions — an
expanded allele reaches further than the reference tract has positions. Deriving
`locus_len - positions_covered` from an over-long reach would wrap; with a saturating subtraction it
would silently relabel a right-anchored read as left-anchored, which is the failure mode this step
was warned about.

Blast radius, exactly as the spec predicted: the enum, `num_obs_along_locus`, the STR generator's four
minting sites, its complete/partial tally and its sort key, and seven dump tools. The two mint sites
that passed the variant **as a function value** could not be retyped — a struct variant is not a
function — so `partial` now takes an `AnchoredBorder` marker plus the locus, which it needs anyway to
place a right-flush run.

### B2 — the read group and `placed_left` (`5054515`)

`read_group: ReadGroupId` joins the observation's **identity**: the STR tally's bucket key becomes
`(bases, read_coverage, read_group)`, so an allele seen from two read groups is two rows. It also
joins the **sort** key, and that is load-bearing rather than tidy — without it two cells differing
only by group would tie, and `HashMap` iteration order is seeded per process, so the output would be
non-deterministic run to run on any multi-group sample.

`placed_left: u32` is folded per read with production's own rule (`alignment_start < anchor`, strictly
left). `placed_start` is deliberately not added.

### B3 — the doc fold-ins (`e76fe2a`)

Docs only. Four fold-ins, written as **dated notes that preserve the original arguments** rather than
overwrite them, so the displaced reasoning stays readable: `locus_generation.md` §3 (the type block,
the `Observed` run, the cell-table consequence, the `reads_without_observation` lower-bound caveat),
`arch/locus_generation.md` §1 and its reconciliation row, and `read_preparation.md` §3 — whose
"reuse production's `PreparedRead` as-is" decision is **reversed**, marked in all three places a
reader could land.

### Extra — the dump fixture made discriminating (`45dab6b`)

Not in the plan; see §6.3.

## 4. Tests added, and what they pin

| test | what it proves |
|---|---|
| `one_allele_from_two_read_groups_is_two_rows_that_sum_back` (`ssr.rs`) | the group split is **computed**, not defaulted — with `read_group` constant the two reads would merge into one row of two — and collapsing the axis recovers the single-group total |
| `placed_left_counts_only_reads_starting_strictly_left_of_the_anchor` (`ssr.rs`) | the read starting *exactly* on the anchor is the discriminating case: it separates `placed_left` from the `placed_start` ng does not carry, and an implementation using `<=` passes everything else |
| `a_read_running_off_the_left_is_a_right_partial` (`ssr.rs`, strengthened) | asserts `offset_in_locus == 6 - 4` on a 6-base tract — the one number a mint that forgot the locus length could not produce |
| `two_read_groups_split_the_rows_and_the_counts_sum_back` (dump) | end-to-end: 4 rows where the single-group fixture gives 3, run-level per-read totals unmoved, and collapsing the group axis reproduces the single-group dump exactly |
| `a_single_read_group_fixture_is_unchanged_by_the_group_axis` (dump) | the other direction — no `(bases, coverage)` cell appears twice when there is nothing to split |
| `the_fixtures_partials_are_asymmetric_and_so_can_catch_a_side_swap` (dump) | the fixture keeps the property that makes it able to fail (§6.3) |

## 5. Validation

Run in the container (`./scripts/dev.sh`), per commit:

- `cargo fmt --all --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — no diagnostics.
- `cargo test --all-features` — **2495 passed / 0 failed / 4 ignored** at B3 (2493 at B1).
- `cargo test --examples` — green (10 for `ng_ssr_loci_dump`).
- Host-native: 2495.

Two standard commands excepted by hand, red independently of this work — PROJECT_STATUS *Standing
project-wide items*.

**The oracles, both halves.**

*B1.* The STR dump over its committed fixture is **byte-identical** before and after, on all three
delimiters (unit-robust, unit-slip, flat-gap) — stronger than the "equivalent modulo the encoding"
the plan asked for, because the labels are now derived from flushness and come out the same. Run
through the real binary over a materialised copy of the fixture, not through the test harness.

*B2.* Single read group: byte-identical to B1's output, again on all three delimiters — nothing
split, because there was nothing to split, which is what "free at one read group" has to mean. Two
read groups: the same six reads dealt across two groups give the row split and the sum-back above.

## 6. Deviations, and one defect found in the oracle itself

1. **`AnchoredBorder`, a new marker enum in `ssr.rs`.** The plan says the two function-value mint
   sites must be "restructured, not retyped" but not into what. A marker plus the locus is the
   smallest shape that works, and it makes the locus-length dependency explicit at the call site
   instead of hiding it in a closure.
2. **`from_left` / `from_right` / `is_flush_left` / `is_flush_right` are new API** the plan does not
   mention. Without them every mint and every label site would repeat the clamp-then-derive
   arithmetic, which is precisely where a wrong depth would come from.
3. **B1's oracle could not fail, and fixing it is its own commit (`45dab6b`).** The dump's fixture had
   two *symmetric* partials — `pl` and `pr` witnessed the same 20 tract positions, so their rows
   carried identical bases and identical counts, differing only in the label; sorting put the
   left-flush run first either way. **Swapping left for right at the mint site left the dump
   byte-identical**, checked by applying the swap. The acceptance anchor was blind to the one
   property the reshape most affects.

   The property was not uncovered — both `classify` side tests fail under that swap — but the
   *anchor* was, so a third partial (`pl2`, reaching 12 positions instead of 20) breaks the symmetry
   and the swap now moves the output. Verified by re-applying it.

   It is a separate commit rather than folded into B2 because **B2 turned out to need no rebaseline
   at all** — the single-group dump came through byte-identical — so there was no shared rebaseline to
   ride on, and the plan's "one rebaseline, not four" principle did not apply. The fixture change is
   better isolated than buried.

## 7. Open

- The reshape's one recorded consequence: **a run covering the whole locus is flush with both
  borders**, so a right-anchored read whose reach reaches the locus length is not distinguishable
  from a left-anchored one. Inherent to one-run-plus-offset, not to this implementation, and it reads
  correctly — a read that witnessed every position is unconstrained from either side. Documented on
  `from_right`.
- **The STR path's rows now split by read group**, which is a change to that path's output on any
  multi-group sample. Whether the STR cohort work then *replaces* its inferred sample groups with
  declared ones is that work's call; this milestone only makes it expressible.
- Milestone B ends the prerequisites plan. **Plan 2** (copy the walker, prove it identical) is next,
  and Checkpoint B is the gate: it cannot start until the type it fills is settled, which it now is.
