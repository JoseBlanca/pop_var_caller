# ng candidate alleles — C1: admission and the remapping

*2026-08-24. Step C1 of [`../../ng/impl_plan/candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md),
branch `ng-candidate-alleles`, on top of `eff6cf16`. Design authority:
[`../../ng/spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §3, §6.1, §6.2, §8 and
[`../../ng/arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §2.3, §3.1, §3.2.*

---

## 1. Plan

The first step that produces the **whole** answer rather than a piece of it. `select_generic`
seeds `CandidateAlleles` with the merge's allele 0, admits every alternative that cleared the
admission rule **in the merge table's own order**, and fills `AlleleRemap` — the map from the
merge's allele indices onto the new dense candidate ids — as it goes.

The plan gives it its own commit for a stated reason: **an off-by-one in the remapping hands the
calling loop a real but wrong allele's evidence, which is a wrong genotype rather than a panic.**

The cap is C2 and the leftover is C3, so this step returns `SelectionVerdict::Selected` always and
one zeroed leftover entry per covering sample — the right length in the right order.

## 2. Assumptions and departures

**The test fixtures moved.** C1's tests need the same hand-built loci as B1's, so the seven fixture
helpers left `mod tests` for a `#[cfg(test)] pub(super) mod fixtures` — the shape
`src/ng/run/cohort_merge/mod.rs` already uses, and for the reason its own doc gives: a locus
fixture is fiddly enough that two copies drift. **The bodies are unchanged**; only their home and
their visibility. The 54 tests that were there before pass unaltered.

**Arch §3.1's sentence about the order cannot be implemented as written**, and writing C2 and C3
in a review worktree is what showed it. It says one pass "admits the survivors in table order,
applies the cap, and fills the leftover". The cap has to run *before* admission, because admission
needs to know what survived it; the leftover has to run *after*, because it is per sample where
admission is per allele and it reads the finished remapping. The real shape is four stages: fold,
cap, admit, leftover. Recorded in `select_generic`'s own doc comment; the arch sentence is raised
at Checkpoint C.

**The merge's distinctness invariant is now held in debug builds.** `CohortObservation::alleles`
documents that each sequence appears once and `AlleleTable` enforces it by interning on the bytes,
but a review probed a table holding one sequence twice and found both copies admitted as separate
candidates — one allele's evidence split across two, and a genotype that looks ordinary. A
`debug_assert!` rather than a release one, because the check is a scan of the table and the
invariant belongs to the producer.

## 3. Changes made

- `src/ng/calling/allele_candidates/generic.rs` — **new**: `select_generic`, and 16 tests.
- `src/ng/calling/allele_candidates/mod.rs` — `pub mod generic;`; the `#[allow(dead_code)]` on
  `summarise_alleles` and on `AlleleSummary::cleared_the_bar` removed, now that `select_generic`
  is a live root for both; the fixtures lifted into their own module; and two stale doc paragraphs
  corrected (see §5).

## 4. Tests added

Sixteen. The plan's oracle is the round trip **with a hole in the middle**: five sequences of which
the middle alternative is the only one no sample earned, so the merge's indices 1, 3 and 4 survive
as dense candidate ids 1, 2 and 3, and index 2 answers `None`. A reviewer confirmed the hole is
load-bearing: a remapping that *counts* rather than remaps satisfies all three of `AlleleRemap`'s
own assertions and is caught only by the fixtures that have a gap.

Also: the evidence hand-off of arch §3.2 reproduced row by row, with one allele shown from two read
groups that must both survive and stay apart; the reference-only outcome; admission in table order
against rank order; the two standing properties from the plan; the leftover's length and its
zeroing; a covering sample whose reads all stopped inside the locus; an empty allele table; and the
scratch reused across loci in both directions.

## 5. What the review changed

Three agents in three isolated worktrees; full account in
[`../reviews/ng_candidate_alleles_c1_2026-08-24.md`](../reviews/ng_candidate_alleles_c1_2026-08-24.md).

**Two Blockers, and for the third step running they were tests that could not fail** rather than
wrong code:

- **the whole admission rule could be replaced by a cohort read total** — `cohort_reads >= 2` in
  place of the per-sample question — and all 65 tests stayed green, because no fixture had two
  samples each lending an alternative *less* than the floor. A cohort term in that rule is the one
  thing this module exists to prevent, and it admits error alleles in proportion to cohort size;
- **the rule's share could be dropped** and the suite stayed green, because every fixture was
  shallow enough that the floor decided. The shipped share binds above 41 compared reads, which is
  exactly where the GIAB trio runs at 30× and 300×.

Three Majors: a covering sample with partials and no rows was in no fixture, so a leftover length
that skipped it survived — and that is a *shifted* leftover, not a short one; nothing asserted the
leftover entries were **zero**, so a non-zero `earned_reads_cut_by_the_cap` on the line C3 rewrites
would emit every covering sample at every locus as missing; and the duplicate-sequence case above.

**And one wrong mechanism in my own prose, proved by deleting the code.** The `# Panics` note said
an empty allele table would otherwise panic in `CandidateAlleles::new` about an empty reference
allele. A reviewer deleted the assertion and got `index out of bounds` two statements earlier —
`CandidateAlleles::new` is never reached. Corrected, with what actually happens.

## 6. Validation

All in the container, on the committed tree:

- `cargo fmt --check` clean; `cargo clippy --lib --tests --all-features -- -D warnings` clean;
  `cargo clippy --lib --all-features -- -D warnings` clean;
- `cargo test --lib` **4,265 passed, 0 failed, 14 ignored** in 45.6 s, against 4,249 at `eff6cf16`.

**Nine mutations, nine killed.** Five before the review — the rule inverted, the remapping off by
one, the reference offered to the rule, admission in reverse table order, and every survivor
admitted with the reference's bases. Four after, each killing a Blocker or Major the review found.

## 7. Tradeoffs and follow-ups

- **C2 must carry a test pinning admission order at a binding cap.** A reviewer deleted the
  sort-back that returns the kept prefix to merge-table order and **no C1 test noticed**, though
  the whole `ALT` column was permuted.
- **`CandidateAlleles::admit` accepts an empty allele where `new` refuses one**, with the same
  unparseable-record consequence one VCF column to the right. Out of this step's scope; the guard
  belongs in `calling/mod.rs`. Raised at Checkpoint C.
- **The `# Panics` note about 65,535 alternatives goes when C2 lands** and makes that width
  unreachable.
