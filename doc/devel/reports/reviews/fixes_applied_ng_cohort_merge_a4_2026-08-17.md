# Fixes applied — ng cohort merge, step A4

*2026-08-17, branch `ng-cohort-merge`. Input:
[the A4 review](ng_cohort_merge_a4_2026-08-17.md) — 1 Blocker, 2 Major, 11 Minor, 7 Nits
over three category checklists. Every finding is accounted for below.*

## Findings table

| ID | Title | Severity | Decision | Status |
|---|---|---|---|---|
| B1 | the chained-width test did not pin what it claimed | Blocker | Apply | **Applied** |
| M1 | `min_alt_obs` never exercised at any value but 2 | Major | Apply | **Applied** |
| M2 | the width test applies to STR loci, against spec §3.1 | Major | Ask | **Raised at Checkpoint A** |
| Mi1 | `Failed`'s "counted" had no owner | Minor | Apply | **Applied** |
| Mi2 | "the runs where it matters most" named the wrong runs | Minor | Apply | **Applied** |
| Mi3 | `span_of`'s "the one spelling" untrue; the file used the other | Minor | Apply | **Applied** |
| Mi4 | the "not counted" rule cited to a section without it | Minor | Apply | **Applied** |
| Mi5 | `LocusCloser::over` documented neither new parameter | Minor | Apply | **Applied** |
| Mi6 | `bound` where the spec says `max_cohort_locus_span` | Minor | Apply | **Applied** |
| Mi7 | the bound never read at two different values | Minor | Apply | **Applied** |
| Mi8 | the ceiling invariant unchecked on the verdict wiring | Minor | Apply | **Applied** |
| Mi9 | spec §15's bystander fixture absent | Minor | Apply | **Applied** |
| Mi10 | wrap the width in a `LocusSpan` newtype | Minor | Defer | **Deferred** |
| Mi11 | allocate `members` only for a `Build` locus | Minor | Defer | **Deferred** |
| Nits | seven | Nit | Apply / Defer | **2 Applied, 5 Deferred** |

## The Blocker

**B1 — the test that named spec §3.2's central claim did not pin it.** Its doc said "an
implementation checking members rather than the closed locus would build it", and its
widest member was 11 bases against a bound of 10 — so a member-wise rule failed the locus
too, and mutating the walk to judge the widest member **passed all 31 tests**.

Fixed by repairing the fixture to spans 10, 10 and 7, closing to the same 21 bases, and by
**asserting** the member widths at `<= 10` rather than describing them. The test's doc now
records the earlier fixture and what it let through, so the reason the numbers are what
they are does not have to be rediscovered.

**Verified after the fix:** the member-wise mutant fails exactly that test — 31 passed, 1
failed, against 34 passed, 0 failed on the real code.

## The Major that was not applied

**M2 — `max_cohort_locus_span` governs generic loci only (spec §3.1), and `judge` applies
it to every locus.** Gating on kind needs a rule the design does not have: §4.1's grouping
chains observations of any kind into one cohort locus, and nothing says how a locus whose
members are an STR tract in one sample and a SNP in another takes a kind. That is design,
not implementation, so it is **raised at Checkpoint A** rather than guessed. It is inert
today — nothing feeds the walk until milestone C — and it is recorded at `judge` with the
question stated and the milestone that must settle it.

## The rest

`min_alt_obs` now has a test at 1, where the threshold changes character rather than
degree — spec §15's `keep_threshold_one_is_variant_filter` row. The bound is asserted at
two different values, so a rule comparing against a number of its own fails. The verdict
wiring is pinned at the coordinate ceiling, where `region.len()` panics in debug and
answers 0 in release — the answer that would let the widest locus expressible pass the
width verdict. Spec §15's bystander fixture is in: one sample's over-wide deletion taking
another sample's ordinary SNP down with it, which is what a failed locus *costs*.

Four doc claims were corrected, three of them because they were wrong rather than vague:
the "counted" obligation now names who counts; the ordering argument no longer claims to
protect the long-read case, which it does not, and names the low-coverage case, which it
does; and `span_of` no longer claims to be a uniqueness it does not have. The policy
parameter is `max_cohort_locus_span` everywhere, matching the spec and its own sibling.

## Deferred, with reasons

- **A `LocusSpan(u64)` newtype.** The comparison now lives in one function whose ordering
  is directly tested; B3 is where the module's type vocabulary next gets a considered pass.
- **Allocating `members` only for a `Build` locus.** The fix makes the field's meaning
  depend on the verdict, which is worse than the allocation until a consumer exists to
  justify it. A4's own doc says where that is decided.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` **34 passed, 0 failed** (31 before). Full-suite
figures are in the commit message.

**Three mutants re-run after the fixes, each now killed:** the member-wise width rule
(B1), `region.len()` in the verdict wiring (Mi8), and the two branches of `judge` swapped
(the step's own ordering guard, which fails exactly the two tests that exist for it).
