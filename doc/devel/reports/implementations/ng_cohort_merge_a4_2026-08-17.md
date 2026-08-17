# ng cohort merge — A4: the two verdicts, width first

*Implementation report, 2026-08-17. Step A4 of
[the plan](../../ng/impl_plan/cohort_merge.md) — the step it marks **"own commit, do not
bundle"**; design authority [spec](../../ng/spec/cohort_merge.md) §3.2, §3.3, §4.1, §4.3
and [arch](../../ng/arch/cohort_merge.md) §3.*

## 1. Plan

Each closed locus gets one of three verdicts, decided **width first**: `Failed` (wider
than `max_cohort_locus_span`), `TooQuiet` (fewer non-reference reads than `min_alt_obs`),
`Build`.

## 2. Why the order is the step

**Both orders drop the locus, so nothing downstream can tell them apart.** A
reference-only chain wider than the bound qualifies for both verdicts: it is too wide to
build and too quiet to be worth building. What changes with the order is the **failed
count** — and spec §3.3 makes that count the only signal an operator has that the bound is
charging more than they expected. Judged width first, the locus is ground the caller
*refused*, and counted. Judged variability first, it is ground the caller *examined and
found empty*, and not counted; the count would then under-report, silently, in exactly the
runs where the bound is doing the most work.

That is why the plan gives this step its own commit: a defect here produces a plausible
wrong number rather than a crash, and `git bisect` needs it isolated.

**Measured, not asserted.** Swapping the two branches of `judge` fails exactly the two
tests that exist to pin the order — `the_verdict_is_decided_width_first` and
`a_reference_only_chain_wider_than_the_bound_is_failed_not_too_quiet` — and nothing else:
`test result: FAILED. 29 passed; 2 failed`.

## 3. Changes made

All in [`close.rs`](../../../../src/ng/run/cohort_merge/close.rs):

- **`Verdict`** — `Failed` / `TooQuiet` / `Build`, with the ordering argument on the type
  rather than buried at the branch.
- **`judge(span, non_reference_reads, bound, min_alt_obs) -> Verdict`** — a free function
  and a pure one, so the ordering can be checked directly rather than only through a
  fixture. It is the only place either comparison is written.
- **`ClosedLocus::verdict`**, filled by the walk. `LocusCloser::over` now takes the two
  parameters.
- **`span_of(region)`** — lifted out of `ClosedLocus::span` so the width a locus is
  *judged* on and the width it *reports* cannot become two different numbers.

**Both boundaries go the inclusive way, and both follow the spec's own words.** A locus
exactly `bound` bases wide is built — the bound is the widest the caller undertakes to
build, not the first width it refuses (§4.1, "a locus at most that wide goes on"). A total
exactly `min_alt_obs` is built — a locus is kept when its non-reference reads *reach* the
threshold (§4.3).

**Every closed locus still comes out, whatever its verdict.** A failed locus is not
silence: it owns its ground and displaces the loci that overlap it (§3.2, §6.1), which the
organiser cannot do for a locus it was never told about.

## 3a. What the review changed

[Review](../reviews/ng_cohort_merge_a4_2026-08-17.md),
[fixes](../reviews/fixes_applied_ng_cohort_merge_a4_2026-08-17.md). Three things are worth
naming here:

- **A Blocker in a test's claim about itself.** `a_chain_of_narrow_observations_fails_on_its_closed_width`
  said in its own comment that a member-wise rule would build the locus, and its widest
  member was 11 bases against a bound of 10 — so a member-wise rule failed it too.
  Mutating the walk to judge the widest member passed all 31 tests. The fixture is
  repaired and the member widths are now asserted, not described.
- **A Major that is a design question, raised at Checkpoint A**: spec §3.1 says the bound
  governs *generic* loci only, and `judge` applies it to every locus. Gating on kind needs
  a rule the design does not have.
- Four doc claims corrected, three of them wrong rather than vague — including one that
  claimed the verdict order protects the long-read case, which it does not.

## 4. Tests added

**Nine** — six written with the step and three added by the review, taking `close.rs` from
23 tests to 32. The review also strengthened three existing tests in place rather than
adding new ones (the repaired chained-width fixture, a second bound value, and the member
widths asserted rather than described).

The six written with the step:

- `a_reference_only_chain_wider_than_the_bound_is_failed_not_too_quiet` — **spec §15's
  fourth new test**: 21 bases against a bound of 10, and not one non-reference read, so it
  qualifies for both verdicts and the width one must stand.
- `the_verdict_is_decided_width_first` — the same rule on `judge` alone, all four cells of
  the two tests crossed, with the both-qualify cell asserted explicitly.
- `a_locus_exactly_at_the_bound_is_built` and `a_locus_exactly_at_the_keep_threshold_is_built`
  — the two inclusive boundaries, each with the neighbouring value that must go the other
  way.
- `a_failed_locus_is_emitted_with_its_ground_and_its_neighbours_are_untouched` — §15's
  first new test: the failed locus comes out with its span, and the loci either side are
  judged on their own terms.
- `a_chain_of_narrow_observations_fails_on_its_closed_width` — §15's second: three
  observations of spans 10, 10 and 7, none reaching the bound of 10, closing to 21. An
  implementation checking members rather than the closed locus builds it. **The member
  widths are asserted, and the review is why**: the first version of this fixture had a
  member of 11 — over the bound — so the member-wise rule it names failed the locus too,
  and the mutant passed the whole suite.

And the three the review added:

- `at_a_threshold_of_one_any_non_reference_read_builds` — spec §15's
  `keep_threshold_one_is_variant_filter` row. Every other fixture uses a threshold of 2, so
  substituting the default at the call site had survived the suite.
- `a_failed_locus_suppresses_the_bystander_variants_chained_into_it` — §15's first
  fixture, and what a failed locus costs: one sample's over-wide deletion takes another
  sample's ordinary SNP down with it.
- `a_locus_at_the_coordinate_ceiling_is_judged_on_its_true_width` — the verdict is wired to
  this module's span rather than `GenomeRegion::len()`, which answers 0 in release where
  the widest locus expressible would then pass the width verdict.

## 5. Validation

In the container: `cargo fmt --check` clean; `cargo clippy --lib --all-features --
-D warnings` clean; `cargo test --lib ng::run::cohort_merge` **34 passed, 0 failed** (25
before this step). Full-suite figures are in the commit message.

## 6. Tradeoffs and follow-ups

- **The failed count is not summed here.** Spec §3.3 requires it to reach the run summary;
  the organiser sums it (E1) and the emission step owns the surface (§13). A4 produces the
  verdict that makes the count possible.
- **Nothing consumes a verdict yet.** B1 assembles the `Build` loci; E2 uses the failed
  spans to displace overlapping loci.
- **The two scheduling-invariance tests spec §15 asks for** — the failed set identical at
  1, 2 and 8 builders and at two partitions — belong to **E4**, which is where builders
  exist to vary.
