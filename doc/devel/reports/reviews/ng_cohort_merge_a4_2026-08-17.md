# Review — ng cohort merge, step A4 (the two verdicts, width first)

*2026-08-17, branch `ng-cohort-merge`, working tree at stash commit `1b6150a8`. Three
category checklists, two sub-agents, one mutation-testing in an isolated worktree.
Per-category audit trail: `tmp/review_2026-08-17_ng-cohort-merge-a4/`.*

## 1. Scope

- **Reviewed:** `Verdict`, `judge`, `ClosedLocus::verdict`, the two new parameters on
  `LocusCloser::over`, `span_of`, and their tests.
- **Categories dispatched:** `reliability` (mutation, in a worktree), `naming`, `smells`
  — the last carrying a spec-fidelity pass against §3.1, §3.2, §3.3, §4.1, §4.3 and §15.
  **Not dispatched:** the rest, for the reasons A3's review gives; nothing new here has an
  error path, concurrency, or a module change.

## 2. Verdict

**Approve-with-changes.** One Blocker — in a *test's own claim about itself* — plus two
Major and eleven Minor. All applied but one, which is a design question for the owner.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | 31 passed, 0 failed (at review time) |

**Mutation testing: 8 run, 4 survived, 0 changed-no-behaviour.** Every survivor was shown
to answer differently from the real code on a constructed input, and each demonstration
was re-run against the unmutated tree before being recorded.

## 4. Findings

### Blocker

**B1 — `a_chain_of_narrow_observations_fails_on_its_closed_width` did not distinguish the
implementation its own doc comment said it distinguished.**
**Confidence:** High, measured. *(reliability)*
The test's comment claimed "an implementation checking members rather than the closed
locus would build it". Its widest member was `region(10, 20)` — **11 bases against a bound
of 10**, so a member-wise rule fails that locus too, for the wrong reason. Mutating the
walk to judge the widest member instead of the closure **passed all 31 tests**. Spec §3.2's
central claim — "what matters is how wide the locus ended up, not how it got there" — was
unpinned, in the test named for it.

This is the failure the plan's own skill warns about most loudly: **a wrong number in the
author's claim about the author's own fixture**, which survives because a reader re-reads
the sentence instead of re-running the test.

**Fixed:** the fixture is now spans 10, 10 and 7 — no member reaching the bound — closing
to the same 21 bases, and the member widths are **asserted** rather than described, at
`<= 10` rather than the old `<= 11`. **Verified after the fix:** the member-wise mutant now
fails exactly that test (31 passed, 1 failed).

### Major

**M1 — `min_alt_obs` was never exercised at any value but 2**, so substituting
`MinAltObs::DEFAULT` at the call site survived the whole suite, and spec §15's
`keep_threshold_one_is_variant_filter` row was unrun. *(reliability)*
**Fixed** by `at_a_threshold_of_one_any_non_reference_read_builds`, where the threshold
stops being a degree and becomes a plain "did anyone vary at all".

**M2 — the width test is applied to every locus, where spec §3.1 says
`max_cohort_locus_span` "governs **generic** loci only".** *(smells, fidelity)*
An STR locus's span is its reference tract — a fact about the reference rather than a claim
about the reads — and §3.1 expects the true ceiling on an emitted observation to be "the
larger of `max_cohort_locus_span` and the widest STR tract the segmentation admits", 100
rather than 50 at the two defaults. As written, a 60-base tract would be failed and
counted, inflating the one number §3.3 says an operator reads.
**Not fixed, and deliberately: this is a design question, raised at Checkpoint A.** Gating
on kind needs a rule the design does not have — §4.1's grouping chains observations of any
kind into one locus, and neither the spec nor the architecture says how a *cohort* locus
whose members are an STR tract in one sample and a SNP in another takes a kind. Inventing
that rule is not this step's to do. The divergence is inert until milestone C, since
nothing feeds the walk yet, and it is now recorded at `judge` with the question stated.

### Minor — applied

- **The `Failed` doc claimed "counted" with no owner in the code and no test.** Now names
  who counts: the builder that consumes this walk, which returns its region's failed-locus
  count (spec §3.3, §6.3), with the walk's part being to yield the locus rather than
  swallow it. *(smells)*
- **"in exactly the runs where it matters most" was wrong, and pointed at the wrong runs.**
  §3.3's named case is long-read data, where a genuine mid-size deletion carries far more
  than `min_alt_obs` non-reference reads and is `Failed` under either order. The loci the
  order actually moves are wide *and* below the threshold — the low-coverage corner. The
  doc now says which. *(smells)*
- **`span_of`'s "the one spelling" was not true**, and this file used the other three
  times, including inside the assertion of the width-rule test. Claim narrowed; the test
  now uses `span_of`. *(smells)*
- **`TooQuiet`'s "not counted" was cited to spec §4.3**, which never mentions counting.
  Recited to §4.3 for the threshold and §1.3/§3.3 for why only a failure is counted.
  *(smells)*
- **`LocusCloser::over` gained two parameters and documented neither.** Both now
  documented, with the property that matters at a call site: neither affects *closing*.
  *(smells)*
- **The parameter was named `bound` while the spec, and the file's own doc comments, call
  it `max_cohort_locus_span`** — where the sibling carries `min_alt_obs` verbatim. Renamed
  throughout; the fixture helper is `max_span` beside `keep_at`. *(naming)*
- **The bound was never read at two different values**, so a hardcoded `span > 10` passed
  all 31 tests. Now the same 11-base locus is asserted to build under a bound of 11.
  *(reliability)*
- **The coordinate-ceiling invariant was checked on `span()` but not on the verdict
  wiring**, so `region.len()` at the call site survived. Now pinned by
  `a_locus_at_the_coordinate_ceiling_is_judged_on_its_true_width` — built by hand, because
  the test helper sizes its reference bases with `region.len()` and panics there.
  *(reliability)*
- **Spec §15's bystander fixture was absent**: one sample's over-wide deletion suppressing
  another sample's ordinary SNP chained into it. Added — it is what a failed locus *costs*.
  *(reliability)*
- `Build`'s doc regained its cross-reference. *(naming)*

### Minor — not applied, with reasons

- **Wrap the locus width in a `LocusSpan(u64)` newtype.** Reasonable, and the threshold it
  is compared against is already a newtype. Deferred: the comparison now lives in exactly
  one function whose ordering is directly tested, so the newtype would buy less here than
  the churn costs, and B3's `SampleSupport` is where the module's type vocabulary next
  gets a considered pass. *(naming)*
- **Allocate `members` only for a `Build` locus.** Real — every dropped locus pays an
  allocation only a built locus's consumer reads. Deferred because the fix makes the
  field's meaning depend on the verdict, which is worse than the allocation until a
  consumer exists to justify it; A4 and B1 are where a real caller decides.
- Three nits on fixture-helper naming.

## 5. What's good

- The mutation agent proved all four survivors changed behaviour before recording them,
  and re-ran each demonstration against unmutated code.
- It separated the two ordering guards and showed they are not redundant: an
  implementation leaving `judge` correct but inlining the wrong order in `next()` is caught
  only by the walk-level test, so both are load-bearing.
- The fidelity pass enumerated seven distinct claims spec §3.2 makes about a failed locus
  and checked each against the code, rather than checking the summary.
