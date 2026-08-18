# ng cohort merge — E2: overlap resolution, and why it is a safety net

*Implementation report, 2026-08-18. Step E2 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §6.1, with §3.2 and §3.3 for what a failed locus is.*

> **This is the first draft, and [the review](../reviews/ng_cohort_merge_e2_2026-08-18.md)
> changed the evidence in it.** The argument below survived 120,000 randomised layouts through
> the real merge. What did not survive is what §3 offered as proof of it: the disjointness guard
> sat *after* the byte-identity comparison, where the cached output is already the oracle's and
> nothing can overlap — so it could never fire, and eight tests failing on a deliberately broken
> driver included none of it. Nothing fed a real builder's output to a real organiser either.
> Both are fixed, and four of this report's own claims were wrong; the corrected account is in
> [the fix report](ng_cohort_merge_e2_fixes_2026-08-18.md).

## 1. Plan

Implement spec §6.1's rule at the organiser's release point: **of two loci that overlap, the
one whose first position is earlier stands, whether it was emitted or failed — one rule, no
special case.** Then settle, before writing a single test, whether the rule is a live one or a
safety net, because the plan prescribes a fixture that turns out to demonstrate the opposite of
what it was written for.

## 2. The rule, in code

Two fields on `Organiser` and two methods, in
[`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs):

- `owned_through: Option<GenomePosition>` — the last base a resolved locus owns;
- `displaced_locus_count: u64` — how often the rule fired;
- `resolve_and_release` — one region's loci and failed spans, taken as **one** genome-ordered
  sequence;
- `claim` — the rule itself: a locus whose first base is at or before `owned_through` loses its
  ground; anything else takes it and moves the frontier to its last base.

**The two vectors are merged rather than handled apart** because the rule is about ground, and
a failed locus owns ground exactly as an emitted one does. Two separate passes would let a
failed span be judged against a frontier that had already run past it — which is
`loci_and_failed_spans_are_resolved_in_one_genome_order`, and the mutation that pins it.

**The test is the first base, not the intersection.** A locus belongs to the builder whose
region holds its first position, so the earlier *start* is the earlier owner. On any pair the
walk can produce the two tests agree; they differ on pairs it cannot, and the first-base form
is the one spec §6.1 states.

## 3. The rule is a safety net, and this is the argument

**Under `build_region`'s input contract, two loci owned by different building regions cannot
overlap.** C1's review established it for the whole-stretch driver and left open whether D2's
division changed it; it does not. The argument, which is now in the `Organiser` doc comment
because a branch nobody can explain is worse than no branch:

1. **A builder is handed every observation overlapping its own ground, including ones that
   opened before it.** The cache evicts only what *ends before* the region's first base
   (`ObservationCache::evict_before`), so everything reaching into the region survives.
2. **So a locus L owned by an earlier region, reaching into this one, arrives connected.** Take
   any member of L that starts inside this region: it reaches at least its own start, so it is
   in the window; so is the member it overlaps, and so on backwards. The sub-chain is unbroken
   and it begins before this region's first base.
3. **So this builder closes that ground as one locus starting before its own first base**, and
   skips it as an earlier region's (`build_region`'s ownership rule). Every locus it *does* own
   begins after that chain ends — which is after L's own reach.

So the rule is kept for three reasons and not because anything reaches it: it is spec §6.1's;
the contract it rests on belongs to whoever feeds the builders rather than to `organise.rs`;
and it costs one comparison per locus.

**Tested as a safety net, which means from both ends.** The organiser's eleven new tests are
all **fabricated** — they say what the rule does if it is ever reached. The other end is
`serial.rs`'s `refuse_overlapping_ground`. *(The review corrected this paragraph twice over: the
guard was placed where it could not fire, and its reach was overstated. What it covers, and the
end-to-end check added beside it, are in the fix report.)*

**What the plan asked for, and what changed.** E2's step says the fixture "must contain" a wide
deletion beginning before a building region and reaching into it. That fixture exists —
`a_locus_reaching_across_a_building_region_boundary_is_built_once_and_whole` in `serial.rs` —
and it demonstrates the *exclusion*, not the rule: the locus is built once, whole, by the
region that owns its first base, and nothing overlaps it. Writing a test that claimed to
exercise displacement through it would have been a test whose fixture cannot reach the branch.
It is instead the fixture the end-to-end organiser check runs on — see the fix report.

**`displaced_locus_count` is how a run would say the argument had failed**, rather than leaving
it to be inferred from a locus quietly missing from the output. It is expected to be zero.

## 4. Deviations and things found

**E2 found five fixtures in E1's own tests that no builder could produce** — failed spans that
overlapped each other, and locus/failed pairs sharing a first base, inside a single region's
outcome. One walk over one region closes disjoint loci and judges each one, so a region's two
lists are ascending and mutually disjoint. The five are corrected, and `outcome_of`'s doc now
states the constraint. That the rule fired on them at all is the clearest evidence it works.

**One dead branch removed after mutation testing.** `claim` first wrote the frontier with a
`max`, "so it can only move forward". A mutation replacing that `max` with plain assignment
survived every test — because a locus reaching `claim` starts past the frontier and ends at or
after its own start, so its last base is past the frontier and the `max` can never bite. It is
now a plain assignment with the reason stated. A branch no input can take reads as a hazard
someone guarded against, and sends the next reader looking for the case.

**Nothing about the design changed**, and nothing outside this step was touched.

## 5. Tests

Eleven in `organise.rs`, two in `serial.rs`, one guard called from a shared helper. The module
went from 189 tests to **202** (`cargo test --lib ng::run::cohort_merge`).

| test | what it pins |
|---|---|
| `a_locus_starting_on_ground_an_earlier_locus_owns_is_dropped` | the rule |
| `the_boundary_is_the_earlier_locus_last_base` | opening *on* the last base loses; on the next base wins |
| `a_failed_locus_displaces_what_starts_inside_its_span` | spec §3.2 — the case the rule exists for |
| `a_failed_locus_displaced_by_an_earlier_one_is_not_counted` | a refusal the run never owned stays out of the total |
| `loci_and_failed_spans_are_resolved_in_one_genome_order` | one merged sequence, not two passes |
| `one_wide_locus_displaces_every_locus_that_opens_inside_it` | three dropped behind one 291-base locus |
| `a_locus_kept_after_another_owns_its_own_ground_too` | the frontier follows the latest owner, not the first |
| `a_locus_on_the_next_contig_is_never_displaced` | ownership does not cross a contig |
| `the_frontier_carries_across_an_empty_region` | the frontier is the organiser's, not one region's |
| `resolution_follows_region_order_rather_than_arrival_order` | why the release waits for the predecessor |
| `nothing_is_displaced_when_no_two_loci_overlap` | the ordinary run |
| `the_disjointness_guard_refuses_two_overlapping_loci` | the guard is not vacuous |
| `the_disjointness_guard_weighs_a_failed_span_as_ground` | it weighs failures as ground |

**Twelve mutations, twelve killed** (`tmp/mutate_e2.sh`). Nine in `organise.rs`: no
displacement at all; the boundary made strict; the frontier set to the region's first base;
failed spans not displacing; a displaced failure still counted; the two lists resolved as two
passes; the frontier freezing at the first locus kept; ownership crossing a contig; the counter
never rising. Three in `serial.rs`: the guard's condition made vacuous, the guard ignoring
failed spans, the guard ignoring the contig.

Two of those twelve survived a first round and are why two things in §4 changed — the frontier
`max`, and the guard having no test of its own.

## 6. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `202 passed; 0 failed`.
- `cargo test --lib` — `3825 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out;
  finished in 614.91s` (3,812 before this step; the thirteen new tests are the difference).
- `tmp/mutate_e2.sh` — 12 mutations, 12 killed.

## 7. Follow-ups

- **The counter has no consumer.** Where `displaced_locus_count` surfaces belongs with the
  failed-locus count, in the run summary the emission step owns (spec §13).
- **The exclusion argument rests on the cache's eviction rule**, which E3 may change when
  several builders read the cache at once. If E3 evicts less eagerly the argument only gets
  stronger; if it ever evicts *more*, the argument has to be re-made.
