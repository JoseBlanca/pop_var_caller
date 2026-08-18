# ng cohort merge — E1: applying the review's fixes

*Fix-application report, 2026-08-18, against
[the E1 review](../reviews/ng_cohort_merge_e1_2026-08-18.md). Every finding was applied; two
items were deliberately not, and both are named in §3.*

## 1. What changed

**Four Majors and every Minor and Nit, in
[`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs).**

- **A gap at the tail of a run is now refused.** `is_finished` and `finish` both take
  `regions_handed_out` — how many building regions the run dealt out, their indexes being
  `0..regions_handed_out`. Without it a run that lost its last regions was indistinguishable
  from one that finished. `finish` also refuses an outcome submitted for an index the run says
  it never handed out, which is the same class of hand-out bug `submit` already refused.
- **The reorder map is written outside the assertion.** `submit` checks
  `!held_outcomes.contains_key(&index)` and then inserts. Writing the map inside the assertion's
  condition put the module's only insertion in an expression that a later edit to
  `debug_assert!` would stop evaluating — and this crate's release profile leaves debug
  assertions off, so every outcome would have been dropped in the shipped binary and in no test.
  Checking first also keeps the *first* outcome rather than replacing it on the way to the
  panic, which is now its own test.
- **`MissingRegionResults` became `RunEndedShort`, a three-variant enum.** The old struct's one
  message stated one cause for two faults — it told a run that had merely stopped draining that
  "a gap stalled the ordered drain" — and `{ 0, 0 }`, a run that ended short by nothing, was
  constructible. The variants are `RegionsNeverReleased`, `LociNeverDrained` and both together;
  the two with a stall carry `first_stalled`, the index whose absence held the drain, which is
  what an operator can map to a building region and to a builder. `RegionIndex` gained
  `Display` for those messages.
- **Both counts are pinned now.** `finish_names_both_counts_when_a_stall_and_undrained_loci_coincide`
  puts two regions never released against two loci never drained, and the two held regions carry
  three loci between them — so a count of *loci* reads 3 where the count of *regions* reads 2.
  `a_region_that_never_submits_holds_back_every_region_after_it` was restrengthened the same way
  (one, two and three loci over three regions). `the_refusal_names_each_count_against_its_own_noun`
  asserts the rendered message with unequal counts, so a swap inside the format string cannot
  pass.

**Everything else the review asked for:** `new()` written field by field and `finish`
destructuring `self`, so a field E2 or E3 adds must be answered for; `RegionOutcome`
destructured where it is consumed, as `build.rs` does; `held`/`released`/`next_expected` →
`held_outcomes`/`released_loci`/`next_expected_region`; `release_ready` →
`release_regions_in_turn`, leaving "ready" to the arch-mandated public name; `#[must_use]` on
the four methods that warrant it, `drain_ready` with a reason; the `PANIC-FREE` marker on the
cursor's `expect` and `saturating_add` on the counter beside it, with the difference explained;
`submit`'s doc reason rewritten around *when* a fault is caught rather than *what* it is about;
`drained` → `drained_regions`; the doc on `a_region_submitted_twice_is_refused` corrected — it
described the other test's fixture.

## 2. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `189 passed; 0 failed` (168 before E1, 184 after
  the first draft, 189 after these fixes: five more tests).
- `cargo test --lib` — `3812 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out;
  finished in 579.59s` (3,807 before the fixes; the five new tests are the difference).
- **`tmp/mutate_e1_v2.sh` — 15 mutations, 15 killed**, including the three the review found
  alive (the swapped message, `finish` reporting one count at a time, regions counted as their
  loci), the two forms of the trailing-gap defect, and the `debug_assert!` hazard reproduced as
  a mutation (insert-then-check, which lets the second outcome displace the first).

## 3. Not applied, and why

- **Splitting `organise.rs` into a cache file and an organiser file.** The review is right that
  the module doc's argument for keeping them together was false — `serial.rs` already calls
  `cover` and `evict_before` from a sibling module, so they can never become file-private, and
  `build.rs` cannot reach the cache at all. **The paragraph is corrected to say so.** The split
  itself moves D1 and D2's code and contradicts the architecture's file tree, so it is raised at
  Checkpoint E rather than taken inside E1.
- **A `RegionTicket` the organiser mints and `submit` consumes**, which would make both of
  `submit`'s panics unrepresentable and hand `finish` its own count. There is no hand-out loop
  to attach it to until E3; recorded there.

## 4. Follow-ups this run leaves

- **Arch §5 writes `RunError::MissingRegionResults { count: usize }`.** E1 ships a three-variant
  `RunEndedShort` instead, for the reasons above. Amending arch §5 is the owner's call.
- **The organiser still does not own the observation cache**, so D2's open item stands: a failed
  merge leaves the cache advanced, and making that unrepresentable waits on E3.
