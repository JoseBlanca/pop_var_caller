# Code review — ng generic locus generator, Milestone D (D1–D3)

**Date:** 2026-07-30 · **Scope:** `e626850..31be8ee` on `ng-generic` (D1 `6993704`, D2 `7bfcd8a`,
D3 `31be8ee`) · **Fixes:** this commit ·
**Impl report:** [ng_locus_generation_pileup_generator_d_2026-07-29.md](../implementations/ng_locus_generation_pileup_generator_d_2026-07-29.md)

## The fan-out did not happen, and that is the first thing to say about this review

Five category agents were dispatched, each meant to run in its own git worktree, as the previous
three milestones did. **Seven of eight launches were killed by API 529 overload** — the five
originals and two retries, over about forty minutes. One agent completed: **reliability**, plus
the test-challenge pass, re-launched on a different model tier as a test of whether the overload
was tier-specific. It was.

So this is a **one-category review**. Four checklists did not run: `errors`/`defaults`,
`idiomatic`/`smells`/`unsafe_concurrency`, `module_structure`/`naming`/`refactor_safety`, and
`extras` (spec conformance and the measurements). Their prompts are ready to re-issue verbatim.
Milestone D should not be treated as fully reviewed on the strength of this.

**One process finding, from the agent that did run: the worktree isolation did not hold.** It
worked in `/Users/jose/devel/pop_var_caller-ng-generic` — the author's own worktree — applying and
reverting ten mutations there, and wrote its report into that tree's `tmp/`. No damage resulted
(the tree is clean at `31be8ee` and the suite re-verified at 2,724), but a mutation-heavy agent
sharing the author's tree is exactly the collision the per-agent worktree rule exists to prevent,
and with several such agents it would have been unrecoverable. Worth diagnosing before the next
fan-out.

## Verdict: 2 Majors, both "a mutation nothing catches", one of them Milestone D's own

The reliability agent re-ran all seven mutations the impl report claims are caught and confirmed
every one, then tried three of its own. Two of those three fail **no test** in either suite.

### Major 1 — `flush_all`'s `ever_contributed` guard was untested. **Fixed.**

`src/ng/locus_generation/pileup/active_read_set.rs`

`reads_silent_over_footprint` is fed by the active set's **two** exits and D2 pinned only one.
`expire_passed` — a read the walker has passed — is covered by
`a_read_silent_at_every_position_is_counted_rather_than_lost`. `flush_all` is the other, and on
the generic path it is **not an edge case**: a region walk stops at `region.end` while the reads
reaching into the halo are still active, so *every bounded walk ends by flushing reads that never
expired*. Deleting the guard there — counting every flushed read as silent — left the **whole
2,724-test suite green**, which the author reproduced before fixing.

The counter would then have over-reported on every region of every real run. That is the failure
mode this milestone exists to eliminate: a number nobody can see being wrong.

**Fixed** by `a_read_still_active_when_the_walk_stops_is_counted_by_what_it_contributed`
(`generator.rs`), and the fix needed two goes:

- **The first draft could not fail either.** With one silent read and one contributing read, the
  correct guard and a guard *inverted* to count the contributors both total 1. Two contributing
  reads make the two answers 1 and 2.
- Mutation-verified in **both** directions: guard deleted → 3 against 1; guard inverted → 2
  against 1; guard restored → green.

That is the thirteenth test on this branch found unable to fail, and the second in two milestones
to be a test written *for a review finding*.

### Major 2 — `refold_live_reads`' contributor-skip has no regression test. **Carried, with the reviewer's sketch.**

`src/ng/locus_generation/pileup/open_record.rs`

Deleting `if contributors.iter().any(|c| c.read_id == read_id) { continue; }` changes nothing
observable in either suite (202 lib + 10 dump), and makes the `contributors` parameter unused —
which confirms the skip is that parameter's only remaining use.

**Carried rather than fixed, for three reasons.** It is **A3's** code, not Milestone D's, and the
function's own doc comment already records the gap deliberately: *"Unpinned, and deliberately so…
Mutating the skip away leaves the whole suite green… Do not read the absence of a failing test
here as the absence of a reason."* The reviewer looked for a correctness counterexample and did
not find one — the carried `contribution` makes a double re-place idempotent — so the exposure is
a future edit breaking that invariant silently, plus a possible change in allele **creation
order**, which is observable in the output. And the test wants the record's internal allele table
at a chosen walker position, which is a fixture shape none of the existing ones have.

The reviewer's sketch, for whoever writes it:
`refold_live_reads_skips_a_read_that_is_also_the_widening_contributor` — build a record where read
`r` folds at one position, then at a later position both `r` has an event (so it is a contributor
there) and a second read anchors a deletion that widens the record; assert the allele order
matches what the fold loop alone would produce.

### Not filed: `apply_events_into`'s `run_end.max(event_end)`

The reviewer also mutated `witnessed = Some((run_start, run_end.max(event_end)))` to drop the
`.max()`; no test fails. Filed as a note rather than a finding, with the reasoning: within one
read's CIGAR-ordered event stream `event_end` is non-decreasing, so the two forms agree on every
input the debug-asserted precondition admits, and the `.max()` reads as defence for a case that
precondition already rules out. **Carried as a note to check against `decompose.rs`'s
events-overlapping construction** before anyone relies on it, since the reviewer could not build a
failing input either way.

## What the review confirmed

All seven prescribed mutations fail a test, and **one is caught more widely than the impl report
claimed**: `coverage_of` always returning `Complete` fails 9 lib tests plus **4** dump fixtures,
where the report said two. Corrected in the impl report.

The agent also audited the permanent anchor and its fixture, and could not break either:

- **The floor assertions bound something real.** `anchored > total * 5` and
  `anchored_multi_base > total` stop the fixture degenerating into single-base loci where the
  equality is vacuous; `widens > 0` stops it losing the one property that makes
  `generate_uniform_events` worth more than `generate_complete`. All three read accumulated
  per-run counts, not constants.
- **It could not construct an input where the fixture leaves a read stale**, and traced why by
  hand: every read on a contig is built from one shared template, so wherever one read's event
  triggers a widen every other live read has an event at that position and is a contributor there
  — not the "live but silent" read a stale fold needs. Both paths that can make a folded read stop
  being a contributor mid-record (the column cap, the mate-overlap collapse) operate on the same
  post-truncation list `refold_live_reads` is driven from, so the read is re-folded from its own
  cursor. The argument depends on no two reads on a contig having different starts or CIGARs,
  which is what the fixture-property assertion checks.

## One documentation contradiction. **Fixed.**

`parity.rs` — the anchor's doc said `generate_uniform_events` gives every read one event set "so
**no record widens at all**", and two paragraphs later said the opposite, correctly: "not by
stopping records widening, which it does not". The first was a leftover from the draft whose
`record_widen_events == 0` assertion had already failed at 7 widens; the test itself agrees with
the second (it asserts `widens > 0`). Corrected to say what is true — the shared event set removes
the *staleness*, not the widen — because the contradiction made the property impossible to verify
by reading.

## Validation after fixes

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` **2,725 passed** (2,724 + the new test); the example's 10 tests green;
`cargo doc --no-deps` still 12 pre-existing unresolved links.
