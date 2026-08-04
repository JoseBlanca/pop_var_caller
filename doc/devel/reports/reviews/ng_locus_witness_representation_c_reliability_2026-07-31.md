# Milestone C — reliability review (ng locus-witness representation)

*Reviewed at `82b13a0` (C1 `ebe3685`, C0 `6805e42`, C2 `761d53e`, C3 `82b13a0`), branch
`ng-pileup-generator`. Method: read the four commits against the spec/arch, then **mutate**.
Every mutation below was compiled and run; every quoted output is real.*

**Baselines, measured here.**

```
ng::locus_generation      test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.55s
--lib --bins --tests --examples --all-features
                          test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.04s
                          (+ 12, 6, 5, 11, 11, 4, 1, 11, 5, 10, 3 across the example/bin targets, all ok)
```

## Mutation ledger

| # | mutation | result |
|---|---|---|
| A | `read_agreed_with_reference`: revert C3's `runs.len() != 1` guard to the pre-C3 enclosing-extent read | **green, whole suite** — finding 1 |
| B | `refold_live_reads`: drop `witnessed.refill_from(...)`, keep the stale witness | 4 failed — sound |
| C | `witness_of`: clamp the enclosing extent instead of each run | 7 failed — sound |
| P | `witness_of`: drop the empty-after-clamp filter (`first < past_last` → `<=`) | 1 failed, at the `expect` — sound but thinly pinned |
| E | `num_obs_along_locus`: delete both `.min(len)` clamps | **green, whole suite** — finding 4 |
| Q | `WitnessedLocusPositions::is_flush_right`: `>=` → `==` | **green, whole suite** — finding 6 |
| F | C0 guard in `ssr::classify_read`: left border only | **green, whole suite** — finding 2 |
| G | C0 guard in `ssr::classify_read`: deleted outright | **green, whole suite** — finding 2 |
| R | probe: assert `refill_from`'s *documented* failure contract | fails (1 test) — finding 5 |

Mutations A, E, Q, F and G were each re-run against the **whole** suite
(`--lib --bins --tests --examples --all-features`) because the generic-dump tests C3 flipped
live in `examples/`, outside the narrow target. All five stayed green there too
(`2752 passed; 0 failed` on the lib target and `0 failed` on every example/bin target). The tree
was restored between every mutation; `git diff --stat HEAD` is empty at the end.

---

## 1. Major — C3's one design decision, `read_agreed_with_reference`, is pinned by nothing

`src/ng/locus_generation/pileup/open_record.rs:520-534`

C3's commit message and the implementation report both single this out as "one decision C1
deferred here": a multi-run witness names no single reference slice, so it answers `false` and
the read keeps its chain id. The code is a three-line guard:

```rust
let mut runs = state.witnessed.runs();
if runs.len() != 1 {
    return false;
}
let (start, end) = runs.next().expect("exactly one run, checked a line above");
```

**Mutation A** put the pre-C3 body back verbatim — the enclosing extent, which is the exact
behaviour C3's 15-line doc block argues against ("the fabrication this milestone removes,
wearing a different hat"):

```rust
let (start, mut end) = runs.next().expect("a witnessed set is never empty");
for (_, run_end) in runs {
    end = run_end;
}
```

Output — `ng::locus_generation`:

```
############ MUTATION A (src/ng/locus_generation/pileup/open_record.rs) ############
test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.60s
```

and the whole suite, examples included:

```
############ FULL-SUITE MUTATION A (src/ng/locus_generation/pileup/open_record.rs) ############
test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.47s
(every other target: ok, 0 failed)
```

The two existing multi-run walk fixtures (`holey` = `ACNTA`, `holed` with an `N` at 11) cannot
discriminate: their concatenated bases are shorter than the enclosing reference slice, so both
readings answer `false` for different reasons. The chain-id output is identical, which is why
nothing moves.

**Why it matters.** `chain_ids` is emitted per observation and the invariant
`chain_ids.len() <= num_obs` rests on this answer. A future edit that "simplifies" the guard
away restores the pre-C3 semantics silently, and on RNA-seq — the only input where holes exist —
it would start claiming a spliced read agreed with reference bases inside its intron.

**Fix.** One unit test in `open_record.rs`'s test module, on a record whose *enclosing* slice
the read's bases do equal. E.g. reference `b"ACG"` at `pos = 100`, a folded read in the `ACG`
bucket with `witnessed` runs `[(100,101), (102,103)]`: pre-C3 compares `reference[0..3] ==
b"ACG"` and answers `true`; C3 answers `false`. That single fixture kills mutation A.

---

## 2. Major — C0's guard catches both borders, but the test that says so cannot fail (and the guard is now dead code)

`src/ng/locus_generation/ssr.rs:757-759` (the guard), `ssr.rs:858-863` (the fallback),
test `classify::tests::a_read_covering_only_a_flank_is_outside_the_tract`.

The C0 commit message states the test's purpose explicitly — "asserts both borders, because they
arrive by mirror-image routes and a fix catching one would halve the population and look like it
worked" — and quotes the output of deleting the guard. **Both claims stopped holding at C2.**

**Mutation F**, guard the left border only:

```rust
// MUTANT F: only the left border is guarded.
RepeatSpan::FromLeft(tract) if tract.is_empty() => {
```

```
############ MUTATION F (src/ng/locus_generation/ssr.rs) ############
test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.76s

############ FULL-SUITE MUTATION F (src/ng/locus_generation/ssr.rs) ############
test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.30s
(every other target: ok, 0 failed)
```

**Mutation G**, delete the guard outright (both arms fall through to `partial`):

```
############ MUTATION G (src/ng/locus_generation/ssr.rs) ############
test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.72s

############ FULL-SUITE MUTATION G (src/ng/locus_generation/ssr.rs) ############
test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.50s
(every other target: ok, 0 failed)
```

Both also green on the whole suite (see the ledger). The reason is C2's own change: `partial()`
now ends with

```rust
let Some(read_witness) = (match border { … }) else {
    return Classified::NoObservation(NoObservationReason::OutsideTract);
};
```

and for an empty tract `reach == 0`, so `from_left(0, len)` and `from_right(0, len)` both answer
`None` and both produce `NoObservation(OutsideTract)` — byte for byte what the guard produces.
The guard is therefore **behaviourally unreachable as a decision**: no input can distinguish the
two code paths, so no test can pin it, and the commit message's quoted "deleting the guard"
output no longer reproduces.

This is not a wrong answer today. It is a guard that survives its own deletion, sitting behind a
test whose doc comment claims it discriminates.

**Fix.** Pick one home for the decision and say so:
- either delete the `classify_read` arm and let `partial()`'s `else` be the single point (the
  test then genuinely pins the surviving code), keeping a comment that records why the earlier
  guard is redundant; or
- keep the guard and give it a test that reaches *it* rather than the classification — e.g. a
  `#[cfg(test)]`-visible assertion that `partial()` is never called with an empty tract.

Either way, correct the test's doc comment and the C0 commit's claim in the implementation
report, because a reader currently believes a mutation is covered that is not.

---

## 3. Major — `ng_ssr_aligner_bakeoff` sums three of C0's four no-observation reasons

`examples/ng_ssr_aligner_bakeoff.rs:372-373`

```rust
counts.reads_without_observation =
    gen_counts.no_border_anchored + gen_counts.low_quality + gen_counts.window_truncated;
```

C0 added `SsrGeneratorCounts::outside_tract` and updated `examples/ng_ssr_loci_dump.rs`, which
even carries a warning against this exact slip:

```rust
// Every reason, named — `outside_tract` is the largest of the four on real data, so
// summing three of them would report a fraction of the reads that yielded nothing.
report.reads_without_observation = counts.no_border_anchored
    + counts.low_quality
    + counts.window_truncated
    + counts.outside_tract;
```

The bake-off example was not updated. The comment immediately above the broken line says "the
header numbers are the authoritative accounting identity, not a re-tally of the emitted rows",
and the struct's doc says "the accounting header, mirroring `ng_ssr_loci_dump`" — so this header
now silently under-reports by the *largest* term. On the numbers C0 measured (tomato chr01,
`SRR7279503`) that is 6,704 reads missing from a total of ~9,265, and `reads_fetched` no longer
balances against `obs_complete + obs_partial + reads_without_observation + reads_capped`.

Not found by mutation — there is no test on this example to mutate, which is the second half of
the finding.

**Fix.** Add `+ gen_counts.outside_tract`, and consider giving `SsrGeneratorCounts` a
`reads_without_observation()` method so the sum has one home and a fifth reason cannot be
forgotten twice.

---

## 4. Minor — `num_obs_along_locus`'s clamp is called "the guard" and nothing pins it

`src/ng/locus_generation/mod.rs:104-110`

The comment above it is emphatic: "**This clamp is the guard, not a second one** … `ReadWitness`
cannot know its own locus, so the invariant is not expressible on the type — it can only be
checked here, against the region actually in hand. Clamping rather than `debug_assert`: … a
debug-only guard compiles out of the release build this repo actually runs (a trap it has
recorded hitting twice)."

**Mutation E** deleted both clamps:

```rust
// MUTANT E: the consumer-side clamp deleted.
let from = start as usize;
let to = (end as usize).max(from);
```

```
############ MUTATION E (src/ng/locus_generation/mod.rs) ############
test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.79s

############ FULL-SUITE MUTATION E (src/ng/locus_generation/mod.rs) ############
test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.68s
(every other target: ok, 0 failed)
```

The nearest test, `a_run_reaching_beyond_the_locus_is_clamped`, does not exercise it: it builds
its witness with `ReadWitness::from_left(9, LocusLen::from_positions(3))`, which the *constructor*
already clamps to `(0,3)`, so the run handed to `num_obs_along_locus` is in range. Without the
clamp the failure is a slice-index panic, i.e. the exact release-build crash the comment says it
exists to prevent — reachable because `ReadWitness::Partial`'s field is public and
`WitnessedLocusPositions::from_half_open_runs` will build any canonical set.

**Fix.** A test that goes round the constructors:
`obs(b"AAA", ReadWitness::Partial { positions: WitnessedLocusPositions::from_half_open_runs([(0, 40)]).unwrap() }, 4)`
on `region(1, 3)`, asserting `vec![4, 4, 4]`. Add a second run wholly past the end
(`[(0,2),(30,40)]`) so both the `from` and the `to` clamp are covered.

*Related nit:* with canonical runs (`start < end`) `.max(from)` in the same line is unreachable —
`min(len)` is monotone, so `to >= from` always. Harmless, but it reads as a live guard.

---

## 5. Minor — `refill_from`'s documented failure contract is not what the code does

`src/ng/locus_generation/pileup/witnessed_ref.rs:87-106`

```rust
/// … `false` — with both this set and the buffer unchanged
/// — when the runs describe nothing.
pub(super) fn refill_from(&mut self, buf: &mut WitnessedRefRuns) -> bool {
    let Some(canonical) = canonicalise_runs(std::mem::take(buf)) else {
        // `take` emptied the buffer; give it back, so "nothing witnessed" leaves both
        // sides exactly as they were.
        return false;
    };
```

**Nothing gives it back.** `std::mem::take` has already emptied `buf` and the early return drops
the taken runs. The comment describes an action that is not in the code, and the doc's "the
buffer unchanged" is false. `take_from`'s doc has the same over-claim ("`None` — and the buffer
untouched").

**Probe R** added the documented assertion to the existing test
`taking_and_refilling_leave_the_callers_buffer_empty_and_reusable`, which today stops one line
short of checking the buffer:

```rust
buf.push((7, 7));
assert!(!set.refill_from(&mut buf), "an empty run is not a set");
assert_eq!(buf.as_slice(), &[(7, 7)][..], "PROBE R: …the buffer is unchanged on failure");
```

Result:

```
############ PROBE R (refill_from buffer contract, lib only) ############
test ng::locus_generation::pileup::witnessed_ref::tests::taking_and_refilling_leave_the_callers_buffer_empty_and_reusable ... FAILED
thread '…::taking_and_refilling_leave_the_callers_buffer_empty_and_reusable' (33647) panicked at src/ng/locus_generation/pileup/witnessed_ref.rs:224:9:
assertion `left == right` failed: PROBE R: refill_from's documented contract says the buffer is unchanged on failure
  left: []
 right: [(7, 7)]
test result: FAILED. 299 passed; 1 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.96s
```

**Severity is Minor because the *behaviour* is the one the fold wants** — an empty buffer is
exactly what the next `apply_events_into` needs, and the only caller (`refold_live_reads`) treats
`false` as unreachable. The defect is that a reader is told a rollback happens; someone relying
on it (a future second caller that wants to retry with the same runs) gets silence.

**Fix.** One-line doc change on both constructors: "the buffer is left **empty** either way". If
the rollback is actually wanted, it is `*buf = taken;` before the `return false`.

*Also worth tightening while there:* in `refold_live_reads` (`open_record.rs:1733-1742`)
`*allele_index = new_index;` runs **before** the refill, so a `false` refill would leave the read
pointing at its new bucket with its old witness — a mismatch caught only by a `debug_assert`.
Assigning the index after a successful refill costs nothing and removes the ordering hazard.

---

## 6. Minor — `is_flush_right`'s `>=` is justified in prose and pinned by nothing

`src/ng/locus_generation/witness.rs:165-169`

```rust
/// `>=` rather than `==` for the same reason the one-run predicate uses it: a producer's
/// reach is measured in read bases, which diverge from locus positions under stutter, so
/// a run may be handed a length the locus does not have.
```

**Mutation Q** replaced it with `==`:

```
############ MUTATION Q (src/ng/locus_generation/witness.rs) ############
test result: ok. 300 passed; 0 failed; 1 ignored; 0 measured; 2456 filtered out; finished in 9.80s

############ FULL-SUITE MUTATION Q (src/ng/locus_generation/witness.rs) ############
test result: ok. 2752 passed; 0 failed; 5 ignored; 0 measured; 0 filtered out; finished in 41.56s
(every other target: ok, 0 failed)
```
 Every fixture in the tree ends its last run exactly at
`locus_len`, because `from_left`/`from_right` clamp — so the over-long case the doc names is
never constructed. It is constructible through the public `from_half_open_runs`.

**Fix.** Extend `flushness_is_read_off_the_first_and_last_run` with a run past the locus end:
`from_half_open_runs([(0, 12)])` against `LocusLen::from_positions(9)` must still be flush right.

---

## 7. Minor — the STR reason-counter test covers three of the four reasons

`src/ng/locus_generation/ssr.rs:1485-1515`. The test is named
`no_observation_reasons_are_counted_by_reason_and_in_total` and its doc says "**Each**
no-observation reason lands in its own run-level counter", but C0's `OutsideTract` was not added
to its fixture, so the tally arm at `ssr.rs:1164`
(`NoObservationReason::OutsideTract => counts.outside_tract += 1`) is never executed by a test
that checks where it lands. Combined with finding 3, C0's counter is the least-covered thing the
milestone added and the one with the largest population on real data.

**Fix.** Add the fourth outcome to the fixture and assert `counts.outside_tract == 1`.

---

## 8. Minor — after C3, `reads_without_observation` is structurally always zero on the generic path, and the field still promises otherwise

`src/ng/locus_generation/mod.rs:50-53`, `pileup/open_record.rs:1988` and `:1491-1517`.

C3's own doc says the no-observation path is "now hard to reach at all". It is in fact
unreachable: `process_position` skips a contributor whose window is empty (`open_record.rs:1988`),
and every event that survives `events_overlapping`'s bounds test contributes a non-empty clipped
run (Match clips to `[lo,hi)`; Insertion requires `lo <= anchor < hi`; Deletion requires
`anchor < hi && anchor+len+1 > lo`). So `apply_events_into` cannot return `false` from either
call site, and `note_no_observation` never fires. Two walk tests now assert exactly that
(`record.reads_without_observation.is_empty()`).

Meanwhile the emitted field still says:

```rust
/// Reads that covered this locus but produced no observation at all. A scalar with
/// no positions: *no coverage* and *coverage that said nothing* are different
/// states, and only one means "look at the mapping" (spec §3).
```

A read whose every base is `N` or adaptor-masked — the class that doc describes — is filtered
before the fold and lands only in the run-level `reads_silent_over_footprint`, never in the
locus. So on generic loci the field is a constant 0 and spec §1 goal 2 is served only by the STR
generator.

**Fix (documentation, at least).** State on the field that the generic path always reports 0
since C3 and name `PileupGeneratorCounts::reads_silent_over_footprint` as where that population
went; or, if the per-locus distinction is wanted, record the silent read against the records its
window covers instead of `continue`ing past it.

---

## 9. Nit — stale doc after C0

`src/ng/locus_generation/ssr.rs:1122-1123`: "`counts` carries the **run-level** totals
(complete/partial observations and the three no-observation reasons)". There are four.

---

## What I checked and found sound

- **`witness_of`'s per-run clamp (C2's core claim).** Mutation C replaced it with a clamp of the
  enclosing extent and **7 tests failed**, including
  `witness_of_clamps_each_run_rather_than_the_extent_enclosing_them`,
  `witness_of_a_witness_with_a_hole_is_not_complete_however_far_its_ends_reach`,
  three walk fixtures and `every_divergence_from_production_is_one_of_the_six_named_classes`.
  Genuinely pinned, and by more than one test.
- **`Complete` decided on total coverage, not on the outermost edges.** Same mutation; the
  hole-at-both-borders fixture fails. Correct as written: the runs are disjoint and inside the
  locus, so `positions_covered() == end - pos` is exactly "covered every position".
- **The `expect` at `open_record.rs:226-227` is unreachable from the fold, and I tried.**
  Every run stored on a read is clipped into `[record_pos, record_end)` by `apply_events_into`
  at fold time (both call sites pass the record's own `pos` and `alleles[0].seq`), a record's
  anchor never moves, `alleles[0]` only grows, and `evict_unsupported_alleles` never touches
  index 0 — so no stored run can lie outside the *final* footprint. The `u16` narrowing cannot
  saturate a run to empty either: `PileupGeneratorConfig::check` rejects
  `max_record_span > MAX_RECORD_SPAN_CEILING (= u16::MAX)`, `to_walker_config` is `pub(super)`
  so the ceiling cannot be bypassed, and `open_new`/`widen` both refuse a wider span — hence
  `past_last - record_pos <= 65_535` is exactly representable. Mutation P (removing the
  empty-after-clamp filter) reaches the `expect` from exactly one test and nothing else:
  ```
  ---- …::witness_of_drops_a_run_that_falls_outside_the_footprint_entirely stdout ----
  panicked at src/ng/locus_generation/pileup/open_record.rs:228:10:
  test result: FAILED. 299 passed; 1 failed; …
  ```
  which is the right shape — the panic exists, and no walk input reaches it.
- **`apply_events_into`'s buffer contract, both directions.**
  `apply_events_into_clears_the_witness_buffer_it_was_handed` asserts the leak case and the
  empty case; the C1 report's own mutation output confirms the leak is caught. The `false`
  branch cannot leave runs behind (the only pushes are inside the loop after `clear()`, and the
  loop only reaches them when `witnessed` is `Some`, which is also what makes the final `match`
  return `true`). No stale run can leak between reads or between records.
- **The hole branch itself.** `event_start > run_end` on half-open runs is the right test;
  `apply_events_adjacent_events_stay_one_run` pins the off-by-one downward and
  `apply_events_a_hole_in_the_middle_is_recorded_as_two_runs` plus two walk fixtures asserting
  `runs().len() == 2` on a **one-position** hole pin it upward (a `> run_end + 1` mutation would
  fail them). The runs pushed are ascending and disjoint under precondition 4, and
  `take_from`/`refill_from` canonicalise anyway.
- **`refold_live_reads`' swap.** Mutation B (keep the old witness, clear the buffer, report
  success) fails four tests including two named for the property and the six-class parity
  anchor:
  ```
  test …::finalise_counts_the_witnesses_against_the_footprint_the_record_ended_with ... FAILED
  test …::widening_updates_the_witnessed_extent_of_a_read_that_is_not_a_contributor ... FAILED
  test …::a_read_whose_witness_splits_at_a_widen_keeps_its_observation ... FAILED
  test …::parity::every_divergence_from_production_is_one_of_the_six_named_classes ... FAILED
  test result: FAILED. 296 passed; 4 failed; …
  ```
  I also traced the skip cases: a contributor skipped by `refold_live_reads` is re-folded by the
  fold loop into the same record (its window over the record is non-empty because it once folded
  there and the window only grows), and an expired read keeps a witness that is still correct
  because it is held in absolute coordinates and re-clamped at `finalise`.
- **`canonicalise_runs` under the two witnessed-set types.** One implementation, generic over
  `(T, T)`, exercised by a proptest that checks representation *and* the position set; the
  reference-axis type has its own fixtures for merge, containment and holes. `start >= end`
  rejects both the empty and the reversed run. Nothing here needed a second look.
- **`sort_key` over a set.** `(u8, &[(u16,u16)])` is a total order, `Complete` sorts first, and
  `witness_order_is_total_over_witnesses_of_several_runs` covers the prefix case that a
  fixed-width key could not separate.
- **`WitnessedLocusPositions`'s `Hash`/`Eq` as an observation-identity component.** `SmallVec`
  hashes as a slice, so an inline set and the same set after a spill compare equal — pinned by
  `a_set_that_merged_down_to_two_runs_equals_the_same_two_built_directly`. The STR tally keys a
  `HashMap` on `ReadWitness` and is safe for the same reason.
- **C0's classification is what the spec asks for.** `locus_generation_ssr.md` §3 says a partial
  must have witnessed at least one tract position and names `NoObservationReason::OutsideTract`;
  `Between` with an empty tract is deliberately left alone. The code matches. (Only the *test*
  for it is the problem — finding 2.)

---

## Suggested order of work

1. Finding 3 (one line, wrong number in a tool that is being used for decisions).
2. Finding 1 (one test; it is the only thing standing between C3's decision and a silent revert).
3. Finding 2 (decide where the decision lives, then fix the test and the report's claim).
4. Findings 4, 6, 7 (three small tests).
5. Findings 5, 8, 9 (doc corrections), plus the `allele_index`-before-refill ordering nit.
