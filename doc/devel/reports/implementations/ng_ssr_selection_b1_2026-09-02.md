# B1 — the repeat ladder: one tract's sequences grouped by how many repeats they carry

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone B step B1, executed on the STR loop plan's branch
([`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md) A1, §Branch).
**Design:** [`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §3;
[`arch/candidate_alleles_ssr.md`](../../ng/arch/candidate_alleles_ssr.md) §2.1.
**Module:** [`src/ng/calling/allele_candidates/ssr.rs`](../../../../src/ng/calling/allele_candidates/ssr.rs)
(new), with the scratch field in
[`mod.rs`](../../../../src/ng/calling/allele_candidates/mod.rs).

---

## What this step is, in one paragraph

At an ordinary locus the merge's alleles are an unordered set, and narrowing them is a support
bar and a cap. At a repeat tract they are ordered: a tract of 11 repeats is adjacent to one of 10
and far from one of 4, and the genotype prior and the stutter model are both written on that
ordering. **The ladder is that ordering, made explicit** — every sequence in the merge's table
placed on the rung named by `bases.len() / motif.period()`, floored, with the rungs ascending and
the cohort's most-supported rung recorded. It is what nomination, admission and the periodicity
test all read, and it computes the repeat count that the genotype prior must agree with to the
integer.

## The step's stated dependency did not bind

The plan makes B1 *depend on A1* — the merge carrying `LocusKind`, which the parallel
[observations plan](../../ng/impl_plan/run_ssr_observations.md) owns and has not merged to `main`.
**It did not block this step**, because arch §1 already settled the shape: the tract's detail is a
separate argument until the merge carries it, "the honest shape today and the shape that survives
when the merge fixes it". So `build_ladder` takes `motif: &Motif` and knows nothing about where it
came from. What still waits for A1 is the call site — `select_ssr` reading the motif off
`CohortObservation::kind` instead of being handed it — which is step D3's and the loop plan's.

## What landed

**`RepeatLadder`, filled by `build_ladder`, living in the shared `SelectionScratch`.** Its
accessors are the rung count, a rung's repeat count, a rung's merge-table indices, a rung's
cohort read total, the lookup from a repeat count to its rung, and the modal repeat count.

Four choices are worth naming, three of them departures the plan's latitude covers:

- **The rungs share one index buffer.** Arch §2.1 sketches `rungs: Vec<Rung>` with each rung
  owning a `Vec<u32>`, which is a heap allocation per rung per locus — the cost a scratch buffer
  exists to remove. The indices instead live in one buffer sorted by `(repeat count, table
  index)` and a rung names its slice of it. Same information, same order, nothing allocated per
  locus after the first.
- **The sort key is the pair, not the repeat count alone.** That makes the order total, so it
  does not depend on the sort's stability or on the order samples were walked in — which is what
  the byte-identical-at-any-worker-count requirement needs
  ([`spec/candidate_alleles.md`](../../ng/spec/candidate_alleles.md) §8). **This guarantee is
  structural and not test-covered**: at three alleles an unstable sort will not reorder anything,
  so no fixture this module can build could tell the two keys apart.
- **A rung's reads are summed from the fold, not from a second walk of the rows.**
  `summarise_alleles` already visits every sample's rows once and records each allele's cohort
  total; a rung's total is a sum over its sequences. So `build_ladder` must run after the fold on
  the same locus, which it asserts on the fold's width and states in its own documentation.
- **The build asserts the ladder is empty instead of emptying it.** See *What the mutations
  found*.

**`modal_repeat_count` is kept although the prior no longer takes it.** Arch §5 records that the
genotype prior's seed was re-indexed on 2026-08-27 by offset from the **reference** tract length,
so the cohort's commonest length is not an input to it. The mode is still needed here: the
periodicity test measures each read's offset from it (spec §7). Whether anything outside this
module should receive it is step D3's question, not this one's.

## Assumptions, where the design left a gap

**The reference tract's rung can exist carrying zero reads, and that is not production's
`occupied` test.** The merge interns the reference at index 0 whether or not a read landed on it,
so every rung here holds a sequence the table holds — but not necessarily one any read reached.
Production's `±1` rescue asks `cohort_support(length) > 0`
([`candidate_set.rs:221`](../../../../src/ssr/cohort/candidate_set.rs)), which is false for such a
rung. `rung_of_repeat_count` therefore answers *which rung*, not *is it occupied*, and a caller
wanting production's question asks it together with `cohort_reads_at`. Pinned by
`the_reference_rung_exists_with_no_reads_on_it` so that step C2, which builds the rescue, meets
the difference as a failing assertion rather than a doc comment.

## Tests — 14 new, and what each is for

The plan names three oracles for this step; they are the first three below.

| test | what it pins |
|---|---|
| `two_sequences_of_one_length_share_a_rung_and_stay_distinct` | two ten-base sequences differing by one interior base share rung 5 as two indices |
| `a_length_that_is_not_a_whole_number_of_units_lands_on_the_floored_rung` | seven bases at a dinucleotide join the six-base sequence on rung 3, counted once |
| `the_mode_is_the_rung_with_the_most_reads_and_not_the_reference_s` | 30 reads at four repeats beat the reference's 2 at five |
| `the_mode_breaks_a_tie_toward_the_shorter_rung` | twelve reads each, shorter rung wins |
| `a_rungs_reads_include_a_sample_that_cleared_no_bar` | the one read that failed the support rule still decides the mode, 7 against 6 |
| `a_repeat_count_the_table_does_not_hold_has_no_rung` | rung 4 absent between 3 and 5 answers `None`, not its neighbour |
| `the_reference_rung_exists_with_no_reads_on_it` | the difference from production's occupancy test, above |
| `a_homopolymer_keys_every_length_to_its_own_rung` | period 1, where floor division is the identity |
| `a_second_locus_leaves_no_rung_of_the_first` | a smaller second locus over the first's buffers |
| `building_a_ladder_twice_without_a_fold_between_is_refused` | the assertion that replaced the redundant clear |
| `an_unbuilt_ladder_refuses_its_mode` | 0 is a legal repeat count, so the answer must be a panic |
| `an_unbuilt_ladder_refuses_an_occupancy_question` | `None` everywhere reads as "the table holds none of them" |
| `a_ladder_refuses_a_rung_it_does_not_hold` | a reader walking another locus's rung count |
| `a_ladder_refuses_a_fold_from_another_locus` | the width check on the scratch's fold |

`resetting_the_scratch_leaves_no_value_from_an_earlier_locus` in `mod.rs` was extended to fill
and then check the ladder, because nothing else can fail if `reset_for` forgets it.

## What the mutations found

Five deliberate defects, each applied to the tree, run, and copied back:

| mutation | outcome |
|---|---|
| `reset_for` drops `ladder.clear()` | caught — `resetting_the_scratch_leaves_no_value_from_an_earlier_locus` |
| the mode's tie goes to the longer rung (`>` → `>=`) | caught — `the_mode_breaks_a_tie_toward_the_shorter_rung` |
| the key stops dividing by the period | caught — 8 tests |
| a merged rung stops adding its second sequence's reads | caught — `a_length_that_is_not_a_whole_number_of_units_lands_on_the_floored_rung` |
| `build_ladder` drops its own `ladder.clear()` | **survived** |

**The survivor was a real finding and the fix is not another test.** `reset_for` already empties
the ladder and the fold calls it, so the `clear()` inside `build_ladder` could not fail: it was a
second owner of one rule, and being the redundant one it could only ever hide the case where the
fold had not run. It is now an assertion that the ladder is empty, which turns a silent append of
two loci onto one ladder into a named panic, and
`building_a_ladder_twice_without_a_fold_between_is_refused` exercises it.

## Validation

All in the container (`./scripts/dev.sh`), on the tree as committed:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::calling::allele_candidates` — **107 passed**, against **93** on the same
  filter with this step stashed.
- `cargo test --all-targets --all-features` — library **5,922 passed, 0 failed, 14 ignored**, and
  every integration and doc target green. The run exits 101 on **one pre-existing failure outside
  this work**: an index-out-of-bounds at `benches/psp_writer_perf.rs:386`, in production's psp
  writer bench, already recorded against the calling loop's block in `PROJECT_STATUS.md`.
- `cargo doc --no-deps` — **26 unresolved-link errors, the same 26 as on the pre-change tree**;
  none of them in these files.

`build_ladder` and every accessor are unreachable outside this module's own tests until
`select_ssr` exists, so the file carries `#![cfg_attr(not(test), expect(dead_code))]` — the
expectation rather than an `allow`, so that the first real caller turns the line into a compile
error and deletes it.
