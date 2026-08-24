# ng read likelihoods — E1: the stutter model says *repeats*, because *frame* meant something else here

*Implementation report, 2026-08-24. Branch `ng-calling-likelihoods`, worktree
`../pop_var_caller-calling-likelihoods`. Step E1 of
[`calling_read_likelihoods.md`](../../ng/impl_plan/calling_read_likelihoods.md), Milestone E, on
top of `bb7a41e9`. **This is the first step of the STR path**, and it changes names and documents
only — no arithmetic moves.*

## 1. What it is

The stutter model — how likely it is that a copy of an allele `L` bases long produced a read
showing `L + Δ` bases — was built during the alignment work and named after HipSTR's fields:
`in_up`, `in_down`, `in_geom`, `out_up`, `out_down`, `out_geom`, `equal`. Those names carry *in
frame* and *out of frame*, which
[`read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §1.3 bans outright: *frame* is
borrowed from coding sequence, it says nothing about repeats to a reader who has not met it, and
in this repository it was read as meaning *inside the tract* against *in the flanks*, which is a
different distinction entirely.

E1 renames the seven to the spec's vocabulary, keeps HipSTR's names in the doc comments beside
them, and makes the two document edits spec §7 asks for.

## 2. The seven names

| was | is | what it is |
|---|---|---|
| `equal` | `same_length_share` | the share of reads showing the allele's own length |
| `in_up` | `whole_repeat_longer_share` | the share a whole repeat longer, at any size |
| `in_down` | `whole_repeat_shorter_share` | the share a whole repeat shorter, at any size |
| `in_geom` | `whole_repeat_one_step_share` | **of the reads that slipped by whole repeats**, the share that moved by exactly one |
| `out_up` | `part_repeat_longer_share` | the share longer by part of a repeat |
| `out_down` | `part_repeat_shorter_share` | the share shorter by part of a repeat |
| `out_geom` | `part_repeat_one_step_share` | **of the reads that changed by part of a repeat**, the share that moved by exactly one base |

Both `StutterRates` (public fields) and `StutterModel` (private fields plus accessors) carry them,
and every doc comment names HipSTR's field beside the new name — `in_up_`, `out_geom_`, and so on
— so someone reading the two side by side is not left translating.

**Two names the rename deliberately leaves alone.** `GEOM_MIN` and `GEOM_MAX` keep their spelling:
they are the clamp on a one-step share, `ng::calling::likelihood` re-exports them and asserts they
are the same two numbers ([`mod.rs:151`](../../../../src/ng/calling/likelihood/mod.rs)), and
renaming them reaches outside this step. *Geometric* is not banned vocabulary — only *frame* is —
so the constants' docs now say what they bound (a one-step share) and note that `geom` is HipSTR's
name for it. `MAX_SLIP` also stays, because **E2 replaces it with two cutoffs** and renaming it
twice would be churn.

## 3. What the rename could have got wrong, and what already stops it

A rename is behaviour-preserving or it is a silent wrong answer, and the six rates are three pairs
that a transposition swaps invisibly. Three mappings could have been crossed:

- **longer against shorter** — the direction split, which is the whole asymmetry stutter has;
- **whole-repeat against part-repeat** — the two regimes;
- **a direction share against its one-step share.**

**None of the three needed a new test, because the existing fixtures already separate them**, and
that is a fact about the fixture rather than about luck. `all_distinct()` gives all six rates
**different** values — 0.03, 0.07, 0.95, 0.004, 0.012, 0.8 — for exactly this reason, recorded on
the fixture since the alignment module's own review. On top of it:

- `the_whole_repeat_branch_reproduces_the_published_formula` asserts `probability(+3n)` is
  `0.03 · 0.95 · 0.05^(n−1)` and `probability(−3n)` is `0.07 · …`, at period 3, for n = 1..5. A
  longer/shorter swap moves both sides and fails; so does routing a multiple of the period into
  the part-repeat shares (0.004/0.012).
- `the_part_repeat_branch_reproduces_the_published_formula_in_both_directions` does the same for
  Δ ∈ {1, 2, 4, 5, 7} against 0.004/0.012 and the re-indexed size.
- `the_same_length_share_is_the_remainder_when_the_floor_does_not_bind` pins the derived share
  against `1 − 0.05 − 0.05 − 0.01 − 0.01` and against `probability(0, ·)`.
- `the_two_hipstr_parameter_sets_are_kept_as_matched_rows` pins each constructor's seven values
  one at a time.

So the tests that hold this step are the ones that were already there, and they hold it because
their fixture refuses to give any two rates the same value.

## 4. Eight test names changed; no assertion did

A rename to the spec's vocabulary that leaves the old words in the test names has not done the
job. Eight of the fifteen carried one — *in frame*, *out of frame*, *unit*, *equal* or
*geometric* — and each names a quantity this step gave a different noun:

| was | is |
|---|---|
| `the_in_frame_branch_reproduces_the_published_formula` | `the_whole_repeat_branch_reproduces_the_published_formula` |
| `the_out_of_frame_branch_reproduces_the_published_formula_in_both_directions` | `the_part_repeat_branch_reproduces_the_published_formula_in_both_directions` |
| `out_of_frame_sizes_compress_onto_consecutive_ranks` | `part_repeat_sizes_compress_onto_consecutive_ranks` |
| `every_change_is_in_frame_at_period_one` | `every_change_is_a_whole_repeat_change_at_period_one` |
| `a_single_unit_slip_outweighs_a_larger_one` | `a_single_repeat_slip_outweighs_a_larger_one` |
| `the_cutoff_counts_units_in_frame_and_base_pairs_out_of_frame` | `the_cutoff_counts_repeats_on_one_branch_and_base_pairs_on_the_other` |
| `equal_is_the_remainder_when_the_floor_does_not_bind` | `the_same_length_share_is_the_remainder_when_the_floor_does_not_bind` |
| `the_geometrics_are_held_strictly_inside_zero_and_one` | `the_one_step_shares_are_held_strictly_inside_zero_and_one` |

The other seven test names were already free of it and did not move.

**No assertion changed.** Every numeric literal in the test module's *code* is identical before
and after — checked by extracting them all and diffing the counts, not by reading:
the only numbers that moved are two `§5.2` citations in doc comments that became `§4.2`. What did
change inside the bodies is identifier spelling (`hostile.equal()` → `hostile.same_length_share()`,
the local `one_unit` → `one_repeat`) and two assertion *messages*. So the plan's "existing tests
green unchanged" holds in the sense that matters: same fixtures, same expected values, same
tolerances.

The alignment module's own review report
([`ng_alignment_a3_2026-07-23.md`](../reviews/ng_alignment_a3_2026-07-23.md)) cites six of these
tests by name, **two of them renamed here** —
`the_out_of_frame_branch_reproduces_the_published_formula_in_both_directions` and
`a_single_unit_slip_outweighs_a_larger_one`. It is a historical record and was not edited; the
table above is the translation.

## 5. The two document edits spec §7 asks for

[`alignment.md`](../../ng/spec/alignment.md) §5.2 stated the whole distribution — the regimes, the
seven parameters, the formula, the re-indexing, the placements, the cutoff — and
[`read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §4.2 now states it too, in full, as its
owner. Two spellings of one distribution are two things that can drift apart, so §5.2 **is now the
pointer**: it says who owns the distribution, lists what to go there for, and keeps only what the
owner does not carry.

**Three things stayed**, and each because §4.2 does not carry it:

- **the second silent trap** — the one-step shares are clamped inside (0, 1) and the same-length
  share is floored, which matters most to a consumer *here*, since an aligner prices slips
  relative to no slip and therefore divides by that share;
- **the two HipSTR parameter rows as matched sets**, and the note that an earlier draft of that
  spec paired one number from each;
- **which grain the parameters belong to** — per locus in HipSTR, per read group per stratum in
  ng, and belonging to a sample group either way because stutter depends on library chemistry.

Its *in frame / out of frame* wording is gone, replaced by §1.3's whole-repeat / part-repeat, with
one parenthesis giving the correspondence. The two in-document references to §5.2 that promised
"the parameters, the formulas and the two silent conversion traps … in full" (§4.2's requirement
list) and the reuse map's row were repointed at §4.2 as well.

**The rest of `alignment.md` still says *in frame* / *out of frame*** — its §4.2, for one. Spec §7
asked for §5.2's wording, and widening the edit to a document this step does not otherwise touch
is not E1's to make. The parenthesis in §5.2 gives a reader the mapping.

## 6. What was touched outside `stutter.rs`, and why it is forced

Six aligners derive their slip costs from the model and one example builds a `StutterRates`
literal, so the rename reaches them by construction:

`ssr_anchor_firm.rs`, `ssr_anchor_robust.rs`, `ssr_best_path_unit_slip.rs`, `ssr_unit_robust.rs`,
`ssr_noise_robust.rs`, `ssr_robust_indel.rs`, `examples/ssr_delimiter_comparison.rs`.

In each, the change is the field or accessor spelling, the doc comments that quote those spellings
in a formula, and the local `ln_equal` → `ln_same_length` (it names the renamed quantity and no
other). **Their own prose keeps its words** — an aligner's *unit slip* stays a unit slip — because
that vocabulary is the alignment module's and not this step's to reword. Two module-doc citations
of "spec §5.2" for the distribution were repointed at `read_likelihoods.md` §4.2, since that is
where the formula they quote now lives.

## 7. Validation

All in the container, from this worktree.

| command | result |
|---|---|
| `./scripts/dev.sh cargo test` | **4,354 passed, 0 failed, 14 ignored** in the library target — identical to `bb7a41e9`, measured on the clean tree before the edit |
| `./scripts/dev.sh cargo test --all-features` | same 4,354 / 0 / 14; 4,448 / 0 / 18 across every target |
| `./scripts/dev.sh cargo clippy --lib --all-features --tests -- -D warnings` | exit 0, no warnings |
| `./scripts/dev.sh cargo check --examples --all-features` | exit 0 |
| `./scripts/dev.sh cargo fmt --check` | exit 0 |

`ng::alignment::stutter::tests` holds **15 tests before and after**, all passing, none added and
none removed. `ng::calling::likelihood` holds **162**, untouched — the generic path does not read
this module.

*(`--all-targets` clippy is red on `main` in `examples/ng_duplicated_class_harness.rs` and
`benches/freebayes_bookkeeping.rs`, unrelated to this branch; `--lib --all-features --tests` is the
gate, as the Checkpoint C/D handoff records.)*

## 8. Deviations from the plan

- **Eight test names changed, where the plan said "existing tests green unchanged".** Assertions
  are untouched; the names carried the vocabulary the step exists to remove. §4 above is the
  translation table.
- **`GEOM_MIN` / `GEOM_MAX` and `MAX_SLIP` keep their names** — the first pair because renaming
  them reaches into `ng::calling::likelihood`, the third because E2 replaces it. §2 gives the
  reasoning.
- **Six aligner files and one example changed**, which the plan's scope line did not name. The
  rename forces it: they read the renamed accessors.

Nothing here changes a design decision, so none of it was escalated.

## 9. What E2 and E3 inherit

`MAX_SLIP` still applies one number to two scales, and its doc comment still records what that
costs — E2 splits it into `MAX_WHOLE_REPEAT_SLIP` and `MAX_PART_REPEAT_SLIP` and makes the
discarded mass reported rather than silent. E3 adds `stutter_rates_for(&Slippage)` and the
sums-to-one tripwire. Both now have the vocabulary to be written in.
