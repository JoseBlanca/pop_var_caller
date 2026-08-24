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
they are the clamp on a one-step share, `ng::calling::likelihood` re-exports them
([`mod.rs:151`](../../../../src/ng/calling/likelihood/mod.rs)), and renaming them reaches outside
this step. The deeper reason the review supplied is better than blast radius: **`GEOM_MIN` bounds
two different quantities** — the one-step-share floor and the derived same-length floor — so an
honest rename splits the constant in two, and that is a decision rather than a substitution. *Geometric* is not banned vocabulary — only *frame* is —
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
- `the_two_hipstr_parameter_sets_are_kept_as_matched_rows` pins **five** of
  `hipstr_shipped`'s seven values and **four** of `hipstr_em_start`'s, one at a time.
  (An earlier draft of this report said "each constructor's seven values"; the step's own
  review measured it. `part_repeat_shorter_share` is asserted for neither constructor and the
  same-length share for neither — the test body is unchanged from `bb7a41e9`, so this is a
  gap in an inherited fixture rather than something the rename introduced.)

So for `probability` — the whole of the arithmetic — the tests that hold this step are the ones
that were already there, and they hold it because their fixture refuses to give any two rates the
same value.

**That was not true of the accessors, and the step's review is what established it.** Every test
that read an *accessor* used a fixture whose longer and shorter shares are equal (`0.05/0.05`,
`0.1/0.1`, `0.01/0.01`), and the seven tests that use `all_distinct()` all go through
`probability`. So making `part_repeat_longer_share()` return the shorter field left the **whole
library green at 4,354 passing tests** while the accessor returned 0.012 where 0.004 is right;
the mirror mutation behaved the same, and so did the whole-repeat pair. The two part-repeat
accessors have no caller outside the module yet — **the genotyping likelihood is the named coming
one**, so a crossed pair would have been waiting for step F2 rather than caught before it.
`every_accessor_returns_its_own_rate` reads all seven on `all_distinct()`; measured here, it
fails on each of those three mutations, and the source was restored from a checksum-verified
copy afterwards.

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

**The rest of `alignment.md` still says *in frame* / *out of frame*** — 14 occurrences, its §4.2
among them. Spec §7 asked for §5.2's wording specifically, and a document-wide vocabulary sweep is
a larger edit than a rename step should carry; the parenthesis in §5.2 gives a reader the mapping
meanwhile. **This is a cost, not a clean boundary**: `alignment.md` now speaks two vocabularies
160 lines apart, and six identifiers in the aligner test modules and the delimiter example
(`Scenario::OutOfFrameIndel`, `an_out_of_frame_change_still_has_a_route`) carry the retired words
by exactly the criterion that renamed eight test names inside `stutter.rs`. **Recorded for the
owner as a step of its own**, not left implicit.

*(The claim is about vocabulary, not about which files were opened: this commit does touch
`alignment.md` outside §5.2, in two hunks — the §4.2 requirement list and the reuse map row — both
of which cite §5.2 and had to be repointed.)*

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

| command | at the rename (`80ecd863`) | after the review's fixes |
|---|---|---|
| `./scripts/dev.sh cargo test` — library target | **4,354 passed, 0 failed, 14 ignored** — identical to `bb7a41e9`, measured on the clean tree before the edit | **4,358 / 0 / 14** — four added tests, one replaced |
| `./scripts/dev.sh cargo test --all-features` — every target | 4,448 / 0 / 18 | 4,452 / 0 / 18 |
| `clippy --lib --all-features --tests -- -D warnings` | exit 0, no warnings | exit 0, no warnings |
| `cargo check --examples --all-features` | exit 0 | exit 0 |
| `cargo fmt --check` | exit 0 | exit 0 |

`ng::alignment::stutter::tests` held **15 tests before and after the rename** — none added, none
removed — and holds **19** after the review: four added
(`every_accessor_returns_its_own_rate`, `sanitizing_maps_each_ill_formed_rate_to_its_documented_value`,
`the_extreme_length_changes_score_zero_rather_than_overflowing`, and the property test
`any_rates_yield_a_distribution_that_is_zero_past_the_cutoff`), and one replaced
(`part_repeat_sizes_compress_onto_consecutive_ranks` → `part_repeat_probabilities_step_down_one_rank_at_a_time`).
`ng::calling::likelihood` holds **162** throughout, untouched: the generic path does not read
this module.

*(One full-suite run failed with a **rustc internal compiler error** in the codegen backend
while building `examples/ng_catalog_window_probe`, a file nothing here touches. It did not
reproduce: the immediate re-run was clean. Recorded because a green log that follows a red one
should say why.)*

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

---

## 10. What the review found, and what it changed

Five agents, one worktree each, against `80ecd863`: naming, refactor safety, reliability,
module structure, and an intent-and-numbers audit. **The rename itself came back clean** — all
seven names match the spec's §4.2 table character for character, every accessor carries HipSTR's
original field name, and no stale identifier survives in `src/ng/` or `examples/`.
Behaviour-preservation was measured rather than argued: **29,787 evaluated cells across three
all-distinct parameter sets, the two named constructor rows and the release sanitizer, maximum
absolute difference 0.0** against both the parent commit and the specification's own formula, with
zero disagreeing cells.

**The one Major was a hole the rename did not open but did put at risk** — the accessor
transposition that the whole library could not see (§3 above). Fixed by
`every_accessor_returns_its_own_rate`, and the fix was mutation-tested three ways.

Also applied:

| what | why |
|---|---|
| `sanitized` destructures `StutterRates` exhaustively, and `new`'s validation array comes from `every_rate` | measured: adding a seventh rate failed to compile at the two constructors **and nowhere else** — the new rate would have been neither validated nor sanitized nor stored |
| the two sanitizers renamed `sanitize_direction_share` / `sanitize_one_step_share`; `the_five_masses_…` → `the_five_shares_…` | undeclared survivors of the retired vocabulary |
| the positional `regime(f64, f64, f64, i64)` helper became a `Regime` literal built at each call site | three same-typed shares in a row is the hazard `StutterRates`'s own doc argues about thirty lines above; a crossed pair now reads as two disagreeing words on one line |
| `ln_same_length` → `ln_same_length_share` at 27 sites | half a name for the quantity the step renamed |
| eight doc formulas rejoined into one code span each | splitting `` `ln(a · b)` `` / `` `− ln(c)` `` across lines made rustdoc emit two adjacent `<code>` elements; a code span carries across a line break |
| five citations in `stutter.rs` repointed | the wholesale "spec §5.2 → `read_likelihoods.md` §4.2" substitution was wrong five times: three of those claims live in the **alignment** spec, one cites §4.2 for the clamps trap that §5.2 deliberately keeps, and one said §4.2 "says to decide" where §4.2 has decided |
| `alignment.md` §5.2 regained the part-repeat-estimator follow-up | **a real loss**: §5.2 was one of only two recorded homes for it, and `read_likelihoods.md` §10 files the item as "Home: unowned, and that is the finding" — so deleting the paragraph orphaned two live citations. The sentence "any comparison involving part-repeat reads inherits that weakness" existed nowhere else |
| `read_likelihoods.md` §4.2 and §7 record that the repointing was made | both still read "the edit is not made here", while §7's own preamble says the three documents "must say the same thing" |
| line anchors in the plan and arch updated | the rename moved every line in `stutter.rs`; **E2's own bullet was sending the next implementer to `stutter.rs:63`**, where `MAX_SLIP` no longer is (it is at 78) |

### The mutation round, and the four holes it found in tests that already passed

The reliability pass ran **50 mutations: 41 killed, 9 survived, and 4 of those nine changed no
behaviour on any input the suite supplies** — so five genuine survivors. It also settled the
question §3 above asks: **all 15 unordered pairs of the six rates, transposed inside the
constructor, are killed**, by between two and eight tests each. The narrowest — the two
part-repeat direction shares — is caught by exactly two, which is what the fixture's 0.004
against 0.012 and the negative-Δ half of the part-repeat formula test were put there for.

What survived, and what now kills it:

| the mutation that lived | why nothing saw it | the test added |
|---|---|---|
| `NaN` in a direction slot sanitized to **1.0** — the *most*-stutter end, where the contract promises the least; and the one-step equivalent sending `NaN` to `GEOM_MAX` | `ill_formed_rates_still_yield_probabilities` reaches that path for all 36 slot-and-value combinations and asserts only that the output is finite and inside `[0, 1]` — **not one of the 36 had its value checked** | `sanitizing_maps_each_ill_formed_rate_to_its_documented_value`, a table of eight inputs against both slot kinds |
| the one-step clamp firing only at the exact endpoints, letting a share of 0.001 through | no fixture supplied a one-step share strictly inside `(0, 0.01)` or `(0.99, 1)`. It matters: 0.001 reaches an aligner as `ln(1 − 0.001) ≈ −0.001`, an almost-free slip extension | the same table's last two rows |
| `hipstr_em_start`'s whole-repeat shorter share edited 0.1 → 0.2, and its part-repeat shorter share 0.01 → 0.02 | the matched-rows test named five of one row's six values and four of the other's | that test extended to **all twelve** |
| a size computed by negation rather than `unsigned_abs`, which panics on `i64::MIN` | the widest change any test passed was ±60 | `the_extreme_length_changes_score_zero_rather_than_overflowing` |

**And one test was doing nothing.** `part_repeat_sizes_compress_onto_consecutive_ranks` built no
`StutterModel` and called nothing in the module — it recomputed `bp_diff - bp_diff / period` in
its own body and compared the answer with `[1, 2, 3, 4, 5]`, which is a statement about Rust's
truncating division. Across all 50 mutations **it never once killed one**, including the two that
corrupt the very re-indexing it was named for. Replaced by
`part_repeat_probabilities_step_down_one_rank_at_a_time`, which asserts the same property as a
ratio between consecutive values *through* `probability`, so it re-spells nothing the
implementation says.

**A property test now covers the rate space**, because the sweep it joins walks three hand-built
models — three corners of six dimensions — and that stops being enough at step E3, where
`stutter_rates_for` starts handing this type **fitted** numbers. No fitted row has ever been
through this suite.

Every one of the six mutations above was re-run here against the added tests: each is killed, and
the source was restored from a checksum-verified copy after each. *(One batch run hit its timeout
mid-loop and left a mutation on disk — caught by the checksum check, restored, and re-run one
mutation at a time. It is exactly the failure mode the project's own commit discipline warns
about, and the reason that check exists.)*

**Two claims in this report were wrong and are corrected above**: the matched-rows test pins five
and four values rather than seven and seven, and §3's "the existing fixtures already separate
them" held for `probability` but not for the accessors.

**One inherited doc sentence was wrong about its own code.** `StutterModel::new` promised that "a
non-finite rate becomes `0` (no stutter), the no-information end of the scale". That holds for
`NaN` and `−∞`; **`+∞` clamps to 1.0** in a direction slot and to `GEOM_MAX` in a one-step slot,
which is the opposite end. The wording came from `bb7a41e9` unchanged, and nothing would have
caught it — writing the table test is what exposed it. Corrected to say what the code does and
why the two cases differ: `NaN` is missing information, an overshoot is a magnitude past the
range.

## 11. Left for the owner, not decided here

Five things the review raised that reach past a rename step, listed with a recommendation each in
the handover rather than settled in this commit.

**The one with a live defect behind it: three of the six aligners cannot tell an expansion from a
contraction.** Swapping `open_expansion` with `open_contraction` in
[`ssr_anchor_robust.rs`](../../../../src/ng/alignment/ssr_anchor_robust.rs) survives the whole
`ng::alignment` suite, and it is **not** a no-op there — measured on that file's own
contraction-biased fixture, the two costs are −3.489 and −2.641, so the swap moves a slip by
0.847 nats. In [`ssr_noise_robust.rs`](../../../../src/ng/alignment/ssr_noise_robust.rs) and
[`ssr_robust_indel.rs`](../../../../src/ng/alignment/ssr_robust_indel.rs) it is worse in kind: every
fixture there is `hipstr_shipped()`, whose two whole-repeat shares are both 0.05, so both opens
are −2.9191921964316565 clean **and** swapped, bit for bit — **no transposition in those files can
ever be caught**. That is the same blindness `all_distinct()` was created to remove, in the module
that fixture lives in. It is pre-existing, it is not this step's to fix, and the fix is small:
give those three files an asymmetric fixture, as `ssr_anchor_firm.rs` already has.

The other four: renaming `StutterRates`, whose noun says *rates* while all six fields say *share*;
a vocabulary sweep of `alignment.md`'s remaining 14 occurrences and six aligner identifiers;
unifying the six aligners' near-identical `SlipCosts::from_model`; and the tautological
`assert_eq!` at `ng::calling::likelihood`'s `mod.rs:2129-2130`, which compares a re-export to
itself and cannot fail.
