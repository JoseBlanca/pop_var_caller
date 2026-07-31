# Code Review: ng locus witness representation — Milestone D (consumers and surfaces)

**Date:** 2026-07-31
**Reviewers:** four category sub-agents, one isolated git worktree each, run **sequentially**
**Scope:** `5baee76..86f60c2` on `ng-pileup-generator` — D3, D4, D5, D6
**Status:** Approve-with-changes — all applied

---

## 1. Scope

- **What:** the four steps that give the witnessed set its constructors, its surfaces and its
  regression anchor — `e46d089` (D3, `from_witnessed_runs`, which replaced the plan's `from_run`
  by owner decision), `e95d9a2` (D4, the dumps' labels and the shared derivation), `c631e3b`
  (D5, the census's hole counters), `86f60c2` (D6, the spliced fixture).
- **In scope:** `src/ng/locus_generation/**`, `examples/ng_*`, `examples/shared/`, the milestone's
  docs, and the one downstream dashboard that parses a moved label.
- **Out of scope:** `src/pileup/` (production, frozen); Milestone E.
- **Categories dispatched:** `reliability` (mutation sweep), `behaviour_safety` (did each step
  change only what it claimed?), `module_structure` + API-fitness (can it carry E and step 7?),
  and **`naming` + documentation accuracy — run as an agent this time**, which the Milestone C
  review recorded as its own protocol gap.

## 2. Verdict

**Approve-with-changes.** No Blocker. One Major that was a false claim in a doc comment, four
Majors of API and documentation fitness, and a tail of stale statements.

The headline is that **the category run by reading last time was the highest-yield one when run
as an agent**: naming and documentation accuracy produced five Majors, every one of them a claim
that was checkable and wrong. The single most valuable finding in the whole round — D5's counters
being documented as "expected to read zero" while the census that prints them reads 400 — came
from an agent running the test with `--nocapture` and reading the output, which is not
proofreading.

**The other pattern worth recording: three findings were guards or claims that C-milestone and
D-milestone work had *invalidated for each other*.** `witness_of`'s clamp and drop filter were
justified by a comment that D3's delegation made false; `from_witnessed_runs`'s start clamp was
dead from birth while its test doc claimed both clamps were load-bearing; and the type doc's
"`witness_of` decides `Complete`" sentence sat three lines above D3's own paragraph saying it no
longer does.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` | clean, at review time and after fixes |
| `cargo test --lib --bins --tests --examples --all-features` | 2,844 → **2,847** passed, 0 failed |
| `ng::locus_generation` | 311 → **312** passed |
| STR dump on the tomato CRAM | byte-identical to the C0 baseline, before and after |

The four agents ran **one at a time**, each in its own worktree, for the reason Milestone C
recorded: four container VMs cannot start at once on this machine. Between them they applied and
reverted **more than fifty mutations** and left no tree dirty.

## 4. Top findings (all applied)

1. **D5's counters were documented as "expected to read zero" and the census that prints them
   reads 400 reads / 528 positions** — the synthetic corpus emits `CigarOp::Skip`, so it contains
   spliced reads. The claim was the stated reason for withholding a floor, so the milestone's own
   deliverable was the one number nothing would notice going to zero.
2. **`from_witnessed_runs` silently dropped a malformed run its sibling constructor rejects** —
   `[(0,4), (9,4)]` gave `Some` with `(0,4)` alone, where `from_half_open_runs` gives `None`. That
   shape is the transposed `(offset, length)` pair, not a state a read can be in.
3. **"The shape is what protects it" was false** — `from_witnessed_runs([(0,9)], LocusLen(6))`,
   the run `from_left(9, LocusLen(6))` builds internally, answers `Complete`. The precondition is
   kept by the call site and nothing else.
4. **The only tests aimed at the generic dump's per-run witness bound never ran in CI**, which
   runs `--lib --tests` and so misses every example suite.
5. **The two research dumps' labels were pinned by nothing** — reverting D4's rename, or labelling
   every partial `complete`, was invisible, because neither binary had a single test.

## 5. Findings

### Major (all applied)

**M1 — `parity.rs`: "expected to read zero" is false on the census that prints it**
- **Category:** behaviour_safety
- **Confidence:** High. `cargo test … every_divergence… -- --nocapture` prints
  `the hole deliverable (D5): 400 reads … blind over 528 positions`.
- **Applied:** the floor D5 withheld — `holed_witness_reads > 0 && hole_positions >=
  holed_witness_reads` — beside the fabrication and stale-widen floors, plus corrected docs in the
  field, both printed lines, the plan and the implementation report. Mutation: disabling the hole
  block now fails the census with `the hole deliverable is 0 reads / 0 positions`, where before it
  failed nothing.

**M2 — `witness.rs`: a malformed run was dropped where the sibling constructor rejects it**
- **Category:** module_structure
- **Confidence:** High. Compiled: `from_half_open_runs([(0,4),(9,4)]) == None` beside
  `from_witnessed_runs([(0,4),(9,4)], len) == Some(Partial { [(0,4)] })`.
- **Applied:** `start >= end` on an *input* run now rejects the whole set; a run merely outside
  the locus is still clamped and dropped, which is a state rather than a mistake. The two cases
  are now different in the code and named apart in the doc. New test
  `a_run_whose_start_is_not_before_its_end_rejects_the_whole_set`; under the old body it reports
  `left: Some(Partial { [(0,4)] }), right: None`.

**M3 — `witness.rs` / arch §1.1: the ruler precondition is not protected by the shape**
- **Category:** naming
- **Confidence:** High. Demonstrated by a compiled test: the run any caller would write for a
  reach — `(0, reach)` — answers `Complete` at `from_witnessed_runs([(0,9)], LocusLen(6))`, which
  is the exact saturating STR case D3's split exists to keep out of `Complete`.
- **Applied:** the claim is replaced by what is true — only the call site keeps the precondition —
  and by what would fix it: a newtype for a run known to be in locus coordinates, recorded against
  spec §6's deferred sealing of `Partial`'s fields.

**M4 — CI never ran the example suites, where D4 put its two new guards**
- **Category:** reliability
- **Confidence:** High. `.github/workflows/ci.yml` ran `cargo test --lib --tests --all-features`;
  examples default to `test = false` for `--tests`.
- **Applied:** `--examples` added, with the reason in the workflow. 33 example tests now run on a
  pull request, at about 0.03 s.

**M5 — the two research dumps' labels were pinned by nothing**
- **Category:** reliability
- **Confidence:** High. `grep -c '#[test]'`: `ng_ssr_cohort_stutter` 0, `ng_ssr_aligner_bakeoff` 0.
  Reverting D4's rename and labelling every partial `complete` both left every suite green.
- **Applied:** one test in each, pinning the four strings, with the reason: the dashboards key on
  this column — one maps it into outcome classes, the other selects `coverage == "complete"` to
  decide which reads carry an exact length, so a partial mislabelled `complete` feeds a censored
  lower bound into a stutter distribution silently.

**M6 — `witness.rs:240`: the type doc still said `witness_of` decides `Complete`**
- **Category:** naming
- **Confidence:** High. Three lines below, D3's own paragraph says the decision moved; the code
  says so at `open_record.rs:255`.
- **Applied:** corrected to name `from_witnessed_runs`.

**M7 — the knife-edge was still 16 in the spec and the arch, and the spec's version was off by one**
- **Category:** naming
- **Confidence:** High. The landed fixture pins the geometry; a 16 bp deletion ends at 44, **two**
  positions short of exon 2 at 46, not one. D6 corrected the plan and left the spec and arch.
- **Applied:** both corrected to 17, with the spec recording where the wrong number came from.

### Minor (all applied)

- **`from_witnessed_runs`'s clamp on a run's *start* was dead code**, and the test doc claimed
  deleting *either* clamp fails an assertion. It survived its own deletion over the whole suite:
  a run starting at or past the locus end has a clamped end no greater than its start and is
  dropped either way. Deleted, with the reason recorded in the test that used to overclaim.
- **`witness_of`'s clamp and drop filter were justified by a comment D3 falsified.** Both are
  still there; what each now buys is different — the filter keeps a wholly-outside run from
  reaching the constructor's new rejection rule as `(x, x)`, and the right-edge clamp keeps the
  `u16::try_from` total. Said so.
- **`witness_of`'s `LocusLen` saturates**, so the delegated comparison differs from the pre-D3 one
  on a 65,536-wide footprint — unreachable, since the config gate caps `max_record_span` at
  `u16::MAX`. The claim is now "for every footprint the config gate admits", which is what was
  measured.
- **D5's two subtractions used `saturating_sub` and `-` two lines apart.** Both are right and the
  asymmetry is real — the hole count's terms come from one set, where span ≥ coverage by
  construction; the fabrication count subtracts from a footprint the witness has never been
  checked against. Written down.
- **The arch's migration table still said `witness_of` keeps the `Complete` short-circuit**; the
  plan's Scope and spec §4 still listed `from_run`; spec §2 still named `rows_observed`. All
  corrected, each pointing at where the decision was made.
- **The generic dump carried D4's comment appended to the pre-D4 one** — the same sentence twice,
  and a forward-looking "D4 owns what this column finally says" six lines above the note that D4
  landed. Merged. Its one surviving `observed:` (in a test comment) is gone, which makes the
  report's "the last user-visible uses" true.
- **`examples/shared/witness_side.rs` said "one body instead of three"** while a fourth site
  remains — `ng_ssr_divergent_reads`, which derives the side from a `(ReadWitness, Vec<u8>)` pair
  with **no locus length in hand**, so it cannot call the shared function at all. Recorded, with
  what folding it in would cost.
- **The implementation report's own mutation numbers were wrong in four places** — "306 passed; 5
  failed" for a mutation that gives 304/7; "C3's three" for five; "same 6" for 7; "7 tests fail"
  reading the *passed* column; "one home instead of two" where D3 moved one home rather than
  merging two; "six trailing lines" for four; one line citation off by four. All corrected against
  re-run output.

### Recorded, not applied — for the owner at Checkpoint D

- **A witness flush at both borders cannot say which border it anchored**, and D4 made
  `(true, _) => Left` a shared, named function — the obvious thing for step 7 to reach for. A
  fourth `WitnessSide` for "flush both, not `Complete`" would force a consumer to decide, and
  would move the STR dump's labels.
- **D5's counters live on a `#[cfg(test)]` census** that only measures where *production's* walker
  also succeeds record-for-record. Getting the RNA-seq number plan E4 wants would need them on
  ng's own `PileupGeneratorCounts`, which any BAM can produce.
- **Nothing relates a run to a slice of `bases`**, which step 7 will need — and it is not merely
  hypothetical: an existing fixture has a `Complete` witness over 5 positions carrying 2 bases, so
  "slice by run length" is wrong under any indel. This is spec §3.2/§6's deferred read-axis
  question, and the review confirms the deferral rather than reopening it.

## 6. What's good

- **The API was checked by writing the code that comes next.** A reviewer wrote Milestone E and
  step 7 consumers against the public API, compiled and ran them; that is what produced three of
  the four structure findings.
- **`witness_of` was proven behaviourally unchanged, not asserted.** A differential of the
  verbatim pre-D3 body against the new one over 98,288 inputs — anchors × widths × all 4,095
  run-masks on 12 positions — agrees on every one, 12,240 of them `Complete`.
- **D6's fixture survived every attempt to fool it.** No reviewer could change what the walk does
  to a spliced read and keep both tests green, and "absent before C3" was confirmed by restoring
  the discard: the record drops from two observations to one for every deletion ≥ 18.
- **The encapsulation is real from inside the crate as well as outside** — bypass attempts from an
  example *and* from `witness`'s own parent module both fail to compile.
- **The D4 rename was not a silent break**: the one downstream consumer maps both spellings behind
  a named assertion, and the other reads only `complete`.

## 7. Commands to re-verify

```
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --lib --bins --tests --examples --all-features
```

Each fix fails under its own mutation:

| test | mutation |
|---|---|
| `every_divergence_from_production_is_one_of_the_six_named_classes` | `if hole_positions > 0` → `if false` |
| `a_run_whose_start_is_not_before_its_end_rejects_the_whole_set` | drop the `start >= end` rejection, clamp instead |
| `the_coverage_column_spells_the_four_cases_…` (both research dumps) | revert either label to `partial_left` |
| `a_run_reaching_past_the_footprint_is_caught_even_when_the_totals_fit` | per-run bound → the enclosing formula |
| `the_census_counts_a_hole_and_the_positions_inside_it` | `span()` → `span() - span()`, or `+= 1` for `num_obs` |
