# ng candidate alleles — A2: the output vocabulary

*Implementation report, 2026-08-24. Branch `ng-candidate-alleles`, worktree
`../pop_var_caller-candidate-alleles`. Step A2 of
[`candidate_alleles.md`](../../ng/impl_plan/candidate_alleles.md), Milestone A, on top of A1
(`52646056`).*

## 1. Plan

What one call to selection returns, and the buffers it returns it from. Five public types and
one private one, no logic:

- `SelectionVerdict` — `Selected`, `Truncated { dropped }`, `NotPeriodic`, `#[non_exhaustive]`.
- `UnmatchedSupport` — one sample's dropped reads and their error mass.
- `AlleleRemap` — the merge table's indices mapped onto the surviving ids.
- `LocusSelection` — the bundle: table, verdict, per-sample leftover, remapping.
- `SelectionScratch` — the fold's per-worker buffers.
- `AlleleSummary`, private — one allele's fold across samples, which B1 fills and B2 ranks.

## 2. Assumptions and departures

**`AlleleRemap` carries four methods where arch §2.3 names one.** The architecture declares
`candidate_for` and the tuple struct. Nothing could build one: the field is private, so
`all_dropped(table_len)` and `admitted(table_index, candidate)` are what C1 needs, and
`table_len` / `num_admitted` are what a test and the truncation count need. The design is
unchanged — the type still holds `Option<AlleleId>` and never a sentinel — and the additions are
the constructor and the writer that declaration implies.

**`SelectionScratch` carries `new`, `reset_for` and `table_len`.** Arch §2.4 declares the two
buffers and says they are "cleared and refilled per locus"; `reset_for` is that sentence as a
method, and it is here rather than in B1 because it is vocabulary — where the clearing happens is
the difference between a correct fold and one that reads a previous locus's values, and putting
it in the type is what keeps every future caller from re-deciding it.

**`LocusSelection::num_alternatives`** is added because two later readers need it and neither
should subtract one by hand: the genotype prior divides its alternative concentration by this
number (`calling_priors.md` §4), and the truncation count is defined against it.

**`AlleleSummary` carries no leftover fields, and this is the one departure worth the owner's
attention.** Arch §2.4 says the summary holds "the largest share of one sample's reads it took,
how many samples cleared the bar, its cohort read total, **and the reads and mass it would
contribute to the leftover**". The first three are here; the last is not, and `reached_the_bar`
is here instead.

The reason: **the leftover is per sample, not per allele.** `LocusSelection::unmatched` runs
parallel to `CohortObservation::per_sample`, so a cohort-wide mass on a per-allele summary has no
reader — and step C3's own oracle requires the pool to be *the bitwise sum of the dropped rows'
`q_sum`*, which means walking the merge's own rows in the second pass rather than reading a total
off this struct. A per-allele cohort mass would be a field nothing reads, and one that invites
exactly the re-derivation C3 forbids. **This is raised at Checkpoint A as an edit arch §2.4 is
owed, not taken unilaterally.**

**What A2 does not do, deliberately:** nothing constructs a `LocusSelection`. Plan step A2 asks
for the two parallelism invariants "in doc comments **and** asserted where they are built", and
at A2 nothing builds one — the construction sites are C1 (the table and the remapping) and C3
(the leftover). Both invariants are in the doc comments; the assertions land with the code that
can violate them. Recorded here so the obligation is not lost.

## 3. Changes made

One file, [`src/ng/calling/allele_candidates/mod.rs`](../../../../src/ng/calling/allele_candidates/mod.rs),
391 added lines. No other file changed: `pub mod allele_candidates;` went in with A1.

Two shapes carry an invariant the type enforces rather than documents:

- **`AlleleRemap::admitted` refuses a second write to one merge allele.** Without it, an
  off-by-one in C1's admission loop would leave one merge allele silently pointing at another's
  evidence, which is a wrong genotype rather than a panic — the failure the plan isolates C1 for.
- **`AlleleRemap::candidate_for` refuses an out-of-range index** rather than returning `None` for
  it. `None` means *dropped*; an allele the table does not hold is a caller bug, and returning
  `None` for both would let a support row naming allele 7 of a 5-allele table read as an
  ordinary drop.

## 4. Tests added

Eight, all in `allele_candidates::tests`:

| test | what it pins |
|---|---|
| `the_remapping_answers_none_for_a_dropped_allele_and_a_dense_id_for_a_survivor` | the round trip of arch's "Test & bench shape" — five alleles, the **middle** one dropped, every survivor's dense id back and `None` for the hole |
| `the_admitted_ids_are_dense_from_the_reference_with_no_gaps` | the property that makes a candidate id usable as an index into `CandidateAlleles`, written over the survivors rather than as four literals |
| `a_row_naming_an_allele_the_table_does_not_hold_is_a_caller_bug` | the out-of-range read asserts rather than answering `None` |
| `one_merge_allele_cannot_be_admitted_twice` | the double-admission refusal |
| `a_sample_with_nothing_dropped_has_a_zero_pool` | spec §12's "a pool of zero and no branch taken to produce it", including that the zero is positive so a later sum cannot inherit a sign |
| `the_verdict_separates_a_full_list_from_a_truncated_one` | `Truncated`'s payload is compared, not just its discriminant |
| `resetting_the_scratch_leaves_no_value_from_an_earlier_locus` | a written summary does not survive a reset **to a larger table** — the case a `resize` without a `clear` gets wrong |
| `resetting_the_scratch_to_a_smaller_table_shrinks_it` | the opposite direction, which a `clear`-less implementation gets right — both are needed, because either alone admits a wrong implementation |

**The middle hole in the first two is the point.** A remapping that subtracted a constant, or
returned the table index unchanged, agrees with the fixture on alleles 0 and 1 and disagrees on 3
and 4.

## 5. What the review changed

Four agents, each in its own worktree: reliability, naming, idiomatic + errors, and a
**design-fidelity** pass written for this step — does the code match arch §2.2–§2.4 field by
field, and can B1 through D1 be built on it. Full account in
[the review report](../reviews/ng_candidate_alleles_a2_2026-08-24.md). Four distinct defects above
Minor:

- **A Blocker:** `num_alternatives` had no test — **no test built a `LocusSelection` at all** — so
  deleting its `- 1` passed all 13. Every answer moves by one, including 0 → 1 at a locus that
  selected down to the reference alone, which is 27.4% of tomato loci; the genotype prior divides
  its concentration by that number.
- **A Major, found twice independently:** `admitted` guarded the merge index three ways and the
  candidate id not at all, so **two merge alleles could be recorded onto one `AlleleId`** — both
  indices in range, each written once, no bounds check able to see it. The hand-off would then
  score two sequences as one allele. Closed by carrying the admission count and asserting the id
  is the next dense one, which is also order-independent and so holds for the repeat-tract path.
- **A Major:** `reset_for` cleared its fields by name, so **adding a buffer would have compiled**
  — and `arch/candidate_alleles_ssr.md` §5 commits to adding one. Destructuring makes it
  `error[E0027]` instead.
- **A Major, raised by three agents:** `LocusSelection` had four `pub` fields and no constructor,
  so the parallelism invariant its doc comment states had nowhere to be checked. `LocusSelection::new`
  now asserts both halves.

**The judgment call this report flagged was checked and confirmed**, with a third argument the
author had not made: a cohort total is a sum in allele-major order where C3's oracle demands the
per-sample rows' own sum, so the bitwise check would fail *by construction* if the total were the
source. Arch §2.4 and plan step B2 both owe an edit; neither was made here.

Also applied: `admitted` → `admit`, `all_dropped` → `with_all_dropped`,
`num_alternatives` → `alternative_allele_count` (the crate's existing spelling),
`reached_the_bar` removed as a second copy of `samples_clearing_the_bar > 0`, `Clone` dropped from
`SelectionScratch`, the allocation-contract sentence corrected, and three panic tests added.

## 6. Validation

All in the container, on the tree as committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean; and `cargo clippy --lib`
  alone, which is where `dead_code` fires and `--tests` hides it — also clean.
- `cargo doc --lib --no-deps` — 23 unresolved intra-doc links, all pre-existing on `main`, none in
  this file.
- `cargo test --lib allele_candidates` — 20 passed.
- `cargo test --lib` — see the commit message.
- **Mutation re-check on the fixed tree:** removing `- 1` from the alternative count fails 1 test;
  removing the dense-id assertion fails 2.

**One attribute needs explaining.** `AlleleSummary::cleared_the_bar` carries
`#[allow(dead_code)]`: its shipping caller is C1's admission pass, and today only the test calls
it. `#[expect]` would have been the better tool and **cannot express this** — it is unfulfilled in
the test build, where the test is a caller, and satisfied in the library build, where nothing is.
Measured both ways; the `allow` names C1 so a reader knows when it goes.

## 7. Tradeoffs and follow-ups

- **`AlleleRemap` allocates per locus** — a `Box<[Option<AlleleId>]>` of the merge table's length,
  outside the scratch. It is part of the output rather than a working buffer, so it is not the
  per-locus allocation arch §1 rules out; arch §6 already records whether a bitset plus a prefix
  sum would be better as an open item, with "measure before changing".
- **`SelectionScratch` stands alone** and becomes a field of `CallingScratch` when that type
  exists (`calling_loop.md` A1). Nothing about the shape changes when it moves.
- **The cap of 0 or 1 is still representable**, carried over from A1 and raised at Checkpoint A.
  Until it is settled, `select_generic` asserts a cap of at least 2 at step C2.
