# ng calling loop — E1: the input edge, and three joins nothing in the types enforces

**Step:** E1 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — evidence shaping.
**Design authority:** [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §2, §7;
[`arch/read_likelihoods.md`](../../ng/arch/read_likelihoods.md) §2;
[`arch/candidate_alleles.md`](../../ng/arch/candidate_alleles.md) §5.1.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. What landed

**`src/ng/calling/evidence_shaping.rs`** — one merged cohort locus, plus candidate selection's
verdict on it, shaped into the `LocusEvidence` the loop reads. Data-shaping only; no arithmetic.

- `GenericEvidenceScratch` — the per-worker buffers, reset and refilled at every locus, one entry
  per run sample.
- `GenericEvidenceScratch::narrow` and `fill_views`, and `shape_generic_locus` as the one-call
  spelling.
- `shape_ssr_locus` — the repeat-tract half, which is almost nothing to do and says why.

**A module of its own rather than inside arm A's file**, which is a deviation from the plan's
scope line. The shaping is arm-independent: a second arm of the seam
([`calling_bakeoffs.md`](../../ng/impl_plan/calling_bakeoffs.md)'s B and C) reads the same
evidence, and D1's review already noted five helpers stranded inside arm A's file for the same
reason.

## 2. Three joins, and none of them is enforced by the types

**Two per-sample lists in different orders.** `CohortObservation::per_sample` holds only the
samples that *covered* the locus, each naming its own index in the run's sample order;
`LocusEvidence`'s list is one entry per run sample. A sample that covered nothing gets
`GenericSampleEvidence::empty()` — which scores every genotype alike and leaves the prior to
decide, the right answer rather than a special case. And `LocusSelection::unmatched` is parallel
to the **merge's** covering samples, so the uncallable ruling has to travel by that index and land
in the run's.

**⚑ The row order the row builder requires is already given, and the first draft of this report
said the opposite.** It said selection admits candidates in ranking order, so the remapping
permutes and the module must sort. `select_generic` does the reverse: it ranks only to decide
*which* alternatives survive a binding cap, then puts the survivor list back into the merge
table's own index order before admitting — and its own
`the_survivors_of_a_binding_cap_are_admitted_in_the_merge_tables_order` pins that, written after a
reviewer deleted the restoring sort and every test still passed. **So the sort here is insurance,
not description**: `AlleleRemap::admit` constrains only the candidate ids, so a rank-order
selection would be type-legal, and `GenericSampleEvidence::new`'s order check is a `debug_assert`
that a release run would not raise. The two tests that catch the sort's removal both build shapes
no shipped producer can emit, and now say so.

**⚑ The pooled leftover is selection's number, and adding the narrowing's own would double it.**
`GenericObservation::fill_from_supported_alleles`' doc says *"a caller adds this to whatever
selection hands it"* — and that would be wrong. Both quantities are Σ `ln P(error)` over exactly
the rows whose allele has no candidate: selection sums them at selection time (its own test calls
that *"the sum of the dropped rows' own `q_sum`, to the bit"*), and the narrowing sums the same
rows again as it walks them. **This module takes selection's and asserts the narrowing's against
it** — a free cross-check between two independently written walks over one set of rows. The
misleading sentence is in `likelihood/`, whose plan has merged; it is left for its owner rather
than edited here, and the check is what stops a caller acting on it.

Double-counting would not have shown up as a wrong genotype: the leftover is the same number under
every genotype and cancels in the comparison. It would have shown up in the data likelihood, which
emission and the site quality read.

## 3. What it allocates

**Per worker for everything it fills** — the narrowed rows (one buffer per run sample), the two
per-sample maps and the remapping as a slice.

**The one per-locus allocation is the caller's list of views**, and two things make it
unavoidable rather than one. A `GenericLocusSample` borrows the rows the scratch holds, so a buffer
holding both would have to name its own lifetime — which is why the list cannot live *inside* the
scratch. And a `Vec` is **invariant in its element type**, so a caller cannot hold one across two
loci either: the element type names the lifetime of the borrow, and reusing the list holds the
first locus's borrow open into the second. The craft review settled this by compiling four callers
— `E0499` on the one-call spelling, `E0502` on the two-step form. The list is a few dozen bytes per
sample, against a called locus's own output — one `Genotype` per sample — which is the same order.

**A row buffer per sample rather than one flat buffer with spans**, because the narrowing this
calls clears what it fills; spans would need a second, appending spelling of it, and two functions
that must stay identical is the worse trade. **What that costs at a large cohort is now stated
rather than implied**: no buffer is ever shrunk, so the rows held at high-water are Σ over samples
of *that sample's widest locus*, where a flat buffer would hold the largest single locus's total.
The first is never smaller, and the two coincide only if every sample's widest locus is the same
one.

## 4. Tests

**Twenty-two**, 4,694 → 4,716 on the library target. **Six of them are the review's**, each
closing something measured to survive the first suite; ten are `#[should_panic]`.

| test | what it pins |
|---|---|
| `a_sample_that_covered_nothing_keeps_its_place_and_shows_nothing` | the run-order join: a covering-list join by position would give sample 1 sample 2's reads, and both are legal evidence |
| `the_narrowed_rows_are_sorted_on_the_candidate_key_and_not_the_merges` | the insurance case: a remapping that reverses two alternatives, which today's selection cannot produce and the test now says so |
| `two_read_groups_of_one_allele_stay_apart_and_in_group_order` | the second half of the sort key |
| `the_leftover_is_selections_number_and_not_twice_it` | −7.5, not −15.0 |
| `a_leftover_that_disagrees_with_the_narrowings_own_sum_is_refused` | the cross-check between the two walks |
| `the_uncallable_ruling_lands_on_the_run_sample_the_merge_entry_names` | the second covering entry is run sample 2 |
| `the_partial_observations_are_borrowed_from_the_merges_own_rows` | borrowed, not copied — asserted on the pointer |
| `a_sample_that_covered_the_last_locus_but_not_this_one_shows_nothing` | the worker's buffers are reused, and a stale row is the previous locus's reads |
| `a_tract_takes_the_generators_observations_and_sets_no_sample_aside` | spec §5.0.1 on the repeat-tract path |
| `the_one_call_spelling_hands_back_evidence_a_caller_can_use` | that the wrapper is callable at all — its signature borrows the buffers and the views for the evidence's own lifetime, and a wrapper nobody can call is worse than none |
| `a_sample_ruled_uncallable_at_one_locus_is_callable_at_the_next` | the review's Blocker: a stale leftover emits a missing genotype at a locus the sample merely did not cover |
| `the_allele_mapping_is_this_locus_own_and_not_the_last_ones` | the same for the remapping buffer, seen in the candidate ids that come back |
| `a_shifted_covering_sample_with_two_alternatives_is_sorted_like_the_first` | that the sort runs for every sample: restricting it to the first covering entry survived every other fixture |
| `a_worker_shapes_one_locus_after_another_on_one_scratch` | the loop shape E3 has to copy — the view list *inside* the body, because a `Vec` is invariant in its element type |
| ten `#[should_panic]` | the remapping's width, the leftovers' length, a run of no samples, a covering sample past the run, covering samples out of order, the two-walk disagreement, views filled for a different locus, from a never-narrowed scratch, and against a rebuilt covering list |

## 5. Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4716 passed; 0 failed; 14 ignored`. Before E1: **4,694**.
- `cargo test --release --lib ng::calling --all-features` — `670 passed; 0 failed; 3 ignored`.
  Before E1: **648**.
- **The release-held checks: E1 adds eight.** Downgraded all eight to `debug_assert` together and
  re-ran under `--release`: **9 failed**, every one reached. **A ninth was deleted rather than
  kept**: `shape_ssr_locus` restated `LocusEvidence::ssr`'s own emptiness refusal, and the first
  battery showed why that is not free — the test passed with the restatement downgraded, because
  the constructor's check fired instead. A release check no test can reach on its own is one the
  suite cannot keep honest.

## 6. What the review found

**Two agents in worktrees: one on tests and mutation plus the diff's own claims, one on six craft
checklists.** Verdict: **1 Blocker, 5 Majors, 12 Minors** and **9 of 30 claims wrong**. All
applied.

**The Blocker and one Major were the same shape: a buffer whose reset nothing tested.** A sample
ruled uncallable at one locus kept the flag at the next locus it did not cover — and emission
writes that as a missing genotype. Deleting the reset left the whole suite green, because
`fill_views` returns an empty view for a non-covering sample whatever the buffer holds, so every
existing assertion passed. The same held for the allele remapping.

**The three joins came back differently from how they were argued.** The leftover join is right,
and the double-count it guards against is real — the reviewer traced both producers to the same
predicate over the same rows. The sort's justification was backwards (§2). And `fill_views`' guard
compared only the *count* of covering samples, so a locus covered by `{0, 1}` narrowed and one
covered by `{2, 3}` filled passed, pairing one locus's rows with another's partials: it now checks
the locus's own region, and then the join itself at every sample.

**One question the craft agent settled by compiling four callers**: both entry points need a fresh
view list at each locus, declared inside the loop body — not because of this module's signatures
but because a `Vec` is invariant in its element type, so a reused list holds the first locus's
borrow open into the second. Nothing said so, and the single test called the function once,
outside a loop. There is now a test in the shape E3 has to copy.

**Nine wrong claims of thirty**, and the substantive ones were mechanisms: three sentences about
selection's admission order, "`GenericSampleEvidence::new` checks the order" (it is a
`debug_assert`), a `calling_em_loop.md` §7 citation that is about cohort size where the claim's
home is `read_likelihoods.md` §3.3, and a comment describing a stale row as *scored* when it is
merely unread. Three were this report's own: the test count, a claim that one test replaced two
others' coverage, and the number of `#[should_panic]` tests.

## 7. What this step owes

- **E2 gathers the pre-pass's outputs**, including the STR substitution rate D1's tract refusal
  names.
- **E3 is where this shaping meets the driver**, over candidates from `select_generic` rather than
  from a fixture.
- **⚑ `GenericObservation::fill_from_supported_alleles`' doc tells a caller to add its return to
  selection's pool.** That doubles the leftover. The function is in `likelihood/`, and its plan has
  merged; recorded for its owner.
