# Code review — ng calling loop E1: the input edge

**Scope:** the working-tree diff of step E1 of
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md), on top of `b10864b7` — one new module,
`src/ng/calling/evidence_shaping.rs`, and one `pub mod` line.
**Date:** 2026-08-26. **Verdict: request changes** — **1 Blocker, 5 Majors, 12 Minors, 5 Nits**,
and **9 of 30 of the diff's own claims wrong**. All applied; see
[the fix report](fixes_applied_2026-08-26_e1.md).

**Two agents, each in its own worktree**: one on reliability and on step 8a's re-derivation
(**23 mutations run, 8 survived, 4 changed no behaviour**), one on six craft checklists, every
proposed refactor compiled in its own tree.

---

## Blocker — a reset nothing tested, and the failure is a missing genotype

`narrow` resets each sample's leftover; nothing checked that it does. **A sample the candidate step
ruled uncallable at one locus kept `genotype_must_be_missing` at the next locus it did not cover**
— and `fill_views` writes that flag for *every* sample, covering or not, so emission would write a
missing `GT` at a locus where the sample merely showed nothing. Deleting the reset left the whole
`ng::calling` suite green, because every existing assertion looked at the *view*, and a
non-covering sample's view is empty whether the buffer was reset or not. **Fixed** by
`a_sample_ruled_uncallable_at_one_locus_is_callable_at_the_next`, and by the reset being a
destructured field-by-field reset that a sixth field cannot slip past.

## Majors

**M1 — the sort's justification stated the opposite of what selection does.** The module said
selection admits candidates in ranking order, so the remapping permutes and the rows must be
re-sorted. `select_generic` puts its ranked survivor list *back* into the merge table's index
order before admitting, and its own
`the_survivors_of_a_binding_cap_are_admitted_in_the_merge_tables_order` pins that — written after
a reviewer deleted the restoring sort and every test still passed. **The sort stays**, because
`AlleleRemap::admit` constrains only the candidate ids and a rank-order selection would be
type-legal, and because `GenericSampleEvidence::new`'s order check is a `debug_assert` a release
run would not raise. **Only the three sentences explaining it were wrong**, and the two tests that
catch the sort's removal now say that their fixtures are shapes no shipped producer can emit.

**M2 — the same untested-reset shape for the allele remapping.** Fixed, with a test that sees it
in the candidate ids that come back.

**M3 — `fill_views`' cross-locus guard compared only the *count* of covering samples.** A locus
covered by `{0, 1}` narrowed and one covered by `{2, 3}` filled passed the guard and paired one
locus's rows with another's partials — both legal evidence, and neither the caller's. It now checks
the locus's own **region**, which also refuses a scratch never narrowed at all, and then the join
itself at every sample: the covering entry a row was narrowed under must still name that run
sample.

**M4 — no fixture had a shifted run index *and* a second alternative**, so restricting the sort to
the first covering entry survived the suite. Fixed by
`a_shifted_covering_sample_with_two_alternatives_is_sorted_like_the_first`.

**M5 — the row-buffer reset's comment named a hazard that cannot occur**, and the test written for
it could not fail on it: `fill_views` never reads a non-covering sample's rows. The comment now
says the buffer is *unread* rather than *scored*, and a test-only `rows_left_for` is what lets the
test see the thing it claims.

## The three joins, as they turned out

**The leftover join is right, and the double-count it guards against is real.** The reviewer
traced both producers to the same predicate over the same rows: selection's pool and the
narrowing's return are Σ `ln P(error)` over the rows whose allele has no candidate. Adding them —
which `fill_from_supported_alleles`' own doc tells a caller to do — would double every covering
sample's pooled error mass.

**Both entry points need a fresh view list at each locus**, settled by compiling four callers:
`E0499` on the one-call spelling with the list hoisted, `E0502` on the two-step form. Not this
module's signatures but `Vec`'s invariance in its element type. Nothing said so, and the single
test called the function once, outside a loop; there is now a test in the shape E3 must copy.

**The run-order join was right and is now checked twice** — at the shaping, and again at every
sample when the views are filled.

## Minors, applied

The module renamed `evidence` → `evidence_shaping` (it owns no evidence *type*) and the scratch
`GenericEvidenceShaping` → `GenericEvidenceScratch` (the crate's word for that role); the three
parallel per-sample `Vec`s folded into one `Vec<NarrowedRunSample>`, inside the module whose stated
job is removing positional joins; the two-walk tolerance made relative and named, with the measured
association error behind it — 200 dropped rows differ by 1.4e-12, 600 by 1.1e-11, 2,000 at −4.0e6
by 3.7e-9, which the old absolute `1e-9` would have reported as a selection defect; a
`calling_em_loop.md` §7 citation corrected to `read_likelihoods.md` §3.3; a dead `truncate` after
`resize_with` deleted; `#[must_use]` on both entry points; a `run_sample_count == 0` refusal where
the count arrives rather than two calls later.

## Nine wrong claims of thirty

Three sentences about selection's admission order; "`GenericSampleEvidence::new` checks the order"
(it is a `debug_assert`); the §7 citation; a comment calling a stale row *scored* when it is
unread; and three of this report's own — the test count (4,710, not 4,709), a claim that one test
replaced two others' coverage, and the number of `#[should_panic]` tests.

## Out of scope

**`GenericObservation::fill_from_supported_alleles`' doc tells a caller to add its return to
selection's pool**, which doubles the leftover. It is in `likelihood/`, whose plan has merged.
Raised for its owner.

## Verification

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4716 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `670 passed; 0 failed; 3 ignored`.
- The eight release-held checks downgraded together: **9 failures**, every check reached.
