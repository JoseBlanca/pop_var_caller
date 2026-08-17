# Code review — ng cohort merge, B2: unification into one allele table

*2026-08-17. Six category checklists, each in its own isolated git worktree detached at the
commit under review. Per-category audit trail in
`tmp/review_2026-08-17_ng-cohort-merge-b2/`.*

## 1. Scope

- **What:** the working-tree diff of step B2, as commit `66cb435a` (parent `4f30e334`,
  branch `ng-cohort-merge`). One file changed, `src/ng/run/cohort_merge/build.rs`,
  +879/−15.
- **In scope:** that file — `AlleleTable`, `AlleleLookup`, `alleles_of_sample`,
  `ReadAlleleScratch`, `MemberPlacement::compose_into`, and the tests.
- **Out of scope:** `close.rs` and `mod.rs` (committed earlier, unchanged); `src/pileup/`
  and `src/var_calling/` (frozen production); step B3, deliberately unbuilt.
- **Categories dispatched:** reliability (always; the mutation pass), errors (a new
  assertion class), idiomatic (always), refactor_safety (`project_into` was restructured),
  smells (always), extras (hot path, stable output, diff-matches-intent).

## 2. Verdict

**Approve with changes** — 2 Blockers, 7 Majors, 12 Minors, 6 Nits. Both Blockers are
missing test coverage rather than wrong behaviour; every Major is applied or answered
below.

## 3. Execution status

Run by the author in the container before dispatch, and independently by four agents in
their own worktrees:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. (`--all-targets` is red on
  this branch and was red identically before this work: 49 pre-existing errors in
  `examples/`, `benches/` and other modules' test code.)
- `cargo test --lib ng::run::cohort_merge` — 80 passed at review time.
- `cargo test --lib` — 3,703 passed, 0 failed, 11 ignored.

No finding is labelled "needs verification": every claim below was reproduced by running
something.

## 4. The two Blockers — both untested classes, both now pinned

**B1 — no multi-record fixture contained an indel**, so nothing pinned that a composed
allele is closed on the *reference* it has consumed rather than on its own length. Every
cross-record fixture used substitutions, where the two counts are equal at every step.
Mutating the close to `composed.len()` left all 80 tests green. Fixed by
`an_allele_composed_across_records_closes_on_the_reference_it_consumed_not_its_length`: one
sample's read composes a six-base allele over a five-base locus and another's a four-base
one, and the mutation now fails it (measured — `left: [.., "ATGTG", "ACTTA"]` against
`right: [.., "ATGTGA", "ACTT"]`).

**B2 — no locus had two samples that each hold several records**, so the per-sample reset of
the shared working buffer was unobservable. Deleting it left all 80 green. Fixed by
`two_samples_that_each_hold_several_records_are_derived_independently` — **and the first
version of that fixture did not kill the mutant either**, because the two samples used
different read ids, which is the easy case: with distinct ids the stale sightings compose
the same allele again and the table is unchanged. The ids space is per *file*, so two
samples sharing id 7 is the ordinary case; with that, a carried-over buffer makes sample 1's
read look sighted at four records out of two and its allele is removed instead of built. The
mutation now fails the test (measured).

## 5. The Majors, and what happened to each

| # | Finding | Resolution |
|---|---|---|
| M1 | **The overlap guard only fired when a read was named at both records.** Where two overlapping records of one sample carry *different* reads, every read fails the presence test, and the sample's whole evidence vanishes with no panic — proved with a fixture returning a table holding the reference alone. (errors) | **Applied.** A structural check now runs once per sample in `LocusReferenceBases::over`: records must be disjoint and ascending. Two tests, one with a shared read and one without. |
| M2 | **The fast path's doc claimed an equivalence the code does not have.** A fragment whose mates overlap and disagree is two sequences of one record under one id; with one record both are emitted, with several the read is removed. Proved by running the same data through both branches. (smells) | **Applied as documentation and a test**, and the behaviour kept: with one record there is nothing to compose across, so each mate's sequence is already complete evidence; composing across records needs the fragment as a unit. `a_read_showing_two_things_at_a_sole_record_still_contributes_both`. |
| M3 | **The docs said a qualifying read "is known to have covered the whole locus"**, where what is required is presence at every record *that sample* minted; ground outside those records is reference by construction. (extras) | **Applied.** The module header and `AlleleTable` now defer to `alleles_of_sample`, which states both things the derivation does *not* claim — including the fragment whose mates flank an unsequenced insert. |
| M4 | **Reads removed as evidence were removed silently**, with no count, and unrecoverable downstream — the shape spec §3.3 argues against for failed loci. (extras) | **Applied.** `alleles_of_sample` answers how many it removed and `AlleleTable::reads_removed_as_evidence()` carries the locus total; where it surfaces is B3's or C1's. Pinned. |
| M5 | **The overlap assertion asserts a producer's guarantee and carried no "must become a `RunError`" note**, unlike its two siblings. (errors) | **Applied** on the structural check that replaced it. |
| M6 | **The cross-record derivation composes and hashes one allele per read**, where the one-record branch does it per distinct sequence: measured at 8.1 ms for one locus at 1,000 samples × 300 reads, against 33 µs for the ordinary shape. A prototype composing once per distinct *read pattern* gave the same table and removed 39%. (extras) | **Recorded, not applied.** The dedup also yields the read multiplicity per allele, which B3 needs, so it belongs in the step that will consume it rather than as a change with no test-visible effect here. |
| M7 | **The derivation discards which reads backed which allele**, and `index_of`'s doc recommended a route to B3 that recovers identity but not moments. (extras) | **Doc applied, signature not.** `index_of` now says it answers identity only, and why the moments cannot come from the bytes. What B3's callback needs depends on the attribution rule, which is the owner's to settle. |

## 6. Minors and Nits applied

- The `(ChainId, u32, u32)` triple became `ReadSighting`, whose field order *is* the sort
  order (three categories converged; a swap at the push site compiles).
- The hand-rolled run-grouping loop became `slice::chunk_by` (two categories).
- `AlleleIndex` → `AlleleLookup`, since the architecture already uses *allele index* for the
  position inside the table.
- `intern`'s return is now used: the reference's index is asserted to be `REFERENCE_ALLELE`.
- `MemberPlacement`s are built once per record rather than once per (read, record).
- The working buffer is cleared on entry rather than after the fast path.
- The depth-cap test asserted a distinction its own fixture could not make
  (`reads_discarded_by_cap` is never read here). It now asserts the *sameness* — the table
  with a capped read equals the table without one — which is a property that can fail.
- The `by_read` doc claimed one allocation "for the whole run" where the type's own doc says
  once per locus; corrected.
- The chain-id assertion's comment now says why it cannot fire on the STR path (two STR
  records of one sample can never chain, segments being the reference's partition) and what
  the `num_obs == 0` escape is for — the reviewer measured that clause as dead.

## 7. Declined, with reasons

- **`Arc<[u8]>` to stop each distinct allele's bytes being held twice.** It changes the type
  arch §4 fixes for `CohortObservation::alleles`, and converting back at B3 would hand the
  saving straight back. Recorded at the field; measurable at B3.
- **Hoisting the working buffer above the locus** (`over_with(locus, &mut scratch)`). The
  builder that owns a region is where it belongs (plan C1), and it is documented there.
- **Adding the psp obligations to arch §5.** Three assertions in this file now instruct a
  later step to convert them, and §5 mentions none of them — a real gap, and a **design-doc
  edit this plan does not make silently**. Offered to the owner under the block's `Open:`.

## 8. What's good

- The refactor was proved rather than argued: `project_into` was restored under a second
  name and compared against the new one over **859 shapes** — seven loci including two at
  the coordinate ceiling, crossed with every offset and width and with insertions,
  deletions, exact matches and empty sequences — **zero disagreements**, with both buffers
  pre-dirtied so a missing `clear()` would have shown.
- Determinism was proved three ways, not asserted: the byte-keyed map is touched at four
  places and none iterates; the sort key is unique so stability cannot matter; and the test
  binary run as two separate processes printed byte-identical allele orders (`ahash` pulls
  `getrandom`, so this was worth checking rather than reasoning about).
- The "presence, not amount" probe that caught B0 did **not** reproduce: truncating each
  sequence's ids to the first kills six tests.

## 9. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh cargo test --lib
```
</content>
