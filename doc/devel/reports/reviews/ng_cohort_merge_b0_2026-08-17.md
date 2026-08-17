# Code Review: ng cohort merge — B0 (every read is named)

**Date:** 2026-08-17
**Reviewer:** rust-code-review skill (orchestrator, two category sub-agents in isolated worktrees)
**Scope:** step B0 of the cohort-merge plan, captured as `d5558ae5` on branch `ng-cohort-merge`
**Status:** Approve-with-changes

---

### 1. Scope

- **What was reviewed:** one step's diff in the generic locus generator — two one-line mint changes, a
  deleted predicate, the differential's normalisation, and a new test helper.
- **Reviewed against:** `d5558ae52b484d99b6a726ffde20afa65f61b9c6`, branch `ng-cohort-merge`.
- **In-scope files:** [open_record.rs](../../../../src/ng/locus_generation/pileup/open_record.rs),
  [fast_column.rs](../../../../src/ng/locus_generation/pileup/fast_column.rs),
  [parity.rs](../../../../src/ng/locus_generation/pileup/parity.rs),
  [pileup/mod.rs](../../../../src/ng/locus_generation/pileup/mod.rs),
  [generator.rs](../../../../src/ng/locus_generation/pileup/generator.rs),
  [genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs), and the
  [impl report](../implementations/ng_cohort_merge_b0_2026-08-17.md).
- **Out of scope:** `src/ng/run/cohort_merge/`; the psp encoding, which does not exist; the STR path,
  which records no ids and needs none; the pre-existing `--all-targets` clippy failures.
- **Categories dispatched:** `reliability` (the change is one line whose absence is silent — mutation
  testing is the whole review) and `smells` carrying the claim-verification pass (the step's claims are
  about a departure from production, and six of them are checkable).

Both agents ran in their own worktrees, detached at the reviewed commit. The isolation failure of the
B1 review did not repeat.

### 2. Verdict

**Approve-with-changes.** No finding says the mint is wrong. The Blocker and the Major are both about
what the tests could not see — which is the shape this loop keeps producing.

### 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::locus_generation` | `ok. 362 passed; 0 failed; 1 ignored` |
| `cargo test --lib` | `ok. 3679 passed; 0 failed; 11 ignored` (570.94 s) |

Findings labelled "Needs verification": 0. Mutation testing: **4 run, 2 survived, 0 changed no
behaviour**.

### 4. Top findings

**B1 (Blocker) — "every read is named" was pinned only as "some read is named."**
`parity.rs`'s new assertion checked `!chain_ids.is_empty()`. Keeping only the *first* id of every
reference-matching observation, in both mint paths, left the whole library suite green:
`3680 passed; 0 failed`. The mutation was proven to change behaviour — at position 45 of the tiling
fixture the ids print `[1, 2]` unmutated and `[1]` mutated — so the property B0 exists to establish was
untested. The two "revert the change" mutations (withhold ids for bucket 0; restore the fast lane's
`base != ref_base`) were both caught, each by that same assertion.

**M1 (Major) — the new helper accepted a walk that renumbers reads mid-region.**
`assert_same_evidence_up_to_chain_renaming` rebuilt its renaming map per locus, so adding 100 to every
id from position 60 on passed. That is the merge's read-linking broken *inside* a segment — the exact
failure the helper's own doc claims to guard. Measured on the fixture: the renaming is the identity
before the join and `[(1,3),(2,4),(3,5)]` after, so one map per **region** is the true property.

**M2 (Major) — a comment justifying the surviving chain-id equality still stated the deleted rule**,
four lines above the new filter it was supposed to explain, and argued for removing it.

### 5. Minor

- `SequenceObservation::chain_ids`, the **public** field the merge will consume, still described the id
  as "what lets a later step chain observations into a haplotype" — the ruling lived only on the
  crate-private `KeyedObservation::chain_ids`.
- The reference-matching filter is now spelled in five places, and the two new ones dropped the
  complete-observation guard `matches_reference`'s own doc requires: a partial's bases stop where its
  read's witness stopped, so comparing them against the whole locus's reference asks about bases the
  read never saw.
- Two `assert_eq!(len, len)` calls left behind after the helper subsumed them, still carrying the strong
  message of the check they replaced.
- Seven comments and one dangling test name in the copied suite still said "the walker drops REF chain
  ids".

### 6. Nits

A per-locus `reference_bases.clone()` in the normalisation; two undifferentiated `HashMap<ChainId,
ChainId>` in the helper; a doc link pointing at the type rather than the field; the new assertion
dumping every observation on failure; an unreachable `num_obs == 0` disjunct.

### 7. Claim verification — six claims, all correct

1. **The straddling read's two identities — CHECKED-CORRECT, exactly.** A probe over the tiling fixture
   confirms `r1` (25–54) is the only read crossing the cut; the whole walk calls it 1 throughout, the
   split walk 1 through position 50 and **3** from 51.
2. **3,680 → 3,679 — CHECKED-CORRECT.** The diff adds no `#[test]` and removes exactly one, the test of
   the deleted predicate.
3. **`read_agreed_with_reference` had no caller left — CHECKED-CORRECT.**
4. **The parity normalisation is symmetric — CHECKED-CORRECT, with the qualification** that it ran over
   all observations rather than the complete ones (see Minor).
5. **The segment rule says what the report says — CHECKED-CORRECT**, with two refinements: §4.3 states
   the rule for *observations* rather than loci, and §6.2 makes segment independence conditional on
   every sample sharing one segmentation — a second dependency of the same kind.
6. **Production's `allele_index == 0` rule — CHECKED-CORRECT**, and its neighbouring comment carries a
   number nobody had cited: **the ids it drops are ~96.6% of all chain ids on real cohorts**, which is
   the measured size of the cost this step accepts.

### 8. What's good

- Both "revert the change" mutations were caught by the differential against production, which is the
  test that had to keep working across a deliberate departure from it.
- The report raised the straddling-read discovery itself rather than leaving it to the review, and
  named the rule it depends on.
- The correction about the depth cap — that a capped read takes the same branch as an uncovered one, so
  no identities are needed — was made against the author's own earlier claim.

### 9. Commands to re-verify

`./scripts/dev.sh cargo fmt --check`; `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`;
`./scripts/dev.sh cargo test --lib`. Per-category findings are in
`tmp/review_2026-08-17_cohort_merge_b0/`.
