# Code Review: ng generic locus generator — the generator, Milestone B

**Date:** 2026-07-29
**Reviewer:** rust-code-review skill (orchestrator) — 5 category sub-agents, **one git worktree each**
**Scope:** the Milestone B diff, `942d59f..54b8eb6`
**Status:** Request-changes → **fixes applied**, see §11

---

## 1. Scope

- **What was reviewed:** a branch diff — B1 `51347ec`, B2 `86b20b5`, B3 `54b8eb6`.
- **In-scope files:** [open_record.rs](../../../../src/ng/locus_generation/pileup/open_record.rs),
  [genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs),
  [mod.rs](../../../../src/ng/locus_generation/pileup/mod.rs),
  [parity.rs](../../../../src/ng/locus_generation/pileup/parity.rs),
  [tests.rs](../../../../src/ng/locus_generation/pileup/tests.rs),
  [copy_fidelity.rs](../../../../src/ng/locus_generation/pileup/copy_fidelity.rs).
- **Out of scope:** `src/pileup/**` (production, frozen); the four still-verbatim copies for
  *content* (their guard was checked); Milestones C and D.
- **Categories:** reliability, refactor_safety, errors, smells+naming, extras.

**Method note — the Milestone A flaw was fixed and the fix worked.** That review ran nine
agents against one shared worktree and they overwrote each other's mutations; five had to
retreat to private worktrees mid-run, and both Blockers needed serial re-verification. This
review gave every agent `isolation: worktree` from the start. Every agent mutation-tested
freely, every result is first-hand, and **nothing required re-verification**. Do this from
now on.

## 2. Verdict

**Request-changes** — three Blockers, all missing *assertions* rather than wrong code, and
all three hidden by the same mechanism. One Major was a genuine defect. All applied.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib` | `2679 passed` at review time; **2684** after fixes |
| `cargo doc --no-deps` | 12 unresolved links, all pre-existing (a 13th, introduced by the fixes, was caught and removed) |
| `cargo audit` | **not run** — not installed in the container; no dependency changes |
| soak, host-native, 20,000 cases | 2,253,903 records; tolerated class 3,073 (0.14 %) |

Findings labelled "Needs verification": **0**.

## 4. Open questions and assumptions

1. **`[profile.release]` does not set `debug-assertions`.** Every `debug_assert` in the walk
   — including the coverage-class invariant added in Milestone A — is compiled out of every
   soak run. The soak proves divergence, not invariants. Is that intended?
2. **The back-projection's cost is now measured.** Three properties of the emitted type were
   invisible through it (§6, B1–B3). D1's forward projection removes the need for it; until
   then the three new native tests are the only cover.
3. **Milestone B costs +15.1 % wall / +24.5 % allocations, B1 alone +12.5 %.** Spec §7 says a
   bad number is a performance problem to solve, not a design to reconsider. D3 decides.

## 5. Top 3 priorities

1. **M1** — `reads_discarded_by_cap` double-counted a read on 240 records in ~506,000.
2. **B1** — rows never merging left the whole suite *and the 20,000-case soak* green.
3. **B3** — the per-read chain-id rule could be replaced wholesale with production's, green.

## 6. Findings

### Blocker

**B1: [open_record.rs:566](../../../../src/ng/locus_generation/pileup/open_record.rs#L566) —
`observation_rows` never merging two reads that share a row identity is undetectable.**
*Categories: reliability.* One row per read, every `num_obs == 1`, all 158 tests green — and
still green at 20,000 soak cases. The three B1 fixtures test only the *splitting* half, and
the differential is blind by construction because `to_pileup_record` merges rows back by
bases before the walkers are compared. The projection undoes exactly the defect.
*Fixed:* `rows_merge_the_reads_that_share_an_identity`, on ng's own type.

**B2: [open_record.rs:696](../../../../src/ng/locus_generation/pileup/open_record.rs#L696) —
the emitted `region` is asserted nowhere.** *Categories: reliability, errors (convergent).*
Dropping the `saturating_sub(1)` is green; replacing the end with `Position(0)` is *also*
green. `to_pileup_record` carries the contig and the start and discards the end, so the 44
inherited tests and the whole oracle are projected through a function that throws the field
away. `region.len()` sizes `num_obs_along_locus`'s depth vector and defines flush-right.
The value itself was verified **correct**; only unpinned.
*Fixed:* `the_emitted_region_covers_the_footprint_inclusively`, one-base and widened.

**B3: [open_record.rs:536](../../../../src/ng/locus_generation/pileup/open_record.rs#L536) —
the per-read chain-id rule's answer is asserted nowhere.** *Categories: errors.* Replacing
the whole body with production's positional `return state.allele_index == 0` → 158/158 green
at 10,000 cases. Making every partial row silently lose its ids → green again. Both branches
are heavily exercised (10,826 ids kept on partial rows, 23,644 dropped), so this is a missing
assertion, not a missing fixture: `classify_record` compares ids by equality only on
`generate_complete`, where ng's rule and production's coincide by construction, and the
fabrication census — the only test that sees partial rows — never looks at ids.
*Fixed:* `only_the_reads_that_departed_from_the_reference_carry_a_chain_id`.

### Major

**M1: [open_record.rs:604](../../../../src/ng/locus_generation/pileup/open_record.rs#L604) —
`reads_discarded_by_cap` over-reports. A real defect.** *Categories: errors.*
`!folded_reads.contains_key(id)` conflates "the cap kept it out" with "its witness was
non-contiguous", so a read that survived the cap, folded, then lost its row to a hole is
counted in **both** `reads_without_observation` and `reads_discarded_by_cap`. Reproduced on
**240 records of ~506,000**. *Fixed* by excluding A5's set, and pinned by
`a_read_that_lost_its_row_to_a_hole_is_not_also_counted_as_capped`.

**M2: `to_pileup_record` silently drops a field added to `ObservedSequence`.**
*Categories: refactor_safety.* `finalise`'s destructure catches a field going *in*; nothing
caught it never coming *out* — the `placed_start` failure mode one type down, on the type
most likely to gain fields. *Fixed:* exhaustive destructure, with the two deliberately
discarded fields named.

**M3: the projection is the exact inverse of `finalise`'s mapping, so a shared error
cancels.** *Categories: refactor_safety.* Swapping `num_fwd`/`placed_left` in **both** leaves
158 green with every emitted locus wrong; swapping in `finalise` alone fails 6. The suite
sees the round trip, not the product. *Mitigated* by B1–B3's native tests, which assert the
emitted values without passing through the projection; fully removed by D1.

**M4: the determinism digest can be disarmed silently.** *Categories: reliability, extras.*
Making it hash the record *count* instead of the locus, with a real hash-order bug live,
leaves 158 green; `records > 1000` counts loci, not sensitivity. Separately, an inherited
`PVC_DETERMINISM_CHILD` makes the top-level run take the child branch and assert nothing
while still reporting `ok`. *Fixed:* a positive control plus a parent-marker guard — and the
**first control was itself inadequate** (removing a read changes the record count, which the
digest includes), so it now rewrites every read's MAPQ, moving the evidence while leaving
which records exist alone.

**M5: comment truth — nine stale comments, four describing machinery B2 deleted.**
*Categories: smells.* The worst is structural: a **76-line doc block containing three
functions' docs concatenated** (with a surviving `/// (continued)` marker) attached to
`read_agreed_with_reference`, leaving `finalise` and `observation_rows` — the two
load-bearing functions — with **no doc at all**, and `observation_rows`' "**this line** is the
determinism guarantee" pointing sixty lines from the sort it names. Also: `AlleleSupportStats`
still described the deleted `PileupRecord` boundary; `RecordWitness` said "two counts … only
until B2"; `FoldedReadState::read_group` said "B1 still owes …" when B1 had done it; the copy
inventory was one file stale in three places. *All fixed.*

### Minor

- **Mi1 — a test that could not fail, introduced by B2 while removing a real assertion.**
  `placed_left_is_per_record` asserted `num_obs - placed_left == 1` immediately after
  asserting `num_obs == 2` and `placed_left == 1`. **Fixed** to check what the deleted
  `placed_start` half actually pinned: that the counter is per record.
- **Mi2 — B3's end-to-end cap test asserted `discarded > 0` summed**, which survives an
  off-by-one in the collected slice (`contributors[cap + 1..]` reports three of four).
  **Fixed** to a per-locus `== 4`; mutation-proven by the reviewer.
- **Mi3 — the hot-path costs named in the code were the wrong ones.** Measured: the
  `Vec<u32>`+sort is 2.2 % *and load-bearing*, the linear `find` 0 %, hash-keying the rows
  *worse*. The real cost is a `bases.clone()` once per **read** rather than once per **row**.
  **Fixed** by borrow-comparing the identity.
- **Mi4 — naming.** `as_pileup_record` allocates, copies and is lossy, which `as_` promises
  the opposite of; `note_no_observation` acquired a second caller in B3 that its name made
  false. **Fixed:** `to_pileup_record`, `note_read_id_once`, plus `record_evidence`,
  `drive_production`, `truncated_read_ids_buf`.
- **Mi5 — the row sort and `coverage_order` have no test**, though B2 is honest that the sort
  is not the determinism guarantee. **Recorded**, not fixed: the property it *is* kept for
  (order follows identity, not arrival) is in-process observable and cheap to add at D1.
- **Mi6 — `add_contribution`/`subtract_contribution` assign field by field**, so a new stats
  field compiles there and never accumulates. **Recorded.**

### Nits

`finalise` computes `coverage_of` a second time per read (already computed in
`observation_rows`) — drift surface as well as cost; `RecordWitness`'s complete/partial split
is read by nothing outside this file's tests.

## 7. Out of scope observations

- `[profile.release]` sets no `debug-assertions` (see §4.1).
- Nothing in `benches/` drives ng's walker; both perf measurements this plan has needed were
  taken with throwaway probes.
- `cargo doc --no-deps` has 12 pre-existing unresolved links in `ssr`/`em`/`sfs`.

## 8. Missing tests to add now

Added: the three Blocker tests, the double-count test, the determinism positive control, the
per-locus cap count, and the replacement for Mi1. Still owed: Mi5's sort test, and a direct
test of `to_pileup_record` itself (D1 retires it, so its value is time-limited).

## 9. What's good

- **Worktree isolation made the whole review first-hand.** Five agents, zero collisions, no
  finding needing re-verification — against Milestone A's five-of-nine retreat.
- **Every mutation-table row in the B1/B2/B3 commit messages reproduced**, including the two
  recorded as *uncaught* and the exact figures (`left: 4 / right: 1`, `15.3 %`).
- **`read_agreed_with_reference` is correct on every reachable input** — all three defensive
  clamps proven unreachable across ~1.6 M records.
- **The census ceiling added at the Milestone A review fired on its first real occasion**,
  turning a one-sided `placed_start` projection into a 15.3 % failure instead of a silently
  restated measurement.
- **`copy_fidelity` still guards what it claims**: the four named files are byte-identical
  bar the sanctioned additions, and the four released ones genuinely diverge.

## 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
PVC_PARITY_CASES=5000 cargo test --release --lib ng::locus_generation::pileup::parity   # host-native
```

## 11. Fixes applied

All three Blockers, M1–M5 and Mi1–Mi4, in `f810547`. Each fix that closes a hole was verified
by re-running the mutation that opened it:

| finding | verification |
|---|---|
| B1 | rows-never-merge now fails `rows_merge_the_reads_that_share_an_identity` alone |
| B2 | the region off-by-one now fails `the_emitted_region_covers_the_footprint_inclusively` alone |
| B3 | the positional rule restored now fails `only_the_reads_..._carry_a_chain_id` alone |
| M1 | removing the A5 exclusion now fails `a_read_that_lost_its_row_to_a_hole_...` |
| M4 | a digest hashing only the record count now fails the positive control |

Mi5, Mi6 and the nits are recorded with reasons. Suite 158 → **163** in the module, **2684**
overall; `cargo doc` back to its 12 pre-existing failures.
