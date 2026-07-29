# ng generic locus generator — the generator, Milestone B: ng's own locus type

**Date:** 2026-07-29 · **Plan:**
[locus_generation_pileup_generator.md](../../ng/impl_plan/locus_generation_pileup_generator.md)
steps B1–B3 · **Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md)
§3, §6, §7, §12 · **Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md) §1.2

Implementation report for Milestone B of plan 3 of 3. Three commits, one per step, plus the
review's fixes.

## 1. Plan

Stop emitting production's `PileupRecord` and emit ng's own `SampleLocusObservations`: rows
whose identity is `(bases, read_coverage, read_group)`, deterministically ordered, with the
two per-record counters a `PileupRecord` has nowhere to put.

## 2. Assumptions

One, and the review priced it: that routing the 44 inherited tests through a back-projection
would cost only the three properties listed on it. It cost three more — §7.

## 3. Changes made

### B1 — a row is a read's bases, its coverage and its lane (`51347ec`)

`ObservationKey { bases, read_coverage, read_group }`, and `finalise` re-derives rows **per
read** instead of reading them off the per-bucket totals. Two of the three parts of a row's
identity are facts about a *read*: a bucket knows its bases, not that one of its reads
witnessed the whole footprint while another saw one position of fourteen, nor that they came
from different lanes.

**Implemented as the identity realised at `finalise`, not as a fold-time bucket key** —
recorded deviation from the step's wording. `read_coverage` cannot be a fold-time key at all
(A4), and arch §1.2 already puts the bucketing at `finalise`. The output was left unchanged,
projected back onto the positional allele list, so the stage-1 differential still proved the
re-derivation faithful — the same "keep the oracle alive across the risky change" that put A0
first.

### B2 — the walk emits ng's own locus type (`86b20b5`)

`finalise` and the walker return `SampleLocusObservations`. Rows sorted; chain ids dropped by
a **per-read** rule replacing production's positional one; `placed_start`,
`to_production_support` and the whole `PileupRecord` conversion deleted.

**The per-read chain-id rule is the substantive change.** Production's `allele_index == 0`
named a unique row only while there was one row per allele; rows now split by coverage and
group, so a reference-matching read can sit in a *partial* row whose bases are a prefix of
the reference bytes and never compare equal to them. "A read that agreed with the reference
across everything it witnessed carries no chain id" is decidable at fold time, stable under
every split, and identical to production's when the rows are one-per-allele.

**The inherited suite goes through one projection, not 67 hand-edits (owner, 2026-07-29).**
Spec §12 sanctions mechanical adaptation; `to_pileup_record` is that adaptation in one
reviewable place. `tests.rs` needed **two** edits and left `copy_fidelity`, which now guards
four files.

### B3 — the reads a cap discarded (`54b8eb6`)

`reads_discarded_by_cap`, resolved at `finalise` against `folded_reads` rather than counted
at the cap, with the truncated read ids plumbed from the walk into the fold.

**The obvious per-record count is the wrong quantity.** Production counts truncated
*positions*, run-wide. A read can be capped at one position of a footprint and survive at
another, and if it folds at all it folds with its **whole window** — so counting truncation
events per record flags records whose support is complete. The quantity is "reads that had
events inside this footprint and were truncated at every position where they did, so folded
nowhere", which is a membership list resolved at the end, not a counter.

## 4. Validation

| | |
|---|---|
| `cargo fmt --check` / `clippy --all-targets --all-features -D warnings` | clean |
| `cargo test --lib` | **2684 passed** (2648 at the start of the plan) |
| `cargo test --lib ng::locus_generation::pileup` | 163 passed |
| `cargo doc --no-deps` | 12 unresolved links, all pre-existing |
| host-native `cargo test --lib` | same result outside the container |

**Soak, host-native, 20,000 cases** — unchanged across all three steps and the fixes:
2,253,903 records compared, 197,380 multi-base, tolerated class **3,073 (0.14 %)**; eviction
census clean over 1,008,679 emitted records; fabrication census **19,703 of 1,010,515
(1.9 %)**.

> **`PVC_PARITY_CASES` does not reach the container.** `scripts/dev.sh` forwards only
> `CARGO_TARGET_DIR` and `HOME`, so a soak invoked through it silently runs the default 1,600
> cases in under a second. Soak runs go host-native.

**Measured cost (review, throwaway probe):** Milestone B is **+15.1 % wall, +24.5 %
allocations**, with **B1 alone +12.5 %**. Not depth-driven — +11.0 % at ~90 reads/case,
because the per-read work is bounded by `max_snp_column_depth`. The costs the code called out
were the wrong ones: the `Vec<u32>`+sort is 2.2 % *and is the determinism guarantee*, the
linear `find` 0 %, hash-keying the rows *worse*. The real one — a `bases.clone()` once per
read rather than once per row — is fixed.

## 5. Deviations from the plan

1. **B1's key is realised at `finalise`, not in the fold** (§3), per arch §1.2.
2. **The inherited suite goes through `to_pileup_record`** rather than being hand-translated
   (owner decision; §7 prices it).
3. **B2 pulled nothing forward from D1.** The differential still compares `PileupRecord`s,
   with ng's side projected back. D1 builds the forward projection and retires this one.
4. **`reads_discarded_by_cap` excludes A5's set** — not in the plan's wording, but required:
   without it a read is counted in both counters (§6).

## 6. The defect this milestone shipped and the review caught

`reads_discarded_by_cap` over-reported. `!folded_reads.contains_key(id)` conflates "the cap
kept it out" with "it folded and then lost its row when its witness turned out
non-contiguous". **240 records in ~506,000** counted one read in *both* per-record counters,
which tell a model different things: one says the support is a subsample of the depth, the
other says a read covered the locus and said nothing usable.

## 7. Review, and what it changed

Five category agents, **each in its own git worktree**
([report](../reviews/ng_locus_generation_pileup_generator_b_2026-07-29.md)): **3 Blockers, 5
Majors, 6 Minors.** The isolation fixed Milestone A's fan-out failure completely — zero
collisions, every result first-hand, nothing needing re-verification.

**All three Blockers were missing assertions, and all three were hidden by the same
mechanism.** `to_pileup_record` merges rows back by bases and discards `region.end`, so:
rows never merging left the suite *and the 20,000-case soak* green; the emitted region could
be anything; and the per-read chain-id rule could be replaced wholesale with production's,
green at 10,000 cases, with both branches heavily exercised (10,826 ids kept on partial rows,
23,644 dropped). **This is the price of deviation 2, now measured rather than assumed.** Three
tests now assert these on ng's own type.

Five Majors: the double-count above; `to_pileup_record` silently dropping a field added to
`ObservedSequence` (the `placed_start` failure mode one type down); the projection being the
exact inverse of `finalise`'s mapping, so a *shared* error cancels; a determinism digest that
could be disarmed; and nine stale comments, four describing machinery B2 deleted.

**Two of the review's findings were about tests I wrote:**

- `placed_left_is_per_record` asserted `num_obs - placed_left == 1` immediately after
  asserting `num_obs == 2` and `placed_left == 1` — arithmetic on the two lines above it,
  introduced *while removing* a real assertion. The "test that cannot fail" pattern again.
- B3's end-to-end cap test asserted `discarded > 0` summed, which survives an off-by-one in
  the collected slice.

**And the first positive control I wrote for the determinism digest was itself inadequate** —
removing a read changes the record count, which the digest includes, so it passed against a
digest that hashed nothing else. It now rewrites every read's MAPQ.

## 8. Checkpoint B — what is open

None blocks Milestone C.

1. **`[profile.release]` sets no `debug-assertions`**, so every `debug_assert` in the walk —
   including A4's coverage-class invariant — is compiled out of **every soak run**. The soak
   proves divergence, not invariants. A decision, not a defect.
2. **The back-projection's cost is measured and time-limited.** D1's forward projection
   removes it; until then the three new native tests are the only cover for the merge, the
   region and the chain-id rule.
3. **+15.1 % wall / +24.5 % allocations for the milestone.** Spec §7: a bad number here is a
   performance problem to solve, not a design to reconsider. **D3 decides**, and this is its
   baseline.
4. **The row sort and `coverage_order` have no test** (review Mi5), and
   `add_contribution`/`subtract_contribution` assign field by field (Mi6).
5. **Nothing in `benches/` drives ng's walker** — both perf measurements this plan has needed
   used throwaway probes. Worth committing one before Milestone C adds a region walk.
6. **Checkpoint A's remaining items** and plan 2's four stand, none blocking.
