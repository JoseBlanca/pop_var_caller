# ng cohort merge — E3 review: the builders in parallel

*Review of step E3 of [the plan](../../ng/impl_plan/cohort_merge.md), 2026-08-18. Two
sub-agents, each in its own worktree detached at the step's working-tree diff (`41a6837c`,
parent `36c0326b`): one on concurrency, reliability and refactor safety, one on the idiomatic,
naming, structure, defaults and documentation checklists.*

## 1. Verdict

**Approve-with-changes.** The concurrency is sound and nobody could break it. What both agents
found instead is that **the one thing E3 adds was invisible to the suite**: the round.

## 2. The concurrency held under attack

The concurrency agent could not make the output depend on how the builders interleave, and
could not make a builder see a cache anyone was drawing forward. It ran:

- **200 repeats** of the same merge — every rendering identical;
- **rayon pools of 1 to 8 threads** — identical to each other and to the default pool;
- **400 random layouts** against the oracle, with the region width drawn 1 to 30 and the
  builder count drawn 1 to 20 — no disagreement, 365 with observations, 371 with failed spans,
  398 with a record straddling an analysed edge;
- **4,096 builders over 30 regions**, and analysed regions shorter than one building region —
  all agreeing.

Its reasoning matched: `with_observations` takes `&self` and reads only `covered_to` and each
sample's held observations; there is no `unsafe`, no interior mutability, no lock; the cover
loop holds `&mut` and must finish before the shared reborrow exists; `collect` blocks until
every builder returns, so the shared borrow is over before the next round's `&mut`; and
`RegionOutcome` has no lifetime parameter, so no window can escape into a result.

## 3. Findings

### M1: the round's size was pinned by an inequality that a doubled round also passes

**Category:** reliability. **Confidence:** High. Measured.

`the_cache_holds_a_whole_round_and_not_one_region` asserted only that eight builders hold more
than one. The agent changed `take(builders)` to `take(builders * 2)` and the whole suite stayed
green — the held-record table moved from 2/4/4/12/29 to 4/4/12/29/62, every count now holding
what the next count up held, and at 16 builders the whole 62-record fixture. **The round is the
design decision this step took**, and its memory claim was covered by an assertion that a 2×,
4× or whole-region round satisfies.

**Fixed:** the test asserts the exact table, `[2, 4, 4, 12, 29]`. All three round-size mutations
now die on it.

### M2: nothing could tell a parallel merge from a sequential one

**Category:** unsafe_concurrency. **Confidence:** High. Measured.

Replacing `.par_iter()` with `.iter()` needed no other change and left the suite green,
including the agent's own six probes. That is by design — the whole claim of the step is that
the answer is unchanged — but it left the one new property unobserved, and the ordering
guarantee resting on an unstated fact (`rayon::slice::Iter` is an `IndexedParallelIterator`,
and `Map` of one is too).

**Fixed at the type level:** the round's builders go through `in_region_order`, whose bound is
`IndexedParallelIterator`. Verified — `.iter()` now fails to compile with *the trait bound
`Map<Iter<'_, GenomeRegion>, …>: IndexedParallelIterator` is not satisfied*, and so would a
`filter` or `flat_map` inserted before the collect. **Whether the builders genuinely occupy
several threads is still untested**, and the module doc now says so rather than leaving the
next reader to look for the test.

### M3: the builder count should be a newtype like the module's other three parameters

**Category:** defaults, naming. **Confidence:** High.

`builders: NonZeroUsize` sat beside three `NonZero*` newtypes with documented defaults. Two
consequences: two same-primitive counts passed positionally, transposable by a reader; and the
parameter that sets this module's resident memory had nowhere to say what a run should use.

**Fixed:** `CohortLocusBuilderRegionsInFlight`, with `one_per_worker_thread()` — a default that
is a rule rather than a constant, because what it should be depends on the machine's cores and
on how much memory the cohort's width leaves. A test pins the rule, since a default ignoring the
pool would pass every other test in the module.

### Mi1–Mi6, all applied

- **`S: Sync` is a constraint on the run's future `ObservationSource`** and nothing said so.
  The agent noted that `organise.rs`'s own `CountingSource` holds `Rc<Cell<usize>>` and so
  cannot be merged in parallel at all. Now stated on the driver, with the escape hatch named.
- **A panicking builder is reachable and undocumented.** Two overlapping records of one sample
  pass the cache — which checks only that starts do not go backwards — and trip `build_region`'s
  disjointness assertion on a rayon worker. The agent confirmed rayon joins every builder before
  re-raising, so the cache is unborrowed when the caller regains it; that it is *not* rewound;
  and that `Organiser::finish` never runs, so the run-ended-short guard is skipped. Documented,
  and now a test.
- **`RegionOutcome` was consumed by field access**, so a field it gains would be silently
  dropped — the failure `organise.rs` guards against twice, with a comment each time. Now
  destructured.
- **`refuse_malformed_analysed_regions` was `pub(super)` in `serial.rs`** though all three
  drivers call it. Moved to `mod.rs`, where the module's other shared vocabulary lives.
- **Three test helpers and the 600-base fixture were duplicated** between the two files whose
  outputs are compared — in the step that created the shared-fixture home for exactly that.
  `source_of`, `width`, `in_flight`, `render`, `refuse_any_difference` and the fixture are now
  shared.
- **The whole-outcome `assert_eq!`s** reintroduced the 26 kB failure message `serial.rs`
  documents avoiding. Replaced by `refuse_any_difference`, which names the first differing entry.

### Nits, applied

`Vec::with_capacity` on a caller-supplied count; `as usize` narrowing a `u64`; `round` →
`regions_in_round`; `merged_in_parallel` → `merge_fixture_in_parallel`; `held_at` →
`records_held_at`; both module docs still calling this step "still to come".

## 4. Documentation truth — one wrong number

| claim | verdict |
|---|---|
| "320 bases at 16 builders on 20-base regions" | **CHECKED-CORRECT** |
| "48,000 open cursors against 3,000" | **CHECKED-CORRECT**; `n` was never defined in the spec |
| **"at most 15 idle builders once"** at a 600-base interval on 20-base regions | **WRONG — it is 2.** 600 ÷ 20 is 30 regions, so at 16 in flight the rounds are 16 and 14. 15 is the bound `builders − 1` over all interval lengths, quoted against an interval that determines the answer exactly |
| the held-record table 2/4/4/12/29 | **CHECKED-CORRECT**, re-measured independently by both agents |
| "no test reaches" the two end assertions | **CHECKED-CORRECT**, with a distinction: every merge evaluates them; no fixture makes either fire |
| "three of them want it" on the shared fixture | **WRONG as written** — the pronoun points at test modules, of which there are two; three *drivers* want it |

## 5. Missing tests, all added

Two intervals on the **same** contig with a gap, and an interval shorter than one region — the
class the two-contig fixture cannot catch, because there the cache refuses the cover outright. A
randomised sweep against the oracle with the count in flight drawn too. A panicking builder. The
overlap half of the shared guard, which the inverted-region fixture does not reach.

## 6. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh bash tmp/mutate_e3.sh    # 13 mutations; 11 killed, 2 documented safety nets
```
