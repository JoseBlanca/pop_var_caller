# ng cohort merge — E1 review: ordered release

*Review of step E1 of [the plan](../../ng/impl_plan/cohort_merge.md), 2026-08-18. Three
sub-agents, each in its own worktree detached at the step's working-tree diff
(`6589d440`, parent `51eec029`). Per-category files kept as an audit trail under
`tmp/review_2026-08-18_e1_ordered_release/`.*

## 1. Scope

- **What:** the working-tree diff of E1 — `RegionIndex`, `Organiser`, `MissingRegionResults`
  and sixteen tests in
  [`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs), plus one doc paragraph in
  [`mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs).
- **Out of scope:** `build.rs`, `close.rs`, `serial.rs`, the `ObservationCache`, and the two
  things E1 deliberately defers — overlap resolution (E2) and cache ownership (E3).
- **Categories, in three agents:** reliability + refactor_safety; errors + idiomatic; naming +
  module_structure + smells.

## 2. Verdict

**Approve-with-changes.** One Major found independently by all three agents, one release-build
hazard, and one pair of unpinned counts. All applied — see
[the fix report](../implementations/ng_cohort_merge_e1_fixes_2026-08-18.md).

## 3. Execution status

Every agent built and mutated in its own worktree. `cargo fmt --check` clean,
`cargo clippy --lib --all-features -- -D warnings` clean, `cargo test --lib` `3807 passed;
0 failed; 11 ignored` — the last measured independently by the reliability agent
(`570.13s`) and by the author (`701.00s`).

**Mutation testing: 13 run by the reliability agent, 3 survived, 1 changed no behaviour**, on
top of the 9 the author had run and killed. The three survivors are M1, M3 and M4 below.

## 4. Top 3

- **M1 — a gap at the *tail* of a run finishes `Ok`.** All three agents. The step's own
  contract, missed in the one shape no caller can notice.
- **M2 — the only insert into the reorder map sat inside an `assert!` condition.** One edit to
  `debug_assert!` — which is what the production code this shape is carried from uses — would
  drop every region's outcome in the shipped binary and in no test.
- **M3 — neither count in the refusal is pinned.** Two behaviour-changing mutations of `finish`
  survived the whole suite.

## 5. Findings

### M1: `organise.rs` — a run whose last regions never submit finishes `Ok`

**Categories:** reliability, errors, naming (cross-category note) — three independent agents.
**Confidence:** High.

The organiser learns of a missing index only when a *later* index has already arrived and is
sitting in `held`. A gap at the tail leaves both buffers empty, so `is_finished()` was `true`
and `finish()` returned `Ok`. Probed on the reviewed code: a run that hands out five regions and
receives outcomes for 0 and 1 printed `is_finished: true`, `finish: Ok(1)` — the output stopped
after region 1 and the failed-locus total was short by everything regions 2–4 refused.

Spec §6.3 and the plan's E1 line both say a gap must be an error, never a truncation. The
trailing gap is the half that got through.

**Fix applied:** `is_finished` and `finish` both take `regions_handed_out`.

### M2: `organise.rs` — the map insert was a side effect inside an `assert!` condition

**Category:** reliability. **Confidence:** High.

`submit` performed the module's only insertion into `held` inside the assertion expression.
`debug_assert!` does not evaluate its condition when debug assertions are off, and the analogous
check in `var_calling/vcf_writer.rs:239-244` **is** a `debug_assert!`, while `Cargo.toml`
records that `[profile.release]` deliberately leaves debug assertions off. The agent changed
that one macro name and ran with `CARGO_PROFILE_DEV_DEBUG_ASSERTIONS=false`: `33 passed; 15
failed` — every region's outcome dropped on the floor, invisible in CI and live only in the
shipped binary.

**Fix applied:** check first, then insert — which also keeps the *first* outcome rather than
replacing it on the way to the panic.

### M3: `organise.rs` — neither count in the refusal is pinned

**Category:** reliability. **Confidence:** High. Two survivors, each proved to change behaviour:

- *finish reports one count at a time.* Nothing exercised both counts non-zero at once, and
  that state is reachable.
- *`regions_never_released` counts held loci rather than held regions.* The only test of that
  count submitted three regions of exactly one locus each — a fixture where two different rules
  give the same number.

A third, Minor, survivor: swapping the two field names inside the `#[error(…)]` string changed
the message and passed every test, because nothing asserted the rendered `Display`.

**Fix applied:** one fixture where the counts differ from each other *and* from the loci held,
one test of both counts at once, one test of the rendered message.

### M4: `organise.rs` — the refusal states a cause that is false in one of its two cases

**Category:** errors. **Confidence:** High.

One `#[error]` string covered two different faults. Against a run that released two loci and
never drained them it printed `0 region result(s) never released and 2 locus/loci released but
never drained — a gap stalled the ordered drain`. No gap stalled anything there. The error also
carried counts but never the one identifier that makes it actionable: *which* index never
arrived, which `finish` knows exactly. And `{ 0, 0 }` — "a run that ended short by nothing" —
was constructible and printable.

**Fix applied:** a three-variant enum, `RunEndedShort`, each variant with its own true message,
carrying the first stalled index.

### Mi1: the module-doc argument for keeping the organiser in the cache's file is false

**Category:** module_structure. **Confidence:** High.

Two of its three claims do not hold: `serial.rs:159,163` already calls `cover`/`evict_before`
from a sibling module inside a driver the plan never removes, so those methods can never become
file-private; and `build.rs` cannot reach the cache at all, because `build_region` takes
`&[&[SampleLocusObservations]]`. The reachable narrowing is `pub(super)`, which a file of its
own would get equally. The agent split the file and compiled it (184 tests green).

**Partly applied:** the paragraph is corrected to say what is true. **The split itself is
raised at Checkpoint E rather than taken**, because it moves D1 and D2's code and the
architecture's file tree names `organise.rs`.

### Mi2–Mi8, all applied

- `Organiser`'s `Default`-derived `new` and field-accessing `finish` would silently swallow a
  field E2 or E3 adds — `new` spelled out, `finish` destructures `self`.
- `release_ready` field-accessed `RegionOutcome` against the convention `build.rs` states in
  its own comment — destructured.
- `held` / `released` / `next_expected` are bare participles — `held_outcomes`,
  `released_loci`, `next_expected_region`.
- "ready" meant two different things in `release_ready` and `drain_ready` — the private one is
  `release_regions_in_turn`.
- `submit`'s doc gave a reason that misdescribes arch §5 (which makes `MissingRegionResults` an
  error, not a panic) — rewritten around *when* the fault is caught rather than *what* it is
  about.
- No `#[must_use]` anywhere; a dropped `drain_ready()` is a silently truncated output.
- The `expect` lacked the project's `PANIC-FREE` marker, and the neighbouring counter was
  unguarded where the cursor was checked.

### Nits, applied

`drained` → `drained_regions`; the doc on `a_region_submitted_twice_is_refused` described the
*other* test's fixture; one doc line over the file's width.

## 6. What's good

- **The `from_fn`/`pop_front` drain is right and should not be "simplified".** The errors agent
  compiled `drain(..)`, `mem::take` and a `Vec`-returning variant: all three lose the partial
  drain, and the `Vec` one allocates per call for nothing.
- **`submit` panicking on a hand-out bug is the right line**, and stricter than the production
  shape it carries — `vcf_writer.rs:243` uses `debug_assert!`, which this crate compiles out.
- **All four citations of `var_calling/vcf_writer.rs` were read and are true** (`:152-158`,
  `:168-176`, `:246`, `:256`).
- **The implementation report's own numbers check out** — 168 → 184 measured at both commits,
  and three mutation rows re-verified exactly. One wording snag: "the module" meant
  `ng::run::cohort_merge` in one sentence and `organise.rs`'s test module in the next.

## 7. Out of scope, for the owner

- **Arch §5 still writes `MissingRegionResults { count: usize }`.** E1 now ships a
  three-variant `RunEndedShort`. Either the arch is amended or the deviation stays a
  report-only note — an owner call, raised at Checkpoint E.
- **Splitting `organise.rs`** (Mi1), also for Checkpoint E.
- **`held` is a `BTreeMap` whose operations are all exact** (`remove`, `insert`, `len`), so a
  `HashMap` would serve; the ordering is load-bearing only for `submit`'s guard. Not a defect —
  a note for whoever revisits the container.
- **A `RegionTicket` the organiser mints and `submit` consumes** would make both panics
  unrepresentable and hand `finish` its own count. Needs E3's hand-out loop to attach to;
  weigh it there against the plain `regions_handed_out` argument.

## 8. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh bash tmp/mutate_e1_v2.sh      # 15 mutations, all must be killed
```
