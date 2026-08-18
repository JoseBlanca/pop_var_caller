# ng cohort merge — E4 review: the milestone assertion

*Review of step E4 of [the plan](../../ng/impl_plan/cohort_merge.md), 2026-08-18. One
sub-agent, in its own worktree detached at the step's working-tree diff (`081a059c`, parent
`572fc76e`), on the reliability and refactor-safety checklists plus a documentation-truth pass.*

## 1. Verdict

**Approve-with-changes.** One Major: the sweep's own fixture-shape assertions could not fail.
Every measured claim in the step's report checked out; two prose claims did not.

## 2. The Major — assertions that cannot fail

`the_parallel_merge_is_the_oracles_at_every_width_and_count` asserted, once and before the
sweep, that the oracle carries the locus at 305–330 and exactly one refused span at 420–510.
**Both are properties of the fixture, not of any width**: the oracle builds one region per
analysed interval, so those assertions hold whatever the width list is changed to. The shape the
sweep exists *for* — a locus reaching across a building-region boundary — was unguarded.

Measured, by counting the building regions each shape touches:

| width | regions in the interval | regions the deletion touches | regions the refused span touches |
|---|---|---|---|
| 1 | 600 | 26 | 91 |
| 3 | 200 | 9 | 31 |
| 20 | 30 | 2 | 6 |
| 47 | 13 | 2 | 3 |
| 600 | 1 | **1** | 1 |

So four widths divide the deletion and one does not. A width list edited to wide widths, or a
fixture whose deletion moved, would leave both assertions green and the sweep comparing five
undivided merges against an undivided oracle — the milestone reported proved while proving
nothing about division.

**Fixed:** each width asserts its own division count, straight from the divider the driver uses,
and the doc names width 600 as the undivided control rather than implying it straddles.

## 3. How widely a parallel-driver defect is caught — measured

**Nine mutations, none survived, none changed no behaviour.** Six of the nine now reach
`serial.rs`; before E4 the same defects were confined to `parallel.rs`.

| mutation | tests failed | of which in `serial.rs` |
|---|---|---|
| evict at the round's last region | 11 | 3 |
| evict at the first region's last base | 13 | 4 |
| cover only the round's first region | 16 | 6 |
| submit the round's outcomes in reverse | 11 | 3 |
| drop the drained loci | 14 | 5 |
| do not gather the failed spans | 5 | 2 |
| never evict | 1 | 0 |
| a round of one region whatever the count | 1 | 0 |
| drop the round's fifth region | 5 | **0** |

Three things this says. **Two defects are invisible to every output comparison and always will
be** — removing eviction, and ignoring the count in flight — because both leave the answer
byte-identical; the memory table is the sole guard and the module doc says so. **The helper's
sample was blind to a class**: a defect confined to a round's fifth region or later failed no
`serial.rs` test, because its rounds never held more than four; sixteen in flight is now in the
sample. And **E4's new test is not a restatement of E3's sweeps**: not gathering the failed spans
fails it and not the earlier count sweep, whose fixture carries no refused locus.

## 4. Documentation truth

Every measured claim in the step's report checked out and was re-derived independently: 11 tests
failing with 3 in `serial.rs` on the eviction break, against 7 all in `parallel.rs` at the
parent; 224 → 225 tests; and the minted-observations fixture really does route through the
shared helper — proved not by grep but by two mutations failing it by name.

Two prose claims were wrong, and both are corrected:

- **"a width of 1 makes 600 regions, so the round is always full"** — false at 16 in flight:
  600 = 37 × 16 + 8, so the last round holds eight.
- **`serial.rs`'s "four is the smallest that puts several builders in the same round on every
  width this file uses"** — false at width 600, where the analysed stretch is one building
  region and every count gives the same one-builder round.

One more, corrected before the review reported: the deletion was said to reach "across whatever
boundary a width puts near it", which is untrue at width 600.

## 5. Minors, all applied

- **`render` read `RegionOutcome`'s fields rather than destructuring**, and it is what "the same
  answer" means for every comparison in the module — a field added to that type would have
  dropped out of all of them silently. The driver beside it already destructures for exactly
  this reason.
- **The sweep pinned no size**: the oracle produces 50 observations on the shared fixture and
  nothing said so, though the fixture is shared by three files.
- **A round larger than the interval was untested**, though the driver's own note explains why
  the round's buffer is grown rather than reserved on a caller-supplied count. Now a test at
  `usize::MAX` in flight.

## 6. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
```
