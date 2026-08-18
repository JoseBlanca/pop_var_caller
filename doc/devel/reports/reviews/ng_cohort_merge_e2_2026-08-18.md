# ng cohort merge — E2 review: overlap resolution

*Review of step E2 of [the plan](../../ng/impl_plan/cohort_merge.md), 2026-08-18. Three
sub-agents, each in its own worktree detached at the step's working-tree diff (`00573f50`,
parent `a915e52a`). Per-category files under `tmp/review_2026-08-18_e2_overlap_resolution/`.*

## 1. Scope

- **What:** `Organiser`'s `owned_through` and `displaced_locus_count` fields,
  `resolve_and_release`, `claim`, `first_base_of`, eleven tests and five corrected E1
  fixtures in [`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs); and
  `refuse_overlapping_ground` with its two tests in
  [`serial.rs`](../../../../src/ng/run/cohort_merge/serial.rs).
- **The review's main charge:** attack the author's claim that the displacement rule is a
  **safety net nothing in a healthy run reaches**, and judge whether the argument written into
  the `Organiser` doc is both sound and legible.
- **Categories, in three agents:** reliability + refactor_safety; idiomatic + errors +
  defaults; naming + smells + module_structure, plus a documentation-truth pass.

## 2. Verdict

**Approve-with-changes.** The argument survived a serious attack. The *evidence* for it did
not: the guard offered as independent proof sat where it could never fail, and nothing
anywhere fed a real builder's output to a real organiser. Both are fixed.

## 3. The safety-net argument held

The reliability agent built a harness running `merge_cohort_through_cache`'s exact loop into a
real `Organiser` and swept **120,000 random layouts** — 1 to 6 samples, 1 to 3 contigs, records
1 to 200 bases, analysed regions with gaps, building widths 1 to 40, `max_cohort_locus_span` 1
to 64, `min_alt_obs` 1 to 4. Of those, 87,650 produced loci and 83,588 produced a failed span:

```
seeds=120000 with_obs=87650 with_failed=83588 displacements=0 overlaps=0
```

It also checked each step against the real functions and found the chain sound: `evict_before`
drops exactly what does not reach the eviction point; `with_observations` hands out the suffix
from the first observation reaching the region's first base; `LocusCloser` never truncates a
chain, since `max_cohort_locus_span` is read only to *label* a closed locus; and a
`ClosedLocus`'s end is a running maximum seeded with its own start, so no locus is inverted.

## 4. Findings

### M1: the disjointness guard was placed where it could not fire

**Category:** reliability. **Confidence:** High. Measured.

`refuse_overlapping_ground` ran at the *end* of `the_outcome_both_drivers_agree_on`, after the
byte-identity comparison had already proved the cached output equal to the oracle's — and the
oracle's loci come from one walk per analysed region, so they are disjoint by construction.
The agent broke the cached driver (eviction at the building region's last base instead of its
first), which really does make it emit overlapping loci: **eight tests failed and none of the
eight was this guard.** Moved above the comparison, the same break trips it by name.

**Fixed:** the call moved above the comparison. Reproduced here — the guard now reports
`the merge produced contig 0:161-338 and contig 0:209-338, which share ground`.

### M2: nothing connected a builder to the organiser

**Category:** reliability. **Confidence:** High.

All eleven new tests build their input from `locus_at`, a `CohortObservation` with an empty
allele table. No test anywhere submitted a real `build_region` outcome to an `Organiser`, so
the counter meant to be the alarm had never been read on a builder's output. The plan's E2
entry asks for the opposite: a fixture with "a wide deletion beginning before a building region
and reaching into it".

**Fixed:** `refuse_displaced_loci` runs the driver's loop — evict, cover, build, submit — into
a real organiser and asserts nothing is displaced. It is called from the shared driver helper,
so six of `serial.rs`'s twenty-eight tests reach it including the two hundred random layouts,
and from the 305–330 deletion fixture at a 20-base width, which is the shape the plan named.
Its discrimination is measured: moving its eviction point one line makes it report *a builder
produced a locus on ground an earlier locus already owned*.

### M3: two `expect`s in the merge, and a std idiom that removes them

**Category:** idiomatic. **Confidence:** High.

`resolve_and_release` peeked two iterators and then took with `.expect("just peeked")`,
unreachable but unannotated, where the same file annotates its other `expect` with
`// PANIC-FREE:`. The agent compiled a `Peekable::next_if` form that consumes the failed spans
in a plain `for` and needs no `expect` at all, preserving the locus-first tie-break.

**Fixed:** taken as offered, with `claim_and_release` split out.

### Mi1–Mi7, all applied

- **An inverted `GenomeRegion` walked the frontier backwards** and silently disarmed the rule,
  with the counter reading zero — probed on the reviewed code. The rest of the file
  (`with_observations`, `cover`, `building_regions_of`) orders its ends and says why; `claim`
  did not. Now `first_base_of` and `last_base_of` order theirs, with a test at both ends.
- **A displaced locus must own nothing**, and a mutation making it extend the frontier survived
  every test — every fixture had the displaced locus ending inside the standing owner's span.
- **The guard's boundary (`<` on inclusive ends) and its contig-first sort key** were both
  unpinned; each mutation survived. Three tests added.
- **`finish` dropped `displaced_locus_count`** — the natural call order lost the number.
  `finish` now returns a `MergeTally` of both counts.
- **"Earlier" was unspecified** where region order and start order disagree. The code means
  *claimed first*; now stated and pinned.
- **The eviction bullet named the cache** where the load-bearing fact is the caller's choice of
  eviction *point* — which E3 rewrites. Reworded to name the discipline and its owner.
- **`Default` was a second construction path** past the longhand `new()`. Now forwards to it.
  `claim` gained `#[must_use]`.

## 5. Documentation truth — four wrong claims out of six

The naming agent checked every claim in the new prose. Four were wrong:

| claim | verdict |
|---|---|
| "the cache evicts only what ends before the region's first base" | **CHECKED-CORRECT** |
| "no fixture reaches the rule through a real merge" | **CHECKED-CORRECT**, and weaker than the truth — no driver builds an organiser at all |
| "`serial`'s **drivers** assert on **every fixture** in that file that the loci are **disjoint and ascending**" | **WRONG on all three counts**: one driver, six of twenty-eight tests, and disjointness only — the guard sorts before it walks |
| "**four** of E1's fixtures did exactly that" | **WRONG** — five |
| `resolve_and_release`'s "both arrive ascending and disjoint — they come from one walk" | **WRONG both ways**: it claims less than the code needs (the two vectors must be jointly non-overlapping, which neither vector says alone) and more than the organiser can know (`submit` takes any outcome from any caller) |
| the `max` comment's "a branch no input can take" | **CHECKED-CORRECT** for well-formed spans — a probe counted zero retreats over the whole suite — and false for an inverted one, which is Mi1 |

All six sentences are rewritten. The wrong ones were all claims about **this work's own
fixtures**, which is where this project's reviews keep finding them.

## 6. What's good

- The `resolve_and_release`/`claim` split was judged right by two agents: a mutating method
  returning "did you get it" is `HashSet::insert`'s shape, and putting the frontier and the
  count in one place is what makes "one rule, no special case" true in code and not only in
  the doc.
- `u64` saturating arithmetic on both counters was judged right — totals, not cursors.
- Two mutations were proved **no-ops rather than survivors**: the merge tie-break `<=` → `<`
  (zero ties over the whole suite) and the frontier `max` (zero retreats). Recorded so nobody
  writes a test against a hazard that does not exist.

## 7. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --lib ng::run::cohort_merge
./scripts/dev.sh bash tmp/mutate_e2.sh          # 21 mutations, all must be killed
```
