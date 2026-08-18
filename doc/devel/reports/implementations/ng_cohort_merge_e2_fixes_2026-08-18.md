# ng cohort merge — E2: applying the review's fixes

*Fix-application report, 2026-08-18, against
[the E2 review](../reviews/ng_cohort_merge_e2_2026-08-18.md). Every finding was applied.*

## 1. The two that mattered

**The disjointness guard could not fire.** `refuse_overlapping_ground` ran at the end of
`the_outcome_both_drivers_agree_on`, after the byte-identity comparison had already proved the
cached output equal to the oracle's — and the oracle's loci come from one walk per analysed
region, so they are disjoint by construction. The reviewer broke the cached driver, which
really does make it emit overlapping loci: eight tests failed and **none of the eight was this
guard**. The call now runs *before* the comparison, and reproducing the same break here trips
it by name: `the merge produced contig 0:161-338 and contig 0:209-338, which share ground`.

**Nothing fed a builder's output to an organiser.** Every one of E2's organiser tests is
fabricated — that is deliberate, since the claim is that no builder can produce an overlapping
pair — but it left the claim itself untested end to end, and the plan's E2 entry asks for a
fixture with a wide deletion beginning before a building region and reaching into it.
`refuse_displaced_loci` now runs the driver's own loop — evict at the building region's first
base, cover, build, submit — into a real `Organiser`, and asserts nothing was displaced. It is
called from the shared driver helper, so **six of `serial.rs`'s twenty-eight tests reach it**,
the two hundred random layouts among them, and from the 305–330 deletion at a 20-base building
width, which is the plan's named shape.

**What it is discriminating against, measured.** It re-runs the driver's loop rather than
calling the driver, because `merge_cohort_through_cache` returns one merged outcome where the
organiser needs them region by region — so a defect in the *driver's* loop is caught by the
byte-identity comparison and not here. What it pins is the **eviction discipline**: move its
eviction from the building region's first base to its last, one line, and the suite reports *at
building regions of 29 bases a builder produced a locus on ground an earlier locus already
owned*. That matters because the review found the safety-net argument reads as though it rests
on the cache, when it rests on where the caller chooses to evict — which is exactly what E3
rewrites.

## 2. The rest, all applied

- **`resolve_and_release` no longer peeks-then-takes.** `Peekable::next_if` consumes the failed
  spans in a plain `for`, so neither `expect` remains; `claim_and_release` is split out. The
  locus-first tie-break is unchanged.
- **An inverted `GenomeRegion` walked the frontier backwards** and silently disarmed the rule
  with the counter reading zero. `first_base_of` and `last_base_of` now order their ends, as
  `with_observations`, `cover` and `building_regions_of` already do and for the same stated
  reason. Both ends are pinned — a locus written back-to-front owning its ground, and one
  arriving on owned ground losing it.
- **A displaced locus owns nothing**, which no fixture separated: every one had the displaced
  locus ending inside the standing owner's span, so a `claim` that extended the frontier on
  displacement survived the whole suite.
- **The guard's boundary and its sort key** are pinned: two spans sharing exactly one base are
  refused, two adjacent ones are not, and three spans interleaved across two contigs are
  refused — the case a contig-less sort key hides.
- **`finish` returns a `MergeTally`** of the failed and displaced counts, so the natural call
  order cannot lose the one whose whole purpose is to be noticed.
- **"Earlier" means claimed first**, not lower-numbered, where region order and start order
  disagree. Stated on `claim` and pinned.
- **The eviction bullet names the discipline and its owner** rather than the cache.
- **`Default` forwards to `new()`** instead of deriving a second construction path past the
  longhand constructor; `claim` gained `#[must_use]`.
- **Six documentation claims rewritten, four of which were wrong** — see the review's §5. The
  worst said `serial`'s *drivers* assert over *every fixture* that loci are *disjoint and
  ascending*; it is one driver, six of twenty-eight tests, and disjointness only.

## 3. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `209 passed; 0 failed` (189 before E2).
- `cargo test --lib` — `3832 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out;
  finished in 599.40s` (3,812 before E2).
- **`tmp/mutate_e2.sh` — 21 mutations, 21 killed.** Nine on the rule itself, three on the
  frontier's handling of inverted spans and the tally, three on the disjointness guard's
  condition, sort key and boundary, two on the eviction point (the driver's and the end-to-end
  check's), and one each on the displaced locus's ownership, the resolution order and the
  meaning of "earlier".

Two mutations the reviewer proved to be **no-ops rather than survivors**, recorded so nobody
writes a test against a hazard that does not exist: the merge tie-break `<=` → `<` (a probe
counted zero ties over the whole suite), and a `max` on the frontier (zero retreats).

## 4. Follow-ups

- **`displaced_locus_count` still has no consumer.** Where it surfaces belongs with the
  failed-locus count, in the run summary the emission step owns (spec §13); `MergeTally` is the
  shape that carries both out of the organiser.
- **The eviction discipline is E3's to keep.** With several builders in flight the safe point
  is the earliest live region's first base, not the latest, and `refuse_displaced_loci` is the
  regression test for getting it wrong.
