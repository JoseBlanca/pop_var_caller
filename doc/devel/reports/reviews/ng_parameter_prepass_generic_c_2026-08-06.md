# ng step 4, generic path — Milestone C review (C1–C3)

**Date:** 2026-08-06. **Scope:** `generic/depth_and_alt_reads.rs` and
`generic/accumulators.rs` at `9d3e090d`, covering commits `1977c465` (C1), `4b01681d`
(C2), `eccccf7c` (the C2 amendment) and `9d3e090d` (C3).
**Verdict:** Request-changes → resolved in `1c7b4490`.

Three agents in isolated worktrees, eleven categories, about sixty mutations. Per-category
files are in the gitignored `tmp/review_2026-08-06_c/`.

| agent | categories | outcome |
|---|---|---|
| 1 | reliability, refactor_safety, extras | **2 Blockers**, 6 Majors, ~4 Minors; 20 mutations survived |
| 2 | errors, naming, defaults, tooling | 2 Majors, ~7 Minors |
| 3 | idiomatic, smells, module_structure, unsafe_concurrency | **1 Blocker**, 3 Majors, ~7 Minors |

## 1. Three Blockers, and all three were tests that could not fail

That is five review rounds out of seven on this plan where the Blocker was a missing test
rather than wrong code. The code keeps being right; what keeps being absent is anything
that would notice if it were not.

1. **`merge` could drop the read-group table and all four counter sums with the suite
   green** — found independently by agents 1 and 3, which is convergent evidence. Five
   separate deletions survived. With the by-group map skipped, a region-sharded walk
   produces an **empty read-group table**, which is every error rate this step exists to
   fit, and nothing panics. With the `loci_overlapping_previous` sum dropped, the counter
   the design says *must be zero* reads zero on a sharded run whatever the shards saw —
   the single guard that the windowed table has not entered a site twice.

   The cause is a fixture: `three_shards_merged_in_every_order_give_the_same_counters`
   compares five zeros against five zeros, because its twelve loci carry no upstream cap,
   no unwitnessed reads, no depth above the cap and no overlap. The two tests that *do*
   exercise counters never merge.

2. **`SelectionWalk`'s resumability across the all-kept fast path was never exercised.**
   Deleting the population decrement in that arm leaves every later group drawing against
   a population one too high — and nothing panics, because the entry's attribution still
   sums to its own total, so `add_attributed_site`'s cross-check passes. The arm fires
   exactly where the design says the damage lands: the deep, allele-rich sites where the
   cap bites. It is the neighbour of the bug C2 caught pre-commit, and unlike that one it
   was invisible.

3. **`merge` checked three of `new`'s four arguments.** The unguarded one was the ploidy
   map: two shards handed different maps merged silently and every cell one of them filled
   was scored against the wrong set of genotype classes — a haploid region has two and a
   diploid three. Adding `Arc::ptr_eq` immediately failed two of the author's own tests,
   because the fixture minted a fresh map per shard, a configuration no real run produces.

## 2. Findings that were wrong prose, and one that was upstream

- **arch §2.2's worked example is on the wrong side of its own argument.**
  `round(1 × 124/500) = 0`, so "a 500-read site with one alternative read becomes
  `(124, 1)`" cannot happen; the sign reversal it illustrates is at depth 248. The design
  document, the code and C2's commit message all carried it. Corrected in the first two.
- **"two groups, which is every sample that has more than one library at all" is false**,
  and it cost allocations rather than only accuracy. `spec/read_groups.md` records 133
  samples with two libraries, 20 with three and four with 7, 16, 16 and 42. At two inline,
  `add_locus` allocated twice per locus from three groups upward — 200 per 100
  steady-state loci, measured, about 47 ns a locus — breaking the architecture's
  no-allocation contract on exactly the multi-library samples the machinery exists for.
- **"binning a site costs no atomic at all"** was contradicted by an `Arc::clone` per
  locus; the reviewer showed `let edges = &self.edges;` compiles, so the claim can be made
  true rather than softened. Deferred with the split-borrow rewrite.
- **"No world reaches a cap — the deepest site is 125, against a cap of 124"** is
  self-refuting as written; the research note hedges it across two caps.
- Smaller: `f·(1−f)` described `f` as "the fraction dropped" when it is the fractional
  part of the rescaled count; "one site in a thousand" for a homozygous non-reference site
  named no organism, where tomato's measurement is 6 per kilobase.

**Nine other quantitative claims were checked and hold**, including the hypergeometric
mean and variance (62.248 and 23.358), the stochastic round's 0.186, `1/r = 4.03`, the
modulo bias at 4.34 × 10⁻¹⁶, the 250/8,000 pileup caps, 1,550 of 1,707, and the 1-based
window arithmetic.

## 3. What the reviewers confirmed rather than found

Worth recording, because each was an argument the author had made without evidence:

- **`GenericAccumulators` is `Send + Sync`**, and works under a real rayon
  `par_iter().fold(…).reduce(…)` over 4,000 loci, matching the sequential walk exactly.
  `merge(&mut self, other: Self)` slots into `reduce` through a one-line wrapper.
- **`previous_end` not being merged is sound** for that shape.
- **`add_locus` allocates nothing per locus at one and two read groups**, and the
  `mem::take` pair has no reachable early return between take and put-back. A panic is
  reachable and costs only the allocation, leaving a valid empty `Vec`.
- **The `covered_spans` linear `find`, the `Arc<dyn PloidyMap>` virtual call and the
  per-locus `Arc::clone` are each 1–2 ns** and not worth acting on. Only the `SmallVec`
  was.

## 4. Deferred, with reasons

- **`AccumulationCounts` is one type with two meanings** — the stored value's
  `shard_spans_overlapping` is always zero, only `adjustments()`'s is real. A
  stored/reported split would say it in the type. `merge` now destructures exhaustively,
  which closes the concrete hazard.
- **`CountedSite` means three slightly different things** across its three producers; on
  the by-group path it is a group's slice of a site. No consumer sums it today.
- **`generic/mod.rs` now does four jobs** and splits `WindowIndex` from `WindowKey` across
  files. Agent 3's answer to "has the whole become incoherent" was: yes, in `mod.rs` only.
- **The split-borrow rewrite of `add_locus`** — built by a reviewer, passes, removes both
  `mem::take`s and the per-locus `Arc::clone`. Worth taking when that function is next
  touched.
- **`shard_spans_overlapping` is a lower bound**, not a count: three shards over one
  stretch report 1. Documented rather than changed, because zero-detection is exact and
  that is the property that matters.
