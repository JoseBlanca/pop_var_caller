# ng step 4, generic path — Milestone B review (B1–B4)

**Date:** 2026-08-06. **Scope:** `src/ng/parameter_estimation/generic/histogram.rs` across
Milestone B, reviewed in three rounds against the commits that produced it.
**Verdict per round:** Request-changes → resolved, all three.

Per-category agent files are in the gitignored `tmp/review_2026-08-06_b1b2/`,
`tmp/review_2026-08-06_b3/` and `tmp/review_2026-08-06_b4/`; this is the synthesis.

## 1. How it was run

| round | commits reviewed | agents | categories | outcome |
|---|---|---|---|---|
| 1 | `b3bc6541` (B1), `2768ded2` (B2) | 3 | reliability, refactor_safety, extras, errors, naming, defaults, tooling, idiomatic, smells, module_structure, unsafe_concurrency | 1 Blocker, 8 Major, 14 Minor, 6 Nit → `367a1450` |
| 2 | `375159b4` (B3) | 2 | the same, minus tooling | 2 Blocker, 4 Major, ~10 Minor → `498a9f53` |
| 3 | `b01b0212` (B4) | 2 | the same | 2 Blocker, 6 Major, ~9 Minor → `698302f1` |

Every agent ran in its own git worktree, detached at the commit under review, with a
two-file existence check and a hard stop if either was missing. **The write to the main
checkout's `tmp/` was refused by worktree isolation for all but one agent**; each wrote
inside its own worktree and named the absolute path, and the files were copied across
before the worktrees were pruned. That is now a known cost of the isolation and the
prompts should keep instructing for it.

## 2. What the rounds actually found

**Every substantive finding came from mutation, not from reading.** Across ~60 mutations
the suite killed most and let a handful through, and the ones it let through are the
findings worth carrying:

### The five Blockers, and four of them are missing tests

1. **`add_attributed_site`'s covered-position accumulation had no test** (round 1).
   Deleting the line left all 2,963 tests green. Every other attributed test asserts on
   `cells()` or `total_loci()`, neither of which reads the covered positions, and the one
   test that does read them used `add_site` twice and the attributed arm never. That arm
   is the entry point for every window of a multi-library sample.
2. **Two same-typed counters, transposable, surviving everything** (round 2). `cells()`
   unpacked the attributed cell's `(count, depth_sum)` positionally. Transposed, a cell
   of three sites at depths 18, 19 and 20 came back holding **57 sites at a mean depth of
   0.05** — one alternative read scored at a twentieth of a read — with all 2,983 tests
   green. Both round-2 agents found it independently. It is the same shape the A5 review
   found in two `SmallVec<[f64; 3]>` fields.
3. **The B3 oracle touched one bin of twenty** (round 2). Recording depth zero for every
   site below depth 18 passed the whole suite. Depths 0–17 are the entire exact-per-depth
   region, where 97 sites in 100 of a three-read tomato cohort sit.
4. **B4's overflow guards, two of three untested** (round 3). Silently dropping an
   attributed cell's overflow, and silently dropping the merge's covered-position
   overflow, both left the module green. Each drops evidence and returns normally, which
   is the module's defining hazard exactly.
5. **The one test that looked like it covered the third could not fail on that half**
   (round 3): it overflows *both* halves at 3 × 10⁹, so the site count's `?`
   short-circuits and the depth-sum guard is never the reason for the `None`. The
   sites-fit-depth-sum-does-not case is the only one that actually happens.

Of B4's 17 mutations, 14 were caught and **all three survivors were the `checked_add`
refusals B4 itself added** — every mutation that changed a value on a reachable path was
caught, including the transposition class that survived everything in B3. That is the
split-at-every-boundary oracle working; what it structurally cannot reach is a guard whose
failure is a silent no-op.

### The findings that were wrong code

- **A guard that fired in one of three regimes** (round 1). `SiteKey::attributing`'s
  documented duplicate-read-group panic sat after the early return that pools, so a
  listing whose duplicated entries summed above the bound pooled quietly —
  `attributing(bin, &[(g4,3),(g4,3)])` returning `alt_reads = 6` with no panic. It
  escapes exactly where the damage is worst, and `add_attributed_site`'s cross-check
  cannot catch it because a caller that double-counted builds both arguments from the
  same reads.
- **Neither entry point checked a site's depth against the cap** (round 1). A 500-read
  site entered the 98–124 bin silently, and B3's per-cell depth sum is what would then
  carry a depth no site in that bin could have had.
- **A supplied `widen` function could narrow** (round 3). `absorb` took
  `impl Fn(N) -> C`; `|counter| counter as u32` compiled and four billion sites came back
  as three at a mean depth of 1.0, the truncation happening before `checked_add` saw the
  value. `N: CellCounter + Into<C>` makes it `error[E0277]`.
- **The read-group table has no fold to widen it** (round 3). `whole_sample_histogram`
  folds *windows*; the read-group histogram is genome-wide and is not keyed by them, so
  at `u32` its busiest cell's depth sum passes 4.29 × 10⁹ about a third of the way
  through a human sample — this module's own widening argument, applied to the table the
  fold cannot reach. The width choice is C3's; the docs and both overflow messages now
  say which width that table needs.
- **The fold dropped the architecture's "for one ploidy"** (round 3), from its signature
  and its prose alike. A table carries no ploidy — `cells()` stamps one on read — so a
  mixed-ploidy fold produces cells all scored against one genotype set, and haploid
  sites, which can never be heterozygous, enter the heterozygosity fit as diploid ones.
- **Nothing ran a merge across threads** (round 3), which is the whole purpose of B4:
  one accumulator per region shard, filled in parallel, merged at the end. `Send + Sync`
  held for both widths and was asserted nowhere.
- **A test that did not test its name** (round 3): the empty-fold test asserted 583
  cells, which a fold building its result from a *fresh* ladder also satisfies — a table
  of the right shape that no later merge could accept. It now asserts pointer identity
  against the ladder it was handed.

### Documentation accuracy, again the highest-yield category

Round 1 found "four is the measured default" (research §2.5 says the bound "is not
currently buying precision"), "0.3% in each genotype frequency" (0.23% and 0.30%), a
`u32`/`u64` justification quoting the number that shows the widening is *unnecessary*,
and two doc claims the file's own test disproves. Round 2 found **a wrong number in the
author's own commit message** — "2,875 cells" where the fixture holds 125 cells and 2,825
sites, overstating the oracle's reach 23-fold, with the test asserting only
`cells.len() > 100`. Round 3 found none: every one of its nine named figures checked out.

That is three rounds of five on this plan with a wrong number in prose, and the pattern is
consistent — **the figures that go unchecked are the ones describing the author's own
test's reach**, not the ones taken from the research note.

## 3. What was declined, and why

- **Renaming `SiteKey::attributing` and `add_attributed_site`** (round 1, single-source).
  The finding's real content — that a reader picks between the two per *site* when the
  question is per *sample* — is what the Major-2 doc fix closed, and `add_site` is the
  architecture's name.
- **Renaming `depth_and_alt_reads.rs`** (round 1). It is named in the plan's scope list
  and in A1; renaming a file the plan names is an owner call, recorded in the
  implementation report.
- **`DepthAndAltReads` as a shape name** (round 1, filed as a Nit by its own reporter).
  The architecture's name.
- **Renaming `whole_sample_histogram` to `fold_windows_of_one_ploidy`** (round 3). The
  reviewer's point stands — the name reads as "all of it" and invites the mixed-ploidy
  fold — but it is the name the plan's B4 and the architecture both use, and C3's method
  of the same name will pass the selected windows. The obligation went into the doc
  instead.

One finding was accepted in a stronger form than proposed: the reviewer asked for a doc
warning about `Arc::new((*edges).clone())`, then observed that `DepthBinEdges: Clone` has
no user and cannot have one — `new()` takes no arguments, so every ladder a run can build
is byte-identical. Dropping the derive turns a documented footgun into a compile error.

## 4. What the reviews are owed next

- **C3 must declare the read-group histogram as `DepthAltHistogram<u64>`.** This is the
  one finding that outlives Milestone B.
- **`merge` is not atomic**, which round 3 verified: a table whose merge panicked has
  already absorbed the cells the fold reached. Nothing recovers from those panics today;
  it is documented so that a later `catch_unwind` is not written on the assumption that
  it could.
- **Two Milestone-A owner items are still open** — the arch module table naming four
  files under `generic/` where five exist, and the "29% covered by runs" figure in
  `spec` §6.5 and `arch` §5.3 against the research note's realised 0.2629.
