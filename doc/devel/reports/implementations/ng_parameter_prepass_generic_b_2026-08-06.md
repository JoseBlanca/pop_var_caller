# ng step 4, generic path — Milestone B: the cell table

**Date:** 2026-08-06. **Branch:** `ng-parameter-estimation`, from `5d13beff`.
**Plan:** `doc/devel/ng/impl_plan/parameter_prepass_generic.md`, Milestone B (B1–B4),
stopping at Checkpoint B. **Design:** `doc/devel/ng/arch/parameter_prepass_generic.md`
§2.2 and §3, `doc/devel/ng/spec/parameter_prepass_generic.md` §4 and §9,
`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §4.

## 1. What was built

One file, `src/ng/parameter_estimation/generic/histogram.rs`, holding the tally every fit
in step 4 reads. Each covered position reduces to a key — a depth bin, an exact
alternative count, and for a multi-library sample which library each of the first few
alternative reads came from — and the table counts how many positions showed each key,
alongside the sum of those positions' exact depths. Eight hundred million positions
become a few hundred counters.

| commit | step | what landed |
|---|---|---|
| `b3bc6541` | B1 | `DepthAndAltReads`, `SiteKey`, `CellCounter`, `MAX_ATTRIBUTED_ALT_READS` |
| `2768ded2` | B2 | `DepthAltHistogram<C>`: storage, `add_site`, `add_attributed_site`, `cells`, `total_loci`, `total_covered_positions` |
| `367a1450` | — | the B1+B2 review applied |
| `375159b4` | B3 | the per-cell depth sums and `mean_depth_in_cell` — **own commit, do not bundle** |
| `498a9f53` | — | the B3 review applied |
| `b01b0212` | B4 | `merge`, `absorb`, `whole_sample_histogram` |
| `698302f1` | — | the B4 review applied |

`src/ng/parameter_estimation/generic/depth_bins.rs` gained `DepthBinEdges::bins()` and a
doc warning about cloning a ladder rather than its handle. Nothing else was touched.

## 2. Deviations from the architecture, and why

The architecture's own rule is that its signatures are illustrative and its measurements
are not. Five signatures were deviated from; every measurement was reproduced.

- **`SiteKey` is a struct with private fields, not the public two-variant enum** of arch
  §2.2. A variant's fields are exactly as public as the enum, and two of this key's
  properties are invariants whose violation is silent: the read-group listing must be
  canonical (sorted, no zero entries), and which arm a site takes must follow from its
  counts rather than from the call site. Neither is expressible on an enum.
  `attribution()` returns a borrowed two-variant *view*, `Attribution`, where the
  objection does not transfer — a view is returned, never constructed by a caller.
- **`attribution()` returns `Attribution`, not `Option<&[…]>`.** An `Option` invites
  `unwrap_or(&[])`, which turns "several libraries contributed and this key forgot which"
  into "no library showed an alternative read" — different arithmetic, silently.
- **`add_site` takes the locus's reference span as a second `Bp` argument.** Arch §2.2
  calls `total_covered_positions` "accumulated alongside" and gives `add_site` one
  argument; a locus's span is not derivable from its reads.
- **`covered_positions` is `u64`, not the generic counter width.** The memory argument
  that forces a width choice on the cell vectors — eight bytes a cell, ~4.7 kB a window,
  ~37 MB a tomato sample — does not reach one scalar per table, and a read-group table
  accumulates genome-wide with no fold to widen it, putting a human genome's 3.1 × 10⁹
  analysable positions at 72% of the `u32` ceiling.
- **`cells()` returns `Vec<Cell>`, not `Vec<(SiteKey, Ploidy, u64)>`.** This vector is the
  seam between `generic/` and `fitting/`, walked 161 times per fit by code that would
  otherwise be written in `.2`. Introduced now because retro-fitting it once the fits
  exist would touch every fit.

Two internal types have no counterpart in the architecture: `CellTally { sites, depth_sum }`
(see §4) and the private `absorb`. `whole_sample_histogram` is a free function taking the
already-selected windows, where arch §3 has it as a method on `GenericAccumulators`
taking a `Ploidy`; C3 writes that method and it calls this. The restriction the parameter
carried is stated on the free function as an obligation, because nothing in it can check
one. And `DepthBinEdges` lost its `Clone` derive, so that cloning a ladder rather than
its handle is a compile error rather than a merge that panics hours later.

## 3. Deviation from the plan: where `depth_sums` lands

The plan puts the per-cell depth sums in **B2** and `mean_depth_in_cell` in **B3**. They
landed together in B3.

A private field written and never read fails `-D warnings`, so a B2 that accumulated
depths with no reader would not have been green. The repo has hit this once before — the
alignment plan's A2/B1/B2 merged into one commit for the same reason — but that
resolution is unavailable here, because B3 is one of the six steps the plan marks *own
commit, do not bundle*, precisely so a `git bisect` over a moved depth can land on it.
Splitting the mechanism across two commits would have put the storage in one and the only
thing that reads it in another; keeping it whole is what the isolation is for.

## 4. What the reviews found

Three loops, eleven agents, ~60 mutations. **Two Blockers, both of them missing tests
rather than wrong code**, and one real bug.

- **B1+B2 — a guard that fired in one of three regimes.** `SiteKey::attributing`
  documented an unconditional panic on a duplicate read group, but the check sat after
  the early return that pools, so `attributing(bin, &[(g4,3),(g4,3)])` returned a pooled
  key with `alt_reads = 6` and no panic — a double-counted read set reaching a cell as a
  heterozygote invented out of nothing. `add_attributed_site`'s own cross-check cannot
  catch it, because a caller that double-counted builds both arguments from the same
  reads and the two agree.
- **B1+B2 — the Blocker was a test that was not there.** Deleting `add_covered_positions`
  from `add_attributed_site` left the whole suite green. That arm is the entry point for
  every window of a multi-library sample; lost, each would report zero covered positions
  and the inbreeding coefficient would come out of a weighting nothing prints.
- **B1+B2 — neither entry point checked a site's depth against the cap.** `bin_for` is
  total by design, so a 500-read site entered the 98–124 bin silently. B3's depth sum is
  exactly the quantity that would then hold a depth no site in that bin could have had.
- **B3 — two same-typed counters, transposable, surviving everything.** `cells()`
  unpacked the attributed cell's `(count, depth_sum)` positionally; transposed, a cell of
  three sites at depths 18, 19 and 20 came back holding **57 sites at a mean depth of
  0.05**, and all 2,983 tests passed. The cause: the oracle table is built entirely
  through `add_site`, so it contains zero attributed cells, and nothing asserted
  `Cell::sites` or `Cell::mean_depth` on that arm at all.
- **B3 — the oracle touched one bin of twenty.** Recording depth zero for every site
  below depth 18 passed the whole suite, and depths 0–17 are the entire exact-per-depth
  region where 97 sites in 100 of a three-read tomato cohort sit.
- **B4 — the `widen` closure could narrow.** `absorb` took `impl Fn(N) -> C`, which the
  type system cannot hold to being a widening: `|counter| counter as u32` compiled, and
  four billion sites came back as three at a mean depth of 1.0, the truncation happening
  before `checked_add` ever saw the value. Replaced by `N: CellCounter + Into<C>`, which
  makes the narrowing `error[E0277]`.
- **B4 — the read-group table has no fold to widen it.** `whole_sample_histogram` folds
  *windows*; the read-group histogram is genome-wide and is not keyed by them, so at
  `u32` its busiest cell's depth sum passes 4.29 × 10⁹ about a third of the way through a
  human sample. The module had already made this argument once, for `covered_positions`,
  and the cell counters did not inherit it. Recorded in the type's doc and in both
  overflow messages; **the width choice itself belongs to C3** (§6).
- **B4 — the fold dropped the architecture's "for one ploidy" restriction** from both its
  signature and its prose. A table carries no ploidy — `cells()` stamps one on read — so
  folding a haploid window into a diploid sample's table lets haploid sites, which can
  never be heterozygous, enter the heterozygosity fit as diploid ones. Restored as a
  stated obligation on the caller, since nothing in the fold can check it.
- **B4 — all three of its overflow guards but one were untested**, and each drops
  evidence and returns normally when broken. Everything else B4 added was caught by the
  split-at-every-boundary oracle, including the transposition class that survived
  everything in B3.
- **B4 — nothing ran a merge across threads**, which is the step's entire purpose.
  `Send + Sync` held for both widths but was asserted nowhere.

Three responses were structural rather than local. `CellTally { sites, depth_sum }`
replaced both the anonymous pair and the two parallel vectors, so a transposition stops
compiling on either arm and the dense arm's two halves can no longer be written at
different indices. `absorb` opens with an exhaustive destructure, so a field added later
stops it compiling rather than being silently left out of every merge. And
`CellCounter::checked_add` is checked rather than wrapping, because `[profile.release]`
leaves `overflow-checks` off.

**Three wrong numbers in prose**, one of them in this work's own commit message: "2,875
cells" where the B3 fixture holds **125 cells and 2,825 sites** — the one figure stating
how wide the oracle is, overstating it 23-fold, with the test asserting only
`cells.len() > 100`. Also "four is the measured default" for `MAX_ATTRIBUTED_ALT_READS`
(research §2.5 says the bound "is not currently buying precision" and two is equally
good) and "0.3% in each genotype frequency" (the adopted row is 0.23% and 0.30%). B4's
prose was the first round of this plan with no wrong number in it.

## 5. The oracles, and what each is worth

| step | oracle | what it caught, or could not |
|---|---|---|
| B2 | one site into all 583 cells; 583 cells of one site each | the **only** test that catches a mis-sized or mis-placed row — both width mutations died by it alone |
| B3 | `alt_reads ≤ mean_depth_in_cell` at every cell of a table of depths 100–124 | kills the bin mean. **One-sided**: passing the alternative count to the depth sum instead of the depth leaves it green, because the wrong value saturates the inequality |
| B3 | the same table's bin mean shown violating it — 112.46, twelve counts 113–124 | kills exactly the mutation the identity cannot see. The plan's insistence on a *pair* was load-bearing |
| B3 | all 7,875 legal `(depth, alt)` pairs over all 583 cells, against a tally the test keeps | added after review; the identity alone reached one bin of twenty |
| B4 | the fixture split at every one of its nine boundaries, merged both ways round | 18 comparisons where the plan asked for one. Caught 14 of 17 mutations; the three it could not reach were the overflow guards, now tested |
| B4 | a `u64` fold equals a `u32` walk, cell for cell | genuine rather than trivial because `Cell` is width-independent |
| B4 | shards filled on separate threads, folded, against the single walk | added after review; nothing had run the merge across threads |

## 6. Open items for the owner

- **The read-group histogram's counter width is C3's to set, and it must be `u64`.**
  Arch §3 sketches `BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram>` — the `u32`
  default — which overflows its depth sums on a human sample with no fold to widen it.
  The type doc and both overflow messages now say so; the declaration is C3's.
- **Two owner items carried from Milestone A remain open**: the arch module table and the
  plan's A1 still name four files under `generic/` where five exist, and
  `spec` §6.5 / `arch` §5.3 still say "29% covered by runs" where the research note's
  realised `F` is 0.2629.
- **`cargo test --all-targets` is red on this branch and it is not step 4's.**
  `benches/psp_writer_perf.rs:386` panics with `index out of bounds: the len is 3300000
  but the index is 3300000`, in its own `flush_block_one` priming loop, which walks
  `phase_records` until a projected block size crosses a threshold and runs off the end
  when it never does. That bench and all of `src/psp/` are byte-identical to `5d13beff`
  and reference nothing under `parameter_estimation`. Everything here was validated with
  `--lib --bins --tests`, which is the whole test suite — benches carry no tests.
- **Declined review findings, recorded rather than silently dropped**: renaming
  `SiteKey::attributing` / `add_attributed_site` (single-source, and the confusion it
  named — per-site versus per-sample — is what the doc fix closed); renaming
  `depth_and_alt_reads.rs`, which the plan names; `DepthAndAltReads` as a shape name,
  which is the architecture's.

## 7. Validation

All via `./scripts/dev.sh`, at every commit:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --all-features` — **2,942 → 2,996 passed, 0 failed,
  5 ignored**. `ng::parameter_estimation` holds **94** tests, 40 → 94 over the milestone.
- `cargo doc --no-deps --lib` — **12 unresolved links, the pre-existing baseline**,
  unchanged. `Cargo.toml` sets `broken_intra_doc_links = "deny"` and the other three
  gates do not cover rustdoc; design-doc references are written as plain code spans
  throughout, never as intra-doc links.
