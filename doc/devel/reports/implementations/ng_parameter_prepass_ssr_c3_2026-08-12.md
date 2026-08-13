# ng step 4, the STR path — C3: the read cap, drawn from the locus's position

*Implementation report, 2026-08-12. Step C3 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md) — the plan's second
**own commit, do not bundle** step, because its failure is silent. With the review that followed
and the fixes applied: two agents, 20 mutations, 6 behaviour-changing survivors. Design authority:
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.1, §2.3 and
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1.*

## What the step is

A locus deeper than twelve reads is entered from a **uniform random subsample** of its reads down
to that cap, never by dropping the locus: dropping deep loci is depth-dependent selection, the
bias this step exists to remove, while a uniform subsample leaves the bucket counts distributed
exactly as at the lower depth, so the cap costs precision and never bias.

**The draw is seeded from the locus's position and from nothing else**, which is the plan's named
silent failure. Seeded so, a region-sharded walk and a single-threaded one keep the same reads at
every locus and their tables merge as an equality; seeded from a counter, a thread or the clock,
they diverge by a few reads at each deep locus — a difference nothing short of comparing two whole
walks would show, and one that would make a fitted level move with the thread count.

## Recorded deviations from the architecture

1. **The sampler is lifted, not copied.** `SelectionWalk`, `seed_at` and `splitmix64` move out of
   the SNP/indel path's `depth_and_alt_reads.rs` into a shared
   `parameter_estimation/subsample.rs`, unchanged. Both caps are the same draw and it is a subtle
   one; the plan named no reuse target, and copying it would have left two copies of an algorithm
   whose one contract is that it produces the *same* stream. The generic path keeps its own
   two-category `hypergeometric_draw` on top of it. The reviewer diffed the moved items against
   the previous commit: the only code differences are visibility, one parameter rename, and
   `as u64` → `u64::from` on the same type — the draw is bit-for-bit the stream it was.
2. **`shape_of` takes the tally and the position**, where arch §2.3 sketches one function from the
   locus. The caller needs what the tally counts *besides* the shape — the reads left out for
   having witnessed only part of the tract — and re-deriving the tally inside would walk each
   locus's observations twice on a walk of millions.

## What the review changed

**Blocker — the test guarding this step's central property could not fail.** The fixture helper
meant to move a locus assigned the new start before deriving the end from it, so the end stayed
put and the *span* varied instead: `(start + old_end) − start` is `old_end`. Every test that
claimed to move a locus was varying two things at once, and a draw seeded from the region's span
rather than its position passed all 32 tests — in production that is a seed depending on tract
length alone, so every same-length tract in the genome would share one stream. With the helper
corrected the mutant draws one distinct shape over 200 positions where the real code draws ten. I
reproduced it: it now fails four tests.

**Blocker — the variance check's tolerance was wider than the whole signal.** The test's doc said
a draw made with replacement "would match the mean and miss the variance". At depth 300 with a cap
of 12 the two variances are 2.5686 (without replacement) and 2.6667 (with) — a difference of
0.0981, against a tolerance of 0.1. The test could not have rejected the wrong model at any sample
size, and with-replacement sampling is the natural wrong implementation of "take twelve of these
reads": it double-counts, so the cap stops being unbiased. Now 100,000 draws and a tolerance of
0.05, which leaves the true model ten times the margin. I wrote the with-replacement draw and it
now fails at 2.6813 against 2.5686.

**Major — the draw is documented as a format and nothing pinned it.** `seed_at`'s own note says it
may not change silently with a dependency bump, because it decides which reads a fit sees; every
other test asserts a *relation*, and reordering the walk or moving where its state advances leaves
them all green while changing every capped locus's draw — including the SNP/indel path's, since
that code is now shared. Four positions are pinned to the shapes the algorithm actually draws.

**Major — no fixture combined a capped locus with partial witnesses**, so widening the draw's
population to include them entered a locus at five reads instead of twelve, silently, since a
short shape is a legal shape. **Major — the shallow path's guard passthrough was untested**, so
dropping it entered ten reads as six and biased the guard share towards zero on the majority of
loci. **Minor — the cap boundary was tested from one side**: depth 13, where the walk first runs,
and depth 1 are now both fixtures.

**On the property the plan states, an honest limit.** The plan asks that the same locus give the
same draw "in every shard layout". That half is **not testable at this step and I did not fake
it**: `shape_of` takes one locus and never sees a shard, so a test looping over the same pure
function in two groupings would pass by construction. What *is* testable is covered — the draw is
reproducible, moves with the position, and carries no state between calls (a reviewer's
thread-local stream standing in for a per-worker RNG is killed by the reproducibility test).
The reviewer's file records what C5's test must then assert, and it is worth carrying forward:
an equality rather than a closeness, a cut that splits a stratum across shards, a fixture
dominated by loci deeper than the cap (below it nothing draws, so any implementation passes), at
least two thread counts, and at least two contigs.

**Two wrong numbers of mine, and they were the same error twice**: the draft commit message said
six tests and a suite of 3,476. Nine test functions land — six here plus two new ones with the
lifted sampler and one that moved with it — and after the review's additions the count is twelve
new tests and a suite of 3,482.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation` | **439 passed**, 5 ignored (both paths, so the lift is exercised) |
| `cargo test --lib --bins --tests --all-features` | **3,482 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` gives **33** in `locus_offsets.rs` (23 before
this step) and **3** in `subsample.rs` (one of them moved), so the step adds twelve; the suite
moved 3,470 → 3,482.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets` fails
in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-c3/` — two per-category files (reliability; module
structure, errors and the numbers check) and the reviewed patch.
