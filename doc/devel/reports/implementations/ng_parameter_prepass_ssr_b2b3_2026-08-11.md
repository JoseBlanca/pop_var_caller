# ng step 4, the STR path — B2+B3: the table of locus shapes, its merge, and its two diagnostics

*Implementation report, 2026-08-11. Steps B2 and B3 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), taken in one loop, with
the review that followed and the fixes applied — four agents, 35 mutations, 5 survivors. Design
authority: [`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.2 and
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.1, §5.*

**Two plan steps in one loop, named rather than silent.** B2 is the table's storage and B3 its
merge and two diagnostics — the same type's surface, and B2 alone would be a table that cannot be
merged or read. The plan marks neither *own commit, do not bundle*.

## What the step is

`StratumTable`: how many loci showed each shape, in a `BTreeMap`, plus the two base counts the
substitution rate is divided from. Sparse because the space of possible shapes at one depth is 220
at three reads a locus and 293,930 at twelve (arch §2.2), of which only a data-dependent corner is
occupied — 1.73 million tomato loci made 70,305 entries, one per 25, and HG002 at 300× made 12,727
for 29,811 loci, one per 2.3 (spec §4.1). A `BTreeMap` rather than a `HashMap` because a fit is a
floating-point sum over the entries and floating-point addition is not associative.

`merge` is element-wise integer addition, so a table split any way and merged in either order is
the table the whole walk would have built — an equality, not a tolerance. `substitution_rate` is
mismatched over compared, and it is the maximum rather than a moment estimate. `not_whole_repeat_share`
is the guard share `GUARD_SHARE_LIMIT` bounds.

## Recorded deviations from the architecture

1. **`BaseComparison`, one type where arch §2.2 sketches two `u32` arguments.** The two counts have
   the same type and one is a subset of the other, so a swapped pair is not obviously wrong at the
   call site — and Milestone C's `composition_of` returns exactly this pair. Refusing
   `mismatched > compared` at construction is also what makes the fitted rate a probability by
   construction rather than by hope.
2. **`entries() -> Vec<StratumEntry>` where the plan and arch say `shapes() -> Vec<(LocusShape, u64)>`.**
   Two changes in one: the sibling path overrode the same sketched tuple for the same reason —
   nothing in a tuple says the second member counts *loci* — and the word the design uses for a
   shape plus its count is *entry* (arch §2.2 is titled "The entry, and the table of entries"), so
   a method returning entries should say so. Nothing outside this file consumes it yet.
   `entry_count()` lands beside it, so asking how many entries a stratum holds does not build the
   list.
3. **`PooledBases`, a private pair, rather than two loose `u64` fields.** The same argument as
   `BaseComparison`, applied to the accumulated counts: they were being passed positionally through
   two helpers, which is where a transposition would live.

## What the review changed

**Major — the acceptance oracle for B3 could not fail.** The plan asks that the substitution rate
be "proven to be the maximum, not merely a ratio". My test pinned the table's answer to the literal
`0.0030` within `1e-12` and *then* checked a grid search agreeing within one grid step of `1e-5` —
ten million times looser, so the grid assertions could never be the ones that fail. The reviewer
measured it: perturbing the rate by 5e-6 and by 2e-5 failed the literal both times, never the grid.
The grid was also computed from the test's own constants rather than from anything the table
returned, at one interior point. It now compares the grid's argmax against `table.substitution_rate()`
at four (compared, mismatched) pairs spanning 1 in 400 to 63 in 64, and the two boundary rates the
grid cannot reach — every base matching, every base mismatching — are asserted separately, because
the binomial score is not a number at either endpoint.

**Major — none of the three overflow guards had a test, and all three are reachable.** Their doc
comments called them unreachable ("4.3 billion loci of one shape, which no genome holds"), but
`merge` takes `&Self`, so merging a table with a clone of itself doubles every count and thirty-two
doublings pass `u32::MAX`. Turning each guard into the wrapping add it exists to prevent left all 61
tests green, and two of the three then reported a *wrong number* rather than failing: the stratum's
commonest shape coming back holding no loci, and a substitution rate computed from a wrapped
remainder. Three `#[should_panic]` tests now reach all three.

**Major — `BaseComparison`'s rejection boundary was untested.** The only fixtures were a gross swap
(3 compared, 400 mismatched) and the equal case. Weakening the guard to `> compared + 1` left
everything green — and then a locus of (400, 401) reaches `substitution_rate`, which panics at the
`expect` whose comment says it cannot.

**Minor — two tests that could not fail.** `three_shards_merge_to_the_same_table_in_every_order`
compared the merge orders only with each other, so deleting the base counters from `merge` passed
it; it now compares against the table one walk would have built, two of its three shards share a
shape, and it walks all six orders rather than four. And the entry-order test could not fail for any
`BTreeMap`-keyed table — it now says so on its face (what it guards is the container being changed)
and its three loci carry different base counts, so the table equality beside it has power.

**Minor — a property test for `merge`, which states an algebraic law.** `proptest` was already a
dev-dependency. Any split of any list of up to 24 loci, over a shape space narrow enough that
collisions are common, merges back to the unsplit table in both orders.

**Seven wrong numbers, all mine.** Six were about my own code and one was a citation:

- The claim that a transposed base pair "reports a substitution rate above one" — written at three
  sites — is wrong. The reviewer bypassed the guard and ran it: `ErrorRate` refuses the rate, so the
  run *dies* inside `substitution_rate`; it is only when the bad locus is pooled with well-formed
  ones that it survives as a plausible wrong rate. The milder-sounding outcome is the dangerous one,
  and all three sites now say that.
- "different numbers by a factor of about two on real data", for loci against entries, is the HG002
  ratio (2.3) stated as the general one; tomato is 25, and both figures were sixty lines above.
- "a grid of 100,000 rates" is 99,999 points at a spacing of one hundred-thousandth.
- "the closed form outscores every point on the grid" — it *ties* the best point exactly, because
  the true maximum was chosen to land on a grid point.
- "Three shards merged in every order" walked four of the six permutations. Now six.
- The shape-space figures (220 and 293,930) were cited to spec §4.1, which contains neither; they
  are arch §2.2's.
- "a tomato sample's whole STR complement is about 1.2e9 bases" appears nowhere in the repo. The
  `f64` conversion's safety no longer rests on it: `u64 → f64` is monotone, so `mismatched <= compared`
  survives the conversion, and a correctly-rounded quotient of `a <= b` cannot exceed one.

**Also applied:** `# Panics` sections on `add_locus` and `merge`, including the sibling path's
"not atomic" note; `bases_compared`/`bases_mismatched` spelled the same way everywhere; the guard
share named as such in the doc that computes it; `table_of` taking an iterator rather than a slice;
and eight prose fixes where a property was asserted without the size the file already held.

**Declined:** turning the overflow guards into `Result`. They match the sibling path's idiom, the
release profile leaves `overflow-checks` off so the alternative to a panic is a silent wrap, and
`add_locus` is the per-locus call on a walk of millions.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | clean |
| `cargo test --lib parameter_estimation::ssr` | **70 passed** (47 after B1) |
| `cargo test --lib --bins --tests --all-features` | **3,447 passed, 0 failed, 10 ignored** |

Counted rather than recalled: `grep -c '#\[test\]'` on `stratum_table.rs` gives **37**, of which 14
are B1's, so this step adds **23**; the suite moved 3,424 → 3,447.

I re-ran two of the reviewer's mutants against the fixed suite: deleting the base counters from
`merge` now fails five tests (it previously passed the three-shard one), and weakening the
`BaseComparison` guard by one fails the new boundary test.

**Two gates are red on this branch and neither is this step's**: `cargo clippy --all-targets` fails
in four `examples/` files, and `cargo doc` reports 13 unresolved intra-doc links. Both reproduce
with this branch's changes reverted.

## Audit trail

`tmp/review_2026-08-11_ng-prepass-ssr-b2b3/` — four per-category files (reliability,
errors+naming, idiomatic+smells, numbers), the reviewed patch, and the reviewers' probes.
