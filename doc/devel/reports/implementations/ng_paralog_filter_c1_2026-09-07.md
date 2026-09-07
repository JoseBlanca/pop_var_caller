# The hidden-duplication filter — C1: the scoring context, and the skip that must not reach a tract

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone C, step C1 — *own commit, do not bundle*
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.1, §3.2, §6 traps 1 and 3
**Branch:** `ng-paralog-filter`

## The answer

**Everything the scorer reads that does not change from record to record is now built once, and
the one rule that had to differ between a repeat tract and everything else is a type rather than
a condition someone has to remember.**

Between the calling pass and the scoring pass there is a gap where the coverage histograms are
finally complete. This is what is built in it: one fitted coverage model per sample, the
inbreeding coefficient the parameters file carries for each, the σ₀ slice the scorer takes beside
them, and the per-pass tables it precomputes from the cohort's size.

**The step's silent failure is production's zero-read skip reaching a repeat tract.** Production
drops a sample whose total reads are zero, which is right where reads report their allele and
wrong at a tract, where the filter hands the scorer no read counts at all because slippage has
smeared the split (spec §3.2). Copied without that distinction, every tract in the run scores on
nothing, is never flagged, and the file looks correct. **The reshaped spill entry is what closes
it**: a tract's rows carry no read counts, so there is no total for the skip to test, and
`observation_of`'s match arm for a tract cannot reach the skip at all.

## What was added

| file | what it is |
|---|---|
| [`scoring_context.rs`](../../../src/ng/run/paralog_filter/scoring_context.rs) | 231 lines: `ParalogScoringContext`, `WhyNoCoverageModel` |
| [`scoring_context/tests.rs`](../../../src/ng/run/paralog_filter/scoring_context/tests.rs) | 325 lines: 10 tests |

[`paralog_filter/mod.rs`](../../../src/ng/run/paralog_filter/mod.rs) gains the declaration and two
re-exports.

## Assumptions and recorded deviations

**1. `new` takes `Vec<SampleHistogram>`, not `&[Option<CoverageByGcHistogram>]`.** The spec's type
block says "one fitted model per sample, `None` where the fit was rejected", which describes the
output. The *input* is what the merged window-coverage work hands back, and it already carries the
three named reasons a sample has no histogram — no window finalised, every window under the floor,
no positive median depth. Flattening those to `None` at the boundary and then reporting "61 of 63"
would throw away exactly the distinction the run report wants.

**2. So `WhyNoCoverageModel` has four variants, not one.** Three are the histogram's own absences,
carried through; the fourth wraps the fit's error, which distinguishes an essentially uncovered
sample from one whose depth ran off the top of the histogram — opposite problems a bare count
would merge. Spec §3.5's run-report line asks for "samples whose coverage model was rejected, each
with its reason", and this is what makes that line writable.

**3. It is not `Clone`,** because `CoverageModelError` is not, and that type is a byte-for-byte
copy of production's that the copy guard forbids editing. Callers match on it instead.

**4. `observation_of` returns `None` for a sample index past the cohort rather than panicking.**
The record and the context both claim the run's sample count, so nothing inside the run can reach
it — but the method takes an index, and indexing is where a wiring error turns into a crash on
someone's whole-genome run rather than a diagnosis.

**5. The length disagreement between the two input slices *does* panic.** They come from the same
run in the same sample order, so a mismatch is a wiring error the caller cannot meet, and the
alternative is a scorer that silently drops the tail of the cohort.

## Tests

**10.**

**The one the step exists for:** `every_sample_of_a_tract_scores_even_where_no_read_supports_an_allele`
— three samples at a tract, all with windows, none with read counts, and all three must produce an
observation. Beside it `a_generic_locus_sample_with_no_reads_is_absent_from_the_score`, which is
the same rule's other half and production's behaviour.

**On what reaches the scorer:** alternatives are summed, so a sample with 4 reference and 8
alternative reads arrives as `alt_reads = 8, total_reads = 12` — the pooling that lets a
multiallelic site take this path at all. A sample with no usable window is absent whatever its
reads say, at both kinds of record. A sample whose coverage model was refused is absent from every
record rather than scored as neutral.

**On the bookkeeping:** each of the four ways a sample can lose its model is kept apart from the
others and the fit's own reason survives; the σ₀ slice is the cohort's length with `NaN` — not
zero — where a sample has none, because the scorer returns a *neutral score* rather than an error
if that slice is short, so a wrong length would quietly score every record as saying nothing.

**On the coverage signal end to end:** doubling a sample's window depth doubles its relative copy
number, which is the footprint the filter hunts — and the fixture's own claim that its
`ONE_COPY_DEPTH` is one copy is asserted rather than left in a comment, since the ratio test alone
would pass on a model whose scale was anywhere at all.

**Two mutations, both killed, both reverted from a byte-compared backup.** Run in the container:

| mutation | what happened |
|---|---|
| the tract arm returns `None` instead of `(0, 0)` — spec §6 trap 1 exactly | 1 of 10 fails: `every_sample_of_a_tract_scores_even_where_no_read_supports_an_allele` |
| the `total == 0` skip deleted from the generic arm | 1 of 10 fails: `a_generic_locus_sample_with_no_reads_is_absent_from_the_score` |

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::run::paralog_filter::scoring_context"
    → test result: ok. 10 passed; 0 failed; 0 ignored; 6559 filtered out

    cargo test --all-features --lib --bins --tests
    → 6,554 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 9 errors: 3 needless_lifetimes in cohort_merge, 6 in src/ng/window_coverage/

    cargo fmt --check
    → dirty on 9 files, none of them this step's

The integration failure and the three `cohort_merge` lints are `main`'s at `a33ada0f`; the six
`window_coverage` lints and the nine files come in with the merged branch and fire there
identically. The gate — no failure and no lint this branch's merge base does not already have —
holds, and the lib count moves 6,544 → 6,554.

## Tradeoffs and follow-ups

- **Nothing calls this yet.** Pass two is step C3, which is where the ratios are computed and
  folded into the histogram; C2 is the sink that fills the spill in the first place.
- **The fit config and model params are taken as arguments and defaulted by the caller.** Spec
  §3.6 leaves the window and histogram parameters as compiled defaults until a measurement says
  one needs to move, so nothing exposes them on the command line yet.
- **`observation_of` looks up the sample's row twice** — once through `samples.window(sample)` and
  once in the match arm. At three thousand samples times the record count that is one extra bounds
  check per sample per record; whether it shows up at all is what step D2's pass shares say.
