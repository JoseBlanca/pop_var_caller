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
| [`scoring_context.rs`](../../../src/ng/run/paralog_filter/scoring_context.rs) | 357 lines: `ParalogScoringContext`, `WhyNoCoverageModel`, `CohortSizeMismatch`, `CoverageFitConfigRefused` |
| [`scoring_context/tests.rs`](../../../src/ng/run/paralog_filter/scoring_context/tests.rs) | 580 lines: 19 tests |

**Those are the sizes after this step's review**, whose fixes are folded into the same commit; as
first written the module was 231 lines with 10 tests.

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
The method takes an index, and indexing is where a wiring error turns into a crash on someone's
whole-genome run rather than a diagnosis.

**Corrected at review:** as first written, *nothing* checked a record's width against the cohort's
at all, so a record narrower than the run would have scored on a prefix and a wider one been
truncated — both silently, since the scorer answers a length mismatch with a *neutral score*
rather than an error. `observations_of` now checks once per record and fills a caller-owned
buffer; `observation_of` stays per sample beneath it.

**5. The length disagreement between the two input slices *does* panic.** They come from the same
run in the same sample order, so a mismatch is a wiring error the caller cannot meet, and the
alternative is a scorer that silently drops the tail of the cohort.

**6. A refused fit configuration fails construction rather than every sample** (added at review).
The copied fit re-validates its configuration on every call, so one wrong knob would have reported
that the coverage model rested on 0 of 63 samples with 63 identical reasons — an operator's
mistake dressed as a verdict about the cohort. `new` returns a `Result` for that one case.

**7. `new` takes `&[InbreedingF]`, the validated newtype the parameters file exposes** (changed at
review), and unwraps it once. The copied scorer wants plain `f64`s because it is production's code
unchanged; this is the one boundary between the two, so it is the one place the newtype comes
off — rather than every caller stripping it.

## Tests

**19** — ten as first written, nine added by the review.

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
number, which is the footprint the filter hunts, and the fixture's claim that its `ONE_COPY_DEPTH`
is one copy is asserted rather than left in a comment. **Tightened at review from 0.1 to 1e-9**:
the fixture is exact — its depth counts are symmetric about the peak bin and its GC curve is a
single `[1.0]` — so a tolerance a hundred million times looser than the fixture's own precision
was letting the fitted scale drift by 9% while the constant's name went on claiming otherwise.

**Added at review, and the first of them closes a Blocker:**

- **`the_relative_copy_number_uses_the_windows_own_gc_and_not_a_constant`.** Every test above uses
  a one-GC-bin histogram, and the copied model returns `gc_bias_curve[0]` for every GC value on a
  one-bin curve — so the GC half of the copy number was unobservable and **replacing the GC
  argument with a literal passed all ten tests**. A wrong GC there does not panic and does not give
  `NaN`; it gives a plausible copy number, wrong by the bias multiplier, on every sample of every
  record. The new four-GC-bin fixture has two samples at the same depth and different GC read as
  different copy numbers.
- **`a_window_with_only_one_field_absent_is_absent`.** Both guard fixtures set *both* fields to
  `NaN`, so either half alone caught them and either could have been deleted without a red test.
  The GC half is load-bearing against a **panic**: a `NaN` GC reaches `gc_multiplier`, both range
  comparisons are false, the floor saturates to zero, and the interpolation indexes past the end of
  the curve — inside a file the copy guard forbids editing.
- **`an_infinite_window_depth_is_absent_rather_than_a_winsorised_paralog`**, which pins the
  deliberate divergence from `WindowCoverage::is_absent`: that predicate accepts `±∞`, and an
  infinite depth would divide to `+∞` and be winsorised into a confident four-copy paralog.
- Zero depth reads as zero copies rather than as an absence; a sample far above its one-copy depth
  is handed over **unwinsorised**, because the cap is the scorer's; a record mixing covered and
  absent samples hands over only the covered ones; a record whose width disagrees with the cohort
  is refused in both directions; the record-level buffer is refilled rather than appended to; and a
  refused fit configuration fails construction.

**Five mutations, all killed, each reverted from a byte-compared backup.** Run in the container:

| mutation | what happened |
|---|---|
| the tract arm returns `None` instead of `(0, 0)` — spec §6 trap 1 exactly | 1 of 19 fails: `every_sample_of_a_tract_scores_even_where_no_read_supports_an_allele` |
| the `total == 0` skip deleted from the generic arm | 1 of 19 fails: `a_generic_locus_sample_with_no_reads_is_absent_from_the_score` |
| the GC argument replaced by a constant | **1 of 19** fails — it survived all ten before the review |
| the GC half of the finiteness guard deleted | **2 of 19** fail — it survived all ten before |
| the depth half deleted | **2 of 19** fail — it survived all ten before |

## Validation

Run in the dev container, on this tree.

    cargo test --all-features --lib "ng::run::paralog_filter::scoring_context"
    → test result: ok. 19 passed; 0 failed; 0 ignored; 6559 filtered out

    cargo test --all-features --lib --bins --tests
    → 6,563 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 9 errors: 3 needless_lifetimes in cohort_merge, 6 in src/ng/window_coverage/

    cargo fmt --check
    → dirty on 9 files, none of them this step's

The integration failure and the three `cohort_merge` lints are `main`'s at `a33ada0f`; the six
`window_coverage` lints and the nine files come in with the merged branch and fire there
identically. The gate — no failure and no lint this branch's merge base does not already have —
holds, and the lib count moves 6,544 → 6,563 (6,554 as first written, plus the review's nine).

## Tradeoffs and follow-ups

- **Nothing calls this yet.** Pass two is step C3, which is where the ratios are computed and
  folded into the histogram; C2 is the sink that fills the spill in the first place.
- **The fit config and model params are taken as arguments and defaulted by the caller.** Spec
  §3.6 leaves the window and histogram parameters as compiled defaults until a measurement says
  one needs to move, so nothing exposes them on the command line yet.
- **~~`observation_of` looks up the sample's row twice~~ — fixed, and the cost was measured rather
  than deferred.** One `match` now yields the window and the counts together, which also makes both
  row types destructured exhaustively. The review measured what the duplication was worth:
  `observation_of` costs **1.82 ns a sample** against **2,828 ns** to score the same sample —
  about 1,550× — so the second lookup was under a ten-thousandth of pass two and could never have
  shown up in step D2's shares. **The reason to remove it was the missing length check, not the
  cost**, and this report's earlier deferral to D2 was answering the wrong question.
