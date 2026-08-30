# Review — C5, the five states of spec §5

*2026-08-30. One agent in an isolated worktree at `3bef74ae`, handed the step's diff as a patch.
**0 Blockers, 3 Majors and 6 Minors applied.** Module tests 134 → 139; five mutations run against
the projection, one a row, and each failed exactly its own row's test.*

The mapping onto §5 came back clean — five tests, five rows, none doubled, none tested through a
proxy, and none vacuous. What the pass judged instead was the brief's real demand: whether each
collapse **changes an answer**.

## Row 5 was asserting through a function calling never runs

`Slippage::read_probabilities` is called by the parameter fit and by that test, and by nothing in
`src/ng/calling/`. What the caller actually does with the two states is take
`StutterModel::hipstr_shipped` where the lookup fails and `stutter_model_for` on the numbers where
it succeeds — so **the number a read is scored against is the model's direction share: 0.05 under
the shipped model, and exactly 0 under a slip rate of zero**, because every direction share is
that rate times a share of it. The test computes both now, through the caller's own conversion.

**And its doc claimed more than the caller does.** It said a read one repeat short becomes
impossible and the genotype that would have explained it is ruled out. The emission does collapse
to exactly zero — but the row then charges the read to the outlier term, one in a hundred by
default, so the genotype survives with a penalty. The honest claim is strong enough: the emission
is nothing, and the only thing left to explain the read is junk.

## Row 3's answer could not fail on the thing it named

It closed by asserting the warrant the run *writes*, and `of_run` copies that straight off the
calibration with no logic between — so the assertion restated, as a `Warrant`, the `Provenance`
two lines above it, and the whole-fixture round trip already covered the write.

The answer that genuinely changes is what a **consumer** does with the warrant. Spec §2's rule is
that consumers combine warrants rather than branching on them, so a call resting on one fitted
parameter and one defaulted one is a defaulted call — and `summarise_condition`'s fold, **the only
place in calling that reads this field at all**, is where that lands. The test folds it now.

## What the pass confirmed about the three report-only rows

**Both ceilings are real.** No score is reachable for row 1: `likelihood::ssr`'s own test already
asserts that the three-term form at a fraction of zero equals the two-term form *bit for bit*. And
`ContaminationView::was_measured` is read in exactly three places, none of which reaches a score, a
genotype or a quality. So the run's report is the strongest thing that moves for both, and it is
what the tests assert.

The header said the three rows "change nothing it computes" while row 1's own doc said the answer
is which formula the read likelihood runs. Both are now one claim: the branch is taken, and it
provably returns the same number — which is exactly why the warrant and the report are the only
things that separate these states, and why a reader who inferred either from the value would have
nothing to notice.

## Six Minors

Two of the fixture indices carried two meanings at once — a position in the file's vector and a
read-group id in the projected view — and agreed only because the fixture happens to list its rows
in read-group order, which is the coincidence that broke eight tests when the slippage rows moved.
Both find their row by its own identity now, and the third uses the named constant.

Row 5's collapsed row said no read slips beside twelve thousand expected slipped reads, which is a
second contradiction on top of the collapse it models; it carries no count now. Row 4's collapsed
file leaves a fallback concentration that is no longer the median of the strata beside it, which
`validate` does not check and this test does not read — said in a comment rather than left. And
two helpers were copied verbatim from the sibling test module; one is shared and the other, which
dropped the evidence counts its own doc said mattered, is gone.

## What the tests still do not catch

Named rather than fixed, since each is another row's or another step's: the sequencing batching
when the contamination table is absent; a swap of the two evidence counts, which `was_measured`
only tests for being above zero; `Supplied` collapsing into `Borrowed`, which rides on the
whole-fixture round trip; the flat rung of the length-spectrum ladder, which this module's fixture
never reaches; and a per-stratum row padded to more slippage groups than the fit was run over,
which is invisible to both the file and the lookup.
