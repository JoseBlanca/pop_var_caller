# ng STR observations — C3: the run report reads truthfully with the slot filled

*2026-09-02. Step C3 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §3.2](../../ng/spec/run_ssr_observations.md)'s second accounting debt. Branch
`ng-ssr-observations`.*

## Plan

Filling the tract slot moved a tract's bases from *not called* to *called* with no counting
change — the ground partition sums to a hundred by construction, since every share is of the
three parts' own sum. What it broke is the **wording**, and in the direction that reads as an
improvement: a run over tract-rich ground now reports almost all of its ground as called, while
nothing scores a tract.

Two changes, and the second is the one that stops the report lying.

## Changes made

**The unbuilt line names what is actually unbuilt.** It read *"not called — repeat tracts this
caller has not built yet"*, and repeat tracts now have a generator. What is left in
`unhandled_not_implemented` is the bundle slot: clusters of repeats too close together for any
of them to have clean flanks. The line says that.

**A new line, in loci rather than bases**: *"repeat tracts built and then not called: N
locus/loci — the evidence is gathered and nothing scores it yet"*. It is printed only when N is
above zero, so a run that met no tract says nothing rather than a zero — the file's own rule
that absence and zero are different claims, applied to the report.

Both facts are needed and a reader acts on each differently: the base lines say what ground no
generator looked at, and this says what was looked at, built, merged across the cohort, and
then not scored.

## Tests added and changed

**Added** `tract_loci_the_run_could_not_score_are_a_line_of_their_own_and_only_when_there_are_some` — a run
with none prints no line, a run with seven says seven, and the same fixture's 950 called bases
are asserted beside it, because *that* is why the line has to exist: without it a run would
report a tract's ground as called when nothing was called there.

**Changed** two report tests that quoted the old wording verbatim. Both compare whole lines
rather than substrings, which is what made the change visible rather than silent.

## Validation

In the dev container. `cargo fmt --check` and `cargo clippy --all-targets --all-features -D
warnings` clean. `cargo test --lib --tests --examples --all-features --no-fail-fast` — **5,997
passed, 0 failed, 14 ignored** in the library suite; every integration target green; the three
known locus-dump failures and the psp writer bench unchanged.

## Tradeoffs and follow-ups

**The bundle line will be wrong again if a second slot is ever left unfilled.**
`unhandled_not_implemented` is one number over every unfilled slot, and today exactly one slot
is unfilled, so naming it is honest. A run report that had to name two would need the count
split per slot first.
