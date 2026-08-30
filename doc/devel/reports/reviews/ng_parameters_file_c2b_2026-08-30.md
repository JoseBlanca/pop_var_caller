# Review — C2's second half, the projection back to `RunParameters`

*2026-08-30. Three agents in isolated worktrees at `dc1e1fd2`, each handed the step's diff as a
patch: a correctness pass, a design-fidelity pass against the spec and the plan, and a reading of
the produced file as a geneticist. **1 Blocker, 7 Majors and 11 Minors applied; 4 items recorded
as owed.** Module tests 114 → 120; library 5,485 → 5,491.*

Every finding but one was about a file **nobody's run wrote** — a hand-edited one, which is the
only kind `validate` exists for. The exception is the one that matters most: **this caller could
write a file its own reader refused.**

## The Blocker, and it was the writer refusing itself

**A run with no repeat tracts writes a file `validate` turns down.** `validate` required
`repeat_tracts.slippage_group_by_read_group` to carry one row per read group over `0..n`. The
writer emits a row only for a read group the run *declared a slippage group for* — its own comment
says so — and **a run with no repeat tracts declares none at all**. That is Milestone E's defaults
run, and it is the single-sample case `CLAUDE.md` puts first among the shapes a design has to have
an answer for. Nothing caught it because no writer test calls `validate` and every projection test
started from a fixture that declares all three.

The reader was wrong, not the writer: spec §5's rule is that absence is data, and a read group
with no slippage group is a real state — `StratumFits::at` answers `UnknownReadGroup` for it and no
slippage number is ever looked up under it. That table and the substitution rate are now checked
for **naming only** the run's read groups rather than covering them, which still refuses the mirror
defect: a row keyed to a library the identity block does not list, which nothing can ever read.

## Six Majors about a hand-edited file, and each had a silent symptom

**The sample list could disagree with the read-group table, three ways.** `fitted_from.samples` is
the order every per-sample axis is indexed in, and nothing tied it to anything. A **repeat** in it
gave two rows one index and panicked several frames from the key. A **different order** did not
panic at all: the projection reads the per-sample tables by name into that list and hands calling a
vector indexed by its position, while the writer writes them in the run's own order — so a file
whose two orders disagreed gave every sample its neighbour's inbreeding coefficient and its
neighbour's sequencing batch. And a read-group row naming a sample the list did not hold — a typo —
**slipped past the new batching check entirely**, because that check looked the sample up and
skipped the row when it found nothing. One refusal closes all three: the list is the read-group
table's own samples, in the order they first appear, once each. The batching check now refuses
rather than skipping when either lookup comes back empty.

**Five repeat-tract tables dropped a duplicated row in silence.** Each becomes a map keyed by the
row's own fields, so a second row for one key replaces the first and which of the two a run scored
under is the order they happen to sit in. `StratumFits::over` carries a release-level assert
against exactly this on the fit's side, on the stated grounds that two levels can differ by a
factor of five — and a file is the one input path where a person copying a row to edit it can
produce one. Refused now, quoting the repeated key back, which the projection could not have done
having already lost one of the two.

**An evidence count's unit was dropped and re-minted, and the comment claiming otherwise cited a
test that does not exist.** `Estimate<T>`'s count is a bare number whose unit follows the quantity;
the file names the unit because a reader cannot be sent to the source to find it. So a calibration
row saying `covered_positions = 812344` validated, projected to 812,344, and was written back as
`reads = 812344` — a key the user typed, silently replaced, in the one direction spec §1.2 goal 1
forbids. The three units differ by orders of magnitude on one cohort. `validate` now checks each
quantity's count against the unit that quantity is fitted over, which is the check the comment had
already claimed existed.

**The fallback concentration's warrant could contradict the file's own rows.** The bottom rung
states the median of the concentrations this run's strata fitted wherever any was fitted, and the
compiled-in flat constant only where none was — so `defaulted` beside a non-empty
`length_spectrum_by_stratum` is a claim the file's own rows refute. It had to be refused rather
than left, because nothing downstream can see it: the projection carries only the number and the
writer re-derives the warrant, so a contradiction was rewritten on the way out rather than
reported.

## The reader's Blocker: a correction that left a number with no description

The fixture correction that moved the fallback concentration from `1.25 defaulted` to
`3.5 fitted_here` was right about the mechanism and **took the key's only prose with it**. That
note fired only for a defaulted value, and it was where the file said what the number *is* — "this
many chromosomes' worth of belief, spread flat over whatever lengths a tract offers" — and when it
is used. A `fitted_here` fallback was left with a bare `#` above it and, since it carries no
`observations`, read against the file's own header as a number that was *not measured*: exactly
the state the correction set out to remove. The description is now in the section's own note, which
every run writes, and the per-key note says only what a *defaulted* one adds.

Four more of that pass's findings were applied because they were cheap and false as written: the
`rung` key described four regimes in prose and spelled none, twenty lines above a key that spells
all four of its own; the substitution-rate note named two of that table's three keys and omitted
the one that collides with the file's top-level `ploidy`; "a pair with no row" was one axis short
of the table it described; the calibration note explained a multiplier above one and at one, with a
row at 0.87 beneath it; and the batching invariant correction 1 enforces was nowhere in the file,
two sections after a note arguing that two lanes of one plant must be kept apart.

## Two fixture numbers that could not be what they said

`strata = 12` on a share curve fitted from 2 to 4 repeats. A curve `curve_fitted_on = "this_period"`
was fitted through that period's own strata, and a period has one stratum per reference repeat
count — three of them between 2 and 4. `centre_repeats` sat at 11.5 on the same curve, and `cells`
at 23 on a level curve over eight repeat counts. All three now follow the range they are printed
beside, so a reader can check them.

## Eleven Minors, and the two worth naming

Two cohort-sized linear scans sat inside per-row loops — 9 million string comparisons at the 3,000
samples `CLAUDE.md` commits to, twice over — and now share one index built once. And the outlier
weight's match had a wildcard arm where every other conversion in the file is exhaustive, so a
fifth warrant would have arrived silently as `supplied`; it is spelled out.

The rest: a duplicated copy of `UNMEASURED_READ_GROUP` whose own doc said it *had* to equal the
original, now imported rather than copied; dead `Result` machinery around infallible code; a header
that listed the wrong newtype refusals (a ploidy has no upper bound, a rise shape and a repeat
count do, and two of the refusals written are unreachable); a banner counting eight conversions
where there are six; a claim that the module hands the identity block back, which it does not; and
a test whose name promised both directions where it asserts one and its sibling asserts the other.

**A test for the bottom of the committed range was added rather than traced.** Every projection
test started from the three-read-group, two-sample, four-stratum fixture; the one-sample,
no-repeat-tract, uncontaminated file now has one of its own, and it is the only test that exercises
an empty repeat-tract section.

## Recorded and owed

- **The spec is silent on two rules this step now enforces.** §5 or §9 should say that a sample's
  libraries all ran in one batch, and that a `defaulted` number carries no evidence count — §2
  currently says the opposite of the second in passing. Both are invariants of the objects the file
  projects onto rather than inventions, but the file's reader should not be the only place they are
  written down.
- **D3 will trip on the fallback concentration.** Spec §2.1 demotes a mismatched file's numbers to
  `Supplied` wholesale; this one's warrant is not carried through `StratumFits` at all, it is
  re-derived by the writer from the same strata. A demoted file written back out re-emerges as
  `fitted_here`, so §13 test 5's "every warrant `Supplied`" cannot pass on it without changing the
  writer or the spec.
- **C5 is not made easier by the three §5 assertions this step added.** They assert a *state* —
  `contamination_is_absent()`, `was_measured()`, `Provenance::Defaulted` — where C5's brief is a
  fixture per row in which collapsing the two states **changes an answer**. Two of the five rows
  have no lookup-level assertion of any kind yet.
- **The reader still cannot tell two curves apart.** The fixture's two period-2 slip curves have
  identical coefficients and different fitted ranges, and its two share curves likewise, while the
  section's prose says a period has one curve. Making that true means restructuring which strata
  the fixture holds, and it belongs with C5's rebuild rather than here.
