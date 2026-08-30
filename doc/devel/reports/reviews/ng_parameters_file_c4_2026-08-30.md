# Review — C4, the north-star round trip

*2026-08-30. One agent in an isolated worktree at `14e0e3bb`, handed the step's diff as a patch.
**1 Blocker, 5 Majors and 10 Minors applied.** Module tests 130 → 134.*

The pass found **no correctness defect in the production code** and re-derived all three of the
step's measured numbers independently, finding all three right. What it found instead was that the
one thing the step's own justification rests on was not true.

## The Blocker: the fixture was not shaped like a fit's output

The owner's ruling of 2026-08-30 accepts a synthetic run in place of a real fit **because it is
built to look like one** — the maps keyed as the pre-pass keys them, through the fit's own doors.
In four places it was not.

**A period's strata contradicted each other about whether that period has a curve.** The level's
origin and the outcome's variant were both keyed on `which_repeats % 3`, so the fitted strata took
`LevelSource::Cell` — which `blend_level` emits only where *"there is no curve at this period"* —
while the derived strata beside them existed only because `derive_thin_strata` found one. And
because the third arm went to a refusal, whose provenance is discarded, **`LevelSource::Blend` —
the fit's ordinary case, a stratum with its own evidence weighed against its period's curve —
appeared nowhere in the file at all.** The refusals were `BelowTheFloor`, which is exactly what a
period's curves turn into a derived stratum, so they could not have survived beside them either.

The fixture now decides by period: **the first has no curves**, so its strata are fitted or refused
for want of a spanning read, its levels come from the cell and its shares from the stratum; **the
other two have both curves**, so a stratum with evidence blends and one without is derived whole.
All three level sources and all three share sources now appear, each in a period that can produce
it.

**Both shares' origins were a state no fit produces**, on every row that had one: one share said
its period had no curve while the other, three keys along, named one. And on a derived stratum the
fixture claimed a blend, which by construction cannot have happened —
`derive_thin_strata` writes `Curve` for both shares with no slipped-read count, always.

**The contamination mix was unreachable.** `TooFewMarkers` stamps every read group of every sample;
the only per-unit refusal is `OwnFrequencyIsItsOwnEcho`, and its own comment says the refusal is
the *sample's* and refuses every library it has. The fixture refused one library of a sample and
estimated its sibling. It now refuses a whole sample, and the `leverage` it writes is below the
`MAX_LEVERAGE` above which the fit refuses rather than emits.

**And spec §5's third row was staged without its distinguishing case.** Every read group's minted
mean was built from its own fitted rate, so every multiplier was 1.0 to within the accumulator's
quantum — the column where a transposed axis would hide best — and no calibration was ever
`Defaulted`. One library's fitted rate is now zero, which is what assembly refuses a scale for, so
the file carries a defaulted 1.0 beside fitted ones that differ from each other.

## The size measurement was right and its explanation was wrong

The two row widths — **157 bytes an inbreeding row and 185 a substitution-rate row** — were
re-derived line by line and both check out, as do `(108, 144)` slippage answers and `(24, 16)`
strata. What was wrong was the attribution: the comment blamed the 2026-08-30 key renames, and
**no key either row carries was touched by them**; both rows are byte-identical in shape to the
day spec §9's 146 bytes was written. The whole difference is how wide the numbers and names inside
them are.

So the figure is now stated as a floor, with the two things that push a real cohort above it: a
`read_group` is one digit here and four at 3,000 libraries, and a `bases_compared` count is nine
digits where the fit's own per-read-group total on HG002 is 172,616,054. **§9's 146 reproduces
nothing in the tree**, including the fixture that existed when it was written, so what the test
pins is an order of magnitude rather than a measurement against a measurement.

## Two claims about what the tests cover

**The reason given for the second test was wrong.** It said a file comparison cannot see two strata
whose rows were exchanged — but a slippage row's key is written *from* the map key, so an exchange
writes B's key beside A's numbers and the first test catches it. What the second test genuinely
adds is the **lookup**: `StratumFits::at`'s four answers for a cell with no numbers, the three
rungs of `length_spectrum_at`, and the bottom rung's median. None of that is written anywhere in
the file. The header says so now, and the rung ladder is counted so that the walk cannot compare
one rung with itself thirty-six times.

**"Every accessor the run exposes is compared" was false**: the sequencing batching was covered
only through the file. It is compared directly now.

## What the trip does and does not reach

The pass traced it: **the only field in the whole file where a value lost on the way *in* is
invisible to the six-stage trip is the tract ladder's fallback concentration's warrant**, which the
writer re-derives from whether any stratum was fitted. `validate` refuses the contradiction
instead, which is the guard that stands in its place. The identity block is fed back into the
second write rather than round-tripped, because `to_run_parameters` carries none of spec §6's four
bindings — they are D2's. Both are now said in the module header rather than left to be inferred.

Everything else round-trips, **including the two things `RunParameters` drops**: the calibration's
evidence count and the inbreeding warrants both come back through the projection's own bundle and
are fed to the second write from there, not from the original.

## Ten Minors

The fixture's number generator claimed every value was a function of every index that keys it, and
was linear with small weights — so the 72 substitution rows carried only 53 distinct values, and
three of the numbers were not functions of all their keys. It mixes its indices now. `borrowed`
was a per-group vector on a field that names repeat counts. The length-spectrum walk ran six times
identically inside the read-group loop. `the_reach_of` re-implemented `SlippageCurve::reach`, which
is the exact defect a previous rename found in the other fixture. The row-size helper would have
averaged comment lines in. The rates helper was duplicated character for character between two test
modules in one file. "Thirty-one calls to `assemble`" is thirty-three. And the tract ladder's middle
rung had **no** rows at all: the run now pools a spectrum for each period it drew curves for, so
all three rungs of the ladder are exercised.
