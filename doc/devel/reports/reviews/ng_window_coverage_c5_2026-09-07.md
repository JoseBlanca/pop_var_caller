# window coverage — C5 review: the measurement survived, the machinery around it did not

**Date:** 2026-09-07
**Reviewed:** the working tree of plan step C5, at `086b26cc wip: C5 for review`
**Implementation report:** [ng_window_coverage_c5_2026-09-07.md](../implementations/ng_window_coverage_c5_2026-09-07.md)
**Fixes applied in:** C5's own commit

Two agents on the step's diff, each detached at it in its own worktree: one on reliability, errors
and every number the step's prose claims, one on naming, idiom, module structure and refactor
safety. **5 Major, 8 Minor.** All applied. Both agents found the same three things independently.

## Verdict

**The measurement is right and neither agent could break it.** `finish` has exactly one call site
in the library and one in the probe, both consuming their receiver, so a sample's histogram is
finished exactly once and only after both of the streaming caller's error returns. The run's
histograms are byte-identical to the whole-store walk's, and direct mode's to psp mode's, on two
independent slices — this session's six accessions over 200 kb and the reviewer's three over
400 kb. `windows_folded + windows_under_the_floor` equals the covered-position count **exactly** in
every sample of both, including the reviewer's third, where the floor silenced 377 of 301,320.

**What was weak was everything around it.** Both failure reports printed a fixed-length prefix of a
line whose cells begin at field 9 and run to field 20,058, so a wrong cell produced two
identical-looking truncations — or, through the shell, nothing at all. A new function landed inside
another function's doc comment. `into_sources` became a caller-less method that finishes every
accumulator and drops the result. And the `.histograms` path was derived in three places, in a
module whose own doc argues that both halves of a format must live together.

## Findings

### 1. Major, CONFIRMED — the failure report cannot point at the failure, in both oracles

At the shipped bin counts a fitted line is eight header fields and 50 × 401 = **20,050 cells**,
about 126 kB. The probe printed the first 160 characters of each side; the oracle diffed the two
files cut to their first 120 columns. Every cell is past both cuts.

**Measured**, by perturbing one cell of a real histogram file: the probe printed two byte-identical
lines under different labels and then panicked; the shell printed `the two modes' coverage
histograms differ` and nothing else. A second gap in the same place: with one side short a line,
the per-line loop never reached the missing samples and the report was empty too.

**Fixed on both sides, and they now name what differs** — a header field by its name, a cell by the
GC bin and depth bin it holds:

```
SRS3394606.psp  histogram-differs-at  the cell at GC bin 22 depth bin 169  run=1  walk=0
sample 0: the cell at GC bin 1 depth bin 1 is 4 from alignments and 5 from psps
```

The short-file case is now caught before the comparison, by a line count taken from the sample
columns of the VCF that same run wrote.

### 2. Major, CONFIRMED — the new function was inserted inside another function's doc comment

`histograms_the_run_recorded` landed between `windows_the_run_recorded`'s three-paragraph doc block
and its `fn`, with no blank line. The result compiles and reads as one comment: the histogram
reader was documented as owning the row format, as checking sample-index ranges and as tolerating
two rows for one locus — none of which it does — and the row reader had no comment at all. A reader
chasing a duplicate-row failure would follow that paragraph to a function that never reads a row.

**The same defect this project's review history already records once**, for `into_sources` on
2026-09-04. Fixed by moving the function below.

### 3. Major — `into_sources` became a caller-less method that pays for the histograms and drops them

Both real callers moved to `into_sources_and_histograms` in this same commit; nothing in `src/`,
`examples/` or `tests/` called the old one. What remained was a `pub` method that finishes every
accumulator — fitting the depth axis, folding up to 10,000 held-back windows, moving an 80.2 kB
counts vector a sample — and then drops all of it, with a doc line telling readers it was for "the
callers that have nothing to do with them yet".

Deleted, its two paragraphs worth keeping merged into the surviving method, which now states
outright that there is deliberately no sources-only form and that a caller writes `.0`. **The plan
says "joined by", so this is a recorded deviation**, and the reason is that the shorter name is the
one the next person adding a mode reaches for.

### 4. Major — the `.histograms` path was derived three times

Twice as the same four lines in two crate targets, once more in the shell. The module's own doc
opens with *"Both halves of the row format live here … because the reader is in another crate
target"*; the histogram half had not been given the same treatment. Drift would be loud rather than
silent — the probe panics — but its message says "a run that did not finish", which would send the
reader after the wrong cause entirely.

One `histograms_beside` now, called from the library and the probe, with a unit test over
`out/x.windows`, `windows.tsv`, `/` and `..`. The script cannot call it and says so, naming the
constant.

### 5. Major — a fitted line's cell count was never checked against its own scheme

`write_the_histogram` writes `gc_bins`, then `depth_bins`, then whatever `counts` holds. **Both
sides render through this one function**, so a `counts` disagreeing with `gc_bins × (depth_bins +
1)` produces a line no reader can index by bin — and two such lines compare *equal* when both are
wrong the same way, which is a pass in both of the comparisons this file exists for. Structurally
right today (the accumulator allocates exactly that and never resizes) and every field is `pub`.
An assertion naming the sample and both numbers, with a `#[should_panic]` test.

### 6. Major — `StoredCohortTallies`'s lists are paired by position, and the reason given was falsified by the field beside it

`per_sample[i]` and `window_coverage_histograms[i]` describe the same sample, and twenty lines
above, the same file asserts the same invariant for two other lists with the comment "a mis-pairing
here does not crash and does not go out of range". The histograms are a third list with exactly
that exposure and had no assertion — and the number they carry is the one the filter divides a
locus's depth by, so a slipped pairing is a wrong copy number for two samples rather than a wrong
line in a report.

The field's doc justified its home by saying it "is not a fact about the stored file … this is what
the *calling* run accumulated as it decoded" — which `StoredSample::read` one screen below
contradicts, its own doc reading "What **this run** drew out of the file, counted as it decoded".
Assertion added; the doc now states the pairing, says the type does not enforce it, and gives the
actual reason for the placement.

### 7. Minor, applied — prose this step made false, and names carrying the wrong meaning

- The probe's module doc still said **"nothing in the merge calls `finish` until plan step C5"** —
  this *is* C5. The expected difference survives for a different reason and now says so.
- The rendering was credited to `record_the_histograms` where the probe renders through
  `write_the_histogram`.
- `first_bytes` counted characters; `run_said` named two different things one file apart, one a
  `String` and one a `HashMap`. Both renamed.
- A stray space in a panic message; a doc link whose text named a method and whose target was the
  type; `fields[2 + 6..]` written as arithmetic the reader has to do.
- The cache's destructure comment claimed a field added to the sample window "has to be answered
  for", but the nested type was reached through `.accumulator`, so a field added *there* would have
  been dropped silently. Both levels destructured.
- The module doc's "resolved once" paragraph described the rows only; two runs in one process
  interleave their rows and *overwrite* their histograms.
- The oracle's new section had a list ranked by an aside, and could not see a histogram file that
  lost a line. Reordered, and the line count is now checked against the run's own sample columns.

## Mutations run

| mutation | outcome |
|---|---|
| the accumulator dropped instead of finished | caught — the new cache test |
| the tail windows folded a second time | caught — 12 tests, including the production-parity differential |
| `record_the_histograms` writes `sample + 1` | **survives the library suite**; caught by the probe on real data, **not** by the mode-equivalence oracle — both modes render the same wrong index |
| the overflow column omitted from every GC row | caught — the rendering's own test |
| the depth bin width written as digits, not bits | caught — the same test |
| the probe compares only the first sample | **survives** — nothing tests the example; with a difference planted in sample 1 the mutant exits 0 |
| the oracle's `fitted`-count guard deleted | **survives the comparison** — a run where no sample has a histogram passes; the guard is load-bearing |
| the cell count not matching the scheme | caught by the assertion added for it |

**What the two survivors say together.** The histogram half has exactly one guard on real data —
one probe run, by hand — and nothing in CI. The library suite covers the rendering and the
finishing; it covers neither the file the run writes nor the comparison that reads it. That is the
same standing the row half has had since C3, so not a regression, but a sample-index defect in
`record_the_histograms` ships unless someone runs the probe.

## Numbers checked

Measured on the reviewer's own slice — tomato1's first four 100 kb regions, three accessions at
26.5×, 9.3× and 12.7× mean head depth, a fact about that corner and not about the caller:

- **Run histogram = whole-store walk histogram, per sample, byte for byte**, three of three.
- **`windows_folded` + `windows_under_the_floor` = covered positions**, exactly: 398,034 + 0;
  392,756 + 0; 300,943 + 377. The third is the floor silencing about 1 window in 800 with the
  counters still summing exactly. They can sum exactly only where no `N` is observed — an `N`
  advances the frontier and becomes no centre — and none of these regions holds one.
- **What `finish` closes**: 250, 250 and 151 centres, taking `windows_folded` from 397,784 to
  398,034, 392,506 to 392,756 and 300,792 to 300,943 — **about 6 windows in 10,000**. Without it
  each histogram is short by its last half-window and the byte-identity fails outright.
- **The other half of `finish`'s job was not exercised**: every sample passed the 10,000-window
  scale sample long before its stream ended, so "fits the depth axis for a pass that ended early"
  is carried by unit tests only. It fires below about 10 kb of covered ground, or at a very sparse
  sample — this caller's own low end.
- **Mode equivalence on that slice**: 3,747 records, VCFs identical; 11,241 window rows identical,
  10,946 a measurement; histograms identical, 3 samples, 3 fitted.
- `0.25` → `4598175219545276416` correct; "two GC rows of three cells each" correct; "600 covered
  positions, floor silenced none" correct; "after the two error returns" correct; "80.2 kB a
  sample" correct (50 × 401 × 4 = 80,200).

## Looked at and found sound

`finish` exactly once per sample and at the right moment, with a failed run finishing nothing; the
three silences surviving the round trip, each with its own word, and the oracle deliberately not
comparing a non-`fitted` line away; the probe's accumulator configuration, which spells the same
six constants in a literal with no `..` so it cannot drift from the run's; the recorder inert when
off, re-confirmed by the byte-identical VCFs; the exhaustive match and the `CoverageByGcHistogram`
destructure, which make a fifth silence and a new field compile errors; and the decision to compare
as text rather than structurally — it hides a defect inside the shared renderer, which the
renderer's own tests carry, and it prevents the class that would report a *pass*.
