# Review — the two names the key revision left orphaned

*2026-08-30. Two agents in isolated worktrees at `d20a6bce`, each handed the step's diff as a
patch: one reading the produced file as a geneticist who has never seen the code, one checking
the rename for completeness and over-reach. **1 Blocker / 4 Majors applied, 12 findings
recorded and left**, all of them older than this change.*

The renames themselves came back clean from both passes. The reader **guessed all three keys
correctly on first read and the rows confirmed all three**; the correctness pass found no
occurrence of an old spelling left where it matters, nothing renamed outside the module, and no
assertion made vacuous. What the rename broke was four comments that had kept the old words.

## What the change was

| was | is |
|---|---|
| `level_origin` | `share_of_reads_that_slip_origin` |
| `shares_origin` | `shorter_share_and_fall_off_origin` |
| `share_of_reads_from_elsewhere` | `share_of_reads_from_another_sample` |

The Rust types `LevelOrigin` and `SharesOrigin` did **not** move. They mirror the pre-pass's
`LevelProvenance` and `SharesProvenance`, which is the separation `Warrant` already keeps from
`Provenance`: the file's spelling is a compatibility surface and the in-memory type's name is
not. Both passes were asked to judge that and both accepted it; the correctness pass checked
that no doc comment left behind names a key that does not exist.

## Applied

**B1 — "a share curve also records `curve_fitted_on`" became false half the time.** Before the
rename the file held a `level_origin` and a `shares_origin`, so *a share curve* was bound: it
meant the curve under the second. The rename deleted the word *level* from the file and put
*share* into the first key's name, and the file holds two kinds of curve — the one under
`share_of_reads_that_slip_origin`, which counts `cells` and has no `curve_fitted_on`, and the
ones under `shorter_share_and_fall_off_origin`, which count `strata` and do. The reader spent
minutes deciding whether the missing key on one row was a defect or a category they had not
understood. The paragraph now names which curve carries it and says the other is a different
fit and has none.

**M1 — the module header still told a reader to write `fraction = 0`.** `mod.rs`'s five-state
table, which is the header's answer to spec §5, spelled the contamination fraction `fraction` —
a key `ContaminationMeasurement`'s `deny_unknown_fields` refuses at parse time. Stale since the
*previous* rename, and the more misleading of that rename's two missed sites, because it is the
row a reader consults when asking how the file says *measured and found clean*. The TOML sketch
53 lines below it was fixed in this step's first draft; the table was missed until this pass.

**M2 — the reworded paragraph named neither of the two keys it was about.** The first draft said
"the two `_origin` keys record where those numbers came from". To learn what they are called a
reader has to go to the row that carries them, which is **871 characters wide**. The paragraph
before it names its three subjects inline; this one now does the same, and says which key covers
the slip share and which covers the other two.

**M3 — "each of the two carries its own `expected_slipped_reads`" is false on one of the three
rows.** Row 1 of the fixture (period 2, 11 repeats) carries one origin key and no count at all,
because both numbers came whole off a curve. The paragraph now states both absences: a row whose
two shares were not fitted here has no `shorter_share_and_fall_off_origin` key, and an origin
carries no count where its number was taken whole from a curve.

**M4 — a mechanically renamed comment described a path that never existed.** `validate.rs`
records why one test compares whole refusal paths rather than last segments: the code once
emitted `shares_origin.shorter_share` where the key is `shorter_share_smoothing`. `perl` rewrote
the first half of that historic string, leaving a path neither the old code nor the new one ever
produced. Reworded to state the defect rather than quote a path.

Also applied: `SlippageRow`'s two field docs still called the number *the level* (the word the
previous rename removed from the file); `mod.rs`'s header pointed at `tests/testdata/every_shape.toml`
where the golden lives under this module; and `PROJECT_STATUS.md` still carried both renames as
open rulings.

## Recorded and left — twelve findings, none of them this change's

The reader was asked to read `every_shape_as_written.toml` as an artefact. It is not one: it is
the golden copy of a fixture built to exercise **every shape**, which is why `mod.rs`'s header
says to read the sibling file "as a record of the key names, not as the file a run will
produce". Four of the findings are that fixture holding combinations no single run can:

- **the fallback concentration is `1.25` marked `defaulted` in a file that fitted a stratum**,
  where the comment beside it says a run that fitted any stratum states the median of its own
  and marks it `fitted_here` — and the median of the one fitted concentration in the file is
  3.5, not 1.25;
- **a read group with 3 reads in 100 from another sample, alone in its sequencing batch**, where
  the batch is defined two sections up as the population a contaminating read is drawn from. The
  old name (*from elsewhere*) covered a contaminant outside the cohort; the new one does not, so
  the rename is what made this checkable;
- **sample TS-1's two read groups in two different batches** while the per-sample table puts
  TS-1 in one of them, with no comment saying whether that is legal;
- **`warrant = "defaulted"` beside `observations = { reads = 4 }`**, which the projection's own
  documented rule forbids — a defaulted value carries no evidence count.

**These belong to step C5**, which owes one fixture per row of spec §5, each built so that
collapsing the two states it separates changes an answer. A fixture a run could actually have
produced is what that step needs and what this one does not have.

The remaining eight are comment gaps older than this change, owned by B3, and are recorded here
rather than fixed because none was made worse by the rename:

- the four warrants are listed in the header and **defined nowhere**, and *borrowed* is the one a
  reader cannot guess — borrowed from what? The file gives no donor;
- `rung = "fitted_curve"` names four regimes in prose and gives the spelling of **none**, twenty
  lines above `curve_fitted_on`, which gives all four in parentheses;
- `fitted_from_reads_of`'s two values are never mentioned in any comment;
- the substitution-rate row is keyed on three things and its comment enumerates two, omitting
  `ploidy` — which also appears at the top of the file, with nothing saying whether the two may
  differ;
- `span` is used to define `shares_by_repeat_offset` and is not a key anywhere;
- the per-period spectrum has no `reference_repeats`, so the anchoring rule stated for that key
  is undefined for one of the two arrays that carries it;
- `built_in_default` names a stated constant with no origin, where spec §8 requires the origin
  beside it — the same gap B3 recorded as unfixable, since a slippage number carries a smoothing
  origin and no warrant;
- the calibration comment explains a multiplier above one and at one, and the file has a row at
  **0.87**.

## Line width

No comment line in the produced file exceeds 80 characters. Worth knowing before the next word
is added: the longest is exactly 80. Counted with `perl -CSD`, because `awk 'length>80'` counts
bytes and reports five false positives on the lines holding an em dash.

## Tests

103 module tests, 0 failed, unchanged in number — this step adds no test. Both golden files were
regenerated by the two `#[ignore]`d `regenerate_*` tests and their diffs read line by line; the
whole of the first is the rename, and the second is the rename plus the repeat-tract section's
comment. The produced file goes from 176 lines to 183 and from 93 comment lines to 100 — the
seven are what saying which origin is which, and which curve carries `curve_fitted_on`, costs.
