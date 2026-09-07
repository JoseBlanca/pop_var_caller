# window coverage — D2 review (naming, structure, prose): the sweep took D1's shape and left D1's paperwork behind

**Date:** 2026-09-07
**Reviewed:** commit `f3228dff wip: D2 for review`, diffed against `38e06a3f`
**Remit:** naming, idiom, code smells, refactor safety, module structure, and the prose the step
ships. Correctness, error handling and the numbers are a separate reviewer's.
**Implementation report:** [ng_window_coverage_d2_2026-09-07.md](../implementations/ng_window_coverage_d2_2026-09-07.md)
**The step before it, and its two reviews:**
[D1's report](../implementations/ng_window_coverage_d1_2026-09-07.md),
[D1 design](ng_window_coverage_d1_design_2026-09-07.md),
[D1 correctness](ng_window_coverage_d1_correctness_2026-09-07.md)
**Fixes:** applied in this worktree, listed at the end.

## Verdict

**The question this review was set — did D2 follow D1's rulings or quietly diverge — has a split
answer, and the split falls exactly along the line between code and prose.**

**In the code, D2 followed every ruling that had been made, and made the same judgement calls
again.** `BinSchemeSweep::close` takes `self`, consumes the arms and hands back the answers, so
there is no half-closed sweep, no `Option` around an accumulator and no `expect` in a hot path —
which is the shape D1's review had to argue the previous sweep into. The sweep rides inside
`WindowRecomputation` rather than becoming a third `RecordMeasurement`, for the same reason and
with a stronger one besides: the arm it checks has to be compared against the walk's own histogram,
which exists only inside the recomputation. The printed rows tag every column they carry
(`median-depth=`, `folded=`, `overflowed=`). The word *arm* is defined before the fourth question
uses it, because the third question defined it.

**In the prose, three of D1's four rulings were re-broken, in the same places and the same
forms.** The module doc's own statement of how the file is organised — the paragraph D1's review
rewrote to say the third question rides inside the second — still says *the third*, and a reader
who takes it at its word does not know a fourth measurement exists. The section that introduces the
fourth question ends "All three of those constants are provisional", while the constants themselves,
edited by this same commit, now say the opposite. And the sentence that carries the ruling —
"the range has 780 times the margin the worst store needs" — carries a ratio of *overflow shares* in
the language of *range width*, which is the same defect D1's review found in "moves by 6 in 100",
one step earlier, in the same three-document set.

**Two numeric claims are wrong and both are contradicted by the step's own evidence.** "Nine of the
ten sample-stores" fitted a shallow median from the prefix — the report's own table two lines below
shows eight, with two above 1.0, and the raw walk output confirms eight. "The safe band is 5 to 40"
reads as measured at both ends: no setting tried, 2.5 included, actually made the fit reject a
sample, and 40 is where the settings stop rather than where a cost appears.

**One place D2 should have diverged from D1 and did not.** `WindowsUnderEachFloor` earns its own
type because it is printed *and* summed across stores. The bin-scheme answers are printed and never
summed — an axis fitted from one sample's depth has no total — and nothing said so. That silence
matters more here than it would have in D1, because `WindowTotals::add` reads its store by field
access, so a per-store answer that ought to be totalled and is not is a zero row, not a compile
error.

**13 findings, all applied. 7 Major, 6 Minor, no Blocker.**

## Findings

### Major 1 — the module doc's "One walk, several measurements" describes a file with three questions

The section is the file's own account of how it is organised. D1's review rewrote its last
paragraph to say the third question rides inside `WindowRecomputation` and why. D2 adds a fourth
that rides in the same place and leaves the paragraph reading:

> **The third question is not one of those, and rides inside [`WindowRecomputation`] instead.**

A reader who has just read the module doc's opening — which D2 *did* update to say there are four
questions — comes to the organisation section and finds only three accounted for. This is D1's
Major 1 recurring one step later.

**Fixed**: the paragraph now covers both, and names the one arm each has that the recomputation can
check — the floor sweep's arm at the shipped floor against the walk's count of absent windows, the
bin-scheme sweep's arm at the run's own configuration against the walk's histogram.

### Major 2 — "All three of those constants are provisional" contradicts the constants this commit rewrote

The module doc's fourth-question section ends with that sentence. In the same commit,
`src/ng/window_coverage/mod.rs` marks all three constants measured and kept, and spec §3.4 says
they stand. Two documents in one commit, opposite claims about the same three numbers.

This is D1's Major 2 in its exact shape: there, a value the step un-provisioned was still called
provisional 45 lines down the same file.

**Fixed**: the section now says the three were starting values until this measurement, that it kept
all three, and that re-running the sweep is how a later change to any of them is checked — which is
what the probe is for and what the present tense should be about.

### Major 3 — spec §3.6 still marks two settled constants soft, one of them D1's

The interface block in §3.6:

> `min_window_positions` (soft, §3.3), `depth_scale_windows` 10,000 (soft, §3.4).

§3.3 settled the first at D1 and §9 records it resolved; §3.4 and §9 settle the second at D2. The
one place a reader looks up what the accumulator is configured with says both are still soft.
Neither step's edits reached it.

**Fixed**: the line names the values, points at the two sections, and says both were the soft ones
and both were measured and kept.

### Major 4 — "780 times the margin the worst store needs" is a share ratio worn as a range width

The number is real: at the shipped range the worst sample-store puts 2.6 windows in every 10,000
past the top of the axis, against a fit that rejects at 2,000 in 10,000 — 780 times under. What
the sentence does with it is attach it to the range:

> the range has about 780 times the margin the worst store needs
> (`DEPTH_SCALE_WINDOWS`'s doc; and in the report, "a range with 780 times the headroom it needs",
> and "the reason is that the range has 780 times more headroom than the worst store needs")

A reader carrying that reads the range as shrinkable by something like that factor. The table
directly beneath it says otherwise: quartering the range, from 10 to 2.5, multiplies the overflow
share by 430. **The margin is in the guard, not in the range**, and the two are not interchangeable
because overflow grows far faster than linearly as the axis shortens.

This is D1's Major 3 — a ratio carried in the document's own unit for a different quantity — in the
same three documents.

**Fixed** in all three places: the claim is now stated in overflow shares (2.6 in 10,000 against a
guard at 2,000 in 10,000), and what the range itself has is shown by the neighbouring setting —
at a range of 5 the same prefix on the same store puts 35 in 10,000 over.

### Major 5 — "nine of the ten" for the prefix bias is eight, and the report's own table shows it

The claim appears in the report, in `DEPTH_SCALE_WINDOWS`'s doc, and in spec §3.4:

> the median fitted from the first 10,000 windows came out **below** the median over every window
> of the same store in nine of the ten

The report's table lists ten ratios; two are above 1.00 (1.02 and 1.08). I re-read the raw walk
output that the table was built from — `scheme_*.txt`, the ten `as-the-run-fits-it` and
`scale-sample=every-window` rows — and it is eight below and two above. "Nine of the ten" is
correct for a *different* claim in the same document, the count of stores that overflow nothing at
all, which is likely where it came from.

The average, 13%, is right as stated: it is the mean shortfall over all ten, the two above included.

**Fixed** in all three, with the two that came out above given their sizes (2% and 8%) rather than
left implicit.

### Major 6 — two `Option`s, an `unreachable!`, and an invariant that belongs in the type

```rust
struct OneBinSchemeAnswer {
    named: BinSchemeName,
    fitted: Option<FittedAxis>,
    silence: Option<String>,
}
```

Exactly one is `Some`; `as_a_row` enforces it with
`unreachable!("an answer is a fitted axis or a named silence")`. D1's review made the neighbouring
call the other way and gave the reason — a sweep and its answer are two things, split them so no
`Option` and no `expect` can be written. The same reasoning applies to one answer that is one of
two shapes, and the cost of the enum is one accessor.

**Fixed**: `WhatOneArmFitted { Axis(FittedAxis), NoHistogram(&'static str) }`, with an
`OneBinSchemeAnswer::axis() -> Option<FittedAxis>` for the callers that only want the axis. The
`unreachable!` is gone and the state it guarded cannot be constructed.

### Major 7 — `DEPTH_BINS`'s headline says the count was confirmed by measurement; its own last line says nothing measured it

> Depth bins … — **confirmed by measurement 2026-09-07**
> …
> Nothing measured forces this count either way; what it decides is the 80.2 kB.

Spec §3.4 makes the same pair, opening "All three measured 2026-09-07" over a bullet that says
"400 bins. Nothing measured forces the count". The body of both is right and honest — what the
measurement settled is the *scaling*, and the count is a memory knob priced later — and the
headline overstates it.

**Fixed**: the constant's headline says it is kept and no longer soft, and that what the
measurement settled is the scaling and not the count; §3.4's opening says "examined against
measurement" rather than "measured".

## Minor findings

### Minor 8 — `BinSchemeName` is a setting, not a name

The type carries the whole of what makes one arm differ from another — `configuration()` is built
from it — and its own doc concedes this ("Which of the two settings an arm varies, and to what").
The *name* is what `as_text()` returns. A reader hitting `named: BinSchemeName` cold expects a
label and finds the thing itself. **Fixed**: `BinSchemeSetting`, field `setting`, throughout.

`FittedAxis`, `OneBinScheme`, `the_bin_schemes_asked_about` and the three variants
(`AsTheRunFitsIt`, `ScaleSampleOf`, `RangeOf`) all say what their value is and are left alone. The
enum's shape — one variant carrying an `f64`, so no `Eq` — is right: the alternative is two fields
where only one is ever set, which is the defect of Major 6 rebuilt in the setting.
`OneBinSchemeAnswer` is the vaguest of the names (D1's equivalent is `WindowsUnderOneFloor`, which
says what it holds), but its content genuinely is *either* an axis *or* a silence, so "answer" is
not a topic word standing in for a quantity here. Left.

### Minor 9 — `THE_FITS_OVERFLOW_GUARD`: the name hides a possessive, and the reason beside it is false

The name reads as "the fits overflow guard"; the value is a share, and the name says neither that
nor whose. Worse, its doc:

> (`paralog::coverage_model::DEFAULT_MAX_OVERFLOW_FRACTION`, which is not `pub` to ng)

is untrue — it is a `pub const` in `pub mod coverage_model` in `pub mod paralog`, reachable from
this example as written. The report's deviation 2 says the opposite of the doc comment ("is `pub`
within production and production is frozen"), and *that* reasoning does not hold either: reading a
public constant is not an edit to frozen production. A copied 0.20 can drift from the number that
would actually judge the sample, which is the one property the doc claims for it.

**Fixed**: renamed `THE_OVERFLOW_SHARE_THE_FIT_ALLOWS` and bound to
`DEFAULT_MAX_OVERFLOW_FRACTION`, so the local name explains the number in this file's words and the
value cannot drift. `src/ng/window_coverage/mod.rs` already names the same constant in prose; that
reference is now an intra-doc link, so the toolchain checks the path exists. Deviation 2 rewritten.

### Minor 10 — `as_a_row()` does double duty, and the safety it advertises is not enforced

> also what the check against the walk compares, so that a field added to [`FittedAxis`] is
> compared as soon as it is printed

Two gaps. The function reached `axis.median_depth` and friends by field access, so **a field added
to `FittedAxis` is silently neither printed nor compared** — the doc's promise had nothing behind
it, and the house style (exhaustive destructuring, as this same file does for
`CoverageByGcHistogram` and `WindowRecomputation`) is exactly the device that would give it. And
what is compared is the *rendered* row, so the tolerance is the printed precision — `{:.6}` on the
bin width, `{:.4}` on the median. That is a defensible choice for comparing `f64`s, but it was
undocumented, so a later reader tightening the format would tighten the assertion without meaning
to, and one loosening it would loosen the assertion the same way.

**Fixed**: `as_a_row` destructures `FittedAxis` with every field named, and `close`'s doc says the
comparison is on the printed row and that the row's precision is its tolerance.

### Minor 11 — a wildcard arm over `SampleHistogram`, and `Debug` output stored as data

```rust
silence => Self { named, fitted: None, silence: Some(format!("{silence:?}")) },
```

Two things at once. A variant added to `SampleHistogram` is silently reported as a fourth silence
rather than a compile error — the one place in the new code that departs from the file's exhaustive
matching, and the departure is not argued in place. And the row's text is then whatever the variant
is *called in its declaration*, so renaming a variant changes the probe's output with nothing to
flag it.

**Fixed**: the three silences are matched by name and given stable labels
(`no-window-finalised`, `every-window-under-the-floor`, `median-depth-not-positive`), carried as
`&'static str`. A fourth variant is now a compile error.

### Minor 12 — `close`'s "does the arm exist" check and its comparison are the same lookup written twice

`assert!(answers.iter().any(|a| a.named == AsTheRunFitsIt))` followed by a `for` loop that
`continue`s past every arm that is not that one — a filter written as a loop, over a list where the
match is unique by construction — and then an `assert_eq!` whose custom message re-prints both rows
that `assert_eq!` prints anyway.

**Fixed**: one `find(...).expect(...)`, one `assert_eq!`, message kept (the `should_panic` test's
expected substring still matches, unchanged).

### Minor 13 — three doc comments that stopped describing what happened

- **`BinSchemeSweep`'s doc states a decision rule the step then did not follow**: "Varying
  `depth_scale_windows` says whether 10,000 windows is enough … *if the width barely moves when the
  whole store is used instead, it is*." The width moved — 34% at worst — and 10,000 was kept anyway,
  on the range's margin. A reader applying the stated rule to the shipped answer concludes the
  constant should have moved. **Fixed**: the doc says what each sweep shows, and adds the sentence
  the ruling actually turns on — that the two answers are read together, because what a short scale
  sample costs is paid in overflow, which is the range's quantity.
- **"the same device D1's floor sweep uses"**, in a Rust doc comment with no link and no gloss —
  a step label doing work a plain phrase does (rule 5; D1's review struck the same construction).
  Same in a code comment, "are D1's measurement", and in the report's "How it is measured".
  **Fixed** to "the floor sweep above" / "the floor measurement".
- **`BinSchemeSweep::new` carries no doc comment** — the constructor that builds ten accumulators,
  which is D1's Minor 5 on `FloorSweep::new` recurring. **Fixed**, pointing at
  `BinSchemeSetting::configuration`, which is where the `..` over the run's own configuration is
  written.
- **The report's "Tests added: Five"** over a table of six, with `21 passed (15 before)` two
  sections down. **Fixed** to six, with the counts.

### Minor 14 — the one per-store answer that is never totalled, unremarked

`WindowTotals::add` reads `WindowRecomputation` by field access rather than destructuring — D1's
review flagged this and left it, correctly, as older than the step. D2 adds two fields to that
struct, and one of them, `each_bin_scheme`, is a per-store answer that `add` does not sum. That is
the right call — an axis fitted from one sample's own depth has no total, and the cross-store number
a reader wants is the worst overflow, which is read off the rows — but the field's doc did not say
so, and the neighbouring `windows_under_each_floor` *is* summed. A reader comparing the two has to
guess whether the omission is a decision or a miss.

**Fixed**: the field's doc says it is printed per store and never summed, and why.

## Refactor safety, on the four changes asked about

| change | what happens | verdict |
|---|---|---|
| a field added to `CoverageByGcHistogram` | compile error in `OneBinSchemeAnswer::of` — the destructure names every field, with `window_bp: _` and `windows_under_the_floor: _` spelled and a comment saying why they are not part of the axis | sound as shipped |
| a variant added to `SampleHistogram` | **was** silently a fourth silence, rendered under its declaration name (Minor 11); now a compile error | fixed |
| a field added to `WindowCoverageConfig` | compile error in `the_configuration_a_run_uses`, which still names all six with no `..`; the `..run` inside `BinSchemeSetting::configuration` then carries the run's value to every arm, which is what an arm differing in one setting wants | sound as shipped |
| a field added to `WindowRecomputation` | compile error in `RecordMeasurement::observe`, whose `let Self { … }` names both new fields; but `WindowTotals::add` reads by field access, so a per-store count that ought to be totalled would be a zero row rather than an error (Minor 14, pre-existing) | documented |
| a field added to `FittedAxis` | **was** silently absent from both the row and the walk-tie check (Minor 10); now a compile error in `as_a_row` | fixed |
| an entry added to `the_bin_schemes_asked_about` | no compile error and none wanted — it returns a `Vec` and claims no length, unlike `CANDIDATE_FLOORS`'s `[u32; 25]`, which claims one. A duplicate entry (`RangeOf(10.0)`, say, or `ScaleSampleOf(10_000)`) would print a second identical row; harmless, and detectable in the output | left |

## Duplicated scaffolding between the two sweeps: the duplication is the honest choice

Both sweeps hold `arms: Vec<…>`, build one accumulator per candidate from
`the_configuration_a_run_uses()` with one field moved, loop the arms in `observe`, and consume
themselves in `close`. A generic `Sweep<A>` could hold the `new`/`observe` skeleton.

**It should not.** What would be shared is about ten lines, and the two `observe` bodies are not
the same ten: the floor sweep *counts every window as it pops* — that count is its whole
measurement — while the bin-scheme sweep drops them, because its measurement is the histogram the
arm ends with. Sharing would mean a trait whose one method is the interesting half of each sweep,
paid for by making the interesting half harder to read. The `close` bodies share nothing at all:
one asserts five properties about a distribution, the other three about an axis. Two plain structs,
each readable end to end, is right.

## Looked at and found sound

- **The split D1's review argued for was made without being asked.** `close(self, …)` consumes the
  arms; no arm's accumulator is an `Option`; there is no `empty()` constructor and no half-closed
  state to document. The answer is a plain `Vec<OneBinSchemeAnswer>` rather than a named type, and
  that is right here — `WindowsUnderEachFloor` earns its type by owning `print` *and* `add`, and
  this answer has no `add` to own.
- **Riding inside `WindowRecomputation` is right, and for a stronger reason than D1's.** The arm
  the sweep checks is checked against the walk's own `SampleHistogram`, which only exists there;
  a separate `RecordMeasurement` would have to be handed that histogram anyway, on top of the
  second reference accessor, contig check and depth-rule call D1's review priced. `close`'s ordering
  — the sweep closed after `self.histogram` is set, with a comment saying why — is the load-bearing
  detail and it is stated in place.
- **`ScaleSampleOf(u32::MAX)` as a match pattern** in `as_text`, giving the row
  `scale-sample=every-window` rather than `scale-sample=4294967295`. A const pattern, legal and
  clear, and the doc says what the sentinel means ("one no store here reaches").
- **The median is recovered, not re-derived.** `median = width × bins / range`, using **the arm's
  own** range and not the run's, with a comment saying so — the one place where reaching for the
  run's constant would have silently mis-scaled every `RangeOf` row.
- **The overflow sum reads the layout the way `coverage_model.rs` does** — one overflow column per
  GC row, `gc_bin * (depth_bins + 1) + depth_bins` — and the comment says which layout and why, so
  the number is the one the guard would see rather than a near relative.
- **The tests.** Six, and each states arithmetic a reader can check rather than a conclusion: the
  one-depth stream (median 8, width `8 × 10 / 400 = 0.2`, axis to 80) pins the formula in both
  directions; the two-contig stream pins the range against the guard at three settings; the
  10,300-then-19,700 stream pins the failure the scale sample trades against, and the short-stream
  test is named as its control and explains what it rules out. "A contig apiece, so a window never
  straddles two depths" is exactly the kind of setup note that stops a reader re-deriving the
  fixture. The two `should_panic` tests name the property, not the constructor.
- **The report's separation of "by construction" from "measured" succeeds where it is hardest.**
  "Fitting the width to the sample puts every median at bin 40 of 400 **by construction** — that is
  the formula, not a measurement — and the measurement is the 62-fold spread that makes it worth
  doing." The section on 400 bins then says plainly that nothing forces the count. "What this step
  does not measure" assigns the fit's accuracy to the filter's own branch and says so twice.
- **The arithmetic in the constants' docs.** 3.97 at production's 0.5× bin is bin 7 and 246.21 is
  bin 492 of 2,000, so three quarters of that axis sits above the deepest median; the fitted axis
  puts every median at bin 40 of 400. All check out. So do 1,972 of 7,666,421 = 2.6 in 10,000,
  0.20 / 0.000257 = 780, 0.00349 / 0.00026 = 13, and 2,000 / 1,103 = 1.8.
- **The mutual reference between `DEPTH_SCALE_WINDOWS` and `DEPTH_RANGE_IN_MEDIANS` (the hazard
  asked about) is accurate in both directions, once Major 4's unit is fixed.** The scale doc's "a
  prefix a third shallow makes the axis a third short" is right — the axis top is the fitted median
  times the range. The range doc's "the whole of that worst store's overflow comes from the prefix,
  and vanishes when the width is fitted from every window" is right and is *measured*: that store's
  `scale-sample=every-window` arm overflows 0 where its shipped arm overflows 1,972. Neither
  overstates its own contribution. What overstated was the size of the margin, and only in the
  scale doc's direction.

## What I changed, so it can be lifted

**`examples/ng_window_coverage_probe.rs`**

1. Module doc: "One walk, several measurements" now covers the third and fourth questions, and
   names the arm each has that the recomputation checks. "The fourth question" no longer calls the
   three constants provisional.
2. `BinSchemeName` → `BinSchemeSetting`; field `named` → `setting`, everywhere including the tests.
3. `OneBinSchemeAnswer`'s two `Option`s → `fitted: WhatOneArmFitted`, a two-variant enum; new
   `axis()` accessor; `unreachable!` deleted.
4. `OneBinSchemeAnswer::of` matches `SampleHistogram`'s three silences by name and labels them
   `&'static str`, instead of a wildcard arm rendering `Debug`.
5. `as_a_row` destructures `FittedAxis` with every field named; `close`'s doc says the comparison is
   on the printed row and that the row's precision is the tolerance.
6. `close`: `find(...).expect(...)` plus one `assert_eq!` in place of an `any` assert, a
   `continue`-filtered loop and a message that re-printed both rows. Panic text unchanged in the
   substring the `should_panic` test matches.
7. `THE_FITS_OVERFLOW_GUARD` → `THE_OVERFLOW_SHARE_THE_FIT_ALLOWS`, bound to
   `pop_var_caller::paralog::coverage_model::DEFAULT_MAX_OVERFLOW_FRACTION`; the false "not `pub`
   to ng" replaced by why the binding is not an edit to frozen production.
8. `BinSchemeSweep::new` given a doc comment; `BinSchemeSweep`'s own doc no longer states a decision
   rule the step did not follow, and says the two answers are read together; the "D1's floor sweep"
   and "D1's measurement" shorthands replaced.
9. `each_bin_scheme`'s field doc says it is printed per store and never summed, and why.
10. Tests: new `the_axis(answers, setting)` helper replacing five `.fitted.expect(…)` sites.

**`src/ng/window_coverage/mod.rs`** (doc comments only — no library behaviour is touched)

- `DEPTH_BINS`: headline says kept and no longer soft, and that the measurement settled the scaling
  rather than the count.
- `DEPTH_SCALE_WINDOWS`: "nine of the ten" → eight, with the two above given their sizes; the
  "780 times the margin" sentence restated in overflow shares, with the range-of-5 figure beside it.
- `DEPTH_RANGE_IN_MEDIANS`: "the safe band is 5 to 40" replaced by what was measured — no setting
  tried made the fit reject a sample, ten was kept for the size of its margin, 40 is where the
  settings stop; "headroom" → "margin under the guard", with both ends of the thirteenfold given.
  The guard is now an intra-doc link to `crate::paralog::coverage_model::DEFAULT_MAX_OVERFLOW_FRACTION`.

**`doc/devel/ng/spec/window_coverage.md`**

- §3.4: the leftover future-tense "The plan measures the overflow fraction … (step D2)" removed —
  the paragraph beneath it is that measurement; "All three measured" → "All three examined against
  measurement"; "nine of the ten" → eight; "the safe band is 5 to 40" replaced as above;
  "the range's headroom" → "the range's margin under the guard".
- §3.6: the interface block no longer marks `min_window_positions` and `depth_scale_windows` soft.
- §9: the prefix bias given its worst case (34%) beside its average (13%), since the worst is what
  the margin has to absorb; "the range has that much margin" rephrased.

**`doc/devel/reports/implementations/ng_window_coverage_d2_2026-09-07.md`**

- "The answer": the 780 stated as an overflow share against the guard, not as range headroom.
- The range section: "no setting in that table made the fit reject a sample, 2.5 included" replaces
  "the safe band is 5 to 40"; the three stores that overflow at 2.5 are all named (1,103, 614 and
  191 in 10,000) with their fitted medians, and the depth-ordering claim they do not support is
  withdrawn.
- The scale-sample section: "nine of the ten" → eight, with the two above given their sizes; the
  headroom sentence restated in overflow shares.
- "How it is measured": the two assertion bullets now describe the checks as they are written; the
  "D1's floor sweep" shorthand replaced.
- Deviation 2 rewritten for the binding rather than the copy.
- "Tests added: Five" → six, with the 21-against-15 counts.

## Validation

Run in the container, on this worktree, on the tree as I am leaving it:

- `./scripts/dev.sh cargo test --all-features --example ng_window_coverage_probe` —
  **`test result: ok. 21 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out`**, the same 21 the
  step ships.
- `./scripts/dev.sh cargo check --lib --tests --all-features` — **clean**: "Finished dev profile
  (unoptimized + debuginfo) target(s) in 13.33s", no warnings.
- `./scripts/dev.sh rustfmt --check --edition 2024 examples/ng_window_coverage_probe.rs
  src/ng/window_coverage/mod.rs` — **0 hunks on both files**, exit 0.
- `./scripts/dev.sh cargo clippy --all-features --example ng_window_coverage_probe` — 3 warnings,
  all `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch; none in the
  example.
- **Not re-run: the five real-store walks.** No `.psp` exists in this worktree. The renames and the
  enum do not touch any printed value, and the one output change — the three silences now labelled
  in this file's own words rather than by `SampleHistogram`'s variant names — cannot appear on a
  store that fits an axis, which every one of the five does. The numeric claims I corrected were
  re-derived from the walk outputs the step recorded (`tmp/d2/scheme_*.txt` in the branch's own
  worktree), not from a re-run.
