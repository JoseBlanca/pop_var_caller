# window coverage — D3 review (naming, structure, prose): the answer is the wrong quantity, and the two tests test arithmetic

**Date:** 2026-09-07
**Reviewed:** commit `7fae9525 wip: D3 for review`, diffed against `bdf80ada`
**Remit:** naming, idiom, module structure, refactor safety, and the prose the step ships.
Correctness, error handling and the numbers are a separate reviewer's.
**Implementation report:** [ng_window_coverage_d3_2026-09-07.md](../implementations/ng_window_coverage_d3_2026-09-07.md)
**The two steps before it, and their design reviews:**
[the floor](ng_window_coverage_d1_design_2026-09-07.md),
[the bin scheme](ng_window_coverage_d2_design_2026-09-07.md)
**Fixes:** applied in this worktree, listed at the end.

## Verdict

**The question this review was set — did D3 repeat the pattern the previous step's review found,
where the code rulings were followed and the prose rulings re-broken — has a cleaner answer than
D2's, and a worse one.** There is almost no code here to follow a ruling with: two tests. The
prose broke three of the four defects the earlier reviews named, and it broke a fourth that is
new and is the most serious thing in the step.

**The report's headline compares the wrong quantity against the budget.** It answers "100.2 kB a
sample, a fifth of the per-open-sample budget", and its own table shows a 262 kB term live at the
same time as the other three — the held-back windows, which sit in memory beside a histogram that
is allocated whole at the pass's start, not when the axis is fitted. A budget is about the peak.
The peak is **362 kB a sample, just under three quarters of the 500 kB**, and no sentence in the
step says so. The table's own "when it is live" column had the histogram wrong, which is what
allowed the four terms to look as though they were never live together.

**"Measured" is used for a number that was not measured, in the one document that outlives the
report.** Spec §4 said "**301 MB, measured 2026-09-07**" and §5 "**100.2 kB a sample in all,
measured**". The report is careful and says the opposite — the figure is priced by construction
because the run cannot see it. This is exactly the recurring defect: a claim about the design
presented as a measurement.

**Both tests assert arithmetic the code has no part in.** `501 × 16` and `250 × 48` and
`16_384 × 16` are three multiplications; only one assertion in either test — `counts.len()` —
touched an object. And the first multiplication is wrong: a `VecDeque` allocates its capacity, so
the sliding buffer is **8,192 bytes and not 8,016**, which the report's own careful treatment of
the held-back `Vec`'s doubling should have caught. Both tests now read the real capacities off a
fed accumulator, and the figure is corrected everywhere it appears.

**One claim is an inference dressed as a finding, in both the report and the spec**: that the
run-charged term "is the reference the cache now reads once per cover". What was measured is that
it grows with ground. The identification does not survive the report's own arithmetic — 200 kb of
ground is 200 kB of reference bases against a term of 33.2 MB, some 166 bytes for every base
walked.

**13 findings, all applied. 1 Blocker, 5 Major, 7 Minor.**

## Findings

### Blocker 1 — the answer is compared against the budget in the wrong quantity

The headline is "**100.2 kB a sample, a fifth of the per-open-sample budget**". The table beneath
it lists four terms and draws the total line after three, leaving the fourth — the 262 kB the
depth axis is fitted from — below the line as though it were not charged.

It is charged. Three things make all four live at once:

- `counts: vec![0; config.cell_count()]` runs in `WindowCoverageAccumulator::new`, so the whole
  80.2 kB histogram exists from the sample's first record, not from the fit. The table said
  "from the depth axis being fitted to the end of the pass", which is the opposite of what the
  constructor does;
- the sliding buffer is live from the first record too;
- the merge holds the look-ahead's extra records throughout.

So a sample's peak is `80,200 + 8,192 + 12,000 + 262,144 = 362,536` — **362 kB, 73% of §7.2's
500 kB an open sample**, not 20%. On a human reference after the psp path's contig-list fix the
open sample is 123 kB, so the peak is 486 kB against a 500 kB budget: it fits, with 14 kB to
spare. That is a materially different message from "223 kB", and it is the number a reader
sizing a run needs.

**Fixed** in the report's answer, in spec §4 and §5, and in the follow-up bullet, which now says
the transient is what puts the peak at 362 kB and is the first thing to cut if the budget binds.
The table's "when it is live" column for the histogram is corrected to "allocated whole at the
pass's start".

### Major 2 — the spec calls a constructed figure "measured", twice

Spec §4: "**301 MB, measured 2026-09-07**". Spec §5: "**100.2 kB a sample in all, measured**".
The report says plainly that the per-sample cost is *priced by construction* and that a whole-run
measurement cannot resolve it. A spec outlives every report that feeds it; leaving "measured"
there means the next reader who needs to know how good the number is gets the wrong answer, and
will not find the qualification unless they follow the link.

This is D1's Major 3 and D2's Major 4 in a new shape — there a ratio in one unit was carried in
the language of another; here evidence of one kind is carried in the language of another.

**Fixed**: "Added up from the shipped constants 2026-09-07 with each per-unit size pinned by a
library test; a whole-run measurement beside it could only bound the per-sample cost at 1.8 MB."

### Major 3 — "every byte is asserted by a test" is not true of two of the three terms

The claim appears in the report's answer and again in spec §4 ("Every one of those bytes is
asserted by a test rather than derived in prose"). Of the three charged terms:

- the histogram's 80,200 **is** pinned — the test reads `counts.len()`;
- the sliding buffer's 8,016 was `(window_bp + 1) × size_of::<CoveredPosition>()`, a
  multiplication, and it is the wrong number (finding 4);
- the held records' 12,000 was `250 × 48`, a multiplication whose second factor is a struct size
  and whose first is an assumption — *one record a base* — that no test can pin and that the
  step's own doc comment admits is coverage-dependent.

**Fixed**: the claim is now "each term's per-unit size is pinned by a library test", with the
record count named as the assumption it is, in both documents.

### Major 4 — the sliding buffer is 8,192 bytes, not 8,016

`positions` is a `VecDeque`, and a `VecDeque` allocates its capacity. Fed a dense stream at the
shipped `WINDOW_BP`, its length reaches 501 and **its capacity reaches 512** — 8,192 bytes. The
report applies exactly this correction to the held-back `Vec` ("10,000 pairs of `f64` is 160 kB
of live data, and the list grows by doubling, so its capacity reaches 16,384") and does not apply
it here.

The correction moves the whole-pass total from 100,216 to 100,392 — 100.2 kB to 100.4 kB — and
leaves every derived figure (a fifth of the budget, 301 MB at three thousand) standing.

**Fixed** in the test, which now measures the capacity instead of multiplying, and in every
document that quotes it.

### Major 5 — "it is the reference the cache reads once per cover" is an inference stated as a finding

The report's section is titled "The one thing the measurement can see, **and where it comes
from**", and ends:

> **So the term is not per sample** … it is the reference the cache now reads once per cover
> ([C1](…)) and the ground each cover spans

What the two ground sizes establish is that the term follows the ground and not the cohort. They
say nothing about which allocation it is. And the identification does not survive the report's
own numbers: at 200 kb of ground the term is 33.2 MB, which is **166 bytes for every base
walked**, where a reference slice is one byte a base.

Spec §5 carried the same claim, and carried it more confidently than the report — "the
*reference* the cache reads once per cover does grow with the ground a run walks, and D3
measured that".

**Fixed** in both. The report now separates what was measured from what was guessed, names the
heap profile that would settle it, and adds an argument the report had in hand and did not use:
at 200 kb the two arms differ by 33.2 MB at one sample and 91.1 MB at 63, where a per-sample term
would have put about two gigabytes between them — which is what actually rules out the per-sample
candidates (the look-ahead's held records, and the arena each sample's psp source never
releases). Spec §5 now says the allocation was not measured and gives the size argument in one
clause.

### Major 6 — a spec sentence the step's own measurement contradicts, left standing

Spec §3.3, on what the look-ahead costs:

> about 250 at three reads a position, **some tens of kilobytes at a thousand samples**

The step measures that term at 12 kB a sample, which is 12 MB at a thousand — three orders of
magnitude from "some tens of kilobytes", and the step walked past it. The same sentence is quoted
in `ObservationCache::cover`'s doc comment, which cites §3.3 for it.

This is the "a section that stopped describing what the file now does" defect both earlier
reviews found, in its cross-document form: the step wrote the new number in §4 and left the old
one in §3.3.

**Fixed** in both places, with the psp-mode qualification and a pointer to the test.

### Minor 7 — the two tests name a document, where the house convention names the property

`one_samples_measurement_costs_what_the_memory_report_says_it_does` and
`a_held_record_costs_what_the_memory_report_says_it_does`.

The instinct behind them is right — the report is the thing that goes stale, and the test exists
to catch that. But the codebase already has this test and names it the other way:
`the_observation_stays_a_handful_of_scalars_wide` in `ng/calling/likelihood/mod.rs`, which pins
`size_of::<GenericObservation>() == 32` for the same reason. A name that says what is pinned
survives the report being renamed, superseded or merged into another; a name that says "what the
report says" tells a reader who hits a failure nothing about what broke.

**The place for the document is the failure message, not the name** — and the messages were where
this step was thinnest: eleven of the thirteen assertions across the two tests carried none at
all, so a field added to `LocusSummary` produced `assertion left == right failed: 40 != 32` and
no hint that three paragraphs in two documents needed rewriting.

**Fixed**: renamed to `one_samples_accumulator_is_88_kb_for_the_pass_and_262_kb_more_until_the_axis_is_fitted`
and `a_held_record_in_psp_mode_is_a_summary_and_a_range_at_48_bytes`, and every assertion given a
message naming the figure it guards and the documents that quote it.

### Minor 8 — `250 * 48 == 12_000` is a test of multiplication, and what would make it a test

Taken literally the assertion cannot fail unless one of the two lines above it fails first:
`half_a_window` is asserted to be 250 two lines up, and `a_held_record` is asserted to be 48 three
lines up. It is a tripwire on a number three documents quote, which is worth having — but it is
not a test of the cache, and the step presented it as one ("asserted here rather than left as
arithmetic in a report").

**What would make it a real test is already in the file, immediately below it**:
`a_cover_reads_half_a_window_past_its_region_with_nothing_chaining_there` pins that the cover
really does draw `WINDOW_BP / 2` past its region and no further. That is the property the 12 kB
rests on, and it is the reason the arithmetic is sound rather than a coincidence.

**Fixed**: the product is kept as the tripwire on the documents' figure, its message says so and
shows both factors, and the doc comment points at the neighbouring test as the thing that makes
the count meaningful rather than claiming it here.

### Minor 9 — `size_of::<CoveredPosition>() == 16`: pin it, but say what a failure means

Pinning a layout the compiler chose is brittle in principle. In practice it is right here, and
the repository has already made the call: `GenericObservation`'s test pins 32 **and** its
alignment **and** `!needs_drop`, with a comment explaining the padding. Three reasons to keep it:

- the whole point of the step is that these figures do not drift silently, and a field addition
  is exactly the drift;
- every type involved is `repr(Rust)` over `u64`/`u32`/`f32`/`f64`, whose sizes are fixed on
  every target this project builds for; the risk is a compiler layout change, not a target
  change;
- a failure is cheap to read *if the message says what to do*, which is finding 7.

**Kept, with messages.** I did not add `align_of` or `needs_drop` assertions: nothing in the
memory argument turns on either.

### Minor 10 — `assert!(histogram + sliding_buffer < 100_000, "… the report calls them about 88 kB")`

A loose bound with the precise number in its message is the worst of both. It passes while the
histogram grows by 11 kB and the report goes quietly wrong, which is the exact failure the step
exists to prevent; and the message asserts in prose what the code declined to assert in code.

**Fixed**: an equality on the sum, `88_392`, with a message that names what the figure is and
what the report adds to it.

### Minor 11 — the test doc's own arithmetic is wrong, and its list of terms disagrees with the report's

Three things in the doc comment on the accumulator test:

- "A window spans `window_bp` bases either side of its centre inclusive" — it spans
  `window_bp / 2` either side, 250, which is what makes the buffer 501 entries and not 1,001. The
  arithmetic below it uses the right number; the sentence explaining it does not.
- "at most `window_bp + 1` of them … it is a ceiling because a window is only ever that wide",
  two lines above an inline comment saying "a repeated position is one entry per observation" —
  and the type's own doc says a stream repeating a position "would buffer without bound". Three
  statements, two of them incompatible. The ceiling holds for the one-record-a-position stream
  the walk produces, which is the condition the type states and the test's doc dropped.
- "**The four terms, and each is charged to the per-open-sample budget**" listing the ready deque
  third, then "the first two are live for the whole pass; the last is live only until the
  sample's ten-thousandth window" — leaving the third unaccounted, and listing a set of four that
  is not the report's set of three.

**Fixed**: the geometry corrected, the ceiling given its condition, and the doc reorganised around
the three terms this type owns, with a paragraph saying why the ready deque is priced with the
merge's held records instead.

### Minor 12 — the held-record figure is psp mode's, and nothing said so

`size_of::<LocusSummary>() + size_of::<Range<usize>>()` is what a psp-mode sample holds per
record. Direct mode holds `SampleLocusObservations` beside them — **120 bytes inline plus its
evidence on the heap** — so its look-ahead term is several times larger. The test's doc gestured
at this and garbled it ("direct mode holds the built record too, which the look-ahead does not
change the *size* of, only the count"), and neither the report's table nor spec §4 said which
mode the 12 kB belongs to. The whole-run arm was psp mode throughout, so the scope is real and
easy to state.

**Fixed** in the test name, the test doc, the report's answer, its deviations list, and spec §4.

### Minor 13 — the report's ± is a standard error worn as an interval, and the follow-up bullets

- "a difference of **0.82 ± 0.49** MB a sample, an interval nineteen times wider than the whole
  effect and consistent with zero". As printed, the interval is 0.98 MB wide, ten times the
  effect, and excludes zero. Both of the report's claims are right only if 0.49 is a **standard
  error**: `2 × 1.96 × 0.49 = 1.92`, which is nineteen times 0.10, and `0.82 < 1.96 × 0.49`, so
  the interval spans zero. This is the same defect as D1's "moves by 6 in 100" and D2's "780
  times the margin": a quantity carried in the language of a different one. **Fixed** by naming
  the standard error and showing the interval it implies, −0.1 to +1.8, in both the slope and
  intercept places. *The other reviewer should confirm 0.49 is the standard error; both of the
  report's own claims require it.*
- "the instrument's job here is to show that nothing large appeared, and nothing did", and "that
  is the honest result rather than a failure of the run". Delete both clauses and nothing
  checkable is lost — rule 6. Worse, they hand the reader a defence in place of the result. The
  run does give a result: **a ceiling of 1.8 MB a sample, eighteen times looser than the figure
  the construction gives.** **Fixed** to say that, in the answer and in the whole-run section.
- "**the 39.5 MB a sample is the number that should worry someone**" — tells the reader how to
  react in place of the fact, which is in the next clause (79 times the per-open-psp budget).
  Placement at the head of the follow-ups is right: it is the largest number in the step, it is
  not this plan's, and a reader who stops after the answer should still meet it. What it lacked
  was the thing a follow-up exists to give — **what would be needed to chase it**. **Fixed**: the
  reaction clause struck, and the bullet now names the instrument (a heap profile by allocation
  site, at two or three cohort sizes on this slice) and says why the one used here cannot do it.
- "it is **a fifth of the budget** for the first ten thousand windows" — 262 kB is over half of
  500 kB. **Fixed.**
- "**D2 settled** that its 400 bins are not forced by any measurement" — a bare step label doing
  work a phrase does, which both earlier reviews struck. **Fixed** to "the bin-scheme measurement
  one step earlier", with the link.

## Should the spec carry these figures at all?

**Yes, and the step put them in the right sections** — but §4 was carrying more of the report's
working than it should.

The convention is visible in how the two previous steps closed their own questions. §3.3 carries
the floor's measured default, one table's worth of consequence, and a link to the report that has
the distribution; §3.4 carries the three bin-scheme constants, what each was measured against, and
a link. In both cases the spec holds **the requirement and the settled answer**, and the report
holds the working. Both also add a line to §9's resolved-decisions list.

§4 and §5 are the right homes here because §1.1's sixth goal is a requirement — "memory per
sample stated and bounded, charged against the calling stage's per-sample budget" — and §4 and §5
are where a spec answers it. What did not belong was the report's provenance sentence inside the
spec bullet. **Trimmed**: §4 now carries the two per-sample figures, the breakdown, the comparison
against the budget, the reconciliation with §3.4's histogram-only 240 MB, and one clause on how
the figure was arrived at.

**Not added: a §9 entry.** §9 records *decisions*, and this step decided nothing — it priced a
design already settled. D1's and D2's entries close questions §9 had open; there is no open
question here to close. Worth the author confirming, since it is the one place D3 departs from
the shape of the two steps before it.

## Refactor safety

The question asked was what happens when someone adds a field to one of the four types.

| change | before this review | after |
|---|---|---|
| a field on `CoveredPosition` | `assert_eq!` fails, `16 != 24`, no message | fails naming the 8.2 kB and the two spec sections and the report that quote it |
| a field on `WindowMeans` | fails, `16 != 24`, no message | fails naming the 262 kB and where it is quoted |
| a field on `LocusSummary` | fails, `32 != 40`, no message | fails naming the 48-byte held record and the documents |
| a field on `WindowCoverage` or `GenomePosition` | fails on the tuple's 24, no message | fails naming spec §5's 24-byte deque entry |
| a field on `WindowCoverageConfig` | compile error in the test's literal (no `..`) — the same guarantee the probe's `the_configuration_a_run_uses` has, but undocumented | unchanged, and the comment now says the re-spelling is deliberate |
| `WINDOW_BP` changed | both tests fail on derived figures | unchanged, with messages |
| a growth-strategy change in `Vec`/`VecDeque` | invisible: the figures were multiplied out | fails, because the capacities are now read off a fed accumulator |

The two tests are in the right modules: the accumulator's terms with the accumulator, the held
record with the cache that holds it. Neither reaches across.

## Looked at and found sound

- **The step touches nothing frozen.** `src/sample_summary/`, `src/paralog/` and
  `src/var_calling/` are untouched, and the two tests add no behaviour — the library gains
  assertions and nothing else.
- **Leaving the plan's D3 box unticked** is this branch's practice at a review commit, not an
  oversight: `f3228dff`, the D2 review commit, has D2 unticked and `bdf80ada`, the commit that
  landed it, has it ticked.
- **The measurement's design, as reported.** Nine cohort sizes rather than the plan's three, with
  the reason given (three points cannot separate a slope from an intercept at this scatter); the
  same 63 stores in the same order at every size, so a step from *N* to *M* changes one thing;
  the same psps for both binaries, with equal VCF record counts at all three sizes as the check
  that the comparison is a comparison; `VmHWM` because the container has no `time(1)`; and the
  harness naming its binary rather than taking the newer of two target directories, with the
  reason spelled out. Each of those is a decision a reader can check and disagree with.
- **The ground experiment is the right instrument for the question it was given.** Changing the
  ground at one sample separates a per-run term from a per-sample one in a way that adding
  samples cannot, and the repeat spreads are quoted beside it (5.3–17.8 against 18.6–19.7)
  rather than left implicit.
- **The deviations list is honest about the one that matters** — that the plan asked for a
  per-sample slope from the run and the run cannot resolve one — and says what was reported
  instead. I added a fifth deviation it was missing: the measurement is tomato-only and psp-only,
  where the plan's verification summary asks Milestone D's measurements for both benchmarks.
- **`a_cover_reads_half_a_window_past_its_region_with_nothing_chaining_there`**, the test the new
  one sits beside, is the model the new tests should have followed and now point at: it takes its
  number from the constant, checks the property in both directions (the base at
  `50 + half_a_window` is readable, the one past it is not), and its message says what a failure
  would mean for a run.

## What I changed

**`src/ng/window_coverage/accumulator.rs`** (0 rustfmt hunks, unchanged)

- Test renamed to
  `one_samples_accumulator_is_88_kb_for_the_pass_and_262_kb_more_until_the_axis_is_fitted`.
- The sliding buffer and the held-back transient are now **measured on a fed accumulator** —
  2,000 positions for the buffer, 12,000 for the transient — reading `capacity()` rather than
  multiplying. The buffer's figure is corrected to 8,192; the transient's 262,144 is confirmed as
  the capacity the list actually reaches, and a `matches!` assertion pins that the axis is fitted
  by the ten-thousandth window, which is what ends it.
- `assert!(… < 100_000)` replaced by `assert_eq!(histogram + sliding_buffer, 88_392, …)`.
- Every assertion given a message naming its figure and the documents that quote it.
- Doc comment: the window geometry corrected, the buffer's ceiling given its condition, the four
  terms reorganised as the three this type owns plus a paragraph on why the ready deque is priced
  elsewhere, and a note that the two growing collections are read rather than multiplied.

**`src/ng/run/cohort_merge/observation_cache.rs`** (4 rustfmt hunks, all predating this branch, none in the touched region)

- Test renamed to `a_held_record_in_psp_mode_is_a_summary_and_a_range_at_48_bytes`; messages added
  to each assertion; the redundant `assert_eq!(a_held_record, 48)` and `assert_eq!(half_a_window,
  250)` folded into the messages of the assertions that use them.
- Doc comment: names the mode, gives direct mode's extra cost (120 bytes inline plus heap), says
  the count is the assumption and the size is what is pinned, and points at the neighbouring
  look-ahead test as what makes the count meaningful.
- `cover`'s doc: "some tens of kilobytes across a thousand samples" replaced by the measured
  12 kB a sample and 12 MB at a thousand, with the mode and a pointer to the test.

**`doc/devel/ng/spec/window_coverage.md`**

- §3.3: the look-ahead's cost restated at the measured figure.
- §3.4: "240 MB at three thousand" marked as the histogram alone, pointing at §4 for the total.
- §4: the three-thousand-sample bullet rewritten — both per-sample figures with their breakdown,
  both compared against §7.2, the reconciliation with §3.4, "measured" replaced by how the figure
  was arrived at, and the report's working left to the report.
- §5: the memory bullet rewritten — the terms that make the 100.4 kB and the 362 kB, the ready
  deque explained and explicitly excluded, the run-charged term stated as measured-in-shape and
  unattributed-in-fact.

**`doc/devel/reports/implementations/ng_window_coverage_d3_2026-09-07.md`**

- "The answer": both figures, the corrected table (8,192; the histogram's true lifetime; a peak
  row), the budget comparison at both, and "what is asserted and what is assumed" in place of
  "every byte above is asserted by a test".
- The whole-run section: the standard error named as one, the interval shown, and the ceiling
  stated in place of the two clauses defending the instrument.
- The run-charged term: renamed section, the cohort-table argument added, the per-sample
  candidates ruled out by it, and the reference named as an unconfirmed candidate with the size
  argument against it.
- A **Changes made** section, which the report was missing and both previous steps' reports have,
  listing the spec edits.
- Deviation 5 added (tomato only, psp mode only). Deviation 1 restated.
- Follow-ups: the reaction clause struck and the instrument named for the 39.5 MB; "a fifth of
  the budget" corrected to over half; the step label replaced.
- Test names and what they pin updated throughout.

## Validation

In the container, on this worktree, on the tree as I am leaving it:

- `./scripts/dev.sh cargo test --lib --all-features` — **`6380 passed; 0 failed; 15 ignored`**.
- `./scripts/dev.sh cargo check --lib --tests --all-features` — clean, no warnings.
- `./scripts/dev.sh rustfmt --check --edition 2024 src/ng/window_coverage/accumulator.rs
  src/ng/run/cohort_merge/observation_cache.rs` — **accumulator.rs 0 hunks**, unchanged;
  **observation_cache.rs 4 hunks**, all at lines 260, 971, 1061 and 1762, none in this step's or
  this review's regions, and the same 4 the branch already had.
- `./scripts/dev.sh cargo clippy --lib --tests --all-features` — the same warnings as before these
  edits: three `needless_lifetimes` in `cohort_merge`, three byte-string suggestions, and three
  `1 * n` index expressions in accumulator tests that predate this branch. None in the new code.
- **Not re-run: the whole-run memory measurement.** No psp exists in this worktree, and nothing I
  changed can move a byte of it — the corrections to the constructed figures are arithmetic over
  struct layouts, which the two tests now check on every run.
