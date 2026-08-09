# The runs model's third refusal: when the starting points score alike and answer differently

**Date:** 2026-08-09. **Not a plan step** — it adopts a measurement
(`doc/devel/ng/research/inbreeding_resolution_2026-08-09.md`) that was made after Milestone E
shipped, and it lands before Milestone F begins because it changes what `fit_inbreeding`
returns.

## What the measurement said, and what it turned out to have missed

Research note §4 found that the across-start spread of `F` — already computed, already in
`RunsModelFit::starts_tried`, and never read — separates a genome whose runs are real from one
where the chain is reading sampling noise. Its recommendation was a threshold near 0.05 on
that spread.

**Two things had to be measured before that could be adopted, and the second changed it.**

### 1. The shape the note called unconfirmed (note §6.1)

The note's with-runs column was *exactly* 0.0000 at every shape it tried, all of them 30% of
the genome in runs of 3 Mb. `examples/ng_inbreeding_resolution.rs` gained a §4 that asks the
same question of genomes whose runs are shorter or rarer — 300 kb runs, and 10% and 5%
coverage — crossed with the four evidence levels, five seeds each.

**All 100 answered, none refused, worst recovery error 0.0100, and every one has a spread of
0.0000.** Over §3 and §4 together — **160 fits at twenty shapes** — the largest spread is
0.0000. The suspicion does not hold up: the starts agree wherever the runs are real.

### 2. The shape nobody had named, and it is the one that matters (note §6.2)

**Every genome in the harness is drawn without a floor of false heterozygotes.** Milestone
E3's own fixtures are not: collapsed paralogs and mismapping lift both states together, and
spec §6.2 requires this estimator to survive a floor of five times the real heterozygote
rate.

Applied as the note recommended, the criterion refused three of `runs.rs`'s tests — two of
them fits that recover their genome:

| fixture | fitted `F` against realised | raw spread over nine starts |
|---|---|---:|
| 30% runs, floor 1× | 0.3150 against 0.3153 | 0.0000 |
| 30% runs, floor 3× | recovers, within 0.05 | 0.3250 |
| 30% runs, floor 5× | 0.3157 against 0.3122 | 0.3157 |

**Dumping the nine starts is what explains it, and the explanation is not "the fit is
marginal".** At the 5× floor, six starts land on `F` = 0.3157 and three collapse to 0.0000 —
and the three that collapsed score **1,473 nats worse**. They are not disagreeing with the
answer; they have been rejected by a likelihood ratio of e^1473.

The failure looks nothing like that. On a genome with no runs at five heterozygotes a window
the nine starts return `F` = 0.0010, 0.8497, 0.6015, 0.5722, 0.0003, 0.0051 and 0.0159 — and
every one is within **0.91 nats** of the best, odds of 2.5 to 1. That is research note §3.1's
proof showing through: at coincident states the likelihood is exactly flat in `F`, so every
answer scores the same.

**So the score separates the two populations where the spread does not** — 0.91 nats against
1,473, a factor of 1,600 — while the spreads themselves overlap, 0.30 to 0.33 legitimate
against 0.0004 to 0.998 not.

## What was built

**`MAX_TIED_START_LOG_LIKELIHOOD_GAP = 10.0` nats.** Which starts count as tied with the
best. Eleven times above the measured tie, 147 times below the measured rejection; an odds
ratio of about 22,000 to one. **An absolute number and not a fraction of the total**, because
a likelihood *ratio* is what compares two fits and it does not grow with the genome — the two
totals behind those gaps differ eightfold and the gaps do not track them.

**`MAX_IDENTIFIED_START_SPREAD = 0.05`.** How far the tied starts may disagree about `F`.
Above every measured legitimate fit by three orders of magnitude and 2.3× below the nearest
measured failure.

**`tied_starts` and `spread_across_tied_starts`,** the second exposed as
`RunsModelFit::spread_across_tied_starts()`. A **method rather than a field**, so the number a
consumer reads and the number the fit was accepted on cannot drift apart; the harness reads it
through the same method for the same reason. `tied_starts` returns a **prefix**, because
`starts_tried` is documented best-first — so a stray start further down cannot be picked up.
A `NaN` likelihood counts as tied, which refuses rather than admits.

**`ParameterEstimationError::InbreedingStartsDisagree`,** carrying the sample, how many starts
tied, how many were tried, the spread and the threshold. **A fifth variant rather than a reuse
of `InbreedingStatesNotSeparated`**: that one says the search never found a second state, and
this one says it found a different one from every start with nothing to choose between them.
A message naming the wrong failure sends the next reader hunting a mechanism that did not
occur.

**The check is third and last in `fit_inbreeding`,** after both separation checks. Deliberate:
a collapsed search also has starts that disagree, and it should be reported as the collapsed
search it is. It also means that *reaching* this error proves the other two passed — which is
why the test asserts the variant rather than merely asserting a refusal.

## What it is proven against

- **`the_spread_is_the_range_over_the_starts_that_scored_alike`** — four written-down outcomes
  in best-score order. The widest pair is not the top two (best-and-second-best would report
  0.02 for 0.30); the fourth start sits at the measured 1,473-nat gap and is excluded
  (including it would report 0.32 and refuse the fit); the three kept sit 0.00, 0.02 and 0.30
  behind, so a rule keeping only exactly-equal starts would report 0.00.
- **`starts_that_land_in_different_places_are_refused_rather_than_averaged`** — a genome with
  no runs at five heterozygotes a window, the widest of eight seeds swept at that shape.
  **Seven of the eight are refused**, spreads 0.1147 to 0.8494, and every one of the seven
  passed both other checks.
- **`a_real_run_at_the_same_evidence_level_is_still_answered_and_the_starts_agree`** — the
  eighth seed, which drew a run covering 0.28% of the genome. `F` = 0.0028 against a realised
  0.0028, all nine starts agreeing. **This is the half that stops the criterion being a way to
  refuse everything hard**: five heterozygotes a window is not too little evidence to fit, it
  is where an *absent* signal is read as a total one.
- **`neither_check_refuses_the_closest_legitimate_fit`** — the 5× floor. Asserts both the tied
  spread (< 0.0001) **and** the raw range over all nine (> 0.05), because it is their
  *difference* that carries the argument. A mutation widening the tie gap past 1,473 nats
  would leave every other test in the file green.
- **`f_recovers_a_drawn_genomes_realised_autozygous_fraction`** gained a spread assertion at
  all four nominal levels — the margin the threshold is set against, on this file's own
  genomes rather than only the harness's.
- **Three `const _: () = assert!` blocks** put the measured bounds on both constants at build
  time, so moving either outside its measurement is an `error[E0080]` rather than a red test.
  The same device guards the error-rate ladder's rung count.

### The harness re-run with it in place (note §6.4)

The unit tests say the criterion fires on the fixtures written for it. The sweep re-run says
what it does to every cell that was measured before it existed.

- **The two catastrophic cells now refuse all five seeds** — 3,000 windows × 5 heterozygotes
  and 12,000 × 5, which returned `F` = 0.9912 and 0.9922 on genomes with no runs at all. So do
  five other cells. Across §1, refusals go from **14 of 60 to 48 of 60**.
- **The largest `F` any no-runs genome now returns anywhere is 0.0152**, against a reported
  resolution of 0.0409 at that window count. It was 0.9922. **Every no-runs fit still answered
  comes back below its own reported resolution** — the property that was missing, and the
  reason `resolution` is emitted at all.
- **§3 and §4 are untouched: 0 refusals across all 160 with-runs fits**, same answers, worst
  recovery error 0.0100.

## Recorded deviations

- **The criterion is not the note's.** The note recommends a bare threshold on the raw spread;
  what is built is a threshold on the spread over the tied starts. Adopting the note's form as
  written would refuse the robustness spec §6.2 requires. Recorded in the note itself as §6.
- **A fifth error variant** where arch §5.4 and the plan's A6 name four. The enum is
  `#[non_exhaustive]` and its doc anticipates growth, so this is not a breaking change — but
  the design docs now understate the count, and that is on the drift list below.

## Design-doc drift this creates, for the owner

1. **arch §5.4 and the plan's A6 name four `ParameterEstimationError` variants**; there are
   five.
2. **arch §5.3's contract** rejects a fit "when no start left posterior mass on both states".
   There are now three refusals, and this is the only one of them the architecture names.
3. **spec §6.5's third emitted quantity** is "the spread across starting points — the best and
   second-best `F` and their scores". Two things about that sentence are now measured: the
   spread is a *criterion* and not only a report, and *best-and-second-best* is exactly the
   pair the measurement rejects — at the 5× floor it reports 0.02 where the range over the
   tied starts is 0.30.

## What this does not establish

- **The tie gap is measured at two points**, 0.91 and 1,473 nats, one fixture each. Three
  orders of magnitude apart, so ten nats is not tuned — but nothing measures where between
  them a real genome sits.
- **No real data.** Every genome behind every number here is drawn from the model the fit
  assumes.
- **What it costs G3.** HG002 is outbred, so it may now be refused for disagreeing starts
  rather than for coincident states, or answered. Either is the right answer and a confident
  0.99 is not; which one a real human genome gives is unmeasured.
- **`resolution_at` and `MIN_WINDOWS_TO_FIT_INBREEDING` are untouched**, and the question the
  note raises about them is sharper now rather than answered: the failure a window-count floor
  was standing in for is the one this refuses directly. It is an owner decision and it is
  recorded in `PROJECT_STATUS.md`, not taken here.
