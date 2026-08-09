# Review of the runs model's third refusal — synthesis and fixes applied

**Date:** 2026-08-09. **Reviewed:** `2ceb00a6`. **Three agents in isolated worktrees**:
mutation testing, numbers-against-sources, and design/API. Verdicts:
Approve-with-changes ×3. **1 Blocker, 4 Majors**, plus four wrong numbers and two wrong
stories.

Every finding below was re-verified before being acted on, and one reviewer suggestion was
**measured and rejected** (§5).

---

## 1. The Blocker: the new error's message and payload were entirely untested

The mutation agent set `tied_starts: 0`, `starts: tied`, `threshold: 0.0` and replaced the
whole `#[error]` format string with a stub keeping none of its wording. **All 291 tests
passed.**

That message is the variant's whole purpose. It ends *"this is a genome the model could not
read, not an inbreeding coefficient of zero; supply F instead"*, and a caller who reads it as
*nothing was found* will divide the cohort's diversity by `1 − F` with a defaulted `F` — the
harm the refusal exists to prevent. The sibling `InbreedingStatesNotSeparated` has two guards
on exactly this point; the new variant had none, because its test destructured `spread` and
never called `to_string()`.

**Fixed** by asserting all five fields and three properties of the message in the refusal
test, plus a new `mod.rs` test that additionally requires the two refusals to **render
differently** — a log that cannot tell them apart cannot be acted on, since one says widen the
separations and the other says there is no signal to widen towards.

## 2. Major: the threshold rested on a censored sample, and a build-time assert pinned the
wrong number

`const _: () = assert!(MAX_IDENTIFIED_START_SPREAD < 0.1147, "the narrowest disagreement
measured on a genome with no runs is 0.1147 …")`.

**0.1147 was the narrowest spread among the fits that this very threshold had refused.** The
sample was censored by the number it was being used to justify — the eight-seed sweep behind
it could only report a spread through the error, and any seed below 0.05 was accepted and
never counted. Both the design agent and the numbers agent caught it from opposite directions;
the design agent found accepted no-runs fits at spreads up to 0.0478 at sixteen seeds.

**Re-measured uncensored**, 24 seeds at the fixture's own shape with no runs at all: seven
refused for never separating, and the rest at **0.0213, 0.0264, 0.0323, 0.0337, 0.0483,
0.0521, 0.0541, 0.0584, 0.0621, 0.0675, 0.0722, 0.0740, 0.1027, 0.1144, 0.1368, 0.1375,
0.2749** — a continuum, with no gap at a twentieth or anywhere else.

**So the threshold is not separating two clusters**, and the doc now says so. What makes any
choice in that interval safe is the *other* side: the legitimate population is not merely
small but **exactly 0.0000** across 160 harness fits, and 3.7 × 10⁻¹¹ on this file's own
fixture. The two build-time bounds are now honest and unequal in strength — a tight lower one
(3.7 × 10⁻¹¹), and a deliberately **weak** upper one that says only that the criterion still
refuses the worst case measured, with a comment explaining why no tight upper bound exists.

**The reviewer's proposed replacement justification was also checked and does not hold.** They
suggested the accepted no-runs fits are harmless because each returns `F` below its own
reported resolution. At this fixture (3,600 windows, resolution 0.0291) the five accepted fits
return 0.0019, 0.0213, 0.0322, 0.0359 and **0.0505** — three of the five are above it. That
claim is true of the harness's §1 cells and false here.

## 3. Two Majors: neither constant was behaviourally constrained anywhere in its legal range

- `MAX_IDENTIFIED_START_SPREAD` survived mutation to **0.001** and to **0.1146**.
- `MAX_TIED_START_LOG_LIKELIHOOD_GAP` survived **0.92** and **1400.0** — the entire span its
  build-time asserts permit.

The cause is the same for both: every fixture sits either nine orders of magnitude below the
threshold or twenty times above it. Nothing lived near either constant, so the only thing
holding them was a literal pin, which says a value was typed twice and nothing about whether
it is right.

**Fixed** with two tests built from written-down outcomes rather than drawn genomes: a slice
straddling the tie gap at 5 and 50 nats, and a pair of spreads at 0.0213 (answered) and 0.0521
(refused), both measured values from the uncensored sweep. The second test's doc states
plainly that it pins where the threshold was put and does not show that a twentieth is where
it belongs.

## 4. Wrong numbers and wrong stories

Four numbers, all of them the author's own claim about their own fixture's reach — which is
the eighth round on this plan with that shape:

1. **"160 fits at twenty shapes"**, in five places. 160 fits are twenty run shapes crossed
   with §3's twelve window-count shapes; twenty shapes alone is 100 fits.
2. **"at the 5× floor it reports 0.02 where the range over the tied starts is 0.30"**. Neither
   number comes from that fixture — its six tied starts all return 0.315740, so both readings
   give 0.0000. The pair is from the hand-written four-outcome unit test.
3. **"two of them fits that recover their genome"** — all three of the refused tests do.
4. **"the nine starts return `F` = "** followed by seven values. Seven *distinct* values across
   nine starts; `F` = 0.0010 occurs three times.

Two stories, which this project rates worse:

5. **"an absolute number rather than a fraction of the total because a likelihood ratio does
   not grow with the genome"** is backwards. A log-likelihood difference between two genuinely
   different fits *does* grow with the data — 1,473 nats here would be about 2,900 on a genome
   twice as long. The correct argument is the mirror image: in the **failure** case the two
   fits are the same distribution, so their gap stays of order one at any genome size, which
   is exactly why an absolute cut separates them at every scale. Also dropped: "an odds ratio
   of about 22,000 to one", which reads a likelihood difference as a Bayes factor that the
   boundary-sitting competitor does not license.
6. **"the two states come out at 0.31 and 0.62 of each other"** attached to the disagreement
   fixture, in three files. Those are the *means of two harness cells*; the `runs.rs` fixture
   has its winning start at a ratio of **0.086**. It sends a reader hunting a fixture that does
   not exist.

## 5. One reviewer suggestion that measurement rejected

Both the design and numbers agents observed that `tied_starts` anchors on `starts_tried[0]`,
which need not be a start that separated its states, while the reported `F` comes from
`best_separated_start`. The design agent proposed carrying `separated` on `StartOutcome` and
taking the spread over **tied and separated** starts.

**The mutation agent had already run exactly that mutation and it is killed** by
`starts_that_land_in_different_places_are_refused_rather_than_averaged`: on the no-runs fixture
the disagreeing starts are the *non-separated* ones, so filtering them makes the spread small
and the fit is accepted. Adopting the suggestion would disable the criterion on the fixture it
was built for.

Recorded rather than applied. The anchoring concern is real but narrow, and no fixture reaches
it.

## 6. Smaller fixes applied

- **A `NaN` fitted `F` was silently dropped from the spread** — `f64::max` ignores `NaN`, so an
  all-`NaN` tie set folded to `−∞`, compared false, and **accepted**. Now the spread is `NaN`
  and the check tests `is_nan()` first, because every comparison against `NaN` is false.
- **`RunsModelFit::tied_starts()` is exposed**, and the tie-set size asserted wherever the
  margin is claimed. A spread of zero over one start is an empty comparison and prints the same
  0.0000 as nine starts agreeing — and it happens: at a realised `F` of 0.95 under a five-fold
  false-heterozygote floor, one start in nine ties.
- **`InbreedingStatesNotSeparated`'s doc** named one of its two conditions; it has refused on
  the state-ratio condition since E3.
- **"Four properties" over three bullets**, and note §6.1's prose implying six run shapes where
  there are five.

## 7. What the reviews confirmed

Worth recording, because it is where most of the effort went. The numbers agent verified **31
claims exactly**, including every number that decides behaviour: the 0.9106 and 1473.2179 nat
gaps, 14 → 48 of 60 refusals, 0.9922 → 0.0152, seven of eight seeds, and both test counts. The
mutation agent killed **eight of eight** mutations to how the criterion is *computed*, one of
them caught by four independent tests.

And the design agent answered the question that mattered most — whether this refuses a selfing
landrace, which is what tomato's actual data is. **It does not**: 108 fits at a realised `F` up
to 0.9719, across three evidence levels and false-heterozygote floors of 0 to 5 times the real
rate, all answered with recovery error ≤ 0.0074 and **every one at a tied-start spread of
exactly 0.0000**. There is a structural reason — disagreement needs the likelihood flat in `F`,
which needs the two states to coincide, and high `F` does the opposite by putting more data
into the inside state.

## 8. Verification of the fixes

All six mutations that survived the reviews were re-run against the fixed code and **all six
now fail tests**. One needed a second attempt: the first patch missed because `rustfmt` had
reflowed the closure, and the run reported "pattern not found" rather than a green suite —
which is the only reason it was not mistaken for a survivor.

Gates after the fixes: `fmt` clean, `clippy --all-targets --all-features -D warnings` clean,
`test --lib --bins --tests --all-features` **3,195 passed / 0 failed / 5 ignored**, `doc
--no-deps --lib` at its 12-unresolved-link baseline with none in this module.
