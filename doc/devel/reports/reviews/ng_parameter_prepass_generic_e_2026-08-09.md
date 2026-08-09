# Milestone E — review synthesis

**Date:** 2026-08-09. **Seven agents over four rounds**, each in its own git worktree against
the step's own commit. **At least 120 mutations run; 41 survived the shipped suite.** Every
survivor now has a test that kills it, or is recorded as unreachable with the reason.

| round | commit | agents | mutations | survived |
|---|---|---:|---:|---:|
| E1 | `25ee9a1c` | 3 | 24 + reach | 9 |
| E2 | `2d9b8095` | 2 | 47 | 19 |
| E3 | `5b9b2d5c` | 1 | 34 | 19 |
| E4 | `e1ad9290` | 1 | 15 | 3 |

**Fifteen numeric or prose claims were wrong. Every one was mine, and every one was a claim
about my own fixture's reach** — never a figure copied from the research note, of which about
forty were checked and all held. That is the eighth round of sixteen on this plan with the
same shape.

---

## 1. The finding that matters most: a fixture that could not see what the step is

**E2's `TwoLibraryWorld` ran at twenty reads a library, and at that depth the two blocks are
not coupled at all.** A heterozygote shows about twenty alternative reads and a sequencing
error nought or one, so the classes never overlap: the frequency climb returns the same
answer whatever error rate it is handed — measured, the frequencies climbed at three times
the true rates and at the true rates agree to ten significant figures — and the read-group
scan returns the right rungs from **every** start on the ladder, rung 0 included, which is a
hundred times the true rate. A whole-sample table claiming thirty heterozygotes in a hundred
still gives rungs 80 and 64.

Four mutations survived on that world, one of them **the estimator arch §5.2 describes and
the research harness rejected** — step 2 re-climbing the frequencies at every rung. That is
the thing E2 exists to build.

**It is the same fault E1 had already hit and I failed to carry across.** E1's own
`the_frequencies_handed_in_move_the_fitted_rate` had to move from twenty reads to six for
exactly this reason, and the lesson did not travel three steps down the same milestone.

`CoupledWorld` — three reads a library at Phred 15 and 17, an individual 7.5% heterozygous —
is where they die. The load-bearing test is two samples with **byte-identical** read-group
tables whose whole-sample tables disagree about the individual; it is the only assertion in
the file that fails when the two blocks are disconnected.

## 2. Two code defects

- **E1: the fixed-frequency scan could fabricate a railed fit.** `scan_ladder` refuses a
  weighted cell no genotype of the *model* can produce, but the mixture the sibling scores is
  the model's likelihoods weighted by the **caller's** frequencies, and a genotype held at
  zero — which the frequency check deliberately allows — can make a cell impossible that the
  model says is fine. Every rung then scores −∞, the tie rule hands back whichever came last,
  and `argmax_at_ladder_end` is set: the exact shape of a read group whose true rate lies past
  the ladder, which is the one thing that flag exists to tell apart. And it converges, so E2's
  stopping rule would accept it. Guarded on the **winner** rather than per rung, so a single
  −∞ rung — an ordinary answer — is not refused.
- **E3: `resolution` was read from the padded chain.** The window floor is checked against the
  windows that hold sites; the resolution two lines later was taken from the chain's whole
  length. On a region-restricted run — 3,601 windows of evidence plus one contig whose only
  site sits at window 8,000 — that reports 0.006768 where the truth is 0.029067, understated
  **4.3-fold**, on the one number a consumer uses to decide whether a small `F` means
  anything.

## 3. The Blocker that is recorded rather than repaired

**A genome drawn with no runs at all returns `F` = 0.998**, converged, unflagged, with a
reported resolution of 0.029. `separated_states` cannot see it: windows do cross 0.5, because
the two states are an arbitrary split of sampling noise and the label "inside" lands on the
majority. The guarded direction is *a failed search returns zero*; this is *an absent signal
returns one*.

The mechanism is exact and general: **`MIN_WINDOWS_TO_FIT_INBREEDING` and `resolution_at` are
functions of the window *count* alone, and the noise floor depends on the evidence *per*
window.** Research note §3.6's eight seeds were drawn at 100,000 sites a window; the fixture
that produces 0.998 carries 400, about 250 times less, at 3,600 windows — above the floor of
3,000. Neither cohort in hand is near that regime.

**The repair is a measurement, not a patch** — the floor as a function of evidence per window
— and that measurement does not exist. Recorded on `resolution_at` and here for the owner.

## 4. Nine things nothing could see

Each survived the whole suite and each now has a test.

1. `ReadGroupErrorRateFit::log_likelihood` could be the constant `LogProb(0.0)` — and it is
   the field E2 picks its best-scoring iterate by.
2. A group **far below `MIN_SITES_TO_FIT`** must still be fitted at E1, because E4's gate
   needs the fit and the count to decide against them. Every other fixture held 40,000 or
   200,000 sites.
3. The read group inside the noise parameters is **inert on the read-group path** — every cell
   there is pooled, and the pooled branch reads a share and a rate and never a label — so a
   ladder hoisted out of the per-group loop fits every group correctly.
4. **No fixture had more than one ploidy**, so `climb_frequencies`' per-ploidy loop ran
   exactly once.
5. **Nothing pinned where the coupled loop starts**; every start reaches the same answer on
   the shallow worlds, so only the round count can see it.
6. **`fit_coupled`, the public door, had zero coverage** — replacing its whole-sample map with
   an empty one left all 3,147 tests green. Its stated justification (that reaching it needs a
   locus stream) was false.
7. **`DepthAltHistogram::total_reads` summing only the pooled arm** survived, because every
   fixture the fits build enters through `add_site` and so has an empty attributed map.
8. **`F`'s coverage weighting had no test that could fail**: every window held 400 sites at one
   reference position each, so weighting by covered positions, by site counts, by an unweighted
   mean, and returning `enter/(enter+leave)` are the same number. The new fixture separates
   them at 0.03846, 0.66667, 0.66667 and 0.66102.
9. **The `supplied` parameter E4 added was never exercised non-empty** — all fifteen call sites
   pass an empty map.

Plus two rules that were unreachable and are now named functions with direct tests: the runs
model's **ordering constraint** (every start builds its inside state below the outside one, so
no drawn genome crosses over) and its **best-start selection**.

## 5. Fifteen wrong numbers, all mine

| claim | truth |
|---|---|
| E1: "a library **ten times** noisier than the ladder's worst rung" | **3.2 times**, five Phred past the end |
| E1: "**three** mutants die here and nowhere else" | they die in five, four and three tests; the right number is zero |
| E1: a big diploid arm "is how the first version passed a scan that read one set for both ploidies" | does not reproduce — at `[500, 300, 200]` the test *fails on correct code*, and the named mutant dies on the width check |
| E2: reads and sites "differ by a factor of **forty**" | **twenty**; forty is the two libraries pooled, which is neither estimate's number |
| E2: "the second scores higher than the first" | **bit-identical**, gap exactly 0.0 — so the test could not fail, and `>` was keeping the earlier, less settled iterate |
| E2: `between * 1.001` discriminates log-space from probability-space nearness | it does not; the band is **4 parts in 10,000** and 1.001 overshoots it |
| E2: `TwoLibraryWorld`'s share-weighted-rate assertion | cannot fail for the reason given — `p_j(ε)` is affine in ε |
| E3: `resolution_at(8_004)` = 0.00999 | **0.009970** |
| E3: `resolution_at(31_000)` = 0.00292 | **0.002914** |
| E3: drawn fractions miss nominal "by as much as **0.05**" | **0.027** |
| E3: "the **same genome** is drawn at each floor level" | four different genomes — realised 0.2056, 0.3153, 0.3253, 0.3122 |
| E3: `MEAN_RUN_WINDOWS_AT_START` is "the **middle** of the range" | the **shortest**; the harness uses 20, 30 and 50 |
| E3: "a few dozen calls rather than one per site" | `binomial(n, p)` consumes `n` uniforms |
| E3: "settle in three or four passes" | four to six |
| E4: "module **275 → 276** tests, plus six in the new file" | **269 → 276**, seven added |

Two more were *stories* rather than numbers, and a wrong story is worse: E1's ploidy-arm
mechanism (above) and E3's claim that a collapsed fit can outscore a real one — measured, the
collapsed starts score **2,165 nats worse**, and a collapsed fit is a one-state model, a
constrained special case, so its likelihood *cannot* exceed the two-state maximum.

## 6. What the reviewers proved could not be tested

Recorded so the next reader does not spend a round re-finding it.

- **The runs model's contig restarts change no answer.** The cross-boundary transition term is
  normalised by the wrong contig's total, so it underflows to exactly zero; and an absent
  window has zero covered positions, so it cannot move `F` whatever the chain believes.
- **Summing the runs model's log-likelihood at every window instead of once per contig** — a
  ~300-fold error — survives, because nothing asserts a log-likelihood value and every start
  is scaled alike.
- **The coupled loop's block *order*** cannot be told from the architecture's order: they share
  a fixed point and differ only in which half of a round an unconverged answer came from. The
  module doc now says the order follows the harness and that only the *not re-climbing* is a
  property the tests hold.
- **Keep-the-best-scoring *start*** cannot be separated from keep-the-first by any fixture, for
  the likelihood reason above.
- `resolve_error_rates`' `lender != group` filter is **dead code** — proven by an assertion
  probe over the whole library — because a group that qualifies never reaches that branch.

## 7. Design-doc drift, for the owner

Fifteen locations across four documents. The load-bearing ones:

1. **arch §5.2 describes a different estimator** in three places — step 1's re-climbing, the
   single-library one-iteration claim, and the `MAX_COUPLED_FIT_ITERATIONS` note that rests on
   it. The plan's **E2 step description** states the reverse order of what was built, and its
   second oracle is false of the alternation.
2. **arch §5.2's `fit_coupled` signature** is the only one in `doc/`, and now lacks `sample`,
   `ladder` and `supplied`. Arch §5.3's `fit_inbreeding` likewise.
3. **arch §5.3 says the transition rates are fitted per base.** They are fitted per window and
   converted for reporting, which is what the harness does too.
4. **spec `parameter_prepass.md` §1043-1050 and `parameter_prepass_cohort.md` §317-319** put
   the read-group borrow at the cohort gather; arch §5.4 puts it within the sample, and the
   code follows arch. Spec and arch disagree.
5. **arch §2.4 says an error rate's `observations` is reads**; E1's `sites` field is a site
   count, and E2 is where the conversion happens. Now consistent, but undocumented.
6. **The arch module tree and the plan's scope list** still name five files under `generic/`;
   there are eleven.
7. **arch §5.4's "ties resolve to the lower error rate"** sits three paragraphs below the
   fallback ladder and is about the *rung* ladder. "Ladder" is triple-loaded in that document:
   depth bins, error-rate rungs, fallback.
8. **`GenericEstimationConfig` (arch §1.1) has no supplied-error-rates field**, so E4's third
   rung has no documented source.
