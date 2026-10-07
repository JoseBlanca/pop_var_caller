# Fit precision, step C3 — the repeat-tract climb stops when every number is within a tenth of its error of the peak

**Date:** 2026-10-01; round limit, review fixes and the oracle run 2026-10-02. **Plan:**
[fit_precision.md](../../implementation_plans/fit_precision.md) step C3. **Spec:**
[fit_precision.md](../../ng/spec/fit_precision.md) §4.3 as amended with the owner at this step, §4.5 item 2.
**Branch:** `fit-precision`. The rule the spec first gave was built and measured first, and stopped climbs short;
that is [fit_precision_c3_stopped_2026-10-01.md](fit_precision_c3_stopped_2026-10-01.md). The owner chose its
replacement (2026-10-01), and this report is about that. Review:
[fit_precision_c3_2026-10-01.md](../reviews/fit_precision_c3_2026-10-01.md), fixes:
[fixes_applied_fit_precision_c3_2026-10-01.md](../reviews/fixes_applied_fit_precision_c3_2026-10-01.md).

## 1. What was built

The repeat-tract fit climbs each stratum's numbers by coordinate ascent. A **round** moves every number once by a
golden-section search, and a **walk** is the rounds from one of three starting points; the best-scoring walk is the
stratum's answer. A walk used to stop when a round raised the log-likelihood **a tract, on average**, by less than
10⁻⁶, or after five rounds.

**Now a walk stops when every number is within a tenth of its own standard error of where the likelihood peaks**,
the rule the SNP/indel fit stops by since step B1.

- **The judgement** (`StratumClimb::judge`). At the walk's point, the curvature and slope of the stratum's total
  log-likelihood, by step C1's central differences. The slope comes from the same evaluations: each number is
  scored one step up and one step down already. The Newton step `I⁻¹ g` over the numbers the curvature identifies
  is each number's distance to the peak. It is carried to the number's own scale by the derivative its error is
  carried by. The walk has settled when every number with an error is within `settled_fraction` (0.1) of it; a
  number without one is settled by definition.
- **When it is judged** (`walk_until_settled`). Once the gain still to come, projected from the last two rounds'
  gains, is below ½ · p · 0.1² (0.08 at thirteen classes). That projection was the spec's stopping rule and is now
  only the trigger. A walk is also judged after a round that loses.
- **A judgement that finds the walk unsettled puts the next one off**, by one round more each time: one, two, then
  three rounds later. A walk crossing a plateau was otherwise judged nearly every round. A judgement is one
  curvature, 513 evaluations at thirteen classes, against 307 for a round (17 golden sections of 18 evaluations, and
  the score where the round ends): about 1.7 rounds.
- **A round that loses** is undone back to the best point it stood at: the round's start, or a point part-way
  through. A point part-way through is judged unless a judgement is put off, and the walk goes on from it unless it
  has settled. The round's start is always judged, put off or not, since a next round would repeat the same round
  exactly; if it has not settled the walk stops there, recorded as having lost a round, which is not convergence.
  A round whose score is not a number counts as a loss.
- **`DEFAULT_MAX_ROUNDS` = 40**, from 5 (§2.3).
- **The winning walk's last judgement gives the stratum's errors** whenever its last round was judged at the point
  it returns, settled or not, so step C1's separate curvature is not paid again. Otherwise they are computed as
  before.
- **The record.** `StratumFit::converged` is replaced by `ending` — settled, lost a round, or out of rounds — and
  `walks`, every walk's ending, rounds and judgements. The run's log has one line counting them.
- `SsrFitConfig::stillness` is replaced by `settled_fraction`. The SNP/indel fit's `SETTLED_FRACTION` is shared
  with the repeat-tract fit.

## 2. What was measured

**The comparison** (`the_climbs_stop_loses_nothing_against_a_longer_climb`, ignored). Each drawn stratum is fitted
by the rule, with room for 60 rounds. The same three walks are then climbed further with no rule, and the best
point any of them reached is the reference, 60 rounds (re-run with the reference set to 60 explicitly,
`tmp/fit_precision/c3/ref60_{3,13}class.log`: every line the same as the logs below). Every number's distance is taken in the reference's standard errors
(step C1's). "Short" is how many total log-likelihood units the fit stopped below the reference. "The old rule" is
the same walks stopped where the 10⁻⁶-a-tract rule would have, read off the same paths. Logs
`tmp/fit_precision/c3/fixes_3class_60.log` and `fixes_13class_60.log`, run with the final code but for the round
limit (60 there, 40 now; §2.3 shows the drawn strata's answers are the same at 40).

### 2.1 Drawn strata

| stratum drawn | numbers further than a tenth of an error | furthest | most short | rounds (judgements), all walks | the old rule's rounds |
|---|---|---|---|---|---|
| 3 classes, 30 tracts × 8 samples × 3 reads, 20 draws | none of 136 | 0.035 errors | 0.0010 | 200 (67) | 295 |
| 3 classes, 300 × 20 × 3, 20 draws | none of 140 | 0.044 | 0.0026 | 210 (83) | 288 |
| 3 classes, 300 × 20 × 30, 20 draws | none of 140 | 0.048 | 0.0017 | 197 (101) | 287 |
| 3 classes, 400 × 63 × 3, 5 draws | none of 35 | 0.058 | 0.0017 | 64 (23) | 72 |
| 13 classes, 200 × 20 × 3, 5 draws | 17 of 84, all in one draw | 0.044 in the other four | 0.0019 in the other four | 246 (56) | 75 |

- **At three classes every number of every fit is within 0.058 errors of the reference**, and the walks that
  settled took 2 to 6 rounds. With the judgements counted (73 evaluations each against 127 a round), the climb
  costs about 830 rounds' worth against the old rule's 942, 12% less. The old rule was within 0.060 errors; the
  projection alone, as first built, left 21 of 195 class shares 0.1 to 0.5 errors short.
- **At thirteen classes four of the five draws are within 0.044 errors**, where five rounds of the old rule were
  0.15 to 4.1 errors away. In the fifth draw the winning walk stopped at a round that lost, 17.4 units short, and
  says so; its other two walks settled, lower still. The old rule was 17.6 short there.
- **The comparison asserts what this section claims**: at the default draws, every fit whose winning walk settled
  must be within a tenth of an error of the reference on every number, and within ½ · p · 0.1² of its
  log-likelihood. Both runs passed.

### 2.2 The four-accession oracle cohort

Four tomato accessions at about three reads a position; fifteen strata fitted on their own tracts, thirteen classes.

- **With the final code** (`tmp/fit_precision/oracle_c3g/fit.log`, limit 40): 45 walks, 36 settled, 4 stopped at a
  round that lost and 5 ran out of their 40 rounds; **394 rounds against the old rule's 217**. The winning walk
  settled in 13 of the 15 strata.
- **Time.** The repeat-tract fit took 5 min 13 s, against 1 min 45 s in step C1's run (`oracle_c1g`): 3.0 times as
  long. That ratio is not firm: the SNP/indel half, whose work this step does not change, took 8 min 24 s in this run
  and 14 min 54 s in C1's, so the host's load differed between the two. The run at a limit of 30 took 5 min 53 s for
  344 rounds. Rounds are the firmer measure: 1.8 times the old rule's, before counting the judgements.
- **Against a 40-round climb** (`tmp/fit_precision/oracle_c3n/fit.log`, limit 40, with the code before the
  judgements were put off and before the review fixes; per-stratum lines in
  `tmp/fit_precision/c3/oracle_newton_lines.txt`). Of the 13 strata whose winning walk settled:
  - six are within 0.1 error of the reference on every number;
  - five differ on some number by 0.37 to 234 errors, while their log-likelihood is 0.000 to 0.054 units short. **In
    these five the rule called a walk settled that was not, by its own measure.** In three the number is a class share
    near zero: 2.5 × 10⁻⁵ against 1.1 × 10⁻³ in the reference (1:12), 1.2 × 10⁻⁴ against 1.3 × 10⁻³ (1:14),
    5.7 × 10⁻⁶ against 6.1 × 10⁻⁴ (3:6). In 2:9 it is a share of 0.031 against 0.028, 0.59 errors. In 2:10 it is a
    share of 0.9997 against 1.0000, where the concentration ran off (36 against 7.4 million in the reference). The
    cost in log-likelihood is at most 0.054 units, under the 0.08 the rule allows in total;
  - two are 5.9 (1:9) and 1.5 (1:15) units short, because some walk's path climbs higher than the winner stopped. In
    1:9 one walk ran out of its 40 rounds; in 1:15 none did, and one stopped at a round that lost. The old rule was
    5.9 and 2.0 short there.

  Of the two strata whose winning walk ran out of rounds, 1:8 is the one the old rule left 21.3 units short. With the
  limit at 40 its winner is the 40-round reference climb itself, so the two agree by construction.

### 2.3 How many rounds

With room for 60 rounds, the walks that settled took 2 to 6 rounds at three classes and 6 to 36 at thirteen. At
thirteen classes one walk of fifteen had not settled by 60 (draw 3, not the winner), and none settled between 37
and 60. **So the limit is 40**: it holds every walk measured to settle. Whether a walk ever settles after 40 is not
known; none did here. At five rounds, 9 of 15 walks at thirteen classes ran out.

At 40 the drawn strata's answers are those of the 60-round runs: the one walk that ran past 40 did not win, and a
walk's score never falls from one round to the next, so cutting it short cannot make it win.

**The put-off lengthens walks.** In the run before the review fixes (`newton_13class_40.log` against
`backoff_13class_40.log`), one walk judged at every trigger settled at round 17 and, with the put-off, at round 29.
On the oracle cohort at a limit of 40, 5 of 45 walks ran out of rounds and the winner settled in 13 of 15 strata, as
at 30. Raising the limit from 30 to 40 moved 19 lines of the parameters file: the length spectra of strata 1:8 and 1:11, whose
own fits moved, and the 17 slippage rows of motif length 1, which the period's curves pass that change to.

### 2.4 What moved

- **The four-accession oracle cohort** (`scripts/promote_ng_oracle.sh`, `tmp/fit_precision/oracle_c3g` against step
  C1's run `oracle_c1g`, which matched the old baseline; `tmp/fit_precision/c3/cmp_c3g.py`): 50 of the parameters
  file's 574 lines moved, all of them repeat-tract rows — 37 slippage rows and the length spectra of 13 of the 15
  strata fitted on their own tracts; the SNP/indel fit's lines did not move. Of the 6,706
  calls made with it, 352 changed printed numbers only: the allele frequency in 324 by at most 0.043, 326 sample
  genotype qualities in 117 records by at most 2, and the site quality in 193 by at most 0.5. No genotype changed,
  and no record appeared, vanished or changed its alleles. The baseline's two fit lines are re-recorded in this
  commit; its other five lines did not move.
- **The cross-platform checksums did not move** (`cli::cross_platform_digests` passes unchanged).

## 3. Tests

| test | what it shows |
|---|---|
| `the_gain_still_to_come_is_projected_from_the_last_two` | the projection, its cap, and a first round |
| `a_walk_is_judged_once_its_projected_gain_is_small_and_stops_when_settled` | judged first where the projection falls below the trigger; stops when a judged point has settled |
| `the_trigger_is_on_the_total_log_likelihood_not_the_mean` | the same mean gains trigger a one-tract walk on its first round and a thousand-tract walk on its fifth |
| `a_round_that_loses_is_undone_and_its_best_point_judged` | back to a better point part-way through and on from it; back to the start and stop; a `NaN` is a loss |
| `an_unsettled_verdict_puts_the_next_judgement_off` | judged at rounds 1, 2, 4 and 7 of ten; a loss with nothing better than its start, while put off, judged all the same and the walk ended |
| `a_walk_still_gaining_at_its_last_round_ran_out` | out of rounds, never judged |
| `a_walk_from_a_start_that_is_not_a_number_lost_a_round` | a `NaN` start loses its first round, is judged, and stops |
| `the_target_counts_the_numbers_the_errors_are_taken_over` | 0.08 at one live group and thirteen classes, 0.095 at two |
| `a_point_one_error_from_the_peak_has_not_settled` | on a drawn stratum, settled at the peak of a 20-round climb; with the level moved by one error, the level 0.7 to 1.3 errors away and not settled |
| `a_walk_at_two_slippage_groups_is_judged_on_both` | the same at two slippage groups, group 1's level moved: group 1 about one error away, group 0 not |
| `newton_distances_are_carried_to_each_numbers_scale_like_the_errors` | `I⁻¹ g` carried to each slippage number, share and the concentration as its error is, on a hand-built coupled inverse, with every coordinate kept and with one dropped |
| `the_furthest_number_counts_every_kind_and_a_distance_not_a_number_is_unsettled` | a share or the concentration can be the furthest; numbers without an error are passed over; a `NaN` distance is unsettled |
| `a_point_that_scores_not_a_number_is_not_judged_settled` | a judgement at a `NaN` point is unsettled |
| `the_held_point_scores_what_it_holds` | the best point a real round passed through, scored afresh, has its held score to the bit |
| `the_settled_fraction_is_read_from_the_config` | at a fraction of 10⁻⁷ no walk settles within eight rounds |
| `the_winning_walks_errors_are_the_errors_at_its_answer_in_both_schedules` | at the default limit, the reused errors equal a fresh curvature's, and both schedules give the same walks and endings |
| `tracts_without_reads_do_not_move_where_a_walk_stops` | 300 tracts no sample read leave every walk's ending, rounds and judgements, and the answer, the same |
| `the_climb_summary_counts_each_ending` | the log line's counts; a hand-built fit with no walks not counted |
| `central_differences_give_a_quadratics_second_derivatives` | now also the slope, to 10⁻⁸ |
| `the_climbs_stop_loses_nothing_against_a_longer_climb` (ignored) | §2.1, with its assertions |

- The two-schedule tests now compare every walk's ending, rounds and judgements as well as the numbers.
- `a_number_the_climb_left_at_the_end_of_its_reach_is_not_placed`: its bounds on where the runaway numbers stop are
  loosened (fall-off below 10⁻³, shorter share above 1 − 10⁻³, concentration above 10⁴), since the walks now stop
  sooner on thin strata; each number is still `NotPlaced`.
- Mutations (`tmp/fit_precision/c3/mut/`): the shares' distances without their cross terms, and the concentration's
  left on the log scale, are each killed by `newton_distances_are_carried_…`; the tract count the walk and the
  curvature share (`tracts_a_mean_is_over`) set to one is killed by three error tests.
- Module (`cargo test --release --lib parameter_estimation::joint::ssr_fit`): 60 passed, 3 ignored.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,010 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 8 ignored.

## 4. Choices, deviations and what is left

1. **The judgement is on each number's own scale**, as B1's and as the comparison measures, not on the climb's logit
   and log scales.
2. **The trigger is the spec's projection.** A lower or higher trigger changes only how often a walk is judged.
3. **The put-off grows by one round** after each unsettled judgement. Measured at thirteen classes before the review
   fixes: 53 judgements where judging at every trigger took 104, for 22 more rounds (226 against 204). Four of the
   five draws' answers were the same; the fifth's moved closer to the reference (0.040 errors against 0.083). One
   walk that ran to the limit took 35 judgements in 40 rounds when judged at every trigger, and 8 with the put-off.
4. **A losing round inside the put-off** is judged only when it has nothing better than its start. With a better
   point part-way through, it goes back there unjudged and goes on.
5. **Not addressed:** strata where another walk climbs higher than the settled winner (§2.2). That is the starting
   points, not the stopping rule. And the 256-point average the likelihood stands on
   ([its measurement](fit_precision_quadrature_average_2026-10-01.md)).
