# Each library's own stutter — implementation plan

*Status: plan, 2026-10-10, for the owner's review alongside its spec.*

This plan turns [`doc/devel/ng/spec/str_library_stutter.md`](../ng/spec/str_library_stutter.md) into
build order. It is not a place for new design: every step cites the spec section it builds. Branch
`stutter-multiplier`, from main at `4b528afe`, worktree `../pop_var_caller-stutter-multiplier`; it is
meant to merge into main on its own, independent of `fit-precision-followup`.

## 1. Principles that fixed the order

- **Types and plumbing first, with no number moving.** Step H1 keeps each library apart in the
  evidence and routes a multiplier of one through the fit; the checksums must pass unchanged. A
  difference there is a defect in the plumbing, found before any number is meant to move.
- **The estimator before the loop that uses it.** H2 builds the multiplier's estimate as functions
  and proves it against drawn truth; H3 only wires it into the fit.
- **Isolate every step that moves a number.** H3 (the fit's numbers) and H5 (calls) each are their own
  commit, with the cross-platform and oracle checksums' movement measured with
  `scripts/promote_ng_oracle.sh` and explained before they are re-recorded.
- **The file between the fit and calling.** H4 writes and reads the multipliers; H5 makes calling use
  them. Measured end to end at checkpoint H.

## 2. Preconditions

- Main at or after `4b528afe`.
- The simulated stuttering benchmark (simulator `examples/ng_str_stutter_benchmark.rs`, driver
  `scripts/ng_str_stutter_benchmark.sh`, scorer `scripts/ng_str_stutter_score.py`) is on branch
  `tract-ceiling`; checkpoint H runs it from that worktree with this branch's binary, without merging
  it here. Its first-pass outputs (`tmp/stutter_bench/first_pass/` there) hold the pooled fit's and the
  true stutter's results to compare against.

## 3. The steps

1. ✅ **H1 — each library kept apart, every multiplier one.** The evidence's per-tract counts are keyed by
   library (read group) rather than summed into the slippage group; a map from library to slippage
   group says which class rates a library reads; a per-library multiplier vector, all ones, enters the
   level each read is scored at (`min(mₗ · Lₛ, LEVEL_CEILING)`). The subset's top-up stays on slippage
   groups. **No number moves: the checksums pass unchanged**, and a test fits one cohort both ways (the
   old summed evidence and the new) and finds every class's numbers equal to the bit. *Source:* spec
   §2.1, §3.1 item 1, §3.5.
2. ☐ **H2 — one library's multiplier, with the classes held.** Per library and class, the expected
   slipped and exactly-read reads under the classes' genotype posteriors; the multiplier maximising
   §3.2's objective plus the pull of §3.3, by Newton inside its bounds; its standard error (§3.4); the
   scaling of §2.2; one library alone in its group defaulted (§2.3). Nothing calls it yet. **Tests**: on
   drawn classes with known rates and multipliers (2 to 50 libraries, 3 and 30 reads), with the
   classes held at their true rates, every multiplier within two standard errors of its truth in about
   95 cases in 100; the objective's slope against a finite difference; the pull's limits (no reads →
   one; many reads → the data's own ratio). *Source:* spec §2.2, §2.3, §3.2–3.4.
3. ☐ **H3 — the alternation.** Round 0 as today; then rounds of multipliers, scaling and warm-started
   class refits until every multiplier moved less than `SETTLED_FRACTION` of its error, or
   `MAX_MULTIPLIER_ROUNDS`; the log says how many rounds, how it stopped, and the multipliers'
   spread. **Own commit; numbers move.** Tests: the drawn cohorts of H2 fitted end to end, classes and
   multipliers together, every multiplier within two errors in about 95 in 100; a one-library cohort
   unchanged after any number of rounds. The oracle cohort's move measured and explained. *Source:*
   spec §3.1, §3.5, §3.6.
4. ☐ **H4 — the parameters file, version 3.** `repeat_tracts.stutter_multiplier_by_read_group`, its
   validation (a positive number; every read group covered; `defaulted` with no count), the reader
   accepting versions 1 to 3 with every multiplier one where absent, the golden files, the explanatory
   text. The file's bytes move (a new table, the version); no number. *Source:* spec §5.
5. ☐ **H5 — calling reads the multipliers.** The frozen level of each cell is the class rate times the
   library's multiplier, capped; `--defaults` and a missing row give one, the latter counted in the run
   report. **Own commit; calls move.** Tests: a two-library locus called with multipliers 0.5 and 2
   scores each library's reads at its own level; the per-locus refit pulls toward the library's own
   level. *Source:* spec §4.

> **Checkpoint H — pause for review.** On the simulated benchmark (branch `tract-ceiling`): the
> recovered multipliers against the simulator's, library by library; calling accuracy with the
> multipliers against the pooled fit and the true stutter, PCR and PCR-free, 50 samples at 5× and 30×,
> and one sample at 15×; the number of rounds and the fit's time against the pooled fit's. **Decide
> `P` and `MAX_MULTIPLIER_ROUNDS`** (spec §9) from it.

## 4. Verification summary

| step | proven by |
|---|---|
| H1 | checksums unchanged; one cohort fitted from summed and from kept-apart evidence, equal to the bit |
| H2 | multipliers within two errors of drawn truth ~95 in 100; slope against finite differences; the pull's limits |
| H3 | the same end to end; a one-library cohort unchanged; the oracle's move explained |
| H4 | golden files; version 1 and 2 files still read; a round trip |
| H5 | two-library scoring at each library's own level; the refit's target; the calls' move explained |
| checkpoint H | the benchmark against the simulator's truth |

## 5. Out of scope

Spec §8: per-library substitution rates, direction and fall-off; slippage groups other than one.
