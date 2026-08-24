# ng — calling prerequisites, D1: one mint for a read's error probability

**2026-08-24**, branch `ng-calling-prerequisites`. Step D1 of
[`calling_prerequisites.md`](../../ng/impl_plan/calling_prerequisites.md), against
[`spec/read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §3.2 and §12 test 10.

**How wrong one read is at one place now has a name, and both of the walk's column paths call
it.** Nothing about any number changes; what changes is that there is one definition of the
number instead of two spellings of it.

---

## 1. What changed and why

When the walk folds a read into an observation it charges it an error probability: the worse of
what the instrument said about the bases and what the aligner said about the placement, in log
space. Two paths mint it. The general fold takes the smallest base quality over the read's events
inside the record's footprint; the ordinary-column path takes the one base it is standing on,
because an ordinary column is one position. Both wrote the arithmetic out inline, in two files.

**The reason that has to stop is the calibration the caller is about to do.** The read likelihood
rescales every read's quality by one number per read group,
`fitted error rate ÷ mean minted per-read error` (§3.2). The mean is summed by the parameter
pre-pass and the per-read numbers are minted by the locus generator, so **the numerator and the
denominator are computed in different modules** — and if they are computed by two definitions of
"how wrong is this read", their ratio is not a calibration of anything. §3.2 states the
requirement outright, and §12's tenth test asks for the function to be called from the pre-pass's
side and the generator's side on one read — a test that cannot be written until D2 supplies the
first of those sides.

## 2. Changes made

**[`src/ng/locus_generation/pileup/open_record.rs`](../../../../src/ng/locus_generation/pileup/open_record.rs)**

- `minted_ln_read_error(base_quality, mq_log_err)` — new, `pub(crate)`. Takes a Phred score and
  the read's mapping-quality log error, returns the larger of the two log-error probabilities.
- `ln_bq_for_read` → `min_bq_for_read`, which now returns the window's smallest **Phred score**
  rather than that score already converted. The conversion and the mapping-quality floor happen
  only inside the new function.
- The general fold's mint site calls both in sequence.

**[`src/ng/locus_generation/pileup/fast_column.rs`](../../../../src/ng/locus_generation/pileup/fast_column.rs)**

- The ordinary-column mint site calls `minted_ln_read_error` with the base's own quality. Its
  field comment used to describe the general path's arithmetic in prose and now names the
  function both paths share.

Three tests: what the mint charges, that both column paths charge it, and how the general fold
picks the quality it hands over.

## 3. Deviations from the plan

**The plan says `pub(crate)`, and `pub(crate)` here does not yet mean crate-wide.**
`mod open_record;` is private, so the function is reachable within `pileup` and no further. The
pre-pass's call needs a re-export from `pileup/mod.rs`, and adding one now would be an unused
import that `-D warnings` refuses. **It lands with its first caller, in D2**, which is the step
that makes the pre-pass sum this quantity.

## 4. What the reviews changed

Four agents, each in its own worktree: what must not have moved, every claim re-measured, test
strength by mutation, and where the function should live given the caller that is coming.

**Nothing moved, and the measure is exact.** The parity module keeps a determinism digest over the
`Debug` rendering of every emitted locus, `q_sum` included; it is `1722:f4cc9f1c132ccf66` before
this change and after it. That is a stronger check than the production comparator beside it, which
compares `q_sum` within a relative tolerance of one part in a billion.

**The reversion test says yes, and here it is the wrong question.** Putting the arithmetic back
inline at both call sites leaves all 368 tests green — and no test could do otherwise, because an
extract-function refactor that changes no value has no input that distinguishes "calls the
function" from "computes the identical expression". What is worth defending is the two spellings
*drifting apart*. All four reviews found that the test as first written was blind to that, and that
only the production-parity fixtures caught it — the ordinary-column site held by **1 test in 368**
and the general fold's by 2. The test now drives both paths and catches both, so the number of
tests holding each call site goes from 1 and 2 to 2 and 3.

**Five claims in the comments were wrong, and each is corrected.**

- **The test was billed as §12's tenth test and is not.** That test's two sides are the pre-pass
  and the locus generator, and its headline assertion is that the scaled probabilities' mean comes
  back as the fitted rate. Both sides here are the generator's, and there is no pre-pass side until
  D2. Left uncorrected, a reader marks test 10 done and stops asking for the half that matters.
- **Production's worse-of is not in `var_calling/`**; it is
  [`pileup/walker/open_record.rs:792`](../../../../src/pileup/walker/open_record.rs), over a
  byte-identical copy of the same table. `var_calling/` only consumes the minted number.
- **The pre-pass's accumulator was described in the present tense** and does not exist.
- **The window assertion claimed to rule out taking the window's *last* quality, and did not** —
  the poor event was last. Measured: a build reducing with `.last()` passed the whole test, and
  only the parity fixtures caught it. A third event now puts the poor one in the middle.
- **The fast path was said to differ by having "no window to take a minimum over".** In every
  column that path accepts, the general fold's window holds exactly the same single event — the
  lane refuses any read carrying an indel and any column an open record overlaps. The two reach
  one number by a shorter and a longer route, not over a narrower and a wider window.

**Three doc comments still named `ln_bq_for_read`,** which the rename removed. None was an
intra-doc link, so `cargo doc` was silent about all three.

**The one test was three tests wearing one name, and two of its assertions were not oracles.**

- **It never entered either column path.** It called the shared function twice — once wrapped in
  the window reduction — so its central equality held however either call site was wired. It now
  **drives two walks**: a clean read over plain reference, which takes the ordinary-column path,
  and a read carrying a deletion, which opens a record and takes the general fold. Both are given
  one base quality and one mapping quality and must come back charged the same number. Confirmed by
  measurement: dropping the mapping-quality floor at either call site now reddens the test named
  for those call sites, where before only the parity fixtures noticed.
- **Its expected value was a character copy of the function's body.** That fails when the code
  changes and not when the code is wrong — written with `min` instead of `max`, it would have
  passed just as well. The expectations are now the five literal numbers, and the two that are not
  round were computed rather than typed.
- **It is split three ways**, so that a wrong window reduction is not reported under a name that
  says the two column paths disagree.

**Three gaps closed, each confirmed by applying the mutation and watching the right test alone
fail.** Excluding insertions from the window minimum used to pass (only the parity fixtures
noticed); the zero-quality corner the table's own doc singles out was untested; and only one of the
three event shapes was ever the window's answer. That corner is load-bearing rather than
decorative — the mate-overlap rule silences a losing mate by giving it exactly Phred 0, and since
`ln 1 = 0` is the top of the table a wrong answer there would move every overlapped contribution
from zero to the read's mapping term.

**What the test defends that nothing else did:** `min_bq_for_read(&[], 27) == 27`. Changing the
empty-window fallback to a constant zero leaves every other test in the crate green. That branch is
one the walk cannot reach — the fold skips empty windows — which is why nothing had pinned it.

### What the reviews found about D2, which is worth more than any of the above

**Three findings that change what D2 has to be**, each measured and each with a file and line:

- **An observation's `q_sum` is the wrong sum for the calibration.** §3.2's mean is over error
  *probabilities*; `q_sum` accumulates their logarithms, so `exp(q_sum / n)` is the geometric mean
  and Σ ε is not recoverable from Σ ln ε. At the very spread §3.2 says is the whole call — one read
  at Phred 20 and one at Phred 40 — the arithmetic mean is 0.00505 and the geometric mean 0.001, so
  a scale built from `q_sum` comes out **five times wrong**. §12's ninth test exists to catch
  exactly this substitution. **D2 needs its own sum of `exp` of the minted value, accumulated where
  the mint happens.**
- **The census route cannot call this function at all.** The plan says it "gains its own pair,
  summed over the census sites its fit reads, calling D1's function". `fit_jointly` reads a
  `CohortCensusEvidence` whose per-position unit
  ([`joint/fit.rs:467`](../../../../src/ng/parameter_estimation/joint/fit.rs)) is depth and allele
  counts — **no read, no base quality, no mapping quality**. The last place on that route that sees
  a `SequenceObservation` is the census *writer*. So the accumulator there is a census file-format
  change: a new field on `GenericEvidence`, new bytes in the section encoding, and a bump of
  `census_file.rs`'s `VERSION`, which makes every existing census file a rebuild. **That is a
  decision, not a step**, and it is in §6 and in the owner's list rather than taken here. The
  histogram route is the small one: `accumulate_by_read_group` already walks the complete
  observations per read group, and two more slots there is the "no new traversal" §3.2 asks for.
- **Above about 124 reads a position the two routes' read sets diverge.** The histogram route
  subsamples every site to `MAX_BINNED_DEPTH = 124`
  ([`generic/depth_bins.rs:77`](../../../../src/ng/parameter_estimation/generic/depth_bins.rs)),
  while an observation's `num_obs` and `q_sum` count every read. A denominator summed over 300 reads
  against a rate fitted from 124 is §3.2's "different site sets" in its read-set form. It fires on
  the GIAB trio and never on the tomato cohort at about 3×, which is exactly the shape of defect
  that a benchmark chosen for one corner would hide. **D2 must sum over the kept subsample.**

**One place inside ng where "one shared mint" is still not true**, found while checking the
blast radius: [`locus_generation/ssr.rs:955`](../../../../src/ng/locus_generation/ssr.rs)'s
`ln_p_err_sum` is a third spelling of Phred → `ln(P_err)`, summed over a tract's bases with no
mapping-quality floor at all. Two differences, both measured: the definition, and — setting that
aside — the rounding, since `(-q · ln 10) / 10` and `-q · (ln 10 / 10)` disagree by one unit in the
last place for 23 of the 94 Phred values in `0..=93`. It matters only if the calibration ever spans
repeat tracts, and it is recorded in §6 for whoever brings that path through.

## 5. Validation

All in the dev container, on the tree as committed.

| gate | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::locus_generation` | **370 passed, 0 failed, 1 ignored** (368 before this step's tests were split and one added) |
| `cargo test --lib` | **4,168 passed, 0 failed, 14 ignored**, 576.31 s |
| `cargo doc --no-deps` | 24 unresolved-link errors, 12 redundant-explicit-link-target warnings — the recorded baseline, unchanged |

**Every mutation quoted above was applied to this tree and run.** Across the four reviews, **51
mutations, 47 killed, 4 survivors — and every survivor was proved behaviour-neutral rather than
left as a test gap**: the mint's `pub(crate)`/`pub(super)` distinction, which no caller can observe
because no caller outside `pileup` exists; the argument order of `a.max(b)`, checked exhaustively
over all 68,608 reachable `(base quality, mapping quality)` pairs plus every floating-point corner;
and the table's negative-zero spelling at Phred 0, which every accumulation starts from `+0.0` and
so cannot carry through.

## 6. Follow-ups

- **D2 adds the second caller** and with it the re-export §3 defers. It also needs a **second
  accumulated scalar** — Σ `exp` of the minted value — because `q_sum` gives the geometric mean and
  the scale wants the arithmetic one (§4).
- **The census route's accumulator is a file-format change**, not an addition: `GenericEvidence`
  gains a field, the section encoding gains bytes, and `census_file.rs`'s `VERSION` moves, which
  rebuilds every existing census file. The plan assumed it could call the mint over reads it does
  not have.
- **The mint's home should move when its second caller arrives.** A statistics module reaching into
  `pileup::open_record` is reading the walker's file layout; `locus_generation/read_error.rs` is
  where a shared definition belongs, and `open_record.rs` was released from the byte-for-byte copy
  set at A0 so the table can travel with it. `min_bq_for_read` stays — it is the walker's own
  window reduction, not shared.
- **`ssr.rs`'s `ln_p_err_sum` is a third definition of the same quantity** (§4), unfloored and
  differing in the last place for 23 of 94 Phred values. Owed to whoever takes the repeat path
  through the calibration.
