# Code Review: snp_fit_convergence, plan steps A1, A2 and B1
**Date:** 2026-10-16 (measurements finished 2026-10-17)
**Reviewer:** rust-code-review skill (single reviewer; the categories were covered directly, without sub-agents, as the lead asked)
**Scope:** commits 558263c6 (step A1), ec940a7b (step A2), 93f22f6c (step B1) on `snp-fit-convergence`, and two defects the implementer found
**Status:** Request-changes

---

**How to read the labels.** Finding codes (B1, B2, M1, …) are this review's. Plan steps are always
written *step A1*, *step B1*. "The decrement" is the Newton decrement: the log-likelihood the
quadratic model says the remaining Newton step would gain, which a start now stops on below 0.001.
"The box solve" is `bounded_step`, which finds that step while keeping every parameter inside the
interval the fit allows it.

## The short version

The three steps do what the plan and the spec say, on the inputs the plan measured: the 1- and
4-sample oracle fits reproduce the implementer's files byte for byte, the judging pass runs in the
spec's order, the decrement is computed as specified, and the thread count does not change a bit.
Two defects decide whether a start converges on small cohorts, and both were traced:

1. **The box solve often fails, and when it fails it returns a step outside the box.** On the
   two-sample fixture with sequencing errors, 95 of its 100 solves end at their round limit; on the
   fit module's own 3-sample fixture, 75 of 102. The decrement computed from such a step means
   nothing (290, −36, −1,100 were all seen). That is why the two-sample fixtures stopped converging
   at step A2. Replacing the solve with the textbook primal active-set method removes every failure
   measured and leaves the oracle's pass counts unchanged.
2. **The "50 units lost every other cycle" comes from the error-rate update, not from SQUAREM or the
   Newton step, and it is older than this branch.** When one library's mismapped-position error rate
   comes out below its ordinary-position rate, the update swaps the two. That swap is not a harmless
   relabelling: it hands the ordinary positions (99.7% of them here) the rate fitted to the 0.3%
   that are mismapped, and the log-likelihood falls by 50; the next pass undoes it. Deleting the swap
   removes every loss measured. The same numbers appear at main 51b7be99.

With both fixed, the two-sample fixtures still crawl, because along directions the data do not
determine the remaining step is still worth 0.02 to 0.5 units after 200 passes. Measuring the
decrement only over the parameters that have a standard error makes every digest fixture converge,
in 13 to 84 passes, and changes nothing on the oracle. That last change needs the owner's decision.

---

## 1. Scope

- **Reviewed:** the diff 51b7be99..93f22f6c, three commits.
- **Against:** branch `snp-fit-convergence` at 93f22f6c; the plan
  [`snp_fit_convergence.md`](../../ng/impl_plan/snp_fit_convergence.md) steps A1, A2, B1; the spec
  [`fit_precision.md`](../../ng/spec/fit_precision.md) §2 (amendment of 2026-10-16), §3.5, §3.6 and
  [`parameter_prepass_joint_fit.md`](../../ng/spec/parameter_prepass_joint_fit.md) §3.3.1.
- **In-scope files:** `src/parameter_estimation/joint/fit.rs` (`maximise`, `Judgement`,
  `try_the_newton_step`, `newton_target`, `print_the_cycle`, the tests), `fit/standard_errors.rs`
  (`newton_step`, `NewtonStep`, `decrement_of`, and `bounded_step`, which those call),
  `fit/settled.rs`, `fit/information.rs` (re-fixtured tests), `src/cli/cross_platform_digests.rs`
  (notes and checksums), `src/calling/parameters_file/` (the note's text).
- **Out of scope, with reason:** the parameters-file format and the verdict tables (steps C and D,
  not built yet); the repeat-tract fit (untouched).
- **Categories covered:** reliability (always), errors, naming, defaults (a new config field),
  idiomatic, refactor safety, module structure (four files), float portability (new arithmetic
  reaches the fitted numbers), smells, extras (the diff against its stated intent).
  `unsafe_concurrency` does not apply: the diff adds no threads, locks or shared state; its new
  arithmetic is sequential over sums the existing pass makes. `tooling`: no `Cargo.toml` change.

## 2. Verdict

**Request-changes.** The steps are faithful to the plan, but B1 and B2 below leave every start of
every digest fixture (2, 2 and 21 samples) at the pass limit, 1.6 and 63 log-likelihood units below
what the corrected fit reaches. The caller is committed to cohorts from one sample up
(`CLAUDE.md`, "the range, not the example"); a small cohort that never converges, and says so, is
the honest outcome of these defects, but a fit 63 units short is still a wrong answer.

## 3. Execution status

Commands run (macOS arm64, native cargo, scratch build directories under `tmp/review_AB/`):

| command | exit | result |
|---|---|---|
| `cargo fmt --check` (worktree) | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` (clean copy of 93f22f6c) | 0 | `Finished`, no warnings |
| release lib tests `parameter_estimation::joint` (copy of 93f22f6c plus probes, no probe switched on) | 0 | `440 passed; 0 failed; 1 ignored` |
| the three `cli::cross_platform_digests` checksum tests (same build) | 0 | all `ok` (the probes change nothing when off) |
| `estimate-parameters` on 1 and 4 oracle samples, 4 threads | 0 | both parameters files byte-identical to the implementer's final runs (`tmp/impl_AB/runs/final_{s1,o4}_d0.001.toml`) |

Not run: the full suite (`cargo test --all-targets --all-features`) — the implementer's step-B1
commit reports it; `cargo doc`, `cargo audit` — no dependency or public-doc change. The Linux
container was not used. **Needs verification:** 1 finding (M1's cost at thousands of samples).

## 4. Open questions and assumptions

1. **Should the decrement count only directions the data determine?** (M1.) The spec's guarantee
   is stated for "every determined quantity" (§3.6 item 4); the decrement as built also counts
   parameters whose standard error is absent or wider than their range. Recommendation: yes, as
   defined in M1; it needs a spec amendment, and plan step C2's `cohort_covariance` is the natural
   home for the computation.
2. **Must a library's mismapped-position error rate stay above its ordinary-position rate?** (B2.)
   The swap enforces that order; nothing outside `fit.rs` reads the per-library mismapped rate (only
   the mismapped share and each position's posterior leave the fit). Recommendation: drop the order;
   if a reader needs it, enforce it as a constraint on the update *and* on the Newton target, not as
   a swap.
3. **Replace the box solve, or keep it and fall back to the new one only when it fails?** (B1.)
   Replacing moves the oracle's fitted numbers in their eleventh significant figure (pass counts
   unchanged); the fallback keeps the oracle bit-identical. Recommendation: replace — step E1
   re-records the oracle once anyway.

## 5. Top 3 priorities

1. **B1** — replace the box solve; it fails on 95 of 100 solves on the two-sample fixture and its
   failed answer is what the start stops (or fails to stop) on.
2. **B2** — delete the error-rate swap; it makes the alternation lose 50 to 909 units in a pass,
   older than this branch, and SQUAREM accepts the loss unchecked.
3. **M1** — owner decision on counting only determined directions; without it the corrected fit
   still runs the two-sample fixtures to the pass limit.

## 6. Findings

### Blocker

- [src/parameter_estimation/joint/fit/standard_errors.rs:636-704](../../../../src/parameter_estimation/joint/fit/standard_errors.rs#L636-L704) — **[Blocker]** B1: the box solve fails often on small cohorts, and its failed step lies outside the box
- **Confidence:** High
- **Problem:** `bounded_step` fixes *every* parameter that crosses a bound in one round and frees
  every one whose slope points back in, and repeats for at most `3 × side + 3` rounds. Nothing makes
  that sequence converge, and on near-singular information it cycles. When it stops at the round
  limit it returns the step of its last round, computed *before* that round fixed the parameters
  that crossed — so free parameters are carried past their bounds. Measured by counting the solves
  that reached the limit:

  | input | solves at the limit | example of the step returned |
  |---|---|---|
  | two samples, sequencing errors (digest fixture) | 95 of 100 | density shape b at 36.9 given +9.8 × 10⁶ (bound 50) |
  | two samples, no errors (digest fixture) | 164 of 180 | carrier shape b at its bound 50 given +1.8 × 10³ |
  | 3 drawn samples (`a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width`, 60 passes) | 75 of 102 | decrements −36.4, −37.6, −8.2 at passes 20, 27, 36 |
  | 1 oracle sample | 0 of 27 | — |
  | 4 oracle samples | 0 of 69 | — |

  The decrement computed from a failed step is meaningless and has any sign: 290 (the digest
  fixtures), −36 (above), −1,100 (seen once). A negative one is refused by `has_converged`, as the
  spec's trap says; a positive one is not recognised as a failure at all, and one that happened to
  land in [0, 0.001) would end the start as converged. Step B1's Newton step then moves towards this
  step's clamped target.
- **Why it matters:** this is why the digest fixtures went from converging (18 passes) to the pass
  limit at step A2, and why step B1's Newton tries are mostly refused there (26 to 27 of 29 judging
  passes on the fixture without errors, per the implementer's notes). The code predates the branch (plan step B1 of
  `fit_precision`), but this branch made it decide when a start stops.
- **Suggested fix:** the primal active-set method, prototyped as `primal_box_step` in
  `tmp/review_AB/prototype_with_probes.diff`. From the feasible point `d = 0`: solve the free
  parameters' Newton step from the current point; move along it until the first bound it meets; fix
  that one parameter there and solve again; after a whole step, free the one fixed parameter whose
  slope points furthest back into its interval (with a tolerance of 10⁻⁹ of the slope's size, so
  rounding cannot free and re-fix the same parameter forever); stop when none does. The model's gain
  never falls along the way, so **the decrement can no longer be negative**, and every step it
  returns is inside the box. Measured:

  | input | primal active-set |
  |---|---|
  | every digest fixture, every fit test of the module, 1 and 4 oracle samples | 0 solves unfinished |
  | 1 oracle sample | parameters file byte-identical to step B1's |
  | 4 oracle samples | same passes (53 converged, 46 and 37 agreed); multipliers, inbreeding coefficients, concentrations and the fallback share move by at most 1 × 10⁻¹¹ of their value; repeat-tract slip rates by at most 2 × 10⁻⁵ |
  | `parameter_estimation::joint` release tests | 440 passed, 0 failed |

  Also: drop or rewrite the spec's trap "a negative λ means the solve ended at its round limit"
  (`fit_precision.md` §2) and the matching doc on `NewtonStep::decrement` and `has_converged`, which
  imply that a non-negative decrement means the solve finished — false today (290), and moot with
  the new solve. If the owner prefers bit-identity until step E1, call `primal_box_step` only when
  `bounded_step` reports that it reached its round limit (the prototype's `bounded_step` returns
  that flag).

- [src/parameter_estimation/joint/fit.rs:4334-4344](../../../../src/parameter_estimation/joint/fit.rs#L4334-L4344) — **[Blocker, pre-existing]** B2: the error-rate swap makes the alternation lose up to 909 units in a pass
- **Confidence:** High
- **Problem:** after each library's two error rates are fitted, `maximisation` swaps them when the
  mismapped-position rate came out below the ordinary-position one ("Keep them ordered so a swap
  does not read as a fit"). The mixture's share of mismapped positions is not swapped, and neither
  are the other libraries' rates, so this is not a relabelling: the ordinary class, nearly every
  position, receives the rate fitted to the few mismapped ones. On the two-sample fixture with
  errors, library 0's ordinary rate alternates pass by pass between 2.44 × 10⁻³ and 2.64 × 10⁻⁴,
  and the log-likelihood with it, −482.38 and −532.0 (trace, passes 150–155). See section *The two
  defects, traced* below for the full account.
- **Why it matters:** the expectation–maximisation guarantee — a plain pass never lowers the
  log-likelihood — is what SQUAREM relies on when it accepts the plain two-step point without a check
  ([fit.rs:2557-2559](../../../../src/parameter_estimation/joint/fit.rs#L2557-L2559), "`θ₂` never
  needs a check"). With the swap it does not hold, and a fit oscillates between two states forever.
- **Suggested fix:** delete the swap loop (lines 4334-4344). Measured: zero losses of more than one
  unit on all three digest fixtures (from 37 and 53), all 440 joint-fit tests pass, and no swap ever
  fires on the 1- or 4-sample oracle runs, so their output cannot change. Correct the comment at
  2557-2559, which states a guarantee the M-step does not give exactly (the `JUMP_SLACK` doc already
  records a 0.025-unit drift).

### Major

- [src/parameter_estimation/joint/fit/standard_errors.rs:586-601](../../../../src/parameter_estimation/joint/fit/standard_errors.rs#L586-L601) — **[Major]** M1: the decrement counts directions the data do not determine, so tiny inputs crawl to the limit even when B1 and B2 are fixed
- **Confidence:** High for the measurement; Medium for the proposed definition (owner decision)
- **Assumptions:** that the spec's intent is the guarantee of §3.6 item 4, "every determined
  quantity within a tenth of its error"; open question 1.
- **Problem:** with the box solve and the swap both fixed, the two-sample fixtures still reach the
  limit: the decrement falls only to about 0.02 by pass 200 (start 1, fixture with errors), because the parameters it is still
  gaining along — the density's shapes against their bound, the duplicated class — have no standard
  error (12 of 16 parameters on the fixture with errors; all 16 on the one without). The spec
  already calls such parameters "settled by definition" (§2); the decrement does not.
- **Why it matters:** a one- or two-sample cohort is in the caller's range, and the fit spends 200
  passes per start buying hundredths of a unit along directions it then reports as undetermined.
- **Suggested fix:** measure the decrement over the determined parameters only, with the others free
  to adjust:

  `λ_determined = ½ d_Kᵀ (Σ_KK)⁻¹ d_K`

  where `K` is the parameters that have a standard error, `d_K` their part of the (whole) bounded
  step, and `Σ_KK` their block of the inverse information. In words: the largest
  `½ (Δh / SE_h)²` over every quantity `h` built from determined parameters — exactly the spec's
  √(2λ) guarantee, restricted to the quantities it promises. It equals the current decrement when
  every parameter has an error, so nothing changes on cohorts where everything is determined.
  Prototyped as `determined_decrement` (used only when some parameter with a distance lacks an
  error). Measured with B1's and B2's fixes:

  | input | passes by start, before → after | best log-likelihood |
  |---|---|---|
  | two samples, no errors | 200, 200, 200 (limit) → 16, 16, 16 converged | — |
  | two samples, with errors | 198, 199, 200 (limit) → 22, 13, 13 converged | −483.9 → −482.27 (200 passes of the corrected fit: −482.16) |
  | 21 samples | 199, 199, 198 (limit) → 84, 38, 44 converged | −5108.5 → −5045.75 (200 passes: −5045.64) |
  | 1 oracle sample | 22, 20, 21 → same | file byte-identical |
  | 4 oracle samples | 53, 46, 37 → same | identical to the primal solve alone (B1) |
  | joint-fit release tests | 440 pass → 440 pass | — |

  **Rejected alternative, measured:** simply leaving the parameters without an error out of the
  Newton system. It converges the fixtures too (16; 43, 21, 15; 95, 62, 59), but it holds those
  parameters still while the expectation–maximisation keeps moving them, and that drags the
  determined ones: `a_settled_fit_lands_within_a_tenth_of_an_error_of_a_long_one` fails (the
  4-sample case with the duplicated class stops with a parameter 0.41 of its error from the long
  fit), and the 4 oracle samples take 56, 49, 37 passes instead of 53, 46, 37. **Also measured and
  rejected:** treating a failed solve as "cannot judge" (decrement NaN). It answers nothing: no
  start can stop, and every digest start still runs to the limit.

  **Cost — needs verification.** The prototype inverts the dense matrix, fine to 59 samples (185
  parameters) but cubic in the parameters. Without a bound binding the same number is
  `λ − ½ g_Uᵀ I_UU⁻¹ g_U` (the whole decrement less what the undetermined parameters `U` could gain
  alone; block-inverse identity), which needs only the few undetermined parameters' block. Its
  numerical behaviour on a near-singular `I_UU` was not measured.

- [src/parameter_estimation/joint/fit/information.rs:3974](../../../../src/parameter_estimation/joint/fit/information.rs#L3974) and [fit.rs:2413-2538](../../../../src/parameter_estimation/joint/fit.rs#L2413-L2538) — **[Major]** M2: four of five mutations of the new control flow survive all 440 joint-fit tests
- **Confidence:** High
- **Problem:** each mutation below was switched on by an environment variable in the scratch build
  and the module's release tests run (`tmp/review_AB/unit/mut_*.log`). Each one changes behaviour;
  only the first is caught, and only by the helper's own unit test.

  | mutation | caught by |
  |---|---|
  | a start at the pass limit becomes a yardstick (`is_the_new_yardstick` ignores the ending) | `the_yardstick_is_the_best_converged_start` only; **`a_start_at_the_pass_limit_is_no_yardstick` passes** — at `max_passes = 15` no later start reaches a judging pass that could agree |
  | a negative decrement counts as converged | nothing — though the 3-sample census fixture has decrements of −36.4 and −37.6 at passes 20 and 27, a ready drawn case (step A2's plan asked for one "if found") |
  | step B1's addition removed (tries leave fewer than two passes, the cycle runs anyway) | nothing |
  | an accepted try's pass thrown away and re-run | nothing |
  | the first try always the whole step (`first_newton_fraction` bypassed at its call site) | nothing |

- **Why it matters:** the re-fixturing to `max_passes = 15` kept the yardstick test green by making
  it vacuous; the other four are behaviours the spec names (§2 traps; §3.3.1 points 2 and 3, and the
  pass budget) that only the bit-for-bit reproduction of the prototype pins today, and that check is
  not a test.
- **Suggested fix:** see section 8.

### Minor

- [src/parameter_estimation/joint/fit.rs:2801](../../../../src/parameter_estimation/joint/fit.rs#L2801) — **[Minor]** Mi1: the progress line prints the size of the change, not its sign
- **Confidence:** High
- **Problem:** `print_the_cycle` (new in step B1, carrying the old format) prints `gain.abs()` as
  "log-likelihood moved". A cycle that lost 50 units reads exactly like one that gained 50, which is
  how B2 stayed unnoticed in the run logs.
- **Suggested fix:** print the signed gain (`{gain:+.2e}`) and say "rose" or "fell".

- [src/parameter_estimation/joint/fit.rs:2413](../../../../src/parameter_estimation/joint/fit.rs#L2413) — **[Minor]** Mi2: an accepted try's pass can be left unjudged at the limit
- **Confidence:** High (by reading; not hit on any measured run)
- **Problem:** the loop needs three passes in hand to start a cycle, but a cycle that begins with an
  accepted try's pass has already run its first. When the limit falls there, the loop exits with that
  pass — a judging pass — never judged: a start whose step would have converged is recorded at the
  limit, and the final pass is re-run at the same parameters.
- **Suggested fix:** test `passes + 2 <= max_passes` when `accepted_pass` is `Some`, or judge the
  accepted pass before leaving the loop.

### Nits

- `JointFitConfig::decrement_threshold` ([fit.rs:686](../../../../src/parameter_estimation/joint/fit.rs#L686))
  does not say what zero or a negative value does; the tests rely on zero meaning "never converges".
- `parameters.noisy.swap(group, group)` ([fit.rs:4339](../../../../src/parameter_estimation/joint/fit.rs#L4339))
  swaps an element with itself (pre-existing; goes with B2's fix).

## The two defects, traced

### (a) Tiny inputs never converge after step A2

Two causes stack. **First, B1:** on the digest fixtures nearly every box solve reaches its round
limit, and the decrement is read from a step outside the box (290 most of the time). That alone keeps
every start at the limit. **Second, M1:** once the solve is correct, the decrement is honest but
small and slow — 0.27 at start 1's first judging pass on the fixture with errors, still 0.020 at
pass 196 —
because what is left to gain is along directions the data do not determine. Each candidate the lead
listed was built behind a switch in one scratch build and run on the same inputs:

| candidate | 2 samples, no errors | 2 samples, with errors | 21 samples | 1 oracle sample | 4 oracle samples | joint-fit tests |
|---|---|---|---|---|---|---|
| step B1 as committed | limit ×3 | limit ×3 | limit ×3 | reference | reference | 440 pass |
| failed solve → "cannot judge" | limit ×3 | limit ×3 | limit ×3 | unchanged (no solve fails there) | unchanged (no solve fails there) | not run |
| hold parameters without an error out of the step | 16 ×3 | 43; agreed 21, 15 | 95, 62, 59 | identical | 56, 49, 37 passes (moves) | **1 fails** |
| primal active-set only (B1's fix) | 37, 41, 41 | limit ×3 | limit ×3 | identical | same passes; ≤ 10⁻¹¹ | 440 pass |
| primal + determined-direction decrement | 16 ×3 | 42; agreed 19, 13 | limit ×3 | identical | = primal | 440 pass |
| **primal + determined decrement + no swap** (B1 + M1 + B2) | **16 ×3** | **22, 13, 13** | **84, 38, 44** | **identical** | **= primal** (no swap fires) | **440 pass** |

"What stops the start" when a solve fails: with the primal solve nothing fails on any input measured,
so the question does not arise; should it ever, the start runs to `max_passes` and says so, which is
the honest outcome.

What the recommended combination does to the digest files, against step B1's (to be re-recorded at
step E1 with these notes): no-errors fixture, the two SNPs' QUAL +9.1 and the other sample's GQ
74 → 70, no genotype; with-errors fixture, QUAL −5.2 and −12.6, no genotype; 21 samples, one record
replaced by another and 2 of 462 genotypes changed, QUAL −241 to +344. Step B1's numbers there come
from fits that never converged, so they are not the reference to defend.

### (b) "Every other cycle loses about 50 units" on the two-sample fixture with errors

**Where:** in the plain maximisation step, in the error rates, at the swap of B2. Not in SQUAREM's
jump, not in the Newton step, not in a bound or clamp, not in the reach schedule, not in the final
pass.

**How it was found.** A probe evaluated, whenever a plain pass lost more than one unit, the
log-likelihood with each block of the update applied alone (the three shares, the density's shapes,
the carrier shapes, the error rates, the homozygote excesses). Every one of the 37 losses on this
fixture (up to 99 units) and the 53 on the 21-sample fixture (up to 909 units) is reproduced by the
error rates alone; no other block, applied alone, lowered it at all. The swap fired 1,040 times on
this fixture (library 0: 515, library 2: 512, library 1: 13) and 6,126 times on the 21-sample one;
on the 1- and 4-sample oracle runs it never fired.

**What happens, in one cycle.** On two samples the mismapped class holds about 0.3% of positions. Its
rate for library 0, fitted from those few positions, comes out at 2.64 × 10⁻⁴, below the ordinary
rate of 2.44 × 10⁻³. The swap then gives the ordinary positions 2.64 × 10⁻⁴ — far too clean for
them — and the log-likelihood falls from −482.38 to −532.0. On the next pass the ordinary rate is
re-fitted to 2.44 × 10⁻³ and the log-likelihood returns. SQUAREM's cycle is a plain pass, another
plain pass and a jump; with the jump at its shortest (α = 1) it accepts the second plain pass's point
without a check, so the fit lands on the bad half of the oscillation every other cycle.

**Older than this branch.** At main 51b7be99, with the old settled test switched off so the start
runs past pass 49, the log-likelihood at each cycle's start alternates −482.4 and −532.0 from pass 53
to the limit — the same numbers step A2's note recorded. The old test stopped the start at pass 49,
before the oscillation set in.

**Fix and its measurement.** Deleting the swap: no loss of more than a unit on any digest fixture;
on the fixture with errors the cycle-start log-likelihood climbs steadily to −482.3 instead of
alternating; oracle untouched (the swap never fires there). **Rejected alternative, measured:**
replacing the swap by the order-constrained update (both classes' read counts pooled into one rate
when the order is violated). It removes the losses on the fixture with errors but leaves 12 on the
21-sample fixture, up to 35 units: the Newton step can move a library's two rates out of order, and
from such a point the constrained update is no longer guaranteed to climb.

## 7. Out of scope observations

- The homozygote-excess update builds its quadrature from the density's *new* shapes but reads
  genotype counts accumulated at the *old* shapes' nodes
  ([fit.rs, `maximisation`](../../../../src/parameter_estimation/joint/fit.rs#L4346)). Applied alone it
  never lowered the log-likelihood in the measurements above, but it is not an exact conditional
  maximisation. Suggested follow-up: build the excess's quadrature from the shapes the
  counts were taken at, or note why the difference does not matter.

## 8. Missing tests to add now

`try_the_newton_step` and `maximise`:

- `a_start_whose_decrement_is_negative_is_not_stopped` — the 3-sample census cohort of
  `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` at 60 passes judges decrements of
  −36.4 and −37.6 at passes 20 and 27 with today's solve: assert from the trace that the start did
  not stop at either. (With B1's fix no decrement is negative; the test then becomes "every traced
  decrement ≥ 0".)
- `refused_tries_that_leave_fewer_than_two_passes_end_the_start_at_the_judged_point` — a drawn cohort
  and a `max_passes` placed so that four refused tries end two passes short; assert
  `passes ≤ max_passes` and the returned parameters equal the judged ones.
- `an_accepted_newton_step_costs_one_pass` — with the trace on, assert that a judging pass directly
  after an accepted try has the number of the try's pass (record the Newton outcome in
  `PassSummary`, which also lets the census test assert it exercised both acceptance and refusal —
  measured today: start 1 accepts at passes 20, 36 and 40, but nothing asserts it).
- `the_first_newton_try_after_an_accepted_quarter_is_half_the_step` at the `maximise` level, through
  the same `PassSummary` field.

`a_start_at_the_pass_limit_is_no_yardstick`: choose a fixture and limit where a later start reaches
two or more judging passes before the limit, and assert that with the yardstick rule mutated it
*would* agree (run once by hand to choose the fixture).

`bounded_step` / its replacement:

- `the_box_step_stays_inside_the_box_and_is_worth_at_least_nothing` — a property test over random
  positive-definite and near-singular matrices (condition numbers to 10¹²), random slopes and random
  positions inside or on the bounds: every returned step lands inside the box, the decrement is ≥ 0,
  and the solve reports that it finished. Today's solve fails it.

## 9. What's good

- `try_the_newton_step` takes the pass as a closure ([fit.rs](../../../../src/parameter_estimation/joint/fit.rs#L2725)),
  so the halving rules are unit-tested with stubbed log-likelihoods, including the tie refused and
  the bound halved without a pass.
- `has_converged` is a two-line pure function whose test names the spec's exact cases (−2,959, NaN,
  0, just below, at the threshold) ([settled.rs](../../../../src/parameter_estimation/joint/fit/settled.rs)).
- The decrement tests compare against an independent dense solve on both the whole matrix and the
  arrow, with and without a bound binding ([standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs#L2569-L2658)).
- Every re-recorded checksum carries a note saying what moved and why, which is what made defect (b)
  findable from the commit message alone.

**Claims in the diff, re-derived:** 20 samples converge in 14 passes, 0.0058 errors from 600 passes
— CHECKED-CORRECT; 4 samples in 24 passes, 0.0075 — CHECKED-CORRECT; the second-shape-at-50 case in
46 passes — CHECKED-CORRECT; the 1- and 4-sample reproductions — CHECKED-CORRECT (files
byte-identical); "no parameter has an error" on the two-sample fixture without errors —
CHECKED-CORRECT (16 of 16 without an error at every judging pass); "the decrement stays near 290
because the bounded solve ends at its round limit" — mechanism CHECKED-CORRECT (164 of 180 solves at
the limit at step B1); "where in the cycle they are lost was not traced" — now traced, (b) above.
Not re-derived: `SHORT_OF_CONVERGING`'s "25, 27 and 24 passes", step A2's "step of 2.7 × 10⁶ to the
density's second shape" (step A2's trajectory, not re-run).

**Determinism:** the new code adds no parallel sum; the target, the tries and the decrement are
sequential arithmetic on sums the pass already makes in its fixed chunk order. The census test runs
1, 4 and 8 threads over a 60-pass fit that accepts and refuses tries (measured above) and passes;
the digests' pool-width test passes.

## 10. Commands to re-verify

- `cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --release --lib parameter_estimation::joint`
- `cargo test --release --lib cli::cross_platform_digests`
- With the prototype applied (`tmp/review_AB/prototype_with_probes.diff`, every change behind a
  switch): `PROBE_BOX=det PROBE_RATES=noswap` on the commands above and on
  `tmp/review_AB/run.sh <name> {s1|o4}`; `PROBE_FIT=1` prints the probes; `PROBE_MSTEP=1` the
  per-block loss attribution; `PROBE_MUT=<yardstick|negative|nobreak|rerun|whole>` the mutations.

*`PROJECT_STATUS.md` was not updated; the lead commits and owns it.*

---

## Addendum (2026-10-17): does deleting the swap switch the classes' labels?

The lead asked, before taking B2 to the owner, whether the swap is what keeps the mismapped class
meaning "the higher error rate". For each input, the fit was run twice with B1's and M1's fixes in
place: once keeping the swap and once without it. On the 59-sample input it was also run with the
code as committed. The returned values were read from the fit trace (`PVC_JOINT_FIT_TRACE`).
"Inverted" means a library whose mismapped-position rate ended below its ordinary-position rate.

| input | swap fired (all starts, every pass) | libraries ending inverted, swap kept / swap deleted | mismapped share, kept / deleted | the two files |
|---|---|---|---|---|
| 1 oracle sample | 0 | 0 of 1 / 0 of 1 | 0.00879 / 0.00879 | byte-identical |
| 4 oracle samples | 0 | 0 of 4 / 0 of 4 | 0.00518 / 0.00518 | byte-identical |
| 59 tomato2 samples, one chunk in four | 0 (also 0 in the committed code) | 0 of 59 / 0 of 59 | 0.00777 / 0.00777 | byte-identical |
| 59 tomato2 samples, all chunks (committed code only, 993 s at 18 threads) | **0** | 0 of 59 | 0.0106 | — |
| digest fixture, 2 samples, no errors | 0 | 0 of 3 / 0 of 3 | 0.00019 / 0.00019 | identical fit |
| digest fixture, 2 samples, with errors | 133 / — | **3 of 3** / **2 of 3** | 0.0000066 / 0.0026 | differ |
| digest fixture, 21 samples | 5,532 / — | **10 of 22** / **12 of 22** | 0.40 / 0.45 | differ |

1. **Can labels switch?** Not on any real input: the swap never fired on 1, 4 or 59 tomato samples.
   On the full 59 samples (committed code), the closest library's mismapped rate is 2.66 times its
   ordinary rate (median 12.6, largest 24.7); on 4 samples, 7.7 times; on 1 sample, 46.5 times. On
   the synthetic digest fixtures, labels switch **with or without the swap**. With the new
   convergence rule the returned point is the one a judging pass judged. That point was reached by a
   Newton step or a SQUAREM jump, neither of which applies the swap. So keeping the swap no longer
   guarantees the order of what is returned. In the code as committed the order held (0 inverted on
   all three fixtures), but only because every start ended at the pass limit on an M-step output.
2. **Mismapped share and rates, swap kept against swap deleted:** identical to the bit on every
   real input, because the swap never fires there. On the digest fixtures the two runs differ.
   - With errors, swap kept: share 6.6 × 10⁻⁶. Every library's mismapped rate ends at its floor of
     10⁻⁴, below ordinary rates of 2.0 to 2.5 × 10⁻³.
   - With errors, swap deleted: share 0.0026. Library 1's mismapped rate is 0.237; libraries 0 and
     2 sit at the 10⁻⁴ floor. Best log-likelihood −482.27 against −490.08 with the swap kept.
   - 21 samples: on the swap-deleted run the "mismapped" class holds 45% of positions (40% with the
     swap kept, 44% in the committed code). On this fixture the two classes do not describe
     mismapping at all: the fit is splitting ordinary positions into two error levels.
3. **Did the swap ever fire on the 59-sample input in the current code?** No. It fired zero times in
   the committed code on the full 59 samples and on one chunk in four, across every pass of all three
   starts. So deleting it cannot change that input's output.
4. **What correct alternative exists, if the order matters.**
   - *Relabelling the whole class* (the share becomes one minus itself, and every library's two
     rates exchange) is an exact symmetry of this model: the two classes share the density and the
     duplicated class and differ only in share and rates. But it cannot be applied here. The
     mismapped share is bounded to at most 0.5, so the bound already defines the mismapped class as
     the minority one, and the relabelled share would fall outside it. Nor does any relabelling fix
     the 21-sample fixture, where libraries disagree: 10 or 12 of 22 are inverted, not all.
   - *The order as a constraint everywhere* is the correct way to impose it. When the order is
     violated, the update fits one rate from both classes' pooled read counts. In addition, any point
     a SQUAREM jump or a Newton step proposes is put back in order (a library with inverted rates
     takes their mean for both). This was prototyped (`PROBE_RATES=pool_proj`) and measured:
     - no pass loses more than a unit, and no library ends inverted;
     - it costs likelihood where the data want the inversion. With errors: converged in 17, 16 and 13
       passes, best −482.39 (−482.27 without the order). 21 samples: start 1 at the pass limit at
       −5095.1, against −5045.75 without the order (starts 2 and 3 converge in 38 and 41 passes);
     - the earlier pooled-only variant failed precisely because jumps and Newton steps left the
       ordered region.
   - **Recommendation unchanged:** delete the swap and let the bound on the share define the
     mismapped class. Nothing outside the fit reads a library's mismapped rate. If a reader of the
     parameters file should be told, add a note naming any library whose mismapped rate ended below
     its ordinary one. That never happened on any real input measured. Impose the order as a
     constraint (pooled update plus the jump and Newton points put back in order) only if a consumer
     turns out to need it.

**Side observation on B1's fix, at 59 samples (one chunk in four).** The primal box solve leaves pass
counts unchanged (45 converged; 46 and 37 agreed) and moves SNP/indel numbers by at most 3 parts in a
million (contamination shares). One repeat-tract number moves a lot: a slippage group's chosen
smoothing strength goes from 5.1 × 10⁹ to its 10¹⁰ cap. That is the penalised repeat-tract fit
choosing between two nearly equal options, not the SNP/indel fit. Step E1 should expect it.

Files: runs `tmp/review_AB/runs/{s1,o4,t59q}_{swap,noswap}.*`, `t59q_committed.*`,
`t59_committed.*`; digest runs `tmp/review_AB/digests/labels/`; summaries `tmp/review_AB/labels.py`
and `ratio.py`.
