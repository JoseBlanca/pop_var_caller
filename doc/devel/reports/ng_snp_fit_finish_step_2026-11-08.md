# The SNP/indel fit's Newton step: length from the curvature, blended with the last step

**Date:** 2026-11-07 to 2026-11-08. **Branch:** `fit-finish`. **Spec:**
[`parameter_prepass_joint_fit.md`](../ng/spec/parameter_prepass_joint_fit.md) §3.3.1, amended.

## 1. The answer

On 400 tomato accessions at about 3 reads a position, the SNP/indel fit took 22 minutes instead of
66. Its finish on every position needed 40 passes instead of 132, and the whole
`estimate-parameters` run took 36 minutes instead of 80. Two changes to the Newton step did it:

- **The length of each try comes from the log-likelihood's curvature along the step.** That
  curvature is measured from the try just made, but the next first try is never more than twice
  the fraction last kept.
- **Each step is blended with the last kept one** (preconditioned nonlinear conjugate gradients).

The fit reaches the same answer everywhere it was measured. On small cohorts, where the finish was
already short, it costs a few passes:

| input (`estimate-parameters`) | before: passes, fit time | after: passes, fit time |
|---|---|---|
| kimura, 400 samples, from scratch, 48 threads | 50 + 132, 66 min | 48 + 40, 22 min 10 s |
| tomato2, 59 samples, 6 threads | 44/39/42 + 15, 5 min 36 s | 33/40/44 + 14, 5 min 20 s |
| one tomato2 accession | 59/23/19 + 8, 11 s | 65/21/19 + 10, 13 s |
| HG002 chr20, 5×, one sample | 23/16/14 + 12, 9 s | 26/15/15 + 18, 13 s |
| HG002 chr20, 30×, one sample | 14/10/11 + 7, 6 s | the same passes and numbers |

The passes are those of the starts on the thinned census, then of the winner's finish on every
position (§3.3.2). Every fitted number that has a standard error moves by less than 0.06 of it:

- 0.034 at most on 59 samples;
- 0.06 on the one accession;
- 0.006 on HG002 at 5×.

At 400 samples every number agrees within 0.015 standard errors with the run with the blend and no
cap, and that run agreed within 0.06 with the old rule.

A re-fit after adding 20 samples to a 380-sample fit took 15 minutes instead of 29. That was
measured with an earlier form of the rule (no cap, a sixteenth floor), so it is indicative only.

## 2. Why the old finish was slow

In the old rule a judging pass tried its Newton step first at twice the fraction last kept, and
halved a try that fell. On 400 samples, 64 of the finish's 133 passes were such falling tries,
thrown away. Their drops were not arithmetic noise: summing the log-likelihood over 2 million
positions is exact to about 10⁻⁵ units, and the tries fell by far more. The doubled first try
overshot again and again, because the information summed from per-position slopes misstates the
curvature (`fit_precision.md` §3.2).

The kept steps also undid one another. Along the mismapped share and its error rates, each kept
step closed only a fixed share of the remaining gap. The blend carries the last direction forward,
so the next step does not zig-zag back.

## 3. What was built

All of it is in `maximise` and `try_the_newton_step` (`src/parameter_estimation/joint/fit.rs`).
The judging pass now also hands back each parameter's slope (`NewtonStep::slopes`), which it
already summed.

1. **Length.** A try at fraction `t` that gained `Δ` fixes the parabola `s·t − ½c·t²` through the
   judged point. Here `s` is the slope along the whole step, and `c = 2(s·t − Δ)/t²`. The parabola
   peaks at `t* = s/c` (`peak_along_the_step`). A refused try is followed by one at `t*`, held
   between `t/4` and `t/2`. A kept try hands `t*` on as the next first try, capped at `2t`.
2. **Blend.** `p = d + β·p_prev`, with `β = max(0, g·(d − d_prev)/(g_prev·d_prev))`
   (Polak–Ribière), and the plain Newton step `d` whenever `p` would not point uphill
   (`PreviousStep`). After a step no try kept, the blend restarts.

Nothing new is summed. A try still costs one pass and is kept only if it rises, and the result
still does not depend on the thread count.

## 4. Measured and not adopted

The same inputs, under rules that differ from the adopted one:

| rule | 400 samples: finish, fit | 59 samples | one accession | HG002 5× finish |
|---|---|---|---|---|
| **adopted:** quarter floor, capped, blended | 40, 22 min 10 s | 5 min 20 s | 65/21/19 + 10 | 18 |
| sixteenth floor, capped, blended | 47, 26 min 27 s | 6 min 21 s | 55/21/19 + 10 | 17 |
| sixteenth floor, blended, **not capped** | 53, 29 min 30 s | **9 min 06 s** | 72/21/19 + 10 | 15 |
| sixteenth floor, capped, no blend | 65, 34 min 21 s | 5 min 42 s | 51/20/20 + 9 | 12 |
| sixteenth floor, not capped, no blend | 83, 47 min* | 6 min 15 s | 51/20/19 + 9 | 12 |

\*That run shared the machine with another for part of its time, so its time is high; its pass
count is not affected.

**The cap is needed.** Without it, a short kept try sees too little of the curvature and puts the
peak far beyond itself. The next try overshoots, is cut short, is kept, and the cycle repeats. On 59
samples one thinned start then took 110 passes. On the test suite's four-sample fixture, the start
converged after 46 passes, against 25 before.

**The blend is most of the gain at 400 samples.** At the sixteenth floor, the finish took 47
passes with the blend and 65 without. Below about 60 samples it costs a few passes.

## 5. What the range covers, and what it does not

The inputs measured are:

- one sample at 5× and at 30× (HG002 chr20);
- one sample at 3× (a tomato2 accession);
- 59 samples at 3× (tomato2);
- 400 samples at 3× (kimura).

Not measured: several hundred samples at high depth, and thousands of samples. Both are where the
finish is longest, so they are where the change should help most, but that is an expectation, not
a measurement.

## 6. Tests

- **Rewritten**: the two tests that pinned the old schedule. Now the first try is the whole step,
  then at most twice the fraction kept and at least 1/1024, and unchanged after a refused step. A
  kept try's pass is still judged at the pass limit.
- **New**: the parabola's peak, the shortening of a refused try, and the blend.
- **Redrawn**: the agreeing-later-start fixture, now seed `+ 1002`. Its old seed no longer had a
  second start that would win alone.
- **Restated**: one bound. The shapes-at-rest test required the remaining step in the first density
  shape to be below a hundredth of that shape. The shape has no standard error, and where it stops
  along its ridge depends on the path:
  - main stopped at a shape of 0.232 with 0.00217 left, 94% of that bound;
  - the adopted rule stopped at 0.151 with 0.00164 left, nearer rest in absolute terms, yet past a
    hundredth.

  The bound is now an absolute 0.005, still 58 times below the 0.29 the defect it guards against
  left.
- **Re-recorded**: the 21-sample checksums. The fit converges after 42 passes instead of 49, no
  fitted number with a standard error moves by more than 0.057 of it, and the calls keep the same
  records and genotypes. QUAL moves by at most 4.5, and a GQ by at most 1.

## 7. A side finding, not addressed here

Every fit from scratch on the 400 samples stops at the same log-likelihood, −254,339,651.89.
Restarts from a finished answer climb higher, to −649.11 from the 380-sample fit, so the stop is
about 2.8 units short of the maximum.

The likely cause is the duplicated share, which sits on its upper bound of 0.05. There the
carrier-frequency shapes get no standard error. The stopping rule counts only directions that have
one (`fit_precision.md` §2), so it declares convergence while those shapes could still gain about
3 units. The shapes-at-rest test in §6 touches the same gap.

The answer moves only a little. The carrier shapes move by about 2% (0.121 → 0.119, 7.45 → 7.35),
and the frequency shape `b` by 0.2 standard errors. The error rates at mismapped positions move by
up to 1.3 standard errors in one read group. The data seem to want more than 5% of positions
duplicated, which may deserve its own look.

## 8. Where things are

- **Code and tests:** branch `fit-finish`.
- **Experiment code with environment switches:** branch `fit-finish-kimura` (from tag
  `fit-subset-experiment`), kept as tag `fit-finish-experiment`.
- **Raw runs:**
  - kimura runs (logs, parameters files, per-pass traces and scripts): copied back to
    `tmp/fit_finish/kimura_runs/` in the main checkout; the kimura checkout is deleted;
  - laptop runs (logs, parameters files and scripts): `tmp/fit_finish/check/` in the main checkout.
