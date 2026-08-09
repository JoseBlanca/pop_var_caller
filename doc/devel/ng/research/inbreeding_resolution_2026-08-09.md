# What the runs model's noise floor is a function of — and the one diagnostic that sees it

**Date:** 2026-08-09. **Harness:**
[`examples/ng_inbreeding_resolution.rs`](../../../../examples/ng_inbreeding_resolution.rs),
`cargo run --release --example ng_inbreeding_resolution`. Twelve grid points — three genome
sizes crossed with four levels of evidence per window — five seeds each, three questions.

**The question.** Milestone E3 shipped `resolution_at(windows)`, an interpolation of research
note §3.6's four measured points, and `MAX_IDENTIFIED_STATE_RATIO = 0.9`, calibrated on ten
fitted ratios from one fixture shape. Both are constants standing in for a measurement that
had never been made across shapes. §3.6 varied the **window count** and nothing else, but a
two-state chain classifies a window on how many heterozygotes it holds — `sites ×
heterozygosity`, whose sampling noise is `√h` on a mean of `h`. §3.6's genomes carried
100,000 sites a window at one heterozygote per kilobase: **100 heterozygotes a window**. A
shallow sample or a region-restricted run carries far fewer, and the shipped floor has never
been checked there.

---

## 1. The floor is not a floor. It is a rare catastrophic mode

Genomes drawn with **no runs at all**, so `F` must come back at nothing.

| windows | het/window | refused | mean `F` | worst `F` | shipped | worst ÷ shipped | starts spread |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3,000 | 2 | 1 | 0.0510 | 0.1669 | 0.0409 | 4.08 | 0.8379 |
| 3,000 | **5** | 0 | 0.2008 | **0.9912** | 0.0409 | **24.3** | 0.9979 |
| 3,000 | 20 | 1 | 0.0075 | 0.0152 | 0.0409 | 0.37 | 0.1249 |
| 3,000 | 100 | 2 | 0.0034 | 0.0063 | 0.0409 | 0.15 | 0.0063 |
| 6,000 | 2 | 0 | 0.0145 | 0.0342 | 0.0135 | 2.54 | 0.6615 |
| 6,000 | 5 | 0 | 0.0020 | 0.0048 | 0.0135 | 0.36 | 0.4037 |
| 6,000 | 20 | 1 | 0.0193 | 0.0647 | 0.0135 | 4.81 | 0.9980 |
| 6,000 | 100 | 2 | 0.0006 | 0.0008 | 0.0135 | 0.06 | 0.0138 |
| 12,000 | 2 | 2 | 0.0083 | 0.0145 | 0.0065 | 2.22 | 0.7054 |
| 12,000 | **5** | 2 | 0.3888 | **0.9922** | 0.0065 | **152** | 0.9919 |
| 12,000 | 20 | 1 | 0.0019 | 0.0045 | 0.0065 | 0.69 | 0.1304 |
| 12,000 | 100 | 2 | 0.0009 | 0.0018 | 0.0065 | 0.27 | 0.0092 |

**Three findings.**

- **At 20 and 100 heterozygotes a window the shipped floor is conservative** — worst `F` is
  0.06 to 0.69 of it. That is §3.6's own regime, and tomato's at 100 kb windows. Nothing is
  wrong there.
- **Below that it is wrong, and not by a factor.** Two of sixty fits returned **`F` ≈ 0.99 on
  a genome with no runs at all** — 24 and 152 times the reported floor.
- **It is not a smooth function of either variable.** 3,000 × 5 fails and 6,000 × 5 does not;
  3,000 × 2 gives 0.167 and 12,000 × 2 gives 0.015. Most seeds behave and the occasional one
  is catastrophic. **So no two-argument `resolution_at` can express it** — the thing to
  detect is a failure mode, not a magnitude.

## 2. The state-ratio check cannot see it, because the fit manufactures a separation

The fitted inside heterozygote rate over the fitted outside one, reported whether or not the
fit was accepted. Genomes with no runs have no second state to find; genomes 30% covered by
runs have one, drawn at 0.04.

| windows | het/window | no runs: mean (worst) | with runs: mean (worst) |
|---:|---:|---:|---:|
| 3,000 | 2 | 0.546 (0.818) | 0.040 (0.046) |
| 3,000 | 5 | 0.313 (0.788) | 0.041 (0.044) |
| 3,000 | 20 | 0.692 (0.841) | 0.041 (0.042) |
| 3,000 | 100 | 0.806 (0.837) | 0.040 (0.041) |
| 6,000 | 5 | 0.209 (0.453) | 0.039 (0.040) |
| 12,000 | 5 | 0.620 (0.893) | 0.040 (0.040) |

On the two catastrophic cells the no-runs ratio is **0.31 and 0.62 on average** — far below
the 0.9 threshold, so the fit is accepted. **The chain invents a state separation out of
sampling noise and it looks entirely real.** `MAX_IDENTIFIED_STATE_RATIO` therefore catches
the *coincident-states* failure it was measured on and not this one; they are different
failures with the same symptom.

## 3. Where runs actually exist, the estimator is excellent at every shape

Genomes drawn 30% covered by runs.

| windows | het/window | refused | realised `F` | mean fitted | worst error | starts spread |
|---:|---:|---:|---:|---:|---:|---:|
| 3,000 | 2 | 0 | 0.3035 | 0.3019 | 0.0053 | **0.0000** |
| 3,000 | 5 | 0 | 0.2789 | 0.2794 | 0.0016 | **0.0000** |
| 3,000 | 20 | 0 | 0.2700 | 0.2700 | 0.0001 | **0.0000** |
| 3,000 | 100 | 0 | 0.2985 | 0.2985 | 0.0000 | **0.0000** |
| 6,000–12,000 | 2–100 | 0 | 0.298–0.310 | same | ≤ 0.0020 | **0.0000** |

Sixty fits, **zero refusals**, worst recovery error **0.0053** and zero to four decimals at
100 heterozygotes a window. This is not a generally unreliable estimator. It is specifically
an *absent* signal being read as a total one.

## 4. The diagnostic that separates them is already computed and thrown away

Compare the last column of §1 with the last column of §3.

| | across-start spread |
|---|---|
| genomes **with** runs, all 12 shapes, 60 fits | **0.0000** — every one of the nine starts lands on the same `F` to four decimals |
| genomes **with no** runs | **0.0063 to 0.9980**, and **0.99 on both catastrophic cells** |

**Spec §6.5 already names this**: *"A run where every start returned the same `F` at the same
score has not measured zero autozygosity; it has failed to find anything, and must say so."*
The inverse is what the table shows and what nobody had measured: **when the starts agree,
the answer is real.** `RunsModelFit::starts_tried` carries every start's `F` and E3 reports it
without ever reading it.

**A threshold near 0.05 separates the two populations at every shape measured.** It would
accept all sixty with-runs fits and refuse every catastrophic no-runs fit, along with most
benign no-runs fits — which is the design's own preference for this parameter (*fail rather
than emit*), since on a genome with no runs `F` is not identified and the honest answer is
that nothing was found.

## 5. What this does not settle

- **The threshold is measured on twelve shapes at five seeds.** The two populations are
  separated by two orders of magnitude here, so it is not a tuned constant — but the
  with-runs column is *exactly* 0.0000, which is suspicious enough to be worth confirming on
  a shape whose runs are short or sparse rather than 30% at 3 Mb.
- **No real data.** Every genome here is drawn from the model the fit assumes.
- **`MIN_WINDOWS_TO_FIT_INBREEDING` is untouched by this.** If the across-start spread becomes
  the criterion, a window-count floor may no longer be the right gate at all — 3,000 windows
  at 100 heterozygotes each behaves perfectly and 12,000 at five does not.
