# E6 correctness review: branch `snp-fit-convergence`, HEAD b73be771 (2026-10-17)

*Written by the reviewing agent; filed by the lead (the harness refused the agent's file write).*

**Verdict: approve with changes.** The stopping rule, the Newton step, the box solve, the bound fix
10d2a24c and the determined decrement fit together correctly; the full suite passes; the cross-platform
checksums agree on macOS and in the Linux container at HEAD. Two things need attention before merging:
M1 (a defect reaching only the report and the parameters file; calling reads neither) and M2 (the merge
with main).

**Scope:** `git diff main...HEAD` (merge base 51b7be99), 20 commits, 36 files; findings already fixed
from the AB, B and CD reviews not repeated.

**Commands:** release suite 5,288 passed, 3 failed (the debug-only panic tests), 7 ignored; clippy
`-D warnings` and fmt pass; `cargo doc` fails only on broken links already on the base (none in a
touched file); `cross_platform_digests` in the Linux container: 5 pass at HEAD.

## Major

### M1 — `determined.rs:339-351, 356-389, 417`: impossible errors and a misleading verdict at two samples

A drawn cohort (2 samples, 3 reads, 10,000 positions, duplicated class on, seed 77): the fit returns no
standard error for any of the density's four numbers, yet the report says the share above 1% is
0.073 ± 11.0, Hexp 0.031 ± 4.8, E[f] 0.037 ± 5.8, and `separately_determined = true` with trade cost +∞.
Validation accepts it; the run log calls the pair separately determined. Cause: a number whose error is
wider than its range keeps its covariance (as the spec intends), but the derived numbers go through
`determined::error_from_variance`, which has no range rule, and σ and S move together although neither
is determined. At one sample the group is dropped as not identified, so nothing goes wrong; the band is
where numbers are kept but their errors exceed their ranges. **Fix:** apply the range rule to S, Hexp and
E[f] (each in [0, 1]); in `DensityReport::of` and `DuplicatedReport::of`, mark the pair when neither (or
any) number of its group has an error, recorded in spec §2.1.3; a test on this fixture.

### M2 — the branch does not merge cleanly into main

`git merge-tree main HEAD` conflicts in `scripts/promote_ng_oracle.baseline`, `src/cli/calling_run.rs`
(main removed `beside_the_vcf`) and `src/cli/cross_platform_digests.rs` (4 checksum hunks). Main's
763b1066 ("a calling run records its parameters in its log, not beside its VCF") changes the path D2
relies on to carry the two tables back. **Fix:** merge main into the branch; check the parameters block a
calling run writes into its log carries `fitted_density` and `duplicated_positions` (spec §2.1.3); re-record
the oracle and the five checksums, natively and in the container.

## Minor

- **Mi1 — `standard_errors.rs:1060-1155`:** the box solve's cost grows with the square of the samples
  when many parameters cross an end at different fractions (each round fixes one). One library a sample,
  each excess pushed past 0 at a different fraction: 10 ms at 59 samples, 5.4–5.9 s at 2,000. **Fix:**
  fix every parameter within a tolerance of the shortest room in one round, or a gradient-projection
  solve (Moré–Toraldo); minimum: measure a judging pass's share of wall time on a real 2,000-sample run.
- **Mi2:** the E3 defect has no whole-fit test. A drawn cohort (10 or 20 samples, 20,000 positions, 3
  reads, excess drawn at 1.0, under a second) reproduces it with the bound fix off (8 of 10, 15 of 20
  samples lose their coefficient) and passes at HEAD. Add
  `fit_jointly_returns_a_coefficient_for_every_sample_when_excesses_sit_at_one`.
- **Mi3:** docs describing removed behaviour: `standard_errors.rs:10-13` (the old stopping rule);
  `settled.rs:20-21` ("one number bounds every parameter" — since B4 only determined ones);
  `information.rs:5274-5277` (the ignored measurement's thresholds measured under the old rule).
- **Mi4:** specs lag the code: §3.3.1 step 1 still says the Newton target is "clamped into its interval"
  (since 10d2a24c the golden-section numbers stop half its width inside); §2.1.3's formula says
  "V_trade ≤ 0: no trade" where the code uses ≤ 10⁻¹²·Var(σ).
- **Mi5 — `validate.rs:1233-1240`:** a constant sits between `refuse`'s doc comment and the function.

## Nits

- A whole-step Newton try `v + 1·(t − v)` can land one ulp outside a bound (13% of random starts for a
  target of 0.02, 80% for 10⁻¹²): take the target exactly at fraction 1, or clamp each try.
- `information.rs:4087-4088` and `3983-3985`: two test docs misstate counts (8 accepted steps; 199 passes).
- Two redundant explicit link targets: `determined.rs:137`, `standard_errors.rs:854`.

## Across the range

| input | starts | report |
|---|---|---|
| 1 sample, 3 reads | converged at 17–18; others agreed | both pairs marked, no cost; excess held at 0 |
| 1 sample, 30 / 300 reads | converged at 12 / 8–10 | as above |
| 2 samples, 3 reads | converged at 42 | M1 when the duplicated class is on |
| 4 samples, 300 reads | converged at 21–25 | sensible verdicts (costs 1.41 and 3.90) |

At 2,000 samples a judging pass costs 160 µs for the determined decrement and 1.3–1.9 ms for the box
solve with no parameter near an end, up to about 6 s in Mi1's case.

Scratch: `tmp/review_E6/` (exports, probes `probe_e6_two_samples_report`, `probe_e6_box_rounds`,
`probe_e6_inbred_excess`, logs).
