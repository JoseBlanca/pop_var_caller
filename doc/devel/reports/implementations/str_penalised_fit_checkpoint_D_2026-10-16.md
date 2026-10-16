# Checkpoint D: the penalised repeat-tract fit, checked (2026-10-16)

Plan [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md), milestone D, on branch
`str-penalised-fit` (head `c4a425a6`, not pushed, not merged). Every D step is done; each has its own
report, linked below. This page collects them and the decisions they need before milestone E (deleting
the old fit, called **the climb** below) starts.

## 1. Decisions for the owner

1. **Do D2 and D5 pass, so that the climb can be deleted?** *Recommended: yes.*
   - On real reads (59 tomato2 accessions, about 3 reads a position), the new fit predicts held-out
     loci better than the counting fit, by 0.136 ± 0.053 log-likelihood units a locus over 2,349 loci
     ([D2](str_penalised_fit_D2_2026-10-16.md)).
   - In calling (HG002's tandem-repeat truth and the GIAB trio, at 30× and 50×), the new fit's
     parameters call as the climb's did: right genotypes within 10 tracts of about 2,850 or 3,850 in
     every class and depth, SNP calls identical, one true indel fewer at each depth
     ([D5](str_penalised_fit_D5_2026-10-16.md)). On one sample the new fit took 46 s where the climb
     took 49 minutes.
2. **Give the repeat-tract errors the allele-length prior's curvature** (the thin-data defect).
   *Recommended: yes, as a bug fix, before milestone E.* On thin data the fit had been losing every
   stutter error of a motif length, thick classes included: 89 of D1's 480 coverage runs. The cause is
   measured: the matrix the errors are read from leaves out the curvature of the prior on the
   allele-length shares, so at a thin class it is not invertible
   ([diagnosis](str_penalised_fit_thin_errors_diagnosis_2026-10-15.md)). Adding that curvature (about
   20 lines) inverts the matrix in all 384 runs checked, and returns every lost error. The returned
   errors are about as honest as the ones already reported. Fitted values do not change; reported
   errors do, so the oracle's fitted file would be re-recorded again.
3. **Where errors start: `LOCI_FOR_AN_ERROR` and `LOCI_FOR_A_SHARE_ERROR`.** *Recommended* by D1
   ([report](str_penalised_fit_D1_2026-10-15.md) §1):
   - the three stutter numbers get an error from 3 loci with reads instead of 8;
   - F and the allele-length shares keep 8;
   - a share's rule drops from 3 expected loci to 1.

   On classes of 3 to 7 loci, the shorter share and fall-off errors hold the truth within two errors
   90 to 92 times in 100, the slip rate 84 to 86 (thick classes: 92), and F only 72 to 80.
4. **The SNP/indel fit is where `estimate-parameters` spends its time.** *Recommended: a separate job
   after this branch merges, starting with taking the Newton step it already computes.* On 59 tomato
   samples the command took 61 minutes, 60 of them the SNP/indel fit, whose three starting points each
   stop at the 200-pass limit without converging
   ([performance review](../reviews/perf_estimate-parameters_2026-10-15.md), finding L5). That fit is
   outside this branch's scope.

## 2. What milestone D found

| step | what was checked | result |
|---|---|---|
| [D1](str_penalised_fit_D1_2026-10-15.md) | the bench's held-out cells, thin-class error coverage | the new fit predicts better than the counting fit in all 48 cells (+0.009 to +0.162 a locus), 167 of 168 seed-runs |
| [D2](str_penalised_fit_D2_2026-10-16.md) | tomato2, real reads, held out by whole locus | better by 0.136 ± 0.053 a locus; the margin is in the two- and three-base classes |
| [D3](str_penalised_fit_D3_2026-10-16.md) | kimura's size, drawn: 2,169 samples, 66.2 million pairs | 83 to 173 s at 18 threads, 151 to 252 s at 6; peak 10.5 GiB; budget one hour |
| [D4](str_penalised_fit_D4_2026-10-16.md) | the oracle and cross-platform digests | every move traced to a commit and re-recorded; full suite passes but for the 3 debug-only panic tests |
| [D5](str_penalised_fit_D5_2026-10-16.md) | calling on HG002 and the GIAB trio, new fit against the climb | the same to within 10 tracts and one indel |
| [D6](str_penalised_fit_D6_2026-10-16.md) | does the slip rate's 0.5% low bias move calls | about 4 calls in 100,000, in both directions, against 5 to 6 in 100 already wrong |

**One effect D4 found that its commit (C8) did not describe.** Six three-base classes of 1 to 5 loci
on the oracle lost their own F: once the stutter prior raised their slip rates, their F estimate ran
to the end of its search, and the existing rule writes no F then, so calling uses the cohort's
fallback F (0.613) at those tracts. It is the designed rule; whether a class should keep an F at the
end of its search would be a question for spec §2.10.

## 3. Performance, tonight

`pvc estimate-parameters`, the repeat-tract fit. Every change keeps every fitted number bit-identical:
checked on the reference fixture (153 runs), end to end on the oracle's four samples, and by a review
that ran the old and new code side by side on 180 combinations of edge inputs and thread counts.

| commit | change | measured on the bench, 2,000 samples |
|---|---|---|
| `ae239e26` | a round's read-length logarithms taken once, not once per distinct row | CPU halved at 500 × 30 |
| `e0dc3cf2` | the standard errors reuse the rows scored at the fitted stutter | one-base 30 reads: 12.2 → 9.1 s |
| `8a554b08` | the table of distinct rows built without a heap row per sample | → 8.2 s; 3 reads: 6.9 → 6.0 s |
| `0084426a` | a chunk keeps partial sums only for the rows its samples hold | → 7.7 s; a class of 1,000 loci: 6.3 GB → 2.9 GB peak |

Not done, because D3 shows the fit far inside its budget: keeping every class's table of rows
across rounds, and batching small classes; both would break the arch's one-class-at-a-time memory
rule.

## 4. Commits since the handoff (head `53a3e4f9`)

`ae239e26` perf (logs once) · `7265efd0` combined review of C8–C10 (a fall-off's error now follows
the owner's slipped-reads rule; the stated cause of a failed inversion corrected) · `ef307745`,
`e0dc3cf2`, `8a554b08`, `0084426a` perf · `8c6ce6f5` review of the perf commits · `830049c1` D2a
(the bench's real-psp mode) · `2c0bdf99` D1 · `e463eeee` D3 · `fb2e0225` D6 code · `3f0c5b1e` D5 ·
`932fc9b0` D6 · `2b3416bd` D2 · `c4a425a6` D4.

## 5. Left for later

- Five clippy errors from the newer clippy (1.98) in code older than tonight (`fit_bench.rs`,
  `edge_deletion_gathering.rs`, `expected_counts.rs`, `penalised_fit.rs` ×2); they also stop clippy
  from reaching the bench example.
- The cross-platform digests were re-recorded on macOS only; the Linux container did not run.
- Scratch under the worktree's `tmp/` from tonight (builds of several GB): `perf_ep/`, `perf_review/`,
  `perf_impl/`, `perf_logtable/`, `review_perf/`, `review_C_followups/`, `C_followups_fixes/`,
  `diag_errors/`, `e2e/`, `d1/head_src`, `d1/head_target`, `d2a/`, `review_d2a/`, `d3/`, `d4/`,
  `d5/`, `d6/`, `review_d6/`, `oracle_D4/`.
