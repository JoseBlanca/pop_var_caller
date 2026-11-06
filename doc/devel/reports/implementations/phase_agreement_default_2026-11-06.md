# The phase's path agreement defaults to 0.8 (2026-11-06)

Branch `phase-default`, from main ba43f52b. Specs changed: `neighbour_evidence_phase_output.md` §9,
`neighbour_evidence_read_phase.md` (its glossary). User page: `doc/user/neighbour_evidence.md`.

## 1. The change

The refinement pass writes two neighbouring heterozygotes of a sample as phased together when
enough of its 36 drawn paths agree on their phase and lean at both. How many is
`--phase-min-path-agreement`; its default was 1, all of them (spec §2.3). It is now **0.8, 29 of the
36** (`PathAgreement::DEFAULT`, `phase.rs`; `decide`, `stage.rs`; both subcommands' `default_value`).
The option stays, as a compromise only the user can settle: more phase against more wrong links.
The VCF header's decision line now records `path_agreement=0.8`.

## 2. Why 0.8

The owner asked for the default to be optimized for inbreeding 0 to 0.8 and samples at 5× to 30×.
Three 40-sample cohorts were simulated for it (msprime, 4 Mb over two contigs, 10 samples each at
30×, 20×, 10× and 5×), at inbreeding 0, 0.5 and 0.8, and each was called whole and as subsets of 12
samples (three a depth) and 4 (one a depth), each subset with its own fitted parameters, with main
at ba43f52b. Phase scored by `benchmarks/suite/lib/refine_score.py`; wrong links per 100 links in
phase sets:

| cohort | 1 (before) | 0.95 | 0.9 | 0.85 | **0.8** | 0.75 | 0.6 |
|---|---|---|---|---|---|---|---|
| F 0, 40 samples | 0.05 | 0.04 | 0.06 | 0.08 | **0.10** | 0.14 | 0.40 |
| F 0, 12 | 0.26 | 0.27 | 0.43 | 0.44 | **0.45** | 0.99 | 1.8 |
| F 0, 4 | 2.5 | 3.0 | 2.9 | 3.5 | **4.0** | 8.0 | 9.6 |
| F 0.5, 40 | 0.02 | 0.03 | 0.04 | 0.10 | **0.11** | 0.10 | 0.22 |
| F 0.5, 12 | 0.23 | 0.35 | 0.33 | 0.39 | **0.45** | 0.76 | 1.1 |
| F 0.5, 4 | 3.2 | 2.9 | 3.1 | 4.0 | **4.8** | 8.5 | 9.1 |
| F 0.8, 40 | 0.03 | 0.03 | 0.06 | 0.09 | **0.09** | 0.08 | 0.16 |
| F 0.8, 12 | 0 | 0 | 0 | 0 | **0.11** | 0.10 | 0.20 |
| F 0.8, 4 | 1.5 | 1.4 | 1.8 | 2.3 | **3.3** | 4.3 | 6.8 |

Heterozygotes phased right, 1 → 0.8:

| cohort | 40 samples | 12 samples | 4 samples |
|---|---|---|---|
| F 0 | 13,831 → 16,033 (+16%) | 3,107 → 3,932 (+27%) | 458 → 698 |
| F 0.5 | 8,137 → 8,471 (+4%) | 1,809 → 1,916 (+6%) | 348 → 466 |
| F 0.8 | 3,753 → 3,850 (+3%) | 1,015 → 1,057 (+4%) | 276 → 314 |

- **From 12 samples up, 0.8 keeps wrong links under half a link in 100**, against under 0.3 with all
  paths, and phases 3% to 27% more heterozygotes right — most where the phase was shortest, in the
  outbred cohorts.
- **0.75 is past the edge**: links and wrong links both jump between 0.8 and 0.75, and the
  12-sample outbred cohort reaches 1 wrong link in 100.
- **At 400 samples** (rick's refine_sim, main before the landing change, the 5×, 10× and 30×
  groups): 0.9 and 0.75 both kept wrong links under 0.05 in 100 at every depth, outbred and
  selfing; 0.8 lies between them.
- **At 4 samples no share keeps wrong links under 1 in 100** — all paths left them at 1.5 to 3.2 in
  100, and 0.8 raises them to 3.3 to 4.8. This is the failure spec §6 already states: the paths
  measure how much the draw varies, not whether the model is right, and with 3 other samples the
  templates often do not hold the sample's haplotypes.

Outside the range: on 40 outbred samples all at 2×, all paths phased 144 of 20,754 heterozygotes
right with 1 wrong link of 48, and 0.9 phased 407 with 3 of 164; on 10 of them, 0.9 gave 5 wrong of
87 (all paths 1 of 50). Such cohorts phase almost nothing at any share.

## 3. Tests

- `phase::tests` — `PathAgreement::default()` is `PathAgreement::DEFAULT`.
- `stage::tests::the_header_line_states_the_decision` — the decision line records
  `path_agreement=0.8`.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 4. Reproducing

The cohorts are built with `doc/devel/ng/research/neighbour_evidence/refine_sources/natural_cohort.sh`
(`DEPTHS=30:10,20:10,10:10,5:10`, inbreeding 0, 0.5, 0.8), called with
`call-from-psps --refine-genotypes on --phase-min-path-agreement <share>`, subsets with
`estimate-parameters` on the subset's psps first, and scored with `refine_score.py` (the phase
needs no coverage files).
