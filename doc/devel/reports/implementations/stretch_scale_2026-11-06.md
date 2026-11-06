# Refinement stretches that lengthen with the cohort (2026-11-06)

Branch `stretch-scale`, from main c0917a4c. Spec changed: `neighbour_evidence_pass.md` §3.3.
Follows `local_sources_2026-11-05.md`, whose §4.1 found that choosing templates stretch by stretch
cost the 400-sample outbred cohort about a tenth of its phase.

## 1. What was wrong

The refinement pass chooses, for each sample, the 16 other samples (from phasing round 2, the 16
haplotypes) its two copies are most likely copying. Since 2026-11-05 it chooses them stretch by
stretch: the window's compared records are cut into stretches each holding one unit of expected
disagreement — about two SNPs at frequency ½ — and in each stretch the two candidates that disagree
least with the sample win. Candidates rank by stretches won.

One unit is little evidence. An unrelated sample matches a stretch that short as well as a relative
does fairly often, and the more samples compete, the more of them do. Diagnostic counters on the
refine_sim benchmark (branch `phase-diag`, not merged) measured it in phasing round 1: a stretch had
**2.45 winning samples at 40 samples and 6.5 at 400**, where two were asked for — the extra winners
being ties at the second place. Each chance winner takes stretches from the true relatives, so at
400 samples the templates are diluted with unrelated samples, and the phase suffers first:
heterozygotes phased right in the outbred cohort fell from 108,656 to 98,871 of 145,951.

Making the winners strict — exactly two a stretch, ties broken by sample order — recovered only
1,442 of the 9,785 phased heterozygotes lost (98,871 → 100,313), so the ties were not the cause; the
chance matches themselves were. Longer stretches remove them.

## 2. The change

| what | where | spec |
|---|---|---|
| A stretch holds `I = max(1, ln(N / 8))` units of expected disagreement, `N` the other samples — 1 unit up to 21 of them, 1.6 at 39, 3.9 at 399, 4.8 at 999, 6.4 at 4,999 | `stretch_information`, `TOLERATED_CHANCE_MATCHES`, `MIN_STRETCH_INFORMATION` in `templates.rs`; `choose_templates`, `choose_haplotypes` pass it to `Disagreements::stretches` and `StretchSums::new` | pass §3.3 |

The form follows a crude count: if an unrelated candidate matches a stretch of `I` units as well as
a relative with chance about `e^−I`, the cohort's chance matches number `N·e^−I`, and `I = ln(N/8)`
holds them to about eight. The 8 is measured, not derived (§3). `N` counts samples in every round,
also from round 2 where each offers two haplotypes, because a sample takes at most one of a
stretch's two places and the scale was measured that way.

Nothing else changes: the winners (two, plus ties), the ranking and the landing weights are as
before.

## 3. How it was chosen

The length was set by an environment switch in the diagnostic build and measured on the three
refine_sim cohorts on rick: the F2 (two inbred parents at 30×, 20 offspring at 1×), and the two
400-sample natural cohorts, outbred (F = 0) and selfing (F = 0.9), 100 samples each at 30×, 10×, 5×
and 2×. SNP genotypes, refinement on, all depths together:

| cohort | units | wrong, truth variant, with reads | no call, without reads | heterozygotes phased right |
|---|---|---|---|---|
| outbred 400 | 1 | 2,771 | 4,641 | 98,871 |
| | 3 | 2,558 | 4,568 | 102,936 |
| | 4 | **2,201** | 4,327 | 108,600 |
| | 6 | 2,145 | 4,466 | 111,110 |
| | 10 | 2,247 | 5,494 | 111,771 |
| selfing 400 | 1 | 588 | 897 | 10,059 |
| | 3 | 562 | 895 | 10,054 |
| | 6 | 594 | 1,006 | 10,111 |
| F2 | 1 | **266** | 966 | **3,087** |
| | 2 | 298 | 1,090 | 2,817 |
| | 3 | 337 | 1,153 | 2,678 |
| | 6 | 338 | 1,113 | 2,770 |
| | 10 | 303 | 1,190 | 2,726 |

- **The outbred cohort wants long stretches**: at 4 units it phases 9,729 more heterozygotes right
  than at 1 and calls 570 fewer true variants wrong; past 6 units the no-calls start to rise.
- **The selfing cohort does not care**: every length is within about 1 in 100 of every class. Its
  samples share long runs of one haplotype, so a relative leads by many units whatever the length.
- **The F2 wants one unit.** Every longer length phased fewer heterozygotes and called more wrong.
  Locally, with crossovers ten times rarer than the benchmark's, two units still called 3 to 13 more
  heterozygotes in 1,000 wrong or missing than one did, so it is not the dense crossovers. Why the
  F2 wants short stretches is not traced.
- On the 40-sample natural cohorts (local, 1×–30×), lengths 1 to 6 measured level, within 1 in 100
  of every genotype class.

A fixed length cannot serve both the F2 and the outbred 400; a length growing with the cohort can,
because the F2 has 22 samples and the outbred cohort 400. The first version used `ln(N/7)`, which
gave the F2 1.1 units. A finer sweep of the F2 on the diagnostic build put its best at exactly one
unit:

| units | 0.8 | 0.9 | **1.0** | 1.05 | 1.1 | 1.15 | 1.2 | 1.3 | 1.5 |
|---|---|---|---|---|---|---|---|---|---|
| wrong, truth variant, with reads | 305 | 297 | **266** | 275 | 296 | 302 | 333 | 299 | 312 |
| heterozygotes phased right | 3,043 | 3,038 | **3,087** | 2,924 | 2,906 | 2,968 | 2,889 | 2,908 | 2,837 |

Neighbouring lengths move the wrong calls by 20 to 40 either way, so that row is mostly noise; the
phase is the clearer signal, about 140 fewer phased right at every length above one. `ln(N/8)` keeps
cohorts of up to 22 samples at one unit and moves the 400-sample cohorts from 4.0 to 3.9 units.

## 4. Measurements: the branch against main on rick

The branch (696dd542) against main c0917a4c, both built on rick, `call-from-psps` with refinement on
from the suite's psps (`benchmark_results/pvc/7120e12d/refine_sim/`), scored as the suite scores.
SNP genotypes of samples:

| cohort, depth | wrong, truth variant, with reads: main → branch | wrong, truth reference, with reads | no call, without reads | heterozygotes phased right |
|---|---|---|---|---|
| outbred 400, 2× | 2,016 → **1,495** | 412 → 309 | 4,488 → **3,561** | 16,365 → 19,099 |
| outbred 400, 5× | 821 → 628 | 172 → 145 | 227 → 159 | 24,677 → 27,087 |
| outbred 400, 10× | 223 → 197 | 38 → 27 | 5 → 8 | 27,291 → 29,219 |
| outbred 400, 30× | 1 → 1 | 0 → 0 | 1 → 1 | 28,027 → 30,362 |
| **outbred 400, all** | **3,061 → 2,321** | 622 → 481 | 4,721 → 3,729 | **96,360 → 105,767** of 145,951 |
| selfing 400, 2× | 412 → 388 | 194 → 156 | 769 → 678 | 2,005 → 2,033 |
| selfing 400, 5× | 196 → 182 | 103 → 86 | 48 → 42 | 2,391 → 2,533 |
| selfing 400, 10× | 35 → 29 | 17 → 23 | 1 → 1 | 2,592 → 2,694 |
| **selfing 400, all** | 643 → 599 | 314 → 265 | 818 → 721 | 9,847 → 10,231 of 14,206 |
| F2, 22 samples | 284 → 284 | 13 → 13 | 1,050 → 1,050 | 2,775 → 2,775 of 7,107 |

- **Outbred, 400 samples: the phase lost on 2026-11-05 comes back**, 9,407 more heterozygotes
  phased right (about 6 in 100 of them), at every depth. Wrong genotypes fall at every depth below
  30×, most at 2×, where a quarter of the missed true variants are now called (2,016 → 1,495) and a
  fifth of the genotypes left uncalled for want of reads are now filled. Wrong links: 3 before and
  after.
- **Selfing, 400 samples: small gains** in every class, 384 more phased right, 93 fewer wrong
  genotypes with reads. At 10× false heterozygotes rise by 6 (17 → 23).
- **F2: unchanged, to the bit** — its 21 other samples keep it at one unit, so it runs the same
  arithmetic as main.
- Wrong genotypes without reads (outbred) fall 181 → 109 and false variants there 9 → 5.
- The first version, `ln(N/7)`, gave the F2 1.1 units on the same build pair and cost it 63 phased
  right (2,775 → 2,712) and 32 more wrong with reads (284 → 316); at 400 samples it measured within
  about 1 in 100 of this version in every class (outbred phased right 105,312, wrong with reads
  2,872).
- Speed was not measured cleanly: rick was shared, and the same main build took 681 s on the
  outbred cohort in one run and 893 s in another. The rule adds no work beyond comparing more
  records a stretch; the stretch count falls as the stretches lengthen.

These are main's numbers of today, not those of §3's diagnostic build, which branched before other
merges (the F2 phases 3,087 right there, 2,775 on today's main); each table compares within one
build pair.

## 5. What is left

- **Why the F2 wants one unit** is not traced; it is not crossover density (§3).
- **Above 400 samples nothing is measured.** The rule gives 4.8 units at 1,000 samples and 6.4 at
  5,000; the outbred cohort was level from 4 to 6 units and lost no-calls at 10.
- **The length follows the cohort, not the samples with reads.** In a region where most of a large
  cohort has no reads, few samples really compete, yet stretches stay as long as the whole cohort
  asks. Not measured.
- Between 22 and 400 samples only the 40-sample cohorts were measured, at 1 unit (they would now get
  1.6 and measured level across that range).

## 6. Code review

A review of the first version (2026-11-06). Fixed:

- The docs said the F2 got one unit; `ln(21/7)` gave it 1.1. (Superseded: the scale is now 8 and
  the F2 does get one, §3.)
- A rounding slip in the docs, 4.9 for 5.0 units at 999 samples.
- The function's parameter was named `candidates`, which from round 2 are haplotypes; it is
  `other_samples`, and the docs say samples are counted in every round.

Checked, no change needed: one sample gives `ln(0)`, which `max(1)` turns into one unit, not NaN.

## 7. Tests

- `templates::tests::stretches_lengthen_with_the_log_of_the_cohort` — one unit up to 21 other
  samples, more from 22, and 1.58, 3.91, 4.83 and 6.44 units at 39, 399, 999 and 4,999.
- `templates::tests::a_stretch_ends_once_it_holds_enough_to_judge_by` — now passes the length,
  and also cuts the same records at 0.8 units (three stretches) and 1.5 (one, the leftover joining
  it).
- The two bit-fingerprint tests of a phased window (`refine::tests`) re-recorded: the test window's
  cohort is large enough that its stretches lengthen. Its phased columns rise from 2,822 to 2,856.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 8. Reproducing

On rick, a scratch script not in the repository, `tmp/phase_drop/rick_variant.sh <label> <env>
<cohort>`, re-runs `call-from-psps` with
refinement on from the suite's psps of `benchmark_results/pvc/7120e12d/refine_sim/<cohort>`, and
scores it with the suite's own `pvc_eval` commands and `benchmarks/suite/lib/refine_score.py`;
`BUILD` picks the pvc build. The diagnostic build is branch `phase-diag` (`PVC_DIAG_STRETCH_INFO`
sets the length, `PVC_DIAG_STRICT_WINNERS` makes the winners strict).
