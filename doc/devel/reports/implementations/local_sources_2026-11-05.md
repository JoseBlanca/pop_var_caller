# Choosing the samples to copy stretch by stretch (2026-11-05)

Branch `local-sources`, from main 11a97683. Specs changed: `neighbour_evidence_pass.md` §3.3 and
§3.5, `neighbour_evidence_phasing.md` §2.3. Follows
`refine_sources_2026-11-05.md`, which left the F2's heterozygous offspring at 70 right in 100
against 90 with the two parents as the only sources.

## 1. What was wrong

The refinement pass re-calls each sample by working out which other samples' haplotypes its two
copies are copying along the window. It copies at most 16 other samples (`K`) and, from the second
phasing round, 16 haplotypes. Both were chosen by one score over the whole window: how little the
candidate disagrees with the sample across its 500 compared records.

A sample's two copies change source at every crossover in its ancestry. In an F2 plant, each copy
follows one parent up to a crossover and the other after it. Over the whole window each parent is
therefore half right and half wrong, while a sibling whose crossovers fall near the plant's own
looks closer than either. The whole-window score ranks highest whoever is fairly close everywhere,
not the samples that are exactly right over a part each — and the walk needs the second kind.

## 2. The change

| what | where | spec |
|---|---|---|
| The compared records are cut, in order, into **stretches holding one unit of expected disagreement** — about two SNPs at frequency ½ where the sample is certain | `Disagreements::stretches`, `templates.rs` | pass §3.3 |
| In each stretch the **two candidates with the lowest pessimistic score there win**, with any tied with the second; candidates **rank by stretches won**, then by whole-window score, then sample order | `Disagreements::judged`, `choose_templates`, `choose_haplotypes` | pass §3.3, phasing §2.3 |
| A copy that changes template **lands on each template in proportion to its stretches won plus ½**; the background keeps its 10⁻⁵ | `ChosenTemplate::landing_weight`; `WalkParameters::landing_weights`, `model.rs`; `refine::solve_round` | pass §3.5 |

The disagreement is the one the pass already used (`refine_sources_2026-11-05.md`): observed
`1 − Σ q_self·q_other` against expected `1 − Σ q_self·p`, the cohort's frequency `p`, scored by the
upper bound of its Beta posterior. From round 2 the haplotype choice now uses the same expected
disagreement, `Σ p(a)·P(genotype without a)`, in place of the count of shared records it divided by
before, so both rounds judge alike.

The walk takes the landing weights as a new argument (`walk`, `walk_and_draw`,
`walk_and_draw_paths`, `SampleFit::fit_window`); `None` weighs the templates alike, as before.

## 3. How it was chosen

Five ways of choosing were built behind diagnostic switches (branch `copy-choice-diag`, not merged)
and run on four 40-sample cohorts and the F2, then the two best on fresh simulations:

| approach | what it does | F2 heterozygotes right |
|---|---|---|
| fewer templates (2, 4, 8) | the whole-window ranking, cut shorter | worse everywhere |
| **A, stretch wins** | as built | 80.7% |
| B, stretches grown until candidates separate | a stretch ends only once the third-best has fallen behind the second | worse than A |
| C, Beagle's composite templates | each template re-chosen in every stretch | 18% |
| D, IMPUTE5's exact match | from round 2, the haplotypes sharing the longest run of likeliest alleles around each 25th record | 78.1% |
| A then D | A in round 1, D after | 83.0% |

- **C fails** because a template re-picked every few records from the sample's own reads carries
  those reads back to it: the walk then mostly confirms its first guess.
- **The settings of A** — one unit a stretch, two winners, landing weight wins + ½ — are the best of
  those tried: longer stretches, four winners and fixed stretches of 50 records all measured worse.
  Landing every template alike instead of by wins cost the F2 7 heterozygotes in 100 (74 against 81)
  and made no difference on the natural cohorts.
- **A against A then D, on fresh simulations** (§4): A then D was 1.6 points better on one F2 copy,
  0.1 worse on the other, level on the outbred cohort and 1 to 4 points worse on the selfing 1×
  heterozygotes. D adds a second mechanism for no measurable gain, so A alone was built.

## 4. Measurements

Two fresh copies of each cohort, simulated with new seeds by the benchmark's simulator — 4 Mb,
haplotypes from msprime (`scripts/msprime_haplotypes.py`) for the natural cohorts:

- **F2**: two inbred parents at 30×, 20 offspring at 1×, one crossover per contig per meiosis; scored
  at the SNPs where the parents are different homozygotes;
- **outbred** and **selfing (F = 0.9)**: 40 samples, 10 each at 30×, 5×, 2× and 1×; scored at the
  truth's biallelic SNPs.

Genotypes called right, with called wrong in brackets, "main" the pass on 11a97683 and "branch" this
change; copy 1 / copy 2:

| cohort, depth, true genotype | genotypes (copy 1 / 2) | main | branch |
|---|---|---|---|
| F2 offspring 1×, heterozygous | 6,416 / 7,458 | 71.6 (6.2) / 67.4 (7.3) | **82.1 (3.1) / 80.1 (3.4)** |
| F2 offspring 1×, homozygous | 7,364 / 7,082 | 94.4 (0.2) / 93.8 (0.1) | 94.5 (0.2) / 93.9 (0.1) |
| outbred 1×, heterozygous | 5,552 / 5,333 | 53.3 (13.6) / 53.4 (14.9) | **58.0 (11.3) / 58.1 (12.0)** |
| outbred 1×, homozygous alternative | 2,363 / 2,558 | 81.6 (1.5) / 82.1 (1.9) | 83.3 (1.9) / 83.5 (2.2) |
| outbred 1×, homozygous reference | 18,995 / 19,209 | 75.4 (0.4) / 73.4 (0.6) | 75.8 (0.5) / 73.8 (0.6) |
| outbred 2×, heterozygous | 5,282 / 5,547 | 75.4 (11.1) / 75.2 (11.1) | **78.8 (8.8) / 79.4 (8.4)** |
| outbred 5×, heterozygous | 5,447 / 5,678 | 93.0 (5.0) / 93.6 (4.3) | 94.6 (3.5) / 94.8 (3.2) |
| selfing 1×, heterozygous | 335 / 495 | 83.0 (5.4) / 85.1 (6.9) | **86.6 (1.8) / 90.7 (2.2)** |
| selfing 2×, heterozygous | 523 / 795 | 93.3 (3.1) / 87.4 (3.5) | 93.3 (2.7) / 90.1 (2.1) |
| selfing 5×, heterozygous | 557 / 579 | 93.0 (3.4) / 95.2 (2.6) | 95.5 (1.1) / 97.8 (0.3) |
| any cohort, 30×, any genotype | | | the same to 0.1 in 100 |

- **The gains are in heterozygotes and in fewer wrong calls.** F2 heterozygotes called wrong fall
  from about 7 in 100 to about 3; selfing 1× heterozygotes called wrong from about 6 in 100 to about
  2.
- **No genotype class loses** beyond the difference between the two copies: outbred 1× homozygous
  alternatives called wrong move 1.5 → 1.9 and 1.9 → 2.2 in 100 while more of them are called
  right.
- **Selfing heterozygotes are few** — 335 to 795 a group — so their figures move by a point for every
  3 to 8 genotypes.
- **The branch column is the build after the code review** (§6), which changed two rules: from
  round 2 a stretch's two winning places go to two samples, so a homozygous sample's equal copies
  do not take both, and what is left after a window's last full stretch joins that stretch rather
  than being judged on its own. Against the build before them (approach A of §3) every figure
  above moves by under 1 in 100 — the largest, selfing 1× heterozygotes, by 0.6 and 0.8, 2 to 4
  genotypes.
- Before the review, the clean build was checked against the diagnostic switch it replaces on three
  of these cohorts:
  the F2 call set is identical; the outbred one differs at 7 of 40 samples on the first contig only,
  where the per-sample walk parameters are still being fitted and the landing weights' rounding,
  computed in a different order, tips some near-tied fits. Its scores agree to 0.1 in 100.
- These figures were measured at 40 samples, 1× to 30×, on 4 Mb. The refine_sim benchmark on rick
  repeats them at 400 samples, beside GATK, Beagle 4.1 and STITCH.

**Cost.** `call-from-psps` with the pass on, on the outbred copy 1 (40 samples, 4 Mb), 6 threads,
three runs of each alternating: median wall 10.31 s on the branch against 10.29 s on main, CPU
45.7 s against 45.6 s. The laptop carried another session's load (load average about 20), so the
runs are noisy, but the difference is smaller than the scatter. In the synthetic timing test
(`diag_source_choice_timing` on `copy-choice-diag`, not merged), the choice itself costs 6.6 ms a
sample at 400 samples against 4.8 ms before, and grows with the cohort as the old choice did — every
candidate against every compared record; the refine_sim benchmark on rick measures it at 400
samples.

### 4.1 The refine_sim benchmark on rick, 400 samples (2026-11-06)

`run_suite.sh --only refine_sim --caller pvc`, results in `benchmark_results/pvc/7120e12d/`, against
the pvc run before (`8f635ece`, which predates `refine_sources_2026-11-05.md` too, so the change
measured is both) and the competitors' run of 2026-11-05. The F2 is the same design as §4; the
natural cohorts have 400 samples, 100 each at 30×, 10×, 5× and 2×. SNPs, genotypes of samples, with
the refinement pass on:

| cohort, depth, reads | wrong, truth variant: before → now | wrong, truth reference | no call | right |
|---|---|---|---|---|
| F2 1×, with reads | 2,208 → **266** | 56 → 61 | 5 → 2 | 11,077 → 13,028 |
| F2 1×, without reads | 14 → 45 | 0 → 2 | 3,449 → **966** | 4,531 → 6,990 |
| outbred 2×, with reads | 1,842 → 1,838 | 1,510 → 1,549 | 140 → 27 | |
| outbred 2×, without reads | 40 → 132 | 7 → 9 | 25,682 → **4,407** | 81,040 → 102,173 |
| outbred 5×, without reads | 5 → 9 | 0 → 0 | 1,067 → 223 | 4,714 → 5,550 |
| selfing 2×, with reads | 708 → **383** | 308 → 398 | 5 → 2 | |
| selfing 2×, without reads | 34 → 116 | 7 → 23 | 2,700 → 840 | 35,931 → 37,589 |
| selfing 5×, with reads | 266 → 171 | 177 → 233 | 2 → 0 | |

- **The F2 gain repeats at full size**: heterozygotes called homozygous fall from 2,208 to 266, and
  where an F2 plant has no read, the pass now fills 7 in 8 genotypes rather than 4 in 7.
- **The natural cohorts fill far more genotypes without reads** (outbred 2×: no calls 25,682 →
  4,407), at about 1 wrong for every 160 more filled right.
- **Selfing**: false heterozygotes (wrong, truth reference) rise by 90 at 2× and 56 at 5×, against
  325 and 95 fewer missed heterozygotes.
- **Against the competitors** (same run's tables): at 2× outbred, with reads, pvc misses 1,838 true
  variant genotypes, STITCH 1,878, Beagle 7,917, GATK 11,874; pvc leaves 27 without a call, STITCH
  28,997. Without reads, Beagle calls everything and gets 1,886 wrong; pvc leaves 4,407 uncalled and
  gets 141 wrong. On the F2, pvc gets 266 wrong with reads, STITCH 471, Beagle 4,144.
- **Phase moved two ways.** Outbred: heterozygotes phased right fall from 108,656 to 98,871 of
  145,951 (links 77,025 → 60,231; wrong links 13 → 1). Selfing: 8,875 → 10,059 of 14,206. F2: 122 →
  3,087 of 7,107. Which of the two changes cost the outbred phase is not measured.
- **Speed**: the pass (call-from-psps with it on, less with it off) takes 26 s at 400 outbred
  samples against 39 s before, 24 s against 33 s selfing. Peak memory of call-from-psps with the
  pass rises 1,707 → 1,782 MB outbred and 1,058 → 1,368 MB selfing.
- Not this change: pvc's outbred 5× and 10× groups carry about 2,200 to 2,900 false variant
  genotypes with reads (wrong, truth reference), refinement off and on, before and now, against
  under 300 for every competitor. Not traced.

## 5. What is left

- **A shallow close relative can rank below a deep distant one**: it cannot win a stretch it has no
  reads in (spec §3.3). Not measured in a cohort built to show it.
- With the parents as the only sources, the F2's heterozygotes reach 90 in 100; the branch reaches 81.
  The remaining gap is not traced.

## 6. Code review

A review of the branch (2026-11-05, ten findings). Fixed:

- **A homozygous sample took both winning places** in rounds 2 and later, where candidates are
  haplotypes and its two copies score alike — the rule of two winners, one a copy, became one sample
  a stretch in an inbred cohort. The places now go to two samples.
- **The window's leftover stretch** won like a full one however little it held; it now joins the
  stretch before it.
- **Memory**: every candidate's disagreement at every compared record was kept (`2·N·500` numbers a
  sample in a phasing round, about 80 MB at 5,000 samples); now only the sums a stretch are, a pair
  a candidate a stretch (`StretchSums`).
- Docs: the module doc's F2 figures disagreed with the spec's; the constants now say which cohorts
  their tuning was measured on; `BACKGROUND_LANDING_SHARE` says its figures predate landing by wins.

Measured, no change needed:

- **The fill check and the hidden records.** The fill check (pass §4.3) hides records from the
  sample's own reads after its templates were chosen from them; with stretches of about two SNPs a
  hidden record weighs more in the choice than in a 500-record score, which could make the check
  report fewer wrong fills than no-read calls really have. Against the truth it does not: where the
  check has enough genotypes (30× and 5×) its error rate is at or above the true error of the no-read
  calls, before and after the change (outbred copy 1, 5×: check 6.7 → 5.2 in 100, true 1.7 → 0.8;
  selfing copy 2, 5×: check 1.0 → 0.7, true 0 both).

Stated, not fixed: at thousands of samples chance matches over a one-unit stretch may swamp a true
relative's lead (spec §3.3, not measured); a shallow close relative ranks below deep samples that
win the stretches it has no reads in (spec §3.3, §5 here).

## 7. Tests

- `templates::tests::the_samples_a_copy_copies_win_the_stretches_they_match` — a sample copying one
  candidate over each half of the window, against three candidates fairly close everywhere: the two
  it copies are its templates, as the whole-window score would not choose them.
- `templates::tests::a_shallow_relative_wins_the_stretches_it_has_reads_in`,
  `a_near_tie_goes_to_the_sample_with_more_data` — rewritten for the stretch rule, which changes
  their order: a shallow relative wins every stretch it has a read in, and ranks below deep samples
  that win the rest.
- `templates::tests::a_stretch_ends_once_it_holds_enough_to_judge_by`,
  `ties_with_the_last_winner_win_too`, `a_homozygous_sample_does_not_take_both_places`.
- `model::tests::the_walk_with_rare_background_landings_is_a_plain_forward_backward` — now also with
  uneven landing weights on the sample templates, against the full forward-backward: agree to 1e-9.
- The two bit-fingerprint tests of a phased window (`refine::tests`) re-recorded: the change is to
  the pass's arithmetic on purpose. The phased columns of the test window rise from 2,648 to 2,840.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 8. Reproducing

The cohorts are built with `doc/devel/ng/research/neighbour_evidence/refine_sources/natural_cohort.sh`
and scored with `score_natural.py` and `score_f2.py` beside it, as for the change before.
