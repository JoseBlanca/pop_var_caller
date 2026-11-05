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
| F2 offspring 1×, heterozygous | 6,416 / 7,458 | 71.6 (6.2) / 67.4 (7.3) | **82.3 (2.6) / 80.6 (2.9)** |
| F2 offspring 1×, homozygous | 7,364 / 7,082 | 94.4 (0.2) / 93.8 (0.1) | 94.2 (0.2) / 93.7 (0.2) |
| outbred 1×, heterozygous | 5,552 / 5,333 | 53.3 (13.6) / 53.4 (14.9) | **58.1 (11.5) / 58.2 (12.0)** |
| outbred 1×, homozygous alternative | 2,363 / 2,558 | 81.6 (1.5) / 82.1 (1.9) | 84.3 (2.0) / 83.9 (2.2) |
| outbred 1×, homozygous reference | 18,995 / 19,209 | 75.4 (0.4) / 73.4 (0.6) | 75.8 (0.5) / 73.6 (0.6) |
| outbred 2×, heterozygous | 5,282 / 5,547 | 75.4 (11.1) / 75.2 (11.1) | **79.1 (8.2) / 79.5 (8.2)** |
| outbred 5×, heterozygous | 5,447 / 5,678 | 93.0 (5.0) / 93.6 (4.3) | 94.4 (3.7) / 94.6 (3.4) |
| selfing 1×, heterozygous | 335 / 495 | 83.0 (5.4) / 85.1 (6.9) | **87.2 (1.8) / 91.5 (1.6)** |
| selfing 2×, heterozygous | 523 / 795 | 93.3 (3.1) / 87.4 (3.5) | 93.3 (2.7) / 90.2 (2.3) |
| selfing 5×, heterozygous | 557 / 579 | 93.0 (3.4) / 95.2 (2.6) | 95.7 (0.7) / 97.6 (0.5) |
| any cohort, 30×, any genotype | | | moves by at most 0.1 in 100 |

- **The gains are in heterozygotes and in fewer wrong calls.** F2 heterozygotes called wrong fall
  from about 7 in 100 to under 3; selfing 1× heterozygotes called wrong from about 6 in 100 to under
  2.
- **The cost**: homozygous alternatives at 1× in the outbred cohort are called wrong slightly more
  often (1.5 → 2.0 and 1.9 → 2.2 in 100) while more of them are called right; F2 homozygous offspring
  lose 0.1 to 0.2 points. Neither is larger than the difference between the two copies.
- **Selfing heterozygotes are few** — 335 to 795 a group — so their figures move by a point for every
  3 to 8 genotypes.
- The clean build was checked against the diagnostic switch it replaces on three of these cohorts:
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

## 5. What is left

- **A shallow close relative can rank below a deep distant one**: it cannot win a stretch it has no
  reads in (spec §3.3). Not measured in a cohort built to show it.
- With the parents as the only sources, the F2's heterozygotes reach 90 in 100; the branch reaches 81.
  The remaining gap is not traced.

## 6. Tests

- `templates::tests::the_samples_a_copy_copies_win_the_stretches_they_match` — a sample copying one
  candidate over each half of the window, against three candidates fairly close everywhere: the two
  it copies are its templates, as the whole-window score would not choose them.
- `templates::tests::a_shallow_relative_wins_the_stretches_it_has_reads_in`,
  `a_near_tie_goes_to_the_sample_with_more_data` — rewritten for the stretch rule, which changes
  their order: a shallow relative wins every stretch it has a read in, and ranks below deep samples
  that win the rest.
- `templates::tests::a_stretch_ends_once_it_holds_enough_to_judge_by`,
  `ties_with_the_last_winner_win_too`.
- `model::tests::the_walk_with_rare_background_landings_is_a_plain_forward_backward` — now also with
  uneven landing weights on the sample templates, against the full forward-backward: agree to 1e-9.
- The two bit-fingerprint tests of a phased window (`refine::tests`) re-recorded: the change is to
  the pass's arithmetic on purpose. The phased columns of the test window rise from 2,648 to 2,872.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 7. Reproducing

The cohorts are built with `doc/devel/ng/research/neighbour_evidence/refine_sources/natural_cohort.sh`
and scored with `score_natural.py` and `score_f2.py` beside it, as for the change before.
