# Which samples the refinement copies: three changes for low-coverage samples (2026-11-05)

Branch `refine-sources`, from main 1edce6ca. Specs changed: `neighbour_evidence_pass.md` §3.3 and
§3.5, `neighbour_evidence_phasing.md` §2.2, §2.3 and §5.

## 1. What was wrong

The refinement pass (`--refine-genotypes`) re-calls each sample's genotypes by working out, along
the chromosome, which other samples' haplotypes its two copies are copying. On the `refine_sim`
benchmark's F2 — two inbred parents at 30× and 20 offspring at 1×, one crossover per contig per
meiosis — it got the homozygous offspring genotypes right (54 → 89 in 100) and left the
heterozygous ones almost where calling left them (27 → 30 in 100), at the 716 SNPs where the
parents are different homozygotes.

The cause is that the pass preferred, at three points, the samples that carry least information:

1. **Choosing the samples to copy.** The 16 candidates were ranked by their plain disagreement with
   the sample. A 1× sample's genotype probabilities are soft, and two soft samples differ little —
   so 1× siblings ranked above the 30× parents, whose probabilities are 0 or 1. Records where every
   sample carries the same allele counted as agreement too, and those are most of a 1× sample's
   confident records (357 of the F2's 1,073 SNP sites do not separate the parents). For most
   offspring neither parent was among the 16.
2. **The background template.** Beside the 16, each walk copies a template whose allele
   probabilities are the cohort's frequencies. In an F2 every allele stands near 50%, and at one
   read a site a source that predicts either allele with even odds explains a heterozygote exactly
   as well as the two parents do. As one template among 17 alike it took the heterozygous
   stretches.
3. **Soft templates in the walk's emissions.** A 1× sibling's shares sit at the frequencies wherever
   it has no read, so it imitates the background and is never plainly wrong either.

With the two parents as the only sources and no background (a diagnostic switch, not a feature),
heterozygous offspring reached 90 in 100 right, so the walk itself can do it.

## 2. The changes

| change | where | spec |
|---|---|---|
| Rank candidates by disagreement **against what the cohort's frequencies predict for an unrelated sample**: observed `1 − Σ q_self·q_other` over expected `1 − Σ q_self·p`, summed over shared records, scored pessimistically as before | `templates::choose_templates` | pass §3.3 |
| A copy that changes template **lands on the background with probability 10⁻⁵**, the 16 sharing the rest | `WalkParameters::background_landing_share`, `model.rs`; `refine::BACKGROUND_LANDING_SHARE` | pass §3.3, §3.5 |
| From round 2, the **emissions read one allele drawn from each template's shares** (the background's included); the genotype prediction and the drawn haplotypes keep reading the shares | `WalkSite::emission_shares`, `model::site_terms`; drawn in `refine::solve_round` from `StreamPurpose::DefiniteAlleles` | phasing §2.2, §2.3, §5 |

The third is the one that needed care. Drawing definite alleles everywhere — prediction included —
gave the F2 the same gain and cost the outbred 1× samples 12 homozygous-reference genotypes in 100:
at a site without reads, the call is what the copied template carries, and a drawn guess turns an
unsure template into a confident wrong call or, averaged over rounds, into a call too unsure to
write. Keeping the shares for the prediction removed that. Leaving the background soft while every
other template was definite made it the only source never plainly wrong; the walk copied it more,
and the no-read calls pass §4.2 withholds when the background carries more than a tenth of the
weight rose from 8 to 21 in 100 of the outbred 1× homozygous-reference genotypes. Drawing its
allele too removed that.

## 3. Measurements

Three cohorts besides the F2, all simulated on the laptop with the benchmark's simulator, 4 Mb over
two contigs, haplotypes from msprime (`scripts/msprime_haplotypes.py`), 40 samples each:

- **outbred, 30–5×** and **selfing (F = 0.9), 30–5×**: 10 samples at each of 30×, 20×, 10×, 5×;
- **outbred, 30–1×** and **selfing, 30–1×**: 10 at each of 30×, 5×, 2×, 1×.

Genotypes right, at the truth's biallelic SNPs, by the sample's true genotype; the F2 at the 716
SNPs where the parents are different homozygotes. "Off" is calling without refinement; "main" the
pass before these changes; the last column the branch as built.

| cohort, depth, true genotype | count | off | main | ranking + background last resort | + definite alleles in the emissions (branch) |
|---|---|---|---|---|---|
| F2 offspring 1×, heterozygous | 7,107 | 27.0% | 29.6% | 50.8% | **69.6%** |
| F2 offspring 1×, homozygous | 7,213 | 51.6% | 89.1% | 91.8% | **92.4%** |
| outbred 1×, heterozygous | 5,257 | 33.1% | 45.1% | 51.0% | **52.7%** |
| outbred 1×, homozygous alternative | 2,315 | 52.5% | 70.7% | 77.7% | **79.6%** |
| outbred 1×, homozygous reference | 18,688 | 50.2% | 65.1% | 72.0% | **72.7%** |
| outbred 2×, heterozygous | 5,215 | 54.1% | 68.0% | 72.0% | **74.1%** |
| selfing 1×, heterozygous | 838 | 10.4% | 64.3% | 67.4% | **72.8%** |
| selfing 1×, homozygous alternative | 3,996 | 60.3% | 91.8% | 91.7% | **92.1%** |
| outbred 5× (of 30–5×), heterozygous | 5,257 | 86.3% | 96.1% | — | **96.5%** |
| selfing 5× (of 30–5×), heterozygous | 838 | 79.8% | 92.7% | — | **93.6%** |

- **Wrong calls do not rise to pay for the gains.** Heterozygous F2 offspring called wrong: 31.2 →
  6.7 in 100. Outbred 1× heterozygotes called wrong: 18.3 → 13.9; homozygous alternative 1.7 → 2.3;
  homozygous reference 0.4 → 0.6. Selfing 1× homozygous reference wrong: 0.4 → 0.7.
- **The samples at 10× and above** move by at most 0.2 in 100 in any cohort.
- "Main" was measured with main at e6565a37; the later main the branch starts from changed the
  walk's memory layout only, its output bit for bit the same (824ec349). The column before the last
  was measured with the diagnostic build that preceded the branch; the branch's draws come from the
  pass's own random streams rather than that build's, and reproduce its full combination within
  1.6 in 100 everywhere.
- These figures were measured at 40 samples, at most 4 Mb, 1× to 30×; the refine_sim benchmark on
  rick repeats them at 400 samples.

**Cost.** Wall time of `call-from-psps` with the pass on, median of three alternating runs, 6
threads: outbred 30–1× cohort 11.38 s on the branch against 11.30 s on main, the pass itself about
2.9 s of it (8.4 s with the pass off); F2 2.87 s against 2.85 s. The scatter between runs is larger
than either difference.

## 4. Tried and not adopted

- **Definite alleles in the prediction too** — the F2 gain, and a 12-in-100 loss on outbred 1×
  homozygous references (§2).
- **Definite alleles only where a template is at least 90% sure** — safe for the outbred cohort,
  but the F2 kept only 54 in 100 of its heterozygotes: the gain comes from exactly the unsure
  alleles.
- **Hidden records passed on with the sample's own reads.** The fill check (pass §4.3) hides up to
  50 records of each sample a window from its own reads, and other samples copy the guessed
  haplotype there. With the parents as the only sources this cost 10 heterozygotes in 100 (90 → 80);
  with all samples as sources, nothing measurable (70.7 against 68.9 in 100 on the F2, within a
  point on the natural cohorts). Hiding stays as it is.

## 5. What is left

With the parents as the only sources the F2 reaches 90 heterozygotes in 100; the branch reaches 70.
The 1× siblings are still among the 16 copied beside the parents. How to choose sources well
without being told who the parents are is the next investigation.

## 6. Tests

- `model::tests::the_walk_with_rare_background_landings_is_a_plain_forward_backward` — the walk's
  posteriors against a forward-backward over the full transition matrix, written out state by
  state, background landing at 0.05: agree to 1e-9.
- `model::tests::a_sites_draw_follows_forward_times_transition` — now at even and uneven landing
  weights.
- `model::tests::an_even_landing_share_is_no_landing_share`,
  `model::tests::the_emissions_and_the_prediction_read_their_own_shares`.
- `templates::tests::an_unsure_candidate_is_not_taken_for_a_close_one`; `an_absent_sample_ranks_last`
  kept by ranking a candidate with nothing expected to disagree on after every measured one.
- The two bit-fingerprint tests of a phased window (`refine::tests`) re-recorded: the change is to
  the pass's arithmetic on purpose.
- The whole crate's tests pass in release, except two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 7. Reproducing

`doc/devel/ng/research/neighbour_evidence/refine_sources/`: `natural_cohort.sh` builds a cohort
(simulation, mapping, psps, fitted parameters; inside the dev container), `score_natural.py` and
`score_f2.py` score a call set as the tables above, `fill_check_against_truth.py` sets the pass's
fill check beside the true error of its no-read calls.
