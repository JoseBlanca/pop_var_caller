# Landing on templates by stretches won to the fourth power (2026-11-06)

Branch `landing-power`, from main ae21942a. Spec changed: `neighbour_evidence_pass.md` §3.3. Takes
up what `local_sources_2026-11-05.md` §5 left: with its two parents as its only sources, an F2
plant's heterozygotes are called right more often than when the pass chooses its sources itself.

## 1. What was wrong

The refinement pass walks each sample's two copies along the window, each copying one of its 16
templates and changing template now and then. Where a copy changes template, it lands on each
template in proportion to that template's weight — until now the number of stretches it won plus
½ (§3.3: the window is cut into stretches, and in each the two candidates closest to the sample
win).

On the local F2 (two inbred parents at 30×, 20 offspring at 1×, 4 Mb), with today's main:

| templates of each offspring | heterozygotes right | wrong | uncalled |
|---|---|---|---|
| chosen by the pass (main) | 85.2 in 100 | 2.2 | 8.1 |
| the two parents only (oracle) | 89.7 | 1.6 | 4.2 |
| the parents first, then the pass's choice, all landed by wins | 84.0 | 2.2 | 9.3 |
| the parents first, the others landed on a tenth as often | 90.1 | 1.5 | 3.9 |
| the pass's top two only | 53.9 | 15.1 | 26.5 |

So the parents are not missing from the templates — putting them first changes nothing — and the
siblings are needed beside them (the top two alone collapse). What costs the 4.5 heterozygotes in
100 is **how often the walk lands on the siblings**. A log of every round-1 stretch, classified by
the offspring's true state there, shows why the siblings win so much:

| offspring's true state over the stretch | stretches | right parent | wrong parent | sibling in the same state | sibling in another state |
|---|---|---|---|---|---|
| homozygous | 3,830 | 0.69 | 0.16 | 0.62 | 0.53 |
| heterozygous | 3,226 | 0.51 (either parent) | — | 0.53 | 0.97 |

(winners a stretch, two places, more with ties.) At one read a site, a heterozygous offspring reads
parent A's allele at half the sites, exactly as it would an unrelated sample's: no single candidate
can show that it is half of the pair of parents, and the stretch's places go to siblings about three
times in four. The parents still rank first, from the homozygous stretches, but landed on in
proportion to their wins the siblings took copies the parents explain with certainty, and a 1×
sibling's shares, near 50/50 where it has no read, leave a heterozygote uncalled.

## 2. The change

| what | where | spec |
|---|---|---|
| A template's landing weight is `(stretches won + ½)⁴`, not `stretches won + ½`: one that won about twice as many stretches is landed on about 16 times as often | `LANDING_WEIGHT_POWER`, `ChosenTemplate::landing_weight`, `templates.rs` | pass §3.3, §3.5 |

The walk normalises the weights (`WalkParameters::landing_weights`), so only their ratios matter.
A template that won no stretch now lands at an 81st of the weight of one that won a single
stretch, and against a leader that won 45, at about one hundred-millionth of its weight — less than
the background's 10⁻⁵. In practice the walk copies the few templates that won the most.

## 3. How it was chosen

The power was a diagnostic switch (branch `f2-templates`, not merged). On the local F2,
heterozygous offspring right in 100:

| power | 1 | 1.5 | 2 | 2.5 | 3 | 3.5 | **4** | 5 | 6 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|
| right | 85.2 | 84.2 | 85.0 | 87.0 | 88.4 | 88.4 | **89.5** | 90.2 | 89.5 | 89.0 |

The gain starts at 2.5 and is flat from 4 to 6, the oracle's level. On rick (§4) 5 gained a little
more than 4 on the F2 and the outbred cohort, and called more false variants on both (F2 with reads
39 against 32, without 11 against 2; outbred 572 against 544). 4 is the lower end of the plateau.

No principled form is claimed: the power is measured. A form that would explain it — landing by how
likely the template is to be the copy's source, which the stretch wins only estimate — is not built.

## 4. Measurements

**refine_sim on rick**, the diagnostic build at power 4 against main's build, `call-from-psps`
with refinement on from the suite's psps (`benchmark_results/pvc/7120e12d/refine_sim/`), scored as
the suite scores. SNP genotypes, all depths:

| cohort | missed, with reads | false, with reads | uncalled, without reads | wrong, without reads | heterozygotes phased right | wrong links |
|---|---|---|---|---|---|---|
| F2, 22 samples | 284 → **208** | 13 → 32 | 1,050 → **795** | 33 → 48 | 2,775 → **4,607** of 7,107 | 1 → 2 |
| outbred, 400 | 2,321 → 2,260 | 481 → 544 | 3,729 → 3,801 | 114 → 113 | 105,767 → 106,177 of 145,951 | 3 → 2 |
| selfing, 400 | 599 → 604 | 265 → 270 | 721 → 762 | 100 → 101 | 10,231 → 10,207 of 14,206 | 1 → 1 |

- **The F2 gains in every column but false variants**: phase up by two thirds, missed variants down
  by a quarter, a quarter of the genotypes without reads filled. It calls 19 more false variants
  with reads and 15 more wrong genotypes without them.
- **The natural cohorts are level**: outbred 61 fewer missed against 63 more false, out of about
  3.6 million genotypes; selfing within 41 genotypes in every class.

**Locally**, the branch's build against main's:

- the F2 (above): 89.5 heterozygotes in 100 right, the diagnostic's figure exactly; homozygous
  offspring wrong 5 → 16 of 7,213;
- the 40-sample outbred cohort (10 samples each at 30×, 10×, 5×, 2×): every depth and genotype
  class within half a point;
- a 40-sample selfing cohort (F = 0.9, same depths, built for this): within 0.5 points, the largest
  move one heterozygote of 202 at 5×;
- **small cohorts**, subsets of the outbred one with their own fitted parameters: 4 samples (one a
  depth), no class moves by more than 3 genotypes; 12 samples (three a depth), heterozygotes at 5×
  and 2× gain 13 and 10 right, homozygous references there lose 8 and 15.

## 5. What is left

- **The F2's false variants**: 19 more with reads, 15 more wrong without. Where they come from is not
  traced.
- **The rule is fitted, not derived** (§3). It was measured on an F2, two 400-sample cohorts, two
  40-sample cohorts and two small ones; a cohort with a structure unlike all of these — a
  backcross, a collection of a few families — is not measured.
- Phase on the natural cohorts is limited by the output rule, not by this: at 400 outbred samples,
  `--phase-min-path-agreement 0.75` phases 128,612 heterozygotes right with 29 wrong links of
  109,436, against 105,767 with 3 of 69,392 at the default 1 (Beagle: 135,247, 5,443 wrong of
  133,654). Not changed here.

## 6. Code review

The change is one line of arithmetic; checked by hand:

- **Overflow and precision**: a window has a few hundred stretches at most, so the largest weight is
  about 500⁴ ≈ 6·10¹⁰, far inside `f64`; the walk divides by the sum, so the ratios are exact to
  rounding.
- **A template that won nothing** keeps a weight above zero (½⁴), so the walk can still land on it;
  it is now below the background's share, so the walk reaches it about as rarely as the background.
- **Round 2 and later** take the same weights through the same function; no other caller of
  `landing_weight` exists.

## 7. Tests

- `templates::tests::ties_with_the_last_winner_win_too` — the landing weight is `(3 + ½)⁴` for 3
  wins, and one win weighs 81 times no win.
- The two bit-fingerprint tests of a phased window (`refine::tests`) re-recorded; the test window's
  phased columns rise from 2,856 to 2,871.
- The whole crate's tests pass in release, except the two `should_panic` tests outside the pass
  (`alignment::left_align_structured`, `genetics::lgamma`) that panic only through debug
  assertions, which a release build leaves out.

## 8. Reproducing

The local cohorts are built with `doc/devel/ng/research/neighbour_evidence/refine_sources/natural_cohort.sh`
(`DEPTHS=30:10,10:10,5:10,2:10`, inbreeding 0 and 0.9) and the F2 with the simulator's `cross=f2`;
they are scored with `score_f2.py` and `score_natural.py` beside the cohort script. The rick runs
use the scratch script described in `stretch_scale_2026-11-06.md` §8, with
`PVC_DIAG_LANDING_POWER` on the diagnostic build of branch `f2-templates`.
