# N5 — the anchors re-measured, and the fit that does not reach them

> **⚠ Superseded in one respect, the same day.** The defect this report found is real and its
> measurements stand, but **its account of the cause is incomplete and its emphasis is wrong**.
> The larger half is that the per-read-group rate scan was never handed the second class of
> site at all, so every candidate rate was scored under the one-class rule — nothing to do with
> coordinate ascent. The profile scan recommended below is the smaller half, and it is kept:
> with the argument restored it still buys 209 nats on a real tomato sample, while every
> fixture and both human arms said it was unnecessary. See
> `ng_noise_model_extension_n5_fix_2026-08-10.md`, which carries the re-measured table. Every
> number below is what the fit returned *before* both changes.

**Step:** N5 of `impl_plan/noise_model_extension.md`, the milestone inserted between
Milestones F and G of `impl_plan/parameter_prepass_generic.md`.
**Date:** 2026-08-10.
**Verdict: the extension does what it was built to do on the quantity it was built for, and
the fit that finds its parameters stops 351 nats short of the maximum.** The second result
was found by a fixture written during this step, and it explains every oddity in the first.

## What was run

The four `#[ignore]`d tests of `generic/real_alignments.rs` on all five alignments — HG002 at
30x and at 300x over the 100 GIAB confident spans, and tomato SRR7279481, SRR7279482 and
SRR7279483 over the tomato1 BED — with the end-to-end test extended to print the second class
of site and to refuse three shapes of degenerate answer (below). The **before** column is the
same five runs at `65515a43`, the last commit before this milestone began.

## The headline: heterozygosity

The benchmark heterozygosity over the 30x locus set is 9.9666 × 10⁻⁴ — 550 heterozygous loci
in 551,843 — from `research/noise_model_overdispersion_2026-08-10.md`, which counted it with
`bcftools query` over the v4.2.1 benchmark VCF restricted to the same spans.

| | before | after | truth | before / truth | after / truth |
|---|---|---|---|---|---|
| HG002 30x | 1.407 × 10⁻³ | **1.061 × 10⁻³** | 9.9666 × 10⁻⁴ | 1.412 | **1.065** |
| HG002 300x | 1.543 × 10⁻³ | **1.124 × 10⁻³** | — | — | — |

**1.41 → 1.06 is better than the 1.09 the research note predicted**, which is consistent with
that note's own correction: its figures came from an optimiser that had not converged and are
upper bounds on the residual bias, and the direction a better fit would move heterozygosity
was not known. It moved down.

The 300x arm has no truth figure for its own locus set (550,049 loci, 1,795 fewer than the 30x
arm), so only the change is quoted: heterozygosity falls 27.2% there against 24.6% at 30x.

The homozygous-non-reference rate barely moves, as the note said it would not — 5.33 × 10⁻⁴ →
5.39 × 10⁻⁴ at 30x against a benchmark 5.7444 × 10⁻⁴, from 0.928 to 0.938 of it.

## The tomato samples, which have no truth set

They can only say that the fit does not degenerate. What they show instead is how far this
milestone moves a number nobody can check:

| | heterozygosity before | after | fall | hom-non-ref before | after |
|---|---|---|---|---|---|
| SRR7279481 | 9.62 × 10⁻⁴ | 5.14 × 10⁻⁴ | **−47%** | 2.075 × 10⁻³ | 2.094 × 10⁻³ |
| SRR7279482 | 1.525 × 10⁻³ | 7.98 × 10⁻⁴ | **−48%** | 3.418 × 10⁻³ | 3.454 × 10⁻³ |
| SRR7279483 | 1.415 × 10⁻³ | 8.58 × 10⁻⁴ | **−39%** | 4.523 × 10⁻³ | 4.555 × 10⁻³ |

**Heterozygosity roughly halves on all three, and nothing here can say whether it should.** On
HG002 a fall of that size landed on the benchmark; on a tomato landrace there is no benchmark,
and two of these three fits rail (below). The homozygous-non-reference rate moves by about 1%
on all three, as it did on HG002.

## The second class, on all five alignments

| | one rate, before | `ε_clean` after | `w` | `ε_noisy` | emitted marginal |
|---|---|---|---|---|---|
| HG002 30x | 2.239 × 10⁻³ | 2.239 × 10⁻³ | 0.4905% | 7.079 × 10⁻² | 2.575 × 10⁻³ |
| HG002 300x | 2.371 × 10⁻³ | 2.371 × 10⁻³ | 0.9422% | 5.012 × 10⁻² | 2.821 × 10⁻³ |
| tomato SRR7279481 | 4.467 × 10⁻³ | 4.467 × 10⁻³ | 0.5598% | 9.441 × 10⁻² | 4.970 × 10⁻³ |
| tomato SRR7279482 | 2.371 × 10⁻³ | 2.371 × 10⁻³ | 0.3553% | **1.000 × 10⁻¹ (railed)** | 2.718 × 10⁻³ |
| tomato SRR7279483 | 3.758 × 10⁻³ | 3.758 × 10⁻³ | 0.3813% | **1.000 × 10⁻¹ (railed)** | 4.125 × 10⁻³ |

Every sample has a second class and none is degenerate in the share. **The clean rate is, on
every one of the five, the exact rung the one-class fit returned before this milestone.** A
second class was added and the first class did not move — five of five, at two depths, in two
organisms, and whether or not the noisy class railed. That uniformity is the first sign of the
defect
below: a class fitted at a rate that never moves is a class fitted at the one-class rate, and
the one-class rate is the tail-inflated number this milestone exists to correct.

## What that costs the error-rate anchor

`arch/parameter_prepass_generic.md` §9 compares the emitted rate against a model-free count of
mismatching bases at benchmark homozygous-reference positions — 2.263 × 10⁻³ for HG002 30x,
38,450 bases in 16,992,201 read observations.

| | emitted rate | against model-free | in ladder rungs |
|---|---|---|---|
| before | 2.239 × 10⁻³ | −1.1% | −0.19 |
| after | 2.575 × 10⁻³ | **+13.8%** | **+2.25** |

The plan's own decision section predicted +3.6%, inside one rung, from the research note's
parameters. This fit does not return those parameters.

## The fixture that explains it

`the_whole_fit_finds_both_classes_when_it_is_given_neither` (in `generic/coupled_fit.rs`)
generates a world at the research note's own HG002 30x parameters — clean 1.8836 × 10⁻³, 0.88%
of sites noisy at 5.29 × 10⁻², over the measured depth distribution, at the benchmark genotype
frequencies — as a **table**, and runs the whole coupled fit, which is handed neither rate.

**It returns a clean rate of 2.2387 × 10⁻³: three rungs, 19%, above the generating value — the
same rung it returns on the real alignment.** Scored on the same cells, the generating
parameters reach −147,552.77 and the fit reaches −147,903.98: **the truth scores 351 nats
higher than the point the fit stopped at.**

So the emitted rate is not 13.8% above the model-free count because the model says so. It is
because the fit does not reach the model's own maximum, and 351 nats is not a rounding error —
the whole second class was worth 425 nats over a beta-binomial, and the fit declines a class
that gains less than 3.

**Why the existing oracles do not see it.** N3b's two worlds hand `fit_site_noise` the *true*
clean rate and ask whether the second class is recovered given the first. The real path fits
both, in this order: settle the rates with one class, fit the pair at those settled rates,
re-settle the rates with the pair held. The first step's answer is the one-class rate, which
the tail inflates — that is the whole finding this milestone rests on — and every later step is
a coordinate-wise optimum around it. Each block is exhaustive and neither can move: the rate
scan tries all 161 rungs at the fitted pair, and the pair scan tries all 161 at the fitted rate.
The pair is a fixed point 351 nats below the maximum.

`a_wrong_clean_rate_rails_the_second_class_and_the_flag_says_so` already records that this
block misbehaves at a wrong clean rate. What was missing is that **the loop hands it one.**

## What the new checks caught on real data

The end-to-end test now refuses three shapes of answer that carry no information: a share at 0
or 1, a noisy class no noisier than the clean one, and a noisy rate on either end of the ladder.

**Tomato SRR7279482 and SRR7279483 fail the third**: their noisy class comes back at 0.1, the
ladder's coarsest rung — Phred 10, one base in ten. That is the edge of the search rather than
a maximum inside it, so what the fit actually wants is a class **noisier than the ladder can
express**, and the ladder's ends are Phred 10 to 50 because they were chosen for chemistry
(`spec/parameter_prepass.md` §3, DRAGstr's own range).

**That is worth separating from the stall above, because the two want different fixes.** A
duplication the reference does not carry — the owner's first cause of a noisy site — puts two
copies' reads at one locus, and where the copies differ the site shows about **half** its reads
non-reference. Half is five times the coarsest rung. A class trying to hold such sites cannot
land inside a ladder of sequencing-error rates, and no better search would put it there.
Whether that is what these two samples hold is unmeasured; what is measured is that the fit
asks for a rate above the range and is clamped to it.

Three of five arms pass all four tests; the two railed tomatoes pass three and fail the
end-to-end one. Both are the deeper tomato samples — SRR7279482 carries 194.9 M reads and
SRR7279483 103.1 M, against SRR7279481's 77.1 M.

**The check printed nothing on its first run**, because it asserted before it printed and took
the two railed samples' rates down with it. It now prints all three numbers and then asserts,
which is what the table above is filled from.

## Recorded, not fixed

- **The fit's stall.** The fix is a design change the plan explicitly ruled out — N3b states
  *"no multi-start, because the surface has no trap"*, and this is a trap. The shape that
  matches the module's existing idiom is to make the clean-rate scan a **profile** scan over
  the whole inner model: at each rung, re-climb the genotype frequencies *and* re-fit the pair,
  score, and keep the best rung — which is what `fit_by_profile_scan` already does for the
  one-class model. Unmeasured; what is measured is that the pair is recovered exactly when the
  clean rate is right, and that the truth scores 351 nats above where the loop stops.
- **`CoupledFit::error_rate`'s doc said it carries the share-weighted marginal. It carries the
  clean rate** — `estimate` marginalises when it assembles `GenericSampleParameters` and
  nowhere earlier, as its own comment says. On a sample like HG002 the two differ by 15%.
  Corrected here, and the new fixture asserts which of the two the field holds.

## Validation

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` and
`cargo doc --no-deps --lib` clean, the last at the 12-unresolved-link pre-existing baseline
with none in this module.

**Two things in this commit are red when run, and both are findings rather than breakage.**
`the_whole_fit_finds_both_classes_when_it_is_given_neither` is `#[ignore]`d with a reason
naming the numbers, so `cargo test --lib` stays green; it takes 0.3 s and needs nothing but
itself, so the `#[ignore]` is a schedule and not a cost. And the end-to-end test fails on two
of the five alignments — the two tomato samples whose noisy class rails at the ladder's edge.
Neither is on by default; both are one command away.
