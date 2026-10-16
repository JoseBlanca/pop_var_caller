# D6: does a slip rate 0.5% of its value low change calls? (2026-10-16)

Plan step D6 of [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md); the question is spec
§5.1: the production fit's slip rate comes out about half a percent of its value low, and the step asks
whether that bias moves genotype calls.

## 1. The answer

**No, not by an amount that matters.** A slip rate exactly 0.5% of its value low changes about 4 calls
in 100,000, where 5 to 6 calls in 100 are already wrong under the true numbers, and the calls it
changes go as often to right as to wrong. The fitted stutter as a whole — whose slip rates are off by
up to 2–5% class by class — changes 2 to 5 calls in 10,000, again in both directions about equally.

Every figure is from the bench's drawn cohorts (the tomato-like PCR truth, at most 100 loci a class,
seeds 1 to 3 summed). The call is a proxy: each sample called alone from its genotype weights.

| motif length | samples × reads | pairs with reads | wrong under the truth | slip rate exactly 0.5% low: changed (right→wrong / wrong→right) | fitted slip rate only | fitted stutter (all three numbers) |
|---|---|---|---|---|---|---|
| one-base | 500 × 3 | 1,383,966 | 86,753 (6.3%) | 52 (26 / 21) | 162 (78 / 68) | 299 (114 / 163) |
| two-base | 500 × 3 | 1,251,394 | 59,430 (4.7%) | 55 (16 / 31) | 70 (25 / 35) | 154 (63 / 72) |
| one-base | 2,000 × 3 | 5,536,053 | 353,344 (6.4%) | 240 (85 / 130) | 387 (152 / 201) | 2,153 (1,040 / 1,010) |
| two-base | 2,000 × 3 | 5,005,561 | 241,727 (4.8%) | 206 (90 / 88) | 196 (72 / 99) | 2,466 (1,035 / 1,266) |

A changed call whose old and new genotypes are both wrong counts in "changed" but in neither
direction, which is why the two directions do not add up to it.

**Where the fitted stutter's changes sit.** At 2,000 × 3 two-base, 1,887 of the 2,466 changes are in
seed 3, and most of those in the long classes of 22 to 27 repeats (233 to 518 changes a class, up to
0.8% of the class's calls at 27 repeats), where the fitted shorter share and fall-off move calls more
than the slip rate does. They still go both ways (seed 3: 814 right→wrong, 955 wrong→right). At one
base the 18-repeat class holds most of the changes at 500 samples (43 of the 52 in the 0.5% case).

**How far the fitted slip rate was off**, fitted / true − 1, per class: median −0.67% and −0.63%
(one- and two-base) at 500 samples, −0.44% and −0.53% at 2,000, ranging from −5.0% to +5.3% at 500
and from −2.4% to +2.9% at 2,000. That spread, not the 0.5% bias, is what moves the fitted-stutter
calls.

## 2. What a call is, and what was held fixed

- **A call** is the genotype with the largest weight among the 190 genotypes over the fit's window
  (±9 repeats; a copy further out counts as ±9): the weight the fit's own expectation step gives a
  sample's genotype — the prior from the class's allele-length shares, F and the sample's homozygote
  excess, times the probability of its reads under the stutter. Calling proper also uses the other
  samples at the locus; this does not, and nothing in the project simulates reads for this many
  samples.
- **A changed call** is a (locus, sample) pair with reads whose call differs from its call under the
  truth's numbers; **a wrong call** is one that is not the genotype the sample drew.
- **The allele-length shares and F are held at the truth** in every column above, so a changed call
  has one cause. With the fitted shares and F as well, 50 to 130 times more calls move (60,189 and
  78,953 at 2,000 × 3), more of them from wrong to right than back (35,394 against 16,822 one-base);
  that comparison is about the shares, not the stutter, and is recorded in the bench's output only.

## 3. How it was run

`ng_str_fit_bench --calls-out` (commit `fb2e0225`), `--methods D_tract_prior`, at 500 × 3 (4 threads)
and 2,000 × 3 (all threads, each run alone under the memory guard), both motif lengths, seeds 1–3.
Scripts and per-class tables: `tmp/d6/` (`run_500.sh`, `run_2000.sh`, `n500_summary.md`,
`run_2000.log`). Review: `reports/reviews/str_penalised_fit_D6_2026-10-16.md`.
