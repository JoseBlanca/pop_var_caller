# The fit stage — Milestone D: the four commands end to end, and what the fitted numbers change

**Date:** 2026-09-05
**Plan steps:** [parameter_prepass_runs.md](../../ng/impl_plan/parameter_prepass_runs.md) Milestone D, steps D1 and D2
**Branch:** `ng-psp-mode`
**Script:** [`scripts/ng_fit_stage_end_to_end.sh`](../../../scripts/ng_fit_stage_end_to_end.sh)

## The answer

**The four commands compose on real reads**, and calling a cohort with numbers fitted from its own
data rather than with the compiled-in constants **removes 82 of 599 records and changes 115
genotypes in 3,102**.

The 82 it removes are the marginal ones: **median QUAL 5.5, against 125.3 for the 517 records both
runs called**. The two distributions barely meet — the highest QUAL among the dropped is 88.9,
below the median of the kept.

## What was run

Six tomato accessions (`SRR7279481`, `488`, `501`, `533`, `536`, `537`) over the two 100 kb
intervals of `tmp/c1_two_regions.bed`, at about three reads a position:

    generate-psps        6 samples: 3,589,737 bytes of psp and 1,306,395 of census
    generate-census      the same censuses, rebuilt from the stored psps
    estimate-parameters  6 samples fitted, 36,889 bytes of parameters file
    call-from-psps       twice: once --defaults, once --parameters

**The two routes to a census still agree byte for byte on real reads** — all six identical — and
that now covers the minted read-error totals as well, since both producers accumulate them
through the same writer.

## What the fit actually changed, and the honest reading of it

**The comparison is between a tomato cohort called with a human heterozygosity and the same
cohort called with its own fitted curve.** That is what the two files say:

| | reference concentration | alternative total | where it came from |
|---|---|---|---|
| `--defaults` | 1.0 | 0.001 | `stated_heterozygosity` — human data, and the file's own words are *"the one that rests on nothing this run measured"* |
| fitted | 2.2533 | 0.0026 | `fitted_curve` — this cohort's own allele-frequency density |

So the size of the change is partly a fact about **calling tomato with a human default**. On a
human cohort the two would sit much closer, and nothing here says how much closer.

**What the direction says.** The fitted pair is both more concentrated (2.256 against 1.001 in
total) and slightly richer in alternative belief (about 1 part in 858 against 1 in 1,000). A more
concentrated prior takes more evidence to move, which is consistent with the records it drops
being the low-QUAL ones — but this run measured the outcome, not the mechanism, and the two
should not be confused.

## What the file could not fit, and says so

**Contamination.** The file writes no `[contamination]` section, and its own note says what that
means: *"nobody identified any … that is not the same as measured and found clean"*. Six samples
is below what the estimator needs — it wants about a dozen, because the allele frequencies it
judges each sample against are fitted from the run itself. **The reads are then scored by the read
likelihood's plain formula rather than by a correction of size zero**, which is the right
behaviour and is distinguishable in the output from a measured zero.

**The inbreeding coefficients**, which are declared on this route and recorded as `supplied`.

So the file reports **5 of its 7 groups of numbers as fitted from reads**, and names the two that
were not.

## One figure that must not be read as a general one

**On this ground the census keeps 198,182 of the 200,000 analysed bases, and 151 repeat tracts.**
The budget is two million positions and this BED is 200 kb, so the selection keeps very nearly
everything; on a whole tomato genome the same budget keeps about 1 base in 400. Every count here —
the census sizes, the tract count, and how much of the fit rests on tracts — is a fact about a
200 kb BED.

**And 599 records over 200 kb at three reads a position is a small cohort's worth of calls.** The
82 and the 115 are counts from this corner. What they establish is that the fitted numbers reach
the calls and move them in the direction the QUAL distribution predicts, not how much they would
move on a genome.

## What the census cost when it began carrying the minted totals

**1,306,395 bytes against 1,305,915** for the same six accessions over the same intervals — 480
bytes across six samples, or 80 a sample, for two numbers a read group and the table's own length.

## Validation

`cargo test --lib` in the container: 6,271 passed, 0 failed, 15 ignored. `cargo fmt` and
`cargo clippy --all-targets --all-features -D warnings` are clean.
