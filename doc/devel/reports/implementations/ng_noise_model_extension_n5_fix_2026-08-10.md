# The rate scan never saw the second class — the defect N5 found, and the two changes that fix it

**Follows:** `ng_noise_model_extension_n5_2026-08-10.md`, which measured the anchors and found
the fit three ladder rungs and 351 nats away from the maximum. **That report's diagnosis of the
cause was wrong** and is corrected here; its measurements stand, and the remedy it recommended
turns out to be half of the answer.

## What was actually wrong

`fit_read_group_error_rates` — the block that picks each read group's error rate, and the only
place the fit's error rate ever comes from — **was never handed the sample's second class of
site.** It built one `SampleLibraryNoise` per rung with `SampleLibraryNoise::single`, which
carries no site noise, so every candidate rate was scored under the one-class rule against a
table whose tail belongs to the other class.

The consequence is not subtle and it is not a search problem: the scan returns the
tail-inflated rate whatever pair sits beside it, so the outer loop's second half — *settle the
rates again with the second class held fixed* — was re-deriving the number it already had. The
clean rate came back on the one-class rung on **all five real alignments and on a generated table**, which is exactly what a block scoring under the wrong model looks like.

**The first change is the argument.** `site_noise: Option<SiteNoise>` is now a parameter of
that function, and each rung's noise is built with the pair when the sample has a second class.
That alone takes the fixture below from failing by 351 nats to passing.

**The second change is the profile N5's report recommended, and deciding whether it was still
needed took a measurement.** With the argument fixed, the fixture passes in 0.5 s with no
profile at all, and so does a generated table shaped like a tomato sample — so every fixture
available said the profile could be deleted, and it was. **Real tomato SRR7279481 said
otherwise**: scored on its own cells, the argument alone reaches −1,504,289.10 and the profile
reaches −1,504,079.98, **209 nats higher**, reporting 1.42% of sites noisy at 6.310 × 10⁻²
against 1.07% at 7.079 × 10⁻². Both HG002 arms return the identical answer either way. So the
profile is kept, and what justifies it is a real alignment rather than a fixture.

**What the profile is.** The noisy rate is held at each of the 161 rungs in turn while
everything else — the per-read-group clean rates, the genotype frequencies and the noisy share
— is refitted around it, and the best-scoring rung wins. That is a profile likelihood in one
parameter, the shape `fit_by_profile_scan` already uses for the one-class model. One dimension
and not two, because the noisy rate is per sample where `ε` is per read group.

**It pins the rate through `fit_site_noise` handed a ladder of one rung**, rather than through
a second copy of that function's expectation-maximisation. The first draft bypassed it
entirely, which left the function N3b built — and both of its oracles — testing code no caller
ran.

What it costs: 23 s → 37 s on the tomato run, and 0.5 s → **17.2 s** on the fixture. What that
costs the whole suite is inside its own run-to-run spread — thirteen runs of the same command
ranged 79 to 90 seconds — so there is no number to quote. **An earlier draft of this report said
39.5 s**, which timed the draft of the profile that did not route its share climb through
`fit_site_noise`.

## What the fix recovers

`the_whole_fit_finds_both_classes_when_it_is_given_neither` generates a table at the research
note's own HG002 30x parameters and runs the whole fit, which is handed neither rate. Before
the fix it returned a clean rate 19% high and scored 351 nats below the generating parameters.
After it, every rate is the generating rung, the share is within 10⁻³ of the generating 0.88%,
and the fit's score is no lower than the truth's.

## On the five real alignments

| | one rate, before | `ε_clean` | `w` | `ε_noisy` | emitted marginal | heterozygosity |
|---|---|---|---|---|---|---|
| HG002 30x | 2.239 × 10⁻³ | 1.884 × 10⁻³ | 0.8770% | 5.309 × 10⁻² | 2.333 × 10⁻³ | 1.081 × 10⁻³ |
| HG002 300x | 2.371 × 10⁻³ | 1.995 × 10⁻³ | 1.2662% | 4.217 × 10⁻² | 2.504 × 10⁻³ | 1.146 × 10⁻³ |
| tomato SRR7279481 | 4.467 × 10⁻³ | 3.758 × 10⁻³ | 1.4194% | 6.310 × 10⁻² | 4.601 × 10⁻³ | 5.52 × 10⁻⁴ |
| tomato SRR7279482 | 2.371 × 10⁻³ | 2.113 × 10⁻³ | 0.4175% | **1.000 × 10⁻¹ (railed)** | 2.522 × 10⁻³ | 7.85 × 10⁻⁴ |
| tomato SRR7279483 | 3.758 × 10⁻³ | 3.350 × 10⁻³ | 0.4921% | **1.000 × 10⁻¹ (railed)** | 3.825 × 10⁻³ | 8.35 × 10⁻⁴ |

**The clean rate now moves on every one of the five**, by two to three rungs, where before the
fix it sat on the one-class rung on all five.

### Two implementations that share no code now agree

The research note fitted HG002 in Python, from the same cells but through its own
expectation-maximisation. Against this fit:

| | `ε_clean` | `ε_noisy` | `w` |
|---|---|---|---|
| Python, 30x (research note) | 1.895 × 10⁻³ | 5.29 × 10⁻² | 0.88% |
| **this fit, 30x** | **1.884 × 10⁻³** | **5.309 × 10⁻²** | **0.877%** |
| Python, 300x (research note) | 1.952 × 10⁻³ | 4.24 × 10⁻² | 1.28% |
| **this fit, 300x** | **1.995 × 10⁻³** | **4.217 × 10⁻²** | **1.266%** |

Every pair is within a rung of the ladder, at both depths. The note's own caveat — that its
figures came from an optimiser stopped after 150 iterations and are upper bounds — is why this
is worth stating: the two disagree by less than the note's own convergence error.

## The two anchors

**Heterozygosity, against the benchmark's 9.9666 × 10⁻⁴ over the 30x locus set:**

| | heterozygosity | / benchmark |
|---|---|---|
| one rate (before the milestone) | 1.407 × 10⁻³ | 1.412 |
| two classes, rate scan blind to them | 1.061 × 10⁻³ | 1.065 |
| **two classes, fixed** | **1.081 × 10⁻³** | **1.085** |

The research note predicted 1.091 and this fit gives 1.085. **The stalled fit's 1.065 was
closer to the benchmark and was not a better answer** — it was 351 nats short of the maximum,
and a wrong fit landing nearer the truth on one of its four parameters is luck, not accuracy.

**The error rate, against the model-free count of 2.263 × 10⁻³** at benchmark
homozygous-reference positions:

| | emitted rate | against model-free | in rungs |
|---|---|---|---|
| one rate (before the milestone) | 2.239 × 10⁻³ | −1.1% | −0.19 |
| two classes, rate scan blind | 2.575 × 10⁻³ | +13.8% | +2.25 |
| **two classes, fixed** | **2.333 × 10⁻³** | **+3.1%** | **+0.53** |

The plan predicted +3.6% and inside one rung, from the note's parameters. Measured: +3.1%, half
a rung. **The design decision to emit the share-weighted marginal is vindicated by measurement
rather than by the note it was argued from.**

## What does not survive the measurement

**The milestone's claim that the clean rate is more depth-invariant than the single rate is not
reproducible on this ladder.** The research note reported the drift across the four-fold depth
change falling from +6.1% to +3.0%. Measured here:

| | 30x | 300x | drift |
|---|---|---|---|
| single error rate | 2.239 × 10⁻³ | 2.371 × 10⁻³ | +5.90% = **one rung** |
| `ε_clean` | 1.884 × 10⁻³ | 1.995 × 10⁻³ | +5.89% = **one rung** |

Both move by exactly one rung, and a rung is 5.925%. The note's Python fit had a continuous
rate and could express a 3.0% drift; this ladder cannot express anything between zero and one
rung, so the comparison it made is not available here. **What does reproduce is the
heterozygosity drift**: +9.7% on one rate (1.407 → 1.543 × 10⁻³) against +6.0% on two classes
(1.081 → 1.146 × 10⁻³), where the note reported +9.0% and +5.3%.

**And the two deep tomato samples still rail.** SRR7279482 and SRR7279483 still put their noisy
class on the ladder's coarsest rung, 0.1, so the fit still wants a class noisier than a ladder
chosen for chemistry can express. The fix moved their clean rates — 2.371 → 2.113 × 10⁻³ and
3.758 → 3.350 × 10⁻³ — and did not change what they are asking for. That is the separate
decision N5's report set out, and it is still open.

## Validation

`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
`cargo test --lib --bins --tests --all-features` (3,234 passed, 0 failed, 9 ignored) and
`cargo doc --no-deps --lib` at the 12-unresolved-link pre-existing baseline, none in this
module. The suite's `#[ignore]` count is back to F3's nine: the fixture that recorded the
defect now passes and runs by default, in 17.2 s.

Three of the five real-alignment arms pass all four tests; the two railed tomatoes fail the
end-to-end one on the ladder-end check, which is the finding above and not a regression.
