# ng step 4, the STR path — E4: the monotonicity walk

*Implementation report, 2026-08-12. Step E4 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md) — its own commit, as the
plan requires. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.3,
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §4.1 rule 4.*

## What the step is

`merge_until_monotone` holds each period's fitted slippage levels rising with the repeat count.
Slippage genuinely rises with repeat count, so a sequence that dips — tracts of 7 copies coming
out less slippery than tracts of 6 — is reporting the noise in one stratum rather than a fact
about repeats. Where it dips, the two tables are pooled, refitted as one, and both strata report
the pooled answer with `merged_over` naming the set. Pooling can itself dip against the run below
it, so it repeats until the sequence rises — the pool-adjacent-violators shape.

The pooled model is built from the **lowest** repeat count in the set, which is the intersection
of the two supports and not their union: a tract of four copies cannot carry an allele six copies
shorter, so scoring the pooled table against the longer tract's alleles would let the fit place
mass on lengths half its loci cannot have.

## ⛦ A finding for the checkpoint: the search does not recover a level of 0.2

Building the fixtures turned this up, and it is the most important thing in the step.

A stratum whose every locus shows **two short reads in ten** has a slippage level of 0.2 with
every slip in the losing direction. The search returns **1e-5** — the bottom of its range — with
all four starts agreeing to 1.00, at both the coarse and the fine setting. The likelihood at the
truth is about 2,300 nats better than what it returns: 1,200 identical loci score −1,436 at a
level of 0.2 against the −3,750 the returned answer scores, and −3,750 is exactly what a
heterozygous locus with no slippage gives (`1200 × ln(C(10,2)/2¹⁰)`), so the search has settled on
explaining every locus by its genotype rather than by slippage.

**One short read in ten is recovered correctly** (0.0998 against a truth of 0.1), and so are 0.02,
0.04 and 0.06. So the failure is confined to a regime above where real strata sit: the measured
range is 0.0009 to 0.15, and the fixture that fails is at 0.2 with the direction share on its
boundary.

**The diagnostic cannot see it.** `start_spread` is 1.00, which is the value that means "all four
starts agree" — and this is precisely the reading `fitting/multistart.rs` warns about in its own
doc: on a surface where each axis is line-searched over its whole range, agreement is what a
search that never found anything also produces.

The likely mechanism is the coordinate sweep: the level and the direction share are coupled, and
three of the four starts begin with a direction share of 0.35 or more where the truth is ~0. At
those shares any level that puts 20% of reads a copy short also puts several percent a copy long,
which the data has none of, so the level axis is driven to zero before the direction axis is ever
searched. **Not diagnosed further and not fixed here**: the search is E2's and the harness that
could settle it is the exact-bias one, so this is a question for the owner rather than a change to
make inside a merge step.

The E4 fixtures were moved to the band real strata occupy, and the reason is written on the
fixture rather than left implicit.

## An ordering the plan and the composition disagree about

The plan numbers borrowing (E3) before the merge (E4) and says the merge is last "because it reads
the fitted sequence". But a borrowed stratum has no fitted level of its own to be out of order, so
the only sequence the merge can read is the *fits* — and a stratum that borrowed before the merge
holds a value its lender has since changed. So `merge_until_monotone` takes the map of fits and
returns a map of fits, which composes either way, and the natural order is fit → merge → borrow.
**Which order the entry point wires is a question for the checkpoint**; nothing is blocked, because
that wiring is F1's.

## What the review changed

Two agents; the reliability agent ran 21 mutations of which 6 survived.

**The locus floor was not applied here**, only in the walk that searches. A 999-locus stratum
fitting near zero pulls a correctly fitted neighbour from 0.0599 to 0.0327 — 1.83-fold, against
the 15-to-25% a merge is priced at — and stamps the pair as merged. It now applies the same
`thick_enough_to_fit` the borrowing does, and for the same reason: the function is public and may
be handed a map it did not build.

**A pooled refit's starts were never checked.** "Every reported level came from a search that
settled" is an invariant of the whole step, and a refit is a new search, so the function now
returns a `Result` and re-applies the check. No pooled table has yet been found that trips it.

**Four doc claims did not survive checking**, three of them mine and one an inherited slip: a
field named `fitted_over` where the code has `merged_over`; "it is under the locus floor for
exactly that reason", which was false twice over; "close to the loci-weighted mean of their
levels", where what the refit returns is the maximum over the union of the two tables and is 4.3
times away from that mean in the collapsed case above; and a fixture described as slipping "one
read in ten where the one below it slips at two", when the fixture always shows one short read in
ten and the level comes from how many loci slip — which the helper's own doc comment says nine
lines earlier.

**And the reviewer corrected my account of the search failure.** The truth is not scored at
−1,436: that is the analytic ideal, and the model cannot reach it because the direction share
floors at 0.005. A grid search maximises at (0.200, 0.005, 0.002) scoring −1,453.7, so the gap is
**2,296 nats** — "about 2,300" stands, but the mechanism I gave does not. A single start placed
**exactly at the truth** still returns 1e-5. The cause is that the level axis is genuinely
bimodal — *every locus heterozygous and nothing slipping* is a real local maximum at the bottom —
and the golden section's first two probes both land in that mode with the function falling
between them, so the whole upper half of the axis is discarded on step one, at every direction
share tried.

**The mutation run's six survivors were one pattern again**: every fixture was one library, one
ploidy, one motif period and 1,200 loci, so three legs of the grouping filter and the
no-fit branch cost nothing to delete; no fixture cascaded, so pooling twice was untestable; and
the one fixture with two equal levels — which is what a pooled-on-equality rule needs — called
`fit_slippage` by hand instead of the function under test.

## Tests

Eight.

| test | what it pins |
|---|---|
| `a_rising_sequence_of_levels_is_left_alone` | the control: 0.02, 0.04, 0.06 pass through with no fit moved and nothing pooled |
| `a_dip_in_the_sequence_is_merged_and_refitted` | 0.04, 0.02, 0.06 pools the first two, both name the set, the pooled level lands between them over 2,400 loci, and the stratum above is untouched |
| `two_strata_that_agree_are_not_merged` | equal is not a dip |
| `a_pooling_that_still_dips_is_pooled_again` | the cascade: 0.04, 0.06, 0.005 ends as one set of three over 3,600 loci |
| `the_walk_merges_inside_one_library_ploidy_and_period` | each of the three legs, in a fixture where dropping any one puts a dip in front of the walk |
| `a_stratum_under_the_locus_floor_is_not_merged_with_its_neighbour` | the floor, applied to the map it is handed |
| `merging_two_strata_that_agree_costs_exactly_nothing` | the plan's own control — two strata at the same level pool to a **bit-identical** model, so what a real merge costs is the distance between its parts and not the pooling |
| `a_merged_fit_is_scored_against_the_shorter_tracts_alleles` | eleven allele lengths and 66 genotypes from a merge of a four-copy and a six-copy stratum |

Re-run after the fixes: pooling on equality, pooling once instead of until it rises, and dropping
all three legs of the grouping — all three now fail, at three different tests.

## Validation

`cargo fmt --check`, `cargo clippy --lib --all-features -- -D warnings` clean, and
`cargo test --lib --bins --tests --all-features` in the container: 3,470 → **3,478** lib tests, 0
failed, 11 ignored.
