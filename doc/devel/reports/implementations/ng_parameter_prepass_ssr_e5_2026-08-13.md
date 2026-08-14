# ng step 4, the STR path — E5: the record, and the summary a person reads

*Implementation report, 2026-08-13. Step E5 of
[`parameter_prepass_ssr.md`](../../ng/impl_plan/parameter_prepass_ssr.md), the last step of
Milestone E. Design authority:
[`spec/parameter_prepass_ssr.md`](../../ng/spec/parameter_prepass_ssr.md) §4.4, §4.5, §4.6,
[`arch/parameter_prepass_ssr.md`](../../ng/arch/parameter_prepass_ssr.md) §2.4, §4.3.*

## What the step is

The four fits each answer for one stratum at a time. This step puts their answers together into
the record a reader opens when a number looks wrong — `StratumFit`, one per
`(library, stratum, ploidy)` — and aggregates them into the one a reader opens when nothing does:
`StratumFitSummary`, one per library. A sample carries several hundred of the first against the
SNP/indel path's four fits in total, so the second is not a convenience: **a flag nobody reads is
how a badly-fitted parameter reaches a caller**, which is spec §4.4's reason for the summary
existing at all.

`assemble_sample_parameters(accumulators, substitutions, slippage, merges)` is the whole surface,
plus `unexplained_locus_share`, which is the one piece of new mathematics.

## The one piece with no code before this step: `unexplained_locus_share`

Spec §4.6's diagnostic — the share of a stratum's loci whose shape its **own fitted model** calls
very unlikely. Each entry is scored under the fitted slippage parameters, mixed over the fitted
genotype frequencies, divided by the entry's whole-repeat depth, and compared with
`UNEXPLAINED_SHAPE_LN_LIMIT`; the shares are weighed by loci.

**Per read, and that is the whole of why it is not a floor on the shape's likelihood.** A shape's
likelihood is a product over its reads, so a floor on the total is a floor on depth. Measured on a
stratum whose loci come in four genotype classes at twelve reads each and which the fit explains
entirely: **every one of its 1,000 loci sits below the −3.00 limit on the undivided total**, from
−4.37 to −12.52, while per read the same four run from −0.36 to −1.04.

**One term of the fit's own score is taken out: how many orders the reads could have arrived in.**
A shape's likelihood carries `ln(n! / Π n_b!)` beside the reads' own probabilities. It cancels out
of a mixture — which is why the scoring rule keeps it, and why it makes the sums-to-one gate an
identity — but it does not cancel out of a comparison between one locus and a threshold, and it
grows with how many distinct lengths a locus shows, which is the very thing this diagnostic looks
for. Measured, it pays **+0.23 a read to an ordinary locus of one length and +1.07 to a locus
showing four**, so the canonical collapsed duplication would have arrived with nine tenths of the
evidence against it cancelled. That was the review's finding, not mine; `UNEXPLAINED_SHAPE_LN_LIMIT`
is now set against the convention without it, and Milestone H calibrates it there.

**What the denominator is.** Every locus in the stratum, including those no read of which sits on
the whole-repeat grid. Those score an empty product — a likelihood of one — because the model
scores only whole-repeat movement, so they are loci with nothing to say about which alleles they
carry rather than loci that refute them all. Counting them as unexplained would report the guard
share a second time under another name.

**Reported, never acted on**, and the plan refuses in advance to drop the loci that score badly:
that is threshold-then-count, the bias the whole step exists to remove, and it would take real long
alleles before it took artefacts.

## ⛦ What the diagnostic can and cannot see, measured

Five loci in 1,005 whose reads sit at three lengths four motif copies apart — the shape a
duplication the reference does not carry produces — are seen: the share comes back at exactly
5/1,005, those loci score −5.23 a read against the ordinary loci's −0.33, and **the guard share
stays at exactly zero**, because every one of those reads moved by a whole number of copies. That
is why spec §4.6 asks for a second diagnostic rather than a threshold on the first.

**Plant twenty times as many of exactly the same locus and the diagnostic goes quiet.** At 100 in
1,100 the share is **zero**: the planted locus now scores −2.48 a read.

**What absorbs them is the noise parameters, and not — as it first appears — the genotype
frequencies.** The fit hands the planted shape a genotype at *both* densities, a (0, +6) pair at
0.00498 with 5 planted and 0.0910 with 100, each of them the planted share. So of the 24.7 nats the
planted locus gains between the two runs, its genotype's frequency supplies 2.9, about one eighth;
the rest comes from the fall-off going 0.043 to 0.467 and the level from 0.10118 to 0.12198 — and
that fall-off is the whole stratum's, fitted over all 1,100 loci. **The contamination has moved
into the parameters instead of standing out from them, which is the damage spec §4.6 describes, and
it is why Milestone H measures what such loci cost rather than trusting this number to find them.**

## What the record carries that the resolved slippage does not

`StratumFit::fitted_over` and `shares_fitted_over` **name the stratum itself where a number is its
own**, where `StratumSlippage`'s two lists are empty in that case. The resolved answer is read with
the stratum in hand; a record is read beside several hundred others, where an empty list leaves a
reader to work out whether it means *its own* or *not recorded*. Neither list is ever empty in a
record, and the three claims a stratum can make stay distinguishable:

| what happened | `fitted_over` | `shares_fitted_over` | provenance |
|---|---|---|---|
| fitted from its own loci | itself | itself | `FittedHere` |
| kept its level, borrowed its shares | itself | the lender | `FittedHere` |
| borrowed everything | the lenders | the lenders | `Borrowed` |
| merged | the whole set | the whole set | `FittedHere` |

The fourth row is why `merges` is an argument: a pooling is the one thing a resolved slippage
cannot carry, since what it holds is one stratum's model and not the set that model was measured
over.

## Three counters that existed with nowhere to report

`every_climb_settled`, `FitTermination` and the third state of the shares floor were all built in
E2 and E3 and read by nothing. They are now fields on the summary, and none is in the
architecture's sketch of the type, which the doc comments say:

- `strata_with_unsettled_climbs` — fits whose climb over the genotype frequencies ran out of
  passes. Not an error: the surface is concave, so a climb that ran out did not find a wrong
  summit. Expected to be most strata on real data.
- `strata_with_unsettled_searches` — fits whose **outer** search ran out of sweeps. A different
  question, and the one `fitting/` says is the consumer's business: that search has no concavity
  proof at all, so whether a level is where a search settled or where it was stopped is not
  derivable from anything else in the emitted parameters. The review found this one dropped.
- `strata_with_unmeasured_shares` — strata reporting the two shares they measured on fewer moved
  reads than `MIN_SLIPPED_READS_TO_FIT_SHARES` asks, because nobody in their period had enough to
  lend. Spec §4.5 expects every stratum at the bottom of the repeat range to be one of these.
  Derived through `StratumSlippage::keeps_unmeasured_shares()`.

## Where the code went its own way, and why

- **`worst_start_disagreement` names the widest spread of any fit, not the widest among those over
  the limit** as the architecture sketches. `fit_slippage_by_stratum` and `merge_until_monotone`
  both raise `SlippageNotIdentified` on the first fit whose starts land further apart than
  `START_AGREEMENT_LIMIT`, so a sample that reached this summary through this module's own walk has
  none above it — a field restricted to crossers would be permanently empty. A run at 1.002 and a
  run at 1.059 both pass, and only one is comfortable.
- **The guard-share counters are folded over the strata thick enough to speak**, not over every
  stratum. The guard share is a ratio over the reads that moved, so a one-locus stratum with one
  moved read reports 1.0, and most of a genome's 338 strata per read group are thin: a maximum over
  all of them is a near-certain 1.0 from a stratum nobody would act on. Thickness and not *fitted*,
  because the stratum that most needs flagging — one no read of which sits on the grid — is exactly
  the one `fit_slippage` refuses.
- **`low_slippage_substitution` is the rate of the least-slippery stratum that fitted its own
  level**, and not a pool over several. Spec §4.5 says the two paths' rates must meet *where a
  stratum barely slips*; a pool needs a threshold on the level that nothing in the design fixes,
  and the comparison G5 makes is per stratum. A stratum that borrowed its level is not eligible:
  its level describes a neighbour's tracts.
- **`assemble_sample_parameters` takes the resolved map and the merge sets separately**, and where
  both answer for a stratum it checks that they agree about the model — which is what catches a
  caller handing in a resolved map built before the pooling.

## ⛦ Two things for the checkpoint

**Which order borrowing and merging run in.** Unchanged from E4's report and still F1's to wire.
Measured then: borrow-then-merge reports a borrowed stratum at the geometric mean of its lenders'
*pre-merge* levels while both lenders have since moved; merge-then-borrow is flat. The types allow
only one composition anyway — `merge_until_monotone` takes fits and returns fits, `resolve_slippage`
takes fits and returns something else — so the recommendation is to merge first.

**A stratum with no substitution rate aborts the sample, and the case is reachable.** A stratum
every locus of which is witnessed only by reads showing the tract entirely deleted files its shapes
and compares no bases, so `substitution_rate_of` answers `None` — and `StratumFit::substitution` is
not an `Option`, which E1's own doc records as a case the design has not ruled on. E5 fails loudly
rather than inventing a zero that would later be compared against the SNP/indel path's rate as
though it had been measured. Deciding it is the owner's; the recommendation is below.

## Tests

Eighteen. Seven on the diagnostic and eleven on the record and the summary.

| test | what it pins |
|---|---|
| `a_stratum_its_own_fit_explains_reports_nothing_unexplained` | the control: exactly zero, not nearly zero |
| `loci_the_fitted_model_cannot_explain_are_seen_where_the_guard_share_is_blind` | 5/1,005 seen, guard share exactly 0 |
| `a_class_of_loci_the_fit_absorbs_stops_being_unexplained` | the same locus twenty times over is absorbed, and the fall-off is what absorbs it |
| `the_unexplained_share_is_per_read_so_depth_alone_never_flags_a_locus` | all 1,000 loci below the limit on the total, none per read |
| `a_locus_is_divided_by_the_reads_the_model_scored_and_not_by_its_depth` | a locus mixing guard reads with scored ones; over its depth it would clear the limit |
| `a_locus_the_model_scores_no_read_of_counts_in_the_denominator_alone` | 5/1,020 and not 5/1,005; an empty table is 0 |
| `a_stratum_the_fits_frequencies_cannot_reach_is_wholly_unexplained_and_is_named` | a share of 1.0, and the summary naming the worse of two strata |
| `a_record_names_the_strata_its_numbers_came_from` | itself, its lenders, 78 genotypes, bases as the rate's warrant, and a whole borrow not counted as a borrowed share |
| `a_stratum_that_kept_its_level_and_borrowed_its_shares_says_so_in_two_lists` | the third claim, and the counter that keeps it apart |
| `a_merge_is_named_once_and_every_stratum_in_it_names_the_set` | one entry per merge, and 2,400 loci behind a fit no single table holds |
| `a_merge_at_two_ploidies_is_two_claims` | the ploidy half of the dedup key, and a haploid record's twelve genotypes against a diploid's 78 |
| `the_summary_reports_the_least_slippery_fitted_stratums_substitution_rate` | three strata with the quiet one in the middle |
| `the_summary_holds_the_loci_behind_the_thinnest_and_the_thickest_fit` | 1,200 and 5,000, with the borrower behind neither |
| `the_summary_counts_a_climb_that_ran_out_and_shares_nobody_could_lend` | two of three climbs, one capped search, one 1.5-fold start spread named |
| `the_summary_counts_the_strata_this_model_does_not_describe` | the thick ragged stratum counted, the four-locus one and the one exactly at the limit not |
| `two_libraries_get_two_summaries` | nothing crosses between libraries |
| `a_stratum_nothing_resolved_is_refused_rather_than_dropped` | the first panic |
| `a_stratum_with_no_compared_bases_is_refused_rather_than_given_a_zero` | the second |

## What the review changed

Three agents: a mutation run, a numbers check and a design review.

**The mutation run put 47 mutations through the suite and 15 survived, 14 of them provably
changing behaviour** — the same shape as every earlier milestone of this plan: *every fixture was
uniform in the dimension the assertion was about*. Ploidy 2 everywhere, period 2 everywhere, four
reads a locus, `start_spread` 1.0 on every fit, and `fit.loci == table.loci` for every stratum. So
scoring at a hardcoded diploid, dropping the ploidy from the merge dedup key, inverting both
start-spread rules, and reading the stratum's locus count where the fit's was meant all left the
suite green. Two more were reachable-input gaps: no locus anywhere mixed guard reads with
whole-repeat reads, which is the only place the two depths differ, and no test scored a table
against a fit whose mass sits where the reads are not, which is the only place the genotype
*frequencies* do work. Three new tests and eleven added assertions close all fourteen.

**The design review found the arrangement-count defect above** — the one finding that changed a
number rather than a test — **plus the dropped `FitTermination`, the guard counters folded over
thin strata, and a contract that claimed two orderings where the types allow one.** All four are
fixed. Its blocker, the missing substitution rate, is the stop-and-ask above.

**The numbers check found six wrong claims, all mine, and two of them were wrong mechanisms:** why
two locus classes fell below the limit on the undivided total (the two heterozygous classes, and
because a heterozygote spreads its reads over two lengths — not because their genotypes are rare),
and what absorbs planted loci (the fall-off and the level, not the genotype frequencies). Both are
now asserted by the tests rather than only stated. The other four: "one mismatched base in every
read" for one read in four, a genotype ceiling of 91 quoted without its ploidy, a fixture whose
"thickest" alternative was a tie, and an invariant claimed for every run that holds only for this
module's own walk.

## Validation

`cargo fmt --check` and `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib --bins --tests --all-features` green in the container — 3,483 → **3,501** library
tests, 0 failed, 11 ignored, and the 80 integration tests unchanged. The library suite is also
green **natively on the host**, at the same 3,501, which the previous host run predates by four
commits.

**Every mutation the review found surviving was re-run against the fixed tree by hand rather than
taken on the agent's word**, sixteen of them, and all sixteen now fail — including one the first
round of fixes did not close: reading the stratum's own locus count where the fit's was meant
survived because the merge fixture's third stratum held exactly as many loci as the pooled fit's
smaller half.

**Checkpoint E's own condition, which is run by hand:**
`the_search_recovers_a_known_truth_and_no_start_beats_it`, the thirteen-minute accuracy control E2
had invalidated, passes — 589.84 s in the container.
