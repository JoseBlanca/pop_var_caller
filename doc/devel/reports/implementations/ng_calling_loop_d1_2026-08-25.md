# ng calling loop — D1: the driver, and a sample that leaves before the first pass

**Step:** D1 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the table built once,
the outer rounds structurally off.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §2, §4.1, §5,
§5.0, §5.1, §9; [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md) §1, §5, §6.
**Date:** 2026-08-25. **Branch:** `ng-calling-loop`.

---

## 1. What landed

**`SummariseConditionLoop`** — arm A of the step-9 seam, and the first implementation of
`LocusGenotyper`. It is the whole of spec §2's pseudocode: evidence and frozen parameters in,
a `LocusInference` out.

```
candidate alleles := selection's                       (given)
repeat                                  the discovery round — §4.1, Off by default
  repeat                                the slippage round — §5.1, 0 rounds by default
    build the samples × genotypes likelihood table      (once)
    run the frequency loop                              (C1, C2)
  ...
finally
  the final pass                                        (C3b)
```

**The two outer rounds are written as loops whose bodies run once**, rather than left out and
added later — what they buy is a place for
[`calling_bakeoffs.md`](../../ng/impl_plan/calling_bakeoffs.md) to re-enter from. **Each body ends
in an unconditional `break`, and that — not the configuration — is what makes it run once**:
neither loop reads its round's settings, so a token that somehow held a non-default one would
change nothing here. Validation is the second lock rather than the first. **The bodies are not
empty; they hold the whole of the work**, and what a filled-in round adds is a reason to go round
again. `clippy::never_loop` is expected once, with that reason — and measured: turning either
`break` into a `continue` is an `E0384` rather than a second round, because `outcome` is assigned
inside the inner body, and making the binding `mut` to get past that spins forever, because
nothing counts rounds.

**The genotype-likelihood table is built once**, before the frequency loop, and read by every
pass of it. `EmissionCost` — a new instrument on `CallingScratch` — records the table builds, the
row builds, and the emission evaluations those builds were asked for. D2 is what asserts on it.

### 1.1 The ruling B2 and C3b both left open: a set-aside sample gets no row at all

Spec §5.0 says a sample the candidate step ruled uncallable *leaves the loop entirely, before the
first pass*. Two shapes implement that, and B2's M-step recorded the choice as D1's: **give it a
row and skip it everywhere**, or **never give it one**.

**It never gets one.** The scratch is prepared for the samples the locus is called on, and the
rows are the run's sample order with the gaps closed up. What that buys is that *nothing has to be
told to skip anything*: the M-step sums the rows there are, the convergence delta divides by the
chromosomes those rows carry, and the site quality's count axis runs over the same samples. The
alternative needs a "which samples count" list threaded through **five** places — the E-step's
walk, the M-step, the convergence denominator, the site quality's count axis, and the table build,
which must not score the reads of a row nobody will read — and the one that would have been
easiest to forget is the convergence denominator, spec §6's division by cohort chromosomes, which
is load-bearing across the whole cohort range.

Two things follow. `CallingScratch` gains a **row map** (`claim_row_for`,
`run_sample_of_each_row`, `inbreeding_coefficient_by_row`), and the loop's inbreeding coefficients
now live there rather than being passed in — so `run_frequency_loop` lost an argument and the check
that went with it moved onto the scratch, where it fires at the map's first read. And **the final
pass walks the run's samples against a row cursor**, writing `SampleGenotypeCall::Missing` where a
sample has no row, and asking the map at each step whose row this is: a cursor checked only by its
count at the end is satisfied by a *permutation*.

**The row-map check covers two directions that fail differently**, and neither the way the pre-D1
prose said. Measured with the check downgraded: a map one entry **short** panics at the walk's last
read with `index out of bounds` — loud, but naming a slice where the reader needs the cohort. A map
one entry **long** is the silent one: the walk is over the prepared row count, so the surplus is
never read and a scratch claimed for a different locus runs to completion.

**Measured on the fixture:** at a locus with one callable sample and one set aside, the cohort's
expected copies sum to **2.0** — the called sample's two chromosomes — rather than to the run's
four, and the uncallable sample's row is never built (`row_builds` is 1, not 2).

**The review found the join itself untested, and it was the step's one Blocker.** Every set-aside
fixture put the uncallable sample *last*, where a row's index and its sample's index are the same
number — so replacing the table build's `per_sample[run_sample_of_each_row()[row]]` with
`per_sample[row]` left the whole suite green. With the gap **first**, that mutant calls the
surviving sample `0/0` on reads it never saw, copies `[2.0000, 3.9e-6]` against the right
`[0.0077, 1.9923]`: a systematic permutation and no panic. The fixture now exists, and the final
pass no longer re-derives the map — it reads the one the table build read, and asserts row by row
that the sample the map names is the sample the candidate step left callable.

### 1.1a The locus's warrant is derived, not stamped

The first draft passed `Provenance::FittedHere` at the one call site, with no comment. The review
measured what that costs: the fixtures' calibration is `ReadGroupCalibration::defaulted()`, whose
own provenance is `Defaulted`, and the record claimed `FittedHere` — **the exact failure
`LocusInference::weakest_provenance` exists to prevent**, shipped by the field meant to prevent it.

It is now derived: the weakest warrant of the calibrations of the read groups whose reads reached
the locus, folded with `Provenance::weaker_of`. **The ordering that field's doc called undecided is
not** — `parameter_estimation` states the ladder (fitted here, borrowed, supplied, defaulted) and
implements it; that doc comment is corrected here too. What the fold does *not* yet include is the
prior's seed, which carries no provenance at all, and a repeat tract's slippage warrants, which
travel on the scoring contexts step E2 gathers.

### 1.2 Two refusals rather than two approximations

**A repeat tract is refused, and the message says what it is waiting for.** The repeat-tract row
exists and is shipped, but what it takes is a scoring context per `(read group, candidate)` holding
two fitted parameters beside the motif — the stutter model and the STR substitution rate — and one
of the two is neither on `FrozenParameters` nor in `StratumFits`. The pre-pass emits it as a map of its own
(`parameter_estimation::ssr::substitution_rates`), and **gathering the pre-pass's outputs is step
E2's**. Assembling it here would mean inventing that field, the outlier weight's source and the
reachable-length buffer's shape, each of which E2 and E2a settle against real inputs. This is the
one part of D1's contract not delivered, and §4 says what it costs.

**A run that fitted a contamination fraction is refused too.** The mixture has two halves: the
per-read-group fractions, which `FrozenParameters` carries, and the per-locus contaminant allele
frequencies, which come from the loop's own current expected copies and are **step E2a's**.
Scoring against the fractions alone would be a different model rather than a smaller one.

Both follow the rule the two unbuilt loop settings already follow: **refused loudly, never
half-honoured** — and both refuse at the seam's front door, before the worker's scratch is touched.
The release profile aborts on a panic, so one such locus will end a whole cohort run once this is
wired in; it should do that from the arm the run selected rather than from three frames down, and
it should not first leave a shared scratch prepared for a locus nobody scored.

### 1.3 What else the driver assembles

- **The locus's seed concentration**, once per locus, from the run's fitted spectrum and the
  locus's variant class (`calling_priors.md` §2.3). The class is read off the candidate lengths —
  every alternative the reference's own length is a substitution, anything else an insertion or a
  deletion. The two classes take the same seed today, which is exactly why it has a test: nothing
  downstream would notice it being wrong until the split arrives.
- **The error-spread table**, once per locus rather than once per sample: how far an allele's own
  error mass spreads across the locus's others is a property of the candidate sequences and the
  genotype, not of anything a sample showed. It is `genotypes × alleles` and lives on the scratch.
- **The genotype prior as a value**, beside the emission model. Both are seams the design exists
  to compare across — the two priors disagree by 11 points of genotype accuracy on GIAB at 5×
  (`calling_priors.md` §2.2) — so an arm holds both and `name()` says which arm it is.

## 2. Deviations from the plan

- **The STR half of the table build is not here** (§1.2). The plan's D1 entry names "contexts per
  `(read group, candidate)` looked up from `StratumFits`" as this step's; the lookup that exists
  covers three of the four numbers a context needs.
- **`run_frequency_loop` lost its `inbreeding_by_sample` argument** and `summarise_final_pass`
  changed its walk, both consequences of §1.1's ruling. Four `#[should_panic]` tests changed the
  message they expect, because the check moved rather than went.
- **The arm is `pub`, not `pub(crate)`.** A run selects it, and the seam it implements is public.

## 3. Tests

**Nineteen**, all on the driver: 4,672 → 4,691 on the library target. **Eleven of the nineteen are
the review's**, and §6 says what each of those closed.

| test | what it pins |
|---|---|
| `the_driver_calls_genotypes_from_reads_and_builds_the_table_once` | reads in, `1/1` and `0/0` out, copies `[2.0151, 1.9849]`, and `EmissionCost { table_builds: 1, row_builds: 2, emission_evaluations: 4 }` at two passes |
| `a_sample_the_candidate_step_ruled_uncallable_gets_no_row_and_no_vote` | one row for one callable sample; its 40 loud reads never scored; the expected copies sum to 2.0, not 4.0 |
| `a_locus_where_every_sample_was_ruled_uncallable_is_refused` | the case selection's truncation ruling does not cover |
| `a_repeat_tract_is_refused_by_the_driver_rather_than_scored_against_invented_parameters` | §1.2's first refusal, by message |
| `a_run_with_a_fitted_contamination_fraction_is_refused_until_its_other_half_exists` | §1.2's second |
| `a_scratch_prepared_for_more_rows_than_the_locus_can_call_is_refused` | the row-by-row join, at a locus whose rows and samples disagree |
| `the_shipped_arm_names_itself_and_is_object_safe` | a result that cannot name its arm is not auditable |
| `the_variant_class_is_read_off_the_candidate_lengths` | a number that moves nothing today and would be silently wrong when the split arrives |

**The review's eleven**, each closing something measured to survive the first suite:

| test | what it pins |
|---|---|
| `call_locus_scores_each_row_against_the_sample_that_claimed_it_when_the_gap_comes_first` | the Blocker: the uncallable sample **first**, so a row's index and its sample's differ |
| `call_locus_reports_the_weakest_warrant_of_the_parameters_that_reached_the_locus` | a defaulted calibration may not come back as `FittedHere` |
| `a_locus_no_read_reached_has_no_weaker_warrant_to_report` | the fold's identity, and that it is an identity rather than a claim |
| `call_locus_honours_the_runs_pass_cap` | the run's config reaching the loop: `passes = 1, converged = false` against the default's `2, true` |
| `call_locus_claims_each_row_with_its_own_samples_inbreeding_coefficient` | F = 0 beside F = 0.9 on identical reads — qualities 32.32 and 20.75 |
| `emission_evaluations_sum_each_samples_own_observation_count` | one observation beside two: 6, where charging the first row's count for every row gives 4 |
| `emission_evaluations_charge_the_partial_observations_too` | a partial read is an emission the builder was asked for |
| `the_table_is_built_once_at_a_locus_that_takes_four_passes` | one build at four passes, three samples over three alleles |
| `call_locus_calls_a_cohort_of_one` | the hardest corner of the range, end to end |
| `is_callable_rules_no_sample_out_on_a_repeat_tract` | spec §5.0.1's ruling, as a unit, because the tract path is refused before the count is observable |
| `rows_claimed_past_the_end_of_the_run_are_refused` | the one disagreement the per-sample join cannot see |

## 4. What this step owes, and to whom

- **The repeat-tract table build is E2's to unblock and E3's to use.** Until then a tract panics
  with a message naming the missing input. E3's tract half was already fixture-supplied for its
  candidates; this adds that its *parameters* are E2's too.
- **The contamination wiring is E2a's**, on both paths.
- **The locus's warrant is derived from the calibrations and from nothing else.** The prior's
  fitted spectrum carries no provenance at all, and a repeat tract's slippage warrants travel on
  the scoring contexts E2 gathers — so a locus's `weakest_provenance` is today the weakest
  *calibration* that reached it. Whoever gives `SpectrumSeed` a warrant should fold it in here.
- **⚑ A locus where every sample was ruled uncallable is refused, and the owner may want another
  answer.** Candidate selection prefers cutting an allele over refusing a locus, on the ground
  that most samples stay callable; where none does, that argument has nothing to rest on. Emitting
  the locus with every call missing and no expected copies is the other defensible answer. Nothing
  can reach either today — the loop is not wired into a run — so the refusal is a placeholder for
  a decision, not a decision.

## 5. What the review found

**Three agents in worktrees: one on tests and mutation, one on six craft checklists, one
re-deriving the diff's own claims.** Verdict: **1 Blocker, 8 Majors, and 10 of 39 quantitative or
mechanism claims wrong**. Every finding was applied.

**The Blocker was the join this step exists to get right** (§1.1). Every fixture put the
uncallable sample last, so the row index and the sample index coincided and a table filled by row
passed all of them.

**Of the eight Majors, three were about telling one number from another.** Three
`sample_count()` accessors, two meaning the run and — after this step — one meaning the callable
subset, appearing forty lines apart in one function as deliberately different numbers: the
scratch's is now `row_count()`. The final pass re-deriving the row map instead of reading it. And
the warrant stamped rather than derived (§1.1a).

**Two more were tests that could not fail.** The driver's `config` argument was unpinned — dropping
it turns a capped locus (`passes = 1, converged = false`) into a converged one — and every driver
fixture was outbred, so the per-sample inbreeding lookup was invisible; a fixture at F = 0 beside
one at F = 0.9 separates their qualities by 32.32 against 20.75 on identical reads. The emission
counter was asserted only on the fixture shape its own documentation names as the one that hides
the bug: two samples of one observation each. At one observation beside two, `Σ_s obs_s ×
candidates` is 6 where a version charging the first row's count for every row reports 4.

**Ten claims were wrong and all ten were mechanisms, not numbers** — every counted figure in the
first draft checked out. Four doc comments told the same stale story about what a mis-sized row map
costs, and the two directions were each other's answer, swapped (§1.1). The claim that the outer
loops' bodies "cannot be reached" was wrong twice over: they run at every locus and hold all the
work, and what makes them run *once* is the `break` rather than validation. And "prepared first,
claimed second … would hand this locus the last one's samples" describes a carry-over that cannot
happen: the wrong order is caught immediately, with zero claimed rows.

## 6. Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4691 passed; 0 failed; 14 ignored`. Before D1: **4,672**.
- `cargo test --release --lib ng::calling --all-features` — `645 passed; 0 failed; 3 ignored`.
  Before D1: **626**.
- **The release-held checks: D1 adds six.** Downgraded all six to `debug_assert` together and
  re-ran under `--release`: `638 passed; 7 failed`, and every one of the six is reached — the
  repeat tract, the contaminated run, the locus with nobody to call, the per-sample row join, rows
  claimed past the end of the run, and the scratch's row map being one entry per prepared row (two
  tests, one check). **Two more were downgraded to `debug_assert` on the review's finding**: the
  restatements of the first two inside `build_genotype_likelihood_table`, which no test can reach
  because the front door refuses first — a release check the suite cannot keep honest is not one.
