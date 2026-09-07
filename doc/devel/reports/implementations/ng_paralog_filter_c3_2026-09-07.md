# The hidden-duplication filter — C3: scoring the parked records, and where the cut falls

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone C, step C3
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.3, §3.7, §6 trap 4
**Branch:** `ng-paralog-filter`

## The answer

**Pass two exists: it reads the spill once, scores every parked record, fits how common hidden
duplications are in this run, and turns the operator's target false-discovery rate into a
likelihood ratio to cut at.** What it hands pass three is spec §3.7's `ParalogVerdicts` — the
calibration, and one ratio per record in spill order.

**The whole step is about one value.** A record that no sample could speak for is *unscored*: it
is kept, never flagged, and it must not enter the fit. The scorer's own answer there is a ratio of
`0.0`, which is also the number it returns when the two stories are exactly balanced — so passing
it straight through would turn "we could not weigh this" into "this is not a duplication". A run
whose coverage models were all rejected would then fit a duplication rate out of nothing at all and
write a file that looks right (spec §6 trap 4). Two things stand against that here, and they are
both structural rather than a rule someone has to remember:

- **`ParalogScoringContext::score` returns `NaN` when the scorer weighed nothing**, reading the
  scorer's own `samples_used` count rather than re-deriving the condition;
- **the ratio folded into the histogram and the ratio recorded for pass three are the same
  binding.** The histogram drops a `NaN` and the vector keeps it — that is the difference between
  "did not enter the fit" and "kept, never flagged" — but no later edit can fold one value and
  record another.

## What was added

| file | what it is |
|---|---|
| [`pass_two.rs`](../../../src/ng/run/paralog_filter/pass_two.rs) | `ParalogVerdicts`, `PassTwoError`, `score_the_parked_records_and_resolve_the_cut`, and `calibrate_from_the_ratio_histogram` |
| [`pass_two/tests.rs`](../../../src/ng/run/paralog_filter/pass_two/tests.rs) | 13 tests |
| [`scoring_context.rs`](../../../src/ng/run/paralog_filter/scoring_context.rs) | `ParalogScoringContext::score`, and four accessors narrowed to private |

## Assumptions and recorded deviations

### 1. `score` is the only way in, and four accessors became private

C1's review filed this as Mi11 at Medium confidence and deferred it to this step, because it is a
recommendation about C3's shape and C3 did not exist. It holds. `score_locus_for_paralogy` takes
three arguments — the observations, the σ₀ slice, and the precomputed tables — and answers a
mismatch between them with a **neutral score rather than an error** (spec §6 trap 3), so a caller
that paired the wrong two would produce a run that quietly flags nothing. Handing all three out
separately is what made that pairing spellable. `score` builds all three from one place, and
`observation_of`, `observations_of`, `single_copy_depth_sd` and `score_tables` are now private to
the module. This is production's shape: it keeps the same three behind one `score` method.

`sample_count`, `how_many_samples_have_a_coverage_model` and `why_no_model` stay public — they are
what spec §3.5's run-report line reads, and C4 writes that line.

### 2. Unscored is "the scorer weighed nothing", not "no observation was built"

Production screens on whether any observation was constructed
([`score_spilled_locus`, `calibrate.rs:236`](../../../src/var_calling/paralog_filter/calibrate.rs)).
That is one step short: the scorer *also* drops a sample whose one-copy depth spread σ₀ is not
positive and finite, so a record can hand it six observations and have it weigh none of them, and
production folds a `0.0` for that record. Reading `samples_used` — the scorer's own count of what
it weighed — closes the gap and is what spec §3.2 actually says ("a record with no usable sample is
unscored"). The case is reachable and is tested:
`a_record_whose_samples_have_models_but_no_usable_spread_is_unscored`.

**This is the one place where ng does not do exactly what production does**, and it is deliberate.
It costs nothing at the shipped configuration, where σ₀ comes from the fit rather than an override.

### 3. The fallback lives in pass two, not beside the copied statistics

An estimate that did not settle is replaced by the documented rate rather than used, because an
unconverged iterate is not distinguishable from a real estimate by its value alone. Production does
this in `calibrate_from_histogram`, which sits in its file *below* four items ng deliberately does
not copy — so it is not in the span `copy_fidelity.rs` guards, and the file it would join is
compared byte for byte against production's and may not gain a line. It is four lines, and it lives
beside its only caller with the oracle beside it:
`the_fallback_and_the_cut_agree_with_productions_bit_for_bit` compares ng's against production's
across empty, small and 500-ratio histograms, five targets and both convergence outcomes.

### 4. The warning is words, not a print

The plan asks pass two to "warn and fall back". It falls back, and
`ParalogVerdicts::why_the_paralog_rate_is_not_fitted` is the sentence — naming the fallback rate
and how many records there were to fit from. **Printing it is C4's**, because the run report is
where this project puts operator-facing words and pass two has no other output. Spec §3.5 already
asks the run report for π, the cut and whether the estimate converged.

### 5. `records_scored` is one field beyond spec §3.7

The number of ratios that were folded. A run over a million records whose coverage models were all
rejected produces a million `NaN`s, falls back, and writes a VCF indistinguishable from a run where
the filter worked; this is the number that tells the two apart, and it is what the trap-4 test
asserts against the finite count.

### 6. A record that is not the run's cohort stops the pass

Production returns `NaN` and carries on. C1 made the mismatch an error because it is a wiring
error between the sink that filled the spill and the context that reads it, not a property of the
data — and because the scorer answers it with a neutral score, so carrying on means every short
record folds a plausible zero. Pass two propagates it.

## Tests

**13**, all in [`pass_two/tests.rs`](../../../src/ng/run/paralog_filter/pass_two/tests.rs).

Four are about the one value the step exists for — that an unscored record carries `NaN` and never
a zero: `a_record_no_sample_can_speak_for_is_unscored`,
`a_record_whose_samples_have_models_but_no_usable_spread_is_unscored`,
`the_histogram_folds_exactly_the_finite_ratios` (the plan's own test — four records, one of them
unscorable and **sitting in the middle of the file**, because at either end a pass that stopped
early or started late gives the same counts), and
`an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing`.

**Which dimensions the fixtures vary.** Whether a record can be scored; where the unscorable one
sits in the file; repeat tract against generic locus — including a spill of *nothing but* tracts,
so a pass that silently dropped them would fit on nothing rather than merely fitting on less; the
cohort size, at one sample and at six; and whether the rate estimate settles. They hold GC content
constant, at one bin per histogram, and the module doc says so: GC enters through the coverage
model and nothing in this pass reads it, and C1's own review found that a one-GC-bin fixture is
exactly what hides a wrong GC argument — in `scoring_context/tests.rs`, where the GC-varying
fixture now lives.

**Two tests were rewritten because their premise was wrong, and the measurement is what said so.**
The first draft asserted that a target false-discovery rate of zero is unreachable and flags
nothing. On a fixture of twenty records where every sample of a record looks identical, the two
classes separate completely — measured ratios of `+156.26` and `−31.15` and nothing between — so
every target from 0 to 0.5 gives the same cut and removes the same five records. The target was
doing nothing, and a test asserting "some records are flagged" would have passed while proving it.
The fixture is now seven records graded by how many of the six samples look duplicated, none to
all: ratios `−31.15, −3.43, 24.29, 52.01, 79.74, 115.19, 156.26`, and 4, 5 and 7 records removed at
targets of 0, one in a hundred and one in two. Both tests assert the *order* — ratios rising with
the evidence, records removed rising with the target — so neither goes stale if the scorer's
arithmetic changes.

### The mutation ledger

Eight deliberate defects, each rebuilt and re-run against the module's tests in the container, each
file restored afterwards and its checksum compared with the original. **All eight killed.**

| # | the defect | killed by |
|---|---|---|
| 1 | a record with nothing weighed is never called unscored | the three trap-4 tests |
| 2 | an unscored record's ratio is `0.0` instead of `NaN` | the three trap-4 tests |
| 3 | an estimate that did not settle is used rather than replaced | the fallback test, and the parity test |
| 4 | nothing is folded into the histogram | six of the thirteen |
| 5 | only a finite ratio keeps a slot in the vector | the three trap-4 tests |
| 6 | the scored count is the record count | the three trap-4 tests |
| 7 | the cut ignores the operator's target | the parity test |
| 8 | an unfitted rate is never warned about | the two fallback tests |

Defect 7 is worth a line: it changes only `lr_threshold`, the number the header quotes, and **not
which records are removed** — the removal decision reads the curve rather than the threshold. So
none of the eleven behavioural tests can see it, and the only thing that does is the comparison
against production's cut, bit for bit.

## The gate

`main` is red on four checks, so `--all-targets` is unavailable to this branch
(`examples/ng_candidate_selection_probe.rs` does not compile against the current `ClosedLocus`).
What was run, in the container:

| command | at this step | at the merge base `546845b2` |
|---|---|---|
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6595 passed; 0 failed; 15 ignored` | `6582 passed; 0 failed; 15 ignored` |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, three kinds | same 9, same three kinds |
| `cargo fmt --check` | dirty on 9 files, none this step's | the same 9 |

6,595 is 6,582 plus exactly the 13 tests this step adds. The one integration failure —
`ng_calling_loop_calls_genotypes::a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`
— is `main`'s and fails at the merge base too.

**The clippy count is over error *kinds*, not over the three the baseline happened to have.** This
step added a tenth error of a kind the baseline does not carry — `useless use of format!`, in a
test assertion — and a check that grepped for the baseline's three kinds could not have seen it by
construction. That is C2's own reported failure, and the command is the one its review prescribed:

```
cargo clippy --lib --bins --tests --all-features -- -D warnings 2>&1 \
  | grep -E "^error: " | grep -v "could not compile" | sort | uniq -c
```

It is fixed — the assertion now destructures the error rather than formatting it — and the counts
above are from the re-run after the fix.

**`cargo fmt` reformats the nine baseline files**, so they were reverted with `git checkout --`
before staging and the fmt set was re-counted by *unique file*, not by line: 12 files before, 9
after, and the three that went away are this step's.

## What this leaves for C4

- **Pass three**: read the spill again in step with `ratios`, drop or tag, append the two INFO
  fields, write through `write_line`, delete the spill, count.
- **The run report and the header line**, including printing
  `why_the_paralog_rate_is_not_fitted` and the samples-with-a-model count.
- **The two reversals C2 recorded**: `--paralog-fdr` back to spec §3.6's `0.01`, and the refusal
  of a non-zero target deleted.
