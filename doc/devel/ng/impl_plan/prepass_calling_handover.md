# The pre-pass → calling handover — implementation plan

**Status:** draft, 2026-08-27. Scope set by the owner the same day: **the seam only**, and **an
unfitted inbreeding coefficient refuses the run**.

The parameter pre-pass measures what a run's data is like — how often a base is read wrong, how
often a repeat tract slips, how variable the cohort is. The calling loop needs every one of those
numbers, and it reads them from one object, `RunParameters`. **Nothing today builds that object
from what the pre-pass produced.** `RunParameters::assemble` has ten callers and all ten are inside
its own test module; each hands it values written by hand.

This plan builds the one function in between: it takes what the pre-pass returned and returns a
`RunParameters`. It does no input/output, drives no walk, and calls no caller.

---

## Scope

**In:**

1. **The seam** — one function that takes the per-sample generic results, the per-sample repeat-tract
   results, the joint cohort fit, the read groups and the run's ploidy, and returns a
   `RunParameters`.
2. **Two quantities the walk accumulates and the result type drops**, without which the seam cannot
   be written:
   - the per-read-group **minted read errors**, the denominator the read likelihood's error-rate
     scale divides by;
   - the per-stratum **substitution rate**, the rate at which a repeat-tract read shows a wrong
     base rather than a slipped one.
3. **The rule for a sample whose inbreeding coefficient was not fitted**: refuse, naming the sample.

**Out, with owners:**

- **Any driver that walks reads.** No program here runs the pre-pass or the caller. Proving the
  seam on real data is the next plan, and it needs the per-sample walk orchestration that
  `src/ng/run/` does not have yet ([`run_streaming.md`](../spec/run_streaming.md) owns that shape).
- **Wiring `call_locus` into the merge's builder** —
  [`calling_loop.md`](calling_loop.md)'s out-of-scope list, and it additionally waits on the
  unwritten repeat-tract half of candidate selection.
- **Which read groups share a chemistry.** The seam takes the batches it is given. Fitting the
  split is the cohort gather's
  ([`parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md) §7); until it exists the
  honest input is `SequencingBatches::all_together`, which is what the seam will be handed.
- **Relatedness.** Nothing downstream reads it yet.

## What the seam has to bridge

`RunParameters::assemble` takes nine arguments. Six are a straight read or a reshape of something
that already exists:

| what `assemble` wants | where it comes from | the work |
| --- | --- | --- |
| the error rate of each read group | each sample's `GenericSampleParameters::error_rate` | union the per-sample maps |
| each read group's contamination | `JointFit::contamination` | flatten: the per-sample value is **already** `Vec<(ReadGroupId, ContaminationEstimate)>`, so this is not a change of grain |
| which read groups share a chemistry | `SequencingBatches::all_together` | pass through |
| the repeat-tract slippage fits | the `StratumFits` gather | pass through |
| the genotype prior's seed | `RunParameters::seed_from_moments`, on the joint fit's two population moments | pass through |
| the ploidy | the run's configuration | pass through |

Two more are computed during the walk, exposed on the **accumulator**, and dropped by the **result
type** the walk returns. That is the same gap twice:

- `GenericAccumulators::minted_errors()` is public, but `estimate_generic_parameters` consumes the
  accumulator and returns a `GenericSampleParameters` that does not carry the totals.
- `SsrAccumulators::strata()` hands out `(StratumKey, &StratumTable)` pairs and
  `substitution_rate_of(&StratumTable)` already turns one table into an
  `Estimate<ErrorRate>` — but `SsrSampleParameters` carries the *fits* and drops the tables, so the
  two halves never meet.

The ninth is a decision, now taken. Each sample's inbreeding coefficient is an `Option`, because
the fit can fail: a sample can hold too few windows, or the two states the fit separates can fail to
separate. **The seam refuses the run and names the sample.** The pre-pass already
refuses on the same grounds one level up: the `# Errors` note on `GenericAccumulators::estimate`
says of its three inbreeding failures that **none of them has a default**, because a cohort's
diversity divides by `1 − F`, so a coefficient invented rather than measured is amplified rather
than absorbed. The seam repeats that rule at its own edge rather than softening it. The joint fit's
homozygote excess is available and was rejected as the fallback: it is the same fit whose diversity
the coefficient exists to correct, and the census-moments report already carries a warning saying
so.

## Principles (how the order was chosen)

- **The two dropped quantities first.** The seam cannot be written without them, and each is a
  small, separately testable change to a module the seam does not otherwise touch.
- **Types before implementation**, within every milestone (project rule).
- **No new numbers.** Every value the seam produces is one the pre-pass already computed. A step
  that finds itself choosing a constant has found a design question, and it stops and asks.
- **One test per way of being wrong.** The seam's failure mode is a value carried to the wrong key —
  a read group's error rate landing under a neighbour's id, a sample's coefficient landing on
  another sample. Fixtures where every input differs from every other are the only ones that can
  see it.

## Preconditions (already in place)

- `RunParameters::assemble`, with its nine arguments and its own assertions.
- `RunParameters::seed_from_moments`, taking the joint fit's mean alternative-allele frequency and
  heterozygosity.
- The `StratumFits` gather, `SequencingBatches::all_together`, `substitution_rate_of`.
- `ParameterEstimationError`, which already carries the two inbreeding failures the seam's new
  refusal sits beside.

## Branch and merge

Branch `ng-prepass-handover`, cut from `main` at `a6e8472b`, in the existing worktree
`../pop_var_caller-prior-moments`. Sequential — no second worktree.

## The steps

### Milestone A — the minted read errors reach the pre-pass's result

**A1.** Add `minted_errors: BTreeMap<ReadGroupId, MintedReadErrors>` to `GenericSampleParameters`
and fill it where the result is built (`generic/estimate.rs`), from the accumulator's existing
public accessor.

*Fails silently if wrong:* a map filled from the wrong accumulator, or left empty, gives every read
group the defaulted calibration — scale one, every read charged the error floor — and the run still
completes. The test therefore checks a **non-default, per-group-distinct** value survives the trip,
not merely that the map is non-empty.

*Test:* two read groups with different totals, walked, then read off the result and compared
against the accumulator's own.

### Milestone B — every stratum's substitution rate, gathered

**B1.** Add a gather that folds `SsrAccumulators::strata()` through `substitution_rate_of` into a
`BTreeMap<StratumKey, Estimate<ErrorRate>>`, and carry it on `SsrSampleParameters`.

A stratum whose table compared no bases has no rate. `substitution_rate_of` returns `None` there,
and the gather **omits the key** rather than inventing a zero: a zero substitution rate says every
mismatch is a slip, which is the one direction that biases the parameter the repeat-tract design
exists to protect. The doc comment says so, and a test walks a stratum with no compared bases and
asserts the key is absent.

*Test:* three strata, two with different rates and one with none; the map holds two entries, keyed
correctly, with the rates apart.

### Milestone C — the refusal for an unfitted inbreeding coefficient

**C1.** A `ParameterEstimationError::InbreedingNotFittedForSample { sample: String }` — beside the
three inbreeding failures already there — and the sentence that says what to do about it: supply the
coefficient.

*Test:* the error's message names the sample, and a `#[should_panic]`-free path returns it as an
`Err` rather than defaulting.

### Milestone D — the seam

**D1 (types).** `RunParameters::from_prepass(...) -> Result<RunParameters, ParameterEstimationError>`,
taking the per-sample generic results in **sample order**, the per-sample repeat-tract results, the
joint fit, the read groups, the sequencing batches and the ploidy.

Sample order is an argument and not a `BTreeMap` iteration order, because
`inbreeding_coefficient_by_sample` is a `Vec` indexed by the run's sample order and the joint fit
keys its own results by sample **name**. Getting that correspondence wrong is the seam's most
damaging failure and the least visible: every sample gets a coefficient, just not its own.

**D2 (implementation).** The six reads and reshapes, the refusal, and the call to `assemble`.

*Fails silently if wrong:* the sample-order correspondence, and the read-group union. Both get a
fixture in which **no two samples and no two read groups share a value**, so a permuted or
duplicated assignment cannot pass.

*Test:* a three-sample, four-read-group fixture where every number differs; each field of the
returned `RunParameters` is traced back to the input it came from.

### Milestone E — the two ways the run is smallest and largest

The project requires every decision to have an answer at both ends of the cohort and depth ranges.
For a seam that only reshapes, the ends that matter are the **cohort** ones:

- **One sample.** The union over one map, one coefficient, one set of tract fits. The joint fit's
  cohort quantities still exist at one sample; what a one-sample run gives up is stated where the
  seam reads them, not discovered downstream.
- **A thousand samples.** The seam holds every sample's results at once because `assemble` does.
  This milestone measures what that costs — bytes a sample, times a thousand — and writes the number
  down. If it is large, that is a finding for the run driver's plan, not something this seam fixes.

*Test:* the one-sample case as a test; the thousand-sample cost as a measured number in the report.

## Verification summary

| milestone | what proves it |
| --- | --- |
| A | a per-group-distinct total read back off the result |
| B | three strata, one without a rate, keys checked |
| C | the refusal names the sample |
| D | every field traced to its input, on a fixture where no two values are alike |
| E | a one-sample run assembles; the thousand-sample footprint is measured |

Each milestone also runs the whole suite, `cargo fmt --all -- --check` and
`cargo clippy --all-targets --all-features -- -D warnings`.

## Out of scope (next plans)

- **The run driver** — the program that walks every sample, runs both halves of the pre-pass, and
  hands the result to this seam. It is what proves the seam on real data.
- **Wiring `call_locus`** into the merge's builder.
- **Fitting the sequencing-batch split**, the cohort gather's §7.
